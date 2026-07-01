"""Cross-target aggregation over an amino-acid campaign (Phase 8).

Loads every study's ``sensitivity_case_summary.csv`` from a scan directory and builds
a one-row-per-amino-acid table of composition + accessibility metrics, reusing
``sensitivity_summary.compute_sensitivity_metrics`` so the numbers match the
per-study outputs. Pure pandas/numpy — no Cantera. Single-product only.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

from formula_tools import parse_formula
from sensitivity_summary import compute_sensitivity_metrics

METRIC_COLUMNS = [
    "amino_acid", "target_product", "formula", "n_C", "n_H", "n_N", "n_O",
    "molar_volume_cm3_mol", "max_stoichiometric_yield_mol",
    "inventory_accessible_fraction", "min_NH3_accessible",
    "min_C2H2_over_HCN_accessible", "X_eq_at_reference_inventory",
    "accessible_at_zero_offset", "max_X_eq",
    "peak_case_id", "deltaG_positive_crossing_kJ_mol", "deltaG_negative_crossing_kJ_mol",
    "robust_to_pm20", "robust_to_pm40", "n_failed", "n_suspect_balance", "discriminator",
]


def load_campaign(scan_dir: str | Path) -> Dict[str, pd.DataFrame]:
    """Return ``{amino_acid_key: case_summary_df}`` for a scan directory."""
    scan_dir = Path(scan_dir)
    campaign = {}
    for summary in sorted(scan_dir.glob("*/results/sensitivity_case_summary.csv")):
        key = summary.parent.parent.name
        campaign[key] = pd.read_csv(summary)
    return campaign


def _max_stoichiometric_yield(composition: Dict[str, int], ref_inventory: Dict[str, float]) -> float:
    """Limiting-reagent moles of the target from the available feedstock atoms.

    Atom inventory at the reference point: C from HCN(+1)/C2H2(+2); N from HCN(+1)/
    NH3(+1); O from water; H from water(+2)/HCN(+1)/C2H2(+2)/NH3(+3).
    """
    available = {
        "C": ref_inventory["HCN"] * 1 + ref_inventory["C2H2"] * 2,
        "N": ref_inventory["HCN"] * 1 + ref_inventory["NH3"] * 1,
        "O": ref_inventory["H2O"] * 1,
        "H": ref_inventory["H2O"] * 2 + ref_inventory["HCN"] * 1
             + ref_inventory["C2H2"] * 2 + ref_inventory["NH3"] * 3,
    }
    limits = [available.get(el, 0.0) / n for el, n in composition.items() if n > 0]
    return float(min(limits)) if limits else math.nan


def _reference_inventory(case_summary: pd.DataFrame) -> Dict[str, float]:
    """Feedstock at the ΔG-sweep reference point (auto-adapts to a no-NH3 batch)."""
    dg = case_summary[case_summary["substudy_id"] == "deltaG_sweep"]
    r = dg.iloc[0] if len(dg) else None
    g = lambda col: float(r[col]) if (r is not None and col in r and pd.notna(r[col])) else 0.0
    return {"H2O": g("H2O_mol") or 1.0, "HCN": g("HCN_mol"),
            "C2H2": g("C2H2_mol"), "NH3": g("NH3_mol")}


def build_cross_target_table(campaign: Dict[str, pd.DataFrame], species_csv: str | Path,
                             significant_X_threshold: float = 1e-6) -> pd.DataFrame:
    """One row per amino acid: composition + max yield + accessibility metrics."""
    species = pd.read_csv(species_csv).set_index("cantera_name")
    rows = []
    for key, cs in campaign.items():
        m = compute_sensitivity_metrics(cs, significant_X_threshold=significant_X_threshold)
        g, inv, dg = m.get("general", {}), m.get("inventory_landscape", {}), m.get("deltaG_sweep", {})
        target = str(cs["target_product"].iloc[0])
        meta = species.loc[target] if target in species.index else None
        comp = parse_formula(str(meta["formula"])) if meta is not None else {}
        max_yield = _max_stoichiometric_yield(comp, _reference_inventory(cs))
        rows.append({
            "amino_acid": key,
            "target_product": target,
            "formula": str(meta["formula"]) if meta is not None else "",
            "n_C": comp.get("C", 0), "n_H": comp.get("H", 0),
            "n_N": comp.get("N", 0), "n_O": comp.get("O", 0),
            "molar_volume_cm3_mol": float(meta["molar_volume_cm3_mol"]) if meta is not None else math.nan,
            "max_stoichiometric_yield_mol": max_yield,
            "inventory_accessible_fraction": inv.get("accessible_area_fraction"),
            "min_NH3_accessible": inv.get("min_NH3_accessible"),
            "min_C2H2_over_HCN_accessible": inv.get("min_C2H2_over_HCN_accessible"),
            "X_eq_at_reference_inventory": dg.get("X_eq_at_zero_offset"),
            "accessible_at_zero_offset": dg.get("accessible_at_zero_offset"),
            "max_X_eq": inv.get("max_X_eq"),
            "peak_case_id": inv.get("peak_case_id"),
            "deltaG_positive_crossing_kJ_mol": dg.get("offset_crossing_positive_kJ"),
            "deltaG_negative_crossing_kJ_mol": dg.get("offset_crossing_negative_kJ"),
            "robust_to_pm20": dg.get("robust_to_plus_minus_20_kJ"),
            "robust_to_pm40": dg.get("robust_to_plus_minus_40_kJ"),
            "n_failed": g.get("failed_cases"),
            "n_suspect_balance": g.get("suspect_balance_cases"),
        })
    table = pd.DataFrame(rows)
    return table.sort_values("amino_acid").reset_index(drop=True) if len(table) else table


def classify_discriminators(table: pd.DataFrame) -> pd.DataFrame:
    """Tag each amino acid (data-driven), accessibility-aware:

    - ``not_accessible_in_batch`` — not accessible at the baseline (ΔG=0) reference,
      or never accessible in the swept inventory.
    - ``energetically_fragile`` — accessible at baseline but with a finite ΔG crossing
      in range.
    - ``inventory_gated`` — accessible at baseline, no crossing, but over a smaller
      fraction of the inventory grid than the most-accessible amino acid.
    - ``robust_accessible`` — accessible at baseline, no crossing, top accessibility.
    """
    t = table.copy()
    acc = pd.to_numeric(t["inventory_accessible_fraction"], errors="coerce")
    max_acc = acc.max()

    def _tag(row):
        a = row.get("inventory_accessible_fraction")
        azo = bool(row.get("accessible_at_zero_offset"))
        if not azo or pd.isna(a) or a == 0:
            return "not_accessible_in_batch"
        fragile = (pd.notna(row.get("deltaG_positive_crossing_kJ_mol"))
                   or pd.notna(row.get("deltaG_negative_crossing_kJ_mol")))
        if fragile:
            return "energetically_fragile"
        if pd.notna(a) and pd.notna(max_acc) and a < max_acc:
            return "inventory_gated"
        return "robust_accessible"

    t["discriminator"] = t.apply(_tag, axis=1)
    return t


def combined_nh3_groups(combined: pd.DataFrame, significant_log10x: float = -6.0) -> Dict[str, str]:
    """Per amino acid: accessible_no_nh3 / nh3_unlocked / not_accessible (combined view)."""
    groups = {}
    for aa, d in combined.groupby("amino_acid"):
        a = d[d["source_batch"] == "A_no_nh3"]["log10_X_eq"]
        b = d[d["source_batch"] == "B_nh3"]["log10_X_eq"]
        at0 = bool(a.notna().any() and (a >= significant_log10x).any())
        with_nh3 = bool(b.notna().any() and (b >= significant_log10x).any())
        groups[aa] = ("accessible_no_nh3" if at0
                      else "nh3_unlocked" if with_nh3 else "not_accessible")
    return groups


def _snap_ratio(values, target: float) -> float:
    """Nearest sampled C2H2/HCN ratio column to ``target``."""
    arr = np.asarray(sorted(set(float(v) for v in values)), dtype=float)
    return float(arr[int(np.argmin(np.abs(arr - target)))])


def _combined_row(key, nh3_frac, source, snap, r):
    """One combined row, with yield-relative-to-HCN (paper units)."""
    n_eq = float(r["n_eq_mol"]) if ("n_eq_mol" in r and pd.notna(r["n_eq_mol"])) else math.nan
    hcn = float(r["HCN_mol"]) if ("HCN_mol" in r and pd.notna(r["HCN_mol"])) else math.nan
    yfrac = (n_eq / hcn) if (hcn and not math.isnan(n_eq)) else math.nan
    return {
        "amino_acid": key, "NH3_frac": nh3_frac, "source_batch": source,
        "ratio_snap": snap,
        "X_eq": r.get("X_eq"), "log10_X_eq": r.get("log10_X_eq"),
        "formation_call": r.get("formation_call"),
        "n_eq_mol": n_eq, "HCN_mol": hcn,
        "C2H2_mol": float(r["C2H2_mol"]) if ("C2H2_mol" in r and pd.notna(r["C2H2_mol"])) else math.nan,
        "yield_fraction_HCN": yfrac,
        "yield_pct_HCN": 100.0 * yfrac if not math.isnan(yfrac) else math.nan,
    }


def build_nh3_combined(campaign_no_nh3: Dict[str, pd.DataFrame],
                       campaign_nh3: Dict[str, pd.DataFrame],
                       ratio: float = 2.1) -> pd.DataFrame:
    """Per amino acid at a fixed C2H2/HCN ratio: the NH3=0 point (excluded batch) +
    the NH3 ≥ 0.01 series (present batch). Tidy long table with yield-relative-to-HCN."""
    rows = []
    for key in sorted(set(campaign_no_nh3) | set(campaign_nh3)):
        if key in campaign_no_nh3:
            inv = campaign_no_nh3[key]
            inv = inv[inv["substudy_id"] == "inventory_landscape"]
            snap = _snap_ratio(inv["C2H2_over_HCN"], ratio)
            for _, r in inv[inv["C2H2_over_HCN"] == snap].iterrows():
                rows.append(_combined_row(key, 0.0, "A_no_nh3", snap, r))
        if key in campaign_nh3:
            inv = campaign_nh3[key]
            inv = inv[inv["substudy_id"] == "inventory_landscape"]
            snap = _snap_ratio(inv["C2H2_over_HCN"], ratio)
            for _, r in inv[inv["C2H2_over_HCN"] == snap].sort_values("NH3_mol").iterrows():
                rows.append(_combined_row(key, float(r["NH3_mol"]), "B_nh3", snap, r))
    cols = ["amino_acid", "NH3_frac", "source_batch", "ratio_snap", "X_eq", "log10_X_eq",
            "formation_call", "n_eq_mol", "HCN_mol", "C2H2_mol",
            "yield_fraction_HCN", "yield_pct_HCN"]
    return pd.DataFrame(rows, columns=cols)


def _interpretation(group: str, robust_pm40) -> str:
    if group == "accessible_no_nh3":
        base = "Equilibrium-accessible even with no ammonia"
    elif group == "nh3_unlocked":
        base = "Inaccessible without ammonia; unlocked by NH3 in the feedstock"
    else:
        base = "Not accessible at the fiducial in the swept range"
    return f"{base}; {'robust' if robust_pm40 else 'sensitive'} to ±40 kJ/mol Gibbs uncertainty."


def assemble_bridge_table(campaign_no_nh3: Dict[str, pd.DataFrame],
                          campaign_nh3: Dict[str, pd.DataFrame],
                          ratio: float = 2.1,
                          significant_X_threshold: float = 1e-6) -> pd.DataFrame:
    """Paper-extension bridge: one row per amino acid comparing the NH3-excluded
    fiducial to the NH3 fiducial series (yield % of HCN), with discriminator + prose."""
    level = math.log10(significant_X_threshold)
    combined = build_nh3_combined(campaign_no_nh3, campaign_nh3, ratio)
    groups = combined_nh3_groups(combined, level)

    def _yld(frame, nh3):
        m = frame[np.isclose(frame["NH3_frac"].astype(float), nh3)]
        v = m["yield_pct_HCN"].iloc[0] if len(m) else None
        return float(v) if (v is not None and pd.notna(v)) else None

    rows = []
    for key in sorted(set(campaign_no_nh3) | set(campaign_nh3)):
        d = combined[combined["amino_acid"] == key]
        a = d[d["source_batch"] == "A_no_nh3"]
        b = d[d["source_batch"] == "B_nh3"]
        sig_b = b[b["log10_X_eq"] >= level]
        dgm = (compute_sensitivity_metrics(campaign_nh3[key],
               significant_X_threshold=significant_X_threshold).get("deltaG_sweep", {})
               if key in campaign_nh3 else {})
        group = groups.get(key)
        y_no = (float(a["yield_pct_HCN"].iloc[0])
                if len(a) and pd.notna(a["yield_pct_HCN"].iloc[0]) else None)
        rows.append({
            "amino_acid": key,
            "accessible_no_nh3": bool((a["log10_X_eq"] >= level).any()),
            "accessible_with_nh3": bool((b["log10_X_eq"] >= level).any()),
            "min_NH3_significant": (float(sig_b["NH3_frac"].min()) if len(sig_b) else None),
            "yield_pct_HCN_no_nh3": y_no,
            "yield_pct_HCN_1pct": _yld(b, 0.01),
            "yield_pct_HCN_5pct": _yld(b, 0.05),
            "yield_pct_HCN_10pct": _yld(b, 0.10),
            "deltaG_positive_crossing_kJ": dgm.get("offset_crossing_positive_kJ"),
            "robust_accessible_pm40": dgm.get("robust_to_plus_minus_40_kJ"),
            "paper_group": group,
            "workflow_interpretation": _interpretation(group, dgm.get("robust_to_plus_minus_40_kJ")),
        })
    return pd.DataFrame(rows)


def deltaG_matrix(campaign: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Rows = amino acid, columns = ΔG offset (kJ/mol), cell = log10 X_eq (ΔG sweep)."""
    series = {}
    for key, cs in campaign.items():
        dg = cs[cs["substudy_id"] == "deltaG_sweep"].sort_values("deltaG_offset_kJ_mol")
        if len(dg):
            series[key] = dg.set_index("deltaG_offset_kJ_mol")["log10_X_eq"]
    if not series:
        return pd.DataFrame()
    return pd.DataFrame(series).T.sort_index()
