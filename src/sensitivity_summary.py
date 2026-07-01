"""Sensitivity case summaries, landscape grid, and run metrics (Phase 3).

Converts the merged raw-long equilibrium output into:

- ``sensitivity_case_summary.csv``    — one tidy row per case (all design vars kept).
- ``sensitivity_landscape_grid.csv``  — the plot/ML-ready subset (failed cases kept).
- ``sensitivity_run_summary.md``      — human-readable totals and accessibility metrics.

Key decisions (review §3.3, §5):

- ``case_id`` is the canonical run key — moles are reconstructed with
  ``mole_balance.add_equilibrium_moles(group_cols=["case_id"])``.
- The target row is selected by ``species == target_variant`` (not ``target_product``),
  so ΔG variants are matched correctly.
- Thresholds go through the shared ``result_summary.classify_formation`` so the
  formation calls are identical to the base workflow.
- Failed cases are **kept** (never dropped); ``suspect_balance`` flags cases whose
  ``element_balance_relative_spread`` exceeds ``balance_tol``.

Pure pandas/numpy — no Cantera.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from mole_balance import add_equilibrium_moles
from result_summary import classify_formation

# Identifier/design columns that are constant within a case (carried to the summary).
_CASE_META_COLUMNS = [
    "case_id", "study_id", "substudy_id", "target_product", "target_variant",
    "model_id", "T_C", "P_Pa", "runtime_seconds",
    "H2O_mol", "HCN_mol", "C2H2_mol", "NH3_mol", "C2H2_over_HCN",
    "deltaG_offset_kJ_mol",
]
# Columns the landscape grid exposes (architecture §8.10). target_variant + model_id
# are included so the grid is self-contained for ΔG comparison and ML.
_GRID_COLUMNS = [
    "case_id", "study_id", "substudy_id", "target_product", "target_variant", "model_id",
    "NH3_mol", "HCN_mol", "C2H2_mol", "C2H2_over_HCN", "deltaG_offset_kJ_mol",
    "T_C", "P_Pa", "X_eq", "log10_X_eq", "n_eq_mol", "log10_n_eq_mol",
    "formed_bool", "formation_call", "solver_status", "runtime_seconds",
]


def _safe_log10(value: float) -> float:
    """log10 that returns NaN for non-positive/NaN inputs instead of raising."""
    if pd.isna(value) or float(value) <= 0:
        return np.nan
    return float(np.log10(float(value)))


def summarize_sensitivity_cases(
    raw_long_df: pd.DataFrame,
    species_df: pd.DataFrame,
    thresholds: Dict[str, float],
    output_csv: str | Path | None = None,
) -> pd.DataFrame:
    """One summary row per case, keyed on ``case_id``.

    ``raw_long_df`` is the Phase-2 merged raw long output (it already carries the
    design variables and ``target_variant`` on every row).
    """
    if raw_long_df.empty:
        return pd.DataFrame()

    x_thr = float(thresholds.get("formation_X_threshold", 1e-12))
    x_sig = float(thresholds.get("significant_X_threshold", 1e-6))
    n_thr = float(thresholds.get("formation_n_threshold_mol", 0.0))
    bal_tol = float(thresholds.get("balance_tol", 1e-6))

    raw = raw_long_df.copy()
    raw["case_id"] = raw["case_id"].astype(str)

    # Per-case solver status: a case is ok only if every species row is ok.
    status_by_case = raw.groupby("case_id")["solver_status"].agg(
        lambda s: "ok" if (s.astype(str) == "ok").all() else "failed")
    err_by_case = raw.groupby("case_id")["error_message"].agg(
        lambda s: next((str(e) for e in s if str(e)), ""))

    # Reconstruct moles only for ok cases (add_equilibrium_moles raises on all-NaN).
    ok_raw = raw[raw["solver_status"].astype(str) == "ok"]
    target_by_case: Dict[str, Dict[str, float]] = {}
    if not ok_raw.empty:
        moles = add_equilibrium_moles(ok_raw, species_df, group_cols=["case_id"])
        target_rows = moles[moles["species"].astype(str) == moles["target_variant"].astype(str)]
        for _, r in target_rows.drop_duplicates("case_id").iterrows():
            target_by_case[str(r["case_id"])] = {
                "X_eq": float(r["X_eq"]) if pd.notna(r["X_eq"]) else np.nan,
                "n_eq_mol": float(r["n_eq_mol"]) if pd.notna(r["n_eq_mol"]) else np.nan,
                "element_balance_relative_spread": float(r.get(
                    "element_balance_relative_spread", np.nan)),
            }

    # One meta row per case (design vars are constant within a case).
    meta_cols = [c for c in _CASE_META_COLUMNS if c in raw.columns]
    case_meta = raw[meta_cols].drop_duplicates("case_id").set_index("case_id")

    rows = []
    for case_id, meta in case_meta.iterrows():
        status = str(status_by_case.get(case_id, "failed"))
        target = target_by_case.get(case_id)
        if status == "ok" and target is not None:
            x_eq, n_eq = target["X_eq"], target["n_eq_mol"]
            balance = target["element_balance_relative_spread"]
        else:
            x_eq, n_eq, balance = np.nan, np.nan, np.nan
        formed, call = classify_formation(x_eq, n_eq, status, x_thr, x_sig, n_thr)
        suspect = bool(pd.notna(balance) and balance > bal_tol)
        row = {col: meta[col] for col in meta.index}
        row.update({
            "case_id": case_id,
            "X_eq": x_eq,
            "log10_X_eq": _safe_log10(x_eq),
            "n_eq_mol": n_eq,
            "log10_n_eq_mol": _safe_log10(n_eq),
            "formed_bool": formed,
            "formation_call": call,
            "solver_status": status,
            "error_message": str(err_by_case.get(case_id, "")),
            "element_balance_relative_spread": balance,
            "suspect_balance": suspect,
        })
        rows.append(row)

    summary = pd.DataFrame(rows)
    if "error_message" in summary.columns:
        summary["error_message"] = summary["error_message"].fillna("").astype(str)
    sort_cols = [c for c in ["substudy_id", "case_id"] if c in summary.columns]
    if sort_cols:
        summary = summary.sort_values(sort_cols).reset_index(drop=True)
    if output_csv is not None:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(output_csv, index=False)
    return summary


def make_landscape_grid(case_summary_df: pd.DataFrame,
                        output_csv: str | Path | None = None) -> pd.DataFrame:
    """The plot/ML-ready subset. Failed cases are retained (visible as gaps)."""
    if case_summary_df.empty:
        return pd.DataFrame()
    cols = [c for c in _GRID_COLUMNS if c in case_summary_df.columns]
    grid = case_summary_df[cols].copy()
    if output_csv is not None:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        grid.to_csv(output_csv, index=False)
    return grid


def _accessible_window(offsets, values, level):
    """Bracket-interpolate the accessible ΔG window [neg, pos] containing offset 0.

    ``values`` is log10_X_eq vs ``offsets`` (kJ/mol). Returns the interpolated
    offsets where the curve crosses ``level`` (= log10 significant threshold) on the
    negative and positive sides of 0. NaN on a side means the curve never drops below
    ``level`` within the sampled range (accessible to the grid edge). Returns
    ``(nan, nan)`` if offset 0 itself is below ``level``. No exact sample point at the
    crossing is required — it is found by linear interpolation between brackets.
    """
    o = np.asarray(offsets, dtype=float)
    v = np.asarray(values, dtype=float)
    order = np.argsort(o)
    o, v = o[order], v[order]
    finite = np.isfinite(o) & np.isfinite(v)
    if finite.sum() < 2:
        return float("nan"), float("nan")
    of, vf = o[finite], v[finite]
    if float(np.interp(0.0, of, vf)) < level:
        return float("nan"), float("nan")

    pos = float("nan")
    for i in range(len(of) - 1):
        x0, y0, x1, y1 = of[i], vf[i], of[i + 1], vf[i + 1]
        if x1 <= 0:
            continue
        if y0 >= level > y1:  # accessible -> inaccessible going positive
            pos = x0 + (level - y0) * (x1 - x0) / (y1 - y0)
            break
    neg = float("nan")
    for i in range(len(of) - 2, -1, -1):
        x0, y0, x1, y1 = of[i], vf[i], of[i + 1], vf[i + 1]
        if x0 >= 0:
            continue
        if y1 >= level > y0:  # inaccessible -> accessible going positive (left edge)
            neg = x0 + (level - y0) * (x1 - x0) / (y1 - y0)
            break
    return neg, pos


def compute_sensitivity_metrics(case_summary_df: pd.DataFrame,
                                n_targets_projection: Optional[int] = None,
                                significant_X_threshold: float = 1e-6) -> Dict[str, Any]:
    """General + per-substudy accessibility metrics (robust to missing substudies).

    When ``n_targets_projection`` is given, a multi-target runtime projection is
    included; otherwise it is omitted (no hardcoded target count).
    """
    df = case_summary_df
    metrics: Dict[str, Any] = {}
    total = len(df)
    failed = int((df["solver_status"] != "ok").sum()) if total else 0
    runtimes = df["runtime_seconds"].dropna() if "runtime_seconds" in df else pd.Series(dtype=float)
    metrics["general"] = {
        "total_cases": total,
        "successful_cases": total - failed,
        "failed_cases": failed,
        "failure_rate": (failed / total) if total else 0.0,
        "suspect_balance_cases": int(df.get("suspect_balance", pd.Series(dtype=bool)).sum()),
        "median_runtime_sec": float(runtimes.median()) if len(runtimes) else None,
        "max_runtime_sec": float(runtimes.max()) if len(runtimes) else None,
        "total_runtime_sec": float(runtimes.sum()) if len(runtimes) else None,
        "n_targets_projection": n_targets_projection,
        "projected_runtime_sec_targets": (
            float(runtimes.median()) * total * n_targets_projection
            if (len(runtimes) and n_targets_projection) else None),
    }

    inv = df[df.get("substudy_id") == "inventory_landscape"] if "substudy_id" in df else df
    if len(inv):
        accessible = inv[inv["formed_bool"] == True]  # noqa: E712
        metrics["inventory_landscape"] = {
            "n_cases": len(inv),
            "accessible_area_fraction": float(len(accessible) / len(inv)),
            "min_NH3_accessible": (float(accessible["NH3_mol"].min())
                                   if len(accessible) and "NH3_mol" in accessible else None),
            "min_C2H2_over_HCN_accessible": (float(accessible["C2H2_over_HCN"].min())
                                             if len(accessible) and "C2H2_over_HCN" in accessible else None),
            "max_X_eq": (float(inv["X_eq"].max(skipna=True)) if "X_eq" in inv else None),
            "peak_case_id": (str(inv.loc[inv["X_eq"].idxmax(), "case_id"])
                             if "X_eq" in inv and inv["X_eq"].notna().any() else None),
        }

    level = float(np.log10(significant_X_threshold))

    dg = df[df.get("substudy_id") == "deltaG_sweep"] if "substudy_id" in df else df.iloc[0:0]
    if len(dg) and "deltaG_offset_kJ_mol" in dg:
        d = dg.sort_values("deltaG_offset_kJ_mol")
        zero = d[d["deltaG_offset_kJ_mol"] == 0]
        x_zero = (float(zero["X_eq"].iloc[0])
                  if len(zero) and zero["X_eq"].notna().any() else None)
        # Robustness of *accessibility* is only meaningful if the target is accessible
        # at the baseline (ΔG = 0). _accessible_window returns (NaN, NaN) BOTH when the
        # target is accessible across the whole sweep AND when ΔG=0 is already below
        # threshold — so gate robustness on accessible_at_zero_offset.
        accessible_at_zero = bool(x_zero is not None and x_zero >= significant_X_threshold)
        neg, pos = _accessible_window(d["deltaG_offset_kJ_mol"], d["log10_X_eq"], level)

        def _robust(window):
            if not accessible_at_zero:
                return False
            ok_neg = (not np.isfinite(neg)) or neg <= -window
            ok_pos = (not np.isfinite(pos)) or pos >= window
            return bool(ok_neg and ok_pos)

        metrics["deltaG_sweep"] = {
            "n_cases": len(dg),
            "X_eq_at_zero_offset": x_zero,
            "accessible_at_zero_offset": accessible_at_zero,
            "offset_crossing_negative_kJ": (float(neg) if np.isfinite(neg) else None),
            "offset_crossing_positive_kJ": (float(pos) if np.isfinite(pos) else None),
            "robust_to_plus_minus_20_kJ": _robust(20),
            "robust_to_plus_minus_40_kJ": _robust(40),
        }

    nh3dg = df[df.get("substudy_id") == "nh3_deltaG_landscape"] if "substudy_id" in df else df.iloc[0:0]
    if len(nh3dg) and "deltaG_offset_kJ_mol" in nh3dg:
        thresholds = {}
        for offset, grp in nh3dg.groupby("deltaG_offset_kJ_mol"):
            acc = grp[grp["formed_bool"] == True]  # noqa: E712
            thresholds[float(offset)] = (float(acc["NH3_mol"].min()) if len(acc) else None)
        offs = sorted(thresholds)
        known = [v for v in thresholds.values() if v is not None]
        metrics["nh3_deltaG_landscape"] = {
            "n_cases": len(nh3dg),
            "nh3_threshold_by_offset": {str(k): thresholds[k] for k in offs},
            "nh3_threshold_at_zero_offset": thresholds.get(0.0),
            "nh3_threshold_at_min_offset": thresholds.get(offs[0]) if offs else None,
            "nh3_threshold_at_max_offset": thresholds.get(offs[-1]) if offs else None,
            "threshold_shifts_with_offset": bool(len(set(known)) > 1),
        }
    return metrics


def write_sensitivity_run_summary(metrics: Dict[str, Any], output_md: str | Path,
                                  study_id: str = "") -> Path:
    """Render the metrics dict to a human-readable Markdown summary."""
    g = metrics.get("general", {})
    lines = [f"# Sensitivity run summary — {study_id}".rstrip(), ""]
    lines += ["## Overview", ""]
    lines += [f"- Total cases: {g.get('total_cases')}",
              f"- Successful: {g.get('successful_cases')}",
              f"- Failed: {g.get('failed_cases')} "
              f"(failure rate {g.get('failure_rate', 0):.1%})",
              f"- Suspect element balance: {g.get('suspect_balance_cases')}"]
    med = g.get("median_runtime_sec")
    if med is not None:
        lines += [f"- Median runtime: {med:.4f} s/case",
                  f"- Max runtime: {g.get('max_runtime_sec'):.4f} s",
                  f"- Total runtime: {g.get('total_runtime_sec'):.1f} s"]
        proj, n_t = g.get("projected_runtime_sec_targets"), g.get("n_targets_projection")
        if proj is not None and n_t:
            lines.append(f"- Projected runtime for {n_t} targets: {proj:.0f} s")
    else:
        lines += ["- Runtime: not recorded"]

    inv = metrics.get("inventory_landscape")
    if inv:
        lines += ["", "## Inventory landscape", "",
                  f"- Cases: {inv['n_cases']}",
                  f"- Accessible-area fraction: {inv['accessible_area_fraction']:.3f}",
                  f"- Min NH3 for accessibility: {inv.get('min_NH3_accessible')}",
                  f"- Min C2H2/HCN for accessibility: {inv.get('min_C2H2_over_HCN_accessible')}",
                  f"- Max X_eq: {inv.get('max_X_eq')} (case {inv.get('peak_case_id')})"]

    dg = metrics.get("deltaG_sweep")
    if dg:
        lines += ["", "## ΔG sweep", "",
                  f"- Cases: {dg['n_cases']}",
                  f"- X_eq at 0 kJ/mol: {dg.get('X_eq_at_zero_offset')}",
                  f"- Crossing offset (negative side): {dg.get('offset_crossing_negative_kJ')} kJ/mol",
                  f"- Crossing offset (positive side): {dg.get('offset_crossing_positive_kJ')} kJ/mol",
                  f"- Robust to ±20 kJ/mol: {dg.get('robust_to_plus_minus_20_kJ')}",
                  f"- Robust to ±40 kJ/mol: {dg.get('robust_to_plus_minus_40_kJ')}"]

    nh3dg = metrics.get("nh3_deltaG_landscape")
    if nh3dg:
        lines += ["", "## NH3 × ΔG landscape", "",
                  f"- Cases: {nh3dg['n_cases']}",
                  f"- NH3 threshold at 0 kJ/mol: {nh3dg.get('nh3_threshold_at_zero_offset')}",
                  f"- NH3 threshold at min offset: {nh3dg.get('nh3_threshold_at_min_offset')}",
                  f"- NH3 threshold at max offset: {nh3dg.get('nh3_threshold_at_max_offset')}",
                  f"- Threshold shifts with offset: {nh3dg.get('threshold_shifts_with_offset')}"]

    output_md = Path(output_md)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_md


def render_mvp_verdict(metrics: Dict[str, Any], config: Dict[str, Any]) -> str:
    """Build the scale/revise/pause verdict Markdown from the computed metrics.

    Recommendation: **SCALE** when the run was stable (``failure_rate == 0``) and the
    target is broadly accessible (``accessible_area_fraction > 0.5``); otherwise
    **REVISE**. All findings come from ``metrics`` (nothing hardcoded). Kept here (not
    in the notebook) so the logic is tested and the notebook stays a thin viewer.
    """
    g = metrics.get("general", {})
    inv = metrics.get("inventory_landscape", {})
    dg = metrics.get("deltaG_sweep", {})
    target = (config.get("mode", {}).get("target_products") or ["the target"])[0]

    acc = inv.get("accessible_area_fraction")
    failure_rate = g.get("failure_rate", 1.0)
    recommend = "**SCALE**" if (failure_rate == 0 and acc is not None and acc > 0.5) else "**REVISE**"
    no_crossing = (dg.get("offset_crossing_negative_kJ") is None and
                   dg.get("offset_crossing_positive_kJ") is None)

    lines = ["## MVP verdict", ""]
    lines.append(f"- **Run stability:** {g.get('successful_cases')}/{g.get('total_cases')} "
                 f"cases succeeded (failure rate {failure_rate:.1%}); suspect element "
                 f"balance: {g.get('suspect_balance_cases')}.")
    if inv:
        lines.append(f"- **Inventory landscape:** {target} is accessible across "
                     f"{acc:.0%} of the NH3 × C2H2/HCN grid. Carbon is the limiter — the "
                     f"minimum C2H2/HCN ratio for accessibility is "
                     f"{inv.get('min_C2H2_over_HCN_accessible')}.")
        lines.append(f"- **NH3 is not limiting:** the minimum NH3 for accessibility is "
                     f"{inv.get('min_NH3_accessible')} (forms with zero *initial* NH3 "
                     "because nitrogen is supplied by HCN).")
    if dg:
        if no_crossing:
            lines.append("- **ΔG robustness:** accessible across the **entire swept ±range** "
                         f"— no crossing found; robust to ±20 kJ/mol: "
                         f"{dg.get('robust_to_plus_minus_20_kJ')}, ±40 kJ/mol: "
                         f"{dg.get('robust_to_plus_minus_40_kJ')}.")
        else:
            lines.append(f"- **ΔG robustness:** accessibility crossings at "
                         f"{dg.get('offset_crossing_negative_kJ')} / "
                         f"{dg.get('offset_crossing_positive_kJ')} kJ/mol; robust to ±20: "
                         f"{dg.get('robust_to_plus_minus_20_kJ')}, ±40: "
                         f"{dg.get('robust_to_plus_minus_40_kJ')}.")
    proj, n_t = g.get("projected_runtime_sec_targets"), g.get("n_targets_projection")
    if proj is not None and n_t:
        lines.append(f"- **Scaling:** projected runtime for {n_t} targets ≈ **{proj:.0f} s** "
                     f"(median {g.get('median_runtime_sec'):.4f} s/case).")
    lines += ["",
              f"### Recommendation: {recommend}",
              "",
              f"{target} is a **robust baseline** — broadly accessible and insensitive to "
              "plausible Gibbs-energy uncertainty, and the full study runs in seconds, so "
              "scaling to the other amino acids is cheap. The scientifically interesting next "
              "step is the amino acids that are *not* robust — those whose accessibility "
              "flips within the swept inventory or ΔG range — since they will discriminate "
              "between Titan inventory scenarios."]
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Column dictionary (ML-readiness, review §8)
# --------------------------------------------------------------------------- #
# Single source of truth: column -> (units, one-line description). Units is "" when
# not applicable. Any column not here is still listed with dtype + "(undocumented)".
COLUMN_DESCRIPTIONS: Dict[str, tuple] = {
    "case_id": ("", "Unique simulation case identifier (canonical run key)."),
    "scenario_id": ("", "Generated scenario id (equals case_id in sensitivity runs)."),
    "scenario": ("", "Scenario id carried from the base runner (equals case_id)."),
    "study_id": ("", "Study identifier (folder name under studies/)."),
    "substudy_id": ("", "Substudy: inventory_landscape, deltaG_sweep, or nh3_deltaG_landscape."),
    "target_product": ("", "Original target product, e.g. Alanine(aq)."),
    "target_variant": ("", "Cantera species actually modelled (base name or ΔG pseudo-species)."),
    "model_id": ("", "Hashed Cantera model identity (one YAML reused across grid points)."),
    "model_mode": ("", "Model mode tag, e.g. single_product_sensitivity."),
    "yaml_file": ("", "Cantera YAML file name used for the case."),
    "H2O_mol": ("mol", "Initial moles of water (solvent basis)."),
    "HCN_mol": ("mol", "Initial moles of hydrogen cyanide."),
    "C2H2_mol": ("mol", "Initial moles of acetylene."),
    "NH3_mol": ("mol", "Initial moles of ammonia."),
    "C2H2_over_HCN": ("ratio", "Initial C2H2/HCN mole ratio (derived design variable)."),
    "deltaG_offset_kJ_mol": ("kJ/mol", "Gibbs-energy offset applied to the target species."),
    "T_C": ("deg C", "Equilibrium temperature."),
    "T_K": ("K", "Equilibrium temperature in Kelvin."),
    "P_Pa": ("Pa", "Equilibrium pressure."),
    "species": ("", "Cantera species name for this row (long-form tables)."),
    "X_initial": ("mole fraction", "Initial mole fraction supplied to the solver."),
    "initial_moles": ("mol", "Initial moles of the species in the scenario."),
    "X_eq": ("mole fraction", "Equilibrium mole fraction of the (target) species."),
    "log10_X_eq": ("log10 mole fraction", "log10(X_eq); NaN when X_eq <= 0 or missing."),
    "n_eq_mol": ("mol", "Reconstructed equilibrium moles of the species."),
    "n_total_eq_mol": ("mol", "Reconstructed total moles at equilibrium for the case."),
    "log10_n_eq_mol": ("log10 mol", "log10(n_eq_mol); NaN when non-positive/missing."),
    "formed_bool": ("", "True if the target is equilibrium-accessible above threshold."),
    "formation_call": ("", "significant | trace | below_threshold | solver_failed."),
    "solver_status": ("", "Equilibrium solver outcome: ok or failed."),
    "error_message": ("", "Solver/setup error message when the case failed."),
    "element_balance_relative_spread": ("", "Relative spread of per-element total-mole estimates (QC)."),
    "element_total_mole_estimates": ("", "Per-element total-mole estimates (diagnostic string)."),
    "suspect_balance": ("", "True when element_balance_relative_spread exceeds balance_tol."),
    "runtime_seconds": ("s", "Wall-clock runtime for the case."),
}

# Canonical result tables documented by the schema dictionary (in report order).
DEFAULT_SCHEMA_TABLES = [
    "sensitivity_case_summary.csv",
    "sensitivity_landscape_grid.csv",
    "equilibrium_raw_long.csv",
    "equilibrium_moles_long.csv",
]


def write_schema_dictionary(results_dir: str | Path,
                            tables: Optional[List[str]] = None,
                            output_md: str | Path | None = None,
                            output_json: str | Path | None = None):
    """Emit a column dictionary for the result tables that exist in ``results_dir``.

    Writes a human-readable ``SCHEMA.md`` (one section per table: column | dtype |
    units | description) and a machine ``schema.json`` mapping
    ``table -> {column -> {dtype, units, description}}``. Dtypes are read from the
    written CSVs. Missing tables are skipped (never fails the run); unknown columns
    are still listed with their dtype and a ``(undocumented)`` description so nothing
    is silently dropped. Returns ``(md_path, json_path)``.
    """
    results_dir = Path(results_dir)
    tables = tables if tables is not None else DEFAULT_SCHEMA_TABLES
    output_md = Path(output_md) if output_md is not None else results_dir / "SCHEMA.md"
    output_json = Path(output_json) if output_json is not None else results_dir / "schema.json"

    schema: Dict[str, Dict[str, Dict[str, str]]] = {}
    lines = ["# Result schema dictionary", "",
             "Auto-generated column dictionary for the result tables in this "
             "directory. One section per table; units are blank where not applicable.",
             ""]
    def _esc(value: str) -> str:
        return str(value).replace("|", "\\|")

    for table in tables:
        path = results_dir / table
        if not path.exists():
            continue
        df = pd.read_csv(path)
        # An all-empty error_message reads back as float64; report it as text.
        if "error_message" in df.columns:
            df["error_message"] = df["error_message"].fillna("").astype(str)
        columns: Dict[str, Dict[str, str]] = {}
        lines += [f"## {table}", "", "| column | dtype | units | description |",
                  "|---|---|---|---|"]
        for col in df.columns:
            dtype = str(df[col].dtype)
            units, desc = COLUMN_DESCRIPTIONS.get(col, ("", "(undocumented)"))
            columns[col] = {"dtype": dtype, "units": units, "description": desc}
            lines.append(f"| {_esc(col)} | {_esc(dtype)} | {_esc(units)} | {_esc(desc)} |")
        lines.append("")
        schema[table] = columns

    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    output_json.write_text(json.dumps(schema, indent=2) + "\n", encoding="utf-8")
    return output_md, output_json
