"""Aggregate an amino-acid campaign into a cross-target comparison (Phase 8).

For one scan directory (e.g. studies/aa_nh3) writes, under <scan-dir>/aggregate/:
  amino_acid_metrics.csv         one row per amino acid (composition + metrics + tag)
  amino_acid_case_summary.csv    all studies' case summaries concatenated (+ amino_acid)
  SCHEMA.md                      column dictionary for the metrics table
  comparison_summary.md          discriminator groups, prose
  figures/*.png|pdf              ranked bars, target × ΔG heatmap, composition scatter

Usage:
    python scripts/compare_amino_acids.py --scan-dir studies/aa_nh3 \
        --species inputs/amino_acids_species.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from sensitivity_compare import (  # noqa: E402
    build_cross_target_table, classify_discriminators, deltaG_matrix,
    load_campaign,
)

_DESCRIPTIONS = {
    "amino_acid": "Study key (folder name).",
    "target_product": "Cantera target species.",
    "formula": "Molecular formula (CHNOSZ).",
    "n_C": "Carbon atoms.", "n_H": "Hydrogen atoms.", "n_N": "Nitrogen atoms.",
    "n_O": "Oxygen atoms.", "molar_volume_cm3_mol": "Standard partial molar volume (CHNOSZ).",
    "max_stoichiometric_yield_mol": "Limiting-reagent moles formable from the reference feedstock.",
    "inventory_accessible_fraction": "Fraction of the inventory grid where the target is accessible.",
    "min_NH3_accessible": "Smallest NH3 (mol) with accessibility (Batch B only).",
    "min_C2H2_over_HCN_accessible": "Smallest C2H2/HCN ratio with accessibility.",
    "X_eq_at_reference_inventory": "Equilibrium mole fraction at the ΔG-sweep reference inventory.",
    "max_X_eq": "Peak equilibrium mole fraction over the inventory grid.",
    "peak_case_id": "Case id of the peak.",
    "deltaG_positive_crossing_kJ_mol": "ΔG offset (+) where accessibility is lost (None if none in range).",
    "deltaG_negative_crossing_kJ_mol": "ΔG offset (-) where accessibility is lost (None if none in range).",
    "robust_to_pm20": "Accessible across ±20 kJ/mol of Gibbs uncertainty.",
    "robust_to_pm40": "Accessible across ±40 kJ/mol of Gibbs uncertainty.",
    "n_failed": "Solver failures.", "n_suspect_balance": "Cases flagged for element-balance spread.",
    "discriminator": "not_accessible_in_batch | energetically_fragile | inventory_gated | robust_accessible.",
}


def _save(fig, out_dir: Path, stem: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    for fmt in ("png", "pdf"):
        fig.savefig(out_dir / f"{stem}.{fmt}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def _bar(table, value_col, title, xlabel, out_dir, stem):
    d = table.dropna(subset=[value_col]).sort_values(value_col)
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(7, max(4, 0.35 * len(d))))
    ax.barh(d["amino_acid"], d[value_col], color="steelblue")
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    _save(fig, out_dir, stem)


def _heatmap(campaign, out_dir):
    mat = deltaG_matrix(campaign)
    if mat.empty:
        return
    fig, ax = plt.subplots(figsize=(8, max(4, 0.35 * len(mat))))
    offsets = mat.columns.to_numpy(float)
    im = ax.pcolormesh(np.arange(len(offsets) + 1), np.arange(len(mat.index) + 1),
                       mat.to_numpy(float), cmap="viridis")
    ax.set_xticks(np.arange(len(offsets)) + 0.5)
    ax.set_xticklabels([f"{o:g}" for o in offsets], rotation=90, fontsize=6)
    ax.set_yticks(np.arange(len(mat.index)) + 0.5)
    ax.set_yticklabels(mat.index)
    ax.set_xlabel("Delta G offset (kJ/mol)")
    ax.set_title("Amino acid × ΔG offset — log10 equilibrium mole fraction")
    fig.colorbar(im, ax=ax, label="log10 X_eq")
    _save(fig, out_dir, "deltaG_heatmap")


def _composition_scatter(table, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, ycol, ylabel in [
            (axes[0], "max_X_eq", "max X_eq"),
            (axes[1], "inventory_accessible_fraction", "accessible fraction")]:
        d = table.dropna(subset=["n_C", ycol])
        ax.scatter(d["n_C"], d[ycol], color="darkorange")
        for _, r in d.iterrows():
            ax.annotate(r["amino_acid"], (r["n_C"], r[ycol]), fontsize=6,
                        xytext=(2, 2), textcoords="offset points")
        ax.set_xlabel("number of carbons (n_C)")
        ax.set_ylabel(ylabel)
    fig.suptitle("Composition vs accessibility")
    _save(fig, out_dir, "composition_scatter")


def _write_schema(metrics_df, out_path):
    lines = ["# Cross-target metrics — column dictionary", "",
             "| column | dtype | description |", "|---|---|---|"]
    for col in metrics_df.columns:
        desc = _DESCRIPTIONS.get(col, "(undocumented)").replace("|", "\\|")
        lines.append(f"| {col} | {metrics_df[col].dtype} | {desc} |")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_summary(table, out_path, scan_dir):
    groups = {g: table[table["discriminator"] == g]["amino_acid"].tolist()
              for g in ["not_accessible_in_batch", "energetically_fragile",
                        "inventory_gated", "robust_accessible"]}
    lines = [f"# Cross-target comparison — {scan_dir}", "",
             f"{len(table)} amino acids. Discriminator groups:", ""]
    for g, members in groups.items():
        lines.append(f"## {g} ({len(members)})")
        for aa in members:
            r = table[table["amino_acid"] == aa].iloc[0]
            extra = ""
            if g == "energetically_fragile":
                extra = (f" — ΔG crossings −{r['deltaG_negative_crossing_kJ_mol']} / "
                         f"+{r['deltaG_positive_crossing_kJ_mol']} kJ/mol")
            elif g == "inventory_gated":
                extra = (f" — accessible {r['inventory_accessible_fraction']:.0%}, "
                         f"min C2H2/HCN {r['min_C2H2_over_HCN_accessible']}")
            lines.append(f"- {aa}{extra}")
        lines.append("")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def aggregate(scan_dir: Path, species_csv: Path, out_dir: Path) -> pd.DataFrame:
    campaign = load_campaign(scan_dir)
    if not campaign:
        raise FileNotFoundError(f"No case summaries under {scan_dir}/*/results/.")
    table = classify_discriminators(build_cross_target_table(campaign, species_csv))
    out_dir.mkdir(parents=True, exist_ok=True)

    table.to_csv(out_dir / "amino_acid_metrics.csv", index=False)
    concat = pd.concat([cs.assign(amino_acid=k) for k, cs in campaign.items()], ignore_index=True)
    concat.to_csv(out_dir / "amino_acid_case_summary.csv", index=False)
    _write_schema(table, out_dir / "SCHEMA.md")
    _write_summary(table, out_dir / "comparison_summary.md", scan_dir.name)

    figs = out_dir / "figures"
    _bar(table, "inventory_accessible_fraction", "Ranked inventory accessible fraction",
         "accessible fraction", figs, "ranked_accessible_fraction")
    _bar(table, "X_eq_at_reference_inventory", "Ranked X_eq at reference inventory",
         "X_eq (reference inventory)", figs, "ranked_reference_X_eq")
    _heatmap(campaign, figs)
    _composition_scatter(table, figs)
    return table


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-dir", required=True)
    parser.add_argument("--species", default="inputs/amino_acids_species.csv")
    parser.add_argument("--out", default=None, help="Default: <scan-dir>/aggregate")
    args = parser.parse_args(argv)

    scan_dir = PROJECT_ROOT / args.scan_dir
    species_csv = PROJECT_ROOT / args.species
    out_dir = Path(args.out) if args.out else scan_dir / "aggregate"
    if not out_dir.is_absolute():
        out_dir = PROJECT_ROOT / out_dir

    try:
        table = aggregate(scan_dir, species_csv, out_dir)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    counts = table["discriminator"].value_counts().to_dict()
    print(f"Aggregated {len(table)} amino acids → {out_dir}")
    print(f"Discriminator groups: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
