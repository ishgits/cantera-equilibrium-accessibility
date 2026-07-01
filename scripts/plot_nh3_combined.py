"""Combined NH3 figure across both batches (Phase 9, M&P-style).

At a fixed fiducial C2H2/HCN ratio, stitches the NH3 = 0 point (no-ammonia batch) to
the NH3 ≥ 0.01 series (ammonia batch) for each amino acid, and asks which amino acids
are accessible with **no** NH3 vs only **unlocked by** NH3 ≥ 1% of water — the
comparison in Madan & Pearce (2025).

Writes studies/_nh3_combined/: nh3_combined.csv, an amino-acid × NH3% heatmap of
log10 X_eq (analogous to M&P Fig E1), per-amino-acid line plots, and summary.md.

Usage:
    python scripts/plot_nh3_combined.py --no-nh3 studies/aa_no_nh3 --nh3 studies/aa_nh3 --ratio 2.1
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
    build_nh3_combined, combined_nh3_groups, load_campaign,
)

DEFAULT_OUT = "studies/_nh3_combined"
SIGNIFICANT_LOG10X = -6.0   # log10(1e-6); "accessible" threshold for the verdict


def _save(fig, out_dir: Path, stem: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    for fmt in ("png", "pdf"):
        fig.savefig(out_dir / f"{stem}.{fmt}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def _heatmap(combined: pd.DataFrame, out_dir: Path, value_col: str, stem: str,
             title: str, cbar_label: str):
    mat = combined.pivot_table(index="amino_acid", columns="NH3_frac",
                               values=value_col, aggfunc="first").sort_index()
    if mat.empty:
        return
    fig, ax = plt.subplots(figsize=(9, max(4, 0.4 * len(mat))))
    cols = mat.columns.to_numpy(float) * 100.0   # NH3 as % of water
    im = ax.pcolormesh(np.arange(len(cols) + 1), np.arange(len(mat.index) + 1),
                       mat.to_numpy(float), cmap="viridis")
    ax.set_xticks(np.arange(len(cols)) + 0.5)
    ax.set_xticklabels([f"{c:g}" for c in cols], rotation=90, fontsize=6)
    ax.set_yticks(np.arange(len(mat.index)) + 0.5)
    ax.set_yticklabels(mat.index, fontsize=7)
    ax.set_xlabel("NH3 (% of water)  [0% = no-NH3 batch]")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=cbar_label)
    _save(fig, out_dir, stem)


def _line_plots(combined: pd.DataFrame, out_dir: Path, value_col: str, stem: str,
                ylabel: str, threshold=None):
    aas = sorted(combined["amino_acid"].unique())
    ncol = 3
    nrow = int(np.ceil(len(aas) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 2.6 * nrow), squeeze=False)
    for ax, aa in zip(axes.flat, aas):
        d = combined[combined["amino_acid"] == aa].sort_values("NH3_frac")
        ax.plot(d["NH3_frac"] * 100, d[value_col], marker="o", ms=3)
        if threshold is not None:
            ax.axhline(threshold, ls="--", c="gray", lw=0.8)
        ax.set_title(aa, fontsize=8)
        ax.set_xlabel("NH3 (% water)", fontsize=7)
        ax.set_ylabel(ylabel, fontsize=7)
    for ax in axes.flat[len(aas):]:
        ax.axis("off")
    fig.tight_layout()
    _save(fig, out_dir, stem)


def _write_summary(combined: pd.DataFrame, out_path: Path, ratio: float):
    groups = combined_nh3_groups(combined, SIGNIFICANT_LOG10X)
    accessible_no_nh3 = [aa for aa, g in groups.items() if g == "accessible_no_nh3"]
    unlocked_by_nh3 = [aa for aa, g in groups.items() if g == "nh3_unlocked"]
    never = [aa for aa, g in groups.items() if g == "not_accessible"]
    lines = [f"# Combined NH3 comparison (fiducial C2H2/HCN ≈ {ratio})", "",
             "Accessibility threshold: log10 X_eq ≥ −6 (X_eq ≥ 1e-6).", "",
             f"## Accessible with NO NH3 (Batch A) — {len(accessible_no_nh3)}",
             ", ".join(accessible_no_nh3) or "(none)", "",
             f"## Unlocked only by NH3 ≥ 1% (Batch B) — {len(unlocked_by_nh3)}",
             ", ".join(unlocked_by_nh3) or "(none)", "",
             f"## Not accessible in the swept range — {len(never)}",
             ", ".join(never) or "(none)", "",
             "## Comparison to Madan & Pearce (2025)",
             "M&P find only a few amino acids accessible with no NH3 of any kind; "
             "ammonia in the feedstock unlocks the rest. The split above is the "
             "equilibrium-accessibility analogue at the fiducial ratio."]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-nh3", default="studies/aa_no_nh3")
    parser.add_argument("--nh3", default="studies/aa_nh3")
    parser.add_argument("--ratio", type=float, default=2.1,
                        help="Fiducial C2H2/HCN ratio (default 2.1 = 0.042/0.020).")
    parser.add_argument("--out", default=DEFAULT_OUT)
    args = parser.parse_args(argv)

    camp_a = load_campaign(PROJECT_ROOT / args.no_nh3)
    camp_b = load_campaign(PROJECT_ROOT / args.nh3)
    if not camp_a and not camp_b:
        print("No campaigns found — run both batches first.", file=sys.stderr)
        return 1

    combined = build_nh3_combined(camp_a, camp_b, args.ratio)
    out_dir = PROJECT_ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_dir / "nh3_combined.csv", index=False)
    # Two heatmaps: detectability (log10 X_eq) and the direct M&P comparison (yield %).
    _heatmap(combined, out_dir, "log10_X_eq", "nh3_combined_heatmap_log10X",
             "Amino-acid accessibility vs NH3 — log10 equilibrium mole fraction", "log10 X_eq")
    _heatmap(combined, out_dir, "yield_pct_HCN", "nh3_combined_heatmap_yield_pct_HCN",
             "Amino-acid yield vs NH3 — % relative to initial HCN", "yield (% of HCN)")
    _line_plots(combined, out_dir, "log10_X_eq", "nh3_combined_lines_log10X",
                "log10 X_eq", threshold=SIGNIFICANT_LOG10X)
    _line_plots(combined, out_dir, "yield_pct_HCN", "nh3_combined_lines_yield_pct",
                "yield (% HCN)")
    _write_summary(combined, out_dir / "summary.md", args.ratio)

    print(f"Wrote {out_dir} ({combined['amino_acid'].nunique()} amino acids, "
          f"{len(combined)} points).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
