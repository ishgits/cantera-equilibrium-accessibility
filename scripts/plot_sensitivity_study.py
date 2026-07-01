"""Plot a sensitivity study's landscape grid (config-driven styling).

Phase 3. Reads sensitivity_landscape_grid.csv and the `plots:` block of the study
config, and writes figures into studies/<id>/figures/ in each requested format.
Failed cases stay visible (gray cells). matplotlib only.

Usage:
    python scripts/plot_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml
    python scripts/plot_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml --substudy inventory_landscape
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import pandas as pd  # noqa: E402

from sensitivity_design import StudyConfigError, load_study_config  # noqa: E402
from sensitivity_plotting import plot_all  # noqa: E402


def _study_dir(config) -> Path:
    return PROJECT_ROOT / config["study"].get(
        "output_dir", f"studies/{config['study']['study_id']}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--substudy", default=None,
                        help="Plot only one substudy, e.g. inventory_landscape.")
    args = parser.parse_args(argv)

    try:
        config = load_study_config(args.config)
    except StudyConfigError as exc:
        print(f"Config error: {exc}", file=sys.stderr)
        return 1

    if not config.get("outputs", {}).get("make_plots", True):
        print("outputs.make_plots is false — skipping figure generation.")
        return 0

    study_dir = _study_dir(config)
    grid_path = study_dir / "results" / "sensitivity_landscape_grid.csv"
    if not grid_path.exists():
        print(f"No landscape grid at {grid_path}. Run "
              "scripts/summarize_sensitivity_study.py first.", file=sys.stderr)
        return 1

    grid = pd.read_csv(grid_path)
    plots_cfg = config.get("plots", {}) or {}
    written = plot_all(grid, plots_cfg, study_dir / "figures", substudy=args.substudy)

    for name, paths in written.items():
        if paths:
            print(f"{name}: {', '.join(str(p) for p in paths)}")
        else:
            print(f"{name}: no data yet (skipped)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
