"""Summarize a sensitivity study's raw results into case summaries and metrics.

Phase 3. Reads the merged raw-long output produced by the runner and writes:

    studies/<id>/results/equilibrium_moles_long.csv     (ok cases, reconstructed moles)
    studies/<id>/results/sensitivity_case_summary.csv   (one row per case)
    studies/<id>/results/sensitivity_landscape_grid.csv (plot/ML-ready subset)
    studies/<id>/results/sensitivity_run_summary.md     (human-readable metrics)

No Cantera needed — this operates on the raw CSV.

Usage:
    python scripts/summarize_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import pandas as pd  # noqa: E402

from mole_balance import add_equilibrium_moles  # noqa: E402
from sensitivity_design import (  # noqa: E402
    StudyConfigError, load_species_for_config, load_study_config,
)
from sensitivity_thermo import augment_species_metadata_with_variants  # noqa: E402
from sensitivity_design import build_thermo_offsets_table  # noqa: E402
from sensitivity_summary import (  # noqa: E402
    compute_sensitivity_metrics, make_landscape_grid, summarize_sensitivity_cases,
    write_schema_dictionary, write_sensitivity_run_summary,
)


def _study_dir(config) -> Path:
    return PROJECT_ROOT / config["study"].get(
        "output_dir", f"studies/{config['study']['study_id']}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)

    try:
        config = load_study_config(args.config)
        species_df = load_species_for_config(config, PROJECT_ROOT)
    except StudyConfigError as exc:
        print(f"Config error: {exc}", file=sys.stderr)
        return 1

    study_dir = _study_dir(config)
    raw_path = study_dir / "results" / "equilibrium_raw_long.csv"
    if not raw_path.exists():
        print(f"No raw results at {raw_path}. Run scripts/run_sensitivity_study.py "
              "first (in a Cantera-enabled environment).", file=sys.stderr)
        return 1

    raw_long = pd.read_csv(raw_path)
    # Variant species may appear in results; include their metadata for reconstruction.
    species_meta = augment_species_metadata_with_variants(
        species_df, build_thermo_offsets_table(config))
    thresholds = config.get("thresholds", {})

    # Moles long (ok cases only — failed cases have no reconstructable moles).
    ok_raw = raw_long[raw_long["solver_status"].astype(str) == "ok"]
    results_dir = study_dir / "results"
    if not ok_raw.empty:
        moles_long = add_equilibrium_moles(ok_raw, species_meta, group_cols=["case_id"])
        moles_long.to_csv(results_dir / "equilibrium_moles_long.csv", index=False)

    case_summary = summarize_sensitivity_cases(
        raw_long, species_meta, thresholds,
        output_csv=results_dir / "sensitivity_case_summary.csv")
    make_landscape_grid(case_summary,
                        output_csv=results_dir / "sensitivity_landscape_grid.csv")
    metrics = compute_sensitivity_metrics(
        case_summary,
        significant_X_threshold=float(thresholds.get("significant_X_threshold", 1e-6)))
    write_sensitivity_run_summary(
        metrics, results_dir / "sensitivity_run_summary.md",
        study_id=config["study"]["study_id"])

    # Column dictionary, shipped next to the data (review §8, ML-readiness).
    if config.get("outputs", {}).get("write_schema_dictionary"):
        md_path, json_path = write_schema_dictionary(results_dir)
        print(f"Schema dictionary: {md_path}, {json_path}")

    g = metrics["general"]
    print(f"Summarized {g['total_cases']} cases "
          f"({g['failed_cases']} failed, {g['suspect_balance_cases']} suspect balance).")
    print(f"Outputs in {results_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
