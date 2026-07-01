"""Generate a sensitivity study's design matrix and scenario YAML.

Phase 1 of the sensitivity layer. Reads a ``study_config.yaml``, validates it
(plain-English errors), and writes:

    studies/<id>/design_matrix.csv
    studies/<id>/generated_scenarios.yaml
    studies/<id>/thermo_offsets.csv   (only if a ΔG substudy is enabled)

No Cantera is run. Use ``--dry-run`` to preview case counts, grid ranges, and ΔG
variant names *before* writing anything.

Usage:
    python scripts/generate_sensitivity_design.py --config studies/alanine_mvp/study_config.yaml
    python scripts/generate_sensitivity_design.py --config studies/alanine_mvp/study_config.yaml --dry-run
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from sensitivity_design import (  # noqa: E402
    StudyConfigError,
    _SUBSTUDY_BUILDERS,
    _enabled_sweeps,
    build_full_design_matrix,
    build_thermo_offsets_table,
    load_gibbs_seed_for_config,
    load_species_for_config,
    load_study_config,
    validate_study,
    write_design_outputs,
)


def _print_preview(config: dict) -> None:
    enabled = _enabled_sweeps(config)
    print(f"Study:   {config['study']['study_id']}")
    print(f"Target:  {', '.join(config.get('mode', {}).get('target_products') or [])}")
    print(f"Output:  {config['study'].get('output_dir')}")
    print("\nEnabled substudies (run counts):")
    total = 0
    for substudy_id in ("inventory_landscape", "deltaG_sweep", "nh3_deltaG_landscape"):
        if substudy_id not in enabled:
            continue
        n = len(_SUBSTUDY_BUILDERS[substudy_id](config))
        total += n
        print(f"  - {substudy_id:<22} {n:>6} cases")
    print(f"  {'TOTAL':<22} {total:>6} cases")

    matrix = build_full_design_matrix(config)
    print("\nGrid ranges (design variables):")
    for col in ("NH3_mol", "C2H2_mol", "C2H2_over_HCN", "deltaG_offset_kJ_mol"):
        if col in matrix:
            print(f"  - {col:<22} {matrix[col].min():>10.4g} .. {matrix[col].max():<10.4g}")

    offsets = build_thermo_offsets_table(config)
    if not offsets.empty:
        print("\nΔG pseudo-species variants:")
        for _, row in offsets.iterrows():
            print(f"  - {row['variant_species']:<26} "
                  f"({row['deltaG_offset_kJ_mol']:+.0f} kJ/mol)")
    print("\n[dry-run] No files written.")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Path to study_config.yaml")
    parser.add_argument("--dry-run", action="store_true",
                        help="Preview counts/ranges/variants without writing files.")
    args = parser.parse_args(argv)

    try:
        config = load_study_config(args.config)
        species_df = load_species_for_config(config, PROJECT_ROOT)
        gibbs_seed_df = load_gibbs_seed_for_config(config, PROJECT_ROOT)
        validate_study(config, species_df, gibbs_seed_df)
    except StudyConfigError as exc:
        print(f"Config error: {exc}", file=sys.stderr)
        return 1

    if args.dry_run:
        _print_preview(config)
        return 0

    design_df = build_full_design_matrix(config)
    paths = write_design_outputs(design_df, config)
    print(f"Wrote {len(design_df)} cases.")
    for name, path in paths.items():
        print(f"  {name:<22} {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
