"""Scaffold a new sensitivity study folder.

Copies studies/_template/study_config.yaml to studies/<id>/study_config.yaml
(substituting the study id) and creates the study subdirectories, so a user never
has to wonder where files go.

Usage:
    python scripts/new_study.py --id my_study
    python scripts/new_study.py --id my_study --force   # overwrite existing
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = PROJECT_ROOT / "studies" / "_template" / "study_config.yaml"
STUDY_SUBDIRS = ("processed", "models", "results", "figures")
PLACEHOLDER = "__STUDY_ID__"


def scaffold_study(study_id: str, force: bool = False) -> Path:
    """Create studies/<study_id>/ with a runnable study_config.yaml."""
    if not TEMPLATE.exists():
        raise FileNotFoundError(f"Template not found: {TEMPLATE}")

    study_dir = PROJECT_ROOT / "studies" / study_id
    config_path = study_dir / "study_config.yaml"
    if config_path.exists() and not force:
        raise FileExistsError(
            f"{config_path} already exists. Use --force to overwrite."
        )

    study_dir.mkdir(parents=True, exist_ok=True)
    for sub in STUDY_SUBDIRS:
        (study_dir / sub).mkdir(parents=True, exist_ok=True)

    text = TEMPLATE.read_text(encoding="utf-8").replace(PLACEHOLDER, study_id)
    config_path.write_text(text, encoding="utf-8")
    return config_path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--id", required=True, dest="study_id",
                        help="Study id, e.g. my_study (becomes studies/my_study/).")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite an existing study_config.yaml.")
    args = parser.parse_args(argv)

    try:
        config_path = scaffold_study(args.study_id, force=args.force)
    except (FileExistsError, FileNotFoundError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Scaffolded study '{args.study_id}':")
    print(f"  config: {config_path}")
    print("Next: edit that file, then preview with")
    print(f"  python scripts/generate_sensitivity_design.py "
          f"--config {config_path.relative_to(PROJECT_ROOT)} --dry-run")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
