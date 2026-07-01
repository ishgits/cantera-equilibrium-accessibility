"""Scaffold one single-target sensitivity study per amino acid (breadth campaign).

Generates ``studies/amino_acid_scan/<key>/study_config.yaml`` for each amino acid in
``inputs/amino_acids_species.csv`` (rows with ``product_class == amino_acid``), all
cloned from the alanine template with the **same** feedstock, sweeps, and thresholds
so the resulting landscapes are directly comparable. Only the target, the study
identity, the shared input files, and the figure labels differ.

Targets are derived from the species CSV — never hardcoded.

Usage:
    python scripts/new_amino_acid_batch.py
    python scripts/new_amino_acid_batch.py --only glycine,serine --force
"""
from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path
from typing import List, Optional

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

DEFAULT_SPECIES_CSV = "inputs/amino_acids_species.csv"
DEFAULT_GIBBS_SEED = "inputs/amino_acids_gibbs_seed.csv"
DEFAULT_TEMPLATE = "studies/alanine_mvp"
DEFAULT_OUT = "studies/amino_acid_scan"


def amino_acid_targets(species_csv: str | Path) -> List[dict]:
    """Return ``[{key, cantera_name, display}]`` for each amino acid in the CSV."""
    df = pd.read_csv(species_csv)
    aas = df[df["product_class"] == "amino_acid"]
    targets = []
    for _, r in aas.iterrows():
        targets.append({
            "key": str(r["species_key"]),
            "cantera_name": str(r["cantera_name"]),
            "display": str(r["chnosz_name"]),  # natural name, e.g. "aspartic acid"
        })
    return targets


def _substitute(obj, low: str, title: str):
    """Recursively replace the alanine name with the target name in string values."""
    if isinstance(obj, dict):
        return {k: _substitute(v, low, title) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_substitute(v, low, title) for v in obj]
    if isinstance(obj, str):
        return obj.replace("Alanine", title).replace("alanine", low)
    return obj


def _apply_nh3_min(cfg: dict, nh3_min: float) -> None:
    """Force every NH3(aq) sweep axis to start at nh3_min (> 0); never generate 0."""
    for sweep in (cfg.get("sweeps") or {}).values():
        if not isinstance(sweep, dict):
            continue
        spec = sweep.get("variables", {}).get("NH3(aq)")
        if not isinstance(spec, dict):
            continue
        if spec.get("type") == "explicit":
            spec["values"] = [v for v in spec.get("values", []) if float(v) >= nh3_min]
        else:
            spec["min"] = nh3_min


def _exclude_species(cfg: dict, exclude: set) -> None:
    """Remove species from the phase, inventory, and sweep axes; disable sweeps that
    lose a required species axis (e.g. nh3_deltaG_landscape when NH3 is excluded)."""
    if not exclude:
        return
    cfg["model"]["allowed_species"] = [s for s in cfg["model"]["allowed_species"]
                                       if s not in exclude]
    cfg["fixed_inventory"] = {k: v for k, v in (cfg.get("fixed_inventory") or {}).items()
                              if k not in exclude}
    for name, sweep in (cfg.get("sweeps") or {}).items():
        if not isinstance(sweep, dict):
            continue
        if "variables" in sweep:
            sweep["variables"] = {k: v for k, v in sweep["variables"].items()
                                  if k not in exclude}
        if "fixed_inventory" in sweep:
            sweep["fixed_inventory"] = {k: v for k, v in sweep["fixed_inventory"].items()
                                        if k not in exclude}
        if name == "nh3_deltaG_landscape" and "NH3(aq)" in exclude:
            sweep["enabled"] = False   # this substudy inherently needs ammonia


def build_study_config(template_cfg: dict, target: dict, out_dir_rel: str,
                       species_csv: str, gibbs_seed: str,
                       exclude_species: tuple = (), nh3_min: float = 0.01) -> dict:
    """Clone the template config and apply the per-amino-acid substitutions."""
    if nh3_min <= 0:
        raise ValueError("--nh3-min must be > 0 (NH3=0 is the excluded batch only).")
    cfg = copy.deepcopy(template_cfg)
    key, cantera, display = target["key"], target["cantera_name"], target["display"]
    title = display.title()

    cfg["study"]["study_id"] = key
    cfg["study"]["output_dir"] = f"{out_dir_rel}/{key}"
    cfg["study"]["description"] = (
        f"{title} sensitivity — Titan amino-acid accessibility landscape (breadth campaign)")
    cfg["mode"]["target_products"] = [cantera]
    cfg["species_files"]["species_csv"] = species_csv
    cfg["species_files"]["gibbs_seed_wide_csv"] = gibbs_seed
    if "deltaG_sweep" in cfg.get("sweeps", {}):
        cfg["sweeps"]["deltaG_sweep"]["species"] = cantera

    _apply_nh3_min(cfg, nh3_min)
    _exclude_species(cfg, set(exclude_species))

    if "plots" in cfg:
        cfg["plots"] = _substitute(cfg["plots"], display, title)
    return cfg


def scaffold_studies(species_csv: str | Path = DEFAULT_SPECIES_CSV,
                     template_dir: str | Path = DEFAULT_TEMPLATE,
                     out_dir: str | Path = DEFAULT_OUT,
                     only: Optional[List[str]] = None,
                     force: bool = False,
                     gibbs_seed: str = DEFAULT_GIBBS_SEED,
                     exclude_species: tuple = (),
                     nh3_min: float = 0.01) -> dict:
    """Generate the per-amino-acid study configs. Returns {created, skipped} paths."""
    species_csv = Path(species_csv)
    template_cfg = yaml.safe_load(
        (PROJECT_ROOT / template_dir / "study_config.yaml").read_text(encoding="utf-8"))

    all_targets = amino_acid_targets(species_csv)
    if only:
        valid = {t["key"] for t in all_targets}
        unknown = [k for k in only if k not in valid]
        if unknown:
            raise ValueError(
                f"Unknown amino-acid key(s): {unknown}. Valid keys: {sorted(valid)}")
        targets = [t for t in all_targets if t["key"] in only]
    else:
        targets = all_targets

    header = (f"# Auto-generated by scripts/new_amino_acid_batch.py — "
              f"edit the template ({template_dir}) or regenerate.\n")
    out_root = Path(out_dir)
    out_dir_rel = str(out_dir)        # written into each config as the output_dir base
    species_ref = str(species_csv)    # path the generated configs point at
    created, skipped = [], []
    for t in targets:
        study_dir = out_root / t["key"]
        config_path = study_dir / "study_config.yaml"
        if config_path.exists() and not force:
            skipped.append(config_path)
            continue
        cfg = build_study_config(template_cfg, t, out_dir_rel, species_ref, gibbs_seed,
                                 exclude_species=exclude_species, nh3_min=nh3_min)
        study_dir.mkdir(parents=True, exist_ok=True)
        config_path.write_text(
            header + yaml.safe_dump(cfg, sort_keys=False, default_flow_style=False),
            encoding="utf-8")
        created.append(config_path)
    return {"created": created, "skipped": skipped, "targets": targets}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", default=DEFAULT_SPECIES_CSV,
                        help="Species CSV (amino acids derived from product_class).")
    parser.add_argument("--template", default=DEFAULT_TEMPLATE,
                        help="Study folder whose study_config.yaml is the template.")
    parser.add_argument("--out", default=DEFAULT_OUT, help="Output campaign directory.")
    parser.add_argument("--only", default=None,
                        help="Comma-separated subset of amino-acid keys.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing configs.")
    parser.add_argument("--exclude-species", default=None,
                        help="Comma-separated species to exclude from the phase "
                             "(e.g. 'NH3(aq)' for the no-ammonia batch).")
    parser.add_argument("--nh3-min", type=float, default=0.01,
                        help="Minimum NH3(aq) sweep value (>0; default 0.01 = 1%% of water).")
    args = parser.parse_args(argv)

    only = [k.strip() for k in args.only.split(",")] if args.only else None
    exclude = tuple(s.strip() for s in args.exclude_species.split(",")) if args.exclude_species else ()
    try:
        result = scaffold_studies(args.species, args.template, args.out, only, args.force,
                                  exclude_species=exclude, nh3_min=args.nh3_min)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Amino-acid targets: {len(result['targets'])}")
    for p in result["created"]:
        print(f"  created  {p}")
    for p in result["skipped"]:
        print(f"  skipped  {p} (exists; use --force)")
    print(f"Done: {len(result['created'])} created, {len(result['skipped'])} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
