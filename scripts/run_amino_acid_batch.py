"""Run the amino-acid breadth campaign: each scaffolded study through the engine.

Discovers ``studies/amino_acid_scan/<key>/study_config.yaml`` and, for each, calls the
SAME entrypoints the single-study scripts use — ``run_sensitivity_study.main`` →
``summarize_sensitivity_study.main`` → ``plot_sensitivity_study.main`` — in sequence.
Resume-by-default per study; one failing study never aborts the batch.

Usage:
    python scripts/run_amino_acid_batch.py --dry-run
    python scripts/run_amino_acid_batch.py
    python scripts/run_amino_acid_batch.py --only glycine,serine --steps run,summarize
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from sensitivity_design import (  # noqa: E402
    build_full_design_matrix, build_thermo_offsets_table, load_study_config,
)

DEFAULT_OUT = "studies/amino_acid_scan"
ALL_STEPS = ("run", "summarize", "plot")


def discover_studies(out_dir: str | Path, only: Optional[List[str]] = None) -> List[Path]:
    """Return the study_config.yaml paths under ``out_dir`` (optionally a subset)."""
    paths = sorted((PROJECT_ROOT / out_dir).glob("*/study_config.yaml"))
    if only:
        wanted = set(only)
        paths = [p for p in paths if p.parent.name in wanted]
    return paths


def _per_case_estimate() -> float:
    """Median seconds/case from the committed alanine run, or a nominal fallback."""
    rm = PROJECT_ROOT / "studies" / "alanine_mvp" / "run_manifest.csv"
    if rm.exists():
        rt = pd.read_csv(rm)["runtime_seconds"].dropna()
        if len(rt):
            return float(rt.median())
    return 0.0015


def dry_run_summary(config_paths: List[Path]) -> List[dict]:
    """Per-study case count + unique-model count (no execution)."""
    rows = []
    for cfg_path in config_paths:
        cfg = load_study_config(cfg_path)
        design = build_full_design_matrix(cfg)
        offsets = build_thermo_offsets_table(cfg)
        rows.append({
            "key": cfg_path.parent.name,
            "n_cases": len(design),
            "n_models": 1 + len(offsets),   # 1 base model + one per ΔG variant
        })
    return rows


BATCH_SUMMARY_COLUMNS = [
    "study_id", "target_product", "n_cases", "n_ok", "n_failed", "n_suspect_balance",
    "total_runtime_seconds", "inventory_accessible_fraction", "max_X_eq",
    "summary_path", "figures_path",
]


def collect_batch_summary(config_paths: List[Path]) -> pd.DataFrame:
    """One row per study from its written outputs (NaN where an output is missing)."""
    import numpy as np
    rows = []
    for cfg_path in config_paths:
        cfg = load_study_config(cfg_path)
        study_dir = cfg_path.parent
        target = (cfg.get("mode", {}).get("target_products") or [None])[0]
        results = study_dir / "results"
        figures = study_dir / "figures"
        row = {c: np.nan for c in BATCH_SUMMARY_COLUMNS}
        row.update({"study_id": cfg["study"]["study_id"], "target_product": target,
                    "summary_path": str(results / "sensitivity_case_summary.csv"),
                    "figures_path": str(figures)})
        rm_path = study_dir / "run_manifest.csv"
        if rm_path.exists():
            rm = pd.read_csv(rm_path)
            row["n_cases"] = len(rm)
            row["n_ok"] = int((rm["status"] == "ok").sum())
            row["n_failed"] = int((rm["status"] == "failed").sum())
            row["total_runtime_seconds"] = float(rm["runtime_seconds"].dropna().sum())
        cs_path = results / "sensitivity_case_summary.csv"
        if cs_path.exists():
            cs = pd.read_csv(cs_path)
            row["n_suspect_balance"] = int(cs.get("suspect_balance", pd.Series(dtype=bool)).sum())
            inv = cs[cs["substudy_id"] == "inventory_landscape"]
            if len(inv):
                row["inventory_accessible_fraction"] = float((inv["formed_bool"] == True).mean())  # noqa: E712
                row["max_X_eq"] = float(inv["X_eq"].max(skipna=True))
        rows.append(row)
    return pd.DataFrame(rows, columns=BATCH_SUMMARY_COLUMNS)


def run_batch(config_paths: List[Path], steps=ALL_STEPS, force: bool = False,
              limit: Optional[int] = None, substudy: Optional[str] = None,
              run_fn=None, summarize_fn=None, plot_fn=None, progress: bool = True) -> List[dict]:
    """Run each study through the selected steps; never abort on one failure."""
    if run_fn is None:
        from run_sensitivity_study import main as run_fn
    if summarize_fn is None:
        from summarize_sensitivity_study import main as summarize_fn
    if plot_fn is None:
        from plot_sensitivity_study import main as plot_fn

    results = []
    for cfg_path in config_paths:
        key = cfg_path.parent.name
        status, failed_step, message = "ok", None, ""
        if progress:
            print(f"\n=== {key} ===")
        for step in steps:
            try:
                if step == "run":
                    argv = ["--config", str(cfg_path)]
                    if force:
                        argv.append("--force")
                    if limit is not None:
                        argv += ["--limit", str(limit)]
                    if substudy:
                        argv += ["--substudy", substudy]
                    rc = run_fn(argv)
                elif step == "summarize":
                    rc = summarize_fn(["--config", str(cfg_path)])
                elif step == "plot":
                    argv = ["--config", str(cfg_path)]
                    if substudy:
                        argv += ["--substudy", substudy]
                    rc = plot_fn(argv)
                else:
                    continue
                if rc != 0:
                    status, failed_step, message = "failed", step, f"exit code {rc}"
                    break
            except Exception as exc:  # one study's failure must not abort the batch
                status, failed_step, message = "failed", step, repr(exc)
                break
        results.append({"key": key, "status": status,
                        "failed_step": failed_step, "message": message})
    return results


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-dir", default=DEFAULT_OUT,
                        help="Campaign directory of scaffolded studies (e.g. studies/aa_nh3).")
    parser.add_argument("--only", default=None, help="Comma-separated amino-acid keys.")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--substudy", default=None)
    parser.add_argument("--steps", default=",".join(ALL_STEPS),
                        help="Comma-separated subset of run,summarize,plot.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print per-AA case/model counts + projected runtime; run nothing.")
    args = parser.parse_args(argv)

    steps = tuple(s.strip() for s in args.steps.split(",") if s.strip())
    invalid_steps = [s for s in steps if s not in ALL_STEPS]
    if not steps or invalid_steps:
        print(f"Invalid --steps {invalid_steps or '(empty)'}; allowed: {', '.join(ALL_STEPS)}.",
              file=sys.stderr)
        return 1

    only = [k.strip() for k in args.only.split(",")] if args.only else None
    all_paths = discover_studies(args.scan_dir)
    if not all_paths:
        print(f"No studies found under {args.scan_dir}. Run scripts/new_amino_acid_batch.py "
              "first.", file=sys.stderr)
        return 1
    if only:
        valid = {p.parent.name for p in all_paths}
        unknown = [k for k in only if k not in valid]
        if unknown:
            print(f"Unknown study key(s) {unknown}; available: {sorted(valid)}.",
                  file=sys.stderr)
            return 1
    config_paths = discover_studies(args.scan_dir, only)

    if args.dry_run:
        summary = dry_run_summary(config_paths)
        per_case = _per_case_estimate()
        total_cases = sum(r["n_cases"] for r in summary)
        print(f"{'amino acid':<18}{'cases':>8}{'models':>8}")
        for r in summary:
            print(f"{r['key']:<18}{r['n_cases']:>8}{r['n_models']:>8}")
        print(f"{'TOTAL':<18}{total_cases:>8}")
        print(f"\nProjected runtime ≈ {total_cases * per_case:.0f} s "
              f"({len(summary)} studies × ~{summary[0]['n_cases'] if summary else 0} cases, "
              f"~{per_case:.4f} s/case).")
        print("[dry-run] Nothing executed.")
        return 0

    results = run_batch(config_paths, steps=steps, force=args.force,
                        limit=args.limit, substudy=args.substudy)
    n_ok = sum(r["status"] == "ok" for r in results)

    # Campaign-level batch_summary.csv (one row per amino acid) in the scan dir.
    summary_df = collect_batch_summary(config_paths)
    summary_path = PROJECT_ROOT / args.scan_dir / "batch_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    print("\n=== Batch summary ===")
    print(f"{'amino acid':<18}{'status':<10}{'failed step':<12}")
    for r in results:
        print(f"{r['key']:<18}{r['status']:<10}{r['failed_step'] or '':<12}{r['message']}")
    print(f"\n{n_ok}/{len(results)} studies ok. Wrote {summary_path}")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
