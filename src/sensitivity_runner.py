"""Sensitivity study runner (Phase 2).

Turns a Phase-1 design matrix into model/run manifests and executes the cases,
reusing the validated v1.0 engine rather than reimplementing it:

- **Model identity is hashed** over ``(sorted(allowed_species), thermo_variant_id,
  target_variant, phase_name, nasa9_coeff_hash)`` (review §3.1). A constant
  ``allowed_species`` across an inventory sweep therefore collapses to **one** model
  YAML reused by every grid point (625 inventory cases → 1 model).
- **One runner**: :func:`run_sensitivity_manifest` wraps
  ``equilibrium_runner.run_single_yaml_case`` (review §3.4). It adds only the
  orchestration the base lacks — per-case timing, status persistence, design-variable
  merge, ``--only-failed`` / ``--limit``, and never aborting the whole study on one
  failed case. It does **not** duplicate the solver call.
- **YAML generation reuses** ``cantera_yaml._species_block`` / ``_yaml_document``
  (the validated builders), but path resolution is **study-local** (review §3.5);
  the base ``_portable_path`` / ``resolve_yaml_path_from_manifest_row`` helpers are
  deliberately not reused (review §2.4).

Manifest building, YAML generation and provenance need no Cantera; only the actual
``equilibrate`` call does.
"""
from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from cantera_yaml import _species_block, _yaml_document, write_yaml_document
from thermo_fit import COEFF_COLUMNS

MODEL_MANIFEST_COLUMNS = [
    "model_id", "yaml_path", "model_mode", "target_variant", "base_species",
    "thermo_variant_id", "species_set_hash", "thermo_hash",
]
RUN_MANIFEST_COLUMNS = [
    "case_id", "scenario_id", "model_id", "yaml_path", "target_product",
    "target_variant", "T_C", "P_Pa", "case_hash", "status", "runtime_seconds",
    "error_message",
]
# Case-defining fields hashed into case_hash. allowed_species is captured transitively
# via model_id (the model identity hashes sorted allowed_species + thermo); initial
# moles are captured via the inventory *_mol columns. A change to any invalidates a
# prior result, so a matching case_id with a differing case_hash is re-run.
CASE_HASH_FIELDS = [
    "case_id", "target_product", "target_variant", "T_C", "P_Pa",
    "H2O_mol", "HCN_mol", "C2H2_mol", "NH3_mol", "C2H2_over_HCN",
    "deltaG_offset_kJ_mol", "model_id",
]
# Design columns merged onto every raw-long row (review §14.2: preserve design vars).
# NOTE: model_id, target_product and target_variant are intentionally absent — the
# runner stamps those canonical values directly on each raw row, so re-merging them
# here would create *_design duplicates.
_DESIGN_MERGE_COLUMNS = [
    "case_id", "study_id", "substudy_id",
    "H2O_mol", "HCN_mol", "C2H2_mol", "NH3_mol", "C2H2_over_HCN",
    "deltaG_offset_kJ_mol",
]


# --------------------------------------------------------------------------- #
# Model identity
# --------------------------------------------------------------------------- #
def thermo_variant_id_for(target_variant: str) -> str:
    """``Alanine(aq) -> 'base'``; ``Alanine__dG_m040(aq) -> 'dG_m040'``."""
    name = str(target_variant)
    if "__" in name:
        return name.split("__", 1)[1].split("(")[0]
    return "base"


def nasa9_coeff_hash(coeff_df: pd.DataFrame, species: Sequence[str]) -> str:
    """Stable hash of the NASA9 coefficient rows for ``species``.

    Different ``a7`` values (e.g. a ΔG-shifted variant) yield a different hash,
    so the model identity tracks the thermodynamics, not just the species names.
    """
    cols = ["cantera_name", "range_label", "T_low_K", "T_high_K", *COEFF_COLUMNS]
    sub = coeff_df[coeff_df["cantera_name"].isin(list(species))][cols]
    sub = sub.sort_values(["cantera_name", "range_label"]).reset_index(drop=True)
    # Canonical, precision-stable serialization.
    payload = sub.to_csv(index=False, float_format="%.12g")
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


def model_identity(allowed_species: Sequence[str], thermo_variant_id: str,
                   target_variant: str, phase_name: str, coeff_hash: str) -> str:
    """Deterministic ``model_id`` from the full identity tuple."""
    key = "|".join([
        ";".join(sorted(str(s) for s in allowed_species)),
        str(thermo_variant_id),
        str(target_variant),
        str(phase_name),
        str(coeff_hash),
    ])
    return "M_" + hashlib.sha256(key.encode("utf-8")).hexdigest()[:12]


def _allowed_for_row(target_variant: str, allowed_base: Sequence[str]) -> List[str]:
    """Reactant phase set + this case's single target variant (ordered, de-duped)."""
    return list(dict.fromkeys([str(s) for s in allowed_base] + [str(target_variant)]))


def _hash_field(value) -> str:
    """Stable string for a hash field (fixed-precision floats; str otherwise)."""
    try:
        return f"{float(value):.10g}"
    except (TypeError, ValueError):
        return str(value)


def compute_case_hash(row) -> str:
    """Deterministic hash over the case-defining fields (:data:`CASE_HASH_FIELDS`)."""
    key = "|".join(f"{f}={_hash_field(row.get(f, ''))}" for f in CASE_HASH_FIELDS)
    return hashlib.sha256(key.encode("utf-8")).hexdigest()[:16]


def add_case_hash(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of ``df`` with a ``case_hash`` column."""
    df = df.copy()
    df["case_hash"] = [compute_case_hash(r) for _, r in df.iterrows()]
    return df


# --------------------------------------------------------------------------- #
# Manifest building (pure — no Cantera, no writes)
# --------------------------------------------------------------------------- #
def compute_model_table(design_df: pd.DataFrame, allowed_base: Sequence[str],
                        coeff_df: pd.DataFrame, phase_name: str = "aqueous",
                        model_mode: str = "single_product_sensitivity"):
    """Assign a ``model_id`` to every design row and derive the unique-model table.

    Returns ``(design_with_model_id, models_df)``. ``models_df`` has one row per
    unique model (:data:`MODEL_MANIFEST_COLUMNS`) with a planned, study-relative
    ``yaml_path`` of ``models/<model_id>.yaml``. Nothing is written.
    """
    design = design_df.copy()
    model_ids: List[str] = []
    records: Dict[str, Dict[str, Any]] = {}

    for _, row in design.iterrows():
        target_variant = str(row["target_variant"])
        allowed = _allowed_for_row(target_variant, allowed_base)
        tvid = thermo_variant_id_for(target_variant)
        coeff_hash = nasa9_coeff_hash(coeff_df, allowed)
        mid = model_identity(allowed, tvid, target_variant, phase_name, coeff_hash)
        model_ids.append(mid)
        if mid not in records:
            base_species = [s for s in allowed if s != target_variant]
            records[mid] = {
                "model_id": mid,
                "yaml_path": f"models/{mid}.yaml",
                "model_mode": model_mode,
                "target_variant": target_variant,
                "base_species": ";".join(base_species),
                "thermo_variant_id": tvid,
                "species_set_hash": hashlib.sha256(
                    ";".join(sorted(allowed)).encode("utf-8")).hexdigest()[:12],
                "thermo_hash": coeff_hash,
                "_allowed_species": allowed,  # internal: ordered phase list for YAML
            }

    design["model_id"] = model_ids
    models_df = pd.DataFrame(list(records.values()))
    return design, models_df


def model_reuse_stats(design_with_model_id: pd.DataFrame, models_df: pd.DataFrame,
                      expected_models: Optional[int] = None) -> Dict[str, Any]:
    """Report (and optionally assert) the model-reuse ratio (review §5.4)."""
    n_cases = len(design_with_model_id)
    n_models = len(models_df)
    stats = {
        "n_cases": n_cases,
        "n_models": n_models,
        "reuse_ratio": (n_cases / n_models) if n_models else 0.0,
    }
    if expected_models is not None and n_models != expected_models:
        raise AssertionError(
            f"Model-reuse regression: expected {expected_models} unique model(s) "
            f"but built {n_models}. Caching may have broken (per-case YAMLs)."
        )
    return stats


# --------------------------------------------------------------------------- #
# YAML generation (reuses validated cantera_yaml builders)
# --------------------------------------------------------------------------- #
def generate_model_yamls(models_df: pd.DataFrame, species_meta: pd.DataFrame,
                         coeff_df: pd.DataFrame, study_dir: str | Path,
                         phase_name: str = "aqueous") -> pd.DataFrame:
    """Write one Cantera YAML per unique model under ``<study_dir>/models/``.

    Raises a clear error if a species in a model has no metadata or no NASA9
    coefficients (e.g. a ΔG variant whose shifted coefficients were not generated).
    """
    study_dir = Path(study_dir)
    meta = species_meta.set_index("cantera_name", drop=False)
    coeff_species = set(coeff_df["cantera_name"])

    for _, model in models_df.iterrows():
        allowed = list(model["_allowed_species"])
        missing_meta = [s for s in allowed if s not in meta.index]
        if missing_meta:
            raise KeyError(
                f"Model {model['model_id']} needs species metadata for "
                f"{missing_meta} — add rows to the species CSV."
            )
        missing_coeff = [s for s in allowed if s not in coeff_species]
        if missing_coeff:
            raise KeyError(
                f"Model {model['model_id']} has no NASA9 coefficients for "
                f"{missing_coeff}. ΔG variant coefficients are generated by "
                "prepare_study_coefficients (analytic a7 shift) before model building "
                "— run the study to materialize them."
            )
        blocks = [_species_block(meta.loc[s], coeff_df) for s in allowed]
        doc = _yaml_document(allowed, blocks, phase_name=phase_name)
        write_yaml_document(doc, study_dir / model["yaml_path"])
    return models_df


def build_model_manifest(design_df: pd.DataFrame, species_meta: pd.DataFrame,
                         coeff_df: pd.DataFrame, allowed_base: Sequence[str],
                         study_dir: str | Path, phase_name: str = "aqueous",
                         model_mode: str = "single_product_sensitivity",
                         expected_models: Optional[int] = None):
    """Compute identities, generate YAMLs, and write ``model_manifest.csv``.

    Returns ``(design_with_model_id, models_df, reuse_stats)``.
    """
    study_dir = Path(study_dir)
    design, models_df = compute_model_table(
        design_df, allowed_base, coeff_df, phase_name, model_mode)
    stats = model_reuse_stats(design, models_df, expected_models)
    generate_model_yamls(models_df, species_meta, coeff_df, study_dir, phase_name)

    manifest = models_df[MODEL_MANIFEST_COLUMNS].copy()
    manifest.to_csv(study_dir / "model_manifest.csv", index=False)
    return design, models_df, stats


def build_run_manifest(design_with_model_id: pd.DataFrame, models_df: pd.DataFrame,
                       study_dir: str | Path) -> pd.DataFrame:
    """Join each design case to its model and write a **resume-safe** run manifest.

    The design is the canonical case set: one row per case. If ``run_manifest.csv``
    already exists, prior ``status``/``runtime_seconds``/``error_message`` are merged
    back by ``case_id`` so completed work survives a rerun (``--only-failed`` keeps
    working). Cases dropped from the design disappear; new cases start ``pending``.
    Statuses are reset to ``pending`` only when a case is explicitly selected for
    rerun (``--force``, handled in :func:`run_sensitivity_manifest`).
    """
    study_dir = Path(study_dir)
    yaml_by_model = dict(zip(models_df["model_id"], models_df["yaml_path"]))
    rows = []
    for _, case in design_with_model_id.iterrows():
        model_id = str(case["model_id"])
        rows.append({
            "case_id": str(case["case_id"]),
            "scenario_id": str(case["scenario_id"]),
            "model_id": model_id,
            "yaml_path": yaml_by_model.get(model_id, ""),
            "target_product": str(case["target_product"]),
            "target_variant": str(case["target_variant"]),
            "T_C": float(case["T_C"]),
            "P_Pa": float(case["P_Pa"]),
            "case_hash": compute_case_hash(case),
            "status": "pending",
            "runtime_seconds": np.nan,
            "error_message": "",
        })
    fresh = pd.DataFrame(rows, columns=RUN_MANIFEST_COLUMNS)

    manifest_path = study_dir / "run_manifest.csv"
    if manifest_path.exists():
        prior = pd.read_csv(manifest_path)
        if "case_id" in prior.columns:
            prior = prior.drop_duplicates("case_id").set_index("case_id")
            # Carry prior status only when case_id AND case_hash match; a case_id
            # whose case_hash changed (e.g. edited inventory) is reset to pending.
            if "case_hash" in prior.columns:
                prior_hash = fresh["case_id"].map(prior["case_hash"]).astype(str)
                match = pd.Series(fresh["case_hash"].astype(str).values == prior_hash.values,
                                  index=fresh.index)
            else:
                match = pd.Series(False, index=fresh.index)
            for col in ("status", "runtime_seconds", "error_message"):
                if col in prior.columns:
                    mapped = fresh["case_id"].map(prior[col])
                    carried = mapped.where(mapped.notna(), fresh[col])
                    fresh[col] = carried.where(match, fresh[col])
    fresh["status"] = fresh["status"].fillna("pending")
    fresh["error_message"] = fresh["error_message"].fillna("").astype(str)
    fresh.to_csv(manifest_path, index=False)
    return fresh


# --------------------------------------------------------------------------- #
# Study-local YAML path resolution (review §3.5)
# --------------------------------------------------------------------------- #
def resolve_model_yaml_path(study_dir: str | Path, yaml_path: str | Path) -> Path:
    """Resolve a manifest ``yaml_path`` study-locally.

    Tries, in order: ``study_dir/yaml_path``; ``yaml_path`` as-is (absolute);
    ``study_dir/models/<name>``. Does not reuse the base ``models/single_product``
    resolvers.
    """
    study_dir = Path(study_dir)
    yaml_path = Path(yaml_path)
    candidates = [
        study_dir / yaml_path,
        yaml_path,
        study_dir / "models" / yaml_path.name,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


# --------------------------------------------------------------------------- #
# Running
# --------------------------------------------------------------------------- #
def _failed_rows_for_case(scenario_id: str, target_variant: str, temperature_C: float,
                          pressure_Pa: float, scenario_cfg: Dict[str, Any],
                          model_mode: str, error_message: str) -> pd.DataFrame:
    """Synthesize raw-long rows when a case raises before/within the solver.

    Mirrors the base raw-long schema so downstream code is unaffected; X_eq = NaN.
    run_single_yaml_case catches *solver* exceptions internally, but YAML loading
    or mole normalization can raise earlier — this keeps one failure non-fatal.
    """
    T_K = float(temperature_C) + 273.15
    initial = scenario_cfg.get("initial_moles", {}) or {}
    rows = [{
        "scenario": scenario_id,
        "model_mode": model_mode,
        "yaml_file": "",
        "target_product": target_variant,
        "T_C": float(temperature_C),
        "T_K": T_K,
        "P_Pa": float(pressure_Pa),
        "species": str(sp),
        "X_initial": np.nan,
        "X_eq": np.nan,
        "initial_moles": float(amount),
        "solver_status": "failed",
        "error_message": error_message,
    } for sp, amount in initial.items()]
    if not rows:  # no inventory listed; still emit a marker row
        rows = [{
            "scenario": scenario_id, "model_mode": model_mode, "yaml_file": "",
            "target_product": target_variant, "T_C": float(temperature_C), "T_K": T_K,
            "P_Pa": float(pressure_Pa), "species": target_variant, "X_initial": np.nan,
            "X_eq": np.nan, "initial_moles": np.nan, "solver_status": "failed",
            "error_message": error_message,
        }]
    return pd.DataFrame(rows)


def merge_raw_results_with_design(raw_long_df: pd.DataFrame,
                                  design_df: pd.DataFrame) -> pd.DataFrame:
    """Attach design variables + identifiers to every raw-long row (review §14.2).

    Raw rows carry the case id in ``scenario`` (== ``scenario_id`` == ``case_id``).
    """
    if raw_long_df.empty:
        return raw_long_df
    merge_cols = [c for c in _DESIGN_MERGE_COLUMNS if c in design_df.columns]
    design_keyed = design_df[merge_cols].drop_duplicates("case_id")
    raw = raw_long_df.copy()
    raw["case_id"] = raw["scenario"].astype(str)
    merged = raw.merge(design_keyed, on="case_id", how="left", suffixes=("", "_design"))
    return merged


def select_cases_to_run(manifest_df: pd.DataFrame, design_df: pd.DataFrame,
                        only_failed: bool = False, force: bool = False,
                        substudy: Optional[str] = None,
                        limit: Optional[int] = None) -> pd.DataFrame:
    """Return the manifest rows to execute this invocation (single source of truth).

    - resume (default): cases whose status is not ``ok`` (``pending``/``failed``).
    - ``only_failed``: cases whose status is ``failed`` only.
    - ``force``: every selected case, regardless of status.

    ``substudy`` and ``limit`` compose with all three.
    """
    selected = manifest_df
    if substudy is not None and "substudy_id" in design_df.columns:
        sub_cases = set(design_df.loc[design_df["substudy_id"] == substudy, "case_id"].astype(str))
        selected = selected[selected["case_id"].astype(str).isin(sub_cases)]
    status = selected["status"].astype(str)
    if force:
        to_run = selected
    elif only_failed:
        to_run = selected[status == "failed"]
    else:  # resume-by-default
        to_run = selected[status != "ok"]
    if limit is not None:
        to_run = to_run.head(int(limit))
    return to_run


def run_sensitivity_manifest(
    run_manifest_df: pd.DataFrame,
    scenarios: Dict[str, Any],
    design_df: pd.DataFrame,
    study_dir: str | Path,
    output_long_csv: str | Path,
    phase_name: str = "aqueous",
    solver: str = "vcs",
    max_steps: int = 100_000,
    model_mode: str = "single_product_sensitivity",
    only_failed: bool = False,
    limit: Optional[int] = None,
    substudy: Optional[str] = None,
    force: bool = False,
    progress: bool = True,
) -> pd.DataFrame:
    """Execute selected cases, persisting status, and write the merged raw long CSV.

    Resume-by-default: only cases whose status is not ``ok`` run (``--force`` reruns
    all, ``--only-failed`` reruns failures). Reliability (review §14.1): one failing
    case is recorded (``status=failed``) and never aborts the run; ``run_manifest.csv``
    is updated in place after each case. The raw long CSV write is **merge-safe**:
    rerun cases replace their own prior rows; all other completed cases are preserved.
    """
    from equilibrium_runner import run_single_yaml_case  # imports cantera lazily

    study_dir = Path(study_dir)
    manifest = run_manifest_df.copy()
    # A manifest re-read from CSV has an all-empty error_message as float64 (NaN);
    # keep status/error_message object so in-place "" / status writes don't upcast-fail.
    manifest["status"] = manifest["status"].astype(object)
    if "error_message" in manifest.columns:
        manifest["error_message"] = manifest["error_message"].fillna("").astype(str)
    manifest_path = study_dir / "run_manifest.csv"

    selected = select_cases_to_run(manifest, design_df, only_failed=only_failed,
                                   force=force, substudy=substudy, limit=limit)
    if force and len(selected):
        # Reset the rerun set so an interrupted --force run shows pending, not stale ok.
        reset_mask = manifest["case_id"].astype(str).isin(set(selected["case_id"].astype(str)))
        manifest.loc[reset_mask, "status"] = "pending"
        manifest.loc[reset_mask, "runtime_seconds"] = np.nan
        manifest.loc[reset_mask, "error_message"] = ""
        manifest.to_csv(manifest_path, index=False)

    manifest_index = {str(cid): i for i, cid in enumerate(manifest["case_id"])}
    raw_frames: List[pd.DataFrame] = []
    total = len(selected)

    for n, (_, case) in enumerate(selected.iterrows(), start=1):
        case_id = str(case["case_id"])
        scenario_id = str(case["scenario_id"])
        scenario_cfg = scenarios["scenarios"][scenario_id]
        yaml_path = resolve_model_yaml_path(study_dir, case["yaml_path"])
        target_variant = str(case["target_variant"])

        start = time.perf_counter()
        try:
            case_df = run_single_yaml_case(
                yaml_path=yaml_path,
                scenario_id=case_id,
                scenario_cfg=scenario_cfg,
                target_product=target_variant,
                temperature_C=float(case["T_C"]),
                pressure_Pa=float(case["P_Pa"]),
                model_mode=model_mode,
                phase_name=phase_name,
                solver=solver,
                max_steps=max_steps,
            )
            status = "ok" if (case_df["solver_status"] == "ok").all() else "failed"
            error_message = "" if status == "ok" else str(case_df["error_message"].iloc[0])
        except Exception as exc:  # YAML load / normalization / Cantera import, etc.
            error_message = repr(exc)
            status = "failed"
            case_df = _failed_rows_for_case(
                case_id, target_variant, float(case["T_C"]), float(case["P_Pa"]),
                scenario_cfg, model_mode, error_message)
        runtime = time.perf_counter() - start

        case_df = case_df.copy()
        case_df["case_id"] = case_id
        case_df["model_id"] = str(case["model_id"])
        # Restore the stable semantic identity: the solver measured the variant
        # species (target_product was set to it on input), but the output keeps the
        # original product name; only target_variant holds the pseudo-species.
        case_df["target_product"] = str(case["target_product"])
        case_df["target_variant"] = str(case["target_variant"])
        case_df["runtime_seconds"] = runtime
        raw_frames.append(case_df)

        # Persist status in place (restartable; partial run always inspectable).
        idx = manifest_index[case_id]
        manifest.at[idx, "status"] = status
        manifest.at[idx, "runtime_seconds"] = runtime
        manifest.at[idx, "error_message"] = error_message
        manifest.to_csv(manifest_path, index=False)

        if progress:
            ok = sum(f["solver_status"].iloc[0] == "ok" for f in raw_frames)
            print(f"\r  [{n}/{total}] {case_id} {status:<6} "
                  f"({ok} ok) {runtime:.3f}s", end="", flush=True)
    if progress and total:
        print()

    raw_long = pd.concat(raw_frames, ignore_index=True) if raw_frames else pd.DataFrame()
    merged = merge_raw_results_with_design(raw_long, design_df)
    if "error_message" in merged.columns:
        merged["error_message"] = merged["error_message"].fillna("").astype(str)
    output_long_csv = Path(output_long_csv)
    output_long_csv.parent.mkdir(parents=True, exist_ok=True)

    if merged.empty:
        # Nothing ran this invocation — never clobber prior results.
        if output_long_csv.exists():
            return pd.read_csv(output_long_csv)
        return merged

    # Merge-safe write: drop the rerun cases' prior rows, append the fresh ones, so
    # the file stays the union of all completed cases with no duplicates (§14.2).
    selected_case_ids = {str(c) for c in merged["case_id"]}
    if output_long_csv.exists():
        old = pd.read_csv(output_long_csv)
        if "case_id" in old.columns:
            old = old[~old["case_id"].astype(str).isin(selected_case_ids)]
        # Adopt the current canonical schema so a prior file's stale columns (e.g.
        # a pre-fix model_id_design) are not reintroduced via the column union.
        old = old.reindex(columns=merged.columns)
        merged = pd.concat([old, merged], ignore_index=True)
    merged.to_csv(output_long_csv, index=False)
    return merged


# --------------------------------------------------------------------------- #
# Provenance (review §5.2)
# --------------------------------------------------------------------------- #
def _file_sha256(path: str | Path) -> Optional[str]:
    path = Path(path)
    if not path.exists():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _relpath(path: str | Path, repo_root: str | Path) -> str:
    """Repo-relative path string when possible, else the path as given (portable)."""
    try:
        return str(Path(path).resolve().relative_to(Path(repo_root).resolve()))
    except Exception:
        return str(path)


def _git_sha(repo_root: str | Path) -> Optional[str]:
    import subprocess
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(repo_root),
            capture_output=True, text=True, timeout=5)
        return out.stdout.strip() or None
    except Exception:
        return None


def _package_version(name: str) -> Optional[str]:
    try:
        module = __import__(name)
        return getattr(module, "__version__", None)
    except Exception:
        return None


def prior_seed_hash(study_dir: str | Path) -> Optional[str]:
    """Return the gibbs-seed sha256 recorded in a prior ``run_provenance.json``."""
    path = Path(study_dir) / "run_provenance.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text()).get("gibbs_seed_sha256")
    except Exception:
        return None


def seed_is_stale(study_dir: str | Path, seed_path: str | Path) -> bool:
    """True if a prior run recorded a *different* gibbs-seed hash (coeffs may be stale)."""
    prior = prior_seed_hash(study_dir)
    if prior is None:
        return False
    return _file_sha256(seed_path) != prior


def write_run_provenance(study_dir: str | Path, config_path: str | Path,
                         config: Dict[str, Any], repo_root: str | Path = ".") -> Path:
    """Write ``run_provenance.json`` (versions, hashes, timestamp)."""
    study_dir = Path(study_dir)
    repo_root = Path(repo_root)
    species_csv = repo_root / config.get("species_files", {}).get("species_csv", "")
    gibbs_csv = repo_root / config.get("species_files", {}).get("gibbs_seed_wide_csv", "")
    gibbs_seed_sha256 = _file_sha256(gibbs_csv)
    provenance = {
        "study_id": config.get("study", {}).get("study_id"),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_sha": _git_sha(repo_root),
        "versions": {
            "cantera": _package_version("cantera"),
            "numpy": _package_version("numpy"),
            "pandas": _package_version("pandas"),
            "scipy": _package_version("scipy"),
        },
        "config_path": str(config_path),
        "config_relpath": _relpath(config_path, repo_root),
        "config_sha256": _file_sha256(config_path),
        "gibbs_seed_sha256": gibbs_seed_sha256,
        "input_hashes": {
            "species_csv": _file_sha256(species_csv),
            "gibbs_seed_wide_csv": gibbs_seed_sha256,
        },
    }
    out = study_dir / "run_provenance.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    return out
