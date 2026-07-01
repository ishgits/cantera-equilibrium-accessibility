"""Sensitivity study design generation (Phase 1).

Turns a user-edited ``study_config.yaml`` into:

- ``design_matrix.csv``        — one row per intended simulation case.
- ``generated_scenarios.yaml`` — one Cantera scenario per case (with explicit
  ``allowed_species`` phase membership and full inventory keys, zeros retained).
- ``thermo_offsets.csv``       — bookkeeping for any Gibbs-offset (ΔG) substudy.

This module is pure pandas/PyYAML and never imports or runs Cantera. It reuses
the validated v1.0 helpers in ``config_io`` (notably the ``allowed_species``
contract) and does not modify any base-workflow behavior.

Key design invariants:

- Phase membership is explicit via ``allowed_species``; ``0.0`` initial moles are
  meaningful (present-but-zero) and are retained, never tidied away.
- ``case_id`` is the canonical run key and ``scenario_id == case_id``.
- ΔG variants get pseudo-species names here; the analytic a7-shift that gives them
  thermodynamics is applied at run time by
  ``sensitivity_thermo.shift_coeffs_by_gibbs``. This module only needs the
  variant *names*.
"""
from __future__ import annotations

from itertools import product
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import yaml

# Reuse the validated v1.0 loader so the species contract stays single-sourced.
from config_io import load_species_metadata


# Canonical design-matrix column order (architecture §8.1). Inventory columns are
# the alanine-MVP set; they are derived generically from cantera names so the
# layout extends cleanly to other molecules.
DESIGN_COLUMNS = [
    "case_id",
    "study_id",
    "substudy_id",
    "target_product",
    "target_variant",
    "H2O_mol",
    "HCN_mol",
    "C2H2_mol",
    "NH3_mol",
    "C2H2_over_HCN",
    "deltaG_offset_kJ_mol",
    "T_C",
    "P_Pa",
    "scenario_id",
    "model_id",
]

# substudy_id -> short code used in case ids (ALA_INV_000001, ALA_DG_..., ...).
SUBSTUDY_CODES = {
    "inventory_landscape": "INV",
    "deltaG_sweep": "DG",
    "nh3_deltaG_landscape": "NH3DG",
}

# Recognised substudies and the axes each one must declare.
KNOWN_SUBSTUDIES = {"inventory_landscape", "deltaG_sweep", "nh3_deltaG_landscape"}
# C2H2/HCN is the always-required inventory axis; NH3(aq) is optional so a study can
# exclude ammonia entirely (a 1-D C2H2/HCN sweep). nh3_deltaG_landscape inherently
# needs NH3, so it is simply disabled in NH3-excluded studies.
REQUIRED_AXES = {
    "inventory_landscape": ["C2H2_over_HCN"],
    "deltaG_sweep": ["offsets_kJ_mol"],
    "nh3_deltaG_landscape": ["NH3(aq)", "deltaG_offset_kJ_mol"],
}

# Temperature tolerance (K) when checking the run T against the Gibbs-seed range.
# The owner's config runs at 0 C = 273.15 K, a hair below the seed min 273.16 K;
# a small tolerance keeps that valid without weakening the check meaningfully.
_T_RANGE_TOL_K = 1.0


class StudyConfigError(ValueError):
    """Raised for user-input errors in a study config.

    Carries a plain-English message intended to be shown directly to a
    non-coder; the CLI prints ``str(exc)`` instead of a traceback.
    """


# --------------------------------------------------------------------------- #
# Config loading
# --------------------------------------------------------------------------- #
def load_study_config(path: str | Path) -> Dict[str, Any]:
    """Read a study YAML and perform light structural validation.

    Deeper, species-aware checks live in :func:`validate_study`.
    """
    path = Path(path)
    if not path.exists():
        raise StudyConfigError(
            f"Study config not found: {path}. Create it (see "
            f"studies/_template/study_config.yaml) or check the --config path."
        )
    with path.open("r", encoding="utf-8") as f:
        config = yaml.safe_load(f) or {}
    if not isinstance(config, dict):
        raise StudyConfigError(f"Study config {path} did not parse to a mapping.")
    for section in ("study", "mode", "species_files", "base_conditions", "sweeps"):
        if section not in config:
            raise StudyConfigError(
                f"Study config {path} is missing the required top-level "
                f"'{section}:' section."
            )
    if "study_id" not in config["study"]:
        raise StudyConfigError("Study config is missing 'study.study_id'.")
    return config


# --------------------------------------------------------------------------- #
# Small config accessors
# --------------------------------------------------------------------------- #
def _target_product(config: Dict[str, Any]) -> str:
    targets = config.get("mode", {}).get("target_products") or []
    if not targets:
        raise StudyConfigError(
            "mode.target_products is empty — name the product whose "
            "accessibility you are mapping, e.g. 'Alanine(aq)'."
        )
    return str(targets[0])


def _allowed_reactants(config: Dict[str, Any]) -> List[str]:
    allowed = config.get("model", {}).get("allowed_species") or []
    if not allowed:
        raise StudyConfigError(
            "model.allowed_species is empty — list the reactant/solvent species "
            "that make up the Cantera phase (the target is added automatically)."
        )
    return [str(s) for s in allowed]


def _enabled_sweeps(config: Dict[str, Any]) -> Dict[str, Any]:
    sweeps = config.get("sweeps") or {}
    return {name: spec for name, spec in sweeps.items() if spec and spec.get("enabled")}


def _mol_column(cantera_name: str) -> str:
    """Map a Cantera species name to its design-matrix inventory column.

    ``H2O(l) -> H2O_mol``, ``NH3(aq) -> NH3_mol``.
    """
    base = str(cantera_name).split("(")[0]
    return f"{base}_mol"


def _is_inventory_col(column: str) -> bool:
    """True for initial-moles columns (``*_mol``) but not ``deltaG_offset_kJ_mol``."""
    return column.endswith("_mol") and not column.endswith("kJ_mol")


def _study_prefix(study_id: str) -> str:
    """Short uppercase prefix for case ids, e.g. ``alanine_mvp -> ALA``."""
    letters = "".join(ch for ch in str(study_id) if ch.isalpha())
    return (letters[:3] or "STD").upper()


# --------------------------------------------------------------------------- #
# Sweep value generation
# --------------------------------------------------------------------------- #
def make_sweep_values(spec: Dict[str, Any]) -> List[float]:
    """Expand a sweep spec into an explicit list of float values.

    Supported ``type`` values:

    - ``linear``   — ``np.linspace(min, max, points)``.
    - ``logspace`` — ``np.geomspace(min, max, points)`` over the *actual* min/max
      values (not exponents); requires ``min > 0``.
    - ``explicit`` — ``values`` taken verbatim.

    Any axis may also carry ``include_values: [...]`` — exact values merged into the
    generated grid (de-duplicated, sorted), e.g. to hit a paper fiducial like 2.1.
    """
    if not isinstance(spec, dict):
        raise StudyConfigError(f"Sweep spec must be a mapping, got: {spec!r}")
    stype = str(spec.get("type", "")).lower()

    if stype == "explicit":
        values = spec.get("values")
        if not values:
            raise StudyConfigError("Explicit sweep needs a non-empty 'values:' list.")
        base = [float(v) for v in values]
    elif stype in ("linear", "logspace"):
        for key in ("min", "max", "points"):
            if key not in spec:
                raise StudyConfigError(f"{stype} sweep is missing '{key}:'.")
        lo, hi = float(spec["min"]), float(spec["max"])
        points = int(spec["points"])
        if points < 1:
            raise StudyConfigError(f"Sweep 'points' must be >= 1, got {points}.")
        if stype == "linear":
            base = [float(v) for v in np.linspace(lo, hi, points)]
        elif lo <= 0:
            raise StudyConfigError(
                "logspace sweep requires min > 0 (values are actual values, not "
                f"exponents); got min={lo}."
            )
        else:
            base = [float(v) for v in np.geomspace(lo, hi, points)]
    else:
        raise StudyConfigError(
            f"Unknown sweep type {spec.get('type')!r}; use 'linear', 'logspace', or "
            "'explicit'."
        )

    include = spec.get("include_values")
    if include:
        base = base + [float(v) for v in include]
        base = sorted({round(v, 12) for v in base})
    return base


# --------------------------------------------------------------------------- #
# ΔG pseudo-species naming
# --------------------------------------------------------------------------- #
def make_thermo_variant_name(base_cantera_name: str, offset_kJ_mol: float) -> str:
    """Name a Gibbs-offset pseudo-species, preserving the ``(state)`` suffix.

    ``Alanine(aq), -50 -> Alanine__dG_m050(aq)``;
    ``Alanine(aq),   0 -> Alanine__dG_000(aq)``;
    ``Alanine(aq),  20 -> Alanine__dG_p020(aq)``.
    """
    name = str(base_cantera_name)
    if "(" in name:
        stem, state = name.split("(", 1)
        state = "(" + state
    else:
        stem, state = name, ""
    value = float(offset_kJ_mol)
    if abs(value - round(value)) > 1e-9:
        raise ValueError(
            f"ΔG offset must be a whole number of kJ/mol; got {value}. Fractional "
            "offsets would round to a colliding pseudo-species name."
        )
    rounded = int(round(value))
    sign = "m" if rounded < 0 else "p" if rounded > 0 else ""
    tag = f"dG_{sign}{abs(rounded):03d}"
    return f"{stem}__{tag}{state}"


# --------------------------------------------------------------------------- #
# Per-substudy design builders
# --------------------------------------------------------------------------- #
def _base_row(config: Dict[str, Any], substudy_id: str) -> Dict[str, Any]:
    bc = config.get("base_conditions", {})
    return {
        "study_id": str(config["study"]["study_id"]),
        "substudy_id": substudy_id,
        "target_product": _target_product(config),
        "T_C": float(bc.get("temperature_C", 25)),
        "P_Pa": float(bc.get("pressure_Pa", 101325)),
        "model_id": "",  # assigned in Phase 2 after model grouping
    }


def _inventory_row(base_row: Dict[str, Any], inventory: Dict[str, float],
                   c2h2_over_hcn: float, deltaG: float, target_variant: str) -> Dict[str, Any]:
    row = dict(base_row)
    row["target_variant"] = target_variant
    row["C2H2_over_HCN"] = float(c2h2_over_hcn)
    row["deltaG_offset_kJ_mol"] = float(deltaG)
    for species, amount in inventory.items():
        row[_mol_column(species)] = float(amount)
    return row


def build_inventory_landscape_design(config: Dict[str, Any]) -> pd.DataFrame:
    """NH3 × (C2H2/HCN) grid at base thermo (no Gibbs offset)."""
    spec = config["sweeps"]["inventory_landscape"]
    variables = spec.get("variables", {})
    fixed = config.get("fixed_inventory", {})
    hcn = float(fixed.get("HCN(aq)", 0.0))
    h2o = float(fixed.get("H2O(l)", 0.0))

    # NH3 axis is optional: when absent (e.g. an NH3-excluded study) the inventory
    # landscape is a 1-D sweep over C2H2/HCN with no ammonia in the system.
    nh3_values = make_sweep_values(variables["NH3(aq)"]) if "NH3(aq)" in variables else [0.0]
    ratio_values = make_sweep_values(variables["C2H2_over_HCN"])
    target = _target_product(config)
    base = _base_row(config, "inventory_landscape")

    rows = []
    for nh3, ratio in product(nh3_values, ratio_values):
        inventory = {
            "H2O(l)": h2o,
            "HCN(aq)": hcn,
            "C2H2(aq)": ratio * hcn,
            "NH3(aq)": nh3,
        }
        rows.append(_inventory_row(base, inventory, ratio, 0.0, target))
    return pd.DataFrame(rows)


def build_deltaG_sweep_design(config: Dict[str, Any]) -> pd.DataFrame:
    """Fixed-inventory Gibbs-offset sweep on the target species."""
    spec = config["sweeps"]["deltaG_sweep"]
    fixed = spec.get("fixed_inventory", {})
    hcn = float(fixed.get("HCN(aq)", 0.0))
    offsets = make_sweep_values(spec["offsets_kJ_mol"])
    target = _target_product(config)
    base = _base_row(config, "deltaG_sweep")

    ratio = (float(fixed.get("C2H2(aq)", 0.0)) / hcn) if hcn else 0.0
    inventory = {k: float(v) for k, v in fixed.items()}

    rows = []
    for offset in offsets:
        variant = make_thermo_variant_name(target, offset)
        rows.append(_inventory_row(base, inventory, ratio, offset, variant))
    return pd.DataFrame(rows)


def build_nh3_deltaG_landscape_design(config: Dict[str, Any]) -> pd.DataFrame:
    """NH3 × Gibbs-offset grid at fixed HCN/C2H2."""
    spec = config["sweeps"]["nh3_deltaG_landscape"]
    variables = spec.get("variables", {})
    fixed = spec.get("fixed_inventory", {})
    hcn = float(fixed.get("HCN(aq)", 0.0))
    c2h2 = float(fixed.get("C2H2(aq)", 0.0))
    h2o = float(fixed.get("H2O(l)", 0.0))
    ratio = (c2h2 / hcn) if hcn else 0.0

    nh3_values = make_sweep_values(variables["NH3(aq)"])
    offsets = make_sweep_values(variables["deltaG_offset_kJ_mol"])
    target = _target_product(config)
    base = _base_row(config, "nh3_deltaG_landscape")

    rows = []
    for nh3, offset in product(nh3_values, offsets):
        inventory = {
            "H2O(l)": h2o,
            "HCN(aq)": hcn,
            "C2H2(aq)": c2h2,
            "NH3(aq)": nh3,
        }
        variant = make_thermo_variant_name(target, offset)
        rows.append(_inventory_row(base, inventory, ratio, offset, variant))
    return pd.DataFrame(rows)


_SUBSTUDY_BUILDERS = {
    "inventory_landscape": build_inventory_landscape_design,
    "deltaG_sweep": build_deltaG_sweep_design,
    "nh3_deltaG_landscape": build_nh3_deltaG_landscape_design,
}


def build_full_design_matrix(config: Dict[str, Any]) -> pd.DataFrame:
    """Concatenate all enabled substudy designs into one design matrix.

    Assigns ``case_id`` / ``scenario_id`` (``{PREFIX}_{CODE}_{NNNNNN}``) and
    returns columns in the canonical :data:`DESIGN_COLUMNS` order.
    """
    enabled = _enabled_sweeps(config)
    if not enabled:
        raise StudyConfigError(
            "No sweeps are enabled — set 'enabled: true' on at least one substudy "
            "under 'sweeps:'."
        )

    prefix = _study_prefix(config["study"]["study_id"])
    frames = []
    # Iterate in a stable, documented order regardless of YAML key order.
    for substudy_id in ("inventory_landscape", "deltaG_sweep", "nh3_deltaG_landscape"):
        if substudy_id not in enabled:
            continue
        df = _SUBSTUDY_BUILDERS[substudy_id](config)
        code = SUBSTUDY_CODES[substudy_id]
        case_ids = [f"{prefix}_{code}_{i:06d}" for i in range(1, len(df) + 1)]
        df = df.copy()
        df.insert(0, "case_id", case_ids)
        df["scenario_id"] = case_ids
        frames.append(df)

    matrix = pd.concat(frames, ignore_index=True)

    # Ensure every canonical column exists, fill inventory gaps with 0.0, order.
    for col in DESIGN_COLUMNS:
        if col not in matrix.columns:
            matrix[col] = 0.0 if col.endswith("_mol") else ""
    inventory_cols = [c for c in matrix.columns if _is_inventory_col(c)]
    matrix[inventory_cols] = matrix[inventory_cols].fillna(0.0)
    return matrix[DESIGN_COLUMNS].copy()


# --------------------------------------------------------------------------- #
# Scenario YAML generation
# --------------------------------------------------------------------------- #
def design_matrix_to_scenarios_yaml(design_df: pd.DataFrame, config: Dict[str, Any],
                                    output_path: str | Path) -> Dict[str, Any]:
    """Write one Cantera scenario per case and return the scenarios mapping.

    Each scenario carries an explicit ``allowed_species`` (reactant phase set +
    the case's ``target_variant``; single-product means exactly one target variant
    per phase) and a full ``initial_moles`` map with zeros retained.
    """
    reactants = _allowed_reactants(config)
    inventory_cols = [c for c in design_df.columns if _is_inventory_col(c)]
    # Map inventory column -> the reactant cantera name it represents.
    col_to_species = {_mol_column(s): s for s in reactants}

    scenarios: Dict[str, Any] = {}
    for _, row in design_df.iterrows():
        scenario_id = str(row["scenario_id"])
        target_variant = str(row["target_variant"])
        # Phase membership: reactants + this case's single target variant.
        allowed = list(dict.fromkeys(reactants + [target_variant]))
        # Initial moles for every reactant column (zeros retained).
        initial_moles = {}
        for col in inventory_cols:
            species = col_to_species.get(col)
            if species is None:
                continue
            initial_moles[species] = float(row[col])
        scenarios[scenario_id] = {
            "description": (
                f"{row['study_id']} | {row['substudy_id']} | {scenario_id}"
            ),
            "target_products": [target_variant],
            "allowed_species": allowed,
            "initial_moles": initial_moles,
        }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump({"scenarios": scenarios}, f, sort_keys=False, default_flow_style=False)
    return scenarios


def build_thermo_offsets_table(config: Dict[str, Any]) -> pd.DataFrame:
    """Bookkeeping table of ΔG pseudo-species (empty if no ΔG substudy enabled).

    Columns: ``thermo_variant_id, base_species, variant_species,
    deltaG_offset_kJ_mol``. The actual a7 coefficient shift is applied by
    ``sensitivity_thermo.shift_coeffs_by_gibbs``.
    """
    target = _target_product(config)
    enabled = _enabled_sweeps(config)
    offsets: set[float] = set()
    if "deltaG_sweep" in enabled:
        offsets.update(make_sweep_values(enabled["deltaG_sweep"]["offsets_kJ_mol"]))
    if "nh3_deltaG_landscape" in enabled:
        offsets.update(make_sweep_values(
            enabled["nh3_deltaG_landscape"]["variables"]["deltaG_offset_kJ_mol"]
        ))

    rows = []
    for offset in sorted(offsets):
        variant = make_thermo_variant_name(target, offset)
        # variant id is the variant name's "dG_xxx" tag.
        variant_id = variant.split("__", 1)[1].split("(")[0]
        rows.append({
            "thermo_variant_id": variant_id,
            "base_species": target,
            "variant_species": variant,
            "deltaG_offset_kJ_mol": float(offset),
        })
    return pd.DataFrame(rows, columns=[
        "thermo_variant_id", "base_species", "variant_species", "deltaG_offset_kJ_mol",
    ])


# --------------------------------------------------------------------------- #
# Output writer (single funnel; future Parquet export is a one-line addition)
# --------------------------------------------------------------------------- #
def save_table(df: pd.DataFrame, path: str | Path) -> Path:
    """Write a tidy table to CSV. Single writer for all design/result tables."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def _study_dir(config: Dict[str, Any]) -> Path:
    return Path(config["study"].get("output_dir", f"studies/{config['study']['study_id']}"))


def write_design_outputs(design_df: pd.DataFrame, config: Dict[str, Any]) -> Dict[str, Path]:
    """Write design_matrix.csv, generated_scenarios.yaml, and thermo_offsets.csv.

    Returns a mapping of logical name -> written path.
    """
    study_dir = _study_dir(config)
    paths: Dict[str, Path] = {}

    paths["design_matrix"] = save_table(design_df, study_dir / "design_matrix.csv")

    scenarios_path = study_dir / "generated_scenarios.yaml"
    design_matrix_to_scenarios_yaml(design_df, config, scenarios_path)
    paths["generated_scenarios"] = scenarios_path

    offsets_df = build_thermo_offsets_table(config)
    if not offsets_df.empty:
        paths["thermo_offsets"] = save_table(offsets_df, study_dir / "thermo_offsets.csv")
    return paths


# --------------------------------------------------------------------------- #
# Step-0 validator (plain-English failures)
# --------------------------------------------------------------------------- #
def _require_integer_offsets(substudy_name: str, values: List[float]) -> None:
    """Raise if any ΔG offset is not a whole kJ/mol (would collide variant names)."""
    for v in values:
        if abs(v - round(v)) > 1e-9:
            raise StudyConfigError(
                f"{substudy_name}: ΔG offset {v} kJ/mol is not a whole number. "
                "Offsets must be integer kJ/mol so each maps to a unique pseudo-species "
                "name (e.g. dG_p020); fractional values can collide. Use an integer "
                "min/max/points combination or an explicit list of whole numbers."
            )


def _seed_temperature_range(gibbs_seed_df: Optional[pd.DataFrame]):
    if gibbs_seed_df is None or gibbs_seed_df.empty:
        return None
    t_col = gibbs_seed_df.columns[0]
    series = pd.to_numeric(gibbs_seed_df[t_col], errors="coerce").dropna()
    if series.empty:
        return None
    return float(series.min()), float(series.max())


def validate_study(config: Dict[str, Any], species_df: pd.DataFrame,
                   gibbs_seed_df: Optional[pd.DataFrame] = None) -> None:
    """Step-0 validation. Raises :class:`StudyConfigError` with human messages.

    Checks (review §5): target/inventory/swept species and ``allowed_species`` all
    present in the species CSV; thresholds numeric & positive; sweep point counts
    >= 1; no negative generated moles; run temperature within the Gibbs-seed range;
    ``output_dir`` creatable; ``case_id`` uniqueness.
    """
    species_names = set(species_df["cantera_name"])
    species_csv = config.get("species_files", {}).get("species_csv", "the species CSV")

    def _require_species(name: str, where: str) -> None:
        if name not in species_names:
            raise StudyConfigError(
                f"{name} appears in {where} but is not in {species_csv} — check "
                f"spelling or add a row for it."
            )

    # 1. Exactly one target product (one study folder per target; see new_study.py).
    targets = config.get("mode", {}).get("target_products") or []
    if len(targets) != 1:
        raise StudyConfigError(
            "This sensitivity workflow supports exactly one target product per study. "
            "Create one study folder per target (recommended), or implement multi-target "
            "expansion first."
        )
    _require_species(str(targets[0]), "mode.target_products")

    # 1b. Substudy names and required axes (before any builder runs).
    for name, spec in _enabled_sweeps(config).items():
        if name not in KNOWN_SUBSTUDIES:
            raise StudyConfigError(
                f"Unknown substudy '{name}' under sweeps. Valid options are: "
                f"{', '.join(sorted(KNOWN_SUBSTUDIES))}."
            )
        variables = spec.get("variables", {})
        for axis in REQUIRED_AXES.get(name, []):
            present = (axis in spec) or (axis in variables)
            if not present:
                raise StudyConfigError(
                    f"Substudy '{name}' is missing its required axis '{axis}'. "
                    f"It needs: {', '.join(REQUIRED_AXES[name])}."
                )

    # 2. Phase membership.
    for sp in _allowed_reactants(config):
        _require_species(sp, "model.allowed_species")

    # 3. Top-level fixed inventory.
    for sp, amount in (config.get("fixed_inventory") or {}).items():
        _require_species(str(sp), "fixed_inventory")
        if float(amount) < 0:
            raise StudyConfigError(f"fixed_inventory[{sp}] is negative ({amount}).")

    # 4. Thresholds numeric & positive.
    thresholds = config.get("thresholds") or {}
    for key in ("formation_X_threshold", "significant_X_threshold"):
        if key in thresholds and not (float(thresholds[key]) > 0):
            raise StudyConfigError(f"thresholds.{key} must be a positive number.")
    for key in ("formation_n_threshold_mol", "balance_tol"):
        if key in thresholds and float(thresholds[key]) < 0:
            raise StudyConfigError(f"thresholds.{key} must be >= 0.")

    # 5. Per-sweep checks.
    enabled = _enabled_sweeps(config)
    if not enabled:
        raise StudyConfigError(
            "No sweeps are enabled — set 'enabled: true' on at least one substudy."
        )
    for name, spec in enabled.items():
        variables = spec.get("variables", {})
        for var_name, var_spec in variables.items():
            values = make_sweep_values(var_spec)  # also validates point count >= 1
            # A variable with a '(' is a real species axis; check membership + sign.
            if "(" in var_name:
                _require_species(var_name, f"the {name} sweep")
                if min(values) < 0:
                    raise StudyConfigError(
                        f"{name} sweep on {var_name} would generate negative moles "
                        f"(min value {min(values)})."
                    )
            elif var_name == "C2H2_over_HCN" and min(values) < 0:
                raise StudyConfigError(
                    f"{name}: C2H2_over_HCN must be >= 0 (min {min(values)})."
                )
            elif var_name == "deltaG_offset_kJ_mol":
                _require_integer_offsets(name, values)
        if "offsets_kJ_mol" in spec:
            _require_integer_offsets(name, make_sweep_values(spec["offsets_kJ_mol"]))
        for sp, amount in (spec.get("fixed_inventory") or {}).items():
            _require_species(str(sp), f"the {name} fixed_inventory")
            if float(amount) < 0:
                raise StudyConfigError(
                    f"{name}.fixed_inventory[{sp}] is negative ({amount})."
                )
        if name == "deltaG_sweep" and "species" in spec:
            if str(spec["species"]) != _target_product(config):
                raise StudyConfigError(
                    f"deltaG_sweep.species ({spec['species']}) must match the target "
                    f"product ({_target_product(config)})."
                )

    # 5b. No swept or inventory species may sit outside the phase (allowed_species),
    # so an excluded species (e.g. NH3) cannot leak back in via a sweep or inventory.
    allowed = set(_allowed_reactants(config))
    for sp in (config.get("fixed_inventory") or {}):
        if str(sp) not in allowed:
            raise StudyConfigError(
                f"fixed_inventory species {sp} is not in model.allowed_species — "
                "add it to the phase or remove it from the inventory.")
    for name, spec in enabled.items():
        for var_name in spec.get("variables", {}):
            if "(" in var_name and var_name not in allowed:
                raise StudyConfigError(
                    f"{name} sweeps {var_name}, which is not in model.allowed_species.")
        for sp in (spec.get("fixed_inventory") or {}):
            if str(sp) not in allowed:
                raise StudyConfigError(
                    f"{name}.fixed_inventory species {sp} is not in model.allowed_species.")

    # 6. Run temperature within Gibbs-seed range.
    t_range = _seed_temperature_range(gibbs_seed_df)
    if t_range is not None:
        t_k = float(config["base_conditions"].get("temperature_C", 25)) + 273.15
        lo, hi = t_range
        if t_k < lo - _T_RANGE_TOL_K or t_k > hi + _T_RANGE_TOL_K:
            raise StudyConfigError(
                f"Equilibrium temperature {t_k:.2f} K is outside the Gibbs-seed data "
                f"range [{lo:.2f}, {hi:.2f}] K — extend the seed data or change "
                "base_conditions.temperature_C."
            )

    # 7. output_dir creatable.
    try:
        _study_dir(config).mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise StudyConfigError(
            f"Cannot create output_dir {_study_dir(config)}: {exc}."
        ) from exc

    # 8. case_id uniqueness (catches builder regressions).
    matrix = build_full_design_matrix(config)
    if matrix["case_id"].duplicated().any():
        dupes = matrix.loc[matrix["case_id"].duplicated(), "case_id"].tolist()
        raise StudyConfigError(f"Generated duplicate case_ids: {dupes[:5]}.")


# --------------------------------------------------------------------------- #
# Convenience: load species/seed referenced by a config
# --------------------------------------------------------------------------- #
def load_species_for_config(config: Dict[str, Any], repo_root: str | Path = ".") -> pd.DataFrame:
    """Load the species CSV referenced by ``species_files.species_csv``."""
    repo_root = Path(repo_root)
    rel = config.get("species_files", {}).get("species_csv")
    if not rel:
        raise StudyConfigError("species_files.species_csv is not set in the config.")
    return load_species_metadata(repo_root / rel)


def load_gibbs_seed_for_config(config: Dict[str, Any], repo_root: str | Path = ".") -> Optional[pd.DataFrame]:
    """Load the wide Gibbs seed CSV.

    Returns ``None`` only when ``species_files.gibbs_seed_wide_csv`` is absent or
    blank. If it is **set but the file does not exist**, raise immediately so a
    path typo surfaces at validation time, not silently mid-run.
    """
    repo_root = Path(repo_root)
    rel = config.get("species_files", {}).get("gibbs_seed_wide_csv")
    if not rel or not str(rel).strip():
        return None
    path = repo_root / rel
    if not path.exists():
        raise StudyConfigError(
            f"Gibbs seed file not found: {path} — check species_files.gibbs_seed_wide_csv."
        )
    from thermo_fit import read_wide_gibbs_csv  # local import to avoid coupling at module load
    return read_wide_gibbs_csv(path)
