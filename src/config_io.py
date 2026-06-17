"""Configuration and path helpers for the Cantera equilibrium workflow."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List

import pandas as pd
import yaml


REQUIRED_SPECIES_COLUMNS = {
    "species_key",
    "cantera_name",
    "chnosz_name",
    "formula",
    "state",
    "molar_volume_cm3_mol",
    "role",
}


@dataclass(frozen=True)
class ProjectPaths:
    root: Path
    inputs: Path
    src: Path
    data_raw: Path
    data_processed: Path
    data_results: Path
    models_single: Path
    figures_fit: Path
    figures_equilibrium: Path
    figures_diagnostics: Path


def get_project_paths(project_root: str | Path = ".") -> ProjectPaths:
    """Return the canonical project paths and create them if missing."""
    root = Path(project_root).resolve()
    paths = ProjectPaths(
        root=root,
        inputs=root / "inputs",
        src=root / "src",
        data_raw=root / "data" / "raw",
        data_processed=root / "data" / "processed",
        data_results=root / "data" / "results",
        models_single=root / "models" / "single_product",
        figures_fit=root / "figures" / "fit_diagnostics",
        figures_equilibrium=root / "figures" / "equilibrium",
        figures_diagnostics=root / "figures" / "diagnostics",
    )
    ensure_directories(paths)
    return paths


def ensure_directories(paths: ProjectPaths) -> None:
    """Create all directories in a ProjectPaths object."""
    for value in paths.__dict__.values():
        if isinstance(value, Path):
            value.mkdir(parents=True, exist_ok=True)


def load_species_metadata(path: str | Path) -> pd.DataFrame:
    """Load and validate species metadata from CSV.

    Required columns:
        species_key, cantera_name, chnosz_name, formula, state,
        molar_volume_cm3_mol, role
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Species metadata file not found: {path}")
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    missing = REQUIRED_SPECIES_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Species metadata is missing columns: {sorted(missing)}")

    df = df.copy()
    for col in ["species_key", "cantera_name", "chnosz_name", "formula", "state", "role"]:
        df[col] = df[col].astype(str).str.strip()
    df["molar_volume_cm3_mol"] = pd.to_numeric(df["molar_volume_cm3_mol"], errors="coerce")

    if df["species_key"].duplicated().any():
        dupes = df.loc[df["species_key"].duplicated(), "species_key"].tolist()
        raise ValueError(f"Duplicate species_key values found: {dupes}")
    if df["cantera_name"].duplicated().any():
        dupes = df.loc[df["cantera_name"].duplicated(), "cantera_name"].tolist()
        raise ValueError(f"Duplicate cantera_name values found: {dupes}")
    if df["molar_volume_cm3_mol"].isna().any():
        bad = df.loc[df["molar_volume_cm3_mol"].isna(), "species_key"].tolist()
        raise ValueError(f"Missing/non-numeric molar volumes for: {bad}")
    return df


def load_scenarios(path: str | Path) -> Dict[str, Any]:
    """Load scenarios YAML and validate basic structure.

    Each scenario must define ``initial_moles``. Optional keys include:
    - ``allowed_species``: overrides the base species list entirely.
    - ``extra_allowed_species``: appended to the initial_moles species list.
    - ``target_products``: optional list of Cantera species names that restricts
      which target products are modelled for this scenario. When omitted the
      global TARGET_PRODUCTS setting in the notebook is used as the fallback.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Scenario YAML not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    if "scenarios" not in data or not isinstance(data["scenarios"], dict):
        raise ValueError("Scenario YAML must contain a top-level 'scenarios:' mapping.")
    for scenario_id, cfg in data["scenarios"].items():
        if "initial_moles" not in cfg:
            raise ValueError(f"Scenario '{scenario_id}' is missing 'initial_moles'.")
        if not isinstance(cfg["initial_moles"], dict):
            raise ValueError(f"Scenario '{scenario_id}' initial_moles must be a mapping.")
        for species, amount in cfg["initial_moles"].items():
            try:
                float(amount)
            except Exception as exc:
                raise ValueError(
                    f"Scenario '{scenario_id}' has a non-numeric initial mole amount for {species!r}: {amount!r}"
                ) from exc
        if "target_products" in cfg and cfg["target_products"] is not None:
            if not isinstance(cfg["target_products"], list):
                raise ValueError(
                    f"Scenario '{scenario_id}' target_products must be a list of Cantera species names."
                )
    return data


def list_target_products(species_df: pd.DataFrame, target_roles: Iterable[str] = ("product", "target")) -> List[str]:
    """Return Cantera species names whose role marks them as a product/target."""
    roles = {r.lower() for r in target_roles}
    mask = species_df["role"].str.lower().isin(roles)
    return species_df.loc[mask, "cantera_name"].tolist()


def list_base_species(species_df: pd.DataFrame, base_roles: Iterable[str] = ("solvent", "reactant", "additive")) -> List[str]:
    """Return Cantera species names whose role marks them as starting/base species."""
    roles = {r.lower() for r in base_roles}
    mask = species_df["role"].str.lower().isin(roles)
    return species_df.loc[mask, "cantera_name"].tolist()


def scenario_ids(scenarios: Dict[str, Any]) -> List[str]:
    return list(scenarios.get("scenarios", {}).keys())


def target_products_for_scenario(scenario_cfg: Dict[str, Any], fallback_products: List[str]) -> List[str]:
    """Return the target products for one scenario.

    If the scenario defines a ``target_products`` list, that list is used.
    Otherwise the global ``fallback_products`` (from the notebook's
    TARGET_PRODUCTS setting) is returned.

    This lets each scenario restrict which products are modelled without
    requiring separate notebooks or separate species master files.
    """
    if "target_products" in scenario_cfg and scenario_cfg["target_products"]:
        return [str(s) for s in scenario_cfg["target_products"]]
    return list(fallback_products)


def base_species_for_scenario(scenario_cfg: Dict[str, Any]) -> List[str]:
    """Return the base/allowed starting species for a scenario.

    Rules:
    - If ``allowed_species`` is provided, it overrides automatic construction.
    - Otherwise, include species listed in ``initial_moles``.
    - Also include any species listed in ``extra_allowed_species``.

    This function deliberately avoids special-casing NH3 or any other molecule.
    If NH3 belongs in a run, put ``NH3(aq)`` in the scenario's ``initial_moles``
    or ``extra_allowed_species``.
    """
    if "allowed_species" in scenario_cfg and scenario_cfg["allowed_species"] is not None:
        return list(dict.fromkeys([str(s) for s in scenario_cfg["allowed_species"]]))
    species = [str(s) for s in scenario_cfg.get("initial_moles", {}).keys()]
    species.extend([str(s) for s in scenario_cfg.get("extra_allowed_species", []) or []])
    return list(dict.fromkeys(species))
