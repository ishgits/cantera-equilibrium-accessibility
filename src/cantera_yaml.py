"""Cantera YAML generation from species metadata and NASA9 coefficients."""
from __future__ import annotations

from pathlib import Path
from typing import List, Sequence

import pandas as pd
import yaml

from formula_tools import parse_formula
from thermo_fit import COEFF_COLUMNS, load_coefficients, coefficients_for_species


def _float_list(values: Sequence[float]) -> List[float]:
    return [float(v) for v in values]


def _species_block(meta_row: pd.Series, coeff_df: pd.DataFrame) -> dict:
    cantera_name = meta_row["cantera_name"]
    low, high = coefficients_for_species(coeff_df, cantera_name)
    temp_ranges = [float(low["T_low_K"]), float(low["T_high_K"]), float(high["T_high_K"])]
    data_low = _float_list([low[c] for c in COEFF_COLUMNS])
    data_high = _float_list([high[c] for c in COEFF_COLUMNS])
    return {
        "name": cantera_name,
        "composition": parse_formula(meta_row["formula"]),
        "thermo": {
            "model": "NASA9",
            "temperature-ranges": temp_ranges,
            "reference-pressure": 101325.0,
            "data": [data_low, data_high],
        },
        "equation-of-state": {
            "model": "constant-volume",
            "molar-volume": float(meta_row["molar_volume_cm3_mol"]),
        },
    }


def _yaml_document(phase_species: Sequence[str], species_blocks: Sequence[dict], phase_name: str = "aqueous") -> dict:
    return {
        "units": {
            "length": "cm",
            "quantity": "mol",
            "activation-energy": "cal/mol",
            "energy": "J",
            "mass": "g",
            "pressure": "Pa",
            "temperature": "K",
            "time": "s",
        },
        "phases": [
            {
                "name": phase_name,
                "thermo": "ideal-condensed",
                "species": list(phase_species),
            }
        ],
        "species": list(species_blocks),
    }


def write_yaml_document(doc: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(doc, f, sort_keys=False, width=120)
    return path


def _safe_stem(cantera_name: str) -> str:
    return cantera_name.replace("(aq)", "").replace("(l)", "").replace(" ", "_").replace("/", "_")



def _portable_path(path: Path, output_dir: Path) -> str:
    """Return a repo-relative path when possible, otherwise an absolute path."""
    path = Path(path)
    try:
        project_root = output_dir.resolve().parents[1]
        return str(path.resolve().relative_to(project_root))
    except Exception:
        return str(path.resolve())


def generate_single_product_yamls(
    species_df: pd.DataFrame,
    coeffs_csv: str | Path,
    output_dir: str | Path,
    base_species: Sequence[str],
    target_products: Sequence[str],
    phase_name: str = "aqueous",
) -> List[Path]:
    """Generate one YAML per target product."""
    coeff_df = load_coefficients(coeffs_csv)
    meta = species_df.set_index("cantera_name", drop=False)
    output_dir = Path(output_dir)
    written = []
    for target in target_products:
        phase_species = list(dict.fromkeys(list(base_species) + [target]))
        missing_meta = [s for s in phase_species if s not in meta.index]
        if missing_meta:
            raise KeyError(f"Species not found in species metadata: {missing_meta}")
        blocks = [_species_block(meta.loc[s], coeff_df) for s in phase_species]
        doc = _yaml_document(phase_species, blocks, phase_name=phase_name)
        written.append(write_yaml_document(doc, output_dir / f"{_safe_stem(target)}.yaml"))
    return written


def validate_yaml_with_cantera(yaml_path: str | Path, phase_name: str = "aqueous") -> bool:
    """Try loading a YAML file with Cantera.

    Returns True if it loads. Raises ImportError if Cantera is unavailable.
    """
    try:
        import cantera as ct
    except ImportError as exc:
        raise ImportError("Cantera is not installed in this environment.") from exc
    ct.Solution(str(yaml_path), phase_name)
    return True


def generate_single_product_yamls_for_scenarios(
    species_df: pd.DataFrame,
    coeffs_csv: str | Path,
    scenarios: dict,
    output_dir: str | Path,
    target_products: Sequence[str],
    phase_name: str = "aqueous",
) -> pd.DataFrame:
    """Generate scenario-specific single-product YAMLs and return a manifest.

    Scenario-specific YAMLs matter because absent species must be excluded from the
    phase; setting an initial mole amount to zero does not prevent Cantera from
    forming that species if it is present in the YAML.

    If a scenario defines a ``target_products`` list in scenarios.yaml, only those
    products are modelled for that scenario. Otherwise the ``target_products``
    argument (the notebook's global TARGET_PRODUCTS setting) is used as the fallback.
    """
    from config_io import base_species_for_scenario, target_products_for_scenario

    output_dir = Path(output_dir)
    rows = []
    for scenario_id, scenario_cfg in scenarios["scenarios"].items():
        base_species = base_species_for_scenario(scenario_cfg)
        scenario_targets = target_products_for_scenario(scenario_cfg, target_products)
        scenario_dir = output_dir / scenario_id
        written = generate_single_product_yamls(
            species_df=species_df,
            coeffs_csv=coeffs_csv,
            output_dir=scenario_dir,
            base_species=base_species,
            target_products=scenario_targets,
            phase_name=phase_name,
        )
        for path in written:
            target = next((tp for tp in scenario_targets if _safe_stem(tp) == path.stem), path.stem)
            rows.append({
                "scenario": scenario_id,
                "model_mode": "single_product",
                "yaml_path": _portable_path(path, output_dir),
                "yaml_file": path.name,
                "target_product": target,
                "allowed_base_species": ";".join(base_species),
            })
    manifest = pd.DataFrame(rows)
    manifest_path = output_dir / "single_product_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    return manifest
