"""Run Cantera equilibrium sweeps and save raw mole-fraction outputs."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Sequence

import numpy as np
import pandas as pd


def _load_cantera():
    try:
        import cantera as ct
    except ImportError as exc:
        raise ImportError(
            "Cantera is not installed. Install cantera before running equilibrium simulations."
        ) from exc
    return ct


def normalize_moles_to_X(moles: Dict[str, float], allowed_species: Sequence[str]) -> Dict[str, float]:
    """Convert positive initial moles to mole fractions for species present in a YAML."""
    filtered = {k: float(v) for k, v in moles.items() if k in allowed_species and float(v) > 0}
    total = sum(filtered.values())
    if total <= 0:
        raise ValueError("Initial composition has no positive moles after filtering to YAML species.")
    return {k: v / total for k, v in filtered.items()}


def build_initial_moles(scenario_cfg: Dict[str, Any]) -> Dict[str, float]:
    """Return the scenario's initial mole inventory without molecule-specific logic."""
    return {str(k): float(v) for k, v in scenario_cfg["initial_moles"].items()}


def _result_rows_for_solution(
    sol,
    yaml_path: Path,
    scenario_id: str,
    target_product: str,
    temperature_C: float,
    pressure_Pa: float,
    initial_moles: Dict[str, float],
    X0: Dict[str, float],
    status: str,
    error_message: str,
    model_mode: str,
) -> list[dict]:
    T_K = float(temperature_C) + 273.15
    rows = []
    for sp in sol.species_names:
        rows.append({
            "scenario": scenario_id,
            "model_mode": model_mode,
            "yaml_file": yaml_path.name,
            "target_product": target_product,
            "T_C": float(temperature_C),
            "T_K": T_K,
            "P_Pa": float(pressure_Pa),
            "species": sp,
            "X_initial": float(X0.get(sp, 0.0)),
            "X_eq": float(sol.X[sol.species_index(sp)]) if status == "ok" else np.nan,
            "initial_moles": float(initial_moles.get(sp, 0.0)),
            "solver_status": status,
            "error_message": error_message,
        })
    return rows



def resolve_yaml_path_from_manifest_row(case) -> Path:
    """Resolve a manifest YAML path robustly after moving/copying the repo.

    Fresh manifests usually contain valid absolute paths. If the repo was moved
    after generation, this falls back to the canonical scenario/yaml_file layout.
    """
    raw_path = Path(str(case["yaml_path"]))
    if raw_path.exists():
        return raw_path

    scenario = str(case.get("scenario", ""))
    yaml_file = str(case.get("yaml_file", raw_path.name))
    project_root = Path(__file__).resolve().parents[1]
    candidates = [
        project_root / "models" / "single_product" / scenario / yaml_file,
        Path.cwd() / "models" / "single_product" / scenario / yaml_file,
    ]
    if Path.cwd().name == "notebooks":
        candidates.extend([
            Path.cwd().parent / "models" / "single_product" / scenario / yaml_file,
        ])
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return raw_path


def run_single_yaml_case(
    yaml_path: str | Path,
    scenario_id: str,
    scenario_cfg: Dict[str, Any],
    target_product: str,
    temperature_C: float,
    pressure_Pa: float,
    model_mode: str = "single_product",
    phase_name: str = "aqueous",
    solver: str = "vcs",
    max_steps: int = 100_000,
) -> pd.DataFrame:
    """Run one scenario/YAML/temperature equilibrium case."""
    ct = _load_cantera()
    yaml_path = Path(yaml_path)
    sol = ct.Solution(str(yaml_path), phase_name)
    initial_moles = build_initial_moles(scenario_cfg)
    X0 = normalize_moles_to_X(initial_moles, sol.species_names)
    T_K = float(temperature_C) + 273.15
    sol.TPX = T_K, float(pressure_Pa), X0
    status = "ok"
    error_message = ""
    try:
        sol.equilibrate("TP", solver=solver, max_steps=max_steps)
    except Exception as exc:
        status = "failed"
        error_message = repr(exc)

    return pd.DataFrame(
        _result_rows_for_solution(
            sol=sol,
            yaml_path=yaml_path,
            scenario_id=scenario_id,
            target_product=target_product,
            temperature_C=temperature_C,
            pressure_Pa=pressure_Pa,
            initial_moles=initial_moles,
            X0=X0,
            status=status,
            error_message=error_message,
            model_mode=model_mode,
        )
    )


def make_raw_wide(raw_long_df: pd.DataFrame, output_csv: str | Path) -> pd.DataFrame:
    """Create a wide mole-fraction table from the raw long output."""
    index_cols = [
        "scenario", "model_mode", "yaml_file", "target_product", "T_C", "T_K", "P_Pa",
        "solver_status", "error_message",
    ]
    wide = raw_long_df.pivot_table(index=index_cols, columns="species", values="X_eq", aggfunc="first").reset_index()
    wide.columns.name = None
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(output_csv, index=False)
    return wide


def run_equilibrium_manifest(
    manifest_df: pd.DataFrame,
    scenarios: Dict[str, Any],
    temperatures_C: Sequence[float],
    pressure_Pa: float,
    output_long_csv: str | Path,
    phase_name: str = "aqueous",
    solver: str = "vcs",
    max_steps: int = 100_000,
) -> pd.DataFrame:
    """Run equilibrium simulations from a YAML manifest.

    The manifest must contain scenario, model_mode, yaml_path, and target_product.
    This avoids mixing scenario-specific YAMLs with the wrong scenario.
    """
    required = {"scenario", "model_mode", "yaml_path", "target_product"}
    missing = required - set(manifest_df.columns)
    if missing:
        raise ValueError(f"Manifest is missing required columns: {sorted(missing)}")
    all_rows = []
    for _, case in manifest_df.iterrows():
        scenario_id = str(case["scenario"])
        scenario_cfg = scenarios["scenarios"][scenario_id]
        for temp_c in temperatures_C:
            all_rows.append(
                run_single_yaml_case(
                    yaml_path=resolve_yaml_path_from_manifest_row(case),
                    scenario_id=scenario_id,
                    scenario_cfg=scenario_cfg,
                    target_product=str(case["target_product"]),
                    temperature_C=temp_c,
                    pressure_Pa=pressure_Pa,
                    model_mode=str(case["model_mode"]),
                    phase_name=phase_name,
                    solver=solver,
                    max_steps=max_steps,
                )
            )
    df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    output_long_csv = Path(output_long_csv)
    output_long_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_long_csv, index=False)
    return df
