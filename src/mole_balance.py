"""Convert equilibrium mole fractions into moles using elemental conservation."""
from __future__ import annotations

from typing import Dict, Sequence

import numpy as np
import pandas as pd

from formula_tools import parse_formula


DEFAULT_RUN_GROUP_COLS = ("scenario", "model_mode", "yaml_file", "target_product", "T_C")


def run_group_columns(df: pd.DataFrame) -> list[str]:
    """Return the columns that uniquely identify one equilibrium run.

    A run is one (scenario, model_mode, yaml_file, target_product, temperature)
    combination. Only the columns present in ``df`` are returned.
    """
    return [c for c in DEFAULT_RUN_GROUP_COLS if c in df.columns]


def species_compositions(species_df: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    return {row["cantera_name"]: parse_formula(row["formula"]) for _, row in species_df.iterrows()}


def initial_element_moles(initial_moles: Dict[str, float], compositions: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    totals: Dict[str, float] = {}
    for species, n in initial_moles.items():
        if n == 0:
            continue
        if species not in compositions:
            raise KeyError(f"Missing formula/composition for initial species {species!r}")
        for element, count in compositions[species].items():
            totals[element] = totals.get(element, 0.0) + float(n) * float(count)
    return totals


def reconstruct_total_moles_from_elements(
    x_eq: Dict[str, float],
    initial_moles: Dict[str, float],
    compositions: Dict[str, Dict[str, float]],
    min_element_moles: float = 1e-30,
) -> tuple[float, Dict[str, float]]:
    """Infer final total moles from elemental conservation.

    For each conserved element E:
        N_total = initial_element_moles[E] / sum_i(X_i * atoms_E_i)

    Returns the median N_total and per-element estimates. Per-element estimates
    should agree closely. Large disagreement is a diagnostic flag.
    """
    init_elements = initial_element_moles(initial_moles, compositions)
    estimates: Dict[str, float] = {}
    for element, init_amount in init_elements.items():
        if init_amount <= min_element_moles:
            continue
        denom = 0.0
        for species, x in x_eq.items():
            comp = compositions.get(species, {})
            denom += float(x) * float(comp.get(element, 0.0))
        if denom > 0:
            estimates[element] = float(init_amount) / denom
    if not estimates:
        raise ValueError("Could not reconstruct total moles; no positive elemental estimates.")
    values = np.array(list(estimates.values()), dtype=float)
    return float(np.median(values)), estimates


def mole_balance_error(element_estimates: Dict[str, float]) -> float:
    """Return relative spread in element-wise total-mole estimates."""
    vals = np.array(list(element_estimates.values()), dtype=float)
    if vals.size == 0:
        return np.nan
    med = np.median(vals)
    if med == 0:
        return np.nan
    return float((np.max(vals) - np.min(vals)) / abs(med))


def add_equilibrium_moles(
    raw_long_df: pd.DataFrame,
    species_df: pd.DataFrame,
    group_cols: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Add reconstructed equilibrium moles to a raw mole-fraction table."""
    if group_cols is None:
        group_cols = run_group_columns(raw_long_df)
    compositions = species_compositions(species_df)
    rows = []
    for _, group in raw_long_df.groupby(list(group_cols), dropna=False):
        g = group.copy()
        x_eq = dict(zip(g["species"], g["X_eq"]))
        # Initial moles are repeated per species row.
        initial = {}
        for _, r in g.iterrows():
            init = r.get("initial_moles", np.nan)
            if pd.notna(init) and float(init) > 0:
                initial[r["species"]] = float(init)
        n_total, estimates = reconstruct_total_moles_from_elements(x_eq, initial, compositions)
        err = mole_balance_error(estimates)
        for _, r in g.iterrows():
            d = r.to_dict()
            d["n_total_eq_mol"] = n_total
            d["n_eq_mol"] = float(r["X_eq"]) * n_total
            d["element_balance_relative_spread"] = err
            d["element_total_mole_estimates"] = ";".join(f"{k}:{v:.8e}" for k, v in sorted(estimates.items()))
            rows.append(d)
    return pd.DataFrame(rows)
