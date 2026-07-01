"""High-level result summaries for quick workflow inspection."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mole_balance import run_group_columns


def classify_formation(
    x_eq: float,
    n_eq: float,
    solver_status: str,
    x_threshold: float = 1e-12,
    x_significant_threshold: float = 1e-6,
    n_threshold_mol: float = 0.0,
) -> tuple[bool, str]:
    """Return ``(formed_bool, formation_call)`` — the single source of truth for the
    accessibility thresholds.

    - ``"solver_failed"``   — solver status is not ``ok``
    - ``"significant"``     — formed and ``X_eq >= x_significant_threshold``
    - ``"trace"``           — formed and ``x_threshold <= X_eq < x_significant_threshold``
    - ``"below_threshold"`` — solver ok but ``X_eq``/``n_eq`` below threshold

    ``formed_bool`` is ``True`` for both ``significant`` and ``trace``. Both the base
    ``summarize_target_formation`` and ``sensitivity_summary`` call this so their
    classifications stay byte-for-byte identical.
    """
    formed = (
        str(solver_status) == "ok"
        and pd.notna(x_eq)
        and pd.notna(n_eq)
        and float(x_eq) >= float(x_threshold)
        and float(n_eq) >= float(n_threshold_mol)
    )
    if str(solver_status) != "ok":
        call = "solver_failed"
    elif formed and pd.notna(x_eq) and float(x_eq) >= float(x_significant_threshold):
        call = "significant"
    elif formed:
        call = "trace"
    else:
        call = "below_threshold"
    return bool(formed), call


def summarize_target_formation(
    moles_long_df: pd.DataFrame,
    x_threshold: float = 1e-12,
    n_threshold_mol: float = 0.0,
    x_significant_threshold: float = 1e-6,
    output_csv: str | Path | None = None,
) -> pd.DataFrame:
    """Summarize whether each target product appears at equilibrium.

    A target is marked as formed (``formed_bool = True``) when:
    - the solver status is ``ok``;
    - ``X_eq >= x_threshold``; and
    - ``n_eq_mol >= n_threshold_mol``.

    ``formation_call`` uses a two-tier system above the detection threshold:

    - ``"significant"``     — ``X_eq >= x_significant_threshold`` (default 1e-6)
    - ``"trace"``           — ``x_threshold <= X_eq < x_significant_threshold``
    - ``"below_threshold"`` — ``X_eq < x_threshold``
    - ``"solver_failed"``   — solver did not converge

    ``formed_bool`` is ``True`` for both ``"significant"`` and ``"trace"``.

    This is an equilibrium-accessibility diagnostic, not a kinetic/pathway claim.
    """
    required = {"target_product", "species", "X_eq", "n_eq_mol"}
    missing = required - set(moles_long_df.columns)
    if missing:
        raise ValueError(f"moles_long_df is missing required columns: {sorted(missing)}")

    rows = []
    group_cols = run_group_columns(moles_long_df)
    for key, group in moles_long_df.groupby(group_cols, dropna=False):
        key_tuple = key if isinstance(key, tuple) else (key,)
        group_dict = {col: val for col, val in zip(group_cols, key_tuple)}
        target = group_dict["target_product"]
        target_rows = group[group["species"] == target]
        if target_rows.empty:
            rows.append({
                **group_dict,
                "target_present_in_yaml": False,
                "X_eq": np.nan,
                "n_eq_mol": np.nan,
                "formed_bool": False,
                "formation_call": "target_not_present_in_yaml",
                "solver_status": "",
                "element_balance_relative_spread": np.nan,
            })
            continue
        r = target_rows.iloc[0]
        x_eq = float(r["X_eq"]) if pd.notna(r["X_eq"]) else np.nan
        n_eq = float(r["n_eq_mol"]) if pd.notna(r["n_eq_mol"]) else np.nan
        solver_status = str(r.get("solver_status", ""))
        formed, call = classify_formation(
            x_eq, n_eq, solver_status,
            x_threshold=x_threshold,
            x_significant_threshold=x_significant_threshold,
            n_threshold_mol=n_threshold_mol,
        )
        rows.append({
            **group_dict,
            "target_present_in_yaml": True,
            "X_eq": x_eq,
            "log10_X_eq": np.log10(x_eq) if pd.notna(x_eq) and x_eq > 0 else np.nan,
            "n_eq_mol": n_eq,
            "log10_n_eq_mol": np.log10(n_eq) if pd.notna(n_eq) and n_eq > 0 else np.nan,
            "x_threshold": float(x_threshold),
            "n_threshold_mol": float(n_threshold_mol),
            "formed_bool": bool(formed),
            "formation_call": call,
            "solver_status": solver_status,
            "element_balance_relative_spread": r.get("element_balance_relative_spread", np.nan),
        })
    out = pd.DataFrame(rows)
    sort_cols = [c for c in ["scenario", "target_product", "T_C"] if c in out.columns]
    if sort_cols:
        out = out.sort_values(sort_cols)
    if output_csv is not None:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(output_csv, index=False)
    return out


def target_formation_pivot(
    formation_df: pd.DataFrame,
    scenario: str | None = None,
    value_col: str = "log10_X_eq",
) -> pd.DataFrame:
    """Create a target-product by temperature diagnostic pivot."""
    df = formation_df.copy()
    if scenario is not None:
        df = df[df["scenario"] == scenario]
    if df.empty:
        return pd.DataFrame()
    return df.pivot_table(index="target_product", columns="T_C", values=value_col, aggfunc="first").sort_index()
