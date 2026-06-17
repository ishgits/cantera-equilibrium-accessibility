"""NASA9-style Gibbs-energy fitting utilities."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

R_GAS = 8.3144621
COEFF_COLUMNS = [f"a{i}" for i in range(9)]


def nasa9_properties(a: Sequence[float], T: Sequence[float] | float) -> tuple[np.ndarray, np.ndarray]:
    """Compute H(T) and S(T) from NASA9 coefficients.

    Returns H in J/mol and S in J/(mol K). This mirrors the user's previous
    polynomial-fitting notebook.
    """
    a0, a1, a2, a3, a4, a5, a6, a7, a8 = np.asarray(a, dtype=float)
    T = np.asarray(T, dtype=float)
    lnT = np.log(T)
    T_inv = 1.0 / T
    T_inv2 = T_inv * T_inv

    H_RT = (
        -a0 * T_inv2
        + a1 * lnT * T_inv
        + a2
        + (a3 / 2.0) * T
        + (a4 / 3.0) * T**2
        + (a5 / 4.0) * T**3
        + (a6 / 5.0) * T**4
        + a7 * T_inv
    )
    S_R = (
        -0.5 * a0 * T_inv2
        - a1 * T_inv
        + a2 * lnT
        + a3 * T
        + (a4 / 2.0) * T**2
        + (a5 / 3.0) * T**3
        + (a6 / 4.0) * T**4
        + a8
    )
    H = H_RT * R_GAS * T
    S = S_R * R_GAS
    return H, S


def compute_gibbs(a: Sequence[float], T: Sequence[float] | float) -> np.ndarray:
    H, S = nasa9_properties(a, T)
    return H - np.asarray(T, dtype=float) * S


def _residuals(a: Sequence[float], T: np.ndarray, G_obs: np.ndarray) -> np.ndarray:
    return compute_gibbs(a, T) - G_obs


def fit_nasa9_segment(T: np.ndarray, G: np.ndarray, initial_guess: Optional[Sequence[float]] = None) -> np.ndarray:
    """Fit one NASA9 segment to G(T)."""
    try:
        from scipy.optimize import least_squares
    except ImportError as exc:
        raise ImportError("scipy is required for NASA9 fitting. Install scipy before running this step.") from exc

    T = np.asarray(T, dtype=float)
    G = np.asarray(G, dtype=float)
    if T.size < 3:
        raise ValueError("Need at least 3 temperature points for a segment fit; more is strongly recommended.")
    if initial_guess is None:
        initial_guess = np.zeros(9)
        # Give the optimizer a rough free-energy scale through entropy/enthalpy constants.
        initial_guess[7] = np.nanmean(G) / R_GAS
    result = least_squares(_residuals, x0=np.asarray(initial_guess, dtype=float), args=(T, G), max_nfev=200_000)
    if not result.success:
        raise RuntimeError(f"NASA9 fit failed: {result.message}")
    return result.x


def fit_nasa9_two_range(T: np.ndarray, G: np.ndarray, T_split: float = 500.0) -> tuple[np.ndarray, np.ndarray]:
    """Fit low/high NASA9 coefficient sets split at T_split."""
    T = np.asarray(T, dtype=float)
    G = np.asarray(G, dtype=float)
    order = np.argsort(T)
    T, G = T[order], G[order]
    low_mask = T <= T_split
    high_mask = T > T_split
    if low_mask.sum() < 3 or high_mask.sum() < 3:
        raise ValueError(
            f"Need at least 3 points on each side of T_split={T_split}. "
            f"Got low={low_mask.sum()}, high={high_mask.sum()}."
        )
    a_low = fit_nasa9_segment(T[low_mask], G[low_mask])
    # Warm start high with low coefficients.
    a_high = fit_nasa9_segment(T[high_mask], G[high_mask], initial_guess=a_low)
    return a_low, a_high


def _fit_diagnostics(T: np.ndarray, G: np.ndarray, a_low: np.ndarray, a_high: np.ndarray, T_split: float) -> dict:
    pred = np.where(T <= T_split, compute_gibbs(a_low, T), compute_gibbs(a_high, T))
    residual = pred - G
    return {
        "rmse_J_mol": float(np.sqrt(np.mean(residual**2))),
        "max_abs_residual_J_mol": float(np.max(np.abs(residual))),
        "mean_abs_residual_J_mol": float(np.mean(np.abs(residual))),
    }


def fit_all_species(
    gibbs_wide_csv: str | Path,
    coefficients_csv: str | Path,
    diagnostics_csv: str | Path,
    figures_dir: str | Path | None = None,
    T_split: float = 500.0,
    T_low_min: Optional[float] = None,
    T_high_max: Optional[float] = None,
    make_plots: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit NASA9 coefficients for every non-temperature column in a wide G(T) CSV."""
    df = pd.read_csv(gibbs_wide_csv)
    temp_col = None
    for c in df.columns:
        if c.lower().replace(" ", "") in {"t_k", "t(k)", "temperature_k", "temperature(k)"}:
            temp_col = c
            break
    if temp_col is None and "T_K" in df.columns:
        temp_col = "T_K"
    if temp_col is None:
        raise ValueError(f"Could not find a temperature column in {gibbs_wide_csv}")

    T_all = pd.to_numeric(df[temp_col], errors="coerce").to_numpy(float)
    range_mask = np.isfinite(T_all)
    if T_low_min is not None:
        range_mask &= T_all >= T_low_min
    if T_high_max is not None:
        range_mask &= T_all <= T_high_max

    coeff_rows = []
    diag_rows = []
    for species_col in [c for c in df.columns if c != temp_col]:
        G_all = pd.to_numeric(df[species_col], errors="coerce").to_numpy(float)
        mask = range_mask & np.isfinite(G_all)
        T = T_all[mask]
        G = G_all[mask]
        a_low, a_high = fit_nasa9_two_range(T, G, T_split=T_split)
        diag = _fit_diagnostics(T, G, a_low, a_high, T_split=T_split)
        t_min = float(np.min(T))
        t_max = float(np.max(T))
        coeff_rows.append({
            "cantera_name": species_col,
            "range_label": "low",
            "T_low_K": t_min,
            "T_high_K": float(T_split),
            **{f"a{i}": float(v) for i, v in enumerate(a_low)},
        })
        coeff_rows.append({
            "cantera_name": species_col,
            "range_label": "high",
            "T_low_K": float(T_split),
            "T_high_K": t_max,
            **{f"a{i}": float(v) for i, v in enumerate(a_high)},
        })
        diag_rows.append({"cantera_name": species_col, "T_split_K": T_split, **diag})
        if make_plots and figures_dir is not None:
            plot_fit(species_col, T, G, a_low, a_high, T_split, Path(figures_dir))

    coeff_df = pd.DataFrame(coeff_rows)
    diag_df = pd.DataFrame(diag_rows)
    coefficients_csv = Path(coefficients_csv)
    diagnostics_csv = Path(diagnostics_csv)
    coefficients_csv.parent.mkdir(parents=True, exist_ok=True)
    diagnostics_csv.parent.mkdir(parents=True, exist_ok=True)
    coeff_df.to_csv(coefficients_csv, index=False)
    diag_df.to_csv(diagnostics_csv, index=False)
    return coeff_df, diag_df


def plot_fit(species_name: str, T: np.ndarray, G: np.ndarray, a_low: np.ndarray, a_high: np.ndarray, T_split: float, figures_dir: Path) -> None:
    """Save observed-vs-fit and residual plots for one species."""
    import matplotlib.pyplot as plt

    figures_dir.mkdir(parents=True, exist_ok=True)
    order = np.argsort(T)
    T = T[order]
    G = G[order]
    G_fit = np.where(T <= T_split, compute_gibbs(a_low, T), compute_gibbs(a_high, T))
    residual = G_fit - G

    safe = species_name.replace("/", "_").replace("(", "").replace(")", "")

    fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
    ax.plot(T, G, marker="o", linestyle="", label="Observed")
    ax.plot(T, G_fit, linestyle="-", label="NASA9 fit")
    ax.axvline(T_split, linestyle="--", linewidth=1)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Gibbs free energy (J/mol)")
    ax.set_title(species_name)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / f"{safe}_fit.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 3), dpi=200)
    ax.axhline(0, linewidth=1)
    ax.plot(T, residual, marker="o", linestyle="-")
    ax.axvline(T_split, linestyle="--", linewidth=1)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Fit − observed (J/mol)")
    ax.set_title(f"{species_name} residuals")
    fig.tight_layout()
    fig.savefig(figures_dir / f"{safe}_residuals.png")
    plt.close(fig)


def load_coefficients(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"cantera_name", "range_label", "T_low_K", "T_high_K", *COEFF_COLUMNS}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Coefficient table is missing columns: {sorted(missing)}")
    return df


def coefficients_for_species(coeff_df: pd.DataFrame, cantera_name: str) -> tuple[dict, dict]:
    sub = coeff_df[coeff_df["cantera_name"] == cantera_name]
    if sub.empty:
        raise KeyError(f"No NASA9 coefficients found for species {cantera_name!r}")
    low = sub[sub["range_label"].str.lower() == "low"]
    high = sub[sub["range_label"].str.lower() == "high"]
    if low.empty or high.empty:
        raise KeyError(f"Need both low and high coefficient rows for {cantera_name!r}")
    return low.iloc[0].to_dict(), high.iloc[0].to_dict()
