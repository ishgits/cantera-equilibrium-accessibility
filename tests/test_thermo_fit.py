"""Tests for NASA9 fitting: a smooth G(T) curve must be recovered to small residual."""
import numpy as np

from thermo_fit import compute_gibbs, fit_nasa9_two_range


def _smooth_gibbs(T):
    # A smooth, physically-plausible G(T) (J/mol), curved but well-behaved.
    T = np.asarray(T, dtype=float)
    return -560_000.0 + 120.0 * (T - 300.0) - 0.05 * (T - 300.0) ** 2


def test_two_range_fit_recovers_curve():
    T = np.arange(273.0, 644.0, 10.0)
    G = _smooth_gibbs(T)
    a_low, a_high = fit_nasa9_two_range(T, G, T_split=500.0)
    pred = np.where(T <= 500.0, compute_gibbs(a_low, T), compute_gibbs(a_high, T))
    max_abs = np.max(np.abs(pred - G))
    # Tiny relative to |G| ~ 5.6e5 J/mol.
    assert max_abs < 50.0


def test_fit_requires_enough_points_per_segment():
    T = np.array([300.0, 310.0, 320.0, 520.0, 530.0])  # only 2 points above split
    G = _smooth_gibbs(T)
    try:
        fit_nasa9_two_range(T, G, T_split=500.0)
    except ValueError:
        return
    raise AssertionError("expected ValueError for too few high-segment points")
