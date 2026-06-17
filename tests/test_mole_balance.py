"""Tests for elemental-conservation mole reconstruction."""
import pytest

from mole_balance import (
    initial_element_moles,
    reconstruct_total_moles_from_elements,
    mole_balance_error,
)


def test_initial_element_moles():
    comps = {"H2O": {"H": 2, "O": 1}}
    totals = initial_element_moles({"H2O": 1.5}, comps)
    assert totals == {"H": 3.0, "O": 1.5}


def test_single_element_conservation_exact():
    # Two interconverting 1-carbon species; carbon is conserved.
    comps = {"R": {"C": 1}, "P": {"C": 1}}
    initial = {"R": 2.0}
    x_eq = {"R": 0.5, "P": 0.5}
    n_total, estimates = reconstruct_total_moles_from_elements(x_eq, initial, comps)
    assert n_total == pytest.approx(2.0)
    assert estimates["C"] == pytest.approx(2.0)
    assert mole_balance_error(estimates) == pytest.approx(0.0)


def test_per_element_estimates_are_deterministic():
    comps = {"H2O": {"H": 2, "O": 1}, "H2": {"H": 2}, "O2": {"O": 2}}
    initial = {"H2O": 1.0}  # H=2, O=1
    x_eq = {"H2O": 0.5, "H2": 0.25, "O2": 0.25}
    n_total, estimates = reconstruct_total_moles_from_elements(x_eq, initial, comps)
    # H: 2 / (0.5*2 + 0.25*2) = 2/1.5 ; O: 1 / (0.5*1 + 0.25*2) = 1/1.0
    assert estimates["H"] == pytest.approx(2.0 / 1.5)
    assert estimates["O"] == pytest.approx(1.0)
    assert n_total == pytest.approx((2.0 / 1.5 + 1.0) / 2.0)  # median of two values
    assert mole_balance_error(estimates) > 0.0


def test_missing_composition_raises():
    with pytest.raises(KeyError):
        initial_element_moles({"Unknown": 1.0}, {})
