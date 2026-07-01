"""Tests for Phase 4 ΔG support: the exact analytic a7-shift (no refit)."""
from pathlib import Path

import pandas as pd
import pytest

from thermo_fit import COEFF_COLUMNS, compute_gibbs, coefficients_for_species
from sensitivity_thermo import (
    augment_species_metadata_with_variants, shift_coeffs_by_gibbs,
)
from sensitivity_design import build_thermo_offsets_table

PROJECT_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def coeff_df():
    return pd.read_csv(PROJECT_ROOT / "tests" / "fixtures" / "nasa9_coefficients.csv")


@pytest.fixture(scope="module")
def species_df():
    return pd.read_csv(PROJECT_ROOT / "inputs" / "species_example.csv")


def _gibbs(coeffs, name, T):
    """G(T) for a species, picking the NASA9 segment that covers T."""
    low, high = coefficients_for_species(coeffs, name)
    seg = low if T <= float(low["T_high_K"]) else high
    a = [float(seg[c]) for c in COEFF_COLUMNS]
    return float(compute_gibbs(a, T))


# Sampled across both segments, incl. 0 C (273.15 K) and either side of T_split=500.
@pytest.mark.parametrize("T", [273.15, 300.0, 480.0, 520.0, 640.0])
@pytest.mark.parametrize("dG", [-200, -40, 0, 40, 200])
def test_gibbs_shift_is_exact(coeff_df, dG, T):
    variant = shift_coeffs_by_gibbs(coeff_df, "Alanine(aq)", "Alanine__dG_x(aq)", dG)
    combined = pd.concat([coeff_df, variant], ignore_index=True)
    delta = _gibbs(combined, "Alanine__dG_x(aq)", T) - _gibbs(combined, "Alanine(aq)", T)
    assert delta == pytest.approx(dG * 1000.0, abs=1e-3)   # J/mol, to floating point


def test_base_coeffs_unmutated(coeff_df):
    before = coeff_df.loc[coeff_df["cantera_name"] == "Alanine(aq)", "a7"].tolist()
    shift_coeffs_by_gibbs(coeff_df, "Alanine(aq)", "Alanine__dG_p040(aq)", 40)
    after = coeff_df.loc[coeff_df["cantera_name"] == "Alanine(aq)", "a7"].tolist()
    assert before == after


def test_shift_returns_two_rows_with_variant_name(coeff_df):
    v = shift_coeffs_by_gibbs(coeff_df, "Alanine(aq)", "Alanine__dG_m040(aq)", -40)
    assert len(v) == 2                                       # low + high segments
    assert set(v["cantera_name"]) == {"Alanine__dG_m040(aq)"}
    assert set(v["range_label"].str.lower()) == {"low", "high"}


def test_shift_missing_base_raises(coeff_df):
    with pytest.raises(KeyError):
        shift_coeffs_by_gibbs(coeff_df, "NotASpecies(aq)", "V(aq)", 10)


def test_variant_metadata_preserves_properties(species_df):
    config = {
        "study": {"study_id": "t"}, "mode": {"target_products": ["Alanine(aq)"]},
        "sweeps": {"deltaG_sweep": {"enabled": True, "species": "Alanine(aq)",
                   "offsets_kJ_mol": {"type": "explicit", "values": [-40, 0, 40]},
                   "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.02}}},
    }
    offsets = build_thermo_offsets_table(config)
    aug = augment_species_metadata_with_variants(species_df, offsets)
    base = species_df[species_df["cantera_name"] == "Alanine(aq)"].iloc[0]
    for variant in offsets["variant_species"]:
        r = aug[aug["cantera_name"] == variant].iloc[0]
        assert r["formula"] == base["formula"]
        assert r["state"] == base["state"]
        assert r["molar_volume_cm3_mol"] == base["molar_volume_cm3_mol"]
        assert r["product_class"] == base["product_class"]
