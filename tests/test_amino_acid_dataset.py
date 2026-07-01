"""Tests for the shared amino-acid dataset (offline: no Cantera, no pyCHNOSZ)."""
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from config_io import load_species_metadata
from formula_tools import parse_formula
import validate_amino_acid_dataset as vad

SPECIES_CSV = PROJECT_ROOT / "inputs" / "amino_acids_species.csv"
SEED_CSV = PROJECT_ROOT / "inputs" / "amino_acids_gibbs_seed.csv"

# The 18 C/H/N/O amino acids and their expected formulas (per the Phase-6 scope).
EXPECTED_AA = {
    "Glycine(aq)": "C2H5NO2", "Alanine(aq)": "C3H7NO2", "Serine(aq)": "C3H7NO3",
    "Proline(aq)": "C5H9NO2", "Valine(aq)": "C5H11NO2", "Threonine(aq)": "C4H9NO3",
    "Leucine(aq)": "C6H13NO2", "Isoleucine(aq)": "C6H13NO2", "Asparagine(aq)": "C4H8N2O3",
    "AsparticAcid(aq)": "C4H7NO4", "Glutamine(aq)": "C5H10N2O3", "Lysine(aq)": "C6H14N2O2",
    "GlutamicAcid(aq)": "C5H9NO4", "Arginine(aq)": "C6H14N4O2", "Histidine(aq)": "C6H9N3O2",
    "Phenylalanine(aq)": "C9H11NO2", "Tyrosine(aq)": "C9H11NO3", "Tryptophan(aq)": "C11H12N2O2",
}
REACTANTS = {"H2O(l)", "HCN(aq)", "C2H2(aq)", "NH3(aq)"}


@pytest.fixture(scope="module")
def species_df():
    return load_species_metadata(SPECIES_CSV)   # raises on dup keys / bad schema


@pytest.fixture(scope="module")
def seed_df():
    return pd.read_csv(SEED_CSV)


def test_dataset_loads_and_counts(species_df):
    assert len(species_df) == 22
    aas = species_df[species_df["product_class"] == "amino_acid"]
    assert len(aas) == 18
    assert REACTANTS <= set(species_df["cantera_name"])


def test_amino_acid_formulas_are_chno_only(species_df):
    aas = species_df[species_df["product_class"] == "amino_acid"]
    for _, r in aas.iterrows():
        elements = set(parse_formula(r["formula"]).keys())
        assert elements <= {"C", "H", "N", "O"}, f"{r['cantera_name']} has {elements}"
        assert "S" not in elements and "Se" not in elements


def test_formulas_match_expected(species_df):
    by_name = dict(zip(species_df["cantera_name"], species_df["formula"]))
    for name, formula in EXPECTED_AA.items():
        assert name in by_name, f"{name} missing from species CSV"
        assert parse_formula(by_name[name]) == parse_formula(formula), name


def test_all_species_present_in_both_files(species_df, seed_df):
    seed_cols = set(seed_df.columns)
    for name in species_df["cantera_name"]:
        assert name in seed_cols, f"{name} missing a G(T) column in the seed"
    assert len(seed_df) == 38            # 0..370 C in 10 C steps
    assert not seed_df.isna().any().any()  # no gaps


def test_molar_volumes_finite_and_positive(species_df):
    v = species_df["molar_volume_cm3_mol"]
    assert v.notna().all() and (v > 0).all()


def test_validation_report_and_nasa9_fit():
    report = vad.run_validation(SPECIES_CSV, SEED_CSV)
    assert report["ok"], report["errors"]
    assert len(report["fit_rmse_J_mol"]) == 22         # every species fit
    assert report["max_rmse_J_mol"] is not None and report["max_rmse_J_mol"] < 1e4
    assert report["fit_outliers"] == []                # no >10x-median outliers
