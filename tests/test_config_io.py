"""Tests for species/scenario loading and validation."""
import textwrap

import pytest

from config_io import load_species_metadata, load_scenarios, list_target_products


VALID_SPECIES = (
    "species_key,cantera_name,chnosz_name,formula,state,molar_volume_cm3_mol,role,notes,product_class\n"
    "h2o,H2O(l),water,H2O,liq,18.015,solvent,,\n"
    "benzene,Benzene(aq),benzene,C6H6,aq,89.13,reactant,,\n"
    "uracil,Uracil(aq),uracil,C4H4N2O2,aq,112.09,product,,nucleobase\n"
)


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def test_load_valid_species(tmp_path):
    df = load_species_metadata(_write(tmp_path, "sp.csv", VALID_SPECIES))
    assert len(df) == 3
    assert list_target_products(df) == ["Uracil(aq)"]


def test_missing_required_column_raises(tmp_path):
    bad = VALID_SPECIES.replace("molar_volume_cm3_mol,", "")
    bad = bad.replace("18.015,", "").replace("89.13,", "").replace("112.09,", "")
    with pytest.raises(ValueError):
        load_species_metadata(_write(tmp_path, "bad.csv", bad))


def test_duplicate_species_key_raises(tmp_path):
    dup = VALID_SPECIES + "h2o,Water2(aq),water,H2O,aq,18.0,solvent,,\n"
    with pytest.raises(ValueError):
        load_species_metadata(_write(tmp_path, "dup.csv", dup))


def test_non_numeric_molar_volume_raises(tmp_path):
    bad = VALID_SPECIES.replace("18.015", "not_a_number")
    with pytest.raises(ValueError):
        load_species_metadata(_write(tmp_path, "badvol.csv", bad))


def test_load_valid_scenarios(tmp_path):
    yaml_text = textwrap.dedent(
        """
        scenarios:
          cond1:
            description: "test"
            initial_moles:
              H2O(l): 1.0
              Benzene(aq): 0.001
        """
    )
    data = load_scenarios(_write(tmp_path, "sc.yaml", yaml_text))
    assert "cond1" in data["scenarios"]


def test_scenario_missing_initial_moles_raises(tmp_path):
    yaml_text = "scenarios:\n  cond1:\n    description: nope\n"
    with pytest.raises(ValueError):
        load_scenarios(_write(tmp_path, "bad.yaml", yaml_text))


def test_scenario_non_numeric_moles_raises(tmp_path):
    yaml_text = "scenarios:\n  cond1:\n    initial_moles:\n      H2O(l): a_lot\n"
    with pytest.raises(ValueError):
        load_scenarios(_write(tmp_path, "bad2.yaml", yaml_text))


def test_missing_top_level_scenarios_raises(tmp_path):
    with pytest.raises(ValueError):
        load_scenarios(_write(tmp_path, "empty.yaml", "something_else: 1\n"))
