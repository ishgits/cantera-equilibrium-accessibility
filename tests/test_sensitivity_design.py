"""Tests for the Phase 1 sensitivity design layer."""
import copy
from pathlib import Path

import pandas as pd
import pytest
import yaml

from sensitivity_design import (
    StudyConfigError,
    build_full_design_matrix,
    build_inventory_landscape_design,
    design_matrix_to_scenarios_yaml,
    make_sweep_values,
    make_thermo_variant_name,
    validate_study,
)


# A complete alanine-MVP-like config used across tests.
BASE_CONFIG = {
    "study": {"study_id": "alanine_mvp", "output_dir": None},  # output_dir set per-test
    "mode": {"type": "single_product_sensitivity", "target_products": ["Alanine(aq)"]},
    "species_files": {"species_csv": "inputs/species_example.csv"},
    "base_conditions": {"temperature_C": 0, "pressure_Pa": 101325},
    "model": {"allowed_species": ["H2O(l)", "HCN(aq)", "C2H2(aq)", "NH3(aq)"]},
    "thresholds": {
        "formation_X_threshold": 1.0e-12,
        "significant_X_threshold": 1.0e-6,
        "formation_n_threshold_mol": 0.0,
        "balance_tol": 1.0e-6,
    },
    "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.020},
    "sweeps": {
        "inventory_landscape": {
            "enabled": True,
            "variables": {
                "NH3(aq)": {"type": "linear", "min": 0.0, "max": 0.15, "points": 25},
                "C2H2_over_HCN": {"type": "linear", "min": 0.0, "max": 5.0, "points": 25},
            },
        },
        "deltaG_sweep": {
            "enabled": True,
            "species": "Alanine(aq)",
            "offsets_kJ_mol": {"type": "linear", "min": -50, "max": 50, "points": 21},
            "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.020, "C2H2(aq)": 0.042, "NH3(aq)": 0.05},
        },
        "nh3_deltaG_landscape": {
            "enabled": True,
            "variables": {
                "NH3(aq)": {"type": "linear", "min": 0.0, "max": 0.15, "points": 25},
                "deltaG_offset_kJ_mol": {"type": "linear", "min": -50, "max": 50, "points": 21},
            },
            "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.020, "C2H2(aq)": 0.042},
        },
    },
}


def _config(tmp_path, **overrides):
    cfg = copy.deepcopy(BASE_CONFIG)
    cfg["study"]["output_dir"] = str(tmp_path / "study")
    for key, value in overrides.items():
        cfg[key] = value
    return cfg


SPECIES_DF = pd.DataFrame({
    "cantera_name": ["H2O(l)", "HCN(aq)", "C2H2(aq)", "NH3(aq)", "Alanine(aq)"],
})


# --------------------------------------------------------------------------- #
# Sweep values
# --------------------------------------------------------------------------- #
def test_linear_sweep_count_and_endpoints():
    vals = make_sweep_values({"type": "linear", "min": 0.0, "max": 0.15, "points": 25})
    assert len(vals) == 25
    assert vals[0] == pytest.approx(0.0)
    assert vals[-1] == pytest.approx(0.15)


def test_explicit_sweep_preserves_values():
    vals = make_sweep_values({"type": "explicit", "values": [-50, -20, 0, 20, 50]})
    assert vals == [-50.0, -20.0, 0.0, 20.0, 50.0]


def test_logspace_uses_actual_values():
    vals = make_sweep_values({"type": "logspace", "min": 1e-3, "max": 1.0, "points": 4})
    assert vals[0] == pytest.approx(1e-3)
    assert vals[-1] == pytest.approx(1.0)


def test_logspace_rejects_nonpositive_min():
    with pytest.raises(StudyConfigError):
        make_sweep_values({"type": "logspace", "min": 0.0, "max": 1.0, "points": 4})


def test_zero_points_rejected():
    with pytest.raises(StudyConfigError):
        make_sweep_values({"type": "linear", "min": 0.0, "max": 1.0, "points": 0})


def test_include_values_merges_sorted_deduped():
    vals = make_sweep_values({"type": "linear", "min": 0.0, "max": 0.15, "points": 4,
                              "include_values": [0.02, 2.1]})
    assert 0.02 in vals and 2.1 in vals
    assert vals == sorted(vals) and len(vals) == len(set(vals))
    # explicit axis + include_values
    v2 = make_sweep_values({"type": "explicit", "values": [0.01, 0.05],
                            "include_values": [0.10, 0.05]})  # 0.05 duplicate dropped
    assert set(v2) == {0.01, 0.05, 0.10}


# --------------------------------------------------------------------------- #
# Variant naming
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("offset,expected", [
    (-50, "Alanine__dG_m050(aq)"),
    (-5, "Alanine__dG_m005(aq)"),
    (0, "Alanine__dG_000(aq)"),
    (20, "Alanine__dG_p020(aq)"),
    (50, "Alanine__dG_p050(aq)"),
])
def test_variant_naming(offset, expected):
    assert make_thermo_variant_name("Alanine(aq)", offset) == expected


def test_variant_naming_rejects_noninteger_offset():
    with pytest.raises(ValueError):
        make_thermo_variant_name("Alanine(aq)", 20.5)


# --------------------------------------------------------------------------- #
# Design matrix structure & counts
# --------------------------------------------------------------------------- #
def test_inventory_only_count(tmp_path):
    cfg = _config(tmp_path)
    cfg["sweeps"]["deltaG_sweep"]["enabled"] = False
    cfg["sweeps"]["nh3_deltaG_landscape"]["enabled"] = False
    matrix = build_full_design_matrix(cfg)
    assert len(matrix) == 625


def test_all_substudies_count(tmp_path):
    matrix = build_full_design_matrix(_config(tmp_path))
    assert len(matrix) == 1171  # 625 + 21 + 525


def test_case_ids_unique(tmp_path):
    matrix = build_full_design_matrix(_config(tmp_path))
    assert not matrix["case_id"].duplicated().any()
    assert (matrix["scenario_id"] == matrix["case_id"]).all()


def test_no_negative_moles(tmp_path):
    matrix = build_full_design_matrix(_config(tmp_path))
    mol_cols = [c for c in matrix.columns if c.endswith("_mol") and not c.endswith("kJ_mol")]
    assert (matrix[mol_cols] >= 0).all().all()


def test_c2h2_equals_ratio_times_hcn():
    df = build_inventory_landscape_design(BASE_CONFIG)
    expected = df["C2H2_over_HCN"] * 0.020
    pd.testing.assert_series_equal(df["C2H2_mol"], expected, check_names=False)


def test_deltaG_substudy_uses_variant_target(tmp_path):
    matrix = build_full_design_matrix(_config(tmp_path))
    dg = matrix[matrix["substudy_id"] == "deltaG_sweep"]
    # Variant name encodes the offset; base target unchanged.
    assert (dg["target_product"] == "Alanine(aq)").all()
    assert dg["target_variant"].str.startswith("Alanine__dG_").all()
    inv = matrix[matrix["substudy_id"] == "inventory_landscape"]
    assert (inv["target_variant"] == "Alanine(aq)").all()


# --------------------------------------------------------------------------- #
# Scenario YAML emission
# --------------------------------------------------------------------------- #
def test_scenarios_carry_allowed_species_and_zeros(tmp_path):
    cfg = _config(tmp_path)
    matrix = build_full_design_matrix(cfg)
    out = tmp_path / "generated_scenarios.yaml"
    scenarios = design_matrix_to_scenarios_yaml(matrix, cfg, out)

    # Round-trips as valid YAML.
    loaded = yaml.safe_load(out.read_text())["scenarios"]
    assert len(loaded) == len(matrix)

    # The first inventory case has NH3 = 0 and C2H2 = 0 retained (present-but-zero).
    first = scenarios[matrix.iloc[0]["scenario_id"]]
    assert first["initial_moles"]["NH3(aq)"] == 0.0
    assert first["initial_moles"]["C2H2(aq)"] == 0.0
    # Explicit allowed_species = reactants + the case's target variant.
    assert "H2O(l)" in first["allowed_species"]
    assert first["target_products"][0] in first["allowed_species"]

    # A ΔG scenario allows its variant species, not the bare target.
    dg_id = matrix[matrix["substudy_id"] == "deltaG_sweep"].iloc[0]["scenario_id"]
    dg_scenario = scenarios[dg_id]
    assert dg_scenario["target_products"][0].startswith("Alanine__dG_")
    assert dg_scenario["target_products"][0] in dg_scenario["allowed_species"]
    assert "Alanine(aq)" not in dg_scenario["allowed_species"]


# --------------------------------------------------------------------------- #
# Validator (plain-English failures)
# --------------------------------------------------------------------------- #
def test_validator_passes_on_good_config(tmp_path):
    cfg = _config(tmp_path)
    seed = pd.DataFrame({"T_K": [273.16, 643.15]})
    validate_study(cfg, SPECIES_DF, seed)  # should not raise


def test_validator_flags_misspelled_species(tmp_path):
    cfg = _config(tmp_path)
    # A fixed-inventory species not in the CSV is named in a plain-English error.
    cfg["fixed_inventory"] = {"H2O(l)": 1.0, "HCNN(aq)": 0.02}
    with pytest.raises(StudyConfigError) as exc:
        validate_study(cfg, SPECIES_DF)
    assert "HCNN(aq)" in str(exc.value)


def test_validator_flags_temperature_out_of_range(tmp_path):
    cfg = _config(tmp_path)
    cfg["base_conditions"]["temperature_C"] = 1000
    seed = pd.DataFrame({"T_K": [273.16, 643.15]})
    with pytest.raises(StudyConfigError):
        validate_study(cfg, SPECIES_DF, seed)


def test_validator_flags_nonpositive_threshold(tmp_path):
    cfg = _config(tmp_path)
    cfg["thresholds"]["formation_X_threshold"] = 0
    with pytest.raises(StudyConfigError):
        validate_study(cfg, SPECIES_DF)


def test_validator_accepts_integer_offsets(tmp_path):
    # Default config: -200..200 over 11 points = step 40 (all integers).
    validate_study(_config(tmp_path), SPECIES_DF)  # should not raise


def test_validator_flags_noninteger_deltaG_sweep_offsets(tmp_path):
    cfg = _config(tmp_path)
    cfg["sweeps"]["deltaG_sweep"]["offsets_kJ_mol"] = {
        "type": "explicit", "values": [-40, 0, 20.5]}
    with pytest.raises(StudyConfigError, match="whole number"):
        validate_study(cfg, SPECIES_DF)


def test_validator_flags_noninteger_nh3_deltaG_axis(tmp_path):
    cfg = _config(tmp_path)
    # linspace(-50, 50, 4) -> ±16.667, ±50 (non-integer).
    cfg["sweeps"]["nh3_deltaG_landscape"]["variables"]["deltaG_offset_kJ_mol"] = {
        "type": "linear", "min": -50, "max": 50, "points": 4}
    with pytest.raises(StudyConfigError, match="whole number"):
        validate_study(cfg, SPECIES_DF)


def test_load_gibbs_seed_raises_on_set_but_missing(tmp_path):
    from sensitivity_design import load_gibbs_seed_for_config
    cfg = _config(tmp_path)
    cfg["species_files"] = {"species_csv": "inputs/species_example.csv",
                            "gibbs_seed_wide_csv": "data/processed/does_not_exist.csv"}
    with pytest.raises(StudyConfigError, match="Gibbs seed file not found"):
        load_gibbs_seed_for_config(cfg, repo_root=tmp_path)


def test_load_gibbs_seed_returns_none_when_key_absent(tmp_path):
    from sensitivity_design import load_gibbs_seed_for_config
    cfg = _config(tmp_path)
    cfg["species_files"] = {"species_csv": "inputs/species_example.csv"}
    assert load_gibbs_seed_for_config(cfg, repo_root=tmp_path) is None


def test_validator_rejects_multi_target(tmp_path):
    cfg = _config(tmp_path)
    cfg["mode"]["target_products"] = ["Alanine(aq)", "Adenine(aq)"]
    with pytest.raises(StudyConfigError, match="exactly one target"):
        validate_study(cfg, SPECIES_DF)


def test_validator_rejects_unknown_substudy(tmp_path):
    cfg = _config(tmp_path)
    cfg["sweeps"]["bogus_sweep"] = {"enabled": True, "variables": {}}
    with pytest.raises(StudyConfigError, match="Unknown substudy"):
        validate_study(cfg, SPECIES_DF)


def test_validator_rejects_missing_required_axis(tmp_path):
    cfg = _config(tmp_path)
    cfg["sweeps"]["deltaG_sweep"].pop("offsets_kJ_mol")     # required axis for deltaG_sweep
    with pytest.raises(StudyConfigError, match="required axis"):
        validate_study(cfg, SPECIES_DF)


def test_validator_rejects_swept_species_outside_allowed(tmp_path):
    cfg = _config(tmp_path)
    # Sweep a real species that is not in the phase (allowed_species).
    cfg["model"]["allowed_species"] = ["H2O(l)", "HCN(aq)", "C2H2(aq)"]   # NH3 excluded
    with pytest.raises(StudyConfigError, match="not in model.allowed_species"):
        validate_study(cfg, SPECIES_DF)


def test_comment_aware_gibbs_seed_loads(tmp_path):
    # The bundled example seed carries a '#' comment header; it must load + fit.
    from sensitivity_design import load_gibbs_seed_for_config
    from thermo_fit import read_wide_gibbs_csv
    repo_root = Path(__file__).resolve().parents[1]
    df = read_wide_gibbs_csv(repo_root / "inputs" / "example_validation_gibbs.csv")
    assert not df.empty and df.shape[1] > 1
    cfg = _config(tmp_path)
    cfg["species_files"] = {"species_csv": "inputs/species_example.csv",
                            "gibbs_seed_wide_csv": "inputs/example_validation_gibbs.csv"}
    seed = load_gibbs_seed_for_config(cfg, repo_root=repo_root)
    assert seed is not None and not seed.empty
