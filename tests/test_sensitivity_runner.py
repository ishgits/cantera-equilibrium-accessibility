"""Tests for the Phase 2 sensitivity runner (manifests, model reuse, YAML, provenance).

These exercise everything that does not require Cantera; the actual equilibrate()
call is covered by the base engine and is not re-tested here.
"""
import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from sensitivity_design import (
    build_full_design_matrix,
    build_inventory_landscape_design,
    build_thermo_offsets_table,
    write_design_outputs,
)
from sensitivity_design import StudyConfigError
from sensitivity_thermo import (
    augment_species_metadata_with_variants, require_gibbs_seed, shift_coeffs_by_gibbs,
)
from sensitivity_runner import (
    build_model_manifest,
    build_run_manifest,
    compute_model_table,
    merge_raw_results_with_design,
    model_identity,
    model_reuse_stats,
    nasa9_coeff_hash,
    resolve_model_yaml_path,
    run_sensitivity_manifest,
    seed_is_stale,
    select_cases_to_run,
    thermo_variant_id_for,
    write_run_provenance,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ALLOWED_BASE = ["H2O(l)", "HCN(aq)", "C2H2(aq)", "NH3(aq)"]


@pytest.fixture(scope="module")
def coeff_df():
    return pd.read_csv(PROJECT_ROOT / "tests" / "fixtures" / "nasa9_coefficients.csv")


@pytest.fixture(scope="module")
def species_df():
    return pd.read_csv(PROJECT_ROOT / "inputs" / "species_example.csv")


def _small_inventory_config(points=3):
    return {
        "study": {"study_id": "test_study"},
        "mode": {"target_products": ["Alanine(aq)"]},
        "base_conditions": {"temperature_C": 0, "pressure_Pa": 101325},
        "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.020},
        "sweeps": {
            "inventory_landscape": {
                "enabled": True,
                "variables": {
                    "NH3(aq)": {"type": "linear", "min": 0.0, "max": 0.15, "points": points},
                    "C2H2_over_HCN": {"type": "linear", "min": 0.0, "max": 5.0, "points": points},
                },
            },
        },
    }


# --------------------------------------------------------------------------- #
# Model identity
# --------------------------------------------------------------------------- #
def test_thermo_variant_id():
    assert thermo_variant_id_for("Alanine(aq)") == "base"
    assert thermo_variant_id_for("Alanine__dG_m040(aq)") == "dG_m040"


def test_coeff_hash_changes_with_a7(coeff_df):
    base_hash = nasa9_coeff_hash(coeff_df, ["Alanine(aq)"])
    shifted = coeff_df.copy()
    mask = shifted["cantera_name"] == "Alanine(aq)"
    shifted.loc[mask, "a7"] = shifted.loc[mask, "a7"] + 1000.0
    assert nasa9_coeff_hash(shifted, ["Alanine(aq)"]) != base_hash


def test_model_identity_distinguishes_variants():
    base = model_identity(["A", "B", "X(aq)"], "base", "X(aq)", "aqueous", "h1")
    variant = model_identity(["A", "B", "X__dG_p020(aq)"], "dG_p020",
                             "X__dG_p020(aq)", "aqueous", "h2")
    same = model_identity(["B", "A", "X(aq)"], "base", "X(aq)", "aqueous", "h1")
    assert base != variant
    assert base == same  # order-independent over allowed_species


# --------------------------------------------------------------------------- #
# Model reuse / manifests
# --------------------------------------------------------------------------- #
def test_inventory_collapses_to_one_model(coeff_df):
    design = build_inventory_landscape_design(_small_inventory_config(points=5))
    design["case_id"] = [f"C{i}" for i in range(len(design))]
    design["scenario_id"] = design["case_id"]
    design_mid, models_df = compute_model_table(design, ALLOWED_BASE, coeff_df)
    stats = model_reuse_stats(design_mid, models_df, expected_models=1)
    assert stats["n_models"] == 1
    assert stats["n_cases"] == 25
    assert stats["reuse_ratio"] == 25
    assert design_mid["model_id"].nunique() == 1


def test_distinct_variants_make_distinct_models(coeff_df):
    # Synthesize a ΔG variant coefficient set so two target variants coexist.
    variant_rows = coeff_df[coeff_df["cantera_name"] == "Alanine(aq)"].copy()
    variant_rows["cantera_name"] = "Alanine__dG_p020(aq)"
    variant_rows["a7"] = variant_rows["a7"] + 2406.0
    coeffs = pd.concat([coeff_df, variant_rows], ignore_index=True)
    design = pd.DataFrame({
        "case_id": ["c1", "c2"],
        "scenario_id": ["c1", "c2"],
        "target_product": ["Alanine(aq)", "Alanine(aq)"],
        "target_variant": ["Alanine(aq)", "Alanine__dG_p020(aq)"],
        "T_C": [0.0, 0.0], "P_Pa": [101325.0, 101325.0],
    })
    _, models_df = compute_model_table(design, ALLOWED_BASE, coeffs)
    assert len(models_df) == 2


def _full_alanine_config():
    """The real alanine MVP grid dimensions (911 cases, 11 ΔG offsets)."""
    return {
        "study": {"study_id": "alanine_mvp"},
        "mode": {"type": "single_product_sensitivity", "target_products": ["Alanine(aq)"]},
        "base_conditions": {"temperature_C": 0, "pressure_Pa": 101325},
        "model": {"allowed_species": ALLOWED_BASE},
        "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.02},
        "sweeps": {
            "inventory_landscape": {"enabled": True, "variables": {
                "NH3(aq)": {"type": "linear", "min": 0, "max": 0.15, "points": 25},
                "C2H2_over_HCN": {"type": "linear", "min": 0, "max": 5, "points": 25}}},
            "deltaG_sweep": {"enabled": True, "species": "Alanine(aq)",
                "offsets_kJ_mol": {"type": "linear", "min": -200, "max": 200, "points": 11},
                "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.02,
                                    "C2H2(aq)": 0.042, "NH3(aq)": 0.05}},
            "nh3_deltaG_landscape": {"enabled": True, "variables": {
                "NH3(aq)": {"type": "linear", "min": 0, "max": 0.15, "points": 25},
                "deltaG_offset_kJ_mol": {"type": "linear", "min": -200, "max": 200, "points": 11}},
                "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.02, "C2H2(aq)": 0.042}},
        },
    }


def test_model_reuse_911_to_12(coeff_df):
    cfg = _full_alanine_config()
    full = build_full_design_matrix(cfg)
    assert len(full) == 911

    # Build base + 11 analytic ΔG-variant coefficient sets.
    offsets = build_thermo_offsets_table(cfg)
    variants = [shift_coeffs_by_gibbs(coeff_df, o.base_species, o.variant_species,
                                      o.deltaG_offset_kJ_mol)
                for o in offsets.itertuples()]
    coeffs = pd.concat([coeff_df, *variants], ignore_index=True)

    design_mid, models = compute_model_table(full, ALLOWED_BASE, coeffs)
    assert len(models) == 12                                   # 1 base + 11 ΔG variants
    inv = design_mid[design_mid["substudy_id"] == "inventory_landscape"]
    assert inv["model_id"].nunique() == 1                      # inventory reuses one model
    dg_cases = design_mid[design_mid["substudy_id"].isin(
        ["deltaG_sweep", "nh3_deltaG_landscape"])]
    assert dg_cases["model_id"].nunique() == 11                # 11 shared variant models

    # offset-0 variant is a DISTINCT model from the base (different species name).
    base_model = inv["model_id"].iloc[0]
    v0 = design_mid[design_mid["target_variant"] == "Alanine__dG_000(aq)"]
    assert v0["model_id"].iloc[0] != base_model

    # The same variant is reused across deltaG_sweep and nh3_deltaG_landscape.
    m_dg = design_mid[(design_mid["substudy_id"] == "deltaG_sweep") &
                      (design_mid["target_variant"] == "Alanine__dG_p040(aq)")]["model_id"].iloc[0]
    m_nh = design_mid[(design_mid["substudy_id"] == "nh3_deltaG_landscape") &
                      (design_mid["target_variant"] == "Alanine__dG_p040(aq)")]["model_id"].iloc[0]
    assert m_dg == m_nh


def test_reuse_assertion_raises_on_mismatch(coeff_df):
    design = build_inventory_landscape_design(_small_inventory_config(points=3))
    design["case_id"] = [f"C{i}" for i in range(len(design))]
    design["scenario_id"] = design["case_id"]
    design_mid, models_df = compute_model_table(design, ALLOWED_BASE, coeff_df)
    with pytest.raises(AssertionError):
        model_reuse_stats(design_mid, models_df, expected_models=99)


def test_full_manifest_and_yaml(tmp_path, species_df, coeff_df):
    design = build_inventory_landscape_design(_small_inventory_config(points=4))
    design["case_id"] = [f"C{i}" for i in range(len(design))]
    design["scenario_id"] = design["case_id"]
    design["study_id"] = "test_study"
    design["substudy_id"] = "inventory_landscape"

    design_mid, models_df, stats = build_model_manifest(
        design, species_df, coeff_df, ALLOWED_BASE, tmp_path, expected_models=1)
    assert stats["n_models"] == 1

    # model_manifest.csv written with the contract columns.
    mm = pd.read_csv(tmp_path / "model_manifest.csv")
    assert {"model_id", "yaml_path", "thermo_variant_id", "species_set_hash",
            "thermo_hash"} <= set(mm.columns)
    assert (mm["thermo_variant_id"] == "base").all()

    # One YAML written; phase species == allowed_base + target.
    yaml_path = tmp_path / mm.iloc[0]["yaml_path"]
    assert yaml_path.exists()
    doc = yaml.safe_load(yaml_path.read_text())
    assert doc["phases"][0]["species"] == ALLOWED_BASE + ["Alanine(aq)"]

    # Run manifest: exactly one row per design case, each with a valid model_id.
    run_manifest = build_run_manifest(design_mid, models_df, tmp_path)
    assert len(run_manifest) == len(design)
    assert run_manifest["case_id"].is_unique
    assert set(run_manifest["model_id"]) <= set(models_df["model_id"])
    assert (run_manifest["status"] == "pending").all()


def test_missing_variant_coeffs_raises_clear_error(tmp_path, species_df, coeff_df):
    # As in the real pipeline, the variant has metadata (augmented) but no coeffs
    # yet (Phase 4 builds the a7-shifted coeffs) — the coeff check should fire.
    config = {
        "study": {"study_id": "t"}, "mode": {"target_products": ["Alanine(aq)"]},
        "sweeps": {"deltaG_sweep": {"enabled": True, "species": "Alanine(aq)",
                   "offsets_kJ_mol": {"type": "explicit", "values": [-40]},
                   "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.02}}},
    }
    species_meta = augment_species_metadata_with_variants(
        species_df, build_thermo_offsets_table(config))
    design = pd.DataFrame({
        "case_id": ["c1"], "scenario_id": ["c1"], "target_product": ["Alanine(aq)"],
        "target_variant": ["Alanine__dG_m040(aq)"], "T_C": [0.0], "P_Pa": [101325.0],
    })
    with pytest.raises(KeyError, match="Phase 4|coefficients"):
        build_model_manifest(design, species_meta, coeff_df, ALLOWED_BASE, tmp_path)


# --------------------------------------------------------------------------- #
# Path resolution
# --------------------------------------------------------------------------- #
def test_resolve_model_yaml_path_order(tmp_path):
    (tmp_path / "models").mkdir()
    target = tmp_path / "models" / "M_abc.yaml"
    target.write_text("x")
    # study-relative path resolves
    assert resolve_model_yaml_path(tmp_path, "models/M_abc.yaml") == tmp_path / "models/M_abc.yaml"
    # bare name falls back to study_dir/models/<name>
    assert resolve_model_yaml_path(tmp_path, "M_abc.yaml") == target


# --------------------------------------------------------------------------- #
# Raw-long design merge
# --------------------------------------------------------------------------- #
def test_merge_raw_results_with_design():
    design = pd.DataFrame({
        "case_id": ["c1"], "study_id": ["s"], "substudy_id": ["inventory_landscape"],
        "target_variant": ["Alanine(aq)"], "model_id": ["M_x"], "H2O_mol": [1.0],
        "HCN_mol": [0.02], "C2H2_mol": [0.0], "NH3_mol": [0.0], "C2H2_over_HCN": [0.0],
        "deltaG_offset_kJ_mol": [0.0],
    })
    raw = pd.DataFrame({"scenario": ["c1", "c1"], "species": ["H2O(l)", "Alanine(aq)"],
                        "X_eq": [0.9, 0.1]})
    merged = merge_raw_results_with_design(raw, design)
    assert (merged["substudy_id"] == "inventory_landscape").all()
    assert (merged["NH3_mol"] == 0.0).all()
    assert set(merged["case_id"]) == {"c1"}


# --------------------------------------------------------------------------- #
# Variant species metadata
# --------------------------------------------------------------------------- #
def test_augment_species_preserves_base_and_properties(species_df):
    config = {
        "study": {"study_id": "t"}, "mode": {"target_products": ["Alanine(aq)"]},
        "sweeps": {"deltaG_sweep": {"enabled": True, "species": "Alanine(aq)",
                   "offsets_kJ_mol": {"type": "explicit", "values": [-40, 0, 40]},
                   "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.02}}},
    }
    offsets = build_thermo_offsets_table(config)
    augmented = augment_species_metadata_with_variants(species_df, offsets)
    base_row = species_df[species_df["cantera_name"] == "Alanine(aq)"].iloc[0]
    for variant in offsets["variant_species"]:
        vrow = augmented[augmented["cantera_name"] == variant].iloc[0]
        assert vrow["formula"] == base_row["formula"]
        assert vrow["state"] == base_row["state"]
        assert vrow["molar_volume_cm3_mol"] == base_row["molar_volume_cm3_mol"]
    # Base frame untouched.
    assert "Alanine__dG_p040(aq)" not in set(species_df["cantera_name"])


def test_augment_noop_when_no_offsets(species_df):
    empty = pd.DataFrame(columns=["thermo_variant_id", "base_species",
                                  "variant_species", "deltaG_offset_kJ_mol"])
    out = augment_species_metadata_with_variants(species_df, empty)
    assert len(out) == len(species_df)


def test_require_gibbs_seed_present_and_missing():
    assert require_gibbs_seed(
        {"species_files": {"gibbs_seed_wide_csv": "tests/fixtures/gibbs_for_fitting.csv"}}
    ) == "tests/fixtures/gibbs_for_fitting.csv"
    with pytest.raises(StudyConfigError, match="Gibbs seed CSV"):
        require_gibbs_seed({"species_files": {}})


# --------------------------------------------------------------------------- #
# Run loop (fake solver — exercises orchestration without Cantera)
# --------------------------------------------------------------------------- #
def _fake_run_single_yaml_case(fail_case=None):
    def _fake(yaml_path, scenario_id, scenario_cfg, target_product, temperature_C,
              pressure_Pa, model_mode, phase_name, solver, max_steps):
        if scenario_id == fail_case:
            raise ValueError("synthetic solver blowup")
        species = list(scenario_cfg["initial_moles"].keys()) + [target_product]
        return pd.DataFrame([{
            "scenario": scenario_id, "model_mode": model_mode, "yaml_file": "m.yaml",
            "target_product": target_product, "T_C": temperature_C,
            "T_K": temperature_C + 273.15, "P_Pa": pressure_Pa, "species": sp,
            "X_initial": 0.0, "X_eq": 0.5, "initial_moles": 0.0,
            "solver_status": "ok", "error_message": "",
        } for sp in species])
    return _fake


def test_run_loop_persists_status_and_merges_design(tmp_path, coeff_df, monkeypatch):
    import equilibrium_runner
    design = build_inventory_landscape_design(_small_inventory_config(points=2))
    design["case_id"] = [f"C{i}" for i in range(len(design))]
    design["scenario_id"] = design["case_id"]
    design["study_id"] = "test_study"
    design["substudy_id"] = "inventory_landscape"

    design_mid, models_df = compute_model_table(design, ALLOWED_BASE, coeff_df)
    run_manifest = build_run_manifest(design_mid, models_df, tmp_path)
    scenarios = {"scenarios": {
        row["case_id"]: {"initial_moles": {"H2O(l)": 1.0, "NH3(aq)": float(row["NH3_mol"])},
                         "target_products": [row["target_variant"]]}
        for _, row in design_mid.iterrows()}}

    fail_case = design_mid.iloc[1]["case_id"]
    monkeypatch.setattr(equilibrium_runner, "run_single_yaml_case",
                        _fake_run_single_yaml_case(fail_case))

    merged = run_sensitivity_manifest(
        run_manifest, scenarios, design_mid, tmp_path,
        output_long_csv=tmp_path / "results" / "raw_long.csv", progress=False)

    # run_manifest.csv persisted with ok/failed + runtimes.
    final = pd.read_csv(tmp_path / "run_manifest.csv")
    assert int((final["status"] == "ok").sum()) == 3
    assert int((final["status"] == "failed").sum()) == 1
    assert final["runtime_seconds"].notna().all()
    assert final.loc[final["case_id"] == fail_case, "error_message"].iloc[0]

    # Raw long carries design variables + identifiers; failed case retained.
    assert {"case_id", "substudy_id", "NH3_mol", "model_id", "runtime_seconds"} <= set(merged.columns)
    assert (merged["substudy_id"] == "inventory_landscape").all()
    assert fail_case in set(merged["case_id"])


def test_run_loop_limit_and_only_failed(tmp_path, coeff_df, monkeypatch):
    import equilibrium_runner
    design = build_inventory_landscape_design(_small_inventory_config(points=3))
    design["case_id"] = [f"C{i}" for i in range(len(design))]
    design["scenario_id"] = design["case_id"]
    design["study_id"] = "test_study"
    design["substudy_id"] = "inventory_landscape"
    design_mid, models_df = compute_model_table(design, ALLOWED_BASE, coeff_df)
    run_manifest = build_run_manifest(design_mid, models_df, tmp_path)
    scenarios = {"scenarios": {
        row["case_id"]: {"initial_moles": {"H2O(l)": 1.0}, "target_products": [row["target_variant"]]}
        for _, row in design_mid.iterrows()}}
    monkeypatch.setattr(equilibrium_runner, "run_single_yaml_case",
                        _fake_run_single_yaml_case())

    run_sensitivity_manifest(run_manifest, scenarios, design_mid, tmp_path,
                             output_long_csv=tmp_path / "r.csv", limit=4, progress=False)
    final = pd.read_csv(tmp_path / "run_manifest.csv")
    assert int((final["status"] == "ok").sum()) == 4
    assert int((final["status"] == "pending").sum()) == 5  # 9 - 4 limited


# --------------------------------------------------------------------------- #
# Resume-safety: manifest merge, raw union, force/only-failed, full design matrix
# --------------------------------------------------------------------------- #
def _recording_fake(calls, fail_case=None):
    """Fake solver that records which case_ids it was asked to run."""
    def _fake(yaml_path, scenario_id, scenario_cfg, target_product, temperature_C,
              pressure_Pa, model_mode, phase_name, solver, max_steps):
        calls.append(scenario_id)
        if scenario_id == fail_case:
            raise ValueError("synthetic failure")
        species = list(scenario_cfg["initial_moles"].keys()) + [target_product]
        return pd.DataFrame([{
            "scenario": scenario_id, "model_mode": model_mode, "yaml_file": "m.yaml",
            "target_product": target_product, "T_C": temperature_C,
            "T_K": temperature_C + 273.15, "P_Pa": pressure_Pa, "species": sp,
            "X_initial": 0.0, "X_eq": 0.5, "initial_moles": 0.0,
            "solver_status": "ok", "error_message": "",
        } for sp in species])
    return _fake


def _setup_run(tmp_path, coeff_df, points=5):
    design = build_inventory_landscape_design(_small_inventory_config(points))
    design["case_id"] = [f"C{i}" for i in range(len(design))]
    design["scenario_id"] = design["case_id"]
    design["study_id"] = "test_study"
    design["substudy_id"] = "inventory_landscape"
    design_mid, models_df = compute_model_table(design, ALLOWED_BASE, coeff_df)
    build_run_manifest(design_mid, models_df, tmp_path)
    scenarios = {"scenarios": {
        r["case_id"]: {"initial_moles": {"H2O(l)": 1.0}, "target_products": [r["target_variant"]]}
        for _, r in design_mid.iterrows()}}
    return design_mid, models_df, scenarios


def test_select_cases_modes():
    m = pd.DataFrame({"case_id": ["a", "b", "c"],
                      "status": ["ok", "failed", "pending"]})
    d = pd.DataFrame({"case_id": ["a", "b", "c"], "substudy_id": ["x", "x", "x"]})
    assert set(select_cases_to_run(m, d)["case_id"]) == {"b", "c"}                 # resume
    assert set(select_cases_to_run(m, d, only_failed=True)["case_id"]) == {"b"}     # failed only
    assert set(select_cases_to_run(m, d, force=True)["case_id"]) == {"a", "b", "c"}  # force


def test_partial_then_resume_unions_raw(tmp_path, coeff_df, monkeypatch):
    import equilibrium_runner
    monkeypatch.setattr(equilibrium_runner, "run_single_yaml_case", _recording_fake([]))
    design_mid, models_df, scenarios = _setup_run(tmp_path, coeff_df, points=5)  # 25 cases
    out = tmp_path / "results" / "raw.csv"
    rm = pd.read_csv(tmp_path / "run_manifest.csv")

    run_sensitivity_manifest(rm, scenarios, design_mid, tmp_path, out, limit=10, progress=False)
    assert pd.read_csv(out)["case_id"].nunique() == 10

    rm2 = build_run_manifest(design_mid, models_df, tmp_path)   # resume-safe merge
    assert int((rm2["status"] == "ok").sum()) == 10             # prior ok preserved
    run_sensitivity_manifest(rm2, scenarios, design_mid, tmp_path, out, limit=10, progress=False)
    assert pd.read_csv(out)["case_id"].nunique() == 20          # union, not just batch 2


def test_rerun_replaces_not_duplicates(tmp_path, coeff_df, monkeypatch):
    import equilibrium_runner
    monkeypatch.setattr(equilibrium_runner, "run_single_yaml_case", _recording_fake([]))
    design_mid, models_df, scenarios = _setup_run(tmp_path, coeff_df, points=4)  # 16 cases
    out = tmp_path / "results" / "raw.csv"
    rm = pd.read_csv(tmp_path / "run_manifest.csv")
    run_sensitivity_manifest(rm, scenarios, design_mid, tmp_path, out, progress=False)
    raw1 = pd.read_csv(out)
    assert raw1["case_id"].nunique() == 16

    rm2 = build_run_manifest(design_mid, models_df, tmp_path)
    run_sensitivity_manifest(rm2, scenarios, design_mid, tmp_path, out,
                             force=True, limit=1, progress=False)
    raw2 = pd.read_csv(out)
    assert raw2["case_id"].nunique() == 16                      # unchanged
    assert len(raw2) == len(raw1)                               # no duplicate rows
    assert (raw2["case_id"] == "C0").sum() == (raw1["case_id"] == "C0").sum()


def test_only_failed_reruns_only_failures(tmp_path, coeff_df, monkeypatch):
    import equilibrium_runner
    calls = []
    monkeypatch.setattr(equilibrium_runner, "run_single_yaml_case",
                        _recording_fake(calls, fail_case="C3"))
    design_mid, models_df, scenarios = _setup_run(tmp_path, coeff_df, points=3)  # 9 cases
    out = tmp_path / "results" / "raw.csv"
    rm = pd.read_csv(tmp_path / "run_manifest.csv")
    run_sensitivity_manifest(rm, scenarios, design_mid, tmp_path, out, progress=False)
    assert int((pd.read_csv(tmp_path / "run_manifest.csv")["status"] == "failed").sum()) == 1

    calls.clear()
    monkeypatch.setattr(equilibrium_runner, "run_single_yaml_case", _recording_fake(calls))
    rm2 = build_run_manifest(design_mid, models_df, tmp_path)
    run_sensitivity_manifest(rm2, scenarios, design_mid, tmp_path, out,
                             only_failed=True, progress=False)
    assert calls == ["C3"]                                      # only the failed case reran
    assert int((pd.read_csv(tmp_path / "run_manifest.csv")["status"] == "ok").sum()) == 9


def test_plain_resume_skips_ok_and_force_reruns_all(tmp_path, coeff_df, monkeypatch):
    import equilibrium_runner
    calls = []
    monkeypatch.setattr(equilibrium_runner, "run_single_yaml_case", _recording_fake(calls))
    design_mid, models_df, scenarios = _setup_run(tmp_path, coeff_df, points=3)  # 9 cases
    out = tmp_path / "results" / "raw.csv"
    run_sensitivity_manifest(pd.read_csv(tmp_path / "run_manifest.csv"),
                             scenarios, design_mid, tmp_path, out, progress=False)
    assert len(calls) == 9

    calls.clear()
    rm2 = build_run_manifest(design_mid, models_df, tmp_path)
    run_sensitivity_manifest(rm2, scenarios, design_mid, tmp_path, out, progress=False)
    assert calls == []                                          # resume skips all-ok

    calls.clear()
    rm3 = build_run_manifest(design_mid, models_df, tmp_path)
    run_sensitivity_manifest(rm3, scenarios, design_mid, tmp_path, out,
                             force=True, progress=False)
    assert len(calls) == 9                                      # force reruns everything


def test_raw_long_has_single_model_id_column(tmp_path, coeff_df, monkeypatch):
    import equilibrium_runner
    monkeypatch.setattr(equilibrium_runner, "run_single_yaml_case", _recording_fake([]))
    design_mid, models_df, scenarios = _setup_run(tmp_path, coeff_df, points=2)
    out = tmp_path / "results" / "raw.csv"
    merged = run_sensitivity_manifest(pd.read_csv(tmp_path / "run_manifest.csv"),
                                      scenarios, design_mid, tmp_path, out, progress=False)
    assert list(merged.columns).count("model_id") == 1
    assert "model_id_design" not in merged.columns
    assert "model_id_design" not in pd.read_csv(out).columns


def test_merge_drops_stale_columns_from_prior_file(tmp_path, coeff_df, monkeypatch):
    # A raw_long written by older code may carry a stale model_id_design column; a
    # rerun must not reintroduce it via the column union.
    import equilibrium_runner
    monkeypatch.setattr(equilibrium_runner, "run_single_yaml_case", _recording_fake([]))
    design_mid, models_df, scenarios = _setup_run(tmp_path, coeff_df, points=2)  # 4 cases
    out = tmp_path / "results" / "raw.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"case_id": ["C0", "C1"], "species": ["H2O(l)", "H2O(l)"],
                  "model_id": ["M", "M"], "model_id_design": ["M", "M"],
                  "X_eq": [0.1, 0.2]}).to_csv(out, index=False)

    run_sensitivity_manifest(pd.read_csv(tmp_path / "run_manifest.csv"),
                             scenarios, design_mid, tmp_path, out, force=True, progress=False)
    cols = pd.read_csv(out).columns
    assert "model_id_design" not in cols
    assert list(cols).count("model_id") == 1


def test_design_matrix_full_while_manifest_is_subset(tmp_path, species_df, coeff_df):
    # Full design has all 3 substudies; only the inventory (base-thermo) cases are runnable.
    cfg = {
        "study": {"study_id": "alanine_mvp", "output_dir": str(tmp_path / "study")},
        "mode": {"type": "single_product_sensitivity", "target_products": ["Alanine(aq)"]},
        "base_conditions": {"temperature_C": 0, "pressure_Pa": 101325},
        "model": {"allowed_species": ALLOWED_BASE},
        "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.02},
        "sweeps": {
            "inventory_landscape": {"enabled": True, "variables": {
                "NH3(aq)": {"type": "linear", "min": 0, "max": 0.15, "points": 3},
                "C2H2_over_HCN": {"type": "linear", "min": 0, "max": 5, "points": 3}}},
            "deltaG_sweep": {"enabled": True, "species": "Alanine(aq)",
                "offsets_kJ_mol": {"type": "explicit", "values": [-40, 0, 40]},
                "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.02,
                                    "C2H2(aq)": 0.042, "NH3(aq)": 0.05}},
            "nh3_deltaG_landscape": {"enabled": True, "variables": {
                "NH3(aq)": {"type": "linear", "min": 0, "max": 0.15, "points": 3},
                "deltaG_offset_kJ_mol": {"type": "explicit", "values": [-40, 0, 40]}},
                "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.02, "C2H2(aq)": 0.042}},
        },
    }
    study_dir = tmp_path / "study"
    full = build_full_design_matrix(cfg)
    write_design_outputs(full, cfg)                            # always the FULL study
    dm = pd.read_csv(study_dir / "design_matrix.csv")
    assert len(dm) == 9 + 3 + 9
    assert set(dm["substudy_id"].unique()) == {
        "inventory_landscape", "deltaG_sweep", "nh3_deltaG_landscape"}

    runnable = full[full["target_variant"].isin(set(coeff_df["cantera_name"]))]
    design_mid, models_df, _ = build_model_manifest(
        runnable, species_df, coeff_df, ALLOWED_BASE, study_dir)
    rm = build_run_manifest(design_mid, models_df, study_dir)
    assert len(rm) == 9                                        # manifest = runnable subset
    assert set(rm["case_id"]) == set(runnable["case_id"])
    assert len(dm) > len(rm)                                   # design stays full


# --------------------------------------------------------------------------- #
# Phase 4 fix-ups: target_product identity, case_hash staleness, seed hashing
# --------------------------------------------------------------------------- #
def test_target_product_preserved_for_dg_cases(tmp_path, coeff_df, monkeypatch):
    import equilibrium_runner
    monkeypatch.setattr(equilibrium_runner, "run_single_yaml_case", _recording_fake([]))
    design = pd.DataFrame({
        "case_id": ["B0", "D0"], "scenario_id": ["B0", "D0"], "study_id": "s",
        "substudy_id": ["inventory_landscape", "deltaG_sweep"],
        "target_product": ["Alanine(aq)", "Alanine(aq)"],
        "target_variant": ["Alanine(aq)", "Alanine__dG_p040(aq)"],
        "model_id": ["M1", "M2"], "T_C": [0.0, 0.0], "P_Pa": [101325.0, 101325.0],
        "H2O_mol": [1.0, 1.0], "HCN_mol": [0.02, 0.02], "C2H2_mol": [0.0, 0.042],
        "NH3_mol": [0.0, 0.05], "C2H2_over_HCN": [0.0, 2.1],
        "deltaG_offset_kJ_mol": [0.0, 40.0],
    })
    models_df = pd.DataFrame({"model_id": ["M1", "M2"],
                              "yaml_path": ["models/M1.yaml", "models/M2.yaml"]})
    rm = build_run_manifest(design, models_df, tmp_path)
    scenarios = {"scenarios": {
        r["case_id"]: {"initial_moles": {"H2O(l)": 1.0}, "target_products": [r["target_variant"]]}
        for _, r in design.iterrows()}}
    merged = run_sensitivity_manifest(rm, scenarios, design, tmp_path,
                                      tmp_path / "r.csv", progress=False)
    b = merged[merged["case_id"] == "B0"]
    d = merged[merged["case_id"] == "D0"]
    # ΔG case: product stays the original; only the variant holds the pseudo-species.
    assert (d["target_product"] == "Alanine(aq)").all()
    assert d["target_variant"].str.startswith("Alanine__dG_").all()
    # base case: product and variant are equal.
    assert (b["target_product"] == "Alanine(aq)").all()
    assert (b["target_variant"] == "Alanine(aq)").all()


def test_stale_case_hash_resets_to_pending(tmp_path, coeff_df):
    design = build_inventory_landscape_design(_small_inventory_config(points=2))
    design["case_id"] = [f"C{i}" for i in range(len(design))]
    design["scenario_id"] = design["case_id"]
    design["study_id"] = "s"
    design["substudy_id"] = "inventory_landscape"
    design_mid, models_df = compute_model_table(design, ALLOWED_BASE, coeff_df)
    build_run_manifest(design_mid, models_df, tmp_path)
    rm = pd.read_csv(tmp_path / "run_manifest.csv")
    rm["status"] = "ok"
    rm.to_csv(tmp_path / "run_manifest.csv", index=False)

    # Edit one case's inventory (same case_id) -> its case_hash changes.
    changed = design_mid.copy()
    changed.loc[changed["case_id"] == "C0", "NH3_mol"] = 99.0
    rm2 = build_run_manifest(changed, models_df, tmp_path)
    assert rm2.loc[rm2["case_id"] == "C0", "status"].iloc[0] == "pending"   # stale reset
    assert (rm2.loc[rm2["case_id"] != "C0", "status"] == "ok").all()        # others kept


def test_provenance_records_gibbs_seed_hash(tmp_path):
    config = {"study": {"study_id": "t"},
              "species_files": {"species_csv": "inputs/species_example.csv",
                                "gibbs_seed_wide_csv": "tests/fixtures/gibbs_for_fitting.csv"}}
    prov = json.loads(write_run_provenance(tmp_path, "c.yaml", config, PROJECT_ROOT).read_text())
    assert prov["gibbs_seed_sha256"] is not None


def test_seed_is_stale(tmp_path):
    from sensitivity_runner import _file_sha256
    seed = PROJECT_ROOT / "tests" / "fixtures" / "gibbs_for_fitting.csv"
    assert seed_is_stale(tmp_path, seed) is False              # no prior provenance
    (tmp_path / "run_provenance.json").write_text(json.dumps({"gibbs_seed_sha256": "deadbeef"}))
    assert seed_is_stale(tmp_path, seed) is True               # hash differs
    (tmp_path / "run_provenance.json").write_text(
        json.dumps({"gibbs_seed_sha256": _file_sha256(seed)}))
    assert seed_is_stale(tmp_path, seed) is False              # hash matches


# --------------------------------------------------------------------------- #
# Provenance
# --------------------------------------------------------------------------- #
def test_provenance_has_required_keys(tmp_path):
    config = {
        "study": {"study_id": "test_study"},
        "species_files": {"species_csv": "inputs/species_example.csv",
                          "gibbs_seed_wide_csv": "tests/fixtures/gibbs_for_fitting.csv"},
    }
    out = write_run_provenance(tmp_path, "cfg.yaml", config, PROJECT_ROOT)
    prov = json.loads(out.read_text())
    assert prov["study_id"] == "test_study"
    assert "timestamp_utc" in prov
    assert "cantera" in prov["versions"]
    assert prov["input_hashes"]["species_csv"] is not None
