"""Tests for the schema dictionary writer and its CLI gate (no Cantera)."""
import json
import sys
from pathlib import Path

import pandas as pd
import yaml

from sensitivity_summary import COLUMN_DESCRIPTIONS, write_schema_dictionary

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _write_result_tables(results_dir):
    results_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "case_id": ["c0", "c1"], "substudy_id": ["inventory_landscape"] * 2,
        "X_eq": [0.1, 0.2], "log10_X_eq": [-1.0, -0.7], "formed_bool": [True, True],
        "formation_call": ["significant", "significant"],
        "a_made_up_column": [1, 2],   # not in COLUMN_DESCRIPTIONS
    }).to_csv(results_dir / "sensitivity_case_summary.csv", index=False)
    pd.DataFrame({
        "case_id": ["c0", "c1"], "NH3_mol": [0.0, 0.15], "C2H2_over_HCN": [0.0, 5.0],
        "log10_X_eq": [-1.0, -0.7], "solver_status": ["ok", "ok"],
    }).to_csv(results_dir / "sensitivity_landscape_grid.csv", index=False)


def test_writer_documents_all_present_columns(tmp_path):
    _write_result_tables(tmp_path)
    md, js = write_schema_dictionary(tmp_path)
    assert md.exists() and js.exists()

    schema = json.loads(js.read_text())
    assert set(schema) == {"sensitivity_case_summary.csv", "sensitivity_landscape_grid.csv"}
    # Every column present in each CSV is documented (nothing silently dropped).
    for table in schema:
        csv_cols = list(pd.read_csv(tmp_path / table).columns)
        assert list(schema[table].keys()) == csv_cols
        for col, meta in schema[table].items():
            assert set(meta) == {"dtype", "units", "description"}

    # Markdown has a section + table per CSV.
    text = md.read_text()
    assert "## sensitivity_case_summary.csv" in text
    assert "| column | dtype | units | description |" in text


def test_undocumented_column_is_listed(tmp_path):
    _write_result_tables(tmp_path)
    _, js = write_schema_dictionary(tmp_path)
    schema = json.loads(js.read_text())
    made_up = schema["sensitivity_case_summary.csv"]["a_made_up_column"]
    assert made_up["description"] == "(undocumented)"
    assert made_up["dtype"]                      # dtype still captured
    # A documented column carries its curated units/description.
    nh3 = schema["sensitivity_landscape_grid.csv"]["NH3_mol"]
    assert nh3["units"] == "mol"
    assert nh3["description"] == COLUMN_DESCRIPTIONS["NH3_mol"][1]


def test_schema_escapes_pipes_and_types_error_message(tmp_path):
    # formation_call's curated description contains '|' chars; error_message is all-empty.
    pd.DataFrame({"case_id": ["c0"], "formation_call": ["significant"],
                  "error_message": [""]}).to_csv(
        tmp_path / "sensitivity_case_summary.csv", index=False)
    md, js = write_schema_dictionary(tmp_path)

    text = md.read_text()
    assert "\\|" in text                          # pipes inside cells are escaped
    # Every markdown table row keeps exactly the 4-column shape (5 unescaped pipes).
    for line in text.splitlines():
        if line.startswith("| ") and "column" not in line and "---" not in line:
            assert line.replace("\\|", "").count("|") == 5

    schema = json.loads(js.read_text())
    dtype = schema["sensitivity_case_summary.csv"]["error_message"]["dtype"]
    assert dtype in ("str", "string", "object")   # not float64


def test_missing_tables_are_skipped(tmp_path):
    # Only one of the canonical tables exists.
    (tmp_path).mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"case_id": ["c0"], "X_eq": [0.1]}).to_csv(
        tmp_path / "sensitivity_case_summary.csv", index=False)
    _, js = write_schema_dictionary(tmp_path)
    schema = json.loads(js.read_text())
    assert set(schema) == {"sensitivity_case_summary.csv"}   # missing ones skipped, no error


# --------------------------------------------------------------------------- #
# CLI gate: outputs.write_schema_dictionary controls emission
# --------------------------------------------------------------------------- #
def _raw_long():
    rows = []
    for cid, nh3, xa in [("K0", 0.0, 1e-3), ("K1", 0.15, 1e-2)]:
        for sp, x, im in [("H2O(l)", 1 - xa, 1.0), ("Alanine(aq)", xa, 0.0)]:
            rows.append({
                "scenario": cid, "case_id": cid, "study_id": "schema_test",
                "substudy_id": "inventory_landscape", "target_product": "Alanine(aq)",
                "target_variant": "Alanine(aq)", "model_id": "M",
                "model_mode": "single_product_sensitivity", "yaml_file": "m.yaml",
                "T_C": 0.0, "T_K": 273.15, "P_Pa": 101325.0, "species": sp,
                "X_initial": 0.0, "X_eq": x, "initial_moles": im, "solver_status": "ok",
                "error_message": "", "runtime_seconds": 0.01, "H2O_mol": 1.0,
                "HCN_mol": 0.02, "C2H2_mol": 0.0, "NH3_mol": nh3, "C2H2_over_HCN": 0.0,
                "deltaG_offset_kJ_mol": 0.0,
            })
    return pd.DataFrame(rows)


def _make_config(tmp_path, write_schema: bool):
    cfg = {
        "study": {"study_id": "schema_test", "output_dir": str(tmp_path)},
        "mode": {"type": "single_product_sensitivity", "target_products": ["Alanine(aq)"]},
        "species_files": {"species_csv": "inputs/species_example.csv",
                          "gibbs_seed_wide_csv": "tests/fixtures/gibbs_for_fitting.csv"},
        "base_conditions": {"temperature_C": 0, "pressure_Pa": 101325},
        "model": {"allowed_species": ["H2O(l)", "HCN(aq)", "C2H2(aq)", "NH3(aq)"]},
        "thresholds": {"formation_X_threshold": 1e-12, "significant_X_threshold": 1e-6,
                       "formation_n_threshold_mol": 0.0, "balance_tol": 1e-6},
        "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.02},
        "sweeps": {"inventory_landscape": {"enabled": True, "variables": {
            "NH3(aq)": {"type": "linear", "min": 0, "max": 0.15, "points": 2},
            "C2H2_over_HCN": {"type": "linear", "min": 0, "max": 5, "points": 2}}}},
        "outputs": {"write_schema_dictionary": write_schema},
    }
    path = tmp_path / "study_config.yaml"
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    return path


def _run_summarize(config_path):
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    import summarize_sensitivity_study as sss
    return sss.main(["--config", str(config_path)])


def test_cli_emits_schema_when_flag_true(tmp_path):
    (tmp_path / "results").mkdir(parents=True, exist_ok=True)
    _raw_long().to_csv(tmp_path / "results" / "equilibrium_raw_long.csv", index=False)
    rc = _run_summarize(_make_config(tmp_path, write_schema=True))
    assert rc == 0
    assert (tmp_path / "results" / "SCHEMA.md").exists()
    assert (tmp_path / "results" / "schema.json").exists()


def test_cli_skips_schema_when_flag_false(tmp_path):
    (tmp_path / "results").mkdir(parents=True, exist_ok=True)
    _raw_long().to_csv(tmp_path / "results" / "equilibrium_raw_long.csv", index=False)
    rc = _run_summarize(_make_config(tmp_path, write_schema=False))
    assert rc == 0
    assert not (tmp_path / "results" / "SCHEMA.md").exists()
    assert not (tmp_path / "results" / "schema.json").exists()
