"""CLI-level tests for run/plot scripts (Cantera-free: the run path is exercised
only up to the no-op early return)."""
import json
import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import plot_sensitivity_study  # noqa: E402
import run_sensitivity_study  # noqa: E402


def _inventory_config(study_dir, make_plots=True):
    return {
        "study": {"study_id": "cli_test", "output_dir": str(study_dir)},
        "mode": {"type": "single_product_sensitivity", "target_products": ["Alanine(aq)"]},
        "species_files": {"species_csv": "inputs/species_example.csv",
                          "gibbs_seed_wide_csv": "tests/fixtures/gibbs_for_fitting.csv"},
        "base_conditions": {"temperature_C": 0, "pressure_Pa": 101325,
                            "phase_name": "aqueous", "solver": "vcs", "max_steps": 100000},
        "model": {"allowed_species": ["H2O(l)", "HCN(aq)", "C2H2(aq)", "NH3(aq)"]},
        "thresholds": {"formation_X_threshold": 1e-12, "significant_X_threshold": 1e-6,
                       "formation_n_threshold_mol": 0.0, "balance_tol": 1e-6},
        "fixed_inventory": {"H2O(l)": 1.0, "HCN(aq)": 0.02},
        "sweeps": {"inventory_landscape": {"enabled": True, "variables": {
            "NH3(aq)": {"type": "linear", "min": 0, "max": 0.15, "points": 2},
            "C2H2_over_HCN": {"type": "linear", "min": 0, "max": 5, "points": 2}}}},
        "outputs": {"make_plots": make_plots, "write_provenance": True},
        "plots": {"formats": ["png"]},
    }


def _write_config(tmp_path, cfg):
    path = tmp_path / "study_config.yaml"
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    return path


def test_make_plots_false_skips_plotting(tmp_path):
    study_dir = tmp_path / "study"
    cfg_path = _write_config(tmp_path, _inventory_config(study_dir, make_plots=False))
    rc = plot_sensitivity_study.main(["--config", str(cfg_path)])
    assert rc == 0
    # The gate returns before any figure is produced.
    assert not list((study_dir / "figures").glob("*.png")) if (study_dir / "figures").exists() else True


def test_noop_resume_preserves_provenance(tmp_path):
    study_dir = tmp_path / "study"
    cfg_path = _write_config(tmp_path, _inventory_config(study_dir))

    # First invocation builds the design/coeffs/models/manifest (and runs the 4 cases
    # if Cantera is present). Then force a completed manifest + a sentinel provenance.
    run_sensitivity_study.main(["--config", str(cfg_path)])
    rm = pd.read_csv(study_dir / "run_manifest.csv")
    rm["status"] = "ok"
    rm["runtime_seconds"] = 0.001
    rm["error_message"] = ""
    rm.to_csv(study_dir / "run_manifest.csv", index=False)
    (study_dir / "run_provenance.json").write_text(
        json.dumps({"gibbs_seed_sha256": "x", "versions": {"cantera": "SENTINEL"}}))

    # Second invocation is a no-op resume: it must return before write_run_provenance,
    # leaving the sentinel record untouched (P4 reproducibility).
    rc = run_sensitivity_study.main(["--config", str(cfg_path)])
    assert rc == 0
    prov = json.loads((study_dir / "run_provenance.json").read_text())
    assert prov["versions"]["cantera"] == "SENTINEL"
