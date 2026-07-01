"""Tests for the amino-acid batch scaffold + driver (Cantera-free)."""
import sys
from pathlib import Path

import pytest
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from config_io import load_species_metadata
from sensitivity_design import (
    build_full_design_matrix, load_gibbs_seed_for_config, load_species_for_config,
    load_study_config, validate_study,
)
import new_amino_acid_batch as scaffold
import run_amino_acid_batch as driver

SPECIES_CSV = "inputs/amino_acids_species.csv"
TEMPLATE = "studies/alanine_mvp"


# --------------------------------------------------------------------------- #
# Scaffold
# --------------------------------------------------------------------------- #
def test_targets_derived_from_csv_exactly_18():
    targets = scaffold.amino_acid_targets(PROJECT_ROOT / SPECIES_CSV)
    assert len(targets) == 18
    keys = {t["key"] for t in targets}
    assert "alanine" in keys and "glycine" in keys and "tryptophan" in keys
    # Every target is an amino-acid row (not a feedstock species).
    species = load_species_metadata(PROJECT_ROOT / SPECIES_CSV)
    aa_names = set(species[species["product_class"] == "amino_acid"]["cantera_name"])
    assert {t["cantera_name"] for t in targets} == aa_names


def test_generated_configs_are_valid_single_target_and_comparable(tmp_path):
    out = str(tmp_path / "scan")
    result = scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, force=True)
    assert len(result["created"]) == 18

    template = yaml.safe_load((PROJECT_ROOT / TEMPLATE / "study_config.yaml").read_text())
    for cfg_path in result["created"]:
        cfg = load_study_config(cfg_path)
        species_df = load_species_for_config(cfg, PROJECT_ROOT)
        seed_df = load_gibbs_seed_for_config(cfg, PROJECT_ROOT)
        validate_study(cfg, species_df, seed_df)               # single-target guard + checks

        assert len(cfg["mode"]["target_products"]) == 1
        target = cfg["mode"]["target_products"][0]
        assert cfg["sweeps"]["deltaG_sweep"]["species"] == target
        assert cfg["species_files"]["species_csv"] == SPECIES_CSV
        assert cfg["species_files"]["gibbs_seed_wide_csv"] == "inputs/amino_acids_gibbs_seed.csv"
        # Common feedstock / thresholds preserved (comparable landscapes).
        assert cfg["model"]["allowed_species"] == template["model"]["allowed_species"]
        assert cfg["fixed_inventory"] == template["fixed_inventory"]
        assert cfg["thresholds"] == template["thresholds"]
        # Same grids, except every NH3 sweep starts at 0.01 (> 0), never 0.
        inv = cfg["sweeps"]["inventory_landscape"]["variables"]
        assert inv["C2H2_over_HCN"] == template["sweeps"]["inventory_landscape"]["variables"]["C2H2_over_HCN"]
        assert inv["NH3(aq)"]["min"] == 0.01
        assert inv["NH3(aq)"]["max"] == 0.15 and inv["NH3(aq)"]["points"] == 25


def test_only_force_and_noop(tmp_path):
    out = str(tmp_path / "scan")
    r1 = scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, only=["glycine"])
    assert len(r1["created"]) == 1 and r1["created"][0].parent.name == "glycine"

    r2 = scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, only=["glycine"])  # no force
    assert len(r2["created"]) == 0 and len(r2["skipped"]) == 1   # idempotent

    r3 = scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, only=["glycine"], force=True)
    assert len(r3["created"]) == 1                               # overwritten


def test_amino_acid_name_with_space_is_handled(tmp_path):
    out = str(tmp_path / "scan")
    scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, only=["asparticacid"], force=True)
    cfg = yaml.safe_load((tmp_path / "scan" / "asparticacid" / "study_config.yaml").read_text())
    assert cfg["mode"]["target_products"] == ["AsparticAcid(aq)"]
    assert "aspartic acid" in cfg["plots"]["inventory_landscape"]["colorbar_label"]


# --------------------------------------------------------------------------- #
# Part 0/1: --only/--steps validation, NH3 exclusion, NH3-min guard
# --------------------------------------------------------------------------- #
def test_scaffold_rejects_unknown_only_key(tmp_path):
    with pytest.raises(ValueError, match="Unknown amino-acid key"):
        scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, str(tmp_path / "scan"), only=["bogus"])


def test_exclude_species_drops_nh3_and_validates(tmp_path):
    out = str(tmp_path / "scan")
    scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, only=["glycine"], force=True,
                              exclude_species=("NH3(aq)",))
    cfg = load_study_config(tmp_path / "scan" / "glycine" / "study_config.yaml")
    assert "NH3(aq)" not in cfg["model"]["allowed_species"]
    assert "NH3(aq)" not in cfg["sweeps"]["inventory_landscape"]["variables"]
    assert cfg["sweeps"]["nh3_deltaG_landscape"]["enabled"] is False
    assert "NH3(aq)" not in cfg["sweeps"]["deltaG_sweep"]["fixed_inventory"]
    # Must still pass validation (no swept/inventory species outside allowed_species).
    sp = load_species_for_config(cfg, PROJECT_ROOT)
    seed = load_gibbs_seed_for_config(cfg, PROJECT_ROOT)
    validate_study(cfg, sp, seed)
    design = build_full_design_matrix(cfg)
    assert set(design["substudy_id"]) == {"inventory_landscape", "deltaG_sweep"}  # nh3xΔG dropped


def test_nh3_sweep_never_generates_zero(tmp_path):
    out = str(tmp_path / "scan")
    scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, only=["glycine"], force=True)  # Batch B
    cfg = load_study_config(tmp_path / "scan" / "glycine" / "study_config.yaml")
    design = build_full_design_matrix(cfg)
    swept = design[design["substudy_id"].isin(["inventory_landscape", "nh3_deltaG_landscape"])]
    assert (swept["NH3_mol"] > 0).all()        # every NH3 sweep starts at 0.01


def test_nh3_min_must_be_positive():
    with pytest.raises(ValueError, match="nh3-min"):
        scaffold.build_study_config({}, {"key": "x", "cantera_name": "X(aq)", "display": "x"},
                                    "studies/x", "s.csv", "g.csv", nh3_min=0.0)


def test_driver_rejects_bad_steps():
    assert driver.main(["--steps", "bogus"]) == 1
    assert driver.main(["--steps", ""]) == 1


def test_driver_rejects_unknown_only(tmp_path):
    out = str(tmp_path / "scan")
    scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, only=["glycine"], force=True)
    assert driver.main(["--scan-dir", out, "--only", "bogus", "--dry-run"]) == 1


def test_paper_fiducial_config_is_exact_and_valid(tmp_path):
    import run_paper_extension as rpe
    out = str(tmp_path / "px")
    res = scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, only=["glycine"], force=True)
    rpe._fiducialize(res["created"][0], "B")
    cfg = load_study_config(res["created"][0])
    inv = cfg["sweeps"]["inventory_landscape"]["variables"]
    assert inv["C2H2_over_HCN"]["values"] == [2.1]          # exact paper ratio
    assert inv["NH3(aq)"]["values"] == rpe.NH3_SERIES
    assert all(v > 0 for v in inv["NH3(aq)"]["values"])     # NH3 never 0
    assert cfg["sweeps"]["nh3_deltaG_landscape"]["enabled"] is False
    # Passes validation and the NH3 sweep never generates 0.
    sp = load_species_for_config(cfg, PROJECT_ROOT)
    seed = load_gibbs_seed_for_config(cfg, PROJECT_ROOT)
    validate_study(cfg, sp, seed)
    design = build_full_design_matrix(cfg)
    inv_cases = design[design["substudy_id"] == "inventory_landscape"]
    assert len(inv_cases) == len(rpe.NH3_SERIES)            # 6 fiducial NH3 points
    assert (inv_cases["NH3_mol"] > 0).all()


# --------------------------------------------------------------------------- #
# Batch driver
# --------------------------------------------------------------------------- #
def test_dry_run_summary_counts(tmp_path):
    out = str(tmp_path / "scan")
    scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, only=["glycine"], force=True)
    paths = driver.discover_studies(out)
    summary = driver.dry_run_summary(paths)
    assert summary[0]["n_cases"] == 911            # 625 + 11 + 275
    assert summary[0]["n_models"] == 12            # 1 base + 11 ΔG variants


def test_batch_continues_past_a_failing_study(tmp_path):
    out = str(tmp_path / "scan")
    scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, only=["glycine", "serine"], force=True)
    paths = driver.discover_studies(out)
    assert len(paths) == 2

    calls = []

    def fake_run(argv):
        cfg = argv[argv.index("--config") + 1]
        calls.append(("run", cfg))
        return 1 if "serine" in cfg else 0      # serine's run step fails

    def fake_ok(argv):
        calls.append(("step", argv[argv.index("--config") + 1]))
        return 0

    results = driver.run_batch(paths, run_fn=fake_run, summarize_fn=fake_ok,
                               plot_fn=fake_ok, progress=False)
    by_key = {r["key"]: r for r in results}
    assert by_key["glycine"]["status"] == "ok"
    assert by_key["serine"]["status"] == "failed"
    assert by_key["serine"]["failed_step"] == "run"
    # Both studies were visited (the failure did not abort the batch).
    assert len(results) == 2


def test_batch_records_exception_without_aborting(tmp_path):
    out = str(tmp_path / "scan")
    scaffold.scaffold_studies(SPECIES_CSV, TEMPLATE, out, only=["glycine", "serine"], force=True)
    paths = driver.discover_studies(out)

    def fake_run(argv):
        if "glycine" in argv[argv.index("--config") + 1]:
            raise RuntimeError("boom")
        return 0

    results = driver.run_batch(paths, run_fn=fake_run, summarize_fn=lambda a: 0,
                               plot_fn=lambda a: 0, progress=False)
    by_key = {r["key"]: r for r in results}
    assert by_key["glycine"]["status"] == "failed" and "boom" in by_key["glycine"]["message"]
    assert by_key["serine"]["status"] == "ok"
