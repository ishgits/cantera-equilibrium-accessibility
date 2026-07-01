"""Static validation for the generated workflow.

This does not require Cantera or pyCHNOSZ. It checks metadata, scenarios,
cache/coefficients when present, manifest portability, and generated YAML
structure when YAMLs have already been generated.
"""
from pathlib import Path
import sys

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from config_io import get_project_paths, load_species_metadata, load_scenarios, list_target_products, base_species_for_scenario
from chnosz_cache import load_cache
from thermo_fit import load_coefficients

paths = get_project_paths(PROJECT_ROOT)
species_path = paths.inputs / "species_example.csv"
scenario_path = paths.inputs / "scenarios_example.yaml"
cache_path = paths.data_raw / "chnosz_gibbs_cache.csv"
coeff_path = paths.data_processed / "nasa9_coefficients.csv"
manifest_path = paths.models_single / "single_product_manifest.csv"

species = load_species_metadata(species_path)
scenarios = load_scenarios(scenario_path)
target_products = list_target_products(species)

assert not species.empty, "species metadata is empty"
assert scenarios.get("scenarios"), "no scenarios found"
assert target_products, "no target products found; mark rows as role=product or role=target"

species_names = set(species["cantera_name"])
for scenario_id, cfg in scenarios["scenarios"].items():
    base_species = base_species_for_scenario(cfg)
    missing = set(base_species) - species_names
    assert not missing, f"Scenario {scenario_id!r} references species absent from metadata: {sorted(missing)}"
    for sp, amount in cfg["initial_moles"].items():
        assert float(amount) >= 0, f"Scenario {scenario_id!r} has negative initial moles for {sp!r}"

if cache_path.exists():
    cache = load_cache(cache_path)
    assert not cache.empty, "CHNOSZ cache exists but is empty"
    print(f"Cache rows: {len(cache)}")
else:
    print("No CHNOSZ cache found yet; this is OK before first extraction.")

if coeff_path.exists():
    coeffs = load_coefficients(coeff_path)
    assert not coeffs.empty, "coefficient table exists but is empty"
    coeff_species = set(coeffs["cantera_name"].unique())
    missing_coeffs = species_names - coeff_species
    assert not missing_coeffs, f"Coefficient table is missing species: {sorted(missing_coeffs)}"
    print(f"Coefficient rows: {len(coeffs)}")
else:
    print("No NASA9 coefficient table found yet; this is OK before first fit.")

if manifest_path.exists():
    manifest = pd.read_csv(manifest_path)
    assert not manifest.empty, "single-product manifest exists but is empty"
    required = {"scenario", "model_mode", "yaml_path", "yaml_file", "target_product", "allowed_base_species"}
    missing_cols = required - set(manifest.columns)
    assert not missing_cols, f"single-product manifest is missing columns: {sorted(missing_cols)}"

    for _, row in manifest.iterrows():
        scenario_id = row["scenario"]
        yaml_file = row["yaml_file"]
        yaml_path = Path(str(row["yaml_path"]))
        if not yaml_path.exists():
            yaml_path = paths.models_single / scenario_id / yaml_file
        assert yaml_path.exists(), f"Missing YAML for manifest row: {scenario_id}/{yaml_file}"

        doc = yaml.safe_load(yaml_path.read_text())
        phase_species = doc["phases"][0]["species"]
        expected_base = [s for s in str(row["allowed_base_species"]).split(";") if s]
        target = row["target_product"]
        for sp in expected_base + [target]:
            assert sp in phase_species, f"{sp!r} missing from {yaml_path}"

        unexpected_initial_species = set(base_species_for_scenario(scenarios["scenarios"][scenario_id])) - set(phase_species)
        assert not unexpected_initial_species, (
            f"Scenario species missing from YAML {yaml_path}: {sorted(unexpected_initial_species)}"
        )
    print(f"YAML cases: {len(manifest)}")
else:
    print("No single-product manifest found yet; this is OK before YAML generation.")

print("Static workflow validation passed.")
print(f"Species: {len(species)}")
print(f"Scenarios: {len(scenarios['scenarios'])}")
print(f"Target products: {len(target_products)}")
