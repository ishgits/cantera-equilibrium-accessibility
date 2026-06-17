# Run summary — cantera_equilibrium_workflow_v4

## Modeled conditions

- Species file: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/inputs/species_example.csv`
- Scenario file: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/inputs/scenarios_example.yaml`
- Pressure: `100000 Pa` (`1 bar`)
- Equilibrium temperatures (°C): `[20.0]`
- Thermo fitting temperatures (°C): `[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 370.0]`
- Formation threshold: `X_eq >= 1.000e-12`
- Significant threshold: `X_eq >= 1.000e-06`
- Minimum mole threshold: `n_eq_mol >= 0.000e+00`
- Target products modeled: `2`

## Starting inventory

### validation — Validation run: HCN + C2H2 + NH3 in water, 20 C

| Species | Initial moles |
|---|---:|
| `H2O(l)` | 1 |
| `HCN(aq)` | 0.02 |
| `C2H2(aq)` | 0.042 |
| `NH3(aq)` | 0.02 |

## Equilibrium accessibility summary

- Total target runs summarized: `2`
- significant: `2`

| scenario | target_product | T_C | X_eq | n_eq_mol | formation_call | solver_status |
|---|---|---|---|---|---|---|
| validation | Adenine(aq) | 20.0 | 3.752e-03 | 4.000e-03 | significant | ok |
| validation | Cytosine(aq) | 20.0 | 5.650e-03 | 6.000e-03 | significant | ok |

## Reactant depletion diagnostic

| scenario | target_product | T_C | formation_call | most_depleted_species | max_depletion_fraction | most_consumed_species | max_consumed_mol | depletion_call |
|---|---|---|---|---|---|---|---|---|
| validation | Adenine(aq) | 20.0 | significant | HCN(aq) | 1.000e+00 | HCN(aq) | 2.000e-02 | dominant_depletion_detected |
| validation | Cytosine(aq) | 20.0 | significant | HCN(aq) | 1.000e+00 | HCN(aq) | 2.000e-02 | dominant_depletion_detected |

## Saved outputs

- CHNOSZ Gibbs cache: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/data/raw/chnosz_gibbs_cache.csv`
- Gibbs fitting table: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/data/processed/gibbs_for_fitting.csv`
- NASA9 coefficients: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/data/processed/nasa9_coefficients.csv`
- NASA9 fit diagnostics: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/data/processed/nasa9_fit_diagnostics.csv`
- Single-product YAML manifest: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/models/single_product/single_product_manifest.csv`
- Raw equilibrium long table: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/data/results/equilibrium_raw_long.csv`
- Raw equilibrium wide table: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/data/results/equilibrium_raw_wide.csv`
- Equilibrium moles long table: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/data/results/equilibrium_moles_long.csv`
- Raw inspection table: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/data/results/equilibrium_inspection_table.csv`
- Target accessibility summary: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/data/results/target_formation_summary.csv`
- Reactant depletion long table: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/data/results/reactant_depletion_long.csv`
- Reactant depletion summary: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/data/results/reactant_depletion_summary.csv`
- Accessibility figures directory: `/Users/ishaanmadan/Project1_Thermochemical_Modeling/Kendra/cantera_equilibrium_workflow/figures/equilibrium`

## Interpretation note

These outputs report equilibrium accessibility for species explicitly included in each single-product YAML. They do not report kinetic rates, reaction pathways, or formal percent yields. The depletion diagnostic identifies which starting species changed most during equilibrium redistribution.
