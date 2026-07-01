# Result schema dictionary

Auto-generated column dictionary for the result tables in this directory. One section per table; units are blank where not applicable.

## sensitivity_case_summary.csv

| column | dtype | units | description |
|---|---|---|---|
| study_id | object |  | Study identifier (folder name under studies/). |
| substudy_id | object |  | Substudy: inventory_landscape, deltaG_sweep, or nh3_deltaG_landscape. |
| target_product | object |  | Original target product, e.g. Alanine(aq). |
| target_variant | object |  | Cantera species actually modelled (base name or ΔG pseudo-species). |
| model_id | object |  | Hashed Cantera model identity (one YAML reused across grid points). |
| T_C | float64 | deg C | Equilibrium temperature. |
| P_Pa | float64 | Pa | Equilibrium pressure. |
| runtime_seconds | float64 | s | Wall-clock runtime for the case. |
| H2O_mol | float64 | mol | Initial moles of water (solvent basis). |
| HCN_mol | float64 | mol | Initial moles of hydrogen cyanide. |
| C2H2_mol | float64 | mol | Initial moles of acetylene. |
| NH3_mol | float64 | mol | Initial moles of ammonia. |
| C2H2_over_HCN | float64 | ratio | Initial C2H2/HCN mole ratio (derived design variable). |
| deltaG_offset_kJ_mol | float64 | kJ/mol | Gibbs-energy offset applied to the target species. |
| case_id | object |  | Unique simulation case identifier (canonical run key). |
| X_eq | float64 | mole fraction | Equilibrium mole fraction of the (target) species. |
| log10_X_eq | float64 | log10 mole fraction | log10(X_eq); NaN when X_eq <= 0 or missing. |
| n_eq_mol | float64 | mol | Reconstructed equilibrium moles of the species. |
| log10_n_eq_mol | float64 | log10 mol | log10(n_eq_mol); NaN when non-positive/missing. |
| formed_bool | bool |  | True if the target is equilibrium-accessible above threshold. |
| formation_call | object |  | significant \| trace \| below_threshold \| solver_failed. |
| solver_status | object |  | Equilibrium solver outcome: ok or failed. |
| error_message | object |  | Solver/setup error message when the case failed. |
| element_balance_relative_spread | float64 |  | Relative spread of per-element total-mole estimates (QC). |
| suspect_balance | bool |  | True when element_balance_relative_spread exceeds balance_tol. |

## sensitivity_landscape_grid.csv

| column | dtype | units | description |
|---|---|---|---|
| case_id | object |  | Unique simulation case identifier (canonical run key). |
| study_id | object |  | Study identifier (folder name under studies/). |
| substudy_id | object |  | Substudy: inventory_landscape, deltaG_sweep, or nh3_deltaG_landscape. |
| target_product | object |  | Original target product, e.g. Alanine(aq). |
| target_variant | object |  | Cantera species actually modelled (base name or ΔG pseudo-species). |
| model_id | object |  | Hashed Cantera model identity (one YAML reused across grid points). |
| NH3_mol | float64 | mol | Initial moles of ammonia. |
| HCN_mol | float64 | mol | Initial moles of hydrogen cyanide. |
| C2H2_mol | float64 | mol | Initial moles of acetylene. |
| C2H2_over_HCN | float64 | ratio | Initial C2H2/HCN mole ratio (derived design variable). |
| deltaG_offset_kJ_mol | float64 | kJ/mol | Gibbs-energy offset applied to the target species. |
| T_C | float64 | deg C | Equilibrium temperature. |
| P_Pa | float64 | Pa | Equilibrium pressure. |
| X_eq | float64 | mole fraction | Equilibrium mole fraction of the (target) species. |
| log10_X_eq | float64 | log10 mole fraction | log10(X_eq); NaN when X_eq <= 0 or missing. |
| n_eq_mol | float64 | mol | Reconstructed equilibrium moles of the species. |
| log10_n_eq_mol | float64 | log10 mol | log10(n_eq_mol); NaN when non-positive/missing. |
| formed_bool | bool |  | True if the target is equilibrium-accessible above threshold. |
| formation_call | object |  | significant \| trace \| below_threshold \| solver_failed. |
| solver_status | object |  | Equilibrium solver outcome: ok or failed. |
| runtime_seconds | float64 | s | Wall-clock runtime for the case. |

## equilibrium_raw_long.csv

| column | dtype | units | description |
|---|---|---|---|
| scenario | object |  | Scenario id carried from the base runner (equals case_id). |
| model_mode | object |  | Model mode tag, e.g. single_product_sensitivity. |
| yaml_file | object |  | Cantera YAML file name used for the case. |
| target_product | object |  | Original target product, e.g. Alanine(aq). |
| T_C | float64 | deg C | Equilibrium temperature. |
| T_K | float64 | K | Equilibrium temperature in Kelvin. |
| P_Pa | float64 | Pa | Equilibrium pressure. |
| species | object |  | Cantera species name for this row (long-form tables). |
| X_initial | float64 | mole fraction | Initial mole fraction supplied to the solver. |
| X_eq | float64 | mole fraction | Equilibrium mole fraction of the (target) species. |
| initial_moles | float64 | mol | Initial moles of the species in the scenario. |
| solver_status | object |  | Equilibrium solver outcome: ok or failed. |
| error_message | object |  | Solver/setup error message when the case failed. |
| case_id | object |  | Unique simulation case identifier (canonical run key). |
| model_id | object |  | Hashed Cantera model identity (one YAML reused across grid points). |
| target_variant | object |  | Cantera species actually modelled (base name or ΔG pseudo-species). |
| runtime_seconds | float64 | s | Wall-clock runtime for the case. |
| study_id | object |  | Study identifier (folder name under studies/). |
| substudy_id | object |  | Substudy: inventory_landscape, deltaG_sweep, or nh3_deltaG_landscape. |
| H2O_mol | float64 | mol | Initial moles of water (solvent basis). |
| HCN_mol | float64 | mol | Initial moles of hydrogen cyanide. |
| C2H2_mol | float64 | mol | Initial moles of acetylene. |
| NH3_mol | float64 | mol | Initial moles of ammonia. |
| C2H2_over_HCN | float64 | ratio | Initial C2H2/HCN mole ratio (derived design variable). |
| deltaG_offset_kJ_mol | float64 | kJ/mol | Gibbs-energy offset applied to the target species. |

## equilibrium_moles_long.csv

| column | dtype | units | description |
|---|---|---|---|
| scenario | object |  | Scenario id carried from the base runner (equals case_id). |
| model_mode | object |  | Model mode tag, e.g. single_product_sensitivity. |
| yaml_file | object |  | Cantera YAML file name used for the case. |
| target_product | object |  | Original target product, e.g. Alanine(aq). |
| T_C | float64 | deg C | Equilibrium temperature. |
| T_K | float64 | K | Equilibrium temperature in Kelvin. |
| P_Pa | float64 | Pa | Equilibrium pressure. |
| species | object |  | Cantera species name for this row (long-form tables). |
| X_initial | float64 | mole fraction | Initial mole fraction supplied to the solver. |
| X_eq | float64 | mole fraction | Equilibrium mole fraction of the (target) species. |
| initial_moles | float64 | mol | Initial moles of the species in the scenario. |
| solver_status | object |  | Equilibrium solver outcome: ok or failed. |
| error_message | object |  | Solver/setup error message when the case failed. |
| case_id | object |  | Unique simulation case identifier (canonical run key). |
| model_id | object |  | Hashed Cantera model identity (one YAML reused across grid points). |
| target_variant | object |  | Cantera species actually modelled (base name or ΔG pseudo-species). |
| runtime_seconds | float64 | s | Wall-clock runtime for the case. |
| study_id | object |  | Study identifier (folder name under studies/). |
| substudy_id | object |  | Substudy: inventory_landscape, deltaG_sweep, or nh3_deltaG_landscape. |
| H2O_mol | float64 | mol | Initial moles of water (solvent basis). |
| HCN_mol | float64 | mol | Initial moles of hydrogen cyanide. |
| C2H2_mol | float64 | mol | Initial moles of acetylene. |
| NH3_mol | float64 | mol | Initial moles of ammonia. |
| C2H2_over_HCN | float64 | ratio | Initial C2H2/HCN mole ratio (derived design variable). |
| deltaG_offset_kJ_mol | float64 | kJ/mol | Gibbs-energy offset applied to the target species. |
| n_total_eq_mol | float64 | mol | Reconstructed total moles at equilibrium for the case. |
| n_eq_mol | float64 | mol | Reconstructed equilibrium moles of the species. |
| element_balance_relative_spread | float64 |  | Relative spread of per-element total-mole estimates (QC). |
| element_total_mole_estimates | object |  | Per-element total-mole estimates (diagnostic string). |

