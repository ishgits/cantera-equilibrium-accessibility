# Cross-target metrics — column dictionary

| column | dtype | description |
|---|---|---|
| amino_acid | object | Study key (folder name). |
| target_product | object | Cantera target species. |
| formula | object | Molecular formula (CHNOSZ). |
| n_C | int64 | Carbon atoms. |
| n_H | int64 | Hydrogen atoms. |
| n_N | int64 | Nitrogen atoms. |
| n_O | int64 | Oxygen atoms. |
| molar_volume_cm3_mol | float64 | Standard partial molar volume (CHNOSZ). |
| max_stoichiometric_yield_mol | float64 | Limiting-reagent moles formable from the reference feedstock. |
| inventory_accessible_fraction | float64 | Fraction of the inventory grid where the target is accessible. |
| min_NH3_accessible | float64 | Smallest NH3 (mol) with accessibility (Batch B only). |
| min_C2H2_over_HCN_accessible | float64 | Smallest C2H2/HCN ratio with accessibility. |
| X_eq_at_reference_inventory | float64 | Equilibrium mole fraction at the ΔG-sweep reference inventory. |
| accessible_at_zero_offset | bool | (undocumented) |
| max_X_eq | float64 | Peak equilibrium mole fraction over the inventory grid. |
| peak_case_id | object | Case id of the peak. |
| deltaG_positive_crossing_kJ_mol | object | ΔG offset (+) where accessibility is lost (None if none in range). |
| deltaG_negative_crossing_kJ_mol | object | ΔG offset (-) where accessibility is lost (None if none in range). |
| robust_to_pm20 | bool | Accessible across ±20 kJ/mol of Gibbs uncertainty. |
| robust_to_pm40 | bool | Accessible across ±40 kJ/mol of Gibbs uncertainty. |
| n_failed | int64 | Solver failures. |
| n_suspect_balance | int64 | Cases flagged for element-balance spread. |
| discriminator | object | not_accessible_in_batch \| energetically_fragile \| inventory_gated \| robust_accessible. |
