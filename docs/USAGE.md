# Usage guide

This guide covers installation, preparing input files, configuring and running the
notebook, and interpreting every output. If you only want to confirm the workflow
runs, see the Quickstart in the top-level [`README`](../README.md) and run the
bundled example first.

## Contents

1. [Installation](#installation)
2. [How the workflow is organized](#how-the-workflow-is-organized)
3. [Input file 1 — species metadata CSV](#input-file-1--species-metadata-csv)
4. [Input file 2 — scenarios YAML](#input-file-2--scenarios-yaml)
5. [Providing Gibbs-energy data (seeding vs CHNOSZ)](#providing-gibbs-energy-data-seeding-vs-chnosz)
6. [Configuring the notebook (Cell 1)](#configuring-the-notebook-cell-1)
7. [Running](#running)
8. [Outputs](#outputs)
9. [Interpretation](#interpretation)
10. [Troubleshooting](#troubleshooting)

---

## Installation

Core Python packages (pip):

```bash
pip install -r requirements.txt
```

**Cantera** is required to generate/validate model files and run equilibrium. It is
most reliably installed from conda-forge:

```bash
conda install -c conda-forge cantera
```

**pyCHNOSZ** is required *only* if you extract Gibbs energies from the CHNOSZ
database (see [seeding](#providing-gibbs-energy-data-seeding-vs-chnosz)). It depends
on R and the CHNOSZ R package; see the
[pyCHNOSZ project](https://github.com/worm-portal/pyCHNOSZ) for installation. If you
seed Gibbs data from a CSV instead, you do not need pyCHNOSZ.

Verify the non-Cantera plumbing at any time with the offline validator:

```bash
python scripts/validate_static_workflow.py
```

---

## How the workflow is organized

The notebook `notebooks/01_full_cantera_equilibrium_workflow_v4.ipynb` is the
orchestrator. **Only Cell 1 should need editing for normal use.** It imports
reusable functions from `src/` and runs these steps in order:

1. Load species metadata and scenarios.
2. Check/update the CHNOSZ Gibbs-energy cache (extract and/or seed).
3. Build a wide `G(T)` table for fitting.
4. Fit two-segment NASA9 polynomials per species.
5. Generate one single-product Cantera model file per (scenario × candidate product).
6. Validate the model files with Cantera.
7. Run equilibrium at the selected temperatures.
8. Reconstruct equilibrium moles from mole fractions via elemental balance.
9. Write a species-level inspection table.
10. Write the target accessibility summary (the main science table).
11. Write the reactant depletion/source diagnostic.
12. Generate combined accessibility bar charts.
13. Write a Markdown run summary.

---

## Input file 1 — species metadata CSV

A CSV describing every species the workflow may use. Default:
`inputs/species_example.csv`. Copy `inputs/species_template.csv` to start your own.

Required columns:

| Column | Meaning |
|---|---|
| `species_key` | Unique snake_case identifier (your choice). |
| `cantera_name` | Name used inside generated Cantera model files, e.g. `Benzene(aq)`. Must be unique. |
| `chnosz_name` | Name passed to `pyCHNOSZ.subcrt()` for extraction. Ignored if you seed from CSV. |
| `formula` | Expanded chemical formula, **no parentheses or hydrates**, e.g. `C5H5N5`. |
| `state` | `aq` (aqueous) or `liq` (liquid). |
| `molar_volume_cm3_mol` | Molar volume (cm³/mol) for the constant-volume equation of state. |
| `role` | See role taxonomy below. |

Optional columns: `notes` (free text) and `product_class` (a grouping label such as
`nucleobase`; surfaced in displays when present).

### Role taxonomy

`role` determines how a species is treated:

- **Starting inventory** (`solvent`, `reactant`, `additive`): species that can be
  present in a scenario's `initial_moles`.
- **Candidate products** (`product`, `target`): species the workflow tests for
  accessibility. Each gets its own single-product model file.

A species marked `product`/`target` is included as a candidate unless you override
the candidate list in Cell 1 (`TARGET_PRODUCTS`) or per scenario.

### Formula rules

Formulas are parsed by `src/formula_tools.py`. Use expanded element tokens only:
`C5H5N5`, `C4H4N2O2`, `H2O`. Parentheses, charges, and hydrate dots are **not**
supported — expand them first. Element counts may be fractional if needed.

The CSV is validated on load: missing required columns, duplicate `species_key` or
`cantera_name`, and missing/non-numeric molar volumes all raise clear errors.

---

## Input file 2 — scenarios YAML

Defines the modeled conditions. Default: `inputs/scenarios_example.yaml`. Copy
`inputs/scenarios_template.yaml` to start your own. Structure:

```yaml
scenarios:
  my_condition:
    description: "Free-text label (optional)"
    target_products:          # optional; restricts candidates for THIS scenario
      - MyProduct(aq)
    initial_moles:            # required; Cantera name -> initial moles
      H2O(l): 1.0
      MyReactant(aq): 0.001
```

Per-scenario keys:

- **`initial_moles`** (required): mapping of Cantera species name → initial moles.
  Every species named here must exist in the species CSV.
- **`target_products`** (optional): candidate products to model for this scenario.
  If omitted, the global `TARGET_PRODUCTS` (or all `product`/`target` species) is used.
- **`extra_allowed_species`** (optional): adds species to the phase without an
  initial amount.
- **`allowed_species`** (optional): overrides the phase species list entirely.
- **`description`** (optional): echoed into the run summary.

There is no special-casing of any molecule (including water or NH₃). If a species
should be available in a run, put it in `initial_moles` or `extra_allowed_species`.

---

## Providing Gibbs-energy data (seeding vs CHNOSZ)

The workflow needs standard-state Gibbs free energy `G(T)` for every species across
the fitting grid. These values are stored in a tidy cache at
`data/raw/chnosz_gibbs_cache.csv` (regenerated; starts absent). There are two ways
to populate it, and they can be combined in one run:

### Path A — Seed from a wide CSV (no pyCHNOSZ)

Set in Cell 1:

```python
SEED_CACHE_FROM_EXISTING_WIDE_CSV = True
SEED_WIDE_FILENAME = "example_validation_gibbs.csv"
SEED_COLUMN_NAME_MAP = {}   # map CSV column -> Cantera name; empty if headers already match
```

The wide CSV has a temperature column (`T_K`) and one `G_J_mol` column per species.
The bundled `inputs/example_validation_gibbs.csv` demonstrates the format (its values
come from public CHNOSZ). Use `SEED_COLUMN_NAME_MAP` when your column headers differ
from the Cantera names in your species CSV.

### Path B — Extract from CHNOSZ (requires pyCHNOSZ)

Set in Cell 1:

```python
RUN_CHNOSZ_EXTRACTION = True
```

The workflow calls `pyCHNOSZ.subcrt(property="G")` for any species-temperature rows
missing from the cache, using each species' `chnosz_name`. Already-cached rows are
never re-queried, so you can add new species without recomputing old ones. Set
`FORCE_REEXTRACT_CHNOSZ = True` only to rebuild from scratch.

### ⚠ Critical caveat

`pyCHNOSZ` can only return species that exist in the CHNOSZ database. **Estimated or
custom species (e.g. quantum-chemistry-derived Gibbs curves) must be supplied via the
seed-CSV path** — CHNOSZ extraction will fail for them. A typical mixed study seeds
estimated species from CSV *and* extracts CHNOSZ-known species in the same run. Keep
track of provenance: accessibility results inherit the reliability of the underlying
thermodynamic data, which may differ across species.

---

## Configuring the notebook (Cell 1)

Cell 1 groups all user settings. The most commonly edited:

| Setting | Purpose |
|---|---|
| `SPECIES_FILENAME`, `SCENARIO_FILENAME` | Input files in `inputs/`. |
| `THERMO_FIT_TEMPERATURE_GRID_C` | Temperature grid for NASA9 fitting (broad). |
| `T_FIT_SPLIT_K` | Split point between the low/high NASA9 segments. |
| `EQUILIBRIUM_TEMPERATURES_C` | Temperatures to actually run equilibrium at. |
| `PRESSURE_PA` | System pressure. |
| `FORMATION_X_THRESHOLD` | Reporting threshold; below this → `below_threshold`. |
| `SIGNIFICANT_X_THRESHOLD` | Above this → `significant`; between the two → `trace`. |
| `FORMATION_N_THRESHOLD_MOL` | Minimum reconstructed moles to count as formed. |
| Seeding / CHNOSZ flags | See the section above. |
| `TARGET_PRODUCTS` | `None` uses scenario/role-based candidates; or pass an explicit list. |
| `FORCE_*` flags | Force re-extraction, re-fit, or YAML regeneration. |
| `PLOT_X_AXIS_*` | Bar-chart x-axis behavior (`"auto"` or `"fixed"`). |

The shipped values are configured for the validation example and run without
pyCHNOSZ.

---

## Running

Open the notebook and Run All:

```bash
jupyter lab notebooks/01_full_cantera_equilibrium_workflow_v4.ipynb
```

Each cell prints a diagnostic summary. Generated model files, fit figures, and
result tables are cleared and regenerated to match the current run, so the outputs
on disk always correspond to the latest configuration.

---

## Outputs

Written under `data/` and `figures/`:

| File | What it is |
|---|---|
| `data/raw/chnosz_gibbs_cache.csv` | Tidy `G(T)` cache (one row per species/temperature). |
| `data/processed/gibbs_for_fitting.csv` | Wide `G(T)` table used for fitting. |
| `data/processed/nasa9_coefficients.csv` | Fitted two-segment NASA9 coefficients. |
| `data/processed/nasa9_fit_diagnostics.csv` | Per-species fit residual statistics. |
| `models/single_product/…` + `single_product_manifest.csv` | Generated Cantera model files and their manifest. |
| `data/results/equilibrium_raw_long.csv` / `…_wide.csv` | Raw equilibrium mole fractions. |
| `data/results/equilibrium_moles_long.csv` | Mole fractions plus reconstructed moles. |
| `data/results/equilibrium_inspection_table.csv` | Species-level debugging table. |
| `data/results/target_formation_summary.csv` | **Main science table** (see below). |
| `data/results/reactant_depletion_long.csv` / `…_summary.csv` | Depletion/source diagnostic. |
| `data/results/run_summary.md` | Human-readable run summary. |
| `figures/fit_diagnostics/*.png` | Observed-vs-fit and residual plots per species. |
| `figures/equilibrium/*_combined_accessibility_barchart.pdf` (+ `.png`) | Accessibility bar charts. |

### `target_formation_summary.csv` (the main table)

| Column | Meaning |
|---|---|
| `scenario` | Modeled condition. |
| `target_product` | Candidate product for that single-product model. |
| `T_C` | Equilibrium temperature (°C). |
| `X_eq` | Equilibrium mole fraction of the candidate. |
| `n_eq_mol` | Reconstructed equilibrium moles of the candidate. |
| `formed_bool` | `True` when above threshold and the solver succeeded. |
| `formation_call` | `significant`, `trace`, `below_threshold`, or `solver_failed`. |
| `element_balance_relative_spread` | Mole-reconstruction consistency diagnostic. |

### `reactant_depletion_summary.csv`

For each run, reports which starting species changed most (`most_depleted_species`,
`max_depletion_fraction`, `most_consumed_species`, `max_consumed_mol`). The
`depletion_call` avoids inferring a source when the candidate is below threshold.
This is a redistribution diagnostic, **not** a reaction-pathway or limiting-reagent
analysis.

---

## Interpretation

Preferred language for results: **equilibrium-accessible**, **above threshold**,
**trace**, **below threshold**.

- **Accessibility ≠ yield/kinetics.** A `significant` call means the candidate is
  thermodynamically favored when explicitly allowed in the mixture. It says nothing
  about reaction rate, mechanism, or how much would form in a real experiment.
- **Single-product, no competition.** Because each candidate is modeled in
  isolation, candidates are not ranked against one another. For a single-reactant
  scenario the outcome is effectively all-or-nothing, and `X_eq` magnitude largely
  reflects stoichiometric dilution in the solvent rather than relative product
  favorability. Compare a product's call *across conditions*, not products against
  each other within a condition.
- **Mole reconstruction.** Cantera returns mole fractions; absolute moles are
  reconstructed from elemental conservation, taking the median across elements.
  `element_balance_relative_spread` near 0 means the per-element estimates agree;
  a large spread flags a run to inspect before trusting absolute moles.
- **Data provenance matters.** Results are only as reliable as the `G(T)` inputs.
  Note which species came from CHNOSZ versus estimated/seeded sources.

---

## Troubleshooting

- **`pyCHNOSZ is not installed` / extraction fails for a species.** Either install
  pyCHNOSZ, or seed that species from a wide CSV. CHNOSZ cannot return species it
  does not know (custom/estimated ones) — those must be seeded.
- **`Missing Gibbs values for requested rows`.** A species is missing cache entries
  across the fitting grid. Extract or seed the missing rows, or narrow the grid.
- **`NASA9 fit failed` / large residuals.** Ensure ≥3 temperature points on each
  side of `T_FIT_SPLIT_K`. Inspect the plots in `figures/fit_diagnostics/`. Residuals
  of a few hundred J/mol are typically negligible relative to `G` (hundreds of kJ/mol).
- **Solver failures (`solver_status != "ok"`).** Usually an ill-conditioned initial
  composition or a temperature outside the NASA9 fit range. Check the model file and
  the coefficient temperature ranges for the failing species.
- **Empty/odd results.** Run `python scripts/validate_static_workflow.py` to check
  metadata, scenarios, cache/coefficients, and generated model files without Cantera.
