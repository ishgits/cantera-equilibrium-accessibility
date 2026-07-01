# Cantera Equilibrium Accessibility v2.0.0 — Sensitivity Landscapes

## TL;DR

Version **2.0.0** adds a full **sensitivity landscape layer** on top of the original
equilibrium-accessibility workflow.

The original v1 workflow asked:

> Under one chosen inventory and thermochemical dataset, is a target species
> equilibrium-accessible?

The new v2 workflow asks:

> - Is that accessibility prediction stable if the starting inventory changes?
> - Is it stable if the target Gibbs free energy is shifted?
> - Which predictions are robust, which are inventory-gated, and which are
>   thermochemically fragile?

This release adds a study framework, two sensitivity sweep types, a second notebook,
command-line tools, bundled alanine MVP outputs, amino-acid campaign infrastructure,
machine-readable outputs, and expanded documentation.

Importantly: **the v1 results remain byte-for-byte identical.** The sensitivity layer
is additive and does not change the original single-scenario equilibrium workflow.

---

## What's new

### 1. Study-based workflow structure

v2 introduces a self-contained `studies/` framework. Each study lives in its own folder
(`studies/<study_id>/`) and contains its own:

```text
study_config.yaml
design_matrix.csv
generated_scenarios.yaml
thermo_offsets.csv
run_manifest.csv
model_manifest.csv
results/
figures/
processed/
models/
```

The design goal is simple:

> Edit one config file, generate many cases, run/resume safely, summarize, plot, and
> preserve enough provenance for downstream analysis.

The bundled example study is `studies/alanine_mvp/`, committed with completed outputs so
users can inspect the workflow immediately.

### 2. Two sensitivity sweep types

**Inventory landscapes** vary the initial abundances of feedstock species such as
`HCN(aq)`, `C2H2(aq)`, `NH3(aq)`, and `H2O(l)`. These sweeps show whether a target
species is accessible only under narrow starting conditions, or whether accessibility is
robust across broad starting inventories. The main visual output is
`inventory_landscape.png` / `.pdf`; for two-axis designs the plot reads like a chemical
"phase map" of equilibrium accessibility.

**ΔG sweeps** apply exact Gibbs free energy offsets to the target species, testing whether
the accessibility call depends delicately on the target's assumed thermochemistry. The
main visual output is `deltaG_sweep.png` / `.pdf`. These are useful when a target's
thermochemical data come from estimates, analogs, or quantum-chemistry-derived values
with uncertainty.

### 3. New sensitivity notebook

v2 adds `notebooks/02_sensitivity_landscape_workflow.ipynb`, a thin notebook wrapper
around the same tested command-line engine. It lets users:

1. Select a `study_config.yaml`
2. Preview the design matrix
3. Run or resume the study
4. Summarize results
5. Generate plots
6. Read an automated MVP verdict

The notebook deliberately avoids duplicating workflow logic — the real implementation
lives in `src/` and `scripts/`, so notebook and CLI behavior stay aligned.

v2 also adds a **teaching-first entry point**,
`notebooks/00_start_here_sensitivity_landscapes.ipynb`. It reads the committed alanine
MVP outputs and walks through how to interpret inventory landscapes and ΔG sweeps —
`Run All` works even without Cantera installed, since it runs no equilibrium
calculations. This is the recommended starting point for anyone new to the sensitivity
layer.

### 4. Command-line interface for reproducible runs

v2 adds a full CLI pathway for sensitivity studies. Core commands:

```bash
python scripts/new_study.py --id my_study
python scripts/generate_sensitivity_design.py --config studies/my_study/study_config.yaml --dry-run
python scripts/run_sensitivity_study.py --config studies/my_study/study_config.yaml
python scripts/summarize_sensitivity_study.py --config studies/my_study/study_config.yaml
python scripts/plot_sensitivity_study.py --config studies/my_study/study_config.yaml
```

The runner is resume-aware by default; completed cases are skipped unless explicitly
rerun. Useful options include `--dry-run`, `--force`, `--only-failed`, `--limit`,
`--substudy`, and `--rebuild-coefficients`. This makes the workflow practical for large
parameter sweeps and long-running equilibrium campaigns.

### 5. Bundled alanine MVP

The release includes a completed alanine sensitivity MVP at `studies/alanine_mvp/`, with
committed outputs including:

```text
results/sensitivity_case_summary.csv
results/sensitivity_landscape_grid.csv
results/sensitivity_run_summary.md
figures/inventory_landscape.png
figures/deltaG_sweep.png
figures/nh3_deltaG_landscape.png
```

This lets users explore sensitivity outputs without first installing Cantera or running
equilibrium calculations. The alanine MVP demonstrates how to interpret robust
accessibility, inventory-gated accessibility, thermochemical fragility, below-threshold
cases, trace equilibrium formation, and target mole-fraction trends across a design
matrix.

### 6. Amino-acid breadth campaign infrastructure

v2 adds infrastructure for running one sensitivity study per amino acid across a shared
dataset. Core scripts:

```text
scripts/new_amino_acid_batch.py
scripts/run_amino_acid_batch.py
scripts/compare_amino_acids.py
scripts/plot_nh3_combined.py
scripts/run_paper_extension.py
```

The workflow supports two campaign modes: `studies/aa_no_nh3/` (no-ammonia batch, which
excludes `NH3(aq)` from the phase entirely) and `studies/aa_nh3/` (ammonia-present batch,
which starts NH3 sweeps at nonzero NH3 values). The campaign layer writes per-target
summaries and aggregate comparison outputs, including machine-readable metrics and
cross-amino-acid figures. **These campaign outputs are committed as part of the published
dataset.**

### 7. ML-ready outputs

v2 expands the output schema for downstream analysis. Key files:

```text
design_matrix.csv
run_manifest.csv
model_manifest.csv
thermo_offsets.csv
results/equilibrium_raw_long.csv
results/equilibrium_moles_long.csv
results/sensitivity_case_summary.csv
results/sensitivity_landscape_grid.csv
results/schema.json
results/SCHEMA.md
```

These files are intentionally table-oriented and case-indexed so they can be used for
sensitivity analysis, plotting, downstream statistical analysis, machine-learning
datasets, workflow audit trails, and reproducibility checks. Each case has a stable
`case_id`, and model variants are tracked through manifest files.

### 8. Expanded documentation

v2 adds or expands:

```text
docs/USAGE.md
docs/SENSITIVITY_CLI.md
docs/scientific_insights_extension.md
```

These documents explain how to run the original single-scenario workflow, scaffold a new
sensitivity study, run the CLI tools, inspect generated outputs, reproduce the bundled
alanine MVP, run amino-acid campaign studies, and interpret sensitivity classifications.

---

## Compatibility with v1

v2.0.0 is fully backward-compatible with the v1 workflow. The original v1
equilibrium-accessibility pathway remains unchanged:

```text
notebooks/01_full_cantera_equilibrium_workflow_v4.ipynb
scripts/validate_static_workflow.py
inputs/species_example.csv
inputs/scenarios_example.yaml
inputs/example_validation_gibbs.csv
```

The v2 sensitivity layer is additive; it does not change the original single-scenario
model behavior. For the validation example, the v1 outputs remain **byte-for-byte
identical** relative to the v1 workflow. Existing users can continue using the original
workflow exactly as before, while adopting the new sensitivity tools only when needed.

---

## Who should use v2.0.0?

Use v2.0.0 if you want to know whether an equilibrium-accessibility result is robust. This
is especially useful when:

- starting inventories are uncertain
- feedstock ratios are approximate
- target thermochemistry is estimated
- Gibbs free energies come from analogs or quantum chemistry
- you want to compare many related target molecules
- you want machine-readable outputs for larger analysis
- you want to identify which predictions are stable enough to trust

Use the original v1-style notebook if you only need a single equilibrium-accessibility
screen.

---

## Recommended starting points

- **Learn to read sensitivity landscapes (no Cantera needed)** —
  `notebooks/00_start_here_sensitivity_landscapes.ipynb`. A teaching notebook that reads
  the committed alanine MVP outputs; `Run All` to learn interpretation without running
  anything.
- **New users** — `notebooks/01_full_cantera_equilibrium_workflow_v4.ipynb`. Start here to
  understand the original equilibrium-accessibility workflow.
- **Sensitivity studies** — `notebooks/02_sensitivity_landscape_workflow.ipynb`. Use this
  to run or inspect a sensitivity study.
- **Command-line users** — `docs/SENSITIVITY_CLI.md`. The main reference for generating,
  running, summarizing, and plotting sensitivity studies.
- **An already-computed example** — `studies/alanine_mvp/`. Browse the committed alanine
  MVP outputs and figures.

---

## Scientific framing

This workflow evaluates **equilibrium accessibility, not kinetics**. An accessible species
is one that reaches a target mole-fraction threshold under the modeled equilibrium
assumptions. This does **not** imply a reaction pathway, reaction rate, experimental yield,
detectability, or guaranteed formation in a real system.

The sensitivity layer helps answer a narrower but important question:

> Is the equilibrium-accessibility prediction stable across plausible modeling assumptions?

That makes v2 most useful as a screening, prioritization, and robustness-analysis tool.
