# Sensitivity Layer — CLI Cheat Sheet

Quick reference for running and testing the sensitivity workflow from the terminal.
Run everything **from the repo root** (`cantera_equilibrium_workflow/`).
`conftest.py` puts `src/` on the path automatically, so no install/`PYTHONPATH` needed.

---

## One-time setup

```bash
# Core Python deps
pip install -r requirements.txt
pip install pytest                      # for the test suite

# Cantera (only needed once you run equilibrium — Phase 2+). conda-forge is reliable:
conda install -c conda-forge cantera
```

---

## Tests

```bash
# Full suite (base workflow + sensitivity)
pytest -q
#   ...or, if the `pytest` command isn't on PATH:
python -m pytest -q

# Just the Phase-1 design tests, verbose (see every case)
pytest tests/test_sensitivity_design.py -v

# Run a single test by name
pytest tests/test_sensitivity_design.py -k "case_id" -v

# Base-workflow sanity check (no Cantera / pyCHNOSZ needed)
python scripts/validate_static_workflow.py
```

---

## Phase 1 — Design generation  ✅ implemented

```bash
# Preview ONLY: case counts, grid ranges, ΔG variants. Writes nothing.
python scripts/generate_sensitivity_design.py \
  --config studies/alanine_mvp/study_config.yaml --dry-run

# Generate the design outputs for real:
#   studies/alanine_mvp/design_matrix.csv
#   studies/alanine_mvp/generated_scenarios.yaml
#   studies/alanine_mvp/thermo_offsets.csv   (only if a ΔG substudy is enabled)
python scripts/generate_sensitivity_design.py \
  --config studies/alanine_mvp/study_config.yaml
```

### Scaffold a brand-new study

```bash
# Copies studies/_template/ -> studies/<id>/ with a runnable config + subdirs
python scripts/new_study.py --id my_study
python scripts/new_study.py --id my_study --force      # overwrite if it exists
```

---

## Run / summarize / plot  ✅ available

The full pipeline is implemented. Run needs Cantera; summarize/plot work over the
written CSVs without it. Or use the notebook
`notebooks/02_sensitivity_landscape_workflow.ipynb` (same engine, Run All).

```bash
# Run the equilibrium grid (fits base NASA9, a7-shifts ΔG variants, runs cases).
# Resume-by-default: only cases not already 'ok' run.
python scripts/run_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml --dry-run
python scripts/run_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml
python scripts/run_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml --only-failed
python scripts/run_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml --force          # rerun all
python scripts/run_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml --rebuild-coefficients
python scripts/run_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml --limit 50
python scripts/run_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml --substudy inventory_landscape

# Summarize + plot
python scripts/summarize_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml
python scripts/plot_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml
python scripts/plot_sensitivity_study.py --config studies/alanine_mvp/study_config.yaml --substudy inventory_landscape
```

---

## Inspect the generated outputs (handy one-liners)

```bash
# Row count + unique case_ids + per-substudy counts
python3 -c "import pandas as pd; d=pd.read_csv('studies/alanine_mvp/design_matrix.csv'); \
print(len(d),'rows |',d.case_id.nunique(),'unique case_id'); print(d.substudy_id.value_counts())"

# Peek at the design matrix
column -s, -t studies/alanine_mvp/design_matrix.csv | head -20

# Confirm present-but-zero edges are retained (NH3=0 / C2H2=0 kept in initial_moles)
python3 -c "import yaml; s=yaml.safe_load(open('studies/alanine_mvp/generated_scenarios.yaml'))['scenarios']; \
k=list(s)[0]; print(k, s[k]['allowed_species']); print(s[k]['initial_moles'])"

# ΔG variant table
cat studies/alanine_mvp/thermo_offsets.csv
```

---

## Typical workflow order

```text
1. (one-time) put species + G(T) data in inputs/ ; point study_config.yaml at them
2. edit studies/<id>/study_config.yaml          <- the only file you edit
3. generate_sensitivity_design.py --dry-run      <- sanity check counts/runtime
4. generate_sensitivity_design.py                <- write design + scenarios
5. run_sensitivity_study.py --dry-run / then run <- (Phase 2) equilibrium grid
6. summarize_sensitivity_study.py                <- (Phase 3) case summary + metrics
7. plot_sensitivity_study.py                     <- (Phase 3) landscape figures
8. pytest -q                                     <- re-run anytime to confirm nothing broke
```

Tip: a config error (missing species, bad path, typo in `allowed_species`) is reported
as a plain-English message, not a traceback. `--dry-run` is the fastest way to catch them.

---

## Amino-acid breadth campaign (Phases 7–9)

Run one single-target study **per amino acid** over the shared dataset
(`inputs/amino_acids_species.csv`, `inputs/amino_acids_gibbs_seed.csv`), then compare.

**NH3 convention:** `NH3 = 0` is represented **only** by the no-ammonia batch (NH3
excluded from the phase). Every NH3 *sweep* starts at **0.01** (1% of water — M&P's
units); a scaffold guard makes a 0 NH3 sweep value impossible. So we run **two batches**:

```bash
# Batch A — NH3 excluded ("0% NH3"): 1-D C2H2/HCN inventory + ΔG sweep (no NH3 axis)
python scripts/new_amino_acid_batch.py --exclude-species "NH3(aq)" --out studies/aa_no_nh3
# Batch B — NH3 present, every NH3 sweep >= 0.01 (--nh3-min, default 0.01)
python scripts/new_amino_acid_batch.py --out studies/aa_nh3

# Drive each batch through the engine (resume-by-default; one failure never aborts).
python scripts/run_amino_acid_batch.py --scan-dir studies/aa_no_nh3 --dry-run
python scripts/run_amino_acid_batch.py --scan-dir studies/aa_no_nh3
python scripts/run_amino_acid_batch.py --scan-dir studies/aa_nh3
#   subset / steps / resume: --only glycine,serine   --steps run,summarize   --force --limit N

# Cross-target aggregation (Phase 8) per batch -> <scan-dir>/aggregate/
python scripts/compare_amino_acids.py --scan-dir studies/aa_nh3
python scripts/compare_amino_acids.py --scan-dir studies/aa_no_nh3

# Combined M&P-style figure (Phase 9): accessibility vs NH3 % at a fiducial ratio
python scripts/plot_nh3_combined.py --ratio 2.1
```

Each batch run also writes `<scan-dir>/batch_summary.csv` (one row per amino acid).
The aggregator writes `amino_acid_metrics.csv` (composition + accessibility +
fragile/carbon-limited/robust tag), the concatenated case summaries, `SCHEMA.md`, a
`comparison_summary.md`, and ranked/heatmap/scatter figures. The combined script writes
`studies/_nh3_combined/` (tidy CSV + heatmap + per-amino-acid lines + summary).

`studies/aa_nh3/`, `studies/aa_no_nh3/`, and `studies/_nh3_combined/` are
machine-generated, but their outputs are **committed as part of the published dataset**
(re-included explicitly in `.gitignore`), alongside the curated `studies/alanine_mvp/`
demo. Only the transient `studies/amino_acid_scan/` scaffold stays gitignored.

### Direct paper reproduction (Madan & Pearce 2025)

For the *exact* paper fiducial (C2H2/HCN = 2.1; NH3 excluded vs NH3 = 1%..10% of
water), separate from the broad discovery sweep:

```bash
R_HOME=$(R RHOME) python scripts/run_paper_extension.py
#   --only glycine,proline   to subset
```

Writes `studies/_paper_extension/paper_extension_metrics.csv` (one row per amino acid:
accessibility with/without NH3, yield % of HCN at 1/5/10% NH3, ΔG crossing,
robustness, and a paper-group / interpretation) plus `summary.md`. The combined
figure (`scripts/plot_nh3_combined.py`) now also emits a **yield-% heatmap**
(`nh3_combined_heatmap_yield_pct_HCN.*`) for direct M&P comparison alongside the
log10 X_eq detectability heatmap. `studies/_paper_extension/` outputs are committed as
part of the published dataset.
