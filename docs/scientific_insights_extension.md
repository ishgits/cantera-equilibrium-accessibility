# Scientific insights: extending Madan & Pearce (2025)

Synthesis of the two-batch, 18-amino-acid sensitivity campaign, extending **Madan & Pearce 2025,
*PSJ* 6:284**. Same Cantera VCS engine, **single-product weak coupling** (one target pathway at a
time — the deliberate, experiment-aligned modeling choice; product–product competition is *not*
modeled), same Titan feedstock (HCN + C2H2 + NH3 in water), same 0 °C melt-pool condition.

Scope: the 18 C/H/N/O proteinogenic amino acids (the 20 standard minus sulfur-bearing
cysteine/methionine). β-alanine (non-proteinogenic) is also out of scope.

## 1. The campaign reproduces the paper's NH3 result — quantitatively
Two batches, matching the paper's two simulation sets:
- **Batch A — NH3 excluded (the paper's "0% NH3"):** only **alanine and proline** are accessible.
  This is exactly M&P's NH3-free result (alanine, β-alanine, proline) restricted to the CHNO-
  proteinogenic set — β-alanine is the only difference, and only because it is out of scope.
- **Batch B — NH3 present (≥ 1% of water):** all 18 amino acids become accessible at the fiducial
  C2H2/HCN ratio. The 0% → 1% NH3 jump reproduces the paper's sharp unlocking.
- **Yields match the paper.** Expressed as equilibrium moles relative to initial HCN at the fiducial
  inventory, the workflow reproduces M&P's percentages closely (e.g. glycine ≈ 75% vs paper 74.1%,
  phenylalanine ≈ 56.8% vs 56.7%, alanine/proline 100%). This is a quantitative reproduction, not
  just a qualitative pattern match.

## 2. New layer: thermochemical robustness, and what it implies
Perturbing each amino acid's standard Gibbs energy by an exact ±200 kJ/mol (analytic NASA9 a7 shift),
in the NH3-present batch where each amino acid is accessible at zero offset:
- **Only glycine and serine are fragile** — they lose significant accessibility under a large positive
  offset (crossings ≈ +159 and +187 kJ/mol). Every other amino acid is robust across the full
  ±200 kJ/mol range.
- The accessibility answer for most amino acids is therefore governed by **stoichiometry, allowed
  atoms, and feedstock inventory**, not by fine ΔG values. Small/moderate Gibbs uncertainty does not
  flip the binary yes/no for most targets.

**Methodological consequence (connects to the paper's QC-Ochterski method):** the sensitivity layer
identifies *when* high-effort QC thermochemistry is actually worth it — near an accessibility
threshold, where yield magnitude/ranking matters, or for fragile species. Here that points at
glycine and serine specifically; for strongly-accessible or strongly-inaccessible targets, moderate
ΔG uncertainty does not change the conclusion. The workflow is best framed as a **robustness and
prioritization engine**: fiducial simulations answer the first-order yes/no; the sweeps tell you
which conclusions are robust vs. boundary-sensitive and where deeper thermochemistry or experiments
should focus.

## 3. The present-but-zero vs. excluded distinction was the key to the comparison
The explicit phase-membership control (a species in `allowed_species` at zero initial moles can still
*form*; a species not listed cannot) is what makes the paper comparison clean. The paper's 0% NH3 is
NH3 **excluded**; allowing NH3 to *form* (even from zero initial NH3) is a different condition. The
two-batch design enforces this: 0% NH3 is only ever the excluded batch, and NH3 sweeps start at 0.01
(> 0). Because H2O = 1.0, `NH3_mol` is the fraction of water (0.01 = 1% NH3 — the paper's units).

## 4. Caveats
- Single-product weak coupling is the deliberate, experiment-aligned model (one targeted pathway);
  these are not product-competition statements.
- The ΔG sweep is a temperature-independent constant offset (cleanest model of an uncertain reported
  Δ_fG°); an entropy-weighted perturbation would behave differently.
- Glycine and serine are exactly where the input CHNOSZ Gibbs values most affect the conclusion —
  re-verify those two before leaning on their crossing offsets.
- The combined figure should report **yield relative to HCN** (paper units) in addition to
  `log10 X_eq`, and an **exact paper-fiducial grid** (ratio 2.1; NH3 = 0.01–0.10) should be sampled
  for exact reproduction (the broad sweep currently snaps to ratio ≈ 2.083, NH3 ≈ 0.0217).
- Batch-A robustness metrics exclude baseline-inaccessible amino acids before any "robust" claim is
  read from that batch: `classify_discriminators` is accessibility-aware and only tags an amino acid
  `robust_accessible` when it is accessible at the ΔG=0 reference.
