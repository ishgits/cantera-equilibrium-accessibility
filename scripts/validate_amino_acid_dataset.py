"""Offline validation of the shared amino-acid dataset (no Cantera, no pyCHNOSZ).

Checks the committed inputs/amino_acids_species.csv + inputs/amino_acids_gibbs_seed.csv:
- the dataset loads via load_species_metadata (no duplicate keys etc.);
- all amino-acid formulas parse and contain only C/H/N/O;
- every species has a finite molar volume and a full G(T) column over the grid;
- NASA9 fits cleanly for all 22 species, reporting max RMSE and flagging fit outliers.

Usage:
    python scripts/validate_amino_acid_dataset.py
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from config_io import load_species_metadata  # noqa: E402
from formula_tools import parse_formula  # noqa: E402
from thermo_fit import fit_all_species, read_wide_gibbs_csv  # noqa: E402

ALLOWED_ELEMENTS = {"C", "H", "N", "O"}
DEFAULT_SPECIES_CSV = PROJECT_ROOT / "inputs" / "amino_acids_species.csv"
DEFAULT_SEED_CSV = PROJECT_ROOT / "inputs" / "amino_acids_gibbs_seed.csv"


def run_validation(species_csv: str | Path = DEFAULT_SPECIES_CSV,
                   seed_csv: str | Path = DEFAULT_SEED_CSV) -> dict:
    """Run all checks; return a report dict (raises only on file-load failure)."""
    species_df = load_species_metadata(species_csv)        # raises on dup keys etc.
    seed = read_wide_gibbs_csv(seed_csv)
    species_names = list(species_df["cantera_name"])
    amino_acids = species_df[species_df["product_class"] == "amino_acid"]

    errors: list[str] = []

    # 1. Formulas: amino acids must be C/H/N/O only.
    bad_formula = {}
    for _, r in amino_acids.iterrows():
        elements = set(parse_formula(r["formula"]).keys())
        extra = elements - ALLOWED_ELEMENTS
        if extra:
            bad_formula[r["cantera_name"]] = sorted(extra)
    if bad_formula:
        errors.append(f"Non-C/H/N/O formulas: {bad_formula}")

    # 2. Finite molar volumes for every species.
    bad_vol = species_df.loc[~np.isfinite(species_df["molar_volume_cm3_mol"]), "cantera_name"].tolist()
    if bad_vol:
        errors.append(f"Missing/non-finite molar volumes: {bad_vol}")

    # 3. Both files cover all species; the seed has a full G(T) column (no gaps).
    seed_cols = [c for c in seed.columns if c != seed.columns[0]]
    missing_in_seed = [s for s in species_names if s not in seed_cols]
    if missing_in_seed:
        errors.append(f"Species missing a G(T) column in the seed: {missing_in_seed}")
    gappy = [s for s in species_names if s in seed_cols and seed[s].isna().any()]
    if gappy:
        errors.append(f"Species with gaps in their G(T) column: {gappy}")

    # 4. NASA9 fits cleanly for all species; report RMSE + flag outliers (>10x median).
    fit_report = {}
    fit_outliers = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        _, diag = fit_all_species(
            gibbs_wide_csv=seed_csv,
            coefficients_csv=tmp / "coeffs.csv",
            diagnostics_csv=tmp / "diag.csv",
            make_plots=False)
    rmse = dict(zip(diag["cantera_name"], diag["rmse_J_mol"]))
    fit_report = {k: float(v) for k, v in rmse.items()}
    median_rmse = float(np.median(list(rmse.values()))) if rmse else 0.0
    for name, val in rmse.items():
        if median_rmse > 0 and val > 10 * median_rmse:
            fit_outliers.append(name)

    return {
        "n_species": len(species_df),
        "n_amino_acids": len(amino_acids),
        "errors": errors,
        "ok": not errors,
        "fit_rmse_J_mol": fit_report,
        "max_rmse_J_mol": (max(fit_report.values()) if fit_report else None),
        "median_rmse_J_mol": median_rmse,
        "fit_outliers": fit_outliers,
    }


def main() -> int:
    report = run_validation()
    print(f"Species: {report['n_species']} ({report['n_amino_acids']} amino acids)")
    print(f"NASA9 fit RMSE — median {report['median_rmse_J_mol']:.1f} J/mol, "
          f"max {report['max_rmse_J_mol']:.1f} J/mol")
    worst = sorted(report["fit_rmse_J_mol"].items(), key=lambda kv: -kv[1])[:5]
    print("Highest-RMSE species:")
    for name, val in worst:
        print(f"  {name:<20} {val:>10.1f} J/mol")
    if report["fit_outliers"]:
        print(f"Fit outliers (>10x median RMSE): {report['fit_outliers']}")
    if report["errors"]:
        print("\nVALIDATION FAILED:")
        for e in report["errors"]:
            print(f"  - {e}")
        return 1
    print("\nValidation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
