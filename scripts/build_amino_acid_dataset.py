"""Build the shared amino-acid thermodynamic dataset from CHNOSZ (real data only).

Produces, for the breadth campaign (one study per amino acid):

    inputs/amino_acids_species.csv      # 22 rows: 4 feedstock + 18 C/H/N/O amino acids
    inputs/amino_acids_gibbs_seed.csv   # wide G(T): T_K + one column per Cantera name

Scope: the 18 standard amino acids whose formula is C/H/N/O only (the 20 minus
cysteine and methionine — sulfur is not in the Titan feedstock). Formulas and
standard partial molar volumes are pulled from CHNOSZ; G(T) is extracted over the
0–370 C fit grid via the existing ``chnosz_cache`` path. Nothing is fabricated.

Requires pyCHNOSZ/pychnosz (R + the CHNOSZ package). Run:
    R_HOME=$(R RHOME) python scripts/build_amino_acid_dataset.py
"""
from __future__ import annotations

import contextlib
import io
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from chnosz_cache import (  # noqa: E402
    _import_pychnosz, load_cache, make_gibbs_wide, update_cache_with_missing,
)

# 0..370 C in 10 C steps — the same grid used to fit the alanine MVP.
GRID_C = [float(T) for T in range(0, 371, 10)]

SPECIES_COLUMNS = ["species_key", "cantera_name", "chnosz_name", "formula", "state",
                   "molar_volume_cm3_mol", "role", "notes", "product_class"]

# 18 C/H/N/O amino acids: (Cantera name, CHNOSZ name). Neutral aqueous species.
AMINO_ACIDS = [
    ("Glycine(aq)", "glycine"), ("Alanine(aq)", "alanine"), ("Serine(aq)", "serine"),
    ("Proline(aq)", "proline"), ("Valine(aq)", "valine"), ("Threonine(aq)", "threonine"),
    ("Leucine(aq)", "leucine"), ("Isoleucine(aq)", "isoleucine"),
    ("Asparagine(aq)", "asparagine"), ("AsparticAcid(aq)", "aspartic acid"),
    ("Glutamine(aq)", "glutamine"), ("Lysine(aq)", "lysine"),
    ("GlutamicAcid(aq)", "glutamic acid"), ("Arginine(aq)", "arginine"),
    ("Histidine(aq)", "histidine"), ("Phenylalanine(aq)", "phenylalanine"),
    ("Tyrosine(aq)", "tyrosine"), ("Tryptophan(aq)", "tryptophan"),
]

# Feedstock/solvent reused verbatim from inputs/species_example.csv (already CHNOSZ-based).
REACTANTS = [
    dict(species_key="h2o", cantera_name="H2O(l)", chnosz_name="water", formula="H2O",
         state="liq", molar_volume_cm3_mol=18.015, role="solvent",
         notes="liquid water basis species", product_class=""),
    dict(species_key="hcn", cantera_name="HCN(aq)", chnosz_name="HCN", formula="CHN",
         state="aq", molar_volume_cm3_mol=23.9, role="reactant",
         notes="hydrogen cyanide starting inventory", product_class=""),
    dict(species_key="c2h2", cantera_name="C2H2(aq)", chnosz_name="ethyne", formula="C2H2",
         state="aq", molar_volume_cm3_mol=26.038, role="reactant",
         notes="acetylene starting inventory", product_class=""),
    dict(species_key="nh3", cantera_name="NH3(aq)", chnosz_name="NH3", formula="H3N",
         state="aq", molar_volume_cm3_mol=24.5, role="reactant",
         notes="ammonia starting inventory", product_class=""),
]


def chnosz_info(name: str) -> dict:
    """Return CHNOSZ {formula, V, Z, name} for the neutral aqueous species `name`."""
    pcz = _import_pychnosz()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        idx = pcz.info(name, state="aq")
        df = pcz.info(idx)
    r = df.iloc[0]
    return {"formula": str(r["formula"]), "V": float(r["V"]),
            "Z": int(r["Z"]), "name": str(r["name"])}


def build_species_table() -> pd.DataFrame:
    rows = [dict(r) for r in REACTANTS]
    for cantera_name, chnosz_name in AMINO_ACIDS:
        info = chnosz_info(chnosz_name)
        if info["Z"] != 0:
            raise ValueError(f"{chnosz_name}: CHNOSZ species is charged (Z={info['Z']}); "
                             "expected the neutral form.")
        rows.append({
            "species_key": cantera_name.split("(")[0].lower(),
            "cantera_name": cantera_name,
            "chnosz_name": chnosz_name,
            "formula": info["formula"],
            "state": "aq",
            "molar_volume_cm3_mol": info["V"],
            "role": "product",
            "notes": f"CHNOSZ '{chnosz_name}'; standard partial molar volume from CHNOSZ",
            "product_class": "amino_acid",
        })
    return pd.DataFrame(rows, columns=SPECIES_COLUMNS)


def main() -> int:
    inputs = PROJECT_ROOT / "inputs"
    species_path = inputs / "amino_acids_species.csv"
    seed_path = inputs / "amino_acids_gibbs_seed.csv"
    cache_path = PROJECT_ROOT / "data" / "raw" / "amino_acids_gibbs_cache.csv"

    try:
        _import_pychnosz()
    except ImportError as exc:
        print(f"Cannot build dataset: {exc}", file=sys.stderr)
        print("Needed species still to extract: all 22 (4 feedstock + 18 amino acids).",
              file=sys.stderr)
        return 1

    print("Querying CHNOSZ for formulas + molar volumes ...")
    species_df = build_species_table()
    species_df.to_csv(species_path, index=False)
    print(f"Wrote {species_path} ({len(species_df)} species).")

    print(f"Extracting G(T) over {len(GRID_C)} temperatures (0..370 C) for all "
          f"{len(species_df)} species ...")
    update_cache_with_missing(species_df, GRID_C, cache_path)
    cache = load_cache(cache_path)
    make_gibbs_wide(cache, species_df, GRID_C, seed_path)
    wide = pd.read_csv(seed_path)
    print(f"Wrote {seed_path} ({wide.shape[0]} rows x {wide.shape[1]} cols "
          f"= T_K + {wide.shape[1] - 1} species).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
