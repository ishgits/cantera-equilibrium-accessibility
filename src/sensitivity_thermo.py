"""Study-local thermodynamic variants for the sensitivity layer.

Two pieces:

- :func:`augment_species_metadata_with_variants` — for each Gibbs-offset (ΔG)
  pseudo-species, create a study-local species row inheriting the base species'
  formula, state, molar volume, and product class with a new key/name (what the
  Cantera YAML writer and mole reconstruction need). With an empty offsets table
  (e.g. an inventory-only run) it is a no-op pass-through, so the base workflow is
  untouched.
- :func:`shift_coeffs_by_gibbs` — the exact analytic NASA9 coefficient shift
  (``a7 += ΔG·1000/R_GAS``, no refit; review §3.2) that gives a variant its
  thermodynamics.

Nothing here imports or runs Cantera.
"""
from __future__ import annotations

import pandas as pd

from thermo_fit import R_GAS


def shift_coeffs_by_gibbs(coeff_df: pd.DataFrame, base_name: str, variant_name: str,
                          dG_kJ_mol: float) -> pd.DataFrame:
    """Return new low+high NASA9 rows for a Gibbs-shifted pseudo-species.

    Exact analytic shift (no refit, review §3.2)::

        a7 += dG_kJ_mol * 1000 / R_GAS      # on BOTH the low and high segments

    Every other coefficient and the T-ranges are copied unchanged. ``a7`` enters
    only ``H/(RT)`` (as ``a7/T`` → a constant ``R·a7`` in H) and never ``S/R``, so
    shifting it by ``ΔG/R`` adds exactly ``ΔG`` to ``G = H − T·S`` at all
    temperatures, leaving ``Cp(T)`` and ``S(T)`` untouched. The input frame is not
    mutated (``base_name`` rows stay intact).
    """
    base = coeff_df[coeff_df["cantera_name"] == base_name]
    if base.empty:
        raise KeyError(f"No base NASA9 coefficients for {base_name!r} to shift.")
    delta_a7 = float(dG_kJ_mol) * 1000.0 / R_GAS
    out = base.copy()
    out["cantera_name"] = variant_name
    out["a7"] = out["a7"] + delta_a7
    return out.reset_index(drop=True)


def require_gibbs_seed(config: dict) -> str:
    """Return the configured Gibbs seed path, or raise if it is not set.

    NASA9 fitting (Phase 2) cannot proceed without a Gibbs seed CSV, so a missing
    ``species_files.gibbs_seed_wide_csv`` must fail with a clear message rather than
    a downstream ``TypeError``.
    """
    rel = (config.get("species_files", {}) or {}).get("gibbs_seed_wide_csv")
    if not rel or not str(rel).strip():
        from sensitivity_design import StudyConfigError  # local: avoid import cycle
        raise StudyConfigError(
            "species_files.gibbs_seed_wide_csv is not set, but NASA9 coefficient "
            "fitting needs a Gibbs seed CSV. Add it to the study config."
        )
    return str(rel)


def augment_species_metadata_with_variants(species_df: pd.DataFrame,
                                           offsets_df: pd.DataFrame) -> pd.DataFrame:
    """Return ``species_df`` plus one pseudo-species row per ΔG variant.

    Each variant row copies its base species row (formula, state, molar volume,
    role, product class, chnosz name preserved) and changes only ``species_key``,
    ``cantera_name`` and ``notes``. The original frame is never mutated.

    ``offsets_df`` columns: ``thermo_variant_id, base_species, variant_species,
    deltaG_offset_kJ_mol`` (the table written by Phase 1's
    ``sensitivity_design.build_thermo_offsets_table``). Empty ⇒ returns a copy of
    ``species_df`` unchanged.
    """
    if offsets_df is None or offsets_df.empty:
        return species_df.copy()

    by_name = species_df.set_index("cantera_name", drop=False)
    new_rows = []
    existing = set(species_df["cantera_name"])
    for _, off in offsets_df.iterrows():
        base_name = str(off["base_species"])
        variant_name = str(off["variant_species"])
        if variant_name in existing:
            continue  # already present; don't duplicate
        if base_name not in by_name.index:
            raise KeyError(
                f"ΔG base species {base_name!r} is not in the species table; "
                "cannot build its variant metadata."
            )
        row = by_name.loc[base_name].to_dict()
        row["cantera_name"] = variant_name
        if "species_key" in row:
            row["species_key"] = f"{row.get('species_key', base_name)}__{off['thermo_variant_id']}"
        row["notes"] = (
            f"ΔG sensitivity variant of {base_name}: "
            f"{float(off['deltaG_offset_kJ_mol']):+g} kJ/mol Gibbs offset"
        )
        new_rows.append(row)

    if not new_rows:
        return species_df.copy()
    return pd.concat([species_df, pd.DataFrame(new_rows)], ignore_index=True)
