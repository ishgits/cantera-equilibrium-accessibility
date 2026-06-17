"""CHNOSZ Gibbs-energy cache tools.

The cache is long/tidy: one row per species-temperature-state value.
"""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

CACHE_COLUMNS = [
    "species_key",
    "cantera_name",
    "chnosz_name",
    "formula",
    "state",
    "T_C",
    "T_K",
    "G_J_mol",
    "source",
    "extraction_date_utc",
    "notes",
]


def empty_cache() -> pd.DataFrame:
    return pd.DataFrame(columns=CACHE_COLUMNS)


def load_cache(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        return empty_cache()
    df = pd.read_csv(path)
    for col in CACHE_COLUMNS:
        if col not in df.columns:
            df[col] = np.nan
    df = df[CACHE_COLUMNS].copy()
    df["T_C"] = pd.to_numeric(df["T_C"], errors="coerce")
    df["T_K"] = pd.to_numeric(df["T_K"], errors="coerce")
    df["G_J_mol"] = pd.to_numeric(df["G_J_mol"], errors="coerce")
    return df


def save_cache(cache_df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = cache_df.copy()
    out = out[CACHE_COLUMNS]
    out = out.sort_values(["species_key", "state", "T_K"]).reset_index(drop=True)
    out.to_csv(path, index=False)


def requested_grid(species_df: pd.DataFrame, temperatures_C: Sequence[float]) -> pd.DataFrame:
    rows = []
    for _, sp in species_df.iterrows():
        for t_c in temperatures_C:
            rows.append({
                "species_key": sp["species_key"],
                "cantera_name": sp["cantera_name"],
                "chnosz_name": sp["chnosz_name"],
                "formula": sp["formula"],
                "state": sp["state"],
                "T_C": float(t_c),
                "T_K": 273.16 if abs(float(t_c)) < 1e-9 else float(t_c) + 273.15,
            })
    return pd.DataFrame(rows)


def find_missing_rows(cache_df: pd.DataFrame, species_df: pd.DataFrame, temperatures_C: Sequence[float]) -> pd.DataFrame:
    """Return species-temperature rows missing from the current cache."""
    req = requested_grid(species_df, temperatures_C)
    if cache_df.empty:
        return req
    cached = cache_df.dropna(subset=["species_key", "state", "T_C", "G_J_mol"]).copy()
    cached["T_K_round"] = cached["T_K"].round(6)
    req["T_K_round"] = req["T_K"].round(6)
    keys = ["species_key", "state", "T_K_round"]
    merged = req.merge(cached[keys + ["G_J_mol"]], on=keys, how="left", indicator=True)
    missing = merged[merged["_merge"] == "left_only"].drop(columns=["_merge", "G_J_mol", "T_K_round"])
    return missing.reset_index(drop=True)


def _extract_single_species_with_pychNOSZ(chnosz_name: str, state: str, temperatures_C: Sequence[float], exceed_Ttr: bool) -> pd.DataFrame:
    """Call pyCHNOSZ.subcrt and return its output dataframe for one species."""
    try:
        import pyCHNOSZ as pcz
    except ImportError as exc:
        raise ImportError(
            "pyCHNOSZ is not installed in this environment. Install it or pre-populate "
            "data/raw/chnosz_gibbs_cache.csv before running extraction."
        ) from exc

    kwargs = dict(property="G", T=list(temperatures_C), exceed_Ttr=exceed_Ttr)
    state_clean = str(state).strip()
    
    state_clean = str(state).strip().lower()
    state_map = {
        "aqueous": "aq",
        "liquid": "liq",
    }

    if state_clean and state_clean not in {"default", "nan"}:
        kwargs["state"] = state_map.get(state_clean, state_clean)

    data = pcz.subcrt(chnosz_name, **kwargs)
    if data is None or not getattr(data, "out", None):
        raise ValueError(f"No CHNOSZ data returned for {chnosz_name!r} state={state!r}.")

    if chnosz_name in data.out:
        return data.out[chnosz_name].copy()
    # Fallback: take the first returned table.
    return next(iter(data.out.values())).copy()


def extract_missing_rows(
    missing_df: pd.DataFrame,
    exceed_Ttr: bool = True,
) -> pd.DataFrame:
    """Extract missing CHNOSZ rows using pyCHNOSZ.

    This groups missing requests by species so pyCHNOSZ is called once per species,
    not once per temperature.
    """
    if missing_df.empty:
        return empty_cache()

    rows = []
    now = datetime.now(timezone.utc).isoformat()
    group_cols = ["species_key", "cantera_name", "chnosz_name", "formula", "state"]
    for keys, group in missing_df.groupby(group_cols, dropna=False):
        species_key, cantera_name, chnosz_name, formula, state = keys
        temps = sorted(group["T_C"].astype(float).unique().tolist())
        out = _extract_single_species_with_pychNOSZ(chnosz_name, state, temps, exceed_Ttr=exceed_Ttr)
        # pyCHNOSZ output usually has T and G columns. Be forgiving about names.
        cols_lower = {c.lower().strip(): c for c in out.columns}
        t_col = cols_lower.get("t") or cols_lower.get("temp") or cols_lower.get("temperature")
        g_col = cols_lower.get("g") or cols_lower.get("gibbs") or cols_lower.get("g_j_mol")
        if t_col is None or g_col is None:
            raise ValueError(f"Could not identify T/G columns in CHNOSZ output for {chnosz_name}: {out.columns.tolist()}")
        for _, r in out.iterrows():
            t_c = float(r[t_col])
            rows.append({
                "species_key": species_key,
                "cantera_name": cantera_name,
                "chnosz_name": chnosz_name,
                "formula": formula,
                "state": state,
                "T_C": t_c,
                "T_K": 273.16 if abs(t_c) < 1e-9 else t_c + 273.15,
                "G_J_mol": float(r[g_col]),
                "source": "CHNOSZ/pyCHNOSZ",
                "extraction_date_utc": now,
                "notes": "",
            })
    return pd.DataFrame(rows, columns=CACHE_COLUMNS)


def update_cache_with_missing(
    species_df: pd.DataFrame,
    temperatures_C: Sequence[float],
    cache_path: str | Path,
    force_reextract: bool = False,
    exceed_Ttr: bool = True,
) -> pd.DataFrame:
    """Update CHNOSZ cache by extracting only missing species-temperature rows."""
    cache = empty_cache() if force_reextract else load_cache(cache_path)
    missing = find_missing_rows(cache, species_df, temperatures_C)
    if missing.empty:
        save_cache(cache, cache_path)
        return cache
    extracted = extract_missing_rows(missing, exceed_Ttr=exceed_Ttr)
    combined = extracted if cache.empty else pd.concat([cache, extracted], ignore_index=True)
    combined = combined.drop_duplicates(subset=["species_key", "state", "T_C"], keep="last")
    save_cache(combined, cache_path)
    return combined


def make_gibbs_wide(
    cache_df: pd.DataFrame,
    species_df: pd.DataFrame,
    temperatures_C: Sequence[float],
    output_path: str | Path,
    use_column: str = "cantera_name",
) -> pd.DataFrame:
    """Create wide Gibbs table for NASA9 fitting.

    Columns: T_K + one column per species.
    """
    req = requested_grid(species_df, temperatures_C)
    needed_keys = set(req["species_key"])
    df = cache_df[cache_df["species_key"].isin(needed_keys)].copy()
    df["T_K_round"] = df["T_K"].round(6)
    req["T_K_round"] = req["T_K"].round(6)
    merged = req[["species_key", "state", "T_K_round", "T_K", use_column]].merge(df[["species_key", "state", "T_K_round", "G_J_mol"]], on=["species_key", "state", "T_K_round"],how="left")
    
    if merged["G_J_mol"].isna().any():
        bad = merged[merged["G_J_mol"].isna()][["species_key", "T_K"]]
        raise ValueError(f"Missing Gibbs values for requested rows:\n{bad.to_string(index=False)}")
    wide = merged.pivot_table(index="T_K", columns=use_column, values="G_J_mol", aggfunc="first").reset_index()
    wide = wide.sort_values("T_K")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(output_path, index=False)
    return wide


def seed_cache_from_wide_csv(
    wide_csv: str | Path,
    species_df: pd.DataFrame,
    cache_path: str | Path,
    name_map: Optional[dict] = None,
    source: str = "seeded wide CSV",
) -> pd.DataFrame:
    """Seed/update the cache from an existing wide Gibbs-energy CSV.

    This is useful for validating the workflow with old datasets before pyCHNOSZ
    extraction is available. `name_map` maps wide CSV column names to Cantera names.
    """
    wide_csv = Path(wide_csv)
    raw = pd.read_csv(wide_csv, comment="#")
    # Drop empty/unnamed columns and rows without temperatures.
    raw = raw.dropna(axis=1, how="all")
    raw.columns = [str(c).strip() for c in raw.columns]
    # Detect temperature column.
    temp_col = None
    for c in raw.columns:
        lc = c.lower().replace(" ", "")
        if lc in {"t(k)", "t_k", "temperature(k)", "temperature_k"} or lc.startswith("t("):
            temp_col = c
            break
    if temp_col is None:
        # For user files with comment lines, try re-read using the third row as header.
        raw = pd.read_csv(wide_csv, skiprows=2)
        raw = raw.dropna(axis=1, how="all")
        raw.columns = [str(c).strip() for c in raw.columns]
        for c in raw.columns:
            lc = c.lower().replace(" ", "")
            if lc in {"t(k)", "t_k", "temperature(k)", "temperature_k"} or lc.startswith("t("):
                temp_col = c
                break
    if temp_col is None:
        raise ValueError(f"Could not find temperature column in {wide_csv}")

    name_map = name_map or {}
    species_by_cantera = species_df.set_index("cantera_name").to_dict(orient="index")
    rows = []
    now = datetime.now(timezone.utc).isoformat()
    for col in raw.columns:
        if col == temp_col or str(col).startswith("Unnamed"):
            continue
        cantera_name = name_map.get(col, col)
        if cantera_name not in species_by_cantera:
            continue
        meta = species_by_cantera[cantera_name]
        for _, r in raw.iterrows():
            if pd.isna(r[temp_col]) or pd.isna(r[col]):
                continue
            t_k = float(r[temp_col])
            rows.append({
                "species_key": meta["species_key"],
                "cantera_name": cantera_name,
                "chnosz_name": meta["chnosz_name"],
                "formula": meta["formula"],
                "state": meta["state"],
                "T_C": 0.0 if abs(t_k - 273.16) < 1e-6 else t_k - 273.15,
                "T_K": t_k,
                "G_J_mol": float(r[col]),
                "source": source,
                "extraction_date_utc": now,
                "notes": f"seeded from {wide_csv.name}",
            })
    if not rows:
        raise ValueError("No seedable species columns were found in the wide CSV.")
    cache = load_cache(cache_path)
    new_rows = pd.DataFrame(rows, columns=CACHE_COLUMNS)
    combined = new_rows if cache.empty else pd.concat([cache, new_rows], ignore_index=True)
    combined = combined.drop_duplicates(subset=["species_key", "state", "T_C"], keep="last")
    save_cache(combined, cache_path)
    return combined
