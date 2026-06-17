"""Diagnostics and notebook summary helpers for equilibrium accessibility workflows."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from mole_balance import run_group_columns


def _sort_cols(df: pd.DataFrame, candidates: Sequence[str]) -> list[str]:
    return [c for c in candidates if c in df.columns]


def _scenario_description_map(scenarios: Mapping[str, Any] | None) -> dict[str, str]:
    if not scenarios or "scenarios" not in scenarios:
        return {}
    out: dict[str, str] = {}
    for scenario_id, cfg in scenarios["scenarios"].items():
        desc = cfg.get("description", "") if isinstance(cfg, Mapping) else ""
        out[str(scenario_id)] = "" if desc is None else str(desc)
    return out




def _markdown_table(df: pd.DataFrame) -> str:
    """Render a small DataFrame as a Markdown table without optional dependencies."""
    if df.empty:
        return "_No rows._"
    table = df.copy()
    table = table.fillna("")
    cols = [str(c) for c in table.columns]
    lines = []
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("|" + "|".join(["---" for _ in cols]) + "|")
    for _, row in table.iterrows():
        vals = [str(row[c]).replace("\n", " ") for c in table.columns]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def _formation_merge_columns(left: pd.DataFrame, formation_df: pd.DataFrame) -> list[str]:
    base = run_group_columns(left)
    return [c for c in base if c in formation_df.columns]


def make_equilibrium_inspection_table(
    moles_long_df: pd.DataFrame,
    formation_df: pd.DataFrame | None = None,
    scenarios: Mapping[str, Any] | None = None,
    output_csv: str | Path | None = None,
) -> pd.DataFrame:
    """Create a readable species-level raw inspection table.

    This table is intentionally close to the raw Cantera output, but includes
    reconstructed equilibrium moles and the target-level formation call when a
    ``formation_df`` is supplied. It is meant for manual debugging and sanity
    checks rather than publication plotting.
    """
    required = {"scenario", "target_product", "T_C", "species", "X_eq", "initial_moles", "n_eq_mol"}
    missing = required - set(moles_long_df.columns)
    if missing:
        raise ValueError(f"moles_long_df is missing required columns: {sorted(missing)}")

    out = moles_long_df.copy()

    desc_map = _scenario_description_map(scenarios)
    if desc_map:
        out["scenario_description"] = out["scenario"].astype(str).map(desc_map).fillna("")

    if formation_df is not None and not formation_df.empty:
        merge_cols = _formation_merge_columns(out, formation_df)
        formation_cols = [
            c for c in [
                "X_eq", "log10_X_eq", "n_eq_mol", "log10_n_eq_mol",
                "formed_bool", "formation_call",
            ]
            if c in formation_df.columns and c not in merge_cols
        ]
        target_summary = formation_df[merge_cols + formation_cols].drop_duplicates(merge_cols)
        rename = {
            c: f"target_{c}"
            for c in formation_cols
            if c in {"X_eq", "log10_X_eq", "n_eq_mol", "log10_n_eq_mol"}
        }
        target_summary = target_summary.rename(columns=rename)
        out = out.merge(target_summary, on=merge_cols, how="left")

    preferred = [
        "scenario", "scenario_description", "model_mode", "yaml_file", "target_product",
        "T_C", "T_K", "P_Pa", "species", "X_initial", "X_eq", "initial_moles",
        "n_eq_mol", "n_total_eq_mol", "target_X_eq", "target_log10_X_eq",
        "target_n_eq_mol", "target_log10_n_eq_mol", "formed_bool", "formation_call",
        "solver_status", "error_message", "element_balance_relative_spread",
        "element_total_mole_estimates",
    ]
    ordered = [c for c in preferred if c in out.columns]
    extras = [c for c in out.columns if c not in ordered]
    out = out[ordered + extras]

    sort_cols = _sort_cols(out, ["scenario", "target_product", "T_C", "species"])
    if sort_cols:
        out = out.sort_values(sort_cols).reset_index(drop=True)

    if output_csv is not None:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(output_csv, index=False)
    return out


def reactant_depletion_diagnostics(
    moles_long_df: pd.DataFrame,
    starting_species: Sequence[str] | None = None,
    target_products: Sequence[str] | None = None,
    output_csv: str | Path | None = None,
) -> pd.DataFrame:
    """Calculate depletion for starting species in each run.

    ``depletion_fraction = (initial_moles - n_eq_mol) / initial_moles``.
    Positive values mean net consumption; negative values mean net production.

    This is a source/depletion diagnostic only. It does not identify a reaction
    pathway and should not be interpreted as a formal limiting reagent analysis.
    """
    required = {"species", "initial_moles", "n_eq_mol"}
    missing = required - set(moles_long_df.columns)
    if missing:
        raise ValueError(f"moles_long_df is missing required columns: {sorted(missing)}")

    df = moles_long_df.copy()
    target_set = set(target_products or [])

    if starting_species is None:
        mask = pd.to_numeric(df["initial_moles"], errors="coerce").fillna(0) > 0
        if target_set:
            mask &= ~df["species"].isin(target_set)
        else:
            mask &= df["species"] != df["target_product"]
        df = df[mask].copy()
    else:
        df = df[df["species"].isin(list(starting_species))].copy()

    df["n_consumed_mol"] = df["initial_moles"] - df["n_eq_mol"]
    df["depletion_fraction"] = np.where(
        df["initial_moles"] > 0,
        df["n_consumed_mol"] / df["initial_moles"],
        np.nan,
    )

    keep = [
        "scenario", "model_mode", "yaml_file", "target_product", "T_C", "T_K", "P_Pa",
        "species", "X_initial", "X_eq", "initial_moles", "n_eq_mol", "n_consumed_mol",
        "depletion_fraction", "element_balance_relative_spread", "solver_status",
        "error_message",
    ]

    out = df[[c for c in keep if c in df.columns]].sort_values(
        _sort_cols(df, ["scenario", "target_product", "T_C", "species"])
    ).reset_index(drop=True)

    if output_csv is not None:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(output_csv, index=False)
    return out


def summarize_reactant_depletion(
    moles_long_df: pd.DataFrame,
    formation_df: pd.DataFrame,
    starting_species: Sequence[str] | None = None,
    target_products: Sequence[str] | None = None,
    output_long_csv: str | Path | None = None,
    output_summary_csv: str | Path | None = None,
    significant_depletion_fraction: float = 1e-6,
    significant_consumption_mol: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return long and one-row-per-run depletion diagnostics.

    The summary reports two complementary diagnostics:

    - ``most_depleted_species``: largest fractional depletion.
    - ``most_consumed_species``: largest absolute mole consumption.

    For targets below the accessibility threshold, the depletion values are kept
    for debugging, but ``depletion_call`` avoids inferring a source species.
    """
    depletion_long = reactant_depletion_diagnostics(
        moles_long_df=moles_long_df,
        starting_species=starting_species,
        target_products=target_products,
        output_csv=output_long_csv,
    )

    if depletion_long.empty:
        summary = pd.DataFrame()
        if output_summary_csv is not None:
            output_summary_csv = Path(output_summary_csv)
            output_summary_csv.parent.mkdir(parents=True, exist_ok=True)
            summary.to_csv(output_summary_csv, index=False)
        return depletion_long, summary

    group_cols = run_group_columns(depletion_long)
    merge_cols = [c for c in group_cols if c in formation_df.columns]

    formation_keep = [
        c for c in [
            "X_eq", "log10_X_eq", "n_eq_mol", "log10_n_eq_mol", "formed_bool",
            "formation_call", "solver_status", "element_balance_relative_spread",
        ]
        if c in formation_df.columns and c not in merge_cols
    ]
    formation_view = formation_df[merge_cols + formation_keep].drop_duplicates(merge_cols)

    rows: list[dict[str, Any]] = []
    for key, group in depletion_long.groupby(group_cols, dropna=False):
        key_tuple = key if isinstance(key, tuple) else (key,)
        row = {col: val for col, val in zip(group_cols, key_tuple)}
        g = group.copy()
        valid_frac = g.dropna(subset=["depletion_fraction"])
        valid_abs = g.dropna(subset=["n_consumed_mol"])

        if valid_frac.empty:
            row.update({
                "most_depleted_species": "",
                "max_depletion_fraction": np.nan,
            })
        else:
            r_frac = valid_frac.loc[valid_frac["depletion_fraction"].idxmax()]
            row.update({
                "most_depleted_species": r_frac["species"],
                "max_depletion_fraction": float(r_frac["depletion_fraction"]),
            })

        if valid_abs.empty:
            row.update({
                "most_consumed_species": "",
                "max_consumed_mol": np.nan,
            })
        else:
            r_abs = valid_abs.loc[valid_abs["n_consumed_mol"].idxmax()]
            row.update({
                "most_consumed_species": r_abs["species"],
                "max_consumed_mol": float(r_abs["n_consumed_mol"]),
            })

        rows.append(row)

    summary = pd.DataFrame(rows)
    if not formation_view.empty:
        summary = summary.merge(formation_view, on=merge_cols, how="left")

    # Rename target-level quantities so the table is unambiguous.
    summary = summary.rename(columns={
        "X_eq": "target_X_eq",
        "log10_X_eq": "target_log10_X_eq",
        "n_eq_mol": "target_n_eq_mol",
        "log10_n_eq_mol": "target_log10_n_eq_mol",
    })

    def _call(row: pd.Series) -> str:
        solver_status = str(row.get("solver_status", ""))
        formation_call = str(row.get("formation_call", ""))
        formed = bool(row.get("formed_bool", False)) if pd.notna(row.get("formed_bool", np.nan)) else False
        max_frac = row.get("max_depletion_fraction", np.nan)
        max_mol = row.get("max_consumed_mol", np.nan)
        if solver_status and solver_status != "ok":
            return "solver_failed_no_source_inferred"
        if formation_call in {"below_threshold", "target_not_present_in_yaml"} or not formed:
            return "target_below_threshold_no_source_inferred"
        frac_ok = pd.notna(max_frac) and float(max_frac) >= float(significant_depletion_fraction)
        mol_ok = pd.notna(max_mol) and float(max_mol) > float(significant_consumption_mol)
        if frac_ok or mol_ok:
            return "dominant_depletion_detected"
        return "accessible_product_no_significant_depletion_detected"

    if not summary.empty:
        summary["depletion_call"] = summary.apply(_call, axis=1)
        preferred = [
            "scenario", "model_mode", "yaml_file", "target_product", "T_C",
            "target_X_eq", "target_log10_X_eq", "target_n_eq_mol",
            "formed_bool", "formation_call", "most_depleted_species",
            "max_depletion_fraction", "most_consumed_species", "max_consumed_mol",
            "depletion_call", "solver_status", "element_balance_relative_spread",
        ]
        ordered = [c for c in preferred if c in summary.columns]
        extras = [c for c in summary.columns if c not in ordered]
        summary = summary[ordered + extras]
        summary = summary.sort_values(_sort_cols(summary, ["scenario", "target_product", "T_C"])).reset_index(drop=True)

    if output_summary_csv is not None:
        output_summary_csv = Path(output_summary_csv)
        output_summary_csv.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(output_summary_csv, index=False)

    return depletion_long, summary


def write_run_summary(
    output_md: str | Path,
    project_name: str,
    species_file: str | Path,
    scenario_file: str | Path,
    scenarios: Mapping[str, Any],
    target_products: Sequence[str],
    equilibrium_temperatures_C: Sequence[float],
    thermo_fit_temperatures_C: Sequence[float],
    pressure_Pa: float,
    formation_x_threshold: float,
    significant_x_threshold: float,
    formation_n_threshold_mol: float,
    formation_df: pd.DataFrame,
    depletion_summary_df: pd.DataFrame | None = None,
    output_paths: Mapping[str, str | Path] | None = None,
) -> Path:
    """Write a Markdown summary of the completed notebook run."""
    output_md = Path(output_md)
    output_md.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    lines.append(f"# Run summary — {project_name}\n")
    lines.append("## Modeled conditions\n")
    lines.append(f"- Species file: `{Path(species_file).as_posix()}`")
    lines.append(f"- Scenario file: `{Path(scenario_file).as_posix()}`")
    lines.append(f"- Pressure: `{float(pressure_Pa):.6g} Pa` (`{float(pressure_Pa)/1e5:.3g} bar`)")
    lines.append(f"- Equilibrium temperatures (°C): `{list(equilibrium_temperatures_C)}`")
    lines.append(f"- Thermo fitting temperatures (°C): `{list(thermo_fit_temperatures_C)}`")
    lines.append(f"- Formation threshold: `X_eq >= {float(formation_x_threshold):.3e}`")
    lines.append(f"- Significant threshold: `X_eq >= {float(significant_x_threshold):.3e}`")
    lines.append(f"- Minimum mole threshold: `n_eq_mol >= {float(formation_n_threshold_mol):.3e}`")
    lines.append(f"- Target products modeled: `{len(target_products)}`\n")

    lines.append("## Starting inventory\n")
    for scenario_id, cfg in scenarios.get("scenarios", {}).items():
        desc = cfg.get("description", "") if isinstance(cfg, Mapping) else ""
        title = f"### {scenario_id}"
        if desc:
            title += f" — {desc}"
        lines.append(title)
        lines.append("")
        lines.append("| Species | Initial moles |")
        lines.append("|---|---:|")
        for sp, mol in cfg.get("initial_moles", {}).items():
            lines.append(f"| `{sp}` | {float(mol):.8g} |")
        lines.append("")

    lines.append("## Equilibrium accessibility summary\n")
    if formation_df.empty:
        lines.append("No formation summary rows were generated.\n")
    else:
        counts = formation_df["formation_call"].value_counts(dropna=False).to_dict()
        total = len(formation_df)
        lines.append(f"- Total target runs summarized: `{total}`")
        for label in ["significant", "trace", "below_threshold", "solver_failed", "target_not_present_in_yaml"]:
            if label in counts:
                lines.append(f"- {label}: `{counts[label]}`")
        lines.append("")
        cols = [c for c in ["scenario", "target_product", "T_C", "X_eq", "n_eq_mol", "formation_call", "solver_status"] if c in formation_df.columns]
        preview = formation_df[cols].copy()
        if "X_eq" in preview.columns:
            preview["X_eq"] = preview["X_eq"].map(lambda x: "" if pd.isna(x) else f"{float(x):.3e}")
        if "n_eq_mol" in preview.columns:
            preview["n_eq_mol"] = preview["n_eq_mol"].map(lambda x: "" if pd.isna(x) else f"{float(x):.3e}")
        lines.append(_markdown_table(preview))
        lines.append("")

    if depletion_summary_df is not None and not depletion_summary_df.empty:
        lines.append("## Reactant depletion diagnostic\n")
        cols = [c for c in [
            "scenario", "target_product", "T_C", "formation_call", "most_depleted_species",
            "max_depletion_fraction", "most_consumed_species", "max_consumed_mol", "depletion_call",
        ] if c in depletion_summary_df.columns]
        preview = depletion_summary_df[cols].copy()
        if "max_depletion_fraction" in preview.columns:
            preview["max_depletion_fraction"] = preview["max_depletion_fraction"].map(lambda x: "" if pd.isna(x) else f"{float(x):.3e}")
        if "max_consumed_mol" in preview.columns:
            preview["max_consumed_mol"] = preview["max_consumed_mol"].map(lambda x: "" if pd.isna(x) else f"{float(x):.3e}")
        lines.append(_markdown_table(preview))
        lines.append("")

    lines.append("## Saved outputs\n")
    if output_paths:
        for label, path in output_paths.items():
            lines.append(f"- {label}: `{Path(path).as_posix()}`")
    else:
        lines.append("No output path map was supplied.")
    lines.append("")

    lines.append("## Interpretation note\n")
    lines.append(
        "These outputs report equilibrium accessibility for species explicitly included "
        "in each single-product YAML. They do not report kinetic rates, reaction pathways, "
        "or formal percent yields. The depletion diagnostic identifies which starting "
        "species changed most during equilibrium redistribution."
    )
    lines.append("")

    output_md.write_text("\n".join(lines), encoding="utf-8")
    return output_md
