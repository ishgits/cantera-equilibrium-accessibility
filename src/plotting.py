"""Plotting helpers for equilibrium accessibility results."""
from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional

import numpy as np
import pandas as pd


def _resolve_x_limits(
    values: pd.Series,
    formation_threshold: float,
    x_axis_max_mode: Literal["auto", "fixed"] = "auto",
    x_axis_max: float | None = None,
    x_axis_padding_factor: float = 30.0,
    x_axis_floor_factor: float = 5.0,
) -> tuple[float, float]:
    positive_vals = pd.to_numeric(values, errors="coerce")
    positive_vals = positive_vals[positive_vals > 0]
    if positive_vals.empty:
        x_floor = formation_threshold / 100.0
    else:
        x_floor = min(float(positive_vals.min()), float(formation_threshold)) / float(x_axis_floor_factor)
    x_floor = max(float(x_floor), 1e-300)

    observed_max = max(float(positive_vals.max()) if not positive_vals.empty else x_floor, float(formation_threshold))
    if x_axis_max_mode == "fixed":
        if x_axis_max is None:
            raise ValueError("x_axis_max must be supplied when x_axis_max_mode='fixed'.")
        x_right = float(x_axis_max)
    elif x_axis_max_mode == "auto":
        x_right = observed_max * float(x_axis_padding_factor)
    else:
        raise ValueError("x_axis_max_mode must be 'auto' or 'fixed'.")

    if x_right <= x_floor:
        x_right = x_floor * 10.0
    return x_floor, x_right


def plot_combined_accessibility_barchart(
    formation_df: pd.DataFrame,
    scenario: str,
    temperature_C: float,
    output_path: str | Path,
    formation_threshold: float = 1e-12,
    value_col: str = "X_eq",
    title: Optional[str] = None,
    accessible_color: str = "#2CA02C",
    below_threshold_color: str = "#CC3333",
    x_axis_max_mode: Literal["auto", "fixed"] = "auto",
    x_axis_max: float | None = None,
    x_axis_padding_factor: float = 30.0,
    x_axis_floor_factor: float = 5.0,
    label_padding_factor: float = 1.15,
    figure_width: float = 8.8,
    save_png: bool = True,
) -> None:
    """Horizontal bar chart showing accessibility and amount in one figure.

    All target products for the selected scenario/temperature are included.
    Products with ``X_eq >= formation_threshold`` are drawn as filled bars.
    Products below threshold are drawn as hollow outline bars. The dashed line
    marks the selected reporting threshold.

    A display floor is used only to show zero/tiny values on the log-scaled x
    axis; numeric labels always report the actual values from ``value_col``.
    """
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    import textwrap
    from matplotlib.lines import Line2D

    required = {"scenario", "T_C", "target_product", value_col}
    missing = required - set(formation_df.columns)
    if missing:
        raise ValueError(f"formation_df is missing required columns: {sorted(missing)}")

    df = formation_df[
        (formation_df["scenario"] == scenario)
        & (pd.to_numeric(formation_df["T_C"], errors="coerce") == float(temperature_C))
    ].copy()

    if df.empty:
        print(
            f"  [combined_barchart] No data for scenario={scenario!r} "
            f"T={temperature_C} °C — skipping."
        )
        return

    df["_val"] = pd.to_numeric(df[value_col], errors="coerce").fillna(0.0)
    df["_above_threshold"] = df["_val"] >= float(formation_threshold)
    df = df.sort_values("_val", ascending=False).reset_index(drop=True)

    x_floor, x_right = _resolve_x_limits(
        values=df["_val"],
        formation_threshold=float(formation_threshold),
        x_axis_max_mode=x_axis_max_mode,
        x_axis_max=x_axis_max,
        x_axis_padding_factor=x_axis_padding_factor,
        x_axis_floor_factor=x_axis_floor_factor,
    )

    # Plot zeros/tiny missing values as a small visible sliver on the log axis.
    df["_plot_val"] = df["_val"].where(df["_val"] > 0, x_floor * 1.2)
    df["_plot_val"] = df["_plot_val"].clip(lower=x_floor * 1.2)

    n = len(df)
    fig_height = max(3.0, 0.45 * n + 1.4)
    fig, ax = plt.subplots(figsize=(figure_width, fig_height), dpi=600, constrained_layout=True)

    y_pos = np.arange(n)
    bar_height = 0.62

    for y, (_, row) in zip(y_pos, df.iterrows()):
        bar_end = float(row["_plot_val"])
        width = max(bar_end - x_floor, x_floor * 0.2)

        if bool(row["_above_threshold"]):
            ax.barh(
                y,
                width,
                left=x_floor,
                height=bar_height,
                color=accessible_color,
                edgecolor="none",
                zorder=3,
            )
        else:
            ax.barh(
                y,
                width,
                left=x_floor,
                height=bar_height,
                facecolor="none",
                edgecolor=below_threshold_color,
                linewidth=1.2,
                zorder=3,
            )

    ax.axvline(
        formation_threshold,
        color=below_threshold_color,
        linestyle="--",
        linewidth=1.0,
        zorder=2,
    )

    ax.set_xlim(left=x_floor, right=x_right)
    ax.set_xscale("log")
    ax.set_xlabel("Equilibrium mole fraction  ($X_{eq}$)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Target product", fontsize=14, fontweight="bold")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["target_product"].tolist(), fontsize=12)
    ax.invert_yaxis()
    ax.xaxis.set_major_formatter(mtick.LogFormatterSciNotation())
    ax.grid(axis="x", which="both", linestyle=":", linewidth=0.5, alpha=0.6, zorder=1)
    ax.tick_params(axis="x", labelsize=12)

    for y, (_, row) in zip(y_pos, df.iterrows()):
        label_x = max(float(row["_plot_val"]) * float(label_padding_factor), x_floor * 1.35)
        # Keep labels from being clipped when users choose a tight fixed x-axis.
        label_x = min(label_x, x_right / 1.8)
        ax.text(
            label_x,
            y,
            f"{float(row['_val']):.2e}",
            va="center",
            ha="left",
            fontsize=10,
            color="#333333",
            clip_on=False,
        )

    threshold_handle = Line2D(
        [0], [1],
        linestyle="--",
        color=below_threshold_color,
        linewidth=1.0,
        label=f"Threshold ({formation_threshold:.0e})",
    )
    ax.legend(handles=[threshold_handle], fontsize=12, loc="lower right")

    if title is None:
        title = f"{scenario}  |  {temperature_C:g} °C"
    title = "\n".join(
        "\n".join(textwrap.wrap(line, width=72))
        for line in str(title).split("\n")
    )
    ax.set_title(title, fontsize=12, pad=8)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", pad_inches=0.18)
    if save_png and output_path.suffix.lower() == ".pdf":
        fig.savefig(
            output_path.with_suffix(".png"),
            dpi=600,
            bbox_inches="tight",
            pad_inches=0.18,
        )
    plt.close(fig)
    print(f"  [combined_barchart] Saved → {output_path.name}")
