"""Config-driven sensitivity figures (Phase 3).

Reads the ``plots:`` block of ``study_config.yaml`` (every key optional, sensible
defaults) and renders matplotlib-only figures into ``studies/<id>/figures/`` in each
requested format. No seaborn. Failed/missing/inaccessible cases stay **visible**
(review §4, §14.6).

Each heatmap landscape is rendered at **two scales** (config key ``scales``,
default ``[fixed, auto]``): a *fixed* version using the configured ``vmin``/``vmax``
(default −12..0, showing the accessible-vs-inaccessible cliff) and an *autoscaled*
version whose limits come from the robust 2nd/98th percentiles of the finite
``log10_X_eq`` (revealing the within-region gradient). The primary accessibility
boundary is drawn from ``formed_bool`` (not from a fragile ``log10_X_eq`` contour
level), and inaccessible cells are hatched so even a one-cell-wide band is obvious.

The first MVP figure is the NH3 × C2H2/HCN inventory colormap; the ΔG-sweep line and
NH3×ΔG map populate from the ΔG substudy results.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")  # headless/batch safe
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


def _save(fig, output_dir: str | Path, stem: str, plots_cfg: Dict[str, Any]) -> List[Path]:
    """Save a figure to every requested format; return the written paths."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    formats = plots_cfg.get("formats", ["png"]) or ["png"]
    dpi = int(plots_cfg.get("dpi", 300))
    written = []
    for fmt in formats:
        path = output_dir / f"{stem}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        written.append(path)
    plt.close(fig)
    return written


def _edges(centers: np.ndarray) -> np.ndarray:
    """Cell edges from centers (handles a single unique value gracefully)."""
    centers = np.asarray(centers, dtype=float)
    if centers.size == 1:
        c = centers[0]
        step = abs(c) * 0.05 or 0.5
        return np.array([c - step, c + step])
    mids = (centers[:-1] + centers[1:]) / 2.0
    first = centers[0] - (mids[0] - centers[0])
    last = centers[-1] + (centers[-1] - mids[-1])
    return np.concatenate([[first], mids, [last]])


def _pivot(df: pd.DataFrame, x_col: str, y_col: str, value_col: str) -> pd.DataFrame:
    """Pivot a regular grid into a (y × x) table, sorted on both axes.

    ``dropna=False`` keeps fully-inaccessible rows/columns (all-NaN ``log10_X_eq``,
    e.g. the C2H2/HCN = 0 row) so they render as masked gray cells and still anchor
    the accessibility overlay, instead of silently vanishing from the grid.
    """
    pivot = df.pivot_table(index=y_col, columns=x_col, values=value_col,
                           aggfunc="first", dropna=False)
    return pivot.sort_index().sort_index(axis=1)


def _accessible_grid(df: pd.DataFrame, x_col: str, y_col: str,
                     ref_pivot: pd.DataFrame) -> np.ndarray:
    """0/1 accessibility field aligned to ``ref_pivot``.

    Built from ``formed_bool``; any cell with no/failed data counts as 0 (not
    accessible), so the boundary reflects formation, not a log10X level.
    """
    acc = df.copy()
    acc["_acc"] = acc["formed_bool"].astype(str).str.strip().str.lower().isin(
        ["true", "1", "1.0"]).astype(float)
    piv = acc.pivot_table(index=y_col, columns=x_col, values="_acc", aggfunc="max",
                          dropna=False)
    piv = piv.reindex(index=ref_pivot.index, columns=ref_pivot.columns)
    return piv.fillna(0.0).to_numpy(dtype=float)


def _auto_limits(df: pd.DataFrame):
    """Robust (2nd, 98th) percentile color limits from finite ``log10_X_eq``.

    Falls back to min/max when too few finite points; returns ``(None, None)`` when
    there is no finite data (matplotlib then auto-scales).
    """
    vals = pd.to_numeric(df.get("log10_X_eq"), errors="coerce").to_numpy(dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size >= 5:
        vmin, vmax = float(np.percentile(vals, 2)), float(np.percentile(vals, 98))
    elif vals.size > 0:
        vmin, vmax = float(vals.min()), float(vals.max())
    else:
        return None, None
    if vmin == vmax:
        vmin, vmax = vmin - 0.5, vmax + 0.5
    return vmin, vmax


def _render_heatmap(df, x_col, y_col, cfg, vmin, vmax, title, output_dir, stem):
    """Render one heatmap at the given color limits, with the accessibility overlay."""
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    cmap = plt.get_cmap(cfg.get("cmap", "viridis")).with_extremes(
        bad=cfg.get("show_failed_as", "gray"))

    logp = _pivot(df, x_col, y_col, "log10_X_eq")
    x = logp.columns.to_numpy(dtype=float)
    y = logp.index.to_numpy(dtype=float)
    Z = np.ma.masked_invalid(logp.to_numpy(dtype=float))
    mesh = ax.pcolormesh(_edges(x), _edges(y), Z, cmap=cmap, vmin=vmin, vmax=vmax,
                         shading="flat")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(cfg.get("colorbar_label", "log10 equilibrium mole fraction"))

    # --- Primary delineation: accessibility boundary from formed_bool ---
    A = _accessible_grid(df, x_col, y_col, logp)
    boundary_color = cfg.get("boundary_color", "red")
    has_grid = x.size >= 2 and y.size >= 2
    inaccessible = bool((A < 0.5).any())
    accessible = bool((A >= 0.5).any())
    handles = []
    if has_grid and inaccessible and accessible:
        try:
            ax.contour(x, y, A, levels=[0.5], colors=[boundary_color], linewidths=2)
            handles.append(Line2D([0], [0], color=boundary_color, lw=2,
                                  label="accessibility boundary"))
        except Exception:
            pass
    if inaccessible:
        if has_grid:
            try:  # hatch the inaccessible region so a 1-cell band is unmistakable
                ax.contourf(x, y, A, levels=[-0.5, 0.5], colors="none", hatches=["xxx"])
            except Exception:
                pass
        handles.append(Patch(facecolor=cfg.get("show_failed_as", "gray"),
                             hatch="xxx", edgecolor=boundary_color,
                             label="inaccessible (below threshold)"))

    # --- Secondary, optional: log10X threshold contour (only if actually crossed) ---
    thr = cfg.get("threshold_contour_log10X")
    if thr is not None and has_grid and Z.count() > 1:
        zf = Z.filled(np.nan)
        zmin, zmax = np.nanmin(zf), np.nanmax(zf)
        if np.isfinite(zmin) and np.isfinite(zmax) and zmin < float(thr) < zmax:
            try:
                ax.contour(x, y, zf, levels=[float(thr)], colors="white",
                           linewidths=1.2, linestyles=":")
            except Exception:
                pass

    if handles:
        ax.legend(handles=handles, loc="upper left", fontsize=8, framealpha=0.85)
    ax.set_title(title)
    ax.set_xlabel(cfg.get("x_label", x_col))
    ax.set_ylabel(cfg.get("y_label", y_col))
    return _save(fig, output_dir, stem, cfg)


def _heatmap(df, x_col, y_col, cfg, output_dir, stem, default_title):
    """Render the fixed and/or autoscaled versions selected by ``cfg['scales']``."""
    if df.empty:
        return []
    scales = [str(s).lower() for s in (cfg.get("scales") or ["fixed", "auto"])]
    base_title = cfg.get("title", default_title)
    written: List[Path] = []
    if "fixed" in scales:
        written += _render_heatmap(df, x_col, y_col, cfg, cfg.get("vmin"),
                                   cfg.get("vmax"), base_title, output_dir, stem)
    if "auto" in scales:
        vmin, vmax = _auto_limits(df)
        auto_title = cfg.get("title_autoscaled") or f"{base_title} (autoscaled)"
        written += _render_heatmap(df, x_col, y_col, cfg, vmin, vmax, auto_title,
                                   output_dir, f"{stem}_autoscaled")
    return written


def plot_inventory_landscape(grid_df: pd.DataFrame, plots_cfg: Dict[str, Any],
                             output_dir: str | Path) -> List[Path]:
    """NH3 × C2H2/HCN colormap of log10 X_eq. Failed cells shown gray."""
    df = grid_df[grid_df.get("substudy_id") == "inventory_landscape"] \
        if "substudy_id" in grid_df else grid_df
    if df.empty:
        return []
    cfg = (plots_cfg or {}).get("inventory_landscape", {}) or {}
    cfg = {**plots_cfg, **cfg}  # inherit formats/dpi from the top-level plots block
    return _heatmap(df, "NH3_mol", "C2H2_over_HCN", cfg, output_dir,
                    "inventory_landscape", "Inventory accessibility landscape")


def plot_nh3_deltaG_landscape(grid_df: pd.DataFrame, plots_cfg: Dict[str, Any],
                              output_dir: str | Path) -> List[Path]:
    """NH3 × ΔG-offset colormap of log10 X_eq (from the nh3_deltaG_landscape cases)."""
    df = grid_df[grid_df.get("substudy_id") == "nh3_deltaG_landscape"] \
        if "substudy_id" in grid_df else grid_df
    if df.empty:
        return []
    cfg = (plots_cfg or {}).get("nh3_deltaG_landscape", {}) or {}
    cfg = {**plots_cfg, **cfg}
    return _heatmap(df, "NH3_mol", "deltaG_offset_kJ_mol", cfg, output_dir,
                    "nh3_deltaG_landscape", "NH3 threshold vs Gibbs offset")


def plot_deltaG_sweep(grid_df: pd.DataFrame, plots_cfg: Dict[str, Any],
                      output_dir: str | Path) -> List[Path]:
    """ΔG-offset vs log10 X_eq line (from the deltaG_sweep cases)."""
    df = grid_df[grid_df.get("substudy_id") == "deltaG_sweep"] \
        if "substudy_id" in grid_df else grid_df
    if df.empty:
        return []
    cfg = (plots_cfg or {}).get("deltaG_sweep", {}) or {}
    cfg = {**plots_cfg, **cfg}
    d = df.sort_values("deltaG_offset_kJ_mol")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(d["deltaG_offset_kJ_mol"], d["log10_X_eq"], marker="o")
    thr = cfg.get("threshold_line_log10X", -12)
    if thr is not None:
        ax.axhline(thr, linestyle="--", color="gray", linewidth=1)
    ax.axvline(0, linestyle=":", color="black", linewidth=1)
    window = cfg.get("shade_window_kJ")
    if window:
        ax.axvspan(-float(window), float(window), color="gray", alpha=0.12)
    ax.set_title(cfg.get("title", "Accessibility vs Gibbs offset"))
    ax.set_xlabel(cfg.get("x_label", "Delta G offset (kJ/mol)"))
    ax.set_ylabel(cfg.get("y_label", "log10 equilibrium mole fraction"))
    return _save(fig, output_dir, "deltaG_sweep", cfg)


def plot_all(grid_df: pd.DataFrame, plots_cfg: Dict[str, Any], output_dir: str | Path,
             substudy: Optional[str] = None) -> Dict[str, List[Path]]:
    """Render every plot whose substudy has data; return name -> written paths."""
    plotters = {
        "inventory_landscape": plot_inventory_landscape,
        "deltaG_sweep": plot_deltaG_sweep,
        "nh3_deltaG_landscape": plot_nh3_deltaG_landscape,
    }
    out: Dict[str, List[Path]] = {}
    for name, fn in plotters.items():
        if substudy and name != substudy:
            continue
        out[name] = fn(grid_df, plots_cfg, output_dir)
    return out
