"""Tests for Phase 3 config-driven plotting (no Cantera; headless Agg).

Covers the plot fix-ups: dual fixed/autoscaled heatmaps, autoscale limits derived
from finite data, the formed_bool accessibility boundary, and inaccessible cells
staying visible.
"""
import numpy as np
import pandas as pd

from sensitivity_plotting import (
    _accessible_grid, _auto_limits, _pivot, plot_all, plot_inventory_landscape,
)


def _grid(points=4, inaccessible_ratio_rows=(), failed_cells=()):
    """Regular NH3 × C2H2/HCN grid. Rows in ``inaccessible_ratio_rows`` (by index)
    are below threshold; ``failed_cells`` (by flat index) are solver failures."""
    nh3 = np.linspace(0.0, 0.15, points)
    ratio = np.linspace(0.2, 5.0, points)
    rows, k = [], 0
    for ri, r in enumerate(ratio):
        for n in nh3:
            failed = k in failed_cells
            inaccessible = ri in inaccessible_ratio_rows
            if failed:
                log10x, formed, call, status = np.nan, False, "solver_failed", "failed"
            elif inaccessible:
                log10x, formed, call, status = np.nan, False, "below_threshold", "ok"
            else:
                log10x = -2.4 + 0.7 * (r / 5.0)        # ~ -2.4 .. -1.7
                formed, call, status = True, "significant", "ok"
            rows.append({
                "case_id": f"c{k}", "study_id": "s", "substudy_id": "inventory_landscape",
                "target_product": "Alanine(aq)", "NH3_mol": n, "C2H2_over_HCN": r,
                "C2H2_mol": r * 0.02, "HCN_mol": 0.02, "deltaG_offset_kJ_mol": 0.0,
                "T_C": 0.0, "P_Pa": 101325.0,
                "X_eq": np.nan if np.isnan(log10x) else 10 ** log10x,
                "log10_X_eq": log10x, "n_eq_mol": np.nan if np.isnan(log10x) else 0.01,
                "formed_bool": formed, "formation_call": call,
                "solver_status": status, "runtime_seconds": 0.01,
            })
            k += 1
    return pd.DataFrame(rows)


PLOTS_CFG = {
    "formats": ["png"], "dpi": 80,
    "inventory_landscape": {
        "title": "Test", "x_label": "NH3", "y_label": "ratio",
        "cmap": "viridis", "vmin": -12, "vmax": 0, "show_failed_as": "gray",
    },
}


# --------------------------------------------------------------------------- #
# Dual fixed / autoscaled output (Fix 1)
# --------------------------------------------------------------------------- #
def test_default_emits_fixed_and_autoscaled(tmp_path):
    written = plot_inventory_landscape(_grid(inaccessible_ratio_rows=(0,)), PLOTS_CFG, tmp_path)
    names = {p.name for p in written}
    assert names == {"inventory_landscape.png", "inventory_landscape_autoscaled.png"}
    assert all(p.stat().st_size > 0 for p in written)


def test_scales_auto_only(tmp_path):
    cfg = {**PLOTS_CFG, "inventory_landscape": {**PLOTS_CFG["inventory_landscape"],
                                                "scales": ["auto"]}}
    written = plot_inventory_landscape(_grid(), cfg, tmp_path)
    assert [p.name for p in written] == ["inventory_landscape_autoscaled.png"]


def test_scales_fixed_only(tmp_path):
    cfg = {**PLOTS_CFG, "inventory_landscape": {**PLOTS_CFG["inventory_landscape"],
                                                "scales": ["fixed"]}}
    written = plot_inventory_landscape(_grid(), cfg, tmp_path)
    assert [p.name for p in written] == ["inventory_landscape.png"]


def test_autoscale_limits_from_finite_data():
    grid = _grid(points=5, inaccessible_ratio_rows=(0,))
    vmin, vmax = _auto_limits(grid)
    # Real data lives in ~[-2.4, -1.7]; limits must come from it, not -12..0.
    assert -12 < vmin <= vmax < 0
    finite = grid["log10_X_eq"].dropna()
    assert vmin >= finite.min() - 1e-9
    assert vmax <= finite.max() + 1e-9


def test_autoscale_limits_none_when_no_finite_data():
    grid = _grid(points=3, inaccessible_ratio_rows=(0, 1, 2))  # all NaN
    assert _auto_limits(grid) == (None, None)


# --------------------------------------------------------------------------- #
# Accessibility boundary + inaccessible visibility (Fix 2)
# --------------------------------------------------------------------------- #
def test_boundary_drawn_for_mixed_grid(tmp_path):
    # One inaccessible ratio row + the rest accessible -> a real boundary exists.
    written = plot_inventory_landscape(_grid(points=5, inaccessible_ratio_rows=(0,)),
                                       PLOTS_CFG, tmp_path)
    assert all(p.exists() and p.stat().st_size > 0 for p in written)


def test_fully_inaccessible_row_retained_in_grid():
    # Regression: an all-NaN log10X row (e.g. C2H2/HCN = 0) must not be dropped by
    # the pivot, or the accessibility field loses the inaccessible band entirely.
    grid = _grid(points=5, inaccessible_ratio_rows=(0,))
    logp = _pivot(grid, "NH3_mol", "C2H2_over_HCN", "log10_X_eq")
    assert logp.shape[0] == 5                       # the all-NaN row is kept
    A = _accessible_grid(grid, "NH3_mol", "C2H2_over_HCN", logp)
    assert A.shape == (5, 5)
    assert int((A < 0.5).sum()) == 5                # the whole bottom row inaccessible
    assert list(A[0]) == [0.0] * 5


def test_single_inaccessible_row_still_valid(tmp_path):
    # A 1-cell-wide inaccessible band must still render a valid figure.
    grid = _grid(points=6, inaccessible_ratio_rows=(0,))
    written = plot_inventory_landscape(grid, PLOTS_CFG, tmp_path)
    assert all(p.stat().st_size > 0 for p in written)


def test_all_accessible_grid_has_no_boundary_but_renders(tmp_path):
    written = plot_inventory_landscape(_grid(points=4), PLOTS_CFG, tmp_path)
    assert all(p.stat().st_size > 0 for p in written)


def test_failed_cells_render(tmp_path):
    written = plot_inventory_landscape(_grid(points=4, failed_cells=(5,)),
                                       PLOTS_CFG, tmp_path)
    assert all(p.stat().st_size > 0 for p in written)


# --------------------------------------------------------------------------- #
# plot_all + formats
# --------------------------------------------------------------------------- #
def test_plot_all_skips_substudies_without_data(tmp_path):
    out = plot_all(_grid(), PLOTS_CFG, tmp_path)
    assert out["inventory_landscape"]            # produced (2 files)
    assert out["deltaG_sweep"] == []             # no ΔG data yet
    assert out["nh3_deltaG_landscape"] == []


def test_multiple_formats_and_scales(tmp_path):
    cfg = {**PLOTS_CFG, "formats": ["png", "pdf"]}
    written = plot_inventory_landscape(_grid(inaccessible_ratio_rows=(0,)), cfg, tmp_path)
    # 2 scales × 2 formats = 4 files.
    assert len(written) == 4
    assert {p.suffix for p in written} == {".png", ".pdf"}
    assert {p.name for p in written} == {
        "inventory_landscape.png", "inventory_landscape.pdf",
        "inventory_landscape_autoscaled.png", "inventory_landscape_autoscaled.pdf"}
