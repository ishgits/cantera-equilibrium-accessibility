"""Tests for the Phase 3 sensitivity summary (classifier reuse, target_variant
matching, safe log10, failed/suspect handling). No Cantera required."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from result_summary import classify_formation
from sensitivity_summary import (
    _accessible_window, compute_sensitivity_metrics, make_landscape_grid,
    render_mvp_verdict, summarize_sensitivity_cases, write_sensitivity_run_summary,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Species table including a ΔG variant of alanine (same formula as the base).
SPECIES_DF = pd.DataFrame({
    "cantera_name": ["H2O(l)", "HCN(aq)", "C2H2(aq)", "NH3(aq)", "Alanine(aq)",
                     "Alanine__dG_p020(aq)"],
    "formula": ["H2O", "CHN", "C2H2", "H3N", "C3H7NO2", "C3H7NO2"],
})

_RAW_COLS_DEFAULTS = {
    "study_id": "alanine_mvp", "substudy_id": "inventory_landscape",
    "target_product": "Alanine(aq)", "model_id": "M_x", "T_C": 0.0, "P_Pa": 101325.0,
    "runtime_seconds": 0.01, "HCN_mol": 0.02, "C2H2_mol": 0.0, "C2H2_over_HCN": 0.0,
    "deltaG_offset_kJ_mol": 0.0, "error_message": "",
}


def _case_rows(case_id, target_variant, species_x, status="ok", **design):
    """Build raw-long rows for one case (one row per species)."""
    meta = {**_RAW_COLS_DEFAULTS, **design, "case_id": case_id,
            "target_variant": target_variant, "solver_status": status}
    rows = []
    for sp, x in species_x.items():
        init = 1.0 if sp == "H2O(l)" else 0.0
        rows.append({**meta, "species": sp,
                     "X_eq": (x if status == "ok" else np.nan), "initial_moles": init})
    return rows


# --------------------------------------------------------------------------- #
# Classifier (shared with the base workflow)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("x_eq,status,expected", [
    (1e-3, "ok", (True, "significant")),
    (1e-9, "ok", (True, "trace")),
    (1e-15, "ok", (False, "below_threshold")),
    (np.nan, "failed", (False, "solver_failed")),
])
def test_classify_formation_calls(x_eq, status, expected):
    assert classify_formation(x_eq, 1.0, status, 1e-12, 1e-6, 0.0) == expected


def test_classify_formation_safe_with_nan_when_ok():
    formed, call = classify_formation(np.nan, np.nan, "ok", 1e-12, 1e-6, 0.0)
    assert (formed, call) == (False, "below_threshold")


# --------------------------------------------------------------------------- #
# Case summary
# --------------------------------------------------------------------------- #
def test_summary_keys_on_case_and_matches_target_variant():
    # target_product is the bare alanine, but the modelled species is the variant.
    raw = pd.DataFrame(
        _case_rows("c1", "Alanine__dG_p020(aq)",
                   {"H2O(l)": 0.9, "Alanine__dG_p020(aq)": 0.1},
                   substudy_id="deltaG_sweep", deltaG_offset_kJ_mol=20.0))
    summary = summarize_sensitivity_cases(raw, SPECIES_DF, {})
    assert len(summary) == 1
    row = summary.iloc[0]
    assert row["case_id"] == "c1"
    assert row["target_variant"] == "Alanine__dG_p020(aq)"
    assert row["X_eq"] == pytest.approx(0.1)        # matched the variant row, not NaN
    assert row["formation_call"] == "significant"


def test_failed_case_retained_as_solver_failed():
    raw = pd.DataFrame(
        _case_rows("c_fail", "Alanine(aq)", {"H2O(l)": np.nan, "Alanine(aq)": np.nan},
                   status="failed"))
    summary = summarize_sensitivity_cases(raw, SPECIES_DF, {})
    assert len(summary) == 1                          # not dropped
    assert summary.iloc[0]["formation_call"] == "solver_failed"
    assert pd.isna(summary.iloc[0]["X_eq"])
    assert pd.isna(summary.iloc[0]["log10_X_eq"])     # safe log10


def test_suspect_balance_gate():
    # X_eq deliberately inconsistent across elements -> large balance spread.
    raw = pd.DataFrame(
        _case_rows("c_bal", "Alanine(aq)", {"H2O(l)": 0.9, "Alanine(aq)": 0.1}))
    summary = summarize_sensitivity_cases(raw, SPECIES_DF, {"balance_tol": 1e-6})
    assert bool(summary.iloc[0]["suspect_balance"]) is True
    # A very loose tolerance clears the flag.
    loose = summarize_sensitivity_cases(raw, SPECIES_DF, {"balance_tol": 10.0})
    assert bool(loose.iloc[0]["suspect_balance"]) is False


def test_log10_safe_for_zero_x_eq():
    raw = pd.DataFrame(
        _case_rows("c0", "Alanine(aq)", {"H2O(l)": 1.0, "Alanine(aq)": 0.0}))
    summary = summarize_sensitivity_cases(raw, SPECIES_DF, {})
    assert pd.isna(summary.iloc[0]["log10_X_eq"])
    assert summary.iloc[0]["formation_call"] == "below_threshold"


# --------------------------------------------------------------------------- #
# Grid + metrics + md
# --------------------------------------------------------------------------- #
def _inventory_summary():
    raw_rows = []
    for i, (nh3, ratio, x) in enumerate([
        (0.0, 0.0, 1e-15), (0.0, 5.0, 1e-9), (0.15, 0.0, 1e-3), (0.15, 5.0, 1e-2),
    ]):
        raw_rows += _case_rows(
            f"INV_{i}", "Alanine(aq)",
            {"H2O(l)": 1 - x, "Alanine(aq)": x},
            NH3_mol=nh3, C2H2_over_HCN=ratio, C2H2_mol=ratio * 0.02)
    # add one failed case
    raw_rows += _case_rows("INV_fail", "Alanine(aq)",
                           {"H2O(l)": np.nan, "Alanine(aq)": np.nan},
                           status="failed", NH3_mol=0.075, C2H2_over_HCN=2.5)
    return summarize_sensitivity_cases(pd.DataFrame(raw_rows), SPECIES_DF, {})


def test_landscape_grid_keeps_failed_and_has_columns():
    grid = make_landscape_grid(_inventory_summary())
    assert "INV_fail" in set(grid["case_id"])
    # Self-contained for ΔG comparison / ML: target_variant + model_id included.
    assert {"NH3_mol", "C2H2_over_HCN", "log10_X_eq", "formation_call",
            "target_variant", "model_id"} <= set(grid.columns)


# --------------------------------------------------------------------------- #
# Phase 4 — ΔG robustness + NH3×ΔG threshold metrics
# --------------------------------------------------------------------------- #
def _dg_summary(offsets, log10x):
    return pd.DataFrame({
        "case_id": [f"D{i}" for i in range(len(offsets))],
        "substudy_id": "deltaG_sweep",
        "deltaG_offset_kJ_mol": [float(o) for o in offsets],
        "log10_X_eq": [float(v) for v in log10x],
        "X_eq": [10.0 ** v for v in log10x],
        "formed_bool": [v >= -6 for v in log10x],
        "solver_status": "ok", "runtime_seconds": 0.01, "suspect_balance": False,
    })


def test_accessible_window_interpolates():
    import math
    neg, pos = _accessible_window([-40, -20, 0, 20, 40], [-1, -2, -3, -7, -9], -6.0)
    assert pos == pytest.approx(15.0)          # crosses -6 between offsets 0 and 20
    assert math.isnan(neg)                     # never drops below on the negative side


def test_deltaG_robustness_crossing_and_booleans():
    m = compute_sensitivity_metrics(_dg_summary([-40, -20, 0, 20, 40],
                                                [-1, -2, -3, -7, -9]))["deltaG_sweep"]
    assert m["offset_crossing_positive_kJ"] == pytest.approx(15.0)
    assert m["offset_crossing_negative_kJ"] is None
    assert m["robust_to_plus_minus_20_kJ"] is False     # crossing at +15 < 20
    assert m["robust_to_plus_minus_40_kJ"] is False
    assert m["X_eq_at_zero_offset"] == pytest.approx(1e-3)


def test_deltaG_robustness_requires_baseline_accessible():
    # Baseline-INACCESSIBLE (X_eq at ΔG=0 below threshold): robustness must be False.
    inacc = compute_sensitivity_metrics(
        _dg_summary([-40, 0, 40], [-10, -10, -10]))["deltaG_sweep"]
    assert inacc["accessible_at_zero_offset"] is False
    assert inacc["robust_to_plus_minus_20_kJ"] is False
    assert inacc["robust_to_plus_minus_40_kJ"] is False
    # Baseline-accessible across the whole sweep: robust True.
    acc = compute_sensitivity_metrics(
        _dg_summary([-200, 0, 200], [-2, -2, -2]))["deltaG_sweep"]
    assert acc["accessible_at_zero_offset"] is True
    assert acc["robust_to_plus_minus_20_kJ"] is True
    assert acc["robust_to_plus_minus_40_kJ"] is True


def test_deltaG_fully_robust_when_no_crossing():
    m = compute_sensitivity_metrics(_dg_summary([-200, -100, 0, 100, 200],
                                                [-2, -2, -2, -2, -2]))["deltaG_sweep"]
    assert m["offset_crossing_positive_kJ"] is None
    assert m["offset_crossing_negative_kJ"] is None
    assert m["robust_to_plus_minus_20_kJ"] is True
    assert m["robust_to_plus_minus_40_kJ"] is True


def test_nh3_deltaG_threshold_shift():
    rows = []
    spec = {-40: [True, True, True], 0: [False, True, True], 40: [False, False, False]}
    for off, flags in spec.items():
        for nh3, f in zip([0.0, 0.05, 0.1], flags):
            rows.append({
                "case_id": f"N{off}_{nh3}", "substudy_id": "nh3_deltaG_landscape",
                "deltaG_offset_kJ_mol": float(off), "NH3_mol": nh3, "formed_bool": f,
                "X_eq": 1e-3 if f else 1e-15, "log10_X_eq": -3.0 if f else -15.0,
                "solver_status": "ok", "runtime_seconds": 0.01, "suspect_balance": False,
            })
    m = compute_sensitivity_metrics(pd.DataFrame(rows))["nh3_deltaG_landscape"]
    assert m["nh3_threshold_at_min_offset"] == 0.0
    assert m["nh3_threshold_at_zero_offset"] == 0.05
    assert m["nh3_threshold_at_max_offset"] is None      # no accessible NH3 at +40
    assert m["threshold_shifts_with_offset"] is True


def test_render_mvp_verdict_scale_vs_revise():
    cfg = {"mode": {"target_products": ["Alanine(aq)"]}}
    scale = {"general": {"total_cases": 100, "successful_cases": 100, "failed_cases": 0,
                         "failure_rate": 0.0, "suspect_balance_cases": 0,
                         "median_runtime_sec": 0.001, "projected_runtime_sec_targets": 20.0,
                         "n_targets_projection": 18},
             "inventory_landscape": {"accessible_area_fraction": 0.9,
                                     "min_C2H2_over_HCN_accessible": 0.2,
                                     "min_NH3_accessible": 0.0}}
    md = render_mvp_verdict(scale, cfg)
    assert "SCALE" in md and "REVISE" not in md.split("Recommendation")[1]

    revise = {"general": {"total_cases": 100, "successful_cases": 95, "failed_cases": 5,
                          "failure_rate": 0.05, "suspect_balance_cases": 0,
                          "median_runtime_sec": 0.001, "projected_runtime_sec_targets": 20.0,
                          "n_targets_projection": 18},
              "inventory_landscape": {"accessible_area_fraction": 0.2,
                                      "min_C2H2_over_HCN_accessible": 1.0,
                                      "min_NH3_accessible": 0.05}}
    assert "REVISE" in render_mvp_verdict(revise, cfg)


def test_metrics_and_summary_md(tmp_path):
    summary = _inventory_summary()
    metrics = compute_sensitivity_metrics(summary)
    g = metrics["general"]
    assert g["total_cases"] == 5
    assert g["failed_cases"] == 1
    assert g["failure_rate"] == pytest.approx(0.2)
    assert metrics["inventory_landscape"]["n_cases"] == 5
    md = write_sensitivity_run_summary(metrics, tmp_path / "s.md", study_id="alanine_mvp")
    text = md.read_text()
    assert "Total cases: 5" in text
    assert "Inventory landscape" in text
