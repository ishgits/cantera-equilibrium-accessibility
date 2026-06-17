"""Tests for the target-accessibility formation-call logic."""
import numpy as np
import pandas as pd

from result_summary import summarize_target_formation


def _run(target, species_rows, yaml_file):
    """Build long-format rows for one equilibrium run."""
    rows = []
    for sp, x_eq, n_eq, status in species_rows:
        rows.append({
            "scenario": "s1",
            "model_mode": "single_product",
            "yaml_file": yaml_file,
            "target_product": target,
            "T_C": 20.0,
            "species": sp,
            "X_eq": x_eq,
            "n_eq_mol": n_eq,
            "initial_moles": 0.0,
            "solver_status": status,
            "element_balance_relative_spread": 0.0,
        })
    return rows


def _frame():
    rows = []
    # significant: X_eq >= 1e-6
    rows += _run("P1(aq)", [("P1(aq)", 1e-3, 1e-3, "ok")], "p1.yaml")
    # trace: 1e-12 <= X_eq < 1e-6
    rows += _run("P2(aq)", [("P2(aq)", 1e-9, 1e-9, "ok")], "p2.yaml")
    # below_threshold: X_eq < 1e-12
    rows += _run("P3(aq)", [("P3(aq)", 0.0, 0.0, "ok")], "p3.yaml")
    # solver_failed
    rows += _run("P4(aq)", [("P4(aq)", np.nan, np.nan, "failed")], "p4.yaml")
    # target not present in the YAML (no row for the target species)
    rows += _run("P5(aq)", [("H2O(l)", 0.9, 1.0, "ok")], "p5.yaml")
    return pd.DataFrame(rows)


def test_formation_calls():
    out = summarize_target_formation(
        _frame(),
        x_threshold=1e-12,
        n_threshold_mol=0.0,
        x_significant_threshold=1e-6,
    )
    call = dict(zip(out["target_product"], out["formation_call"]))
    assert call["P1(aq)"] == "significant"
    assert call["P2(aq)"] == "trace"
    assert call["P3(aq)"] == "below_threshold"
    assert call["P4(aq)"] == "solver_failed"
    assert call["P5(aq)"] == "target_not_present_in_yaml"


def test_formed_bool_matches_calls():
    out = summarize_target_formation(_frame())
    formed = dict(zip(out["target_product"], out["formed_bool"]))
    assert formed["P1(aq)"] is True or bool(formed["P1(aq)"]) is True
    assert bool(formed["P2(aq)"]) is True   # trace still counts as formed
    assert bool(formed["P3(aq)"]) is False
    assert bool(formed["P4(aq)"]) is False
