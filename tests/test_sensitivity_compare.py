"""Tests for cross-target aggregation + the combined NH3 helper (Cantera-free)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from sensitivity_compare import (
    build_cross_target_table, build_nh3_combined, classify_discriminators, load_campaign,
)

SPECIES_CSV = PROJECT_ROOT / "inputs" / "amino_acids_species.csv"


def _aa_summary(target, accessible_frac, dg_log10x, nh3=0.05):
    """Synthetic case summary: inventory grid + a ΔG sweep."""
    rows = []
    n, n_acc = 10, int(round(accessible_frac * 10))
    for i in range(n):
        formed = i < n_acc
        ratio = (i + 1) * 0.5
        x = 1e-2 if formed else 0.0
        rows.append(dict(substudy_id="inventory_landscape", case_id=f"INV{i}",
                         target_product=target, formed_bool=formed, X_eq=x,
                         log10_X_eq=(-2.0 if formed else np.nan), n_eq_mol=x * 0.04,
                         NH3_mol=nh3, HCN_mol=0.02, C2H2_mol=ratio * 0.02,
                         C2H2_over_HCN=ratio, H2O_mol=1.0, deltaG_offset_kJ_mol=0.0,
                         solver_status="ok", runtime_seconds=0.001, suspect_balance=False,
                         formation_call=("significant" if formed else "below_threshold")))
    for off in range(-200, 201, 40):
        lx = dg_log10x(off)
        formed = lx >= -6
        rows.append(dict(substudy_id="deltaG_sweep", case_id=f"DG{off}",
                         target_product=target, formed_bool=formed, X_eq=10.0 ** lx,
                         log10_X_eq=lx, NH3_mol=nh3, HCN_mol=0.02, C2H2_mol=0.042,
                         C2H2_over_HCN=2.1, H2O_mol=1.0, deltaG_offset_kJ_mol=float(off),
                         solver_status="ok", runtime_seconds=0.001, suspect_balance=False,
                         formation_call=("significant" if formed else "below_threshold")))
    return pd.DataFrame(rows)


def _campaign():
    return {
        "glycine": _aa_summary("Glycine(aq)", 1.0, lambda o: -2.0),               # robust
        "alanine": _aa_summary("Alanine(aq)", 0.8, lambda o: -2.0),               # carbon-limited
        "valine": _aa_summary("Valine(aq)", 0.9, lambda o: -2.0 - max(0, o) * 0.025),  # fragile
    }


def test_one_row_per_study_with_composition_and_metrics():
    table = build_cross_target_table(_campaign(), SPECIES_CSV)
    assert len(table) == 3
    gly = table[table["amino_acid"] == "glycine"].iloc[0]
    assert gly["formula"] == "C2H5NO2" and gly["n_C"] == 2 and gly["n_N"] == 1
    assert gly["inventory_accessible_fraction"] == pytest.approx(1.0)
    assert gly["X_eq_at_reference_inventory"] == pytest.approx(1e-2)   # ΔG at offset 0


def test_max_stoichiometric_yield_limiting_reagent():
    table = build_cross_target_table(_campaign(), SPECIES_CSV)
    gly = table[table["amino_acid"] == "glycine"].iloc[0]
    # C limits: (0.02*1 + 0.042*2) / 2 carbons = 0.104/2 = 0.052 mol.
    assert gly["max_stoichiometric_yield_mol"] == pytest.approx(0.052, abs=1e-6)


def test_classify_discriminators_label_scheme():
    camp = _campaign()
    # A baseline-INACCESSIBLE amino acid: X_eq at ΔG=0 below threshold, nothing formed.
    camp["serine"] = _aa_summary("Serine(aq)", 0.0, lambda o: -10.0)
    table = classify_discriminators(build_cross_target_table(camp, SPECIES_CSV))
    tag = dict(zip(table["amino_acid"], table["discriminator"]))
    assert tag["valine"] == "energetically_fragile"     # accessible at 0, finite crossing
    assert tag["alanine"] == "inventory_gated"          # accessible, below campaign max
    assert tag["glycine"] == "robust_accessible"        # max accessible, no crossing
    assert tag["serine"] == "not_accessible_in_batch"   # inaccessible at baseline


def test_combined_nh3_groups_tags_unlocked():
    from sensitivity_compare import combined_nh3_groups
    combined = pd.DataFrame({
        "amino_acid": ["gly", "gly", "ala", "pro"],
        "source_batch": ["A_no_nh3", "B_nh3", "A_no_nh3", "A_no_nh3"],
        "log10_X_eq": [-10.0, -2.0, -2.0, -10.0]})
    groups = combined_nh3_groups(combined)
    assert groups["gly"] == "nh3_unlocked"          # inaccessible no-NH3, accessible w/ NH3
    assert groups["ala"] == "accessible_no_nh3"
    assert groups["pro"] == "not_accessible"


def test_combine_aligns_zero_point_and_flags_unlocked():
    # Batch A (no NH3): inaccessible at ratio 2.1; Batch B (NH3): accessible.
    a = {"serine": _aa_summary("Serine(aq)", 0.0, lambda o: -2.0, nh3=0.0)}
    # Force the Batch-A inventory point at ratio 2.1 to be below threshold.
    a["serine"].loc[:, "log10_X_eq"] = a["serine"]["log10_X_eq"].where(
        a["serine"]["substudy_id"] != "inventory_landscape", -10.0)
    b = {"serine": _aa_summary("Serine(aq)", 1.0, lambda o: -2.0)}
    # add a ratio == 2.1 inventory column to both so the snap is exact
    for camp, src_lx in [(a, -10.0), (b, -2.0)]:
        cs = camp["serine"]
        for nh3 in ([0.0] if camp is a else [0.01, 0.05, 0.1]):
            cs.loc[len(cs)] = dict(substudy_id="inventory_landscape", case_id=f"R{nh3}",
                target_product="Serine(aq)", formed_bool=(src_lx >= -6), X_eq=10.0 ** src_lx,
                log10_X_eq=src_lx, n_eq_mol=10.0 ** src_lx * 0.04, NH3_mol=nh3, HCN_mol=0.02,
                C2H2_mol=0.042, C2H2_over_HCN=2.1, H2O_mol=1.0, deltaG_offset_kJ_mol=0.0,
                solver_status="ok", runtime_seconds=0.001, suspect_balance=False, formation_call="x")

    combined = build_nh3_combined(a, b, ratio=2.1)
    ser = combined[combined["amino_acid"] == "serine"]
    a_rows = ser[ser["source_batch"] == "A_no_nh3"]
    b_rows = ser[ser["source_batch"] == "B_nh3"]
    assert (a_rows["NH3_frac"] == 0.0).all()
    assert (b_rows["NH3_frac"] > 0).all()
    # unlocked-by-NH3: not accessible at 0%, accessible with NH3.
    assert (a_rows["log10_X_eq"] < -6).any()
    assert (b_rows["log10_X_eq"] >= -6).any()


def _write_campaign(scan_dir, campaign):
    for key, cs in campaign.items():
        (scan_dir / key / "results").mkdir(parents=True, exist_ok=True)
        cs.to_csv(scan_dir / key / "results" / "sensitivity_case_summary.csv", index=False)


def test_aggregate_script_writes_outputs(tmp_path):
    import compare_amino_acids as agg
    scan = tmp_path / "aa_nh3"
    _write_campaign(scan, _campaign())
    out = tmp_path / "agg"
    table = agg.aggregate(scan, SPECIES_CSV, out)
    for f in ["amino_acid_metrics.csv", "amino_acid_case_summary.csv", "SCHEMA.md",
              "comparison_summary.md"]:
        assert (out / f).exists()
    assert (out / "figures" / "ranked_accessible_fraction.png").exists()
    assert len(table) == 3
    # concatenated case summary carries an amino_acid column
    concat = pd.read_csv(out / "amino_acid_case_summary.csv")
    assert set(concat["amino_acid"]) == {"glycine", "alanine", "valine"}


def test_combined_script_writes_and_idempotent(tmp_path):
    import plot_nh3_combined as comb
    a = {"serine": _aa_summary("Serine(aq)", 0.0, lambda o: -2.0, nh3=0.0)}
    b = {"serine": _aa_summary("Serine(aq)", 1.0, lambda o: -2.0)}
    for camp, lx in [(a, -10.0), (b, -2.0)]:
        cs = camp["serine"]
        for nh3 in ([0.0] if camp is a else [0.01, 0.05]):
            cs.loc[len(cs)] = dict(substudy_id="inventory_landscape", case_id=f"R{nh3}",
                target_product="Serine(aq)", formed_bool=(lx >= -6), X_eq=10.0 ** lx,
                log10_X_eq=lx, n_eq_mol=0.006, NH3_mol=nh3, HCN_mol=0.02, C2H2_mol=0.042,
                C2H2_over_HCN=2.1, H2O_mol=1.0, deltaG_offset_kJ_mol=0.0, solver_status="ok",
                runtime_seconds=0.001, suspect_balance=False, formation_call="x")
    _write_campaign(tmp_path / "aa_no_nh3", a)
    _write_campaign(tmp_path / "aa_nh3", b)
    out = tmp_path / "combined"

    argv = ["--no-nh3", str(tmp_path / "aa_no_nh3"), "--nh3", str(tmp_path / "aa_nh3"),
            "--out", str(out), "--ratio", "2.1"]
    assert comb.main(argv) == 0
    assert (out / "nh3_combined.csv").exists()
    assert (out / "nh3_combined_heatmap_log10X.png").exists()
    assert (out / "nh3_combined_heatmap_yield_pct_HCN.png").exists()
    assert (out / "summary.md").exists()
    # yield-relative-to-HCN columns are present (paper units).
    nc = pd.read_csv(out / "nh3_combined.csv")
    assert {"n_eq_mol", "HCN_mol", "ratio_snap", "yield_fraction_HCN", "yield_pct_HCN"} <= set(nc.columns)
    assert comb.main(argv) == 0          # idempotent re-run


def test_combined_yield_columns(tmp_path):
    # Batch A inaccessible, Batch B accessible at ratio 2.1 with known n_eq/HCN.
    a = {"glycine": _aa_summary("Glycine(aq)", 0.0, lambda o: -2.0, nh3=0.0)}
    b = {"glycine": _aa_summary("Glycine(aq)", 1.0, lambda o: -2.0)}
    for camp, lx in [(a, -10.0), (b, -2.0)]:
        cs = camp["glycine"]
        for nh3 in ([0.0] if camp is a else [0.01, 0.05]):
            cs.loc[len(cs)] = dict(substudy_id="inventory_landscape", case_id=f"R{nh3}",
                target_product="Glycine(aq)", formed_bool=(lx >= -6), X_eq=10.0 ** lx,
                log10_X_eq=lx, n_eq_mol=0.006, NH3_mol=nh3, HCN_mol=0.02, C2H2_mol=0.042,
                C2H2_over_HCN=2.1, H2O_mol=1.0, deltaG_offset_kJ_mol=0.0, solver_status="ok",
                runtime_seconds=0.001, suspect_balance=False, formation_call="x")
    combined = build_nh3_combined(a, b, ratio=2.1)
    b_rows = combined[(combined["amino_acid"] == "glycine") & (combined["source_batch"] == "B_nh3")]
    # yield_pct_HCN = 100 * n_eq_mol / HCN_mol = 100 * 0.006 / 0.02 = 30.
    assert b_rows["yield_pct_HCN"].tolist() == pytest.approx([30.0] * len(b_rows))
    assert b_rows["yield_fraction_HCN"].tolist() == pytest.approx([0.3] * len(b_rows))


def test_assemble_bridge_table(tmp_path):
    from sensitivity_compare import assemble_bridge_table
    a = {"glycine": _aa_summary("Glycine(aq)", 0.0, lambda o: -2.0, nh3=0.0)}
    b = {"glycine": _aa_summary("Glycine(aq)", 1.0, lambda o: -2.0)}
    for camp, lx in [(a, -10.0), (b, -2.0)]:
        cs = camp["glycine"]
        for nh3 in ([0.0] if camp is a else [0.01, 0.05, 0.10]):
            cs.loc[len(cs)] = dict(substudy_id="inventory_landscape", case_id=f"R{nh3}",
                target_product="Glycine(aq)", formed_bool=(lx >= -6), X_eq=10.0 ** lx,
                log10_X_eq=lx, n_eq_mol=0.006, NH3_mol=nh3, HCN_mol=0.02, C2H2_mol=0.042,
                C2H2_over_HCN=2.1, H2O_mol=1.0, deltaG_offset_kJ_mol=0.0, solver_status="ok",
                runtime_seconds=0.001, suspect_balance=False, formation_call="x")
    table = assemble_bridge_table(a, b, ratio=2.1)
    r = table[table["amino_acid"] == "glycine"].iloc[0]
    assert not bool(r["accessible_no_nh3"]) and bool(r["accessible_with_nh3"])
    assert r["paper_group"] == "nh3_unlocked"
    assert r["min_NH3_significant"] == pytest.approx(0.01)
    assert r["yield_pct_HCN_1pct"] == pytest.approx(30.0)
    assert "unlocked by NH3" in r["workflow_interpretation"]


def test_load_campaign_reads_from_disk(tmp_path):
    scan = tmp_path / "aa_nh3"
    cs = _aa_summary("Glycine(aq)", 1.0, lambda o: -2.0)
    (scan / "glycine" / "results").mkdir(parents=True)
    cs.to_csv(scan / "glycine" / "results" / "sensitivity_case_summary.csv", index=False)
    campaign = load_campaign(scan)
    assert set(campaign) == {"glycine"}
    assert len(campaign["glycine"]) == len(cs)
