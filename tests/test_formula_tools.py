"""Tests for chemical-formula parsing."""
import pytest

from formula_tools import parse_formula, composition_to_formula


def test_parse_simple():
    assert parse_formula("H2O") == {"H": 2, "O": 1}


def test_parse_multi_element():
    assert parse_formula("C5H5N5") == {"C": 5, "H": 5, "N": 5}
    assert parse_formula("C4H4N2O2") == {"C": 4, "H": 4, "N": 2, "O": 2}


def test_single_atom_counts_as_one():
    assert parse_formula("HCN") == {"H": 1, "C": 1, "N": 1}


def test_fractional_counts_supported():
    comp = parse_formula("C1.5H3")
    assert comp["C"] == pytest.approx(1.5)
    assert comp["H"] == 3


def test_empty_formula_raises():
    with pytest.raises(ValueError):
        parse_formula("")


def test_unparseable_formula_raises():
    # Parentheses are intentionally unsupported and must be expanded first.
    with pytest.raises(ValueError):
        parse_formula("Ca(OH)2")


def test_composition_to_formula_roundtrip():
    formula = "C4H4N2O2"
    assert parse_formula(composition_to_formula(parse_formula(formula))) == parse_formula(formula)
