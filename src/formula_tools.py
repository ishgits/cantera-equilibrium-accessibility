"""Formula parsing and elemental bookkeeping utilities."""
from __future__ import annotations

import re
from typing import Dict

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)([0-9]*\.?[0-9]*)")


def parse_formula(formula: str) -> Dict[str, float]:
    """Parse a simple chemical formula into a composition dictionary.

    Supports formulas like H2O, C5H5N5, C4H5N3O. Parentheses and hydrates are
    intentionally not supported; use already-expanded formulas in species CSV.
    """
    formula = str(formula).strip()
    if not formula:
        raise ValueError("Empty formula.")
    composition: Dict[str, float] = {}
    matched = ""
    for element, count_text in _FORMULA_TOKEN.findall(formula):
        matched += element + count_text
        count = float(count_text) if count_text else 1.0
        composition[element] = composition.get(element, 0.0) + count
    if matched != formula:
        raise ValueError(
            f"Formula '{formula}' could not be fully parsed. Use an expanded formula like C5H5N5."
        )
    # Cantera YAML prefers ints where possible.
    cleaned = {}
    for k, v in composition.items():
        cleaned[k] = int(v) if float(v).is_integer() else float(v)
    return cleaned


def composition_to_formula(composition: Dict[str, float]) -> str:
    """Create a simple formula string from a composition dictionary."""
    parts = []
    for element in sorted(composition):
        count = composition[element]
        if count == 1:
            parts.append(element)
        else:
            parts.append(f"{element}{int(count) if float(count).is_integer() else count}")
    return "".join(parts)
