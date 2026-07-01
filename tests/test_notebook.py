"""Execute the Phase 5 orchestrator notebook end-to-end on the alanine MVP.

Skipped cleanly when Cantera or nbclient/nbformat is unavailable, so CI stays green
without the optional notebook/equilibrium stack.
"""
import importlib.util
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = PROJECT_ROOT / "notebooks" / "02_sensitivity_landscape_workflow.ipynb"


def _missing(*mods):
    return [m for m in mods if importlib.util.find_spec(m) is None]


@pytest.mark.skipif(not NOTEBOOK.exists(), reason="orchestrator notebook not present")
def test_orchestrator_notebook_runs_end_to_end():
    missing = _missing("cantera", "nbclient", "nbformat", "ipykernel")
    if missing:
        pytest.skip(f"notebook execution needs: {', '.join(missing)}")

    import nbformat
    from nbclient import NotebookClient

    nb = nbformat.read(NOTEBOOK, as_version=4)
    client = NotebookClient(
        nb, timeout=600, kernel_name="python3",
        resources={"metadata": {"path": str(NOTEBOOK.parent)}})
    client.execute()  # raises if any cell errors

    results = PROJECT_ROOT / "studies" / "alanine_mvp" / "results"
    assert (results / "sensitivity_case_summary.csv").exists()
    assert (results / "sensitivity_landscape_grid.csv").exists()
    assert (results / "sensitivity_run_summary.md").exists()
    figures = PROJECT_ROOT / "studies" / "alanine_mvp" / "figures"
    assert (figures / "inventory_landscape.png").exists()
