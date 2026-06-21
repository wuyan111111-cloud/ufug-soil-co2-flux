"""Lightweight checks for a completed reproduction run."""

from pathlib import Path

import pandas as pd


def test_expected_outputs() -> None:
    results = Path("results")
    primary = pd.read_csv(results / "primary_csu_clustered_effects.csv")
    paired = pd.read_csv(results / "paired_effects_132.csv")
    integrity = pd.read_csv(results / "data_integrity_checks.csv")

    assert len(primary) == 4
    assert len(paired) == 132
    assert set(primary["Treatment"]) == {"DCD", "Biochar"}
    assert set(primary["Window"]) == {"July", "October"}
    assert (primary["percent_change"] < 0).all()
    assert int(integrity.loc[integrity["Check"] == "Rows", "Observed"].iloc[0]) == 198
