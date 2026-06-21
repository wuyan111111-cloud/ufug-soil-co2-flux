"""Reproduce the paired-effect analysis reported in the UFUG manuscript.

The script accepts either the archived Chinese-header Excel workbook or the
public English-header CSV file. It validates the complete 198-record design,
constructs 132 within-stratum log flux ratios, and writes every statistical
table used by the manuscript and Supplementary Information.

Example
-------
python src/reanalysis.py path/to/original_data.xlsx --outdir results
"""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.formula.api as smf
from scipy import stats


CHINESE_COLUMNS = {
    "窗口": "Window",
    "CSU代码": "CSU",
    "既有植物组成背景": "CSU_label",
    "处理": "Treatment",
    "区组": "Block",
    "永久性collar ID": "Collar",
    "CO₂通量(mg CO₂ m⁻² h⁻¹)": "Flux",
}
PUBLIC_COLUMNS = ["Window", "CSU", "CSU_label", "Treatment", "Block", "Collar", "Flux"]
PUBLIC_OUTPUT_COLUMNS = [
    "Window", "CSU", "Plant_background", "Treatment", "Block", "Collar", "Flux"
]
WINDOW_ORDER = ["July", "October"]
TREATMENT_ORDER = ["CK", "DCD", "Biochar"]
EFFECT_TREATMENTS = ["DCD", "Biochar"]
FORMULA = "LogRatio ~ 0 + C(Treatment):C(Window)"
PLANT_BACKGROUNDS = {
    "Lawn": "Managed lawn (turf species not recorded)",
    "PB": "Pinus bungeana",
    "PG": "Phyllostachys glauca",
    "SJ": "Styphnolobium japonicum",
    "NZE": "Ligustrum lucidum × Iris tectorum",
    "GB×CH": "Ginkgo biloba × Cotoneaster horizontalis",
    "SJ×JN": "Styphnolobium japonicum × Jasminum nudiflorum",
    "OF×RC": "Osmanthus fragrans × Rosa chinensis",
    "AR×ND": "Acer palmatum × Nandina domestica",
    "GYL×EF": "Magnolia grandiflora × Euonymus fortunei",
    "GYL×ND": "Magnolia grandiflora × Nandina domestica",
}


def load_data(path: Path) -> pd.DataFrame:
    """Load and validate the archived Excel workbook or public CSV."""
    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path)
    elif path.suffix.lower() in {".csv", ".txt"}:
        df = pd.read_csv(path)
    else:
        raise ValueError("Input must be .xlsx, .xls, or .csv")

    if set(CHINESE_COLUMNS).issubset(df.columns):
        df = df.rename(columns=CHINESE_COLUMNS)
    if "Plant_background" in df.columns and "CSU_label" not in df.columns:
        df["CSU_label"] = df["Plant_background"]
    if not set(PUBLIC_COLUMNS).issubset(df.columns):
        raise ValueError(f"Required columns not found. Observed columns: {list(df.columns)}")
    df = df[PUBLIC_COLUMNS].copy()

    if len(df) != 198:
        raise ValueError(f"Expected 198 observations, found {len(df)}")
    if df[PUBLIC_COLUMNS].isna().any().any():
        raise ValueError("Required fields contain missing values")
    if df.duplicated(["Window", "Collar"]).any():
        raise ValueError("Duplicate Collar × window record detected")
    if (pd.to_numeric(df["Flux"], errors="coerce") <= 0).any():
        raise ValueError("Flux must be positive")
    if set(df["Window"]) != set(WINDOW_ORDER):
        raise ValueError(f"Unexpected monitoring windows: {sorted(df['Window'].unique())}")
    if set(df["Treatment"]) != set(TREATMENT_ORDER):
        raise ValueError(f"Unexpected treatments: {sorted(df['Treatment'].unique())}")
    if df["CSU"].nunique() != 11 or df["Collar"].nunique() != 99:
        raise ValueError("Expected 11 CSUs and 99 permanent collars")

    expected_per_cell = (
        df.groupby(["Window", "CSU", "Treatment"], observed=True).size().to_numpy()
    )
    if not np.all(expected_per_cell == 3):
        raise ValueError("Each Window × CSU × Treatment cell must contain three blocks")

    df["Window"] = pd.Categorical(df["Window"], WINDOW_ORDER, ordered=True)
    df["Treatment"] = pd.Categorical(df["Treatment"], TREATMENT_ORDER, ordered=True)
    df["Block"] = df["Block"].astype(str)
    df["BlockID"] = df["CSU"].astype(str) + ":B" + df["Block"]
    df["Plant_background"] = df["CSU"].map(PLANT_BACKGROUNDS)
    if df["Plant_background"].isna().any():
        raise ValueError("A CSU is missing from the Latin-name mapping")
    return df


def make_ratios(df: pd.DataFrame) -> pd.DataFrame:
    """Create treatment/CK log ratios within CSU × block × window strata."""
    wide = (
        df.pivot(
            index=["CSU", "Plant_background", "Block", "BlockID", "Window"],
            columns="Treatment",
            values="Flux",
        )
        .reset_index()
    )
    rows: list[dict[str, object]] = []
    for record in wide.itertuples(index=False):
        for treatment in EFFECT_TREATMENTS:
            treatment_flux = float(getattr(record, treatment))
            ck_flux = float(record.CK)
            rows.append(
                {
                    "CSU": str(record.CSU),
                    "Plant_background": str(record.Plant_background),
                    "Block": str(record.Block),
                    "BlockID": str(record.BlockID),
                    "Window": str(record.Window),
                    "Treatment": treatment,
                    "TreatmentFlux": treatment_flux,
                    "CKFlux": ck_flux,
                    "LogRatio": float(np.log(treatment_flux / ck_flux)),
                }
            )
    ratios = pd.DataFrame(rows)
    ratios["Window"] = pd.Categorical(ratios["Window"], WINDOW_ORDER, ordered=True)
    ratios["Treatment"] = pd.Categorical(
        ratios["Treatment"], EFFECT_TREATMENTS, ordered=True
    )
    if len(ratios) != 132:
        raise ValueError(f"Expected 132 paired effects, found {len(ratios)}")
    return ratios


def fit_clustered(ratios: pd.DataFrame, cluster: str, df_t: int):
    result = smf.ols(FORMULA, data=ratios).fit(
        cov_type="cluster",
        cov_kwds={"groups": ratios[cluster], "use_correction": True},
        use_t=True,
    )
    return result, df_t


def contrast(result, terms: dict[str, float], df_t: int) -> dict[str, float]:
    names = list(result.params.index)
    vector = np.zeros(len(names))
    for name, weight in terms.items():
        vector[names.index(name)] = weight
    estimate = float(vector @ result.params.to_numpy())
    se = float(np.sqrt(vector @ result.cov_params().to_numpy() @ vector))
    critical = float(stats.t.ppf(0.975, df_t))
    lower, upper = estimate - critical * se, estimate + critical * se
    return {
        "log_ratio": estimate,
        "se": se,
        "ratio": float(np.exp(estimate)),
        "percent_change": float(100 * (np.exp(estimate) - 1)),
        "ci_low_percent": float(100 * (np.exp(lower) - 1)),
        "ci_high_percent": float(100 * (np.exp(upper) - 1)),
        "p": float(2 * stats.t.sf(abs(estimate / se), df_t)),
    }


def treatment_window_effects(result, df_t: int) -> pd.DataFrame:
    rows = []
    for window in WINDOW_ORDER:
        for treatment in EFFECT_TREATMENTS:
            term = f"C(Treatment)[{treatment}]:C(Window)[{window}]"
            row = contrast(result, {term: 1.0}, df_t)
            row.update(Window=window, Treatment=treatment)
            rows.append(row)
    return pd.DataFrame(rows)


def window_contrasts(result, df_t: int) -> pd.DataFrame:
    rows = []
    for treatment in EFFECT_TREATMENTS:
        july = f"C(Treatment)[{treatment}]:C(Window)[July]"
        october = f"C(Treatment)[{treatment}]:C(Window)[October]"
        row = contrast(result, {october: 1.0, july: -1.0}, df_t)
        row.update(Comparison="October minus July", Treatment=treatment)
        rows.append(row)
    return pd.DataFrame(rows)


def csu_window_effects(ratios: pd.DataFrame) -> pd.DataFrame:
    """Descriptive CSU-specific point estimates for each monitoring window."""
    grouped = (
        ratios.groupby(["CSU", "Plant_background", "Treatment", "Window"], observed=True)[
            "LogRatio"
        ]
        .mean()
        .reset_index()
    )
    grouped["percent_change"] = 100 * (np.exp(grouped["LogRatio"]) - 1)
    return grouped


def exploratory_csu_effects(ratios: pd.DataFrame) -> pd.DataFrame:
    """Across-window CSU summaries based on three independent block values."""
    rows = []
    for treatment in EFFECT_TREATMENTS:
        block = (
            ratios[ratios["Treatment"] == treatment]
            .groupby(["CSU", "Plant_background", "Block"], observed=True)["LogRatio"]
            .mean()
            .reset_index()
        )
        for (csu, label), group in block.groupby(
            ["CSU", "Plant_background"], observed=True
        ):
            values = group["LogRatio"].to_numpy()
            mean = float(values.mean())
            se = float(values.std(ddof=1) / np.sqrt(3))
            critical = float(stats.t.ppf(0.975, 2))
            lower, upper = mean - critical * se, mean + critical * se
            rows.append(
                {
                    "CSU": str(csu),
                    "Plant_background": str(label),
                    "Treatment": treatment,
                    "n_blocks": 3,
                    "log_ratio": mean,
                    "percent_change": float(100 * (np.exp(mean) - 1)),
                    "ci_low_percent": float(100 * (np.exp(lower) - 1)),
                    "ci_high_percent": float(100 * (np.exp(upper) - 1)),
                }
            )
    return pd.DataFrame(rows)


def descriptive_summary(df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        df.groupby(["Window", "Treatment"], observed=True)["Flux"]
        .agg(n="count", mean="mean", sd="std", median="median", minimum="min", maximum="max")
        .reset_index()
    )
    return summary


def integrity_checks(df: pd.DataFrame) -> pd.DataFrame:
    records = [
        ("Rows", len(df), 198),
        ("Unique Collars", df["Collar"].nunique(), 99),
        ("Unique CSUs", df["CSU"].nunique(), 11),
        ("Unique blocks", df["BlockID"].nunique(), 33),
        ("Missing required fields", int(df[PUBLIC_COLUMNS].isna().sum().sum()), 0),
        ("Duplicate Collar × window", int(df.duplicated(["Window", "Collar"]).sum()), 0),
        ("Minimum flux", float(df["Flux"].min()), "positive"),
        ("Maximum flux", float(df["Flux"].max()), "observed range"),
    ]
    return pd.DataFrame(records, columns=["Check", "Observed", "Criterion"])


def write_outputs(input_path: Path, outdir: Path) -> dict[str, object]:
    outdir.mkdir(parents=True, exist_ok=True)
    df = load_data(input_path)
    ratios = make_ratios(df)
    primary, primary_df = fit_clustered(ratios, "CSU", 10)
    sensitivity, sensitivity_df = fit_clustered(ratios, "BlockID", 32)

    outputs = {
        "collar_flux_198_public.csv": df[PUBLIC_OUTPUT_COLUMNS],
        "paired_effects_132.csv": ratios,
        "descriptive_summary.csv": descriptive_summary(df),
        "primary_csu_clustered_effects.csv": treatment_window_effects(
            primary, primary_df
        ),
        "sensitivity_block_clustered_effects.csv": treatment_window_effects(
            sensitivity, sensitivity_df
        ),
        "window_interaction_tests.csv": window_contrasts(primary, primary_df),
        "csu_window_effects.csv": csu_window_effects(ratios),
        "exploratory_csu_effects.csv": exploratory_csu_effects(ratios),
        "data_integrity_checks.csv": integrity_checks(df),
    }
    for filename, frame in outputs.items():
        frame.to_csv(outdir / filename, index=False, encoding="utf-8-sig")

    summary = {
        "input_file": input_path.name,
        "n_original": len(df),
        "n_collars": int(df["Collar"].nunique()),
        "n_csus": int(df["CSU"].nunique()),
        "n_blocks": int(df["BlockID"].nunique()),
        "n_paired_effects": len(ratios),
        "formula": FORMULA,
        "primary_cluster": "CSU",
        "primary_df": 10,
        "sensitivity_cluster": "BlockID",
        "sensitivity_df": 32,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "statsmodels": statsmodels.__version__,
    }
    (outdir / "run_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_data", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("results"))
    args = parser.parse_args()
    summary = write_outputs(args.input_data, args.outdir)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
