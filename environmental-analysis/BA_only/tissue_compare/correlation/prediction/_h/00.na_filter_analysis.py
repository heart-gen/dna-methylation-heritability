"""
NA Filter Analysis for SDOH Variables (African American cohort)
================================================================
Computes individual-level missingness for all candidate SDOH variables,
applies categorical recoding (matching clean_pheno() in 01.corr_pheno.R),
constructs a composite any_trauma_hx variable, then flags variables
exceeding a >15% missingness threshold.

Output: drfe_results/na_filter_summary.tsv
  Columns: variable, pct_missing, n_valid, min_class_n, class_counts,
           excluded, reason

Run from: environmental-analysis/BA_only/tissue_compare/
          correlation/prediction/_m/
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import session_info
from pyhere import here

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PHENO_PATH = Path(here("inputs/phenotypes/_m/phenotypes-AA.tsv"))
OUTDIR     = Path("drfe_results")

NA_THRESHOLD = 0.15   # >15% individual-level missing → exclude

# All SDOH candidates before recoding (individual trauma vars + composite)
SDOH_RAW = [
    "smoking", "codeine", "morphine", "cocaine", "ethanol",
    "antipsychotics", "nicotine", "amphetamines",
    "education", "marital_status",
    "hx_sexual_abuse", "hx_physical_abuse", "hx_other_trauma",
    "hx_military_service", "fsiq",
]

BOOL_COLS = [
    "smoking", "codeine", "morphine", "cocaine",
    "ethanol", "antipsychotics", "nicotine", "amphetamines",
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def recode_pheno(df: pd.DataFrame) -> pd.DataFrame:
    """Apply same recoding as clean_pheno() in 01.corr_pheno.R."""
    df = df.copy()

    # Boolean substance use: True→1, False→0 (preserve NaN)
    for col in BOOL_COLS:
        if col in df.columns:
            df[col] = df[col].map({True: 1, False: 0, "True": 1, "False": 0})

    # Education → 3-level ordinal
    edu_map = {
        "7th": "less_than_hs", "8th": "less_than_hs",
        "Less than 7th": "less_than_hs",
        "9th": "less_than_hs", "10th": "less_than_hs",
        "11th": "less_than_hs", "12th": "less_than_hs",
        "H.S. diploma": "hs", "GED": "hs",
        "1 yr college": "more_than_hs",
        "3 yrs college": "more_than_hs",
        "Associate's or 2 yrs college": "more_than_hs",
        "Bachelor's": "more_than_hs",
        "Master's": "more_than_hs",
        "JD": "more_than_hs",
        "PhD": "more_than_hs",
    }
    df["education"] = df["education"].map(edu_map)

    # Marital status → 3-level
    marital_map = {
        "Single": "single",
        "Married": "married",
        "Divorced": "previously_married",
        "Separated": "previously_married",
        "Widowed": "previously_married",
    }
    df["marital_status"] = df["marital_status"].map(marital_map)

    # Trauma composite: 1 if ANY trauma True, 0 if ALL available False,
    # NA only when ALL four are missing
    trauma_cols = [
        "hx_sexual_abuse", "hx_physical_abuse",
        "hx_other_trauma", "hx_military_service",
    ]
    trauma = df[trauma_cols].copy()
    # Convert to numeric bool (True/False/nan)
    for c in trauma_cols:
        trauma[c] = trauma[c].map(
            {True: 1, False: 0, "True": 1, "False": 0,
             1: 1, 0: 0}
        ).astype("Float64")

    all_missing = trauma.isna().all(axis=1)
    any_true    = (trauma == 1).any(axis=1)

    df["any_trauma_hx"] = np.where(
        all_missing, np.nan,
        np.where(any_true, 1, 0)
    )

    return df


def compute_stats(series: pd.Series, var_name: str) -> dict:
    """Return missingness stats for a single (recoded) variable series."""
    n_total   = len(series)
    n_missing = series.isna().sum()
    pct_miss  = n_missing / n_total
    n_valid   = n_total - n_missing

    valid = series.dropna()
    class_counts = valid.value_counts().to_dict()
    # Convert keys to str for JSON serialisation
    class_counts = {str(k): int(v) for k, v in class_counts.items()}
    min_class_n  = min(class_counts.values()) if class_counts else 0

    excluded = pct_miss > NA_THRESHOLD
    reason   = (
        f"{pct_miss:.1%} missing (>{NA_THRESHOLD:.0%} threshold)"
        if excluded else "passes threshold"
    )

    return {
        "variable":    var_name,
        "pct_missing": round(pct_miss, 4),
        "n_valid":     int(n_valid),
        "min_class_n": int(min_class_n),
        "class_counts": json.dumps(class_counts),
        "excluded":    excluded,
        "reason":      reason,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    OUTDIR.mkdir(exist_ok=True)

    print("Loading phenotype data (African American cohort)...")
    pheno = pd.read_csv(PHENO_PATH, sep="\t")

    # Deduplicate: one row per individual (drop duplicate brnums)
    pheno = pheno.drop_duplicates(subset="brnum")
    n_total = len(pheno)
    print(f"  {n_total} unique individuals")

    # Apply recoding
    pheno = recode_pheno(pheno)

    # Build list of variables to evaluate:
    # replace 4 individual trauma vars with composite
    eval_vars = [v for v in SDOH_RAW if v not in
                 ("hx_sexual_abuse", "hx_physical_abuse",
                  "hx_other_trauma", "hx_military_service")]
    eval_vars.append("any_trauma_hx")

    rows = []
    print(f"\n{'Variable':<20} {'% Missing':>10} {'N Valid':>8} "
          f"{'Min Class N':>12} {'Excluded':>9}")
    print("-" * 65)

    for var in eval_vars:
        if var not in pheno.columns:
            print(f"  WARNING: '{var}' not found in phenotype file — skipping")
            continue
        stats = compute_stats(pheno[var], var)
        rows.append(stats)
        flag = "EXCLUDED" if stats["excluded"] else "keep"
        print(f"  {var:<20} {stats['pct_missing']:>9.1%} "
              f"{stats['n_valid']:>8d} {stats['min_class_n']:>12d}  {flag}")

    summary = pd.DataFrame(rows)
    out_path = OUTDIR / "na_filter_summary.tsv"
    summary.to_csv(out_path, sep="\t", index=False)
    print(f"\nNA filter summary saved to: {out_path}")
    print(f"  Variables kept:    {(~summary['excluded']).sum()}")
    print(f"  Variables excluded: {summary['excluded'].sum()}")
    excluded_vars = summary.loc[summary["excluded"], "variable"].tolist()
    print(f"  Excluded: {excluded_vars}")

    session_info.show()


if __name__ == "__main__":
    main()