import numpy as np
import pandas as pd
import session_info
from pathlib import Path

from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression

from dRFEtools.dRFE import dRFEtools

def load_data(data_path):
    df = pd.read_csv(data_path, sep="\t")
    df = df[df["group"] == "AA"].copy()

    # Encode booleans
    bool_cols = [
        "smoking", "codeine", "morphine", "cocaine",
        "ethanol", "antipsychotics", "nicotine", "amphetamines"
    ]
    for c in bool_cols:
        df[c] = df[c].map({True: 1, False: 0})

    # Education: 3-category
    edu_map = {
        "Less than high school": "less_than_hs",
        "High School": "hs",
        "More than high school": "more_than_hs"
    }
    df["education_3cat"] = df["education"].map(edu_map)

    # Marital status: 3-category
    marital_map = {
        "Single": "single",
        "Married": "married",
        "Divorced": "previously_married",
        "Separated": "previously_married",
        "Widowed": "previously_married"
    }
    df["marital_3cat"] = df["marital_status"].map(marital_map)
    return df


def filter_features(df, frac, var):
    # Pivot to samples × VMRs
    X = df.pivot_table(index="brnum", columns="feature_id", values="meth")
    # Feature filtering
    keep_features = (X.notna().mean() >= frac) & (X.var(skipna=True) >= var)
    return X.loc[:, keep_features]


def filter_target(df, X, yvar):
    y = df.drop_duplicates("brnum")\
          .set_index("brnum")[yvar]\
          .reindex(X.index)
    # Drop missing labels
    valid = y.notna()
    Xy    = X.loc[valid]
    y     = y.loc[valid]
    return Xy, y


def main():
    # Configuration
    DATA_PATH = Path("../../_m/vmr_env_assoc-AA.tsv.gz")
    OUTDIR = Path("drfe_results")
    OUTDIR.mkdir(exist_ok=True)

    MIN_SAMPLES = 40          # minimum samples per task
    MIN_FEATURE_VAR = 1e-5    # variance filter
    FEATURE_MIN_FRAC = 0.6    # feature present in ≥60% samples

    ENV_VARS = [
        "sex", "primarydx", "smoking", "codeine", "morphine",
        "cocaine", "ethanol", "antipsychotics", "nicotine",
        "amphetamines", "education_3cat", "marital_3cat",
    ]

    # Load data
    df = load_data(DATA_PATH)

    # Loop prediction
    results = []
    for region in df["region"].dropna().unique():
        for h2_cat in df["h2_category"].dropna().unique():

            sub = df[(df["region"] == region) & (df["h2_category"] == h2_cat)]
            X = filter_features(sub, FEATURE_MIN_FRAC, MIN_FEATURE_VAR)

            if X.shape[1] < 10:
                continue

            for yvar in ENV_VARS:
                Xy, y = filter_target(sub, X, yvar)

                # Require ≥3 classes
                if y.nunique() < 2:
                    continue

                # Minimum samples per class
                class_counts = y.value_counts()
                if (class_counts < 10).any():
                    continue

                if len(y) < MIN_SAMPLES:
                    continue

                # Set model
                base_estimator = LogisticRegression(
                    penalty="elasticnet", solver="saga", multi_class="multinominal",
                    l1_ratio=0.5, max_iter=5000, class_weight="balanced", n_jobs=1
                )

                pipeline = Pipeline([
                    ("impute", SimpleImputer(strategy="median")),
                    ("scale", StandardScaler()),
                    ("clf", base_estimator)
                ])

                cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=13)
                drfe = dRFEtools(
                    estimator=pipeline, cv=cv, scoring="roc_auc_ovr",
                    min_features_to_select=5, step=0.2, n_jobs=5
                )

                try:
                    drfe.fit(Xy.values, y.values)

                    best_auc = max(drfe.cv_scores_)
                    best_nfeat = drfe.n_features_[np.argmax(drfe.cv_scores_)]

                    results.append({
                        "region": region, "h2_category": h2_cat,
                        "env_var": yvar, "n_samples": len(y),
                        "n_features": X.shape[1], "best_auc": best_auc,
                        "best_n_features": best_nfeat
                    })

                    # Save selected features
                    selected = X.columns[drfe.support_]
                    pd.Series(selected).to_csv(
                        OUTDIR / f"features_{region}_{h2_cat}_{yvar}.txt",
                        index=False
                    )
                except Exception as e:
                    print(f"Failed: {region}, {h2_cat}, {yvar}: {e}")

    # Save summary
    res_df = pd.DataFrame(results)
    res_df.to_csv(OUTDIR / "drfe_summary.tsv", sep="\t", index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()

