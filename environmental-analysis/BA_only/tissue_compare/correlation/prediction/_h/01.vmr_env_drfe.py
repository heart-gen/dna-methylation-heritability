import os
import re
import numpy as np
import pandas as pd
import session_info
from pathlib import Path

from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression

import dRFEtools

def make_file_safe(s):
    """Convert string to file-safe format: lowercase, spaces to underscores."""
    s = str(s).lower().strip()
    s = re.sub(r'\s+', '_', s)  # Replace whitespace with underscores
    s = re.sub(r'[^\w\-]', '', s)  # Remove non-alphanumeric chars except _ and -
    return s


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


def run_drfe_cv(X, y, features, outdir, n_splits=5, elimination_rate=0.2):
    """
    Run dRFE with cross-validation using the official dRFEtools API.

    Returns:
        fold_results: dict of elimination results per fold
        fold_importances: dict of feature importances per fold
    """
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=13)
    fold_results = {}
    fold_importances = {}

    for fold, (train_idx, _) in enumerate(cv.split(X, y)):
        X_train = X[train_idx]
        y_train = y[train_idx]

        # Preprocess data
        imputer = SimpleImputer(strategy="median")
        scaler = StandardScaler()
        X_processed = scaler.fit_transform(imputer.fit_transform(X_train))

        # Set model (multinomial is default in sklearn 1.7+)
        estimator = LogisticRegression(
            penalty="elasticnet", solver="saga",
            l1_ratio=0.5, max_iter=5000, class_weight="balanced"
        )

        # Run dRFE for this fold with ranking enabled
        results, first_pass = dRFEtools.dev_rfe(
            estimator=estimator,
            X=X_processed,
            Y=y_train,
            features=features,
            fold=fold,
            out_dir=outdir,
            elimination_rate=elimination_rate,
            dev_size=0.2,
            RANK=True,
            SEED=True
        )

        fold_results[fold] = results

        # Extract feature importances by refitting on full feature set
        # dRFEtools results don't include importances directly
        estimator.fit(X_processed, y_train)
        # Use L1 norm across classes for multi-class, abs for binary
        coef = estimator.coef_
        if coef.ndim > 1:
            importances = np.linalg.norm(coef, axis=0, ord=1)
        else:
            importances = np.abs(coef).flatten()

        fold_importances[fold] = importances

    return fold_results, fold_importances


def rank_features_across_folds(features, fold_importances):
    """
    Aggregate feature importances across CV folds and create rankings.

    Returns DataFrame with feature rankings for downstream biological analysis.
    """
    n_features = len(features)
    n_folds = len(fold_importances)

    # Initialize importance matrix (folds × features)
    importance_matrix = np.zeros((n_folds, n_features))

    for fold, importances in fold_importances.items():
        if len(importances) == n_features:
            importance_matrix[fold] = importances

    # Aggregate across folds
    mean_importance = np.mean(importance_matrix, axis=0)
    std_importance = np.std(importance_matrix, axis=0)

    # Compute rank per fold, then average ranks
    rank_matrix = np.zeros_like(importance_matrix)
    for fold in range(n_folds):
        # Rank by descending importance (rank 1 = most important)
        rank_matrix[fold] = n_features - np.argsort(np.argsort(importance_matrix[fold]))

    mean_rank = np.mean(rank_matrix, axis=0)

    # Create results DataFrame
    rank_df = pd.DataFrame({
        "feature": features,
        "mean_importance": mean_importance,
        "std_importance": std_importance,
        "mean_rank": mean_rank,
        "min_rank": np.min(rank_matrix, axis=0),
        "max_rank": np.max(rank_matrix, axis=0),
    })

    # Sort by mean importance (descending)
    rank_df = rank_df.sort_values("mean_importance", ascending=False)
    rank_df["final_rank"] = range(1, len(rank_df) + 1)

    return rank_df


def aggregate_fold_results(fold_results):
    """
    Aggregate dRFE results across folds.

    Returns best metric score and optimal feature count.
    """
    # Collect scores across folds for each feature count
    n_features_scores = {}

    for fold, results in fold_results.items():
        for n_feat, result_dict in results.items():
            if n_feat not in n_features_scores:
                n_features_scores[n_feat] = []
            # Use ROC AUC for classification (metrics is nested)
            metrics = result_dict.get("metrics", {})
            score = metrics.get("roc_auc_score", metrics.get("accuracy_score", 0))
            n_features_scores[n_feat].append(score)

    # Average scores across folds
    avg_scores = {
        n: np.mean(scores) for n, scores in n_features_scores.items()
    }

    # Find best feature count using LOWESS on first fold
    first_fold_results = fold_results[0]
    try:
        best_n_features, _ = dRFEtools.extract_max_lowess(
            first_fold_results, frac=0.3, multi=True
        )
    except Exception:
        # Fallback to max score if LOWESS fails
        best_n_features = max(avg_scores, key=avg_scores.get)

    best_score = avg_scores.get(best_n_features, max(avg_scores.values()))

    return best_score, best_n_features


def main():
    # Configuration
    OUTDIR = Path("drfe_results")
    OUTDIR.mkdir(exist_ok=True)

    MIN_SAMPLES = 40          # minimum samples per task
    MIN_FEATURE_VAR = 1e-5    # variance filter
    FEATURE_MIN_FRAC = 0.6    # feature present in ≥60% samples

    na_filter = pd.read_csv(OUTDIR / "na_filter_summary.tsv", sep="\t")
    ENV_VARS = na_filter.loc[~na_filter["excluded"], "variable"].tolist()

    # Get unique regions (sorted for reproducibility)
    all_regions = ["Caudate", "Hippocampus", "DLPFC"]
    h2_categories = ["Heritable", "Non-heritable", "Low Prediction"]

    # SLURM array job support: process single region if SLURM_ARRAY_TASK_ID is set
    slurm_task_id = os.environ.get("SLURM_ARRAY_TASK_ID")
    if slurm_task_id is not None:
        task_idx = int(slurm_task_id)
        if task_idx >= len(all_regions):
            print(f"SLURM_ARRAY_TASK_ID={task_idx} exceeds number of regions ({len(all_regions)})")
            return
        regions_to_process = [all_regions[task_idx]]
        print(f"SLURM array job: processing region index {task_idx} -> {regions_to_process[0]}")
    else:
        regions_to_process = all_regions
        print("Processing all regions (no SLURM_ARRAY_TASK_ID set)")

    # Loop prediction
    results = []
    task_count = 0
    for region in regions_to_process:
        for h2_cat in h2_categories:
            print(f"\nProcessing: {region} / {h2_cat}")

            # Load data
            DATA_PATH = Path(f"../../../../{region.lower()}/correlation/_m/vmr_env_assoc-AA.tsv.gz")
            print("Loading data...")
            df = load_data(DATA_PATH)
            print(f"Loaded {len(df)} rows")

            sub = df[(df["h2_category"] == h2_cat)]
            X = filter_features(sub, FEATURE_MIN_FRAC, MIN_FEATURE_VAR)

            if X.shape[1] < 10:
                continue

            for yvar in ENV_VARS:
                Xy, y = filter_target(sub, X, yvar)

                # Require ≥2 classes
                if y.nunique() < 2:
                    continue

                # Minimum samples per class
                class_counts = y.value_counts()
                if (class_counts < 10).any():
                    continue

                if len(y) < MIN_SAMPLES:
                    continue

                # Create task-specific output directory (file-safe h2_category)
                h2_cat_safe = make_file_safe(h2_cat)
                task_outdir = OUTDIR / f"{region}_{h2_cat_safe}_{yvar}"
                task_outdir.mkdir(exist_ok=True)

                features = np.array(Xy.columns.tolist())

                try:
                    # Run cross-validated dRFE
                    fold_results, fold_importances = run_drfe_cv(
                        Xy.values, y.values, features, task_outdir,
                        n_splits=5, elimination_rate=0.2
                    )

                    # Aggregate results
                    best_score, best_nfeat = aggregate_fold_results(fold_results)

                    # Rank features across folds for biological analysis
                    rank_df = rank_features_across_folds(features, fold_importances)
                    rank_df.to_csv(
                        task_outdir / "ranked_features.tsv",
                        sep="\t", index=False
                    )

                    results.append({
                        "region": region, "h2_category": h2_cat,
                        "h2_category_safe": h2_cat_safe,
                        "env_var": yvar, "n_samples": len(y),
                        "n_features": X.shape[1], "best_score": best_score,
                        "best_n_features": best_nfeat
                    })

                    # Save selected features (top N based on ranking)
                    selected = rank_df.head(best_nfeat)["feature"].tolist()
                    pd.Series(selected).to_csv(
                        task_outdir / "selected_features.txt",
                        index=False
                    )

                    task_count += 1
                    print(f"  [{task_count}] {yvar}: score={best_score:.3f}, "
                          f"n_features={best_nfeat}, n_samples={len(y)}")

                except Exception as e:
                    print(f"  [FAILED] {yvar}: {e}")

    # Save summary (region-specific filename for SLURM array jobs)
    res_df = pd.DataFrame(results)
    if slurm_task_id is not None:
        summary_file = OUTDIR / f"drfe_summary_{regions_to_process[0]}.tsv"
    else:
        summary_file = OUTDIR / "drfe_summary.tsv"
    res_df.to_csv(summary_file, sep="\t", index=False)

    print(f"\nCompleted {task_count} tasks")
    print(f"Results saved to: {summary_file}")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()

