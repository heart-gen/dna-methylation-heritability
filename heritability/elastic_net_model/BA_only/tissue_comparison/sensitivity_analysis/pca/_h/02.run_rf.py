"""
This script calculates heritability classification accuracy 
based on methylation PC values
"""
import pandas as pd
from pathlib import Path
from pyhere import here

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, classification_report
import session_info

def get_pca(tissue):
    fn = here(f"heritability/elastic_net_model/tissue_comparison/pca/_m/{tissue.lower()}_meth_pc.tsv")
    return pd.read_csv(fn, sep='\t')

def run_rf(pca_df):
    # Split data
    X = pca_df.filter(like='PC')
    y = pca_df['h2_category']

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    fold_accuracies = []

    for train_idx, test_idx in cv.split(X, y):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        rf = RandomForestClassifier()
        rf.fit(X_train, y_train)
        y_pred = rf.predict(X_test)

        accuracy = accuracy_score(y_test, y_pred)
        fold_accuracies.append(accuracy)
        print(classification_report(y_test, y_pred))

    print(f'Best fold accuracy: {max(fold_accuracies)}')


def main():
    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
            
            # Read in pcs 
            pc_df = get_pca(tissue)

            # Run random forest
            run_rf(pc_df)

    # Session information
    session_info.show()

if __name__ == "__main__":
    main()