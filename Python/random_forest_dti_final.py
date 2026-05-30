"""
=============================================================================
Random Forest Sensitivity Analysis: Vascular Burden Score + Nigral DTI
=============================================================================

Project: PPMI Vascular Burden Score Analysis (Belnavis et al., Movement
         Disorders 2026)

This script reproduces the manuscript's random forest sensitivity analysis
(Supplementary Analysis S2). The analysis tests whether Vascular Burden Score
components contribute classifier-relevant information beyond candidate nigral
diffusion microstructural biomarkers in the subset of PPMI participants with
baseline DTI data.

-----------------------------------------------------------------------------
DTI metric derivation
-----------------------------------------------------------------------------
Per standard DTI reporting conventions, the random forest uses three derived
diffusion metrics rather than raw tensor eigenvalues:

  FA (fractional anisotropy) - already provided by the PPMI CIND pipeline
  MD (mean diffusivity)      = (L1 + L2 + L3) / 3
  RD (radial diffusivity)    = (L2 + L3) / 2

FA and MD formulas follow the PPMI CIND pre-processing pipeline definitions
(Schuff 2011). RD is computed per standard DTI literature conventions (Basser
and Pierpaoli 1996). The raw eigenvalues (L1, L2, L3) are used only to derive
MD and RD and are not used as features in the classifier; axial diffusivity
(AD = L1) is also not used as a feature given the complex multidirectional
fiber architecture of the substantia nigra.

-----------------------------------------------------------------------------
Design
-----------------------------------------------------------------------------
- Sample:        N = 50 (22 healthy controls, 28 clinically diagnosed PD)
                 Prodromal participants had no baseline DTI in the
                 February 2026 PPMI release.
- Target:        Binary HC vs PD classification.
- Features:      Five total features in the full model:
                   - VBS_micro     (Vascular Burden microvascular sub-score, 0-4)
                   - VBS_macro     (Vascular Burden macrovascular sub-score, 0-3)
                   - Nigral FA     (substantia nigra fractional anisotropy)
                   - Nigral MD     (substantia nigra mean diffusivity)
                   - Nigral RD     (substantia nigra radial diffusivity)
- Classifier:    sklearn RandomForestClassifier
                 - 500 trees
                 - balanced class weights
                 - 5-fold stratified cross-validation
                 - default hyperparameters otherwise
- Metrics:       AUC, balanced accuracy, Cohen's kappa
- Variable
  importance:    permutation importance, 200 iterations per feature
- Feature
  subsets:       Three nested models computed for incremental contribution:
                   1. Full (5 features)
                   2. DTI-only (FA + MD + RD)
                   3. VBS-only (micro + macro sub-scores)

-----------------------------------------------------------------------------
Caveats and interpretation
-----------------------------------------------------------------------------
This analysis is reported in the manuscript as exploratory given the small
imaging-phenotyped subsample (n = 50). The objective is not to demonstrate
strong classification performance but to test whether VBS components capture
variance distinct from nigral DTI microstructure. A weak overall AUC with
non-zero permutation importance for VBS features is consistent with the
orthogonality claim made in the manuscript.

NOTE ON RANDOM SEEDS: With n = 50 and 5-fold CV (10 per held-out fold),
AUC point estimates are sensitive to seed and fold assignment. The qualitative
pattern (DTI features dominating permutation importance, VBS components
contributing measurable additional importance, modest overall classification)
is robust across seeds. Setting a different RANDOM_SEED yields slightly
different point estimates but consistent feature-importance ordering.

-----------------------------------------------------------------------------
Dependencies
-----------------------------------------------------------------------------
Python 3.11+, pandas, numpy, scikit-learn, openpyxl
=============================================================================
"""

import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, cohen_kappa_score


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
ANALYTIC_FILE = Path("F31_PPMI_Vascular_SPSS.xlsx")
DTI_FILE = Path("DTI_Regions_of_Interest_14Jan2026.csv")
OUTPUT_FILE = Path("rf_dti_sensitivity_results.csv")

RANDOM_SEED = 42
N_TREES = 500
N_FOLDS = 5
N_PERMUTATIONS = 200


# -----------------------------------------------------------------------------
# DTI data prep: derive FA, MD, RD from PPMI's per-ROI tensor output
# -----------------------------------------------------------------------------
def build_dti_features(dti_long):
    """Reshape PPMI DTI ROI file and derive standard tensor metrics.

    PPMI's DTI file provides one row per participant per measure (E1, E2, E3, FA)
    with six SN ROIs as columns. This function:

    1. Filters to the substantia nigra rows
    2. Averages across the six SN sub-regions per participant per measure
    3. Renames the three eigenvalues to L1, L2, L3
    4. Computes mean diffusivity   MD = (L1 + L2 + L3) / 3
    5. Computes radial diffusivity RD = (L2 + L3) / 2
    6. Returns per-participant FA, MD, RD (raw eigenvalues are discarded)

    Returns DataFrame with columns: PATNO, FA, MD, RD
    """
    sn = dti_long[dti_long["Tissue"] == "SN"].copy()
    sn["roi_mean"] = sn[["ROI1", "ROI2", "ROI3", "ROI4", "ROI5", "ROI6"]].mean(axis=1)
    wide = sn.pivot_table(index="PATNO", columns="Measure", values="roi_mean", aggfunc="mean")
    wide = wide.rename(columns={"E1": "L1", "E2": "L2", "E3": "L3"})

    # Derive standard DTI metrics from eigenvalues
    wide["MD"] = (wide["L1"] + wide["L2"] + wide["L3"]) / 3
    wide["RD"] = (wide["L2"] + wide["L3"]) / 2

    # Return only the derived features (FA, MD, RD); discard raw eigenvalues
    return wide[["FA", "MD", "RD"]].reset_index()


# -----------------------------------------------------------------------------
# Model fitting helpers
# -----------------------------------------------------------------------------
def fit_and_evaluate(X, y, feature_set_label, seed=RANDOM_SEED):
    """Fit RF and return cross-validated metrics + permutation importance."""
    rf = RandomForestClassifier(
        n_estimators=N_TREES,
        class_weight="balanced",
        random_state=seed,
        n_jobs=-1,
    )
    cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
    y_proba = cross_val_predict(rf, X, y, cv=cv, method="predict_proba", n_jobs=-1)[:, 1]
    y_pred = cross_val_predict(rf, X, y, cv=cv, method="predict", n_jobs=-1)
    auc = roc_auc_score(y, y_proba)
    bal_acc = balanced_accuracy_score(y, y_pred)
    kappa = cohen_kappa_score(y, y_pred)
    rf.fit(X, y)
    perm = permutation_importance(
        rf, X, y, n_repeats=N_PERMUTATIONS, random_state=seed, n_jobs=-1
    )
    importances = dict(zip(X.columns, perm.importances_mean))
    return {
        "feature_set": feature_set_label,
        "n_features": X.shape[1],
        "auc": auc,
        "balanced_accuracy": bal_acc,
        "cohens_kappa": kappa,
        "permutation_importance": importances,
    }


# -----------------------------------------------------------------------------
# Main pipeline
# -----------------------------------------------------------------------------
def main():
    analytic = pd.read_excel(ANALYTIC_FILE)
    print(f"Loaded analytic file: N = {len(analytic)}")

    dti_long = pd.read_csv(DTI_FILE)
    dti = build_dti_features(dti_long)
    print(f"DTI-phenotyped participants (with derived FA, MD, RD): N = {len(dti)}")

    df = analytic.merge(dti, on="PATNO", how="inner")
    df = df[df["Cohort_n"].isin([0, 2])].copy()
    df["target"] = (df["Cohort_n"] == 2).astype(int)
    n_hc = (df["target"] == 0).sum()
    n_pd = (df["target"] == 1).sum()
    print(f"Modeling subset: N = {len(df)} ({n_hc} HC, {n_pd} PD)\n")

    # Feature sets using derived DTI metrics
    feature_sets = {
        "full": ["VBS_micro", "VBS_macro", "FA", "MD", "RD"],
        "dti_only": ["FA", "MD", "RD"],
        "vbs_only": ["VBS_micro", "VBS_macro"],
    }

    results = []
    for label, features in feature_sets.items():
        X = df[features].astype(float)
        y = df["target"].astype(int)
        result = fit_and_evaluate(X, y, label)
        results.append(result)
        print(f"=== {label} (n_features = {len(features)}) ===")
        print(f"  AUC:                  {result['auc']:.3f}")
        print(f"  Balanced accuracy:    {result['balanced_accuracy']:.3f}")
        print(f"  Cohen's kappa:        {result['cohens_kappa']:.3f}")
        print(f"  Permutation importance (mean):")
        for feat, imp in sorted(result["permutation_importance"].items(), key=lambda x: -x[1]):
            print(f"    {feat:<12}: {imp:+.4f}")
        print()

    flat = []
    for r in results:
        for feat, imp in r["permutation_importance"].items():
            flat.append({
                "feature_set": r["feature_set"],
                "n_features": r["n_features"],
                "auc": r["auc"],
                "balanced_accuracy": r["balanced_accuracy"],
                "cohens_kappa": r["cohens_kappa"],
                "feature": feat,
                "permutation_importance_mean": imp,
            })
    pd.DataFrame(flat).to_csv(OUTPUT_FILE, index=False, float_format="%.4f")
    print(f"Results written to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
