"""
=============================================================================
PPMI Vascular Burden Score: Modified Poisson Regression for Adjusted Risk Ratios
=============================================================================

Project: Vascular Phenotyping in Parkinson's Disease: Diabetes Mellitus
         Operationalizes a Microvascular Metabolic Syndrome Cluster
         Across the PPMI Clinical Spectrum

Authors: A. Belnavis, S.Y. Chiu, K. Chen, R. Thorpe, E. Ofori
Journal: Movement Disorders (submitted 2026)
PPMI download date: February 12, 2026
PPMI RRID: SCR_006431

-----------------------------------------------------------------------------
Purpose
-----------------------------------------------------------------------------
Per-factor association testing of diabetes mellitus / prediabetes (DM/Pre-DM)
with each of seven vascular risk factor indicators, using modified Poisson
regression with Huber-White robust variance estimation (Zou 2004).

Modified Poisson regression is preferred over standard logistic regression
when the outcome prevalence is non-rare (here ~17.7%), because the odds ratio
overestimates the relative risk in this regime. The Poisson model with a
log link directly estimates adjusted Risk Ratios (aRRs).

Reference for the method:
    Zou G. A modified poisson regression approach to prospective studies
    with binary data. Am J Epidemiol. 2004;159(7):702-706.

-----------------------------------------------------------------------------
Pipeline
-----------------------------------------------------------------------------
1. Load analytic file (N = 413).
2. For each of 7 vascular risk factors, fit GLM with Poisson family + log
   link, with DM_Group + age + sex + BMI as covariates. (BMI omitted for the
   obesity outcome to avoid collinearity, since obesity is BMI >= 30.)
3. Robust (Huber-White / HC0) standard errors.
4. Compute aRR, 95% CI, p-value.
5. Apply Benjamini-Hochberg FDR correction across the 7 factors.
6. Compose results table.

-----------------------------------------------------------------------------
Output
-----------------------------------------------------------------------------
table2_arr_per_factor.csv  -- the 7-row results table that populates manuscript
                              Table 2 / Figure 2 forest plot.

-----------------------------------------------------------------------------
Dependencies
-----------------------------------------------------------------------------
Python 3.11+
pandas >= 2.0
numpy >= 1.24
statsmodels >= 0.14
openpyxl (for .xlsx reading)

-----------------------------------------------------------------------------
Usage
-----------------------------------------------------------------------------
    python poisson_arr_per_factor.py

Input file path is set as a constant below; adjust if necessary.

=============================================================================
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from pathlib import Path


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
INPUT_FILE = Path("F31_PPMI_Vascular_SPSS.xlsx")  # PPMI analytic file (N=413)
OUTPUT_FILE = Path("table2_arr_per_factor.csv")
ALPHA = 0.05

# Per-factor model specifications.
# Each tuple: (outcome_variable, covariate_list, display_label, category)
# BMI is omitted for the obesity outcome because obesity is defined as BMI>=30,
# which would create perfect collinearity with the BMI covariate.
FACTOR_SPECS = [
    ("HTN",       ["DM_Group", "AGE", "Sex_M", "BMI"], "Hypertension",            "microvascular"),
    ("HLD",       ["DM_Group", "AGE", "Sex_M", "BMI"], "Hyperlipidemia",          "microvascular"),
    ("OSA",       ["DM_Group", "AGE", "Sex_M", "BMI"], "Sleep apnea (OSA)",       "microvascular"),
    ("OBESE_BMI", ["DM_Group", "AGE", "Sex_M"],         "Obesity (BMI >= 30)",     "microvascular"),
    ("STROKE",    ["DM_Group", "AGE", "Sex_M", "BMI"], "Stroke / TIA",            "macrovascular"),
    ("CAD",       ["DM_Group", "AGE", "Sex_M", "BMI"], "Coronary artery disease", "macrovascular"),
    ("AFIB",      ["DM_Group", "AGE", "Sex_M", "BMI"], "Atrial fibrillation",     "macrovascular"),
]


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------
def fit_modified_poisson(df, outcome, covariates):
    """
    Fit modified Poisson regression with robust (Huber-White) standard errors.

    The "modified" Poisson approach uses a Poisson GLM with log link and
    HC0 robust variance estimation, which yields valid risk ratios and
    confidence intervals for binary outcomes (Zou 2004).

    Parameters
    ----------
    df : pd.DataFrame
        Analytic dataset, one row per participant.
    outcome : str
        Name of binary outcome column (0/1).
    covariates : list of str
        Names of predictor columns. First element is treated as the exposure
        of interest for extraction of the aRR.

    Returns
    -------
    dict with keys:
        arr        - exp(beta) for the exposure of interest
        ci_low     - lower 95% confidence limit on the aRR
        ci_high    - upper 95% confidence limit on the aRR
        p_value    - Wald p-value for the exposure coefficient
        n_used     - number of observations entering the model (listwise)
        n_events   - number of outcome events
    """
    sub = df[[outcome] + covariates].dropna()
    y = sub[outcome]
    X = sm.add_constant(sub[covariates])

    model = sm.GLM(y, X, family=sm.families.Poisson()).fit(cov_type="HC0")

    exposure = covariates[0]
    beta = model.params[exposure]
    se = model.bse[exposure]
    return {
        "arr": np.exp(beta),
        "ci_low": np.exp(beta - 1.96 * se),
        "ci_high": np.exp(beta + 1.96 * se),
        "p_value": model.pvalues[exposure],
        "n_used": int(sub.shape[0]),
        "n_events": int(y.sum()),
    }


def benjamini_hochberg(p_values, alpha=ALPHA):
    """Apply Benjamini-Hochberg false discovery rate correction.

    Returns the q-values and a boolean significance flag at the given alpha.
    """
    reject, q_values, _, _ = multipletests(p_values, alpha=alpha, method="fdr_bh")
    return q_values, reject


# -----------------------------------------------------------------------------
# Main pipeline
# -----------------------------------------------------------------------------
def main():
    # 1. Load data
    df = pd.read_excel(INPUT_FILE)
    print(f"Loaded {INPUT_FILE}: N = {len(df)} participants")
    print(f"DM/Pre-DM positive: n = {(df['DM_Group']==1).sum()} "
          f"({100*(df['DM_Group']==1).mean():.1f}%)\n")

    # 2-4. Fit each per-factor model
    rows = []
    for outcome, covars, label, category in FACTOR_SPECS:
        result = fit_modified_poisson(df, outcome, covars)
        rows.append({
            "factor": label,
            "category": category,
            "n_used": result["n_used"],
            "n_events": result["n_events"],
            "aRR": result["arr"],
            "ci_low": result["ci_low"],
            "ci_high": result["ci_high"],
            "p_raw": result["p_value"],
        })

    table = pd.DataFrame(rows)

    # 5. FDR correction across the 7 factors
    table["q_BH"], table["fdr_significant"] = benjamini_hochberg(table["p_raw"].values)

    # 6. Format and report
    print("Per-factor adjusted Risk Ratios (modified Poisson, age/sex/BMI adjusted)")
    print("=" * 90)
    print(f"{'Factor':<28} {'aRR':>6}  {'95% CI':<18}  {'p':>8}  {'q (FDR)':>9}  {'sig'}")
    print("-" * 90)
    for _, r in table.iterrows():
        ci_str = f"({r['ci_low']:.2f}, {r['ci_high']:.2f})"
        sig = "*" if r["fdr_significant"] else " "
        print(f"{r['factor']:<28} {r['aRR']:>6.2f}  {ci_str:<18}  "
              f"{r['p_raw']:>8.3f}  {r['q_BH']:>9.3f}   {sig}")

    # 7. Save
    table.to_csv(OUTPUT_FILE, index=False, float_format="%.4f")
    print(f"\nResults written to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
