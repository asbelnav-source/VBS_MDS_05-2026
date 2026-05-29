"""
=============================================================================
Count-Model Sensitivity Analysis for Composite Vascular Burden Score
=============================================================================

Project: PPMI Vascular Burden Score Analysis (Belnavis et al., Movement
         Disorders 2026)

This script addresses an anticipated reviewer concern: the primary multivariable
regression analyses (Models A, B, C) use ordinary least squares (linear)
regression on the composite Vascular Burden Score, which is technically a
bounded count outcome (range 0-7 for total VBS; 0-4 for microvascular sub-score;
0-3 for macrovascular sub-score).

To verify the robustness of the linear-model inference, we refit each of the
three primary models under two alternative count-model specifications and
compare the diabetes mellitus / prediabetes (DM_Group) coefficient estimates.

-----------------------------------------------------------------------------
Models compared
-----------------------------------------------------------------------------
For each of the three outcomes (VBS, VBS_micro, VBS_macro), the following
three specifications are fit and the DM_Group coefficient is extracted:

  1. Ordinary least squares (linear regression) - the primary manuscript method
  2. Poisson regression with HC0 robust variance - count-model with log link
  3. Negative binomial regression (alpha = 1.0) - handles potential
     overdispersion

For Poisson and negative binomial, the coefficient is also reported as an
incidence rate ratio (IRR = exp(beta)) for clinical interpretability.

-----------------------------------------------------------------------------
Covariate specification
-----------------------------------------------------------------------------
Matches the primary manuscript analysis:
  - Model A (total VBS): adjusted for age, sex, body mass index
  - Model B (microvascular): adjusted for age, sex; BMI omitted because
    obesity (BMI >= 30) is a component of the microvascular sub-score
  - Model C (macrovascular): adjusted for age, sex, body mass index

-----------------------------------------------------------------------------
Interpretation
-----------------------------------------------------------------------------
Substantively concordant findings across linear, Poisson, and negative
binomial specifications support the validity of the primary linear-model
inference. Pearson dispersion statistics from the Poisson fits are also
reported; values close to or below 1.0 indicate no meaningful overdispersion,
in which case Poisson and OLS yield substantively equivalent inference.

-----------------------------------------------------------------------------
Dependencies
-----------------------------------------------------------------------------
Python 3.11+, pandas, numpy, statsmodels, openpyxl

-----------------------------------------------------------------------------
Output
-----------------------------------------------------------------------------
count_model_sensitivity_results.csv - tabular comparison of all three
specifications for each of the three models. Used to populate
Supplementary Table S5.
=============================================================================
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
from pathlib import Path


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
INPUT_FILE = Path("F31_PPMI_Vascular_SPSS.xlsx")
OUTPUT_FILE = Path("count_model_sensitivity_results.csv")

# Three model specifications matching the manuscript Models A, B, C.
# Each: (display label, outcome variable, covariate list)
MODEL_SPECS = [
    ("Model A: Total VBS (0-7)",     "VBS",       ["DM_Group", "AGE", "Sex_M", "BMI"]),
    ("Model B: Microvascular (0-4)", "VBS_micro", ["DM_Group", "AGE", "Sex_M"]),
    ("Model C: Macrovascular (0-3)", "VBS_macro", ["DM_Group", "AGE", "Sex_M", "BMI"]),
]

# Negative binomial dispersion parameter alpha; 1.0 is a standard default that
# allows the model to handle overdispersion. statsmodels does not estimate
# alpha by default through the GLM interface.
NB_ALPHA = 1.0


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def extract_coef(model, predictor):
    """Extract coefficient, SE, 95% CI, p-value for a single predictor."""
    beta = model.params[predictor]
    se = model.bse[predictor]
    ci_low = beta - 1.96 * se
    ci_high = beta + 1.96 * se
    p = model.pvalues[predictor]
    return beta, se, ci_low, ci_high, p


def fit_all_three(df, outcome, covars):
    """Fit OLS, Poisson, and negative binomial; return coefficient summaries."""
    sub = df[[outcome] + covars].dropna()
    y = sub[outcome]
    X = sm.add_constant(sub[covars])
    n = len(sub)

    # 1) Linear regression
    ols = sm.OLS(y, X).fit()
    ols_b, ols_se, ols_lo, ols_hi, ols_p = extract_coef(ols, "DM_Group")
    ols_r2 = ols.rsquared

    # 2) Poisson with HC0 robust SE
    pois = sm.GLM(y, X, family=sm.families.Poisson()).fit(cov_type="HC0")
    p_b, p_se, p_lo, p_hi, p_p = extract_coef(pois, "DM_Group")
    overdisp = pois.pearson_chi2 / pois.df_resid

    # 3) Negative binomial
    try:
        nb = sm.GLM(
            y, X, family=sm.families.NegativeBinomial(alpha=NB_ALPHA)
        ).fit()
        n_b, n_se, n_lo, n_hi, n_p = extract_coef(nb, "DM_Group")
    except Exception as e:
        print(f"    Negative binomial failed: {e}")
        n_b = n_se = n_lo = n_hi = n_p = np.nan

    return {
        "n": n,
        "ols_beta": ols_b, "ols_se": ols_se,
        "ols_ci_low": ols_lo, "ols_ci_high": ols_hi,
        "ols_p": ols_p, "ols_R2": ols_r2,
        "poisson_beta": p_b, "poisson_se": p_se,
        "poisson_ci_low": p_lo, "poisson_ci_high": p_hi,
        "poisson_p": p_p,
        "poisson_IRR": np.exp(p_b),
        "poisson_IRR_low": np.exp(p_lo),
        "poisson_IRR_high": np.exp(p_hi),
        "poisson_pearson_dispersion": overdisp,
        "nb_beta": n_b, "nb_se": n_se,
        "nb_ci_low": n_lo, "nb_ci_high": n_hi,
        "nb_p": n_p,
        "nb_IRR": np.exp(n_b) if not np.isnan(n_b) else np.nan,
        "nb_IRR_low": np.exp(n_lo) if not np.isnan(n_lo) else np.nan,
        "nb_IRR_high": np.exp(n_hi) if not np.isnan(n_hi) else np.nan,
    }


# -----------------------------------------------------------------------------
# Main pipeline
# -----------------------------------------------------------------------------
def main():
    df = pd.read_excel(INPUT_FILE)
    print(f"Loaded {INPUT_FILE}: N = {len(df)}\n")

    results = []
    for label, outcome, covars in MODEL_SPECS:
        print(f"=== {label} ===")
        out = fit_all_three(df, outcome, covars)
        out["model"] = label
        results.append(out)

        print(f"  Linear OLS:        beta = {out['ols_beta']:+.3f} "
              f"({out['ols_ci_low']:+.3f}, {out['ols_ci_high']:+.3f}), "
              f"p = {out['ols_p']:.4f}, R² = {out['ols_R2']:.3f}")
        print(f"  Poisson (HC0):     beta = {out['poisson_beta']:+.3f} "
              f"({out['poisson_ci_low']:+.3f}, {out['poisson_ci_high']:+.3f}), "
              f"p = {out['poisson_p']:.4f}; "
              f"IRR = {out['poisson_IRR']:.3f}; "
              f"dispersion = {out['poisson_pearson_dispersion']:.3f}")
        print(f"  Negative binomial: beta = {out['nb_beta']:+.3f} "
              f"({out['nb_ci_low']:+.3f}, {out['nb_ci_high']:+.3f}), "
              f"p = {out['nb_p']:.4f}; "
              f"IRR = {out['nb_IRR']:.3f}")
        print()

    out_df = pd.DataFrame(results)
    cols = ["model"] + [c for c in out_df.columns if c != "model"]
    out_df = out_df[cols]
    out_df.to_csv(OUTPUT_FILE, index=False, float_format="%.4f")
    print(f"Results written to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
