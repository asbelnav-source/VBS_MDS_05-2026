# Python Analysis Pipeline

Reproducible Python scripts for the modified Poisson regression, sensitivity analyses, and FDR correction reported in:

Belnavis A, Chiu SY, Chen K, Thorpe R, Ofori E. Vascular Phenotyping in Parkinson's Disease: Diabetes Mellitus Operationalizes a Microvascular Metabolic Syndrome Cluster Across the PPMI Clinical Spectrum. Movement Disorders. 2026.

## Contents

| File | Purpose | Reproduces |
|---|---|---|
| `poisson_arr_per_factor.py` | Modified Poisson regression for each of seven vascular risk factors with Huber-White robust variance | Table 2; Figure 2 forest plot |
| `fdr_correction.py` | Benjamini-Hochberg false discovery rate correction utility (importable module + CLI) | FDR-adjusted q-values across the seven per-factor analyses |
| `random_forest_dti.py` | Random forest classification of HC vs PD using VBS sub-scores plus nigral DTI eigenvalues | Figure 3B / 4B sensitivity panel |
| `wmh_sensitivity.py` | Linear regression of log-transformed FreeSurfer 7 white matter hypointensity volume on VBS | Supplementary Analysis S1 |
| `requirements.txt` | Pinned Python dependencies | n/a |

## Dependencies

Python 3.11 or later. Install all required packages with:

​```bash
pip install -r requirements.txt
​```

Core dependencies:

- pandas (data manipulation)
- numpy (numerical operations)
- statsmodels (modified Poisson GLM, FDR correction, OLS)
- scikit-learn (random forest, cross-validation, permutation importance)
- openpyxl (reads the analytic .xlsx file)
- matplotlib (optional, only for figure regeneration)

## Required input files

All scripts assume input files are in the parent directory or current working directory:

- `F31_PPMI_Vascular_SPSS.xlsx` (the merged analytic file, N = 413; produced by `spss/PPMI_VBS_Analysis.sps`)
- `DTI_Regions_of_Interest_14Jan2026.csv` (PPMI source, required for `random_forest_dti.py`)
- `FS7_ASEG_VOL_13May2026.csv` (PPMI source, required for `wmh_sensitivity.py`)

See `data/README.md` for instructions on obtaining the PPMI source files. Raw PPMI data is not redistributed in this repository.

## Usage

Each script is independent and can be run directly:

​```bash
python poisson_arr_per_factor.py
python random_forest_dti.py
python wmh_sensitivity.py
​```

The FDR correction utility can be used as a CLI tool to process any CSV of p-values:

​```bash
python fdr_correction.py input.csv output.csv --pcol p_raw --alpha 0.05
​```

It can also be imported into other Python code:

​```python
from fdr_correction import benjamini_hochberg
q_values, significant = benjamini_hochberg(p_values, alpha=0.05)
​```

## Outputs

Each analytic script writes a CSV results file to the working directory:

- `poisson_arr_per_factor.py` produces `table2_arr_per_factor.csv`
- `random_forest_dti.py` produces `rf_dti_sensitivity_results.csv`
- `wmh_sensitivity.py` produces `wmh_sensitivity_results.csv`

Canonical versions of these output files (matching the manuscript) are deposited in the `outputs/` directory at the repository root.

## Relationship to the SPSS pipeline

The SPSS pipeline (`spss/PPMI_VBS_Analysis.sps`) was the primary analytic workflow for the original F31 fellowship project. It produced the descriptive statistics, multivariable VBS regressions (Models A, B, C), trans-diagnostic ANOVA, and within-cohort Cohen's d effect sizes reported in the manuscript.

The Python pipeline in this folder was added for two reasons:

1. **Modified Poisson aRRs**: The manuscript reports adjusted Risk Ratios from modified Poisson regression (Zou 2004) rather than odds ratios. This method is preferred when the outcome prevalence is non-rare (DM/Pre-DM is 17.7% in this sample). SPSS does not implement modified Poisson natively, so this was re-run in Python using statsmodels with HC0 robust standard errors.

2. **Random forest and WMH sensitivity analyses**: These were exploratory analyses added after Dr. Chiu's review feedback and were implemented in Python from the start.

Both pipelines reference the same merged analytic file and produce results consistent with the published manuscript values, with minor differences in the random forest analysis attributable to seed-dependent cross-validation fold assignment (documented in the docstring of `random_forest_dti.py`).

## Reproducibility notes

- All scripts use deterministic random seeds (`RANDOM_SEED = 42`) where stochastic methods are involved (random forest CV fold assignment, permutation importance).
- The random forest analysis with n = 50 imaging-phenotyped participants is reported in the manuscript as exploratory; AUC point estimates are sensitive to seed choice, but qualitative findings (DTI dominating permutation importance, VBS components contributing measurable additional variance) are robust across seeds.
- Modified Poisson regression results match manuscript values to two decimal places.
- WMH analysis results match manuscript values exactly (β = −0.036, p = 0.357, −3.5% per VBS point).

## Citation

If you use or adapt this code, please cite the manuscript above and acknowledge PPMI per their Data Use Agreement (see `data/README.md` for the verbatim acknowledgment statement).

## Contact

For questions about the analysis, contact the corresponding author: Edward Ofori (Edward.Ofori@asu.edu).
