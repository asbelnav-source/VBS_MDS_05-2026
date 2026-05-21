* =====================================================================
* PPMI Vascular Burden Score Analysis - SPSS Syntax
* =====================================================================
* Project: Vascular Phenotyping in Parkinson's Disease: Diabetes Mellitus
*          Operationalizes a Microvascular Metabolic Syndrome Cluster
*          Across the PPMI Clinical Spectrum
* Author:  Alexander Belnavis (ASU College of Health Solutions)
* Co-authors: S.Y. Chiu, K. Chen, R. Thorpe, E. Ofori
* Data:    F31_PPMI_Vascular.sav (N = 413; PPMI 1.0 enrollment 2010-2019)
* PPMI download date: February 12, 2026
* Output reconstructed from PPMI_Vascular_Phenotyping_Analysis.xlsm
*   (SPSS Statistics v30; output exported 02-MAY-2026)
* =====================================================================
* NOTE: This file documents the SPSS-based primary analysis pipeline that
*       produced the initial logistic regression aORs and the multivariable
*       regression Models A, B, and C of composite VBS. In the final
*       manuscript, the per-factor analyses were re-run as modified Poisson
*       regression (yielding adjusted Risk Ratios with robust variance) in
*       Python (statsmodels) for clinical-epidemiology appropriateness with
*       a common (~17.7%) outcome. The modified Poisson and FDR correction
*       Python code is deposited separately in the same repository.
* =====================================================================



* ---------------------------------------------------------------------
* STEP 1: Data preparation (string-variable cleanup)
* ---------------------------------------------------------------------

DELETE VARIABLES Diabetes_Status_str Cohort_str Sex_str.


* ---------------------------------------------------------------------
* STEP 2: Variable dictionary inspection
* ---------------------------------------------------------------------

DISPLAY DICTIONARY.


* ---------------------------------------------------------------------
* STEP 3: Cohort and sample descriptive statistics
* ---------------------------------------------------------------------

FREQUENCIES VARIABLES=DM_Group Cohort_n Sex_M.

DESCRIPTIVES VARIABLES=AGE BMI VBS VBS_micro VBS_macro years_since_dx
  /STATISTICS=MEAN STDDEV MIN MAX.

DESCRIPTIVES VARIABLES=AGE BMI VBS VBS_micro VBS_macro
  /STATISTICS=MEAN STDDEV MIN MAX.


* ---------------------------------------------------------------------
* STEP 4: Risk-factor and demographic frequencies
* ---------------------------------------------------------------------

FREQUENCIES VARIABLES=Sex_M Cohort_n HTN HLD OSA OBESE_BMI STROKE CAD AFIB
  /ORDER=ANALYSIS.

DESCRIPTIVES VARIABLES=AGE BMI VBS VBS_micro VBS_macro
  /STATISTICS=MEAN STDDEV.

FREQUENCIES VARIABLES=Sex_M Cohort_n HTN HLD OSA OBESE_BMI STROKE CAD AFIB
  /ORDER=ANALYSIS.


* ---------------------------------------------------------------------
* STEP 5: Per-factor logistic regressions (initial aOR analyses)
* ---------------------------------------------------------------------

LOGISTIC REGRESSION VARIABLES HTN
  /METHOD=ENTER AGE Sex_M BMI
  /METHOD=ENTER DM_Group
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20) CUT(.5)
  /PRINT=CI(95) GOODFIT.

LOGISTIC REGRESSION VARIABLES HLD
  /METHOD=ENTER AGE Sex_M BMI
  /METHOD=ENTER DM_Group
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20) CUT(.5)
  /PRINT=CI(95) GOODFIT.

LOGISTIC REGRESSION VARIABLES OSA
  /METHOD=ENTER AGE Sex_M BMI
  /METHOD=ENTER DM_Group
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20) CUT(.5)
  /PRINT=CI(95) GOODFIT.

LOGISTIC REGRESSION VARIABLES OBESE_BMI
  /METHOD=ENTER AGE Sex_M
  /METHOD=ENTER DM_Group
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20) CUT(.5)
  /PRINT=CI(95) GOODFIT.

LOGISTIC REGRESSION VARIABLES STROKE
  /METHOD=ENTER AGE Sex_M BMI
  /METHOD=ENTER DM_Group
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20) CUT(.5)
  /PRINT=CI(95) GOODFIT.

LOGISTIC REGRESSION VARIABLES CAD
  /METHOD=ENTER AGE Sex_M BMI
  /METHOD=ENTER DM_Group
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20) CUT(.5)
  /PRINT=CI(95) GOODFIT.

LOGISTIC REGRESSION VARIABLES AFIB
  /METHOD=ENTER AGE Sex_M BMI
  /METHOD=ENTER DM_Group
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20) CUT(.5)
  /PRINT=CI(95) GOODFIT.


* ---------------------------------------------------------------------
* STEP 6: FDR correction table (Benjamini-Hochberg)
* ---------------------------------------------------------------------

DATASET ACTIVATE FDR_table.

LIST factor p_raw rank fdr_threshold q_value fdr_sig.

DATASET NAME FDR_table.

LIST factor p_raw rank fdr_threshold q_value fdr_sig.


* ---------------------------------------------------------------------
* STEP 7: Composite VBS multivariable regressions (Models A, B, C)
* ---------------------------------------------------------------------

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA CI(95)
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN
  /DEPENDENT VBS
  /METHOD=ENTER DM_Group AGE Sex_M BMI.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA CI(95)
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN
  /DEPENDENT VBS_micro
  /METHOD=ENTER DM_Group AGE Sex_M.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA CI(95)
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN
  /DEPENDENT VBS_macro
  /METHOD=ENTER DM_Group AGE Sex_M BMI.


* ---------------------------------------------------------------------
* STEP 8: Trans-diagnostic ANOVA (cohort x DM_Group)
* ---------------------------------------------------------------------

UNIANOVA VBS BY DM_Group Cohort_n WITH AGE Sex_M BMI
  /METHOD=SSTYPE(3)
  /INTERCEPT=INCLUDE
  /EMMEANS=TABLES(DM_Group*Cohort_n) WITH(AGE=MEAN Sex_M=MEAN BMI=MEAN)
  /EMMEANS=TABLES(DM_Group) WITH(AGE=MEAN Sex_M=MEAN BMI=MEAN) COMPARE ADJ(BONFERRONI)
  /PRINT=DESCRIPTIVE PARAMETER ETASQ
  /CRITERIA=ALPHA(.05)
  /DESIGN=DM_Group Cohort_n AGE Sex_M BMI DM_Group*Cohort_n.


* ---------------------------------------------------------------------
* STEP 9: Within-DM_Group VBS comparisons + Cohen's d
* ---------------------------------------------------------------------

T-TEST GROUPS=DM_Group(0 1)
  /MISSING=ANALYSIS
  /VARIABLES=VBS VBS_micro VBS_macro
  /CRITERIA=CI(.95)
  /ES DISPLAY(TRUE) STANDARDIZER(SD).


* ---------------------------------------------------------------------
* STEP 10: PD-only sensitivity analysis with disease duration
* ---------------------------------------------------------------------

DESCRIPTIVES VARIABLES=VBS DM_Group years_since_dx
  /STATISTICS=MEAN STDDEV.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA CI(95)
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN
  /DEPENDENT VBS
  /METHOD=ENTER DM_Group AGE Sex_M BMI years_since_dx.
