Data Directory
Data availability
This repository does not contain raw PPMI participant data. Per the PPMI Data Use Agreement, raw data cannot be redistributed.
Authorized researchers can access the same data used in this analysis by registering with the Parkinson's Progression Markers Initiative (PPMI) and downloading from the LONI Image and Data Archive.

PPMI portal: https://www.ppmi-info.org
Data download: https://www.ppmi-info.org/access-data-specimens/download-data
RRID: SCR_006431
Data tier used in this analysis: Tier 1 (openly available to all registered PPMI users)

Data snapshot used for the published analysis

PPMI download date: February 12, 2026
Analytic sample: N = 413 (76 healthy controls, 145 prodromal, 192 clinically diagnosed PD)
PPMI 1.0 enrollment era: June 2010 to June 2019

Source files used
The analysis pipeline consumes the following files from a standard PPMI data download. File names reflect the dating convention used in our February 2026 download; reproductions from other download dates will have different timestamps in the filenames but identical content for participants enrolled by June 2019.
FilePurposeDemographics_01May2026.csvAge, sex, race/ethnicitySubject_Cohort_History_01May2026.csvCohort assignment (APPRDX)iu_genetic_consensus_20251025_20May2026.csvLRRK2/GBA/SNCA carrier statusMDS-UPDRS_Part_I_07Mar2026.csvNon-motor symptoms (clinician)MDS_UPDRS_Part_II__Patient_Questionnaire_07Mar2026.csvMotor ADLMDS-UPDRS_Part_III_07Mar2026.csvMotor exam + Hoehn & YahrMDS-UPDRS_Part_IV__Motor_Complications_07Mar2026.csvMotor complicationsMontreal_Cognitive_Assessment__MoCA__08Mar2026.csvCognitive screeningREM_Sleep_Behavior_Disorder_Screening_Questionnaire_20May2026.csvRBDSQUniversity_of_Pennsylvania_Smell_Identification_Test_UPSIT_20May2026.csvUPSIT smell testPD_Diagnosis_History_20May2026.csvPD diagnosis date for years_since_dxDTI_Regions_of_Interest_14Jan2026.csvNigral DTI eigenvaluesFS7_ASEG_VOL_13May2026.csvFreeSurfer 7 segmentation volumes (WMH)
Variable construction
The analytic file F31_PPMI_Vascular_SPSS.xlsx is the merged, derived dataset produced by:

The SPSS pipeline documented in spss/PPMI_VBS_Analysis.sps, which constructs VBS, VBS_micro, VBS_macro, and DM_Group from the PPMI Medical Conditions Log and Vital Signs files.
The Python pipeline in python/ for modified Poisson aRRs, random forest sensitivity analyses, and white matter hypointensity analysis.

Variable definitions:

DM_Group: 1 = diabetes mellitus or prediabetes; 0 = neither
VBS: total Vascular Burden Score (sum of 7 indicators, range 0 to 7)
VBS_micro: microvascular sub-score = HTN + HLD + OSA + OBESE_BMI (range 0 to 4)
VBS_macro: macrovascular sub-score = STROKE + CAD + AFIB (range 0 to 3)
OBESE_BMI: 1 if BMI >= 30 kg/m^2
Cohort_n: 0 = healthy control, 1 = prodromal, 2 = clinically diagnosed PD
years_since_dx: years between PD diagnosis date and baseline assessment (PD cohort only)

Reproducing the analysis

Register with PPMI and download the source files listed above.
Place files in this data/ directory.
Run spss/PPMI_VBS_Analysis.sps in SPSS Statistics 30+ to build F31_PPMI_Vascular_SPSS.xlsx.
Run scripts in python/ to reproduce aRRs and sensitivity analyses.

See the top-level README.md for full reproducibility instructions and dependency information.
Citation
If you use this code or reproduce this analysis, please cite:
Belnavis A, Chiu SY, Chen K, Thorpe R, Ofori E. Vascular Phenotyping in Parkinson's Disease: Diabetes Mellitus Operationalizes a Microvascular Metabolic Syndrome Cluster Across the PPMI Clinical Spectrum. Movement Disorders. 2026.
Acknowledgment
Data used in the preparation of this article were obtained on February 12, 2026 from the Parkinson's Progression Markers Initiative (PPMI) database (www.ppmi-info.org/access-data-specimens/download-data), RRID:SCR_006431. For up-to-date information on the study, visit www.ppmi-info.org.
