# T2D Complications GWAS (CLSA) - Working Phenotype Summary

This repo is set up from the instructions in:

- `0.admin/1.email.pdf`
- `0.admin/T2D-complications_meta-analysis_AnalysisPlan_Phenodefinitions_22122025.docx`

## Requested study concept (from admin docs)

- Cases and controls are **all T2D participants**.
- Controls: T2D without complications and **T2D duration >= 2 years**.
- Main analyses: no strict upper/lower time limit between T2D diagnosis and complication diagnosis (maximize sample size).
- Sensitivity analyses: apply a **2-year timing requirement** between T2D and complication diagnosis for cases (implemented here as complication age >= T2D age + 2 years).
- Requested individual phenotypes in the analysis plan:
  - Ophthalmic
  - Neurological
  - Renal
  - Cardiovascular (includes PAD/PVD in the plan text)
  - Cerebrovascular
- Combined phenotypes:
  - Microvascular
  - Macrovascular
  - Micro+Macro

## Current CLSA implementation status

Implemented (self-report-based, T2D-restricted):

- `OPHTH` (diabetic retinopathy only; FUP1/FUP2)
- `RENAL` (kidney disease/failure and/or dialysis)
- `CARDIO` (heart disease / angina / MI / CAB / PVD / hypertension)
- `CEREBRO` (stroke/CVA / TIA / effects of stroke/TIA)
- `MICRO` (ophthalmic ∪ renal; neurological unavailable in CLSA mapping)
- `MACRO` (cardiovascular ∪ cerebrovascular)
- `MICROMACRO` (microvascular ∪ macrovascular)
- `MAIN` and `SENS2Y`

Not currently available from reviewed CLSA variables:

- `NEURO` (diabetic neuropathy / neurological diabetic complication specific variable not identified)

## Important deviations to report when sharing summary stats

- Microvascular phenotype excludes neurological complications (currently unavailable in CLSA variable mapping).
- Ophthalmic phenotype is conservative and likely under-ascertained (explicit diabetic retinopathy only, not all diabetic eye complications).
- Renal phenotype is broad and may include non-diabetic kidney disease among T2D participants.
- Some timing-sensitive case definitions are unavailable due missing complication/T2D age variables in parts of CLSA (especially FUP1 diabetes age).
- Main analyses drop cases where known complication age precedes known T2D age (to avoid implausible reverse timing), which is slightly stricter than a literal “no time limit” rule.
