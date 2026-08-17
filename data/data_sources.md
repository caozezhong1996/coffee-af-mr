# Data sources

All datasets are publicly available. GWAS summary statistics were accessed
through OpenGWAS (<https://gwas.mrcieu.ac.uk/>) unless noted otherwise.
Access dates: 2026-08-05 to 2026-08-08.

## Exposures

| Role | Phenotype | Source / cohort | OpenGWAS ID | Sample size |
|---|---|---|---|---|
| Primary exposure | Coffee intake, cups/day | UK Biobank | `ukb-b-5237` | 428,860 |
| Sensitivity exposure | Coffee intake (multi-ancestry) | PAGE / Wojcik et al. 2019 | `ebi-a-GCST008028` | 15,837 |
| Subtype exposure | Decaffeinated coffee | Pirastu et al. 2022 | `ebi-a-GCST90096908` | 3–11 SNP instruments |
| Subtype exposure | Ground coffee | Pirastu et al. 2022 | `ebi-a-GCST90096913` | 3–11 SNP instruments |
| Subtype exposure | Instant coffee | Pirastu et al. 2022 | `ebi-a-GCST90096914` | 3–11 SNP instruments |

## Outcomes

| Role | Phenotype | Source | ID / access | Cases / controls |
|---|---|---|---|---|
| Primary outcome | Atrial fibrillation | Nielsen et al. 2018 | `ebi-a-GCST006414` | 60,620 / 970,216 |
| Replication outcome | Atrial fibrillation (`I9_AF`) | FinnGen release R11 | <https://r11.finngen.fi/> (free registration; file `finngen_R11_I9_AF.gz`) | 55,853 / 231,952 |

## Mediators (network MR)

| Mediator | Source | OpenGWAS ID |
|---|---|---|
| Systolic blood pressure | UK Biobank (IEU) | `ieu-b-38` |
| Body mass index | GIANT + UK Biobank | `ieu-b-40` |
| Moderate-to-vigorous physical activity | Klimentidis et al. | `ebi-a-GCST006097` |
| C-reactive protein | — | `ebi-a-GCST90029070` |
| Sleep duration | UK Biobank | `ukb-b-4424` |
| Insomnia (primary) | Million Veteran Program | `ebi-a-GCST90018869` |
| Insomnia (replication) | UK Biobank | `ukb-b-3957` |

## MVMR adjustment covariates (total-effect component)

| Covariate | Source | OpenGWAS ID |
|---|---|---|
| Body mass index | UK Biobank | `ukb-b-19953` |
| Smoking initiation | GSCAN | `ieu-b-4877` |
| Alcohol intake frequency | UK Biobank | `ukb-b-5779` |

## Pharmacovigilance

| Source | Access | Scope |
|---|---|---|
| FDA Adverse Event Reporting System (FAERS) via openFDA | <https://api.fda.gov/drug/event.json> (no key required for low-volume queries) | 11,882,970 serious reports screened; query strings embedded in `code/04_pharmacovigilance/faers_rebuild.py` and Supplementary Table S5 |
