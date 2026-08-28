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

## Comparative beverage-behavior and cross-ancestry module (v1.2)

| Role | Phenotype | Source | ID / access |
|---|---|---|---|
| Corroborating outcome | Atrial fibrillation | AFGen 2018 (Roselli et al.) | `ebi-a-GCST006061` (65,446 cases) |
| Comparative outcome | Atrial fibrillation (`I9_AF`, frozen endpoint) | FinnGen (OpenGWAS freeze) | `finn-b-I9_AF` |
| East Asian outcome | Atrial fibrillation | BioBank Japan | `bbj-a-71` |
| Comparative exposure | Tea intake | UK Biobank | `ukb-b-6066` |
| Comparative exposure | Sugar added to coffee | UK Biobank | `ukb-b-243` |
| Comparative exposure | Sugar added to tea | UK Biobank | `ukb-b-8442` |
| Comparative exposure | Fizzy-drink intake (sugar-sweetened) | UK Biobank | `ukb-b-2832` |
| Comparative exposure | Low-calorie-drink intake (artificially sweetened) | UK Biobank | `ukb-b-19703` |
| Comparative exposure | Artificial sweetener added to tea | UK Biobank | `ukb-b-5867` |
| MVMR covariate | Body mass index | GIANT + UK Biobank | `ieu-b-40` |
| East Asian replication exposure | Coffee consumption | Miyazawa et al. 2024 (BBJ) | harmonized in `revision_outputs/bbj_harmonized.csv` |

## Observational dose-response meta-analysis (v1.2)

| Source | Access | Scope |
|---|---|---|
| PubMed / Europe PMC | search and full-text XML via NCBI E-utilities and Europe PMC (`pmids.txt` in `code/08_observational_meta/`) | Updated Greenland–Longnecker dose-response meta-analysis of prospective cohorts; per-cup trends, high-vs-low contrast, leave-one-out |

## Locus-level colocalization (v1.2)

| Role | Phenotype | Source | ID / access |
|---|---|---|---|
| Exposure | Coffee intake | UK Biobank (IEU OpenGWAS, VCF via OCI) | `ukb-b-5237` (`ukb-b-5237.vcf.gz`, ~300 MB, downloaded by `code/09_colocalization/download_vcf.py`; not stored in this repository) |
| Outcome | Atrial fibrillation | FinnGen R11 summary statistics | `finngen_R11_I9_AF.gz` (<https://r11.finngen.fi/>) |
| Outcome | Atrial fibrillation | Nielsen et al. 2018 summary statistics | `ebi-a-GCST006414` (OpenGWAS `/associations` endpoint, chunked by `fetch`-style queries) |

Windows analysed: ±500 kb around the lead coffee SNP at five loci (FTO, TMEM18, CYP1A1–CYP1A2, AHR, PCMTD2). The PCMTD2 window was re-anchored after a locus-assignment audit (rs6062682, GRCh38 chr20:64,260,467; verified against `finngen_R11_I9_AF.gz`) — see `run_coloc_pcmtd2.R`. Curated per-window extracts used by the coloc runs are versioned under `code/09_colocalization/`.
