# Coffee Consumption and Atrial Fibrillation — Analysis Code

Reproducible analysis code for the manuscript:

> **Coffee consumption and atrial fibrillation: total-effect triangulation and
> network Mendelian randomization suggest opposing blood-pressure and
> adiposity-linked pathways**
>
> Zezhong Cao, Rongrong Chen\*
> Department of Cardiology, Yangzhou University Affiliated Hospital
> (Yangzhou First People's Hospital), Yangzhou, Jiangsu, China
> \*Correspondence: chenrrxh@126.com

The study integrates five components:

1. **Total-effect Mendelian randomization (MR)** — univariable MR of genetically
   proxied coffee intake (cups/day) on atrial fibrillation (AF) in the two
   largest independent AF datasets (Nielsen et al. 2018; FinnGen R11), with
   MR-PRESSO outlier correction, radial MR, multivariable MR (MVMR) adjusting
   for BMI, smoking and alcohol, coffee-subtype instruments, reverse MR, and
   equivalence (TOST) testing.
2. **Network MR mediation** — two-step MR decomposing the total effect through
   six prespecified mediators (systolic blood pressure, BMI, moderate-to-
   vigorous physical activity, C-reactive protein, sleep duration, insomnia),
   with product-of-coefficients indirect effects, Sobel standard errors, MVMR
   direct effects, Steiger directionality filtering, Fisher combination across
   the two AF datasets, and Benjamini–Hochberg FDR control.
3. **Comparative beverage-behavior and cross-ancestry MR** — coffee intake
   alongside six other prespecified UK Biobank beverage behaviors (tea, sugar
   added to coffee/tea, sugar-sweetened and artificially sweetened drinks,
   artificial sweetener in tea) against three AF outcome datasets (AFGen 2018,
   frozen FinnGen endpoint, BioBank Japan), with Steiger filtering, power
   analysis, reverse MR, and a coffee + sugar-in-tea + BMI multivariable model.
4. **FAERS pharmacovigilance** — disproportionality analysis (PRR / ROR / IC)
   of caffeine-containing products and AF reports in the FDA Adverse Event
   Reporting System via the openFDA API, with deduplication, masking
   sensitivity analyses, age/sex-stratified analyses, and a ventricular
   arrhythmia/arrest positive control.
5. **Observational dose-response meta-analysis** — updated Greenland–Longnecker
   random-effects meta-analysis of prospective cohorts (per-cup trends,
   high-vs-low contrast, leave-one-out).

A **locus-level colocalization** module (coloc, ±500 kb windows at the FTO,
TMEM18, CYP1A1–CYP1A2, AHR, and PCMTD2 loci, with prior sensitivity) links
components 1–3 at the variant level.

All analyses use **publicly available data only** (OpenGWAS, FinnGen, openFDA).
No individual-level data were accessed.

## Repository layout

```
code/
  01_total_effect_mr/     Univariable MR, MVMR, radial MR, subtype & reverse MR
  02_network_mr_mediation/  Two-step network MR, mediator instruments, MVMR direct effects
  03_robustness/          Bidirectional MR, FTO-locus exclusion, weak-instrument-robust MVMR
  04_pharmacovigilance/   openFDA/FAERS disproportionality, stratified analyses, positive control, figures
  05_figures/             Manuscript figure scripts (300 dpi TIFF/PNG)
  06_phewas/              Phenome-wide association scan of the 40 instrument SNPs
  07_comparative_beverage/  Seven-beverage comparative MR, cross-ancestry outcomes,
                          Steiger/power/reverse/MVMR upgrades, BBJ replication helpers
  08_observational_meta/  Updated dose-response meta-analysis (Greenland–Longnecker)
  09_colocalization/      Locus-level coloc (FTO, TMEM18, CYP1A1–CYP1A2, AHR, PCMTD2)
data/
  data_sources.md         Every dataset used, with accession IDs and access links
revision_outputs/         Key machine-readable outputs (PheWAS, PRESSO rerun,
                          FAERS stratified + positive control, BBJ replication)
```

## Requirements

**R** (≥ 4.3; developed on R 4.6.0, Windows x86_64):

- CRAN: `TwoSampleMR`, `MRPRESSO`, `ieugwasr`, `quantreg`, `ggplot2`, `patchwork`
- GitHub: `remotes::install_github("WSpiller/MVMR")`,
  `remotes::install_github("WSpiller/RadialMR")`

**Python** (≥ 3.10): `pandas`, `numpy`, `matplotlib` (standard library otherwise).

**OpenGWAS API token** (free): register at <https://api.opengwas.io/profile/>,
then replace the placeholder `jwt <- "YOUR_OPENGWAS_TOKEN"` at the top of each
R script (or set the environment variable `OPENGWAS_JWT`). Tokens expire
~13 days after issue.

## Reproduction order

Run scripts from the repository root so that relative paths resolve.

**1. Total-effect MR (`code/01_total_effect_mr/`)**

| Order | Script | Purpose |
|---|---|---|
| 1 | `mr_main.R` | Primary IVW / weighted-median / MR-Egger / MR-PRESSO vs Nielsen 2018; writes `mr_results/instruments.csv` used by later steps |
| 2 | `mr_step1.R`, `mr_step1b.R`, `mr_p0_fixes.R` | First-pass / exploratory runs |
| 3 | `mr_wave2.R`, `mr_wave2b.R`, `mr_wave3.R`, `mr_wave3b.R` | Robustness extensions (leave-one-out, heterogeneity, plots) |
| 4 | `mr_r11.R` | Replication vs FinnGen R11 (place `finngen_R11_I9_AF.gz` in the working directory; download URL inside the script) |
| 5 | `mr_mvmr_exp.R`, `mr_mvmr_outcome.R`, `mr_mvmr.R` | MVMR vs Nielsen outcome (coffee + BMI + smoking initiation + alcohol frequency) |
| 6 | `mr_mvmr_r11.R` | MVMR vs FinnGen R11 |
| 7 | `condf.R` | Conditional F-statistics for MVMR instruments |
| 8 | `radial.R`, `radial2.R` | Radial MR (modified second-order weights), both outcomes |
| 9 | `wojcik.R`, `wojcik_r11.R` | Independent-instrument sensitivity using the multi-ancestry PAGE/Wojcik 2019 exposure, both outcomes |

Equivalence testing (TOST, bounds OR 0.80–1.25) and power analysis were
computed from the harmonised estimates produced above (Supplementary Table S8).

**2. Network MR mediation (`code/02_network_mr_mediation/`)**

| Order | Script | Purpose |
|---|---|---|
| 1 | `search_mediators.R` | Candidate-mediator screen across OpenGWAS |
| 2 | `med_step1.R` | Step 1: coffee → six mediators; mediator instrument extraction |
| 3 | `med_step1b.R`, `med_step1_insomnia_ukb.R` | Step-1 extensions (UK Biobank insomnia replication) |
| 4 | `med_clump.R`, `med_coffee_chunks.R`, `med_nielsen_chunks.R` | Clumping and chunked extraction of mediator/coffee SNPs in the AF outcome |
| 5 | `med_step2.R` | Step 2: mediator → AF (both datasets); indirect effects, Sobel tests |
| 6 | `med_mvmr.py` | Multivariable IVW of coffee + mediator → AF from local harmonised data (direct effects) |

**3. Robustness (`code/03_robustness/`)**

`revision_part1.R` … `revision_part4b.R`, `revision_analyses.R` — bidirectional
MR (BMI → coffee), FTO-locus-excluded sensitivity analyses, weak-instrument-
robust multivariable estimators (multivariable Egger, multivariable weighted
median), and Steiger directionality checks.

**4. Pharmacovigilance (`code/04_pharmacovigilance/`)**

| Script | Purpose |
|---|---|
| `faers_rebuild.py` | openFDA query of caffeine-containing products and AF reports; PRR / ROR / IC |
| `faers_sens_a.py` | Sensitivity analyses (brand-name phrasing, seriousness restriction, masking) |
| `faers_stratified.py`, `faers_stratified_extra.py` | Age- and sex-stratified disproportionality (post hoc); product-type breakdown |
| `faers_va_control.py` | Ventricular arrhythmia/arrest positive-control analysis |
| `faers_forest.R`, `faers_forest_v2.R` | Forest-plot figures |

Stratified and positive-control results: `revision_outputs/faers_stratified.json`,
`revision_outputs/faers_va_control.json`.

**5. Figures (`code/05_figures/`)**

`mr_figures_final.R`, `mr_figures_r11.R` — MR forest/scatter/funnel/leave-one-out
panels; `make_fig1.py`, `make_figs.py`, `make_integrated_figs.py`, `fig1_v2.py` —
integrated triangulation and mediation-network figures (300 dpi TIFF/PNG).

**6. Comparative beverage-behavior and cross-ancestry MR (`code/07_comparative_beverage/`)**

| Order | Script | Purpose |
|---|---|---|
| 1 | `run1_uvmr.R` | Core univariable MR: seven beverage behaviors → AF (AFGen 2018, `ebi-a-GCST006061`) |
| 2 | `run_batch.R` | Batch rerun for the other outcomes (`finn-b-I9_AF`, `bbj-a-71`); pass the outcome ID as the first argument |
| 3 | `run_reverse_mr.R` | Reverse MR: AF liability → each beverage behavior |
| 4 | `run_mvmr.R` | Multivariable MR: coffee + sugar added to tea + BMI → AF (AFGen 2018) |
| 5 | `run_upgrades.R` | Power / minimum detectable OR, Steiger filtering, radial IVW + RAPS, leave-one-out and single-SNP data |
| 6 | `run_presso.R`, `run_radial.R` | MR-PRESSO and radial-MR sensitivity for the coffee pairs |
| 7 | `mr_miyazawa_replication.py` + `bbj_replication/` | East Asian replication helpers (Miyazawa 2024 BBJ coffee GWAS harmonization; outputs in `revision_outputs/bbj_harmonized.csv`, `revision_outputs/miyazawa_replication.json`) |

Curated run outputs (harmonized data, per-pair estimates, heterogeneity,
pleiotropy, upgrades, European fixed-effect pool) are versioned under
`code/07_comparative_beverage/results/` — these are the numbers behind
manuscript Table 4 and Supplementary Table S16.

**7. Observational dose-response meta-analysis (`code/08_observational_meta/`)**

| Order | Script | Purpose |
|---|---|---|
| 1 | `run_meta_wing.py` | Greenland–Longnecker dose-response pooling from `extractions_meta_wing.csv`; writes `percup_trends.csv`, `high_vs_low.csv`, `loo_percup.csv`, `meta_wing_results.json` |
| 2 | `make_fig.py` | Forest-plot figure (`forest_meta_wing.png`) |

`pmids.txt` lists the included cohorts; `extractions_meta_wing.csv` is the
curated extraction sheet (doses, cases, person-years, adjusted HRs).

**8. Locus-level colocalization (`code/09_colocalization/`)**

| Order | Script | Purpose |
|---|---|---|
| 1 | `subset_windows.py` | Build ±500 kb outcome windows (FinnGen R11, Nielsen) around the lead coffee SNPs |
| 2 | `download_vcf.py` | Download the coffee exposure VCF (`ukb-b-5237.vcf.gz`, ~300 MB, from the IEU OpenGWAS OCI bucket; not stored in this repository) |
| 3 | `extract_coffee.py` | Extract exposure summary statistics for the same windows from the VCF |
| 4 | `run_coloc_v42.R` | coloc (single-causal-variant) at FTO, TMEM18, CYP1A1–CYP1A2, AHR for both outcomes; prior sensitivity (`coloc_v42_prior_sensitivity.csv`) |
| 5 | `run_coloc_pcmtd2.R` | Re-anchored PCMTD2-window run after the locus-assignment audit (rs6062682, GRCh38 chr20:64,260,467; verified against `finngen_R11_I9_AF.gz`) |

All curated per-window extracts needed by steps 4–5 are versioned in this
directory (`fg_pcmtd2_fix.tsv`, `coffee_pcmtd2_fix.csv`), so the coloc runs
reproduce without re-downloading the full summary statistics. The three larger
intermediate extracts (`fg_windows_500kb.tsv`, `nei_windows_500kb.tsv`,
`coffee_vcf_associations.csv`, ~2–5 MB each) are not versioned; they are
regenerated byte-for-byte by running steps 1–3 above. Results:
`coloc_v42_results.csv`, `coloc_v42_prior_sensitivity.csv`,
`coloc_pcmtd2_fix_results.csv`, `coloc_pcmtd2_fix_sensitivity.csv`
(Supplementary Table S20). This v4.2 pipeline supersedes earlier pilot coloc
scripts.

## Data availability

See [`data/data_sources.md`](data/data_sources.md) for every dataset,
accession ID, and access link. Summary statistics were downloaded from
OpenGWAS (<https://gwas.mrcieu.ac.uk/>) and FinnGen R11
(<https://r11.finngen.fi/>, free registration); FAERS data were queried
through the openFDA public endpoint (<https://api.fda.gov/drug/event.json>).

## Notes

- R scripts in this environment occasionally print a harmless segfault after
  the final `DONE` line (Windows/R 4.6.0 shutdown issue); all outputs are
  written before that point.
- `MRlap` (sample-overlap correction) is not installable in this environment
  (requires Rtools + `GenomicSEM` compilation). Overlap bias is addressed
  analytically in the manuscript; the primary exposure GWAS (`ukb-b-5237`) and
  the FinnGen R11 outcome share no samples.
- FAERS query strings are embedded in the Python scripts and in Supplementary
  Table S5.

## License and citation

Code is released under the MIT License (see `LICENSE`). If you use this code,
please cite the manuscript (see `CITATION.cff`).

## Post-hoc confirmation runs (2026-08-16)

- `code/01_total_effect_mr/rerun_presso_nielsen.R` — independent rerun of the
  MR-PRESSO per-SNP outlier test on the Nielsen AF dataset with the official
  MRPRESSO package (v1.0, GitHub sha 3e3c92d, R 4.6.0), mirroring
  `mr_p0_fixes.R` parameters (NbDistribution = 2000, seed 20260816).
  Output: `revision_outputs/presso_rerun_nielsen_nb2000.csv`.
  Result: confirms rs13387939 and rs1421085 as the two outliers (P < 0.05);
  rs2107308 marginal (P = 0.098); distortion test P = 0.10. The RSSobs column
  reproduces the original `p0_presso_outliers_ebi-a-GCST006414.csv` values
  exactly, validating the row-to-SNP mapping used in Supplementary Table S4b.
- `code/06_phewas/run_phewas.R` + `summarize_phewas.py` — phenome-wide
  association scan of the 40 instrument SNPs (OpenGWAS phewas endpoint,
  P < 5e-8; requires an OpenGWAS JWT in the environment variable
  `OPENGWAS_JWT`, never stored in this repository). Outputs in
  `revision_outputs/`: phewas_coffee_snps_raw.csv (3,904 associations across
  1,232 traits), phewas_category_summary.csv, phewas_snp_category_matrix.csv,
  phewas_snp_top2.csv. Summarized in Supplementary Table S15.
