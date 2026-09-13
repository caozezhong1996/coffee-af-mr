#!/usr/bin/env Rscript
# Multivariable MR: coffee intake + sugar added to tea + BMI -> AF (AFGen EUR)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(data.table); library(ieugwasr); library(TwoSampleMR)})
outdir <- "C:/Users/曹泽众/Documents/kimi/workspace/mr_results"

expo_ids <- c("ukb-b-5237", "ukb-b-8442", "ieu-b-40")  # coffee, sugar-in-tea, BMI
labels <- c(`ukb-b-5237` = "Coffee intake", `ukb-b-8442` = "Sugar added to tea", `ieu-b-40` = "BMI")

cat("Extracting exposure instruments (joint) ...\n")
exposure_dat <- mv_extract_exposures(expo_ids, pval_threshold = 5e-8, clump_r2 = 0.001, clump_kb = 10000)
cat("joint instruments:", nrow(exposure_dat), "\n")
print(table(exposure_dat$exposure))

cat("Extracting outcome data ...\n")
outcome_dat <- extract_outcome_data(exposure_dat$SNP, "ebi-a-GCST006061")
outcome_dat$outcome <- "Atrial fibrillation (AFGen 2018)"

mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
res <- mv_multiple(mvdat)
res$result$exposure_label <- labels[res$result$exposure]
print(res$result[, c("exposure_label", "nsnp", "b", "se", "pval")])
fwrite(res$result, file.path(outdir, "mvmr_coffee_sugartea_bmi_AFGen.csv"))
cat("MVMR DONE\n")
