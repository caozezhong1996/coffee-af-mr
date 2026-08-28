#!/usr/bin/env Rscript
# args: outcome_id   -- loops all exposures, skips completed ones
.libPaths("C:/Users/曹泽众/Documents/kimi/workspace/Rlib")
suppressMessages({library(data.table); library(ieugwasr); library(TwoSampleMR)})
outdir <- "C:/Users/曹泽众/Documents/kimi/workspace/mr_results"
args <- commandArgs(trailingOnly = TRUE)
outcome_id <- args[1]

exposures <- list(
  c("ukb-b-8442",  "Sugar added to tea"),
  c("ukb-b-243",   "Sugar added to coffee"),
  c("ukb-b-2832",  "Fizzy drink intake SSB"),
  c("ukb-b-19703", "Low-calorie drink intake ASB"),
  c("ukb-b-5867",  "Artificial sweetener added to tea"),
  c("ukb-b-5237",  "Coffee intake"),
  c("ukb-b-6066",  "Tea intake")
)

run_one <- function(expo_id, expo_label, outcome_id) {
  otag <- gsub("[^A-Za-z0-9]+", "_", outcome_id)
  tag <- gsub("[^A-Za-z0-9]+", "_", paste(expo_label, otag))
  if (file.exists(file.path(outdir, paste0("mr_", tag, ".csv")))) { cat("skip (done):", expo_label, "\n"); return(invisible(NULL)) }
  cat("\n=== ", expo_label, "(", expo_id, ") ->", outcome_id, "\n")
  exp_dat <- tryCatch(extract_instruments(expo_id, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000),
                      error = function(e) { cat("instruments@5e-8 fail:", conditionMessage(e), "\n"); NULL })
  thr <- "5e-8"
  if (is.null(exp_dat) || nrow(exp_dat) < 3) {
    exp_dat <- tryCatch(extract_instruments(expo_id, p1 = 5e-6, clump = TRUE, r2 = 0.001, kb = 10000),
                        error = function(e) { cat("instruments@5e-6 fail:", conditionMessage(e), "\n"); NULL })
    thr <- "5e-6"
  }
  if (is.null(exp_dat) || nrow(exp_dat) < 3) { cat("SKIP: <3 instruments\n"); return(invisible(NULL)) }
  cat("instruments:", nrow(exp_dat), "(p<", thr, ")\n")
  exp_dat$exposure <- expo_label
  out_dat <- tryCatch(extract_outcome_data(exp_dat$SNP, outcome_id),
                      error = function(e) { cat("outcome fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(out_dat)) return(invisible(NULL))
  out_dat$outcome <- outcome_id
  dat <- harmonise_data(exp_dat, out_dat, action = 2)
  dat <- dat[dat$mr_keep, ]
  cat("SNPs after harmonisation:", nrow(dat), "\n")
  if (nrow(dat) < 3) { cat("SKIP: <3 after harmonisation\n"); return(invisible(NULL)) }
  dat$F_stat <- (dat$beta.exposure / dat$se.exposure)^2
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  het <- mr_heterogeneity(dat)
  pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  fwrite(as.data.table(dat), file.path(outdir, paste0("dat_", tag, ".csv")))
  fwrite(res, file.path(outdir, paste0("mr_", tag, ".csv")))
  fwrite(het, file.path(outdir, paste0("het_", tag, ".csv")))
  if (!is.null(pleio)) fwrite(pleio, file.path(outdir, paste0("pleio_", tag, ".csv")))
  ivw <- res[res$method == "Inverse variance weighted", ]
  cat(sprintf("IVW: OR=%.3f (%.3f-%.3f), P=%.3g, nsnp=%d\n",
              exp(ivw$b), exp(ivw$b - 1.96 * ivw$se), exp(ivw$b + 1.96 * ivw$se), ivw$pval, ivw$nsnp))
  invisible(NULL)
}

for (ex in exposures) {
  tryCatch(run_one(ex[1], ex[2], outcome_id), error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
  Sys.sleep(1)
}
cat("BATCH DONE\n")
