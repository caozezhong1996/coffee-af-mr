#!/usr/bin/env Rscript
# Quality-upgrade batch (local data only):
#   1) power / minimum detectable OR at 80% power
#   2) Steiger directionality test + steiger filtering
#   3) radial IVW + RAPS sensitivity estimators
#   4) leave-one-out + single-SNP data (for Figure S1)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(data.table); library(TwoSampleMR)})
outdir <- "C:/Users/曹泽众/Documents/kimi/workspace/mr_results"

exposures <- list(
  c("ukb-b-8442",  "Sugar added to tea"),
  c("ukb-b-243",   "Sugar added to coffee"),
  c("ukb-b-2832",  "Fizzy drink intake SSB"),
  c("ukb-b-19703", "Low-calorie drink intake ASB"),
  c("ukb-b-5867",  "Artificial sweetener added to tea"),
  c("ukb-b-5237",  "Coffee intake"),
  c("ukb-b-6066",  "Tea intake")
)

z_a <- qnorm(0.975); z_b <- qnorm(0.80)

power_rows <- list(); steiger_rows <- list(); sens_rows <- list()

for (e in exposures) {
  expo_id <- e[1]; expo_label <- e[2]
  tag <- gsub("[^A-Za-z0-9]+", "_", expo_label)
  f <- file.path(outdir, paste0("dat_", tag, "_ebi_a_GCST006061.csv"))
  if (!file.exists(f)) f <- file.path(outdir, paste0("dat_", tag, ".csv"))
  dat <- as.data.frame(fread(f))
  cat("\n=== ", expo_label, " (nSNP:", nrow(dat), ")\n")

  ## ---- 1) power from realized SE (AFGen) ----
  ivw <- mr(dat, method_list = "mr_ivw")
  se1 <- ivw$se[1]
  # EUR meta SE (fixed effect with FinnGen)
  ff <- file.path(outdir, paste0("dat_", tag, "_finn_b_I9_AF.csv"))
  se_meta <- NA
  if (file.exists(ff)) {
    dat2 <- as.data.frame(fread(ff))
    ivw2 <- mr(dat2, method_list = "mr_ivw")
    se_meta <- sqrt(1 / (1/se1^2 + 1/ivw2$se[1]^2))
  }
  mde_afgen <- exp((z_a + z_b) * se1)
  pwr_hr11_afgen <- pnorm(log(1.10)/se1 - z_a)
  mde_meta <- if (is.na(se_meta)) NA else exp((z_a + z_b) * se_meta)
  pwr_hr11_meta <- if (is.na(se_meta)) NA else pnorm(log(1.10)/se_meta - z_a)
  power_rows[[length(power_rows)+1]] <- data.frame(
    exposure = expo_label, nsnp = nrow(dat), se_ivw_AFGen = se1,
    MDE80_OR_AFGen = mde_afgen, power_to_detect_OR1.10_AFGen = pwr_hr11_afgen,
    se_ivw_EURmeta = se_meta, MDE80_OR_EURmeta = mde_meta,
    power_to_detect_OR1.10_EURmeta = pwr_hr11_meta)
  cat(sprintf("  AFGen: MDE80 OR=%.2f, power(OR1.10)=%.2f | meta: MDE80 OR=%s, power=%s\n",
      mde_afgen, pwr_hr11_afgen,
      ifelse(is.na(mde_meta),"—",sprintf("%.2f",mde_meta)),
      ifelse(is.na(pwr_hr11_meta),"—",sprintf("%.2f",pwr_hr11_meta))))

  ## ---- 2) Steiger directionality ----
  dt <- tryCatch(directionality_test(dat), error=function(e){cat("  directionality err:",conditionMessage(e),"\n"); NULL})
  sf <- tryCatch(steiger_filtering(dat), error=function(e){cat("  steiger err:",conditionMessage(e),"\n"); NULL})
  n_viol <- NA
  if (!is.null(sf) && "steiger_dir" %in% names(sf)) n_viol <- sum(!sf$steiger_dir, na.rm=TRUE)
  if (!is.null(dt)) {
    steiger_rows[[length(steiger_rows)+1]] <- data.frame(
      exposure = expo_label, nsnp = nrow(dat),
      snp_r2_exposure = dt$snp_r2.exposure, snp_r2_outcome = dt$snp_r2.outcome,
      correct_causal_direction = dt$correct_causal_direction,
      steiger_pval = dt$steiger_pval, steiger_violations = n_viol)
    cat(sprintf("  Steiger: correct direction = %s (P = %.3g), violations = %s\n",
        dt$correct_causal_direction, dt$steiger_pval, ifelse(is.na(n_viol),"—",n_viol)))
  }

  ## ---- 3) radial IVW + RAPS ----
  sens <- tryCatch(mr(dat, method_list = c("mr_ivw_radial", "mr_raps")),
                   error=function(e){cat("  radial/raps err:",conditionMessage(e),"\n"); NULL})
  if (!is.null(sens)) {
    sens$exposure_label <- expo_label
    sens_rows[[length(sens_rows)+1]] <- as.data.frame(sens)
    print(as.data.frame(sens)[, c("method","nsnp","b","se","pval")])
  }

  ## ---- 4) leave-one-out + singlesnp ----
  loo <- tryCatch(mr_leaveoneout(dat), error=function(e) NULL)
  if (!is.null(loo)) fwrite(as.data.table(loo), file.path(outdir, paste0("loo_", tag, "_AFGen.csv")))
  sgl <- tryCatch(mr_singlesnp(dat), error=function(e) NULL)
  if (!is.null(sgl)) fwrite(as.data.table(sgl), file.path(outdir, paste0("singlesnp_", tag, "_AFGen.csv")))
}

fwrite(rbindlist(power_rows), file.path(outdir, "upgrade_power.csv"))
if (length(steiger_rows)) fwrite(rbindlist(steiger_rows), file.path(outdir, "upgrade_steiger.csv"))
if (length(sens_rows)) fwrite(rbindlist(sens_rows), file.path(outdir, "upgrade_radial_raps.csv"))
cat("\nUPGRADES DONE\n")
