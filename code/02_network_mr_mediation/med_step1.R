jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
outdir <- "./mediation_mr"
dir.create(outdir, showWarnings=FALSE)
log <- file(file.path(outdir,"step1_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")

# ---- coffee exposure (existing 40-SNP instrument) ----
inst <- read.csv("mr_results/instruments.csv")
exp <- format_data(inst, type="exposure", snp_col="SNP", beta_col="beta.exposure", se_col="se.exposure",
                   effect_allele_col="effect_allele.exposure", other_allele_col="other_allele.exposure",
                   eaf_col="eaf.exposure", pval_col="pval.exposure")
exp$exposure <- "Coffee intake (cups/day, UKB)"
coffee_snps <- exp$SNP
cat("coffee instruments:", length(coffee_snps), "\n")

# ---- mediator datasets ----
meds <- data.frame(
  key  = c("sleep","insomnia","sbp","mvpa","crp","bmi"),
  id   = c("ukb-b-4424","ebi-a-GCST90018869","ieu-b-38","ebi-a-GCST006097","ebi-a-GCST90029070","ieu-b-40"),
  name = c("Sleep duration (h)","Insomnia (MVP, logOR)","Systolic BP (mmHg)",
           "MVPA (Klimentidis)","C-reactive protein","BMI (GIANT+UKB)"),
  stringsAsFactors=FALSE)

# ---- STEP 1: coffee -> mediator ----
step1_all <- list()
for (i in seq_len(nrow(meds))) {
  m <- meds[i, ]
  cat("\n===== STEP1:", m$key, m$id, "=====\n")
  r <- tryCatch({
    out <- extract_outcome_data(snps=coffee_snps, outcomes=m$id)
    cat("outcome rows:", nrow(out), "\n")
    out$outcome <- m$name
    dat <- harmonise_data(exp, out, action=2)
    d <- dat[dat$mr_keep==TRUE, ]
    cat("kept:", nrow(d), "\n")
    res <- mr(d, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
    res$mediator <- m$key; res$med_id <- m$id
    print(res[, c("method","nsnp","b","se","pval")])
    write.csv(res, file.path(outdir, paste0("step1_", m$key, ".csv")), row.names=FALSE)
    write.csv(d,  file.path(outdir, paste0("step1_harmonised_", m$key, ".csv")), row.names=FALSE)
    het <- mr_heterogeneity(d); write.csv(het, file.path(outdir, paste0("step1_het_", m$key, ".csv")), row.names=FALSE)
    ple <- mr_pleiotropy_test(d); write.csv(ple, file.path(outdir, paste0("step1_pleio_", m$key, ".csv")), row.names=FALSE)
    print(het[, c("method","Q","Q_df","Q_pval")]); print(ple[, c("egger_intercept","se","pval")])
    res
  }, error=function(e){ cat("ERROR:", conditionMessage(e), "\n"); NULL })
  step1_all[[m$key]] <- r
  Sys.sleep(2)
}

# ---- STEP 2a: mediator instruments ----
all_iv <- list()
for (i in seq_len(nrow(meds))) {
  m <- meds[i, ]
  cat("\n===== IV extract:", m$key, "=====\n")
  iv <- tryCatch(extract_instruments(m$id, p1=5e-8, clump=TRUE, r2=0.001, kb=10000),
                 error=function(e){ cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(iv)) {
    cat("instruments:", nrow(iv), "\n")
    iv$mediator <- m$key
    write.csv(iv, file.path(outdir, paste0("iv_", m$key, ".csv")), row.names=FALSE)
    all_iv[[m$key]] <- iv
  }
  Sys.sleep(2)
}
union_snps <- unique(c(coffee_snps, unlist(lapply(all_iv, function(x) x$SNP))))
cat("\nUNION SNPs:", length(union_snps), "\n")
writeLines(union_snps, file.path(outdir, "union_snps.txt"))

# Nielsen outcome for union set (API)
cat("\n===== Nielsen outcome extract (union) =====\n")
out_n <- tryCatch(extract_outcome_data(snps=union_snps, outcomes="ebi-a-GCST006414"),
                  error=function(e){ cat("ERROR:", conditionMessage(e), "\n"); NULL })
if (!is.null(out_n)) { cat("Nielsen rows:", nrow(out_n), "\n")
  write.csv(out_n, file.path(outdir, "outcome_nielsen_union.csv"), row.names=FALSE) }
cat("\nDONE\n")
sink(type="message"); sink(); close(log)
