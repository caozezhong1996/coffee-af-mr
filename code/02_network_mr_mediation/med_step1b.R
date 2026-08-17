jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
outdir <- "./mediation_mr"
log <- file(file.path(outdir,"step1b_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")

# insomnia IV from UKB fallback
iv <- tryCatch(extract_instruments("ukb-b-3957", p1=5e-8, clump=TRUE, r2=0.001, kb=10000),
               error=function(e){ cat("ERROR:", conditionMessage(e), "\n"); NULL })
if (!is.null(iv)) {
  cat("insomnia (ukb-b-3957) instruments:", nrow(iv), "\n")
  iv$mediator <- "insomnia"
  write.csv(iv, file.path(outdir, "iv_insomnia.csv"), row.names=FALSE)
}

# rebuild union
coffee <- read.csv("mr_results/instruments.csv")$SNP
ivfiles <- list.files(outdir, pattern="^iv_.*\\.csv$", full.names=TRUE)
all_iv <- lapply(ivfiles, read.csv)
union_snps <- unique(c(coffee, unlist(lapply(all_iv, function(x) x$SNP))))
cat("UNION SNPs:", length(union_snps), "\n")
writeLines(union_snps, file.path(outdir, "union_snps.txt"))

# Nielsen outcome for union
cat("\n===== Nielsen outcome extract =====\n")
out_n <- extract_outcome_data(snps=union_snps, outcomes="ebi-a-GCST006414")
cat("Nielsen rows:", nrow(out_n), "\n")
write.csv(out_n, file.path(outdir, "outcome_nielsen_union.csv"), row.names=FALSE)
cat("\nDONE\n")
sink(type="message"); sink(); close(log)
