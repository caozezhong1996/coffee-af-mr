jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages(library(ieugwasr))
outdir <- "./mediation_mr"
inst <- read.csv("mr_results/instruments.csv")
coffee <- data.frame(SNP=inst$SNP, pval=inst$pval.exposure, src="coffee", stringsAsFactors=FALSE)
for (k in c("sleep","insomnia","sbp","mvpa","crp","bmi")) {
  iv <- read.csv(file.path(outdir, paste0("iv_", k, ".csv")))
  med <- data.frame(SNP=iv$SNP, pval=iv$pval.exposure, src=k, stringsAsFactors=FALSE)
  both <- rbind(coffee, med)
  both <- both[!duplicated(both$SNP), ]
  cl <- ld_clump(dplyr::tibble(rsid=both$SNP, pval=both$pval, id=both$src), clump_r2=0.001, clump_kb=10000)
  keep <- both[both$SNP %in% cl$rsid, ]
  write.csv(keep, file.path(outdir, paste0("mvmr_snps_", k, ".csv")), row.names=FALSE)
  cat(k, "joint clumped:", nrow(keep), "\n")
  Sys.sleep(1)
}
cat("DONE\n")
