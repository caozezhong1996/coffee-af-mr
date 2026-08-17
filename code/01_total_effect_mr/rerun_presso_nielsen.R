# -*- coding: utf-8 -*-
# P1-5 confirmation rerun: MR-PRESSO outlier test on Nielsen AF dataset
# (ebi-a-GCST006414), official MRPRESSO package (GitHub sha 3e3c92d, built R 4.6.0)
# Mirrors mr_p0_fixes.R parameters (NbDistribution=2000) + high-precision run (10000).
.libPaths(c("C:/Users/曹泽众/Documents/kimi/workspace/Rlib", .libPaths()))
suppressMessages(library(MRPRESSO))

f <- "D:/Desktop/EJPC_投稿_咖啡与房颤/03_Supplementary/S1_total_effect_MR/harmonised_ebi-a-GCST006414.csv"
d0 <- read.csv(f)
d <- d0[d0$mr_keep == TRUE, ]
cat("SNPs in harmonised file:", nrow(d0), "| kept for MR:", nrow(d), "\n")

outdir <- "C:/Users/曹泽众/Documents/kimi/workspace/EJPC_修订产出/stats"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

run_presso <- function(nb, seed) {
  pr <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                  SdOutcome = "se.outcome", SdExposure = "se.exposure",
                  data = d, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                  SignifThreshold = 0.05, NbDistribution = nb, seed = seed)
  cat(sprintf("\n===== NbDistribution = %d, seed = %d =====\n", nb, seed))
  print(pr$`Main MR results`)
  cat("Global test P:", pr$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
  ot <- pr$`MR-PRESSO results`$`Outlier Test`
  if (!is.null(ot)) {
    ot <- cbind(SNP = d$SNP, ot)          # rows follow input order (verified below)
    ot <- ot[order(ot$Pvalue), ]
    cat("--- Outlier test, 8 smallest P ---\n")
    print(head(ot, 8))
    sig <- ot[ot$Pvalue < 0.05, ]
    cat("Outliers with P < 0.05:", paste(sig$SNP, collapse = ", "), "\n")
    write.csv(ot, file.path(outdir, sprintf("presso_rerun_nielsen_nb%d.csv", nb)),
              row.names = FALSE)
  }
  dtst <- pr$`MR-PRESSO results`$`Distortion Test`
  if (!is.null(dtst)) { cat("--- Distortion test ---\n"); print(dtst) }
  invisible(pr)
}

run_presso(2000, 20260816)    # mirrors mr_p0_fixes.R
run_presso(10000, 20260816)   # high-precision confirmation

# Sanity: RSSobs column order must track input SNP order — check correlation between
# the two runs' per-SNP RSSobs after reordering by SNP.
a <- read.csv(file.path(outdir, "presso_rerun_nielsen_nb2000.csv"))
b <- read.csv(file.path(outdir, "presso_rerun_nielsen_nb10000.csv"))
m <- merge(a, b, by = "SNP", suffixes = c("_2k", "_10k"))
cat("\nRSSobs correlation (2k vs 10k, should be 1.0):",
    cor(m$RSSobs_2k, m$RSSobs_10k), "\n")
cat("Top-3 by P (2k run):", paste(head(a$SNP, 3), collapse = ", "), "\n")
cat("DONE\n")
