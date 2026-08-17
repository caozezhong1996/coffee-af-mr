# P1-7 PheWAS: scan 40 coffee instrument SNPs across OpenGWAS
# Token read from environment variable OPENGWAS_JWT (set by caller; never written to disk).
.libPaths(c("C:/Users/曹泽众/Documents/kimi/workspace/Rlib", .libPaths()))
suppressMessages({library(ieugwasr); library(data.table)})

if (nchar(Sys.getenv("OPENGWAS_JWT")) < 100) stop("OPENGWAS_JWT not set")

f <- "D:/Desktop/EJPC_投稿_咖啡与房颤/03_Supplementary/S1_total_effect_MR/harmonised_ebi-a-GCST006414.csv"
snps <- fread(f)$SNP
cat("Instrument SNPs:", length(snps), "\n")

# connectivity probe with 2 SNPs
probe <- tryCatch(ieugwasr::phewas(variants = snps[1:2], pval = 5e-8),
                  error = function(e) { cat("PROBE ERROR:", conditionMessage(e), "\n"); NULL })
if (is.null(probe)) stop("API probe failed")
cat("probe OK, rows:", nrow(probe), "\n")

# full scan in batches of 10 to stay within rate limits
out <- list()
for (i in seq(1, length(snps), by = 10)) {
  v <- snps[i:min(i + 9, length(snps))]
  cat("batch", ceiling(i/10), ":", paste(v, collapse = ","), "\n")
  r <- tryCatch(ieugwasr::phewas(variants = v, pval = 5e-8),
                error = function(e) { cat("  error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r) && nrow(r)) out[[length(out) + 1]] <- r
  Sys.sleep(2)
}
res <- rbindlist(out, fill = TRUE)
cat("TOTAL association rows:", nrow(res), "\n")

od <- "C:/Users/曹泽众/Documents/kimi/workspace/EJPC_修订产出/stats"
fwrite(res, file.path(od, "phewas_coffee_snps_raw.csv"))
cat("saved raw\n")

# quick summary: top traits by SNP count
if (nrow(res)) {
  tab <- res[, .N, by = .(trait)][order(-N)]
  cat("\n--- Top 25 traits by number of associated instrument SNPs ---\n")
  print(head(tab, 25))
  cat("\n--- Per-SNP hit counts (top 10) ---\n")
  print(head(res[, .N, by = snp][order(-N)], 10))
}
cat("DONE\n")
