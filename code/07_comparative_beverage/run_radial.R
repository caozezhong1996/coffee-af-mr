.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(data.table); library(TwoSampleMR)})
outdir <- "C:/Users/曹泽众/Documents/kimi/workspace/mr_results"
exposures <- list(
  c("Sugar added to tea"), c("Sugar added to coffee"), c("Fizzy drink intake SSB"),
  c("Low-calorie drink intake ASB"), c("Artificial sweetener added to tea"),
  c("Coffee intake"), c("Tea intake"))
rows <- list()
for (e in exposures) {
  expo_label <- e[1]; tag <- gsub("[^A-Za-z0-9]+", "_", expo_label)
  f <- file.path(outdir, paste0("dat_", tag, "_ebi_a_GCST006061.csv"))
  if (!file.exists(f)) f <- file.path(outdir, paste0("dat_", tag, ".csv"))
  dat <- as.data.frame(fread(f))
  r <- tryCatch(mr(dat, method_list = "mr_ivw_radial"), error=function(e){cat(expo_label,"ERR:",conditionMessage(e),"\n"); NULL})
  if (!is.null(r)) {
    df <- as.data.frame(r); df$exposure_label <- expo_label
    rows[[length(rows)+1]] <- df
    cat(sprintf("%-32s radial IVW: b=%.4f se=%.4f P=%.4f (nsnp=%d)\n", expo_label, df$b[1], df$se[1], df$pval[1], df$nsnp[1]))
  }
}
if (length(rows)) fwrite(rbindlist(rows, fill=TRUE), file.path(outdir, "upgrade_radial_ivw.csv"))
cat("RADIAL DONE\n")
