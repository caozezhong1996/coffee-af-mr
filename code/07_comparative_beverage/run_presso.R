.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(data.table); library(MRPRESSO)})
outdir <- "C:/Users/曹泽众/Documents/kimi/workspace/mr_results"
for (tag in c("Coffee_intake_ebi_a_GCST006061", "Sugar_added_to_coffee_ebi_a_GCST006061")) {
  cat("\n### ", tag, "\n")
  dat <- tryCatch(as.data.frame(fread(file.path(outdir, paste0("dat_", tag, ".csv")))),
                  error = function(e) { cat("read fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(dat)) next
  pr <- tryCatch(
    mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
              SdOutcome = "se.outcome", SdExposure = "se.exposure",
              OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat,
              NbDistribution = 3000, SignifThreshold = 0.05),
    error = function(e) { cat("presso fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(pr)) next
  print(pr$`Main MR results`)
  cat("Global test P:", pr$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
  ot <- pr$`MR-PRESSO results`$`Outlier Test`
  if (!is.null(ot)) {
    outl <- ot[!is.na(ot$Pvalue) & ot$Pvalue < 0.05, ]
    cat("outliers (P<0.05):", nrow(outl), "\n")
    if (nrow(outl)) print(outl)
  }
  dtst <- pr$`MR-PRESSO results`$`Distortion Test`
  if (!is.null(dtst)) cat("Distortion P:", dtst$Pvalue, "\n")
  capture.output(print(pr), file = file.path(outdir, paste0("presso_", tag, ".txt")))
  cat("saved\n")
}
cat("ALL DONE\n")
