jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
ao <- available_outcomes()
write.csv(ao, "mr_data/opengwas_ao.csv", row.names=FALSE)
cat("total datasets:", nrow(ao), "\n")
pick <- function(pat, n=12) {
  m <- ao[grepl(pat, ao$trait, ignore.case=TRUE), ]
  m <- m[order(-as.numeric(m$sample_size)), ]
  print(m[1:min(n, nrow(m)), c("id","trait","year","sample_size","nsnp","author")], row.names=FALSE)
}
cat("\n### sleep duration ###\n");  pick("^sleep duration")
cat("\n### insomnia/sleepless ###\n");  pick("insomnia|sleepless")
cat("\n### systolic BP ###\n");  pick("systolic blood pressure")
cat("\n### physical activity ###\n");  pick("physical activity")
cat("\n### accelerometer ###\n");  pick("acceler")
cat("\n### CRP ###\n");  pick("C-reactive|CRP")
cat("\n### IL-6 ###\n");  pick("interleukin-6|interleukin 6|IL-6")
