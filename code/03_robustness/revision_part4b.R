# Part 4b: fetch remaining Nielsen SNPs in small sub-chunks
jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
ws <- "."
med <- file.path(ws, "mediation_mr")
mv_exp <- readRDS(file.path(ws, "mr_results", "mv_exp.rds"))
snps <- unique(mv_exp$SNP)
chunks <- split(snps, ceiling(seq_along(snps)/150))
need <- chunks[[3]]
sub <- split(need, ceiling(seq_along(need)/40))
prog_file <- file.path(med, "nielsen_mvmr_sub3.rds")
done <- if (file.exists(prog_file)) readRDS(prog_file) else list()
for (i in names(sub)) {
  if (!is.null(done[[i]])) next
  cat("sub-chunk", i, "/", length(sub), " (", length(sub[[i]]), " SNPs )\n")
  a <- tryCatch(extract_outcome_data(snps=sub[[i]], outcomes="ebi-a-GCST006414"),
                error=function(e){ cat("ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(a)) { done[[i]] <- a; saveRDS(done, prog_file) }
  Sys.sleep(2)
}
cat("done sub-chunks:", length(done), "/", length(sub), "\n")
if (length(done) == length(sub)) {
  main <- readRDS(file.path(med, "nielsen_mvmr_progress.rds"))
  main[["3"]] <- do.call(rbind, done)
  saveRDS(main, file.path(med, "nielsen_mvmr_progress.rds"))
  cat("MERGED chunk 3 into main progress\n")
}
