jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages(library(ieugwasr))
outdir <- "./mediation_mr"
snps <- readLines(file.path(outdir, "union_snps.txt"))
chunks <- split(snps, ceiling(seq_along(snps)/50))
done_file <- file.path(outdir, "nielsen_chunks_progress.rds")
done <- if (file.exists(done_file)) readRDS(done_file) else list()
for (i in names(chunks)) {
  if (!is.null(done[[i]])) next
  cat("chunk", i, "/", length(chunks), "\n")
  a <- tryCatch(associations(chunks[[i]], "ebi-a-GCST006414", proxies=0),
                error=function(e){ cat("ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(a)) { done[[i]] <- a; saveRDS(done, done_file) }
  Sys.sleep(1)
}
alld <- do.call(rbind, done)
cat("rows:", nrow(alld), " unique SNPs:", length(unique(alld$rsid)), "\n")
write.csv(alld, file.path(outdir, "outcome_nielsen_union_raw.csv"), row.names=FALSE)
cat("DONE\n")
