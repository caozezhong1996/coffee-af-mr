jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr); library(MVMR)})
outdir <- "./mediation_mr"
log <- file(file.path(outdir,"step2_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")

meds <- data.frame(
  key  = c("sleep","insomnia","sbp","mvpa","crp","bmi"),
  name = c("Sleep duration","Insomnia","Systolic BP","MVPA","C-reactive protein","BMI"),
  stringsAsFactors=FALSE)

# ---------- load outcome data (local) ----------
niel <- read.csv(file.path(outdir,"outcome_nielsen_union_raw.csv"), check.names=FALSE)
cat("nielsen cols:", paste(names(niel), collapse=","), "\n")
pcol <- intersect(c("p","pval","p_value"), names(niel))[1]
mk_out_nielsen <- function(snps, name){
  x <- niel[niel$rsid %in% snps, ]
  o <- format_data(data.frame(SNP=x$rsid, beta=x$beta, se=x$se, eaf=x$eaf,
                              effect_allele=toupper(x$ea), other_allele=toupper(x$nea),
                              pval=x[[pcol]]),
                   type="outcome", snp_col="SNP", beta_col="beta", se_col="se",
                   effect_allele_col="effect_allele", other_allele_col="other_allele",
                   eaf_col="eaf", pval_col="pval")
  o$outcome <- name; o
}
r11 <- read.delim(file.path(outdir,"r11_union_hits.tsv"), comment.char="", check.names=FALSE)
names(r11)[1] <- "chrom"
mk_out_r11 <- function(snps, name){
  pick <- sapply(strsplit(r11$rsids, ","), function(v){ m <- v[v %in% snps]; if(length(m)) m[1] else NA })
  x <- r11[!is.na(pick), ]; x$SNP <- pick[!is.na(pick)]
  o <- format_data(data.frame(SNP=x$SNP, beta=x$beta, se=x$sebeta, eaf=x$af_alt,
                              effect_allele=toupper(x$alt), other_allele=toupper(x$ref), pval=x$pval),
                   type="outcome", snp_col="SNP", beta_col="beta", se_col="se",
                   effect_allele_col="effect_allele", other_allele_col="other_allele",
                   eaf_col="eaf", pval_col="pval")
  o$outcome <- name; o
}

# ---------- STEP 2: mediator -> AF ----------
for (i in seq_len(nrow(meds))) {
  k <- meds$key[i]
  iv <- read.csv(file.path(outdir, paste0("iv_", k, ".csv")))
  if (!"exposure" %in% names(iv)) iv$exposure <- meds$name[i]
  cat("\n########## STEP2:", k, " nsnp:", nrow(iv), "##########\n")
  for (oc in c("nielsen","r11")) {
    oname <- ifelse(oc=="nielsen", "AF (Nielsen 2018)", "AF/flutter (FinnGen R11)")
    out <- if (oc=="nielsen") mk_out_nielsen(iv$SNP, oname) else mk_out_r11(iv$SNP, oname)
    dat <- harmonise_data(iv, out, action=2)
    d <- dat[dat$mr_keep==TRUE, ]
    cat("\n--- ", k, "->", oname, " kept:", nrow(d), "---\n")
    res <- mr(d, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
    print(res[, c("method","nsnp","b","se","pval")])
    write.csv(res, file.path(outdir, paste0("step2_", k, "_", oc, ".csv")), row.names=FALSE)
    write.csv(d,  file.path(outdir, paste0("step2_harmonised_", k, "_", oc, ".csv")), row.names=FALSE)
    het <- mr_heterogeneity(d); write.csv(het, file.path(outdir, paste0("step2_het_", k, "_", oc, ".csv")), row.names=FALSE)
    ple <- mr_pleiotropy_test(d); write.csv(ple, file.path(outdir, paste0("step2_pleio_", k, "_", oc, ".csv")), row.names=FALSE)
    cat("Egger intercept:", ple$egger_intercept, " p:", ple$pval, "\n")
  }
}

# ---------- INDIRECT EFFECTS (Sobel) ----------
cat("\n########## INDIRECT EFFECTS ##########\n")
getb <- function(f, meth){ x <- read.csv(f); x <- x[x$method==meth, ]; c(b=x$b[1], se=x$se[1], p=x$pval[1]) }
ind <- data.frame()
for (k in meds$key) {
  b1_ivw <- getb(file.path(outdir, paste0("step1_", k, ".csv")), "Inverse variance weighted")
  b1_wm  <- getb(file.path(outdir, paste0("step1_", k, ".csv")), "Weighted median")
  for (oc in c("nielsen","r11")) {
    b2 <- getb(file.path(outdir, paste0("step2_", k, "_", oc, ".csv")), "Inverse variance weighted")
    for (tag in c("ivw","wm")) {
      b1 <- if (tag=="ivw") b1_ivw else b1_wm
      indb <- b1["b"] * b2["b"]
      indse <- sqrt(b2["b"]^2 * b1["se"]^2 + b1["b"]^2 * b2["se"]^2)
      z <- indb/indse; p <- 2*pnorm(-abs(z))
      ind <- rbind(ind, data.frame(mediator=k, outcome=oc, step1_method=tag,
        beta1=b1["b"], se1=b1["se"], p1=b1["p"], beta2=b2["b"], se2=b2["se"], p2=b2["p"],
        indirect=indb, indirect_se=indse, z=z, p=p))
    }
  }
}
print(ind)
write.csv(ind, file.path(outdir, "indirect_effects.csv"), row.names=FALSE)

# ---------- MVMR direct effects: coffee + each mediator ----------
inst <- read.csv("mr_results/instruments.csv")
coffee_exp <- format_data(inst, type="exposure", snp_col="SNP", beta_col="beta.exposure", se_col="se.exposure",
                   effect_allele_col="effect_allele.exposure", other_allele_col="other_allele.exposure",
                   eaf_col="eaf.exposure", pval_col="pval.exposure")
coffee_exp$exposure <- "Coffee intake"
mv_summary <- data.frame()
for (i in seq_len(nrow(meds))) {
  k <- meds$key[i]
  med_exp <- read.csv(file.path(outdir, paste0("iv_", k, ".csv")))
  med_exp$exposure <- meds$name[i]
  both <- rbind(coffee_exp[, intersect(names(coffee_exp), names(med_exp))],
                med_exp[, intersect(names(coffee_exp), names(med_exp))])
  # joint LD clump
  cl <- tryCatch(ld_clump(dplyr::tibble(rsid=both$SNP, pval=both$pval.exposure, id=both$id.exposure),
                          clump_r2=0.001, clump_kb=10000),
                 error=function(e){ cat("clump ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(cl)) next
  both <- both[both$SNP %in% cl$rsid, ]
  cat("\n########## MVMR: coffee +", k, " joint snps:", nrow(both), "##########\n")
  ex <- lapply(split(both, both$exposure), function(x) x)
  mv_exp <- list()
  for (nm in names(ex)) {
    df <- ex[[nm]]
    mv_exp[[nm]] <- df
  }
  # build exposure set in TwoSampleMR mv format
  exposure_dat <- both
  for (oc in c("nielsen","r11")) {
    oname <- ifelse(oc=="nielsen", "AF (Nielsen 2018)", "AF/flutter (FinnGen R11)")
    out <- if (oc=="nielsen") mk_out_nielsen(exposure_dat$SNP, oname) else mk_out_r11(exposure_dat$SNP, oname)
    mvd <- tryCatch(mv_harmonise_data(exposure_dat, out), error=function(e){ cat("mv ERROR:", conditionMessage(e), "\n"); NULL })
    if (is.null(mvd)) next
    res <- tryCatch(mv_multiple(mvd), error=function(e){ cat("mv ERROR:", conditionMessage(e), "\n"); NULL })
    if (!is.null(res)) {
      ror <- generate_odds_ratios(res$result)
      print(ror[, c("exposure","nsnp","b","se","pval","or","or_lci95","or_uci95")])
      ror$mediator <- k; ror$outcome_set <- oc
      mv_summary <- rbind(mv_summary, ror)
    }
    cf <- tryCatch(strength_mvmr(mvd), error=function(e) NULL)
    if (!is.null(cf)) print(cf)
  }
}
write.csv(mv_summary, file.path(outdir, "mvmr_direct_effects.csv"), row.names=FALSE)
cat("\nALL DONE\n")
sink(type="message"); sink(); close(log)
