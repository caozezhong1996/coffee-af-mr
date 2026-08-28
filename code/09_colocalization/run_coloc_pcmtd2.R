suppressMessages({library(data.table); library(coloc)})
base <- "C:/Users/64915/Documents/kimi/workspace/integrated_v4/coloc_v42"
sub  <- fread(file.path(base, "fg_pcmtd2_fix.tsv"))
nei_loc <- fread(file.path(base, "nei_windows_500kb.tsv"))[locus == "PCMTD2"]
cf   <- fread(file.path(base, "coffee_pcmtd2_fix.csv"))
cf_all <- fread(file.path(base, "coffee_vcf_associations.csv"))

comp <- function(a) chartr("ACGT", "TGCA", a)
is_pal <- function(a1, a2) (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C")
harm <- function(eff_a, e1, e2, fr, ref, alt) {
  n <- length(eff_a); beta <- rep(NA_real_, n); eaf <- rep(NA_real_, n); ok <- rep(FALSE, n)
  for (i in seq_len(n)) {
    a1 <- e1[i]; a2 <- e2[i]; f <- fr[i]; rf <- ref[i]; al <- alt[i]
    if (is.na(a1) || is.na(a2) || is.na(rf) || is.na(al)) next
    if (a1 == al && a2 == rf)      { beta[i] <-  eff_a[i]; eaf[i] <- f; ok[i] <- TRUE }
    else if (a1 == rf && a2 == al) { beta[i] <- -eff_a[i]; eaf[i] <- if (is.na(f)) NA else 1 - f; ok[i] <- TRUE }
    else if (comp(a1) == al && comp(a2) == rf) { beta[i] <-  eff_a[i]; eaf[i] <- f; ok[i] <- TRUE }
    else if (comp(a1) == rf && comp(a2) == al) { beta[i] <- -eff_a[i]; eaf[i] <- if (is.na(f)) NA else 1 - f; ok[i] <- TRUE }
  }
  list(beta = beta, eaf = eaf, ok = ok)
}
maf_from <- function(e) ifelse(is.na(e), NA_real_, pmin(e, 1 - e))
cf_all[, m := maf_from(eaf)]
sdy_est <- median(cf_all$se * sqrt(2 * cf_all$n * cf_all$m * (1 - cf_all$m)), na.rm = TRUE)

fg_maf <- maf_from(sub$af_alt)
csub <- cf[match(sub$rsids, rsid)]
h <- harm(csub$beta, csub$ea, csub$nea, csub$eaf, sub$ref, sub$alt)
keep_c <- h$ok & !is.na(csub$beta) & !is.na(sub$beta) & csub$se > 0
pal <- is_pal(sub$ref, sub$alt)
caf_maf <- maf_from(h$eaf); caf_maf[is.na(caf_maf)] <- fg_maf[is.na(caf_maf)]
keep_c <- keep_c & !(pal & caf_maf > 0.42) & !is.na(caf_maf)
d1_all <- list(snp = sub$rsids, beta = h$beta, varbeta = (csub$se)^2,
               type = "quant", N = 428860, sdY = sdy_est, MAF = caf_maf)

nidx <- match(sub$rsids, nei_loc$rs_dbSNP147)
cat("PCMTD2 rsid direct match:", sum(!is.na(nidx)), "/", nrow(sub), "\n")
miss <- is.na(nidx) & !is.na(csub$position)
pos_map <- match(csub$position[miss], nei_loc$POS_GRCh37)
nidx[miss] <- pos_map
cat("position fallback hits:", sum(!is.na(pos_map)), "\n")
nsub <- nei_loc[nidx]
hn <- harm(nsub$Effect_A2, nsub$A2, nsub$A1, suppressWarnings(as.numeric(nsub$Freq_A2)), sub$ref, sub$alt)
n_maf <- maf_from(hn$eaf); n_maf[is.na(n_maf)] <- caf_maf[is.na(n_maf)]

run_one <- function(d1, d2, p12) {
  res <- tryCatch(coloc.abf(d1, d2, p12 = p12), error = function(e) e)
  if (inherits(res, "error")) return(NULL)
  sm <- res$summary; tr <- res$results
  ppcol <- grep("SNP.PP", names(tr), value = TRUE)
  top_snp <- NA_character_; top_pp <- NA_real_
  if (length(ppcol) > 0) { top <- tr[which.max(tr[[ppcol[1]]]), ]; top_snp <- as.character(top$snp); top_pp <- as.numeric(top[[ppcol[1]]]) }
  list(pp = sm[c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")], top_snp = top_snp, top_pp = top_pp, nsnp = sm["nsnps"])
}
res <- list(); sens <- list()
for (tag in c("FinnGenR11", "Nielsen2018")) {
  if (tag == "FinnGenR11") { af_beta <- sub$beta; af_se <- sub$sebeta; af_maf <- fg_maf; N_af <- 287805; s_af <- 55853/287805 }
  else { af_beta <- hn$beta; af_se <- nsub$StdErr; af_maf <- n_maf; N_af <- 1030836; s_af <- 60620/1030836 }
  k <- keep_c & !is.na(af_beta) & !is.na(af_se) & af_se > 0 & !is.na(af_maf) & af_maf > 0 & af_maf < 0.5
  idx <- which(k)
  cat("PCMTD2", tag, "usable SNPs:", length(idx), "\n")
  if (length(idx) < 50) next
  d1 <- lapply(d1_all, function(x) if (is.atomic(x) && length(x) == nrow(sub)) x[idx] else x)
  d2 <- list(snp = sub$rsids[idx], beta = af_beta[idx], varbeta = (af_se[idx])^2, type = "cc", N = N_af, s = s_af, MAF = af_maf[idx])
  r <- run_one(d1, d2, 1e-5)
  res[[tag]] <- data.table(locus = "PCMTD2", outcome = tag, nsnps = r$nsnp,
    PP.H0 = r$pp[1], PP.H1 = r$pp[2], PP.H2 = r$pp[3], PP.H3 = r$pp[4], PP.H4 = r$pp[5],
    top_snp = r$top_snp, top_SNP.PP.H4 = r$top_pp)
  for (p12 in c(5e-6, 5e-5)) {
    r2 <- run_one(d1, d2, p12)
    sens[[paste(tag, p12)]] <- data.table(locus = "PCMTD2", outcome = tag, p12 = p12, nsnps = r2$nsnp, PP.H3 = r2$pp[4], PP.H4 = r2$pp[5])
  }
}
out <- rbindlist(res); sen <- rbindlist(sens)
fwrite(out, file.path(base, "coloc_pcmtd2_fix_results.csv"))
fwrite(sen, file.path(base, "coloc_pcmtd2_fix_sensitivity.csv"))
print(out[, .(locus, outcome, nsnps, PP.H1 = round(PP.H1,3), PP.H3 = round(PP.H3,3), PP.H4 = round(PP.H4,3), top_snp, top_SNP.PP.H4 = round(top_SNP.PP.H4,3))])
print(sen)
