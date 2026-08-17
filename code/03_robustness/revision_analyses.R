# Revision analyses for integrated coffee-AF manuscript
# M2a: bidirectional MR (BMI -> coffee)
# M2b: FTO-excluded sensitivity for coffee -> BMI step1 + indirect effects
# M3: robust MVMR (Egger + weighted median) for confounder MVMR
jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr); library(MVMR); library(quantreg)})
ws <- "."
med <- file.path(ws, "mediation_mr")
out <- file.path(ws, "mr_results")
log <- file(file.path(med, "revision_log.txt"), "w"); sink(log, type="output"); sink(log, type="message")

## ---------------- M2a: bidirectional MR, BMI -> coffee ----------------
cat("\n========== M2a: BMI -> coffee (bidirectional) ==========\n")
ivb <- read.csv(file.path(med, "iv_bmi.csv"))
exp_bmi <- format_data(ivb, type="exposure", snp_col="SNP", beta_col="beta.exposure",
                       se_col="se.exposure", effect_allele_col="effect_allele.exposure",
                       other_allele_col="other_allele.exposure", eaf_col="eaf.exposure",
                       pval_col="pval.exposure")
exp_bmi$exposure <- "BMI (ieu-b-40)"
snps <- unique(exp_bmi$SNP)
cat("BMI instruments:", length(snps), "\n")
chunks <- split(snps, ceiling(seq_along(snps)/150))
out_list <- lapply(chunks, function(s) extract_outcome_data(snps=s, outcomes="ukb-b-5237"))
out_coffee <- do.call(rbind, out_list)
out_coffee$outcome <- "Coffee intake (ukb-b-5237)"
cat("outcome rows:", nrow(out_coffee), "\n")
d_bc <- harmonise_data(exp_bmi, out_coffee, action=2)
d_bc <- d_bc[d_bc$mr_keep==TRUE, ]
cat("harmonised kept:", nrow(d_bc), "\n")
res_bc <- mr(d_bc, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
print(res_bc[, c("method","nsnp","b","se","pval")])
write.csv(res_bc, file.path(med, "rev_bidir_bmi_to_coffee.csv"), row.names=FALSE)
eg_int <- mr_pleiotropy_test(d_bc); print(eg_int)
write.csv(eg_int, file.path(med, "rev_bidir_egger_intercept.csv"), row.names=FALSE)

## ---------------- M2b: FTO-excluded coffee -> BMI step1 + indirect ----------------
cat("\n========== M2b: FTO-excluded sensitivity ==========\n")
h1 <- read.csv(file.path(med, "step1_harmonised_bmi.csv"))
h1 <- h1[h1$mr_keep==TRUE, ]
cat("step1 coffee->BMI SNPs:", nrow(h1), "\n")
print(h1[h1$chr==16, c("SNP","chr","pos")])
fto <- h1$SNP[h1$chr==16 & h1$pos > 52.5e6 & h1$pos < 54.5e6]
presso_out <- c("rs1421085","rs1527961","rs6062682","rs13387939")
cat("FTO-region SNPs in coffee instruments:", paste(fto, collapse=", "), "\n")

sobel <- function(b1, s1, b2, s2){
  est <- b1*b2; se <- sqrt(b2^2*s1^2 + b1^2*s2^2); z <- est/se
  c(est=est, se=se, z=z, p=2*pnorm(-abs(z)))
}
# step2 BMI -> AF IVW estimates
s2n <- read.csv(file.path(med, "step2_bmi_nielsen.csv")); s2r <- read.csv(file.path(med, "step2_bmi_r11.csv"))
get2 <- function(x){ r <- x[x$method=="Inverse variance weighted", ]; c(b=r$b, se=r$se) }
b2n <- get2(s2n); b2r <- get2(s2r)

run_step1 <- function(d){
  m <- mr_ivw(mr_input(bx=d$beta.exposure, bxse=d$se.exposure, by=d$beta.outcome, byse=d$se.outcome))
  c(b=m$Estimate, se=m$StdError)
}
sets <- list(all=h1,
             noFTO=h1[!h1$SNP %in% fto, ],
             noPRESSO=h1[!h1$SNP %in% presso_out, ])
res_tab <- data.frame()
for (nm in names(sets)) {
  d <- sets[[nm]]
  b1 <- run_step1(d)
  in_n <- sobel(b1["b"], b1["se"], b2n["b"], b2n["se"])
  in_r <- sobel(b1["b"], b1["se"], b2r["b"], b2r["se"])
  row <- data.frame(set=nm, nsnp=nrow(d), beta1=b1["b"], se1=b1["se"],
                    ind_nielsen=in_n["est"], p_nielsen=in_n["p"],
                    ind_r11=in_r["est"], p_r11=in_r["p"])
  res_tab <- rbind(res_tab, row); print(row)
}
write.csv(res_tab, file.path(med, "rev_fto_exclusion.csv"), row.names=FALSE)

## ---------------- M3: robust MVMR (confounder model) ----------------
cat("\n========== M3: robust MVMR ==========\n")
mv_exp <- readRDS(file.path(out, "mv_exp.rds"))

read_r11 <- function(f){ x <- read.delim(f, comment.char="", check.names=FALSE); names(x)[1] <- "chrom"; x }
mk_outcome <- function(x, ids, name){
  pick <- sapply(strsplit(x$rsids, ","), function(v){ m <- v[v %in% ids]; if(length(m)) m[1] else NA })
  x$SNP <- pick; x <- x[!is.na(x$SNP), ]
  o <- format_data(data.frame(SNP=x$SNP, beta=x$beta, se=x$sebeta, eaf=x$af_alt,
                              effect_allele=toupper(x$alt), other_allele=toupper(x$ref), pval=x$pval),
                   type="outcome", snp_col="SNP", beta_col="beta", se_col="se",
                   effect_allele_col="effect_allele", other_allele_col="other_allele",
                   eaf_col="eaf", pval_col="pval")
  o$outcome <- name; o
}

robust_mvmr <- function(mvdat, tag){
  cat("\n---- robust MVMR:", tag, "----\n")
  bx <- mvdat$exposure_beta; by <- mvdat$outcome_beta
  w <- 1 / mvdat$outcome_se^2
  expnames <- colnames(bx)
  cat("exposures:", paste(expnames, collapse=", "), "\n")
  # IVW (reference)
  f_ivw <- lm(by ~ bx - 1, weights=w)
  # multivariable weighted median (LAD, tau=0.5)
  f_med <- tryCatch(rq(by ~ bx - 1, tau=0.5, weights=w), error=function(e) NULL)
  # MVMR-Egger oriented by coffee sign
  coffee_col <- grep("[Cc]offee", expnames)[1]
  s <- sign(bx[, coffee_col]); s[s==0] <- 1
  bxo <- bx * s; byo <- by * s
  f_egg <- lm(byo ~ bxo, weights=w)  # with intercept
  res <- data.frame(exposure=expnames,
                    ivw_b=coef(f_ivw), ivw_se=coef(summary(f_ivw))[,2],
                    wm_b=if(is.null(f_med)) NA else coef(f_med),
                    wm_se=if(is.null(f_med)) NA else summary(f_med, se="xy")$coefficients[,2],
                    egger_b=coef(f_egg)[-1], egger_se=coef(summary(f_egg))[-1,2])
  res$ivw_or <- exp(res$ivw_b); res$wm_or <- exp(res$wm_b); res$egger_or <- exp(res$egger_b)
  res$ivw_p <- coef(summary(f_ivw))[,4]
  res$wm_p <- if(is.null(f_med)) NA else summary(f_med, se="xy")$coefficients[,4]
  res$egger_p <- coef(summary(f_egg))[-1,4]
  egg_int <- coef(summary(f_egg))[1, ]; cat("Egger intercept:", egg_int[1], "p:", egg_int[4], "\n")
  print(res[, c("exposure","ivw_or","ivw_p","wm_or","wm_p","egger_or","egger_p")])
  res$dataset <- tag; res$egger_intercept <- egg_int[1]; res$egger_intercept_p <- egg_int[4]
  write.csv(res, file.path(med, paste0("rev_mvmr_robust_", tag, ".csv")), row.names=FALSE)
  cat("-- conditional F --\n"); print(tryCatch(strength_mvmr(mvdat), error=function(e) conditionMessage(e)))
  cat("-- pleiotropy_mvmr (Q_A) --\n"); print(tryCatch(pleiotropy_mvmr(mvdat), error=function(e) conditionMessage(e)))
  invisible(res)
}

# R11 outcome
r11m <- read_r11(file.path(ws, "r11_mvmr_hits.tsv"))
out_r11 <- mk_outcome(r11m, unique(mv_exp$SNP), "AF & flutter (FinnGen R11)")
mvdat_r11 <- mv_harmonise_data(mv_exp, out_r11)
robust_mvmr(mvdat_r11, "r11")

# Nielsen outcome
od_n <- extract_outcome_data(snps=unique(mv_exp$SNP), outcomes="ebi-a-GCST006414")
od_n$outcome <- "AF (Nielsen 2018)"
mvdat_n <- mv_harmonise_data(mv_exp, od_n)
robust_mvmr(mvdat_n, "nielsen")

cat("\nALL DONE\n")
sink(type="message"); sink(); close(log)
