# Revision analyses part 2 (local): M2b fixed (use TwoSampleMR::mr) + M3 R11 robust MVMR
suppressMessages({library(TwoSampleMR); library(MVMR); library(quantreg)})
ws <- "."
med <- file.path(ws, "mediation_mr")
out <- file.path(ws, "mr_results")
log <- file(file.path(med, "revision_log2.txt"), "w"); sink(log, type="output"); sink(log, type="message")

## ---------------- M2b: FTO-excluded coffee -> BMI ----------------
cat("========== M2b: FTO-excluded ==========\n")
h1 <- read.csv(file.path(med, "step1_harmonised_bmi.csv"))
h1 <- h1[h1$mr_keep==TRUE, ]
fto <- h1$SNP[h1$chr==16 & h1$pos > 52.5e6 & h1$pos < 54.5e6]
presso_out <- c("rs1421085","rs1527961","rs6062682","rs13387939")
cat("FTO-region SNPs:", paste(fto, collapse=", "), "\n")
sobel <- function(b1, s1, b2, s2){
  est <- b1*b2; se <- sqrt(b2^2*s1^2 + b1^2*s2^2); z <- est/se
  c(est=unname(est), se=unname(se), z=unname(z), p=2*pnorm(-abs(z)))
}
s2n <- read.csv(file.path(med, "step2_bmi_nielsen.csv")); s2r <- read.csv(file.path(med, "step2_bmi_r11.csv"))
get2 <- function(x){ r <- x[x$method=="Inverse variance weighted", ]; c(b=r$b, se=r$se) }
b2n <- get2(s2n); b2r <- get2(s2r)
sets <- list(all=h1, noFTO=h1[!h1$SNP %in% fto, ], noPRESSO=h1[!h1$SNP %in% presso_out, ])
res_tab <- data.frame()
for (nm in names(sets)) {
  d <- sets[[nm]]
  m <- mr(d, method_list=c("mr_ivw","mr_weighted_median"))
  ivw <- m[m$method=="Inverse variance weighted", ]
  wm  <- m[m$method=="Weighted median", ]
  in_n  <- sobel(ivw$b, ivw$se, b2n["b"], b2n["se"])
  in_r  <- sobel(ivw$b, ivw$se, b2r["b"], b2r["se"])
  in_nw <- sobel(wm$b, wm$se, b2n["b"], b2n["se"])
  in_rw <- sobel(wm$b, wm$se, b2r["b"], b2r["se"])
  row <- data.frame(set=nm, nsnp=ivw$nsnp, beta1_ivw=ivw$b, se1_ivw=ivw$se, p1_ivw=ivw$pval,
                    beta1_wm=wm$b, p1_wm=wm$pval,
                    ind_nielsen=in_n["est"], p_nielsen=in_n["p"],
                    ind_r11=in_r["est"], p_r11=in_r["p"],
                    ind_nielsen_wm=in_nw["est"], p_nielsen_wm=in_nw["p"],
                    ind_r11_wm=in_rw["est"], p_r11_wm=in_rw["p"])
  res_tab <- rbind(res_tab, row); print(row)
}
write.csv(res_tab, file.path(med, "rev_fto_exclusion.csv"), row.names=FALSE)

## ---------------- M3 (R11): robust MVMR ----------------
cat("\n========== M3 R11: robust MVMR ==========\n")
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
  bx <- as.matrix(mvdat$exposure_beta); by <- mvdat$outcome_beta
  w <- 1 / mvdat$outcome_se^2
  expnames <- colnames(bx)
  cat("exposures:", paste(expnames, collapse=", "), " nsnp:", nrow(bx), "\n")
  f_ivw <- lm(by ~ bx - 1, weights=w)
  f_med <- tryCatch(rq(by ~ bx - 1, tau=0.5, weights=w), error=function(e) { cat("rq err:", conditionMessage(e), "\n"); NULL })
  coffee_col <- grep("[Cc]offee", expnames)[1]
  s <- sign(bx[, coffee_col]); s[s==0] <- 1
  bxo <- bx * s; byo <- by * s
  f_egg <- lm(byo ~ bxo, weights=w)
  med_sum <- if (!is.null(f_med)) summary(f_med, se="xy")$coefficients else NULL
  res <- data.frame(exposure=expnames,
                    ivw_b=coef(f_ivw), ivw_se=coef(summary(f_ivw))[,2], ivw_p=coef(summary(f_ivw))[,4],
                    wm_b=if(is.null(med_sum)) NA else med_sum[,1],
                    wm_se=if(is.null(med_sum)) NA else med_sum[,2],
                    wm_p=if(is.null(med_sum)) NA else med_sum[,4],
                    egger_b=coef(f_egg)[-1], egger_se=coef(summary(f_egg))[-1,2], egger_p=coef(summary(f_egg))[-1,4])
  res$ivw_or <- exp(res$ivw_b); res$wm_or <- exp(res$wm_b); res$egger_or <- exp(res$egger_b)
  ei <- coef(summary(f_egg))[1, ]
  cat("Egger intercept:", ei[1], " p:", ei[4], "\n")
  print(res[, c("exposure","ivw_or","ivw_p","wm_or","wm_p","egger_or","egger_p")])
  res$dataset <- tag; res$egger_intercept <- ei[1]; res$egger_intercept_p <- ei[4]
  write.csv(res, file.path(med, paste0("rev_mvmr_robust_", tag, ".csv")), row.names=FALSE)
  cat("-- conditional F --\n"); print(tryCatch(strength_mvmr(mvdat), error=function(e) conditionMessage(e)))
  cat("-- pleiotropy_mvmr --\n"); print(tryCatch(pleiotropy_mvmr(mvdat), error=function(e) conditionMessage(e)))
  invisible(res)
}
r11m <- read_r11(file.path(ws, "r11_mvmr_hits.tsv"))
out_r11 <- mk_outcome(r11m, unique(mv_exp$SNP), "AF & flutter (FinnGen R11)")
mvdat_r11 <- mv_harmonise_data(mv_exp, out_r11)
robust_mvmr(mvdat_r11, "r11")

cat("\nPART2 DONE\n")
sink(type="message"); sink(); close(log)
