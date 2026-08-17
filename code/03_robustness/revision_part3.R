# Part 3: fix M2b p-values + M3 robust MVMR (R11, local)
suppressMessages({library(TwoSampleMR); library(MVMR); library(quantreg)})
ws <- "."
med <- file.path(ws, "mediation_mr")
out <- file.path(ws, "mr_results")

## ---- M2b recompute p ----
h1 <- read.csv(file.path(med, "step1_harmonised_bmi.csv")); h1 <- h1[h1$mr_keep==TRUE, ]
fto <- h1$SNP[h1$chr==16 & h1$pos > 52.5e6 & h1$pos < 54.5e6]
presso_out <- c("rs1421085","rs1527961","rs6062682","rs13387939")
sobel <- function(b1, s1, b2, s2){
  est <- unname(b1*b2); se <- unname(sqrt(b2^2*s1^2 + b1^2*s2^2)); z <- est/se
  list(est=est, se=se, z=z, p=2*pnorm(-abs(z)))
}
get2 <- function(f){ x <- read.csv(file.path(med, f)); r <- x[x$method=="Inverse variance weighted", ]; list(b=r$b, se=r$se) }
b2n <- get2("step2_bmi_nielsen.csv"); b2r <- get2("step2_bmi_r11.csv")
sets <- list(all=h1, noFTO=h1[!h1$SNP %in% fto, ], noPRESSO=h1[!h1$SNP %in% presso_out, ])
res_tab <- data.frame()
for (nm in names(sets)) {
  d <- sets[[nm]]
  m <- mr(d, method_list=c("mr_ivw","mr_weighted_median"))
  ivw <- m[m$method=="Inverse variance weighted", ]; wm <- m[m$method=="Weighted median", ]
  a <- sobel(ivw$b, ivw$se, b2n$b, b2n$se); b <- sobel(ivw$b, ivw$se, b2r$b, b2r$se)
  cw <- sobel(wm$b, wm$se, b2n$b, b2n$se);  dd <- sobel(wm$b, wm$se, b2r$b, b2r$se)
  row <- data.frame(set=nm, nsnp=ivw$nsnp, beta1_ivw=ivw$b, se1_ivw=ivw$se, p1_ivw=ivw$pval,
                    beta1_wm=wm$b, p1_wm=wm$pval,
                    ind_nielsen=a$est, p_nielsen=a$p, ind_r11=b$est, p_r11=b$p,
                    ind_nielsen_wm=cw$est, p_nielsen_wm=cw$p, ind_r11_wm=dd$est, p_r11_wm=dd$p)
  res_tab <- rbind(res_tab, row); print(row, digits=4)
}
write.csv(res_tab, file.path(med, "rev_fto_exclusion.csv"), row.names=FALSE)

## ---- M3 robust MVMR function ----
robust_mvmr <- function(mvdat, tag){
  bx <- mvdat$exposure_beta; by <- mvdat$outcome_beta
  se_out <- mvdat$outcome_se
  w <- 1 / se_out^2
  exp_ids <- colnames(bx)
  lab <- c("ieu-b-4877"="Smoking initiation","ukb-b-19953"="BMI","ukb-b-5237"="Coffee intake","ukb-b-5779"="Alcohol frequency")
  coffee_col <- which(exp_ids=="ukb-b-5237")
  cat("\n---- robust MVMR:", tag, " nsnp:", nrow(bx), " coffee col:", coffee_col, "----\n")
  f_ivw <- lm(by ~ bx - 1, weights=w)
  f_med <- tryCatch(rq(by ~ bx - 1, tau=0.5, weights=w), error=function(e){ cat("rq err:", conditionMessage(e), "\n"); NULL })
  s <- sign(bx[, coffee_col]); s[s==0] <- 1
  bxo <- bx * s; byo <- by * s
  f_egg <- lm(byo ~ bxo, weights=w)
  med_sum <- if (!is.null(f_med)) summary(f_med, se="boot", bsmethod="xy", R=500)$coefficients else NULL
  # rank-se summary lacks p column: compute from t (Value/SE)
  wm_p_val <- NULL
  if (!is.null(med_sum)) {
    if (ncol(med_sum) >= 4) { wm_p_val <- med_sum[1:ncol(bx), 4] }
    else { tv <- med_sum[1:ncol(bx),1] / med_sum[1:ncol(bx),2]; wm_p_val <- 2*pnorm(-abs(tv)) }
  }
  n_e <- ncol(bx)
  res <- data.frame(exposure=lab[exp_ids],
                    ivw_b=coef(f_ivw)[1:n_e], ivw_se=coef(summary(f_ivw))[1:n_e,2], ivw_p=coef(summary(f_ivw))[1:n_e,4],
                    wm_b=if(is.null(med_sum)) NA else med_sum[1:n_e,1],
                    wm_se=if(is.null(med_sum)) NA else med_sum[1:n_e,2],
                    wm_p=if(is.null(wm_p_val)) NA else wm_p_val,
                    egger_b=coef(f_egg)[2:(n_e+1)], egger_se=coef(summary(f_egg))[2:(n_e+1),2], egger_p=coef(summary(f_egg))[2:(n_e+1),4])
  res$ivw_or <- exp(res$ivw_b); res$wm_or <- exp(res$wm_b); res$egger_or <- exp(res$egger_b)
  ei <- coef(summary(f_egg))[1, ]
  cat("Egger intercept:", round(ei[1],5), " p:", signif(ei[4],3), "\n")
  print(res[, c("exposure","ivw_or","ivw_p","wm_or","wm_p","egger_or","egger_p")])
  res$dataset <- tag; res$egger_intercept <- ei[1]; res$egger_intercept_p <- ei[4]
  write.csv(res, file.path(med, paste0("rev_mvmr_robust_", tag, ".csv")), row.names=FALSE)
  # manual Q_A (residual Cochran Q of IVW fit)
  fit_b <- coef(f_ivw)[1:n_e]
  resid <- by - as.vector(bx %*% fit_b)
  Q <- sum(w * resid^2); dfq <- nrow(bx) - n_e
  cat("-- manual Q_A:", round(Q,1), " df:", dfq, " p:", signif(pchisq(Q, dfq, lower.tail=FALSE),3), "\n")
  cat("-- conditional F --\n"); print(tryCatch(strength_mvmr(mvdat), error=function(e) conditionMessage(e)))
  cat("-- pleiotropy_mvmr --\n"); print(tryCatch(pleiotropy_mvmr(mvdat), error=function(e) conditionMessage(e)))
  invisible(res)
}

mv_exp <- readRDS(file.path(out, "mv_exp.rds"))
read_r11 <- function(f){ x <- read.delim(f, comment.char="", check.names=FALSE); names(x)[1] <- "chrom"; x }
r11m <- read_r11(file.path(ws, "r11_mvmr_hits.tsv"))
ids <- unique(mv_exp$SNP)
pick <- sapply(strsplit(r11m$rsids, ","), function(v){ m <- v[v %in% ids]; if(length(m)) m[1] else NA })
r11m$SNP <- pick; r11m <- r11m[!is.na(r11m$SNP), ]
o <- format_data(data.frame(SNP=r11m$SNP, beta=r11m$beta, se=r11m$sebeta, eaf=r11m$af_alt,
                            effect_allele=toupper(r11m$alt), other_allele=toupper(r11m$ref), pval=r11m$pval),
                 type="outcome")
o$outcome <- "AF & flutter (FinnGen R11)"
mvdat_r11 <- mv_harmonise_data(mv_exp, o)
robust_mvmr(mvdat_r11, "r11")
cat("\nPART3 DONE\n")
