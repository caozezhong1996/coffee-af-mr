# Part 4: Nielsen outcome extraction (chunked, resumable) + robust MVMR
jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(MVMR); library(quantreg)})
ws <- "."
med <- file.path(ws, "mediation_mr")
mv_exp <- readRDS(file.path(ws, "mr_results", "mv_exp.rds"))
snps <- unique(mv_exp$SNP)
chunks <- split(snps, ceiling(seq_along(snps)/150))
prog_file <- file.path(med, "nielsen_mvmr_progress.rds")
done <- if (file.exists(prog_file)) readRDS(prog_file) else list()
for (i in names(chunks)) {
  if (!is.null(done[[i]])) next
  cat("chunk", i, "/", length(chunks), "\n")
  a <- tryCatch(extract_outcome_data(snps=chunks[[i]], outcomes="ebi-a-GCST006414"),
                error=function(e){ cat("ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(a)) { done[[i]] <- a; saveRDS(done, prog_file) }
}
if (length(done) < length(chunks)) { cat("INCOMPLETE:", length(done), "/", length(chunks), "\n"); quit(save="no") }
od_n <- do.call(rbind, done); od_n$outcome <- "AF (Nielsen 2018)"
cat("Nielsen outcome rows:", nrow(od_n), "\n")

robust_mvmr <- function(mvdat, tag){
  bx <- mvdat$exposure_beta; by <- mvdat$outcome_beta; se_out <- mvdat$outcome_se
  w <- 1 / se_out^2
  exp_ids <- colnames(bx)
  lab <- c("ieu-b-4877"="Smoking initiation","ukb-b-19953"="BMI","ukb-b-5237"="Coffee intake","ukb-b-5779"="Alcohol frequency")
  coffee_col <- which(exp_ids=="ukb-b-5237")
  cat("\n---- robust MVMR:", tag, " nsnp:", nrow(bx), "----\n")
  f_ivw <- lm(by ~ bx - 1, weights=w)
  f_med <- tryCatch(rq(by ~ bx - 1, tau=0.5, weights=w), error=function(e){ cat("rq err:", conditionMessage(e), "\n"); NULL })
  s <- sign(bx[, coffee_col]); s[s==0] <- 1
  f_egg <- lm(bx*s ~ 0)  # placeholder no-op
  bxo <- bx * s; byo <- by * s
  f_egg <- lm(byo ~ bxo, weights=w)
  med_sum <- if (!is.null(f_med)) summary(f_med, se="boot", bsmethod="xy", R=500)$coefficients else NULL
  n_e <- ncol(bx)
  wm_p_val <- NULL
  if (!is.null(med_sum)) {
    if (ncol(med_sum) >= 4) { wm_p_val <- med_sum[1:n_e, 4] }
    else { tv <- med_sum[1:n_e,1] / med_sum[1:n_e,2]; wm_p_val <- 2*pnorm(-abs(tv)) }
  }
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
  fit_b <- coef(f_ivw)[1:n_e]; resid <- by - as.vector(bx %*% fit_b)
  Q <- sum(w * resid^2); dfq <- nrow(bx) - n_e
  cat("-- manual Q_A:", round(Q,1), " df:", dfq, " p:", signif(pchisq(Q, dfq, lower.tail=FALSE),3), "\n")
  invisible(res)
}
mvdat_n <- mv_harmonise_data(mv_exp, od_n)
robust_mvmr(mvdat_n, "nielsen")
cat("\nPART4 DONE\n")
