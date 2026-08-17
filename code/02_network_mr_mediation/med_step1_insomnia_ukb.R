jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages(library(TwoSampleMR))
outdir <- "./mediation_mr"
inst <- read.csv("mr_results/instruments.csv")
exp <- format_data(inst, type="exposure", snp_col="SNP", beta_col="beta.exposure", se_col="se.exposure",
                   effect_allele_col="effect_allele.exposure", other_allele_col="other_allele.exposure",
                   eaf_col="eaf.exposure", pval_col="pval.exposure")
exp$exposure <- "Coffee intake (cups/day, UKB)"
out <- extract_outcome_data(snps=exp$SNP, outcomes="ukb-b-3957")
out$outcome <- "Sleeplessness/insomnia (UKB)"
dat <- harmonise_data(exp, out, action=2)
d <- dat[dat$mr_keep==TRUE, ]
res <- mr(d, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
print(res[, c("method","nsnp","b","se","pval")])
write.csv(res, file.path(outdir, "step1_insomnia_ukb.csv"), row.names=FALSE)
