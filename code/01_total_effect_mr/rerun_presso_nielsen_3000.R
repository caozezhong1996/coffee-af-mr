# P1-5 precision confirmation: NbDistribution = 3000 (run via script file)
.libPaths(c("C:/Users/曹泽众/Documents/kimi/workspace/Rlib", .libPaths()))
suppressMessages(library(MRPRESSO))
d <- read.csv("D:/Desktop/EJPC_投稿_咖啡与房颤/03_Supplementary/S1_total_effect_MR/harmonised_ebi-a-GCST006414.csv")
d <- d[d$mr_keep == TRUE, ]
pr <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                SdOutcome = "se.outcome", SdExposure = "se.exposure", data = d,
                OUTLIERtest = TRUE, DISTORTIONtest = TRUE, SignifThreshold = 0.05,
                NbDistribution = 3000, seed = 20260816)
ot <- cbind(SNP = d$SNP, pr$`MR-PRESSO results`$`Outlier Test`)
ot <- ot[order(ot$Pvalue), ]
cat("Global P:", pr$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
print(head(ot, 5))
write.csv(ot, "C:/Users/曹泽众/Documents/kimi/workspace/EJPC_修订产出/stats/presso_rerun_nielsen_nb3000.csv",
          row.names = FALSE)
cat("Distortion P:", pr$`MR-PRESSO results`$`Distortion Test`$Pvalue, "\n")
cat("DONE\n")
