suppressMessages(library(ggplot2))
rows <- list(
  c("Caffeine (all products)", "Primary", 131, 20679, 0.807, 0.680, 0.959),
  c("Caffeine (all products)", "Deduplicated", 54, 13040, 0.598, 0.457, 0.781),
  c("Caffeine (all products)", "Masking-corrected", 92, 19897, 0.653, 0.532, 0.801),
  c("Caffeine citrate", "Primary", 11, 2571, 0.544, 0.301, 0.984),
  c("Caffeine citrate", "Deduplicated", 4, 1395, 0.413, 0.155, 1.103),
  c("Caffeine citrate", "Masking-corrected", 11, 2514, 0.618, 0.342, 1.117),
  c("Taurine", "Primary", 6, 591, 1.300, 0.581, 2.905),
  c("Taurine", "Deduplicated", 3, 313, 1.392, 0.446, 4.339),
  c("Taurine", "Masking-corrected", 4, 556, 1.019, 0.381, 2.725))
df <- data.frame(
  group = factor(sapply(rows, `[[`, 1), levels=c("Taurine","Caffeine citrate","Caffeine (all products)")),
  spec  = factor(sapply(rows, `[[`, 2), levels=c("Masking-corrected","Deduplicated","Primary")),
  a = as.numeric(sapply(rows, `[[`, 3)), n = as.numeric(sapply(rows, `[[`, 4)),
  ror = as.numeric(sapply(rows, `[[`, 5)), lo = as.numeric(sapply(rows, `[[`, 6)), hi = as.numeric(sapply(rows, `[[`, 7)))
df$txt <- sprintf("%d / %d   %.2f (%.2f-%.2f)", df$a, df$n, df$ror, df$lo, df$hi)
shapes <- c("Primary"=16, "Deduplicated"=17, "Masking-corrected"=15)
p <- ggplot(df, aes(y=group, x=ror, shape=spec)) +
  geom_vline(xintercept=1, linetype=2, color="grey50") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.25, linewidth=0.5,
                 position=position_dodge(width=0.55)) +
  geom_point(size=2.3, color="#1a5276", position=position_dodge(width=0.55)) +
  scale_shape_manual(values=shapes) +
  scale_x_log10(limits=c(0.1, 6), breaks=c(0.1,0.25,0.5,1,2,4)) +
  labs(x="Reporting odds ratio for atrial fibrillation (log scale, 95% CI)", y=NULL, shape="Analysis",
       caption="Background: 11,882,970 serious FAERS reports (93,053 AF). Deduplicated: case-version deduplication (safetyreportversion >= 2 removed).\nMasking-corrected: reports with concomitant pro-arrhythmic drugs excluded. Energy-drink brands: only 22 serious reports, 0 AF (non-informative, not shown).") +
  theme_bw(base_size=10) + theme(legend.position="bottom", plot.caption=element_text(size=6.5, color="grey40"))
ggsave("./mr_results/Figure5_faers_forest_v2_300dpi.tiff", p, width=8.5, height=5, dpi=300, device="tiff", compression="lzw")
ggsave("./mr_results/Figure5_faers_forest_v2_300dpi.png", p, width=8.5, height=5, dpi=300)
cat("saved\n")
