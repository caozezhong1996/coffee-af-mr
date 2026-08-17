suppressMessages(library(ggplot2))
df <- data.frame(
  label = factor(c("Caffeine (all products)","Caffeine citrate","Taurine","5-hour ENERGY","MONSTER ENERGY","RED BULL"),
                 levels=rev(c("Caffeine (all products)","Caffeine citrate","Taurine","5-hour ENERGY","MONSTER ENERGY","RED BULL"))),
  a = c(134, 48, 6, 103, 53, 50),
  n = c(24546, 32355, 840, 47007, 32991, 32374),
  ror = c(0.67, 0.18, 0.88, 0.27, 0.20, 0.19),
  lo  = c(0.57, 0.14, 0.39, 0.22, 0.15, 0.14),
  hi  = c(0.80, 0.24, 1.97, 0.33, 0.26, 0.25))
df$txt <- sprintf("%d / %d    %.2f (%.2f-%.2f)", df$a, df$n, df$ror, df$lo, df$hi)
p <- ggplot(df, aes(y=label, x=ror)) +
  geom_vline(xintercept=1, linetype=2, color="grey50") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.2, linewidth=0.6) +
  geom_point(size=2.5, color="#1a5276") +
  geom_text(aes(x=4.2, label=txt), size=3, hjust=0) +
  scale_x_log10(limits=c(0.08, 30), breaks=c(0.1,0.25,0.5,1,2,4)) +
  labs(x="Reporting odds ratio (log scale, 95% CI)", y=NULL,
       caption="Counts shown as AF reports / total reports. ROR < 1 = no disproportionality signal.") +
  theme_bw(base_size=10) + theme(plot.caption=element_text(size=7, color="grey40"))
ggsave("./mr_results/Figure5_faers_forest_300dpi.tiff", p, width=8.5, height=4, dpi=300, device="tiff", compression="lzw")
ggsave("./mr_results/Figure5_faers_forest_300dpi.png", p, width=8.5, height=4, dpi=300)
cat("saved\n")
