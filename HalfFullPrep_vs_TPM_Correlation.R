# Compare half or full library prep to TPM with correlation plot
# E. Lamont
# 12/12/24

source("Import_data.R") # to get my_tpm
my_tpm$Gene <- rownames(my_tpm)

# Log10 transform the data
my_tpm_Log10 <- my_tpm %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )

# http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2


###############################################################################
############## orginalTHP1_1e6_3_MtbrRNA vs THP1_1e6_5 ########################
# Neither sample was in a pool
# Both samples are the original THP1 that were MtbrRNA depleted by JA and Melanie

# orginalTHP1_1e6_3_MtbrRNA - This sample was Half library prep done by JA in Nov, AND the hybridization was stopped then restarted
# THP1_1e6_5 - The sample was fully library prepped for the september sequencing AND the hydridization was good 16hr

Sample1 <- "orginalTHP1_1e6_3_MtbrRNA"
Sample2 <- "THP1_1e6_5"
ScatterCorr <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = paste0(Sample1, " vs ", Sample2),
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, ".pdf"),
       path = "HalfFullPrep_Figures",
       width = 7, height = 5, units = "in")










