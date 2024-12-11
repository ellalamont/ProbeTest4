# Compare library prep (full or half) to Mtb # reads
# E. Lamont
# 12/11/24

source("Import_data.R") # to get my_pipeSummary

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
        plot.margin = margin(10, 10, 10, 20))



###########################################################
######## F: THP1_1e6_1 DUALrRNA NUMBER+PERCENT READS ######

# Number Reads
THP1_LibPrepVsReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("originalTHP1_1e6_1_DualrRNA", "originalTHP1_1e6_2_DualrRNA", "originalTHP1_1e6_3_DualrRNA", "originalTHP1_1e6_4_DualrRNA", "originalTHP1_1e6_5_DualrRNA", "originalTHP1_1e6_6_DualrRNA")) %>% 
  ggplot(aes(x = Library_prep, y = N_Genomic)) + 
  geom_point(color = "#FF7F00", size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "Full/Half library prep vs number reads aligned to Mtb",
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra; DualrRNADepletion \nHalf library prep also was re-hybridized") + 
  my_plot_themes
THP1_LibPrepVsReads
ggsave(THP1_LibPrepVsReads,
       file = "THP1_LibPrepVsReads.pdf",
       path = "LibraryPrep_Figures",
       width = 6, height = 4, units = "in")

# Percent Reads
THP1_LibPrepVsPercentReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("originalTHP1_1e6_1_DualrRNA", "originalTHP1_1e6_2_DualrRNA", "originalTHP1_1e6_3_DualrRNA", "originalTHP1_1e6_4_DualrRNA", "originalTHP1_1e6_5_DualrRNA", "originalTHP1_1e6_6_DualrRNA")) %>% 
  ggplot(aes(x = Library_prep, y = P_Genomic)) + 
  geom_point(color = "#FF7F00", size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "Full/Half library prep vs percent reads aligned to Mtb",
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra; DualrRNADepletion \nHalf library prep also was re-hybridized") + 
  my_plot_themes
THP1_LibPrepVsPercentReads
ggsave(THP1_LibPrepVsPercentReads,
       file = "THP1_LibPrepVsPercentReads.pdf",
       path = "LibraryPrep_Figures",
       width = 6, height = 4, units = "in")











