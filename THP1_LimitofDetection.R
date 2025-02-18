# Limit of Detection with the THP1 spiked samples graphs
# E. Lamont
# 1/15/25

source("Import_data.R") # to get my_pipeSummary


# Only looking at the not pooled samples because it didn't work well when I pooled the high and low concentration samples together

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )


###########################################################
############## FILTER JUST THE SAMPLES I WANT #############

LimitofDetect <- my_pipeSummary %>% filter(grepl("new", SampleID)) %>% filter(Pooled_Set == "No") %>% mutate(Ra_cells2 = c(1e2, 1e3, 1e4, 1e5, 1e6))


###########################################################
############ SCATTER: CELL NUMBER VS N_GENOMIC ############

# N_GENOMIC
LimitofDetect_NumReads_Fig1 <- LimitofDetect %>% 
  ggplot(aes(x = Ra_cells2, y = N_Genomic)) + 
  geom_point(size = 6, alpha = 0.8, stroke = 0.8, fill = "#03A9F4", shape = 21) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "THP1 cells spiked with H37Ra",
       subtitle = NULL, 
       x = "# spiked in H37Ra cells", 
       y = "# reads aligning to Mtb genome") + 
  scale_y_continuous(limits = c(0,7000000), breaks = seq(0, 7000000, 1000000)) + 
  scale_x_continuous(trans = "log10") + 
  my_plot_themes
LimitofDetect_NumReads_Fig1
ggsave(LimitofDetect_NumReads_Fig1,
       file = "LimitofDetect_NumReads_1.pdf",
       path = "LimitofDetection",
       width = 6, height = 4, units = "in")

# P_GENOMIC
LimitofDetect_PercentReads_Fig1 <- LimitofDetect %>% 
  ggplot(aes(x = Ra_cells2, y = P_Genomic)) + 
  geom_point(size = 6, alpha = 0.8, stroke = 0.8, fill = "#03A9F4", shape = 21) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  labs(title = "THP1 cells spiked with H37Ra",
       subtitle = "Label is number of reads aligned to Mtb", 
       x = "# spiked in H37Ra cells", 
       y = "% reads aligning to Mtb genome") + 
  scale_x_continuous(trans = "log10") + 
  my_plot_themes
LimitofDetect_PercentReads_Fig1
ggsave(LimitofDetect_PercentReads_Fig1,
       file = "LimitofDetect_PercentReads_1.pdf",
       path = "LimitofDetection",
       width = 6, height = 4, units = "in")


###########################################################
########## COMPARE NUMBER OF GENES (10 OR 100) ############

LimitofDetect_NumGenes1 <- LimitofDetect %>% 
  pivot_longer(cols = c("AtLeast.10.Reads", "AtLeast.100.Reads"), names_to = "AtLeast.X.Reads", values_to = "Num_Genes") %>%
  ggplot(aes(x = Ra_cells2, y = Num_Genes, text = SampleID)) + 
  geom_point(size = 6, alpha = 0.8, stroke = 0.8, fill = "#03A9F4", shape = 21) + 
  facet_grid(~AtLeast.X.Reads, scales = "free") + 
  scale_y_continuous(limits = c(0,4499), breaks = seq(0, 4500, 500)) + 
  labs(title = "THP1 spiked samples # genes aligned to Mtb",
       subtitle = NULL,
       x = "# spiked in H37Ra cells", 
       y = "Number of genes") + 
  scale_x_continuous(trans = "log10") + 
  my_plot_themes
LimitofDetect_NumGenes1
ggsave(LimitofDetect_NumGenes1,
       file = "LimitofDetect_NumGenes1.pdf",
       path = "LimitofDetection",
       width = 6, height = 4, units = "in")

LimitofDetect_10Genes <- LimitofDetect %>% 
  ggplot(aes(x = Ra_cells2, y = AtLeast.10.Reads, text = SampleID)) + 
  geom_point(size = 6, alpha = 0.8, stroke = 0.8, fill = "#03A9F4", shape = 21) + 
  scale_y_continuous(limits = c(0,4499), breaks = seq(0, 4500, 500)) + 
  geom_text_repel(aes(label = format(AtLeast.10.Reads, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  labs(title = "THP1 spiked samples # genes aligned to Mtb with at least 10 TPM",
       subtitle = NULL,
       x = "# spiked in H37Ra cells", 
       y = "# genes with at least 10 reads") + 
  geom_hline(yintercept = 4499/2, linetype = "dashed", alpha = 0.5) + 
  scale_x_continuous(trans = "log10") + 
  my_plot_themes
LimitofDetect_10Genes
ggsave(LimitofDetect_10Genes,
       file = "LimitofDetect_10Genes_v2.pdf",
       path = "LimitofDetection",
       width = 6, height = 4, units = "in")





###########################################################
############ SCATTER: CELL NUMBER VS N_GENOMIC ############

