# Compare pooling to Mtb # reads
# E. Lamont
# 12/10/24

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
######### THP1_1e6_1 to 5 NUMBER+PERCENT READS #########
# Numbers 1-4 were pooled, 5 was not, everything else the same

# Number Reads
THP1_1e6_PoolVsReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("THP1_1e6_1", "THP1_1e6_2", "THP1_1e6_3", "THP1_1e6_4", "THP1_1e6_5")) %>% 
  ggplot(aes(x = Pooled_Set, y = N_Genomic)) + 
  geom_point(color = "#FF7F00", size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  # scale_color_manual(values = c(`B` = "#E31A1C", `C` = "green4", `D` = "#6A3D9A", `F` = "maroon", `No` = "black")) + 
  # scale_shape_manual(values=c(1, 16)) +
  # geom_text_repel(aes(label = Library_prep), size= 2) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(500000,7000000), breaks = seq(1000000, 7000000, 1000000)) + 
  labs(title = "F: Pooled for Capture vs number reads aligned to Mtb",
       subtitle = "All original THP1 spiked with 1e6 cells H37Ra") + 
  my_plot_themes
THP1_1e6_PoolVsReads
ggsave(THP1_1e6_PoolVsReads,
       file = "PoolVsReads_THP1_1e6_F.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")

# Percent Reads
HP1_1e6_PoolVsPercentReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("THP1_1e6_1", "THP1_1e6_2", "THP1_1e6_3", "THP1_1e6_4", "THP1_1e6_5")) %>% 
  ggplot(aes(x = Pooled_Set, y = P_Genomic)) + 
  geom_point(color = "#FF7F00", size = 3, alpha = 0.8) + 
  # scale_color_manual(values = c(`B` = "#E31A1C", `C` = "green4", `D` = "#6A3D9A", `F` = "maroon", `No` = "black")) + 
  # scale_shape_manual(values=c(1, 16)) +
  # geom_text_repel(aes(label = Library_prep), size= 2) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(500000,7000000), breaks = seq(1000000, 7000000, 1000000)) + 
  labs(title = "F: Pooled for Capture vs percent reads aligned to Mtb",
       subtitle = "All original THP1 spiked with 1e6 cells H37Ra") + 
  my_plot_themes
HP1_1e6_PoolVsPercentReads
ggsave(HP1_1e6_PoolVsPercentReads,
       file = "PoolVsPercentReads_THP1_1e6_F.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")



###########################################################
####### A: New THP1 varied Mtb NUMBER+PERCENT READS #######

# Number Reads
newTHP1_PoolVsReads <- my_pipeSummary %>% 
  filter(grepl('new', SampleID)) %>% 
  ggplot(aes(x = Pooled_Set, y = N_Genomic)) + 
  geom_line(aes(group = Ra_cells), linewidth = 0.2) + 
  geom_point(aes(color = Ra_cells), size = 3, alpha = 0.9) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  scale_color_manual(values = Orange_Gradient) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "A: Pooled for Capture vs number reads aligned to Mtb",
       subtitle = "newTHP1 spiked with varied cells H37Ra") + 
  my_plot_themes
newTHP1_PoolVsReads
ggsave(newTHP1_PoolVsReads,
       file = "PoolVsReads_newTHP1_A.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")

# Percent Reads
newTHP1_PoolVsPercentReads <- my_pipeSummary %>% 
  filter(grepl('new', SampleID)) %>% 
  ggplot(aes(x = Pooled_Set, y = P_Genomic)) + 
  geom_line(aes(group = Ra_cells), linewidth = 0.2) + 
  geom_point(aes(color = Ra_cells), size = 3, alpha = 0.9) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  scale_color_manual(values = Orange_Gradient) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "A: Pooled for Capture vs percent reads aligned to Mtb",
       subtitle = "newTHP1 spiked with varied cells H37Ra") + 
  my_plot_themes
newTHP1_PoolVsPercentReads
ggsave(newTHP1_PoolVsPercentReads,
       file = "PoolVsPercentReads_newTHP1_A.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")

###########################################################
###### B: original THP1 1e6 Mtb NUMBER+PERCENT READS ######

# Number Reads
B_THP1_PoolVsReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("originalTHP1_1e6_1_DualrRNA", "originalTHP1_1e6_2_DualrRNA", "originalTHP1_1e6_3_DualrRNA")) %>% 
  ggplot(aes(x = Pooled_Set, y = N_Genomic)) + 
  geom_point(color = "#FF7F00", size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  scale_color_manual(values = c7_2) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "B: Pooled for Capture vs number reads aligned to Mtb",
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra; DualrRNADepletion") + 
  my_plot_themes
B_THP1_PoolVsReads
ggsave(B_THP1_PoolVsReads,
       file = "PoolVsReads_B.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")

# Percent Reads
B_THP1_PoolVsPercentReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("originalTHP1_1e6_1_DualrRNA", "originalTHP1_1e6_2_DualrRNA", "originalTHP1_1e6_3_DualrRNA")) %>% 
  ggplot(aes(x = Pooled_Set, y = P_Genomic)) + 
  geom_point(color = "#FF7F00", size = 3, alpha = 0.8) + 
  scale_color_manual(values = c7_2) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "Pooled for Capture vs percent reads aligned to Mtb",
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra; DualrRNADepletion") + 
  my_plot_themes
B_THP1_PoolVsPercentReads
ggsave(B_THP1_PoolVsPercentReads,
       file = "PoolVsPercentReads_B.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")

###########################################################
###### C: original THP1 1e6 Mtb NUMBER+PERCENT READS ######

# Number Reads
C_THP1_PoolVsReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("originalTHP1_1e6_4_DualrRNA", "originalTHP1_1e6_5_DualrRNA", "originalTHP1_1e6_6_DualrRNA")) %>% 
  ggplot(aes(x = Pooled_Set, y = N_Genomic)) + 
  geom_point(color = "#FF7F00", size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  # geom_text(aes(label = format(N_Genomic, big.mark = ",")), size= 3, nudge_x = 0.1, nudge_y = -150000) + 
  scale_color_manual(values = c7_2) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "C: Pooled for Capture vs number reads aligned to Mtb",
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra; DualrRNADepletion; half library prep") + 
  my_plot_themes
C_THP1_PoolVsReads
ggsave(C_THP1_PoolVsReads,
       file = "PoolVsReads_C.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")

# Percent Reads
C_THP1_PoolVsPercentReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("originalTHP1_1e6_4_DualrRNA", "originalTHP1_1e6_5_DualrRNA", "originalTHP1_1e6_6_DualrRNA")) %>% 
  ggplot(aes(x = Pooled_Set, y = P_Genomic)) + 
  geom_point(color = "#FF7F00", size = 3, alpha = 0.8) + 
  scale_color_manual(values = c7_2) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "C: Pooled for Capture vs percent reads aligned to Mtb",
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra; DualrRNADepletion; half library prep") + 
  my_plot_themes
C_THP1_PoolVsPercentReads
ggsave(C_THP1_PoolVsPercentReads,
       file = "PoolVsPercentReads_C.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")


###########################################################
###### D: original THP1 1e6 Mtb NUMBER+PERCENT READS ######

# Number Reads
D_THP1_PoolVsReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("originalTHP1_1e6_1_MtbrRNA", "originalTHP1_1e6_2_MtbrRNA", "originalTHP1_1e6_3_MtbrRNA")) %>% 
  ggplot(aes(x = Pooled_Set, y = N_Genomic)) + 
  geom_point(color = "#FF7F00", size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  scale_color_manual(values = c7_2) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "D: Pooled for Capture vs number reads aligned to Mtb",
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra; MtbrRNADepletion; half library prep") + 
  my_plot_themes
D_THP1_PoolVsReads
ggsave(D_THP1_PoolVsReads,
       file = "PoolVsReads_D.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")

# Percent Reads
D_THP1_PoolVsPercentReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("originalTHP1_1e6_1_MtbrRNA", "originalTHP1_1e6_2_MtbrRNA", "originalTHP1_1e6_3_MtbrRNA")) %>% 
  ggplot(aes(x = Pooled_Set, y = P_Genomic)) + 
  geom_point(color = "#FF7F00", size = 3, alpha = 0.8) + 
  scale_color_manual(values = c7_2) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "D: Pooled for Capture vs percent reads aligned to Mtb",
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra; MtbrRNADepletion; half library prep") + 
  my_plot_themes
D_THP1_PoolVsPercentReads
ggsave(D_THP1_PoolVsPercentReads,
       file = "PoolVsPercentReads_D.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")



###########################################################
############# E: Sputum NUMBER+PERCENT READS ##############

# Number Reads
E_Sputum_PoolVsReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("S_575533_MtbrRNA", "S_687338_MtbrRNA", "S_575533_DualrRNA", "S_687338_DualrRNA")) %>% 
  ggplot(aes(x = Pooled_Set, y = N_Genomic)) + 
  geom_point(aes(shape = Week), color = "#0072B2", size = 3, alpha = 0.8) + 
  geom_line(aes(group = EukrRNADep_Group), linewidth = 0.2) +
  scale_shape_manual(values = c(`0` = 16, `2` = 17, `4` = 15)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "E: Pooled for Capture vs number reads aligned to Mtb",
       subtitle = "Sputum: pooled samples had MtbrRNA depletion, unpooled had dualrRNA depletion") + 
  my_plot_themes
E_Sputum_PoolVsReads
ggsave(E_Sputum_PoolVsReads,
       file = "PoolVsReads_E.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")

# Percent Reads
E_Sputum_PoolVsPercentReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("S_575533_MtbrRNA", "S_687338_MtbrRNA", "S_575533_DualrRNA", "S_687338_DualrRNA")) %>% 
  ggplot(aes(x = Pooled_Set, y = P_Genomic)) + 
  geom_point(aes(shape = Week), color = "#0072B2", size = 3, alpha = 0.8) + 
  geom_line(aes(group = EukrRNADep_Group), linewidth = 0.2) +
  scale_shape_manual(values = c(`0` = 16, `2` = 17, `4` = 15)) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "E: Pooled for Capture vs percent reads aligned to Mtb",
       subtitle = "Sputum: pooled samples had MtbrRNA depletion, unpooled had dualrRNA depletion") + 
  my_plot_themes
E_Sputum_PoolVsPercentReads
ggsave(E_Sputum_PoolVsPercentReads,
       file = "PoolVsPercentReads_E.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")


###########################################################
############# G: Sputum NUMBER+PERCENT READS ##############

# Number Reads
G_Sputum_PoolVsReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("S_503557", "S_503557_DualrRNA")) %>% 
  ggplot(aes(x = Pooled_Set, y = N_Genomic)) + 
  geom_point(aes(shape = Week), color = "#0072B2", size = 3, alpha = 0.8) + 
  geom_line(aes(group = EukrRNADep_Group), linewidth = 0.2) +
  scale_shape_manual(values = c(`0` = 16, `2` = 17, `4` = 15)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "E: Pooled for Capture vs number reads aligned to Mtb",
       subtitle = "Sputum: pooled samples had MtbrRNA depletion, unpooled had dualrRNA depletion") + 
  my_plot_themes
G_Sputum_PoolVsReads
ggsave(G_Sputum_PoolVsReads,
       file = "PoolVsReads_G.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")

# Percent Reads
G_Sputum_PoolVsPercentReads <- my_pipeSummary %>% 
  filter(SampleID %in% c("S_503557", "S_503557_DualrRNA")) %>% 
  ggplot(aes(x = Pooled_Set, y = P_Genomic)) + 
  geom_point(aes(shape = Week), color = "#0072B2", size = 3, alpha = 0.8) + 
  geom_line(aes(group = EukrRNADep_Group), linewidth = 0.2) +
  scale_shape_manual(values = c(`0` = 16, `2` = 17, `4` = 15)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "E: Pooled for Capture vs percent reads aligned to Mtb",
       subtitle = "Sputum: pooled samples had MtbrRNA depletion, unpooled had dualrRNA depletion") + 
  my_plot_themes
G_Sputum_PoolVsPercentReads
ggsave(G_Sputum_PoolVsPercentReads,
       file = "PoolVsPercentReads_G.pdf",
       path = "Pooling_Figures",
       width = 6, height = 4, units = "in")
