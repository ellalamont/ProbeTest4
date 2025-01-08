# All the Sputum samples from the Sept and Nov sequencing runs 
# Why are some lower read counts than others in the same week? Compare metadata
# E. Lamont 
# 1/8/25

# Compare total_RNA_ng, mRNA_ng, ct, ttd, Hyb_Time, Hyb_Restarted (maybe), 
# Using N_Genomic or P_Genomic
# Facet by Week and just look at W0 and W2!

source("Import_data.R") # to get AllSputum_pipeSummary

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


###########################################################
################### GENOMIC vs Ct #######################
# Week 2 samples do not have any Ct values!

# N_GENOMIC
ctVsReads_sputum <- AllSputum_pipeSummary %>% 
  filter(Week == "0") %>% # Don't want Week 4 and Week 2 does not have Ct values
  ggplot(aes(x = ct, y = N_Genomic)) +
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = setNames(AllSputum_pipeSummary$Colours, AllSputum_pipeSummary$Sputum_Number)) +  # Map Sample to Color
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) +
  facet_grid(~ Week, scales = "free") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,4000000), breaks = seq(0, 4000000, 500000)) + 
  labs(title = "All Sputum: Ct value vs number reads aligned to Mtb",
       subtitle = "All Week 0 (Week 2 samples did not have Ct values") + 
  my_plot_themes
ctVsReads_sputum
ggsave(ctVsReads_sputum,
       file = "ctVsReads_sputum.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")

# P_GENOMIC
ctVsPercent_sputum <- AllSputum_pipeSummary %>% 
  filter(Week == "0") %>% # Don't want Week 4 and Week 2 does not have Ct values
  ggplot(aes(x = ct, y = P_Genomic)) +
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = setNames(AllSputum_pipeSummary$Colours, AllSputum_pipeSummary$Sputum_Number)) +  # Map Sample to Color
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) +
  facet_grid(~ Week, scales = "free") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,4000000), breaks = seq(0, 4000000, 500000)) + 
  labs(title = "All Sputum: Ct value vs Percent reads aligned to Mtb",
       subtitle = "All Week 0 (Week 2 samples did not have Ct values") + 
  my_plot_themes
ctVsPercent_sputum
ggsave(ctVsPercent_sputum,
       file = "ctVsPercent_sputum.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")


###########################################################
################### GENOMIC vs TTD ########################

# N_GENOMIC
ttdVsReads_sputum <- AllSputum_pipeSummary %>% 
  filter(Week != "4") %>% # Don't want Week 4 
  ggplot(aes(x = ttd, y = N_Genomic)) +
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = setNames(AllSputum_pipeSummary$Colours, AllSputum_pipeSummary$Sputum_Number)) +  # Map Sample to Color
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) +
  facet_grid(~ Week, scales = "free") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,4000000), breaks = seq(0, 4000000, 500000)) + 
  labs(title = "All Sputum: TTD value vs number reads aligned to Mtb",
       subtitle = NULL) + 
  my_plot_themes
ttdVsReads_sputum
ggsave(ttdVsReads_sputum,
       file = "ttdVsReads_sputum.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")

# P_GENOMIC
ttdVsPercent_sputum <- AllSputum_pipeSummary %>% 
  filter(Week != "4") %>% # Don't want Week 4 
  ggplot(aes(x = ttd, y = P_Genomic)) +
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = setNames(AllSputum_pipeSummary$Colours, AllSputum_pipeSummary$Sputum_Number)) +  # Map Sample to Color
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) +
  facet_grid(~ Week, scales = "free") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,4000000), breaks = seq(0, 4000000, 500000)) + 
  labs(title = "All Sputum: TTD value vs Percent reads aligned to Mtb",
       subtitle = NULL) + 
  my_plot_themes
ttdVsPercent_sputum
ggsave(ttdVsPercent_sputum,
       file = "ttdVsPercent_sputum.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")

###########################################################
################ GENOMIC vs Total RNA #####################

# N_GENOMIC
TotalRNAVsReads_sputum <- AllSputum_pipeSummary %>% 
  filter(Week != "4") %>% # Don't want Week 4 
  ggplot(aes(x = total_RNA_ng, y = N_Genomic)) +
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = setNames(AllSputum_pipeSummary$Colours, AllSputum_pipeSummary$Sputum_Number)) +  # Map Sample to Color
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) +
  facet_grid(~ Week, scales = "free") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,4000000), breaks = seq(0, 4000000, 500000)) + 
  labs(title = "All Sputum: Total RNA (ng) vs number reads aligned to Mtb",
       subtitle = NULL) + 
  my_plot_themes
TotalRNAVsReads_sputum
ggsave(TotalRNAVsReads_sputum,
       file = "TotalRNAVsReads_sputum.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")

# P_GENOMIC
TotalRNAVsPercent_sputum <- AllSputum_pipeSummary %>% 
  filter(Week != "4") %>% # Don't want Week 4 
  ggplot(aes(x = total_RNA_ng, y = P_Genomic)) +
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = setNames(AllSputum_pipeSummary$Colours, AllSputum_pipeSummary$Sputum_Number)) +  # Map Sample to Color
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) +
  facet_grid(~ Week, scales = "free") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,4000000), breaks = seq(0, 4000000, 500000)) + 
  labs(title = "All Sputum: Total RNA (ng) vs Percent reads aligned to Mtb",
       subtitle = NULL) + 
  my_plot_themes
TotalRNAVsPercent_sputum
ggsave(TotalRNAVsPercent_sputum,
       file = "TotalRNAVsPercent_sputum.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")


###########################################################
################### GENOMIC vs mRNA #######################

# N_GENOMIC
mRNAVsReads_sputum <- AllSputum_pipeSummary %>% 
  filter(Week != "4") %>% # Don't want Week 4 
  ggplot(aes(x = mRNA_ng, y = N_Genomic)) +
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = setNames(AllSputum_pipeSummary$Colours, AllSputum_pipeSummary$Sputum_Number)) +  # Map Sample to Color
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) +
  facet_grid(~ Week, scales = "free") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,4000000), breaks = seq(0, 4000000, 500000)) + 
  labs(title = "All Sputum: mRNA (ng) vs number reads aligned to Mtb",
       subtitle = NULL) + 
  my_plot_themes
mRNAVsReads_sputum
ggsave(mRNAVsReads_sputum,
       file = "mRNAVsReads_sputum.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")

# P_GENOMIC
mRNAVsPercent_sputum <- AllSputum_pipeSummary %>% 
  filter(Week != "4") %>% # Don't want Week 4 
  ggplot(aes(x = mRNA_ng, y = P_Genomic)) +
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = setNames(AllSputum_pipeSummary$Colours, AllSputum_pipeSummary$Sputum_Number)) +  # Map Sample to Color
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) +
  facet_grid(~ Week, scales = "free") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,4000000), breaks = seq(0, 4000000, 500000)) + 
  labs(title = "All Sputum: mRNA (ng) vs Percent reads aligned to Mtb",
       subtitle = NULL) + 
  my_plot_themes
mRNAVsPercent_sputum
ggsave(mRNAVsPercent_sputum,
       file = "mRNAVsPercent_sputum.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")


###########################################################
############## GENOMIC vs Hyb Restarted ###################

# N_GENOMIC
HybRetstartedVsReads_sputum <- AllSputum_pipeSummary %>% 
  filter(Week != "4") %>% # Don't want Week 4 
  ggplot(aes(x = Hyb_Restarted, y = N_Genomic)) +
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = setNames(AllSputum_pipeSummary$Colours, AllSputum_pipeSummary$Sputum_Number)) +  # Map Sample to Color
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) +
  facet_grid(~ Week, scales = "free") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,4000000), breaks = seq(0, 4000000, 500000)) + 
  labs(title = "All Sputum: Hyb restarted vs number reads aligned to Mtb",
       subtitle = NULL) + 
  my_plot_themes
HybRetstartedVsReads_sputum
ggsave(HybRetstartedVsReads_sputum,
       file = "HybRetstartedVsReads_sputum.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")
