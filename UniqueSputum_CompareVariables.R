# Unique Sputum samples from the Sept and Nov sequencing runs 
# Why are some lower read counts than others in the same week? Compare metadata
# E. Lamont 
# 1/20/25

# Unique Sputum: 
# W0 samples: "S_250754_S47", "S_354851_DualrRNA", "S_503557",, "S_503917_DualrRNA"
# W2 samples: "S_349942_DualrRNA_S18", "S_575533_MtbrRNA_S39", "S_349942_DualrRNA", "S_577207_DualrRNA"
# W4 samples: "S_351946_Probe_4A_100", "S_575540_DualrRNA", "S_687338_Probe_4A_100"

# Compare total_RNA_ng, mRNA_ng, ct, ttd
# Using N_Genomic or P_Genomic
# Facet by Week and just look at W0 and W2!

source("Import_data.R") # to get AllSputum_pipeSummary

UniqueSputum_pipeSummary <- AllSputum_pipeSummary %>% 
  filter(SampleID %in% c("S_250754", "S_354851_DualrRNA", "S_503557", "S_503917_DualrRNA", 
                         "S_349942_DualrRNA", "S_575533_MtbrRNA", "S_349942_DualrRNA", "S_577207_DualrRNA", 
                         "S_351946_Probe_4A_100", "S_575540_DualrRNA", "S_687338_Probe_4A_100"))

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right", legend.text=element_text(size=14),
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


my_regression_line <- stat_poly_line(method = "lm", se = TRUE, level = 0.95, color = "black", alpha = 0.3)
my_regression_equations <- stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                                          after_stat(rr.label),
                                                          after_stat(p.value.label), 
                                                          sep = "*\", \"*")))

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default

###########################################################
################### GENOMIC vs Ct #######################
# Week 2 samples do not have any Ct values!

# N_GENOMIC
ctVsReads_sputum <- UniqueSputum_pipeSummary %>% 
  filter(Week != "2") %>% # Week 2 does not have Ct values
  ggplot(aes(x = ct, y = N_Genomic)) +
  geom_point(aes(fill = Week, shape = Week), size = 6, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(-3000000,6000000), breaks = seq(-3000000,6000000, 1000000)) +
  scale_x_continuous(limits = c(16,41), breaks = seq(16,41,2), expand = c(0,0)) +
  labs(title = "Unique Sputum: Ct value vs number reads aligned to Mtb",
       subtitle = NULL,
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = TRUE, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ctVsReads_sputum # + my_regression_line + my_regression_equations
ggsave(ctVsReads_sputum,
       file = "UniqueSputum_ctVsReads_v2.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 7, height = 5, units = "in")

# P_GENOMIC
ctVsPercent_sputum <- UniqueSputum_pipeSummary %>% 
  filter(Week != "2") %>% # Week 2 does not have Ct values
  ggplot(aes(x = ct, y = P_Genomic)) +
  geom_point(aes(fill = Week, shape = Week), size = 6, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,5500000), breaks = seq(0, 5500000, 1000000)) +
  labs(title = "Unique Sputum: Ct value vs percent reads aligned to Mtb",
       subtitle = "Week 2 samples did not have Ct values",
       y = "% reads aligning to Mtb genome") + 
  my_plot_themes
ctVsPercent_sputum + my_regression_line + my_regression_equations
ggsave(ctVsPercent_sputum,
       file = "UniqueSputum_ctVsPercent.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")


###########################################################
################### GENOMIC vs ttd ########################

# N_GENOMIC
ttdVsReads_sputum <- UniqueSputum_pipeSummary %>% 
  ggplot(aes(x = ttd, y = N_Genomic)) +
  geom_point(aes(fill = Week, shape = Week), size = 6, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(-3000000,6000000), breaks = seq(-3000000,6000000, 1000000)) +
  scale_x_continuous(limits = c(0,23), breaks = seq(0,23,2), expand = c(0, 0)) +
  labs(title = "Unique Sputum: TTD value vs number reads aligned to Mtb",
       subtitle = NULL,
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = TRUE, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ttdVsReads_sputum
ggsave(ttdVsReads_sputum,
       file = "UniqueSputum_ttdVsReads_v2.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 7, height = 5, units = "in")




###########################################################
################ GENOMIC vs Total RNA #####################

# N_GENOMIC
TotalRNAVsReads_sputum <- UniqueSputum_pipeSummary %>% 
  ggplot(aes(x = total_RNA_ng, y = N_Genomic)) +
  geom_point(aes(fill = Week, shape = Week), size = 6, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,5500000), breaks = seq(0, 5500000, 1000000)) +
  labs(title = "Unique Sputum: Total RNA (ng) vs number reads aligned to Mtb",
       subtitle = NULL, 
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes
TotalRNAVsReads_sputum
ggsave(ctVsReads_sputum,
       file = "UniqueSputum_ctVsReads.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")


###########################################################
################### GENOMIC vs mRNA #######################

# N_GENOMIC
mRNAVsReads_sputum <- UniqueSputum_pipeSummary %>% 
  ggplot(aes(x = mRNA_ng, y = N_Genomic)) +
  geom_point(aes(fill = Week, shape = Week), size = 6, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,5500000), breaks = seq(0, 5500000, 1000000)) +
  labs(title = "Unique Sputum: mRNA (ng) vs number reads aligned to Mtb",
       subtitle = NULL, 
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes
mRNAVsReads_sputum
ggsave(mRNAVsReads_sputum,
       file = "UniqueSputum_mRNAVsReads.pdf",
       path = "Sputum_ReadsVsVariables_Figures",
       width = 6, height = 4, units = "in")





















