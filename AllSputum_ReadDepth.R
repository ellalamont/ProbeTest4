# Read depth for all the Sputum samples from the Sept and Nov sequencing runs 
# E. Lamont 
# 12/17/24

source("Import_data.R") # to get AllSputum_pipeSummary

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none", legend.text=element_text(size=10),
        legend.title = element_text(size = 10),
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
################### N_GENOMIC vs WEEK #####################

WeekvsReads_sputum1 <- AllSputum_pipeSummary %>% 
  ggplot(aes(x = Week, y = N_Genomic)) + 
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = c14) +  
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,5500000), breaks = seq(0, 5500000, 1000000)) +
  labs(title = "Sputum: Week vs number reads aligned to Mtb",
       subtitle = "Label is number of reads aligned to Mtb", 
       x = "Weeks after start of antibiotics", 
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes
WeekvsReads_sputum1
ggplotly(WeekvsReads_sputum1)
ggsave(WeekvsReads_sputum1,
       file = "WeekvsReads_sputum1.pdf",
       path = "AllSputum_Figures",
       width = 8, height = 5, units = "in")


###########################################################
################### P_GENOMIC vs WEEK #####################

WeekvsPercent_sputum1 <- AllSputum_pipeSummary %>% 
  ggplot(aes(x = Week, y = P_Genomic)) + 
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = c14) +  
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,75), breaks = seq(0, 75, 10)) +
  labs(title = "Sputum: Week vs percent reads aligned to Mtb",
       subtitle = "Label is number of reads aligned to Mtb", 
       x = "Weeks after start of antibiotics", 
       y = "% reads aligning to Mtb genome") + 
  my_plot_themes
WeekvsPercent_sputum1
ggplotly(WeekvsPercent_sputum1)
ggsave(WeekvsPercent_sputum1,
       file = "WeekvsPercent_sputum1.pdf",
       path = "AllSputum_Figures",
       width = 8, height = 5, units = "in")


###########################################################
################### ALL UNIQUE SPUTUM #####################

# W0 samples: "S_250754_S47", "S_354851_DualrRNA", "S_503557",, "S_503917_DualrRNA"
# W2 samples: "S_349942_DualrRNA_S18", "S_575533_MtbrRNA_S39", "S_349942_DualrRNA", "S_577207_DualrRNA"
# W4 samples: "S_351946_Probe_4A_100", "S_575540_DualrRNA", "S_687338_Probe_4A_100"

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default

# N_GENOMIC vs Week
WeekvsReads_Unique_sputum2 <- AllSputum_pipeSummary %>% 
  filter(SampleID %in% c("S_250754", "S_354851_DualrRNA", "S_503557", "S_503917_DualrRNA", 
                         "S_349942_DualrRNA", "S_575533_MtbrRNA", "S_349942_DualrRNA", "S_577207_DualrRNA", 
                         "S_351946_Probe_4A_100", "S_575540_DualrRNA", "S_687338_Probe_4A_100")) %>% 
  ggplot(aes(x = Week, y = N_Genomic)) + 
  # geom_point(aes(shape = SeqRun), fill = "#0072B2", size = 6, alpha = 0.7, stroke = 0.8, color = "black", shape = 21) + # Used for 1
  geom_point(aes(fill = Week, shape = Week), size = 6, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  # guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,5500000), breaks = seq(0, 5500000, 1000000)) +
  labs(title = "Unique Sputum: Week vs number reads aligned to Mtb",
       subtitle = NULL, 
       x = "Weeks after start of antibiotics", 
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes
WeekvsReads_Unique_sputum2
# ggplotly(WeekvsReads_sputum2)
ggsave(WeekvsReads_Unique_sputum2,
       file = "WeekvsReads_UniqueSputum3.pdf",
       path = "AllSputum_Figures",
       width = 6, height = 4, units = "in")

# P_GENOMIC vs Week
WeekvsPercent_Unique_sputum1 <- AllSputum_pipeSummary %>% 
  filter(SampleID %in% c("S_250754", "S_354851_DualrRNA", "S_503557", "S_503917_DualrRNA", 
                         "S_349942_DualrRNA", "S_575533_MtbrRNA", "S_349942_DualrRNA", "S_577207_DualrRNA", 
                         "S_351946_Probe_4A_100", "S_575540_DualrRNA", "S_687338_Probe_4A_100")) %>% 
  ggplot(aes(x = Week, y = P_Genomic)) + 
  geom_point(aes(shape = SeqRun), fill = "#0072B2", size = 6, alpha = 0.7, stroke = 0.8, color = "black", shape = 21) + 
  # scale_fill_manual(values = c14) +  
  # guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  # scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,75), breaks = seq(0, 75, 10)) +
  labs(title = "Unique Sputum: Week vs percent reads aligned to Mtb",
       subtitle = "Label is number of reads aligned to Mtb", 
       x = "Weeks after start of antibiotics", 
       y = "% reads aligning to Mtb genome") + 
  my_plot_themes
WeekvsPercent_Unique_sputum1
ggplotly(WeekvsPercent_Unique_sputum1)
ggsave(WeekvsPercent_Unique_sputum1,
       file = "WeekvsPercent_UniqueSputum1.pdf",
       path = "AllSputum_Figures",
       width = 7, height = 5, units = "in")




