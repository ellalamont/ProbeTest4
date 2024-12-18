# Compare eukaryotic rRNA depletion to Mtb # reads
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
############### SCATTER NUMBER READS THP1 #################
# Looking at all the 1e6 THP1s, should be technical replicates of each other

THP1_1e6_DepVsReads_scatter <- my_pipeSummary %>% 
  filter(Sample_Type == "THP1") %>% 
  filter(Ra_cells == "one_e_6") %>%
  filter(THP1_Set == "Original") %>% # To avoid anything going on with the 2 new ones
  ggplot(aes(x = EukrRNADep, y = N_Genomic)) + 
  geom_point(aes(color = Pooled_Set, shape = Hyb_Restarted), size = 3, alpha = 0.8) + 
  # scale_color_manual(values = c(`B` = "#E31A1C", `C` = "green4", `D` = "#6A3D9A", `F` = "maroon", `No` = "black")) + 
  # scale_shape_manual(values=c(1, 16)) +
  # geom_text_repel(aes(label = Library_prep), size= 2) + 

  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(500000,7000000), breaks = seq(1000000, 7000000, 1000000)) + 
  
  labs(title = "Eukaryotic rRNA depletion vs number reads aligned to Mtb",
       subtitle = "All original THP1 spiked with 1e6 cells H37Ra") + 
  
  my_plot_themes

THP1_1e6_DepVsReads_scatter

ggsave(THP1_1e6_DepVsReads_scatter,
       file = "EukrRNADepVsReads_THP1_1e6_scatter.pdf",
       path = "EukrRNADepletion_Figures",
       width = 6, height = 4, units = "in")


###########################################################
############## SCATTER PERCENT READS THP1 #################
# Looking at all the 1e6 THP1s, should be technical replicates of each other

THP1_1e6_DepVsPercentReads_scatter <- my_pipeSummary %>% 
  filter(Sample_Type == "THP1") %>% 
  filter(Ra_cells == "one_e_6") %>%
  filter(THP1_Set == "Original") %>% # To avoid anything going on with the 2 new ones
  ggplot(aes(x = EukrRNADep, y = P_Genomic)) + 
  geom_point(aes(color = Pooled_Set, shape = Hyb_Restarted), size = 3, alpha = 0.8) + 
  # scale_color_manual(values = c(`High_Low_THP1` = "darkorange4", `Sputum` = "#0072B2", `THP1` = "#FF7F00")) + 
  # scale_shape_manual(values=c(1, 16)) +
  # geom_text_repel(aes(label = Library_prep), size= 2) + 
  
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  
  labs(title = "Eukaryotic rRNA depletion vs percent reads aligned to Mtb",
       subtitle = "All original THP1 spiked with 1e6 cells H37Ra") + 
  
  my_plot_themes

THP1_1e6_DepVsPercentReads_scatter

ggsave(THP1_1e6_DepVsPercentReads_scatter,
       file = "EukrRNADepVsPercentReads_THP1_1e6_scatter.pdf",
       path = "EukrRNADepletion_Figures",
       width = 6, height = 4, units = "in")
  
  
###########################################################
############## SCATTER NUMBER READS SPUTUM ################
# There are a couple matching, but most not

Sputum_DepVsReads_scatter <- my_pipeSummary %>% 
  filter(Sample_Type == "Sputum") %>% 
  
  ggplot(aes(x = EukrRNADep, y = N_Genomic)) + 
  geom_point(aes(color = Pooled_Set, shape = Week), size = 3, alpha = 0.8) + 
  # geom_line(linewidth = 0.2) +
  # scale_color_manual(values = c(`High_Low_THP1` = "darkorange4", `Sputum` = "#0072B2", `THP1` = "#FF7F00")) + 
  # scale_shape_manual(values=c(1, 16)) +
  # geom_text_repel(aes(label = Library_prep), size= 2) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(500000,7000000), breaks = seq(1000000, 7000000, 1000000)) + 
  labs(title = "Eukaryotic rRNA depletion vs number reads aligned to Mtb",
       subtitle = "All Sputum") + 
  my_plot_themes

Sputum_DepVsReads_scatter

ggsave(Sputum_DepVsReads_scatter,
       file = "EukrRNADepVsReads_SputumAll_scatter.pdf",
       path = "EukrRNADepletion_Figures",
       width = 6, height = 4, units = "in")

###########################################################
############# SCATTER PERCENT READS SPUTUM ################
# There are a couple matching, but most not

Sputum_DepVsPercentReads_scatter <- my_pipeSummary %>% 
  filter(Sample_Type == "Sputum") %>% 
  
  ggplot(aes(x = EukrRNADep, y = P_Genomic)) + 
  geom_point(aes(color = Pooled_Set, shape = Week), size = 3, alpha = 0.8) + 
  # geom_line(linewidth = 0.2) +
  # scale_color_manual(values = c(`High_Low_THP1` = "darkorange4", `Sputum` = "#0072B2", `THP1` = "#FF7F00")) + 
  # scale_shape_manual(values=c(1, 16)) +
  # geom_text_repel(aes(label = Library_prep), size= 2) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(500000,7000000), breaks = seq(1000000, 7000000, 1000000)) + 
  labs(title = "Eukaryotic rRNA depletion vs number reads aligned to Mtb",
       subtitle = "All Sputum") + 
  my_plot_themes

Sputum_DepVsPercentReads_scatter

ggsave(Sputum_DepVsPercentReads_scatter,
       file = "EukrRNADepVsPercentReads_SputumAll_scatter.pdf",
       path = "EukrRNADepletion_Figures",
       width = 6, height = 4, units = "in")

###########################################################
########### SCATTER NUMBER READS SPUTUM SUBSET ############
# Just looking at the matching samples

SputumSubset_DepVsReads_scatter <- my_pipeSummary %>% 
  filter(Sample_Type == "Sputum") %>% 
  # filter(EukrRNADep_Group != "NA") %>%
  ggplot(aes(x = EukrRNADep, y = N_Genomic)) + 
  geom_point(aes(fill = Sputum_Number, shape = Week), size = 3, alpha = 0.8) + 
  geom_line(aes(group = Sputum_Number), linewidth = 0.2) +
  scale_shape_manual(values = c(`0` = 21, `2` = 24, `4` = 22)) + 
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  # scale_color_manual(values = c(`High_Low_THP1` = "darkorange4", `Sputum` = "#0072B2", `THP1` = "#FF7F00")) + 
  # geom_text_repel(aes(label = Sputum_Number), size= 2) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(500000,7000000), breaks = seq(1000000, 7000000, 1000000)) + 
  labs(title = "Eukaryotic rRNA depletion vs number reads aligned to Mtb",
       subtitle = "Subset of matching sputum") + 
  my_plot_themes

SputumSubset_DepVsReads_scatter

ggsave(SputumSubset_DepVsReads_scatter,
       file = "EukrRNADepVsReads_SputumSubset_scatter.pdf",
       path = "EukrRNADepletion_Figures",
       width = 6, height = 4, units = "in")

###########################################################
########## SCATTER PERCENT READS SPUTUM SUBSET ############
# Just looking at the matching samples

SputumSubset_DepVsPercentReads_scatter <- my_pipeSummary %>% 
  filter(Sample_Type == "Sputum") %>% 
  # filter(EukrRNADep_Group != "NA") %>%
  ggplot(aes(x = EukrRNADep, y = P_Genomic)) + 
  geom_point(aes(fill = Sputum_Number, shape = Week), size = 3, alpha = 0.8) + 
  geom_line(aes(group = Sputum_Number), linewidth = 0.2) +
  scale_shape_manual(values = c(`0` = 21, `2` = 24, `4` = 22)) + 
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  # scale_color_manual(values = c(`High_Low_THP1` = "darkorange4", `Sputum` = "#0072B2", `THP1` = "#FF7F00")) + 
  # geom_text_repel(aes(label = Sputum_Number), size= 2) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(500000,7000000), breaks = seq(1000000, 7000000, 1000000)) + 
  labs(title = "Eukaryotic rRNA depletion vs percent reads aligned to Mtb",
       subtitle = "Subset of matching sputum") + 
  my_plot_themes

SputumSubset_DepVsPercentReads_scatter

ggsave(SputumSubset_DepVsPercentReads_scatter,
       file = "EukrRNADepVsPercentReads_SputumSubset_scatter.pdf",
       path = "EukrRNADepletion_Figures",
       width = 6, height = 4, units = "in")
