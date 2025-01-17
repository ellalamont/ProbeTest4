# Compare eukaryotic rRNA depletion to # of genes with at least 10 or 100 reads
# E. Lamont
# 12/11/24

source("Import_data.R") # to get my_pipeSummary
my_pipeSummary <- my_pipeSummary %>% mutate(Sputum_Number = str_extract(SampleID, "S_[0-9]+"))

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
################### NUMBER GENES THP1 #####################
# Looking at all the 1e6 THP1s, should be technical replicates of each other

THP1_1e6_DepVsGenes_scatter <- my_pipeSummary %>% 
  filter(Sample_Type == "THP1") %>% 
  filter(Ra_cells == "one_e_6") %>%
  filter(THP1_Set == "Original") %>% # To avoid anything going on with the 2 new ones
  pivot_longer(cols = c("AtLeast.10.Reads", "AtLeast.100.Reads"), names_to = "AtLeast.X.Reads", values_to = "Num_Genes") %>%
  ggplot(aes(x = EukrRNADep, y = Num_Genes, text = SampleID)) + 
  geom_point(aes(color = Pooled_Set, shape = Hyb_Restarted), size = 3, alpha = 0.8) + 
  facet_grid(~AtLeast.X.Reads, scales = "free") + 
  # scale_color_manual(values = c(`B` = "#E31A1C", `C` = "green4", `D` = "#6A3D9A", `F` = "maroon", `No` = "black")) + 
  # scale_shape_manual(values=c(1, 16)) +
  # geom_text_repel(aes(label = Library_prep), size= 2) + 
  scale_y_continuous(limits = c(3000,4500), breaks = seq(3000, 4500, 500)) + 
  labs(title = "Eukaryotic rRNA depletion vs # genes aligned to Mtb",
       subtitle = "All original THP1 spiked with 1e6 cells H37Ra") + 
  my_plot_themes
THP1_1e6_DepVsGenes_scatter
ggplotly(THP1_1e6_DepVsGenes_scatter)
ggsave(THP1_1e6_DepVsGenes_scatter,
       file = "THP1_1e6_DepVsGenes_scatter.pdf",
       path = "EukrRNADepletion_Figures",
       width = 6, height = 4, units = "in")


###########################################################
################ NUMBER GENES SPUTUM ALL ##################

SputumAll_DepVsGenes_scatter <- my_pipeSummary %>% 
  filter(Sample_Type == "Sputum") %>% 
  pivot_longer(cols = c("AtLeast.10.Reads", "AtLeast.100.Reads"), names_to = "AtLeast.X.Reads", values_to = "Num_Genes") %>%
  ggplot(aes(x = EukrRNADep, y = Num_Genes, text = SampleID)) + 
  geom_point(aes(shape = Week, fill = Sputum_Number), color = "black", size = 3, alpha = 0.8) + 
  facet_grid(~AtLeast.X.Reads, scales = "free") + 
  scale_fill_manual(values = c14) +  
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values = c(`0` = 21, `2` = 24, `4` = 22)) + 
  # geom_text_repel(aes(label = Library_prep), size= 2) + 
  # scale_y_continuous(limits = c(3000,4500), breaks = seq(3000, 4500, 500)) + 
  labs(title = "Eukaryotic rRNA depletion vs # genes aligned to Mtb",
       subtitle = "All Sputum") + 
  my_plot_themes
SputumAll_DepVsGenes_scatter

ggsave(SputumAll_DepVsGenes_scatter,
       file = "SputumAll_DepVsGenes_scatter.pdf",
       path = "EukrRNADepletion_Figures",
       width = 6, height = 4, units = "in")
