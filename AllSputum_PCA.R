# PCA plot Sputum from Sept (ProbeTest3) and Nov (ProbeTest4) runs
# 1/8/25


# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

source("Import_data.R") # to get AllSputum_tpm and AllSputum_metadata

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
######################## MAKE PCA #########################

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above


# Think I need to transform the data first
my_tpm_t <- as.data.frame(t(AllSputum_tpm))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 31.0% of variance
summary_PCA[2,1] # PC2 explains 17.9% of variance
summary_PCA[3,1] # PC3 explains 11.1% of variance

###########################################################
################ MAKE PCA PLOT with GGPLOT ################

my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, AllSputum_metadata, by = "SampleID")

fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, fill = Week, shape = SeqRun)) + 
  # ggplot(aes(x = PC1, y = PC2, color = Sample_Type, shape = Strain, text = Replicate)) + 
  geom_point(size = 6, alpha = 0.8) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  # scale_color_manual(values = c(`High_Low_THP1` = "darkorange4", `Sputum` = "#0072B2", `THP1` = "#FF7F00")) + 
  labs(title = "PCA plot AllSputum so far PC1 vs PC2",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
fig_PC1vsPC2
ggplotly(fig_PC1vsPC2)

ggsave(fig_PC1vsPC2,
       file = "PCA AllSputum_PC1vsPC2.pdf",
       path = "PCA_Figures",
       width = 9, height = 6, units = "in")


###########################################################
################### MAKE 3D PCA PLOT ######################

# https://plotly.com/r/pca-visualization/

PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Week# , 
                  # colors = c12,
                  # text = ~Replicate
)
PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")



###########################################################
############## UNIQUE SPUTUM MAKE PCA #####################
# Want a PCA of unique sputum (some samples were sequenced more than once so want to have a PCA of just unique sputum samples so none are repeated)

# W0 samples: "S_250754_S47", "S_354851_DualrRNA_S17", "S_503557_S46"
# W2 samples: "S_349942_DualrRNA_S18", "S_575533_MtbrRNA_S39"


my_tpm_t_UniqueSputum <- my_tpm_t %>% filter(row.names(my_tpm_t2) %in% c("S_250754","S_354851_DualrRNA", "S_503557", "S_349942_DualrRNA", "S_575533_MtbrRNA"))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t_UniqueSputum %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 53.7% of variance
summary_PCA[2,1] # PC2 explains 27.7% of variance
summary_PCA[3,1] # PC3 explains 9.8% of variance

# MAKE PCA PLOT with GGPLOT 

my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, AllSputum_metadata, by = "SampleID")

fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, fill = Week, shape = Week)) + 
  # ggplot(aes(x = PC1, y = PC2, color = Sample_Type, shape = Strain, text = Replicate)) + 
  geom_point(size = 6, alpha = 0.8, stroke = 0.8) +
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  # guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  
  # geom_text_repel(aes(label = Sputum_Number), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  # scale_color_manual(values = c(`High_Low_THP1` = "darkorange4", `Sputum` = "#0072B2", `THP1` = "#FF7F00")) + 
  labs(title = "PCA plot Unique Sputum W0 and W2 only",
       subtitle = "All have >1M reads, all from Nov sequencing run, some dual rRNA depleted",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
fig_PC1vsPC2
ggplotly(fig_PC1vsPC2)

ggsave(fig_PC1vsPC2,
       file = "PCA_UniqueSputum_PC1vsPC2_2.pdf",
       path = "PCA_Figures",
       width = 7, height = 5, units = "in")

