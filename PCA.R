# PCA plot
# 12/10/24

# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

source("Import_data.R") # to get my_tpm

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
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
######################## MAKE PCA #########################

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above


# Think I need to transform the data first
my_tpm_t <- as.data.frame(t(my_tpm))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 16.7% of variance
summary_PCA[2,1] # PC2 explains 11.0% of variance
summary_PCA[3,1] # PC3 explains 9.9% of variance

###########################################################
################ MAKE PCA PLOT with GGPLOT ################

my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, my_metadata, by = "SampleID")

fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, color = Sample_Type, shape = EukrRNADep, label = SampleID, label2 = Week)) + 
  # ggplot(aes(x = PC1, y = PC2, color = Sample_Type, shape = Strain, text = Replicate)) + 
  geom_point(size = 3) +
  geom_text_repel(aes(label = Week), color = "black", size = 2) + 
  scale_color_manual(values = c(`High_Low_THP1` = "darkorange4", `Sputum` = "#0072B2", `THP1` = "#FF7F00")) + 
  labs(title = "PCA plot ProbeTest4: PC1 vs PC2",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
fig_PC1vsPC2
ggplotly(fig_PC1vsPC2)

ggsave(fig_PC1vsPC2,
       file = "PCA_all_PC1vsPC2.pdf",
       path = "PCA_Figures",
       width = 9, height = 6, units = "in")

# NOT DONE
# fig_PC1vsPC3 <- my_PCA_df %>%
#   ggplot(aes(x = PC1, y = PC3, color = Sample_Type, label = SampleID, label2 = Week)) + 
#   geom_point(size = 3) +
#   geom_text_repel(aes(label = Week), color = "black", size = 2) + 
#   scale_color_manual(values = c(`Marmoset` = "#CAB2D6", `Sputum` = "#0072B2", `Saliva` = "#009E73", `THP1` = "#FF7F00")) + 
#   labs(title = "PCA plot ProbeTest3: PC1 vs PC3",
#        x = paste0("PC1: ", summary_PCA[1,1], "%"),
#        y = paste0("PC3: ", summary_PCA[3,1], "%")) +
#   my_plot_themes
# fig_PC1vsPC3
# ggplotly(fig_PC1vsPC3)
# 
# ggsave(fig_PC1vsPC3,
#        file = "PCA_PC1vsPC3.pdf",
#        path = "PCA_Figures",
#        width = 9, height = 6, units = "in")


###########################################################
################### MAKE 3D PCA PLOT ######################

# https://plotly.com/r/pca-visualization/

PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                    type = "scatter3d", mode = "markers",
                    color = ~Sample_Type# , 
                  # colors = c12,
                    # text = ~Replicate
                  )
PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")


###########################################################
################### PCA JUST THP1 #########################

my_PCA_THP1 <- prcomp(my_tpm_t2[grep("THP1_", row.names(my_tpm_t2)),], scale = TRUE)
summary_PCA_THP1 <- format(round(as.data.frame(summary(my_PCA_THP1)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1)

my_PCA_THP1_df <- as.data.frame(my_PCA_THP1$x[, 1:3]) # Extract the first 3 PCs
my_PCA_THP1_df <- data.frame(SampleID = row.names(my_PCA_THP1_df), my_PCA_THP1_df)
my_PCA_THP1_df <- merge(my_PCA_THP1_df, my_metadata, by = "SampleID")

# Add some columns to make the shapes work
# my_PCA_THP1_df$Probe_ng_chr <- as.character(my_PCA_THP1_df$Probe_ng)
# orderd_Probe_ng_chr <- c("0", "10", "15.4", "16.4", "25", "50", "100")
# my_PCA_THP1_df$Probe_ng_chr <- factor(my_PCA_THP1_df$Probe_ng_chr, levels = orderd_Probe_ng_chr)

fig_PC1vsPC2_THP1 <- my_PCA_THP1_df %>%
  ggplot(aes(x = PC1, y = PC2, shape = EukrRNADep, label = Ra_cells)) + 
  geom_point(size = 5, color = "#FF7F00", alpha = 0.8) +
  geom_text_repel(aes(label = Ra_cells), color = "black", size = 3, box.padding = 0.4) + 
  # scale_shape_manual(values=c(1, 16), labels = c("1e6", "1e8")) +
  labs(title = "THP1 PCA plot ProbeTest4: PC1 vs PC2",
       x = paste0("PC1: ", summary_PCA_THP1[1,1], "%"),
       y = paste0("PC2: ", summary_PCA_THP1[2,1], "%"),
       shape = "rRNA Depletion") +
  my_plot_themes
fig_PC1vsPC2_THP1
# ggplotly(fig_PC1vsPC2_THP1)

ggsave(fig_PC1vsPC2_THP1,
       file = "THP1_PCA_PC1vsPC2.pdf",
       path = "PCA_Figures",
       width = 6, height = 4, units = "in")


###########################################################
################### PCA JUST SPUTUM #######################

my_tpm_t2_Sputum <- my_tpm_t2[grep("S_", row.names(my_tpm_t2)),] # Had to make this seperately because need to remove some columns that are all zero here
# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2_Sputum <- my_tpm_t2_Sputum %>% select_if(colSums(.) != 0)

my_PCA_Sputum <- prcomp(my_tpm_t2_Sputum, scale = TRUE)
summary_PCA_Sputum <- format(round(as.data.frame(summary(my_PCA_Sputum)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # Adding round and digits here to set number of digits after the decimal place.

my_PCA_Sputum_df <- as.data.frame(my_PCA_Sputum$x[, 1:3]) # Extract the first 3 PCs
my_PCA_Sputum_df <- data.frame(SampleID = row.names(my_PCA_Sputum_df), my_PCA_Sputum_df)
my_PCA_Sputum_df <- merge(my_PCA_Sputum_df, my_metadata, by = "SampleID")

fig_PC1vsPC2_Sputum <- my_PCA_Sputum_df %>%
  ggplot(aes(x = PC1, y = PC2, shape = EukrRNADep, label = Week)) + 
  geom_point(aes(shape = EukrRNADep), alpha = 0.8, stroke = 0.8, size = 6, color = "#0072B2") +
  geom_text_repel(aes(label = Week), color = "black", size = 3, box.padding = 0.4) + 
  # scale_shape_manual("Sputum \ncollection week", values = c(21,22,24)) +
  # scale_shape_manual(values=c(15, 0, 3)) +
  labs(title = "Sputum PCA plot ProbeTest4: PC1 vs PC2",
       x = paste0("PC1: ", summary_PCA_Sputum[1,1], "%"),
       y = paste0("PC2: ", summary_PCA_Sputum[2,1], "%")) +
  my_plot_themes
fig_PC1vsPC2_Sputum
ggplotly(fig_PC1vsPC2_Sputum)

ggsave(fig_PC1vsPC2_Sputum,
       file = "Sputum_PCA_PC1vsPC2.pdf",
       path = "PCA_Figures",
       width = 6, height = 4, units = "in")

###########################################################
############### PCA JUST PADMINI'S STRAINS ################

my_tpm_t2_HighLowTHP1 <- my_tpm_t2[grep("_T", row.names(my_tpm_t2)),] # Had to make this seperately because need to remove some columns that are all zero here
# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2_HighLowTHP1 <- my_tpm_t2_HighLowTHP1 %>% select_if(colSums(.) != 0)

my_PCA_HighLowTHP1 <- prcomp(my_tpm_t2_HighLowTHP1, scale = TRUE)
summary_PCA_HighLowTHP1 <- format(round(as.data.frame(summary(my_PCA_HighLowTHP1)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # Adding round and digits here to set number of digits after the decimal place.

my_PCA_HighLowTHP1_df <- as.data.frame(my_PCA_HighLowTHP1$x[, 1:3]) # Extract the first 3 PCs
my_PCA_HighLowTHP1_df <- data.frame(SampleID = row.names(my_PCA_HighLowTHP1_df), my_PCA_HighLowTHP1_df)
my_PCA_HighLowTHP1_df <- merge(my_PCA_HighLowTHP1_df, my_metadata, by = "SampleID")

fig_PC1vsPC2_HighLowTHP1 <- my_PCA_HighLowTHP1_df %>%
  ggplot(aes(x = PC1, y = PC2, label = SampleID)) + 
  geom_point(alpha = 0.8, stroke = 0.8, size = 6, color = "darkorange4") +
  geom_text_repel(aes(label = SampleID), color = "black", size = 3, box.padding = 0.4) + 
  labs(title = "High/Low Transmission THP1 PCA plot ProbeTest4: PC1 vs PC2",
       x = paste0("PC1: ", summary_PCA_Sputum[1,1], "%"),
       y = paste0("PC2: ", summary_PCA_Sputum[2,1], "%")) +
  my_plot_themes
fig_PC1vsPC2_HighLowTHP1
ggplotly(fig_PC1vsPC2_HighLowTHP1)

ggsave(fig_PC1vsPC2_HighLowTHP1,
       file = "HighLowTHP1_PCA_PC1vsPC2.pdf",
       path = "PCA_Figures",
       width = 6, height = 4, units = "in")
