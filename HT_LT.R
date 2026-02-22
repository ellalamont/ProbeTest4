# Padmini's HT/LT strains (grown in infected macrophages)
# Don't really trust them because they were prepared with half library prep and half blockers!
# Also the hybridization was with the power cut so it had to be restarted!
# E. Lamont
# 2/12/26


source("Import_data.R")


###########################################################
###################### ORGANIZE DATA ###################### 

# pipeSummary
HTLT_pipeSummary <- my_pipeSummary %>% 
  filter(Sample_Type == "High_Low_THP1") %>% 
  select(-X) %>%
  mutate(Group = case_when(
    str_detect(SampleID, "HI1") ~ "HI1",
    str_detect(SampleID, "LT1") ~ "LT1",
    str_detect(SampleID, "UI") ~ "UI",
    TRUE ~ NA))

# RawReads
ProbTest4_RawReads <- read.csv("Mtb.Expression.Gene.Data.readsM.csv")
HTLT_RawReads <- ProbTest4_RawReads %>% select(X, contains("Thp", ignore.case = FALSE))
names(HTLT_RawReads) <- gsub(x = names(HTLT_RawReads), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it
HTLT_RawReads_f <- HTLT_RawReads %>%
  filter(grepl("^Rv[0-9]+[A-Za-z]?$", X))

# TPM
source("Function_CalculateTPM.R")
HTLT_tpm_f <- CalculateTPM_RvOnly(HTLT_RawReads_f)

stopifnot(all(colnames(HTLT_tpm_f) == HTLT_pipeSummary$SampleID))
stopifnot(all(HTLT_pipeSummary$SampleID == colnames(HTLT_tpm_f)))


###########################################################
########################### PCA ########################### 

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=12),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=8))

# Labelled Colors
my_fav_colors <- c(`HI1` = "red4", `LT1` = "aquamarine2")
# Labelled Shapes
my_fav_shapes <- c(`HI1` = 21, `LT1` = 23)

my_tpm <- HTLT_tpm_f %>% # column_to_rownames(var = "X")
  select(!contains("UI")) # Remove the untreated samples

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm)) # or my_tpm2

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 41.3% of variance
summary_PCA[2,1] # PC2 explains 22.4% of variance
summary_PCA[3,1] # PC3 explains 16.5% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, HTLT_pipeSummary, by = "SampleID", )

PCA_fig <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Group, shape = Group)) + 
  geom_point(aes(fill = Group, shape = Group), size = 5, alpha = 0.8, stroke = 0.8) +
  # geom_text_repel(aes(label = Lineage), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors) +  
  scale_shape_manual(values = my_fav_shapes) + 
  geom_text_repel(aes(label = SampleID), size= 2, box.padding = 0.4, segment.color = "black", max.overlaps = Inf) + 
  labs(title = "HT/LT THP1 infected strains, TPM filtered (Rv genes only)",
       subtitle = "ProbeTest4: Half library prep, half blockers, hybridization restarted after power cut",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_fig
# ggsave(PCA_fig,
#        file = paste0("PCA_v1.pdf"),
#        path = "HTLT_Figures",
#        width = 6, height = 5, units = "in")

# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Group,
                  colors = my_fav_colors)
PCA_3D



