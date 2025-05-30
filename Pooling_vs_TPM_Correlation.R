# Compare pooling to TPM with correlation plot
# E. Lamont
# 12/12/24

source("Import_data.R") # to get my_tpm
my_tpm$Gene <- rownames(my_tpm)
my_tpm_NotScaled$Gene <- rownames(my_tpm_NotScaled)

# Log10 transform the data
my_tpm_Log10 <- my_tpm %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values
my_tpm_NotScaled_Log10 <- my_tpm_NotScaled %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

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

# http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2

##### D #####
 
###########################################################
################### D: LOG10 GGCORRPLOT ###################
# Only going to do Pearson because the number of genes is so high CLT applies and parametric tests can be used

# Samplese are 
# orginalTHP1_1e6_1_MtbrRNA
# orginalTHP1_1e6_2_MtbrRNA
# orginalTHP1_1e6_3_MtbrRNA

# Select the samples of interest
Pool.D_Log10 <- my_tpm_Log10 %>%
  select(orginalTHP1_1e6_1_MtbrRNA, orginalTHP1_1e6_2_MtbrRNA, orginalTHP1_1e6_3_MtbrRNA)

# Make the correlation
my_cor_pearson_Pool.D_Log10 <- cor(Pool.D_Log10, method = "pearson")

min(my_cor_pearson_Pool.D_Log10) # 0.901357

# Plot pearson
Pool.D_PearsonLog10 <- my_cor_pearson_Pool.D_Log10 %>% 
  ggcorrplot(hc.order = FALSE, 
             lab = TRUE, lab_size = 5,
             type = c("full")) + 
  scale_fill_gradient2(limit = c(0.9,1), low = "blue", high =  "red", mid = "white", midpoint = 0.95) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "THP1 Pearson Correlation Log10 transformed", 
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra, All half library prepped, MtbrRNADep \nSamples 1 and 2 were pooled", 
       fill = "Correlation")
Pool.D_PearsonLog10

ggsave(Pool.D_PearsonLog10,
       file = "Pool.D_PearsonLog10.pdf",
       path = "Pooling_Figures",
       width = 7, height = 6, units = "in")

##########################################################################################
############### D: orginalTHP1_1e6_1_MtbrRNA vs orginalTHP1_1e6_2_MtbrRNA ################

orginalTHP1_1e6_1_MtbrRNA_vs_orginalTHP1_1e6_2_MtbrRNA <- my_tpm_Log10 %>% 
  ggplot(aes(x = orginalTHP1_1e6_1_MtbrRNA, y = orginalTHP1_1e6_2_MtbrRNA)) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = "orginalTHP1_1e6_1_MtbrRNA vs orginalTHP1_1e6_2_MtbrRNA",
       x = "orginalTHP1_1e6_1_MtbrRNA Log10(TPM)", y = "orginalTHP1_1e6_2_MtbrRNA Log10(TPM)") + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
orginalTHP1_1e6_1_MtbrRNA_vs_orginalTHP1_1e6_2_MtbrRNA

ggsave(orginalTHP1_1e6_1_MtbrRNA_vs_orginalTHP1_1e6_2_MtbrRNA,
       file = "orginalTHP1_1e6_1_MtbrRNA_vs_orginalTHP1_1e6_2_MtbrRNA.pdf",
       path = "Pooling_Figures",
       width = 7, height = 5, units = "in")

##########################################################################################
############### D: orginalTHP1_1e6_1_MtbrRNA vs orginalTHP1_1e6_3_MtbrRNA ################

orginalTHP1_1e6_1_MtbrRNA_vs_orginalTHP1_1e6_3_MtbrRNA <- my_tpm_Log10 %>% 
  ggplot(aes(x = orginalTHP1_1e6_1_MtbrRNA, y = orginalTHP1_1e6_3_MtbrRNA)) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = "orginalTHP1_1e6_1_MtbrRNA vs orginalTHP1_1e6_3_MtbrRNA",
       x = "orginalTHP1_1e6_1_MtbrRNA Log10(TPM)", y = "orginalTHP1_1e6_3_MtbrRNA Log10(TPM)") + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
orginalTHP1_1e6_1_MtbrRNA_vs_orginalTHP1_1e6_3_MtbrRNA

ggsave(orginalTHP1_1e6_1_MtbrRNA_vs_orginalTHP1_1e6_3_MtbrRNA,
       file = "orginalTHP1_1e6_1_MtbrRNA_vs_orginalTHP1_1e6_3_MtbrRNA.pdf",
       path = "Pooling_Figures",
       width = 7, height = 5, units = "in")





##### F #####
###################################################################
######################## F: LOG10 GGCORRPLOT ######################
# Only going to do Pearson because the number of genes is so high CLT applies and parametric tests can be used

# Samplese are 
# THP1_1e6_1
# THP1_1e6_2
# THP1_1e6_3
# THP1_1e6_4
# THP1_1e6_5

# Select the samples of interest
Pool.F_Log10 <- my_tpm_Log10 %>%
  select(THP1_1e6_1, THP1_1e6_2, THP1_1e6_3, THP1_1e6_4, THP1_1e6_5)

# Make the correlation
my_cor_pearson_Pool.F_Log10 <- cor(Pool.F_Log10, method = "pearson")

min(my_cor_pearson_Pool.F_Log10) # 0.9848202

# Plot pearson
Pool.F_PearsonLog10 <- my_cor_pearson_Pool.F_Log10 %>% 
  ggcorrplot(hc.order = FALSE, 
             lab = TRUE, lab_size = 5,
             type = c("full")) + 
  scale_fill_gradient2(limit = c(0.98,1), low = "blue", high =  "red", mid = "white", midpoint = 0.99) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "F Pool: THP1 Pearson Correlation Log10 transformed", 
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra, All full library prepped, MtbrRNADep \nSamples 1-4 were pooled", 
       fill = "Correlation")
Pool.F_PearsonLog10

ggsave(Pool.F_PearsonLog10,
       file = "Pool.F_PearsonLog10.pdf",
       path = "Pooling_Figures",
       width = 7, height = 6, units = "in")

# F Scatter: THP1_1e6_1 vs THP1_1e6_2
ScatterCorr_THP1_1e6_1_vs_THP1_1e6_2 <- my_tpm_Log10 %>% 
  ggplot(aes(x = THP1_1e6_1, y = THP1_1e6_2)) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = "THP1_1e6_1 vs THP1_1e6_2; Full library prep",
       x = "THP1_1e6_1 Log10(TPM)", y = "THP1_1e6_2 Log10(TPM)") + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr_THP1_1e6_1_vs_THP1_1e6_2
ggsave(ScatterCorr_THP1_1e6_1_vs_THP1_1e6_2,
       file = "ScatterCorr_THP1_1e6_1_vs_THP1_1e6_2.pdf",
       path = "Pooling_Figures",
       width = 7, height = 5, units = "in")

# F Scatter: THP1_1e6_1 vs THP1_1e6_5
ScatterCorr_THP1_1e6_1_vs_THP1_1e6_5 <- my_tpm_Log10 %>% 
  ggplot(aes(x = THP1_1e6_1, y = THP1_1e6_5)) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = "THP1_1e6_1 vs THP1_1e6_5; Full library prep",
       x = "THP1_1e6_1 Log10(TPM)", y = "THP1_1e6_5 Log10(TPM)") + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr_THP1_1e6_1_vs_THP1_1e6_5
ggsave(ScatterCorr_THP1_1e6_1_vs_THP1_1e6_5,
       file = "ScatterCorr_THP1_1e6_1_vs_THP1_1e6_5.pdf",
       path = "Pooling_Figures",
       width = 7, height = 5, units = "in")

###################################################################
################ F: LOG10 CORRELATIONS NOT SCALED  ################
# Want to see what looks different when the data are not scaled

# F Scatter: THP1_1e6_1 vs THP1_1e6_5 NOT SCALED
ScatterCorr_THP1_1e6_1_vs_THP1_1e6_NotScaled <- my_tpm_NotScaled_Log10 %>% 
  ggplot(aes(x = THP1_1e6_1, y = THP1_1e6_5)) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "Not Scaled! THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = "THP1_1e6_1 vs THP1_1e6_5; Full library prep",
       x = "THP1_1e6_1 Log10(TPM)", y = "THP1_1e6_5 Log10(TPM)") + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr_THP1_1e6_1_vs_THP1_1e6_NotScaled
ggsave(ScatterCorr_THP1_1e6_1_vs_THP1_1e6_NotScaled,
       file = "ScatterCorr_THP1_1e6_1_vs_THP1_1e6_NotScaled.pdf",
       path = "Pooling_Figures",
       width = 7, height = 5, units = "in")



##### B
###################################################################
######################## B: LOG10 GGCORRPLOT ######################
# Only going to do Pearson because the number of genes is so high CLT applies and parametric tests can be used

# Samples are 
# orginalTHP1_1e6_1_DualrRNA
# originalTHP1_1e6_2_DualrRNA
# originalTHP1_1e6_3_DualrRNA

# Select the samples of interest
Pool.B_Log10 <- my_tpm_Log10 %>%
  select(orginalTHP1_1e6_1_DualrRNA, originalTHP1_1e6_2_DualrRNA, originalTHP1_1e6_3_DualrRNA)

# Make the correlation
my_cor_pearson_Pool.B_Log10 <- cor(Pool.B_Log10, method = "pearson")

min(my_cor_pearson_Pool.B_Log10) # 0.8326534

# Plot pearson
Pool.B_PearsonLog10 <- my_cor_pearson_Pool.B_Log10 %>% 
  ggcorrplot(hc.order = FALSE, 
             lab = TRUE, lab_size = 5,
             type = c("full")) + 
  scale_fill_gradient2(limit = c(0.83,1), low = "blue", high =  "red", mid = "white", midpoint = 0.91) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "B Pool: THP1 Pearson Correlation Log10 transformed", 
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra, All full library prepped, DualrRNADep \nSamples 2+3 were pooled", 
       fill = "Correlation")
Pool.B_PearsonLog10

ggsave(Pool.B_PearsonLog10,
       file = "Pool.B_PearsonLog10.pdf",
       path = "Pooling_Figures",
       width = 7, height = 6, units = "in")

# B Scatter: orginalTHP1_1e6_1_DualrRNA vs originalTHP1_1e6_2_DualrRNA
Sample1 <- "orginalTHP1_1e6_1_DualrRNA"
Sample2 <- "originalTHP1_1e6_2_DualrRNA"
ScatterCorr <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = paste0(Sample1, " vs ", Sample2),
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, ".pdf"),
       path = "Pooling_Figures",
       width = 7, height = 5, units = "in")

# B Scatter: orginalTHP1_1e6_1_DualrRNA vs originalTHP1_1e6_3_DualrRNA
Sample1 <- "orginalTHP1_1e6_1_DualrRNA"
Sample2 <- "originalTHP1_1e6_3_DualrRNA"
ScatterCorr <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = paste0(Sample1, " vs ", Sample2),
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, ".pdf"),
       path = "Pooling_Figures",
       width = 7, height = 5, units = "in")


# B Scatter: originalTHP1_1e6_2_DualrRNA vs originalTHP1_1e6_3_DualrRNA
Sample1 <- "originalTHP1_1e6_2_DualrRNA"
Sample2 <- "originalTHP1_1e6_3_DualrRNA"
ScatterCorr <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = paste0(Sample1, " vs ", Sample2),
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, ".pdf"),
       path = "Pooling_Figures",
       width = 7, height = 5, units = "in")


##### C
###################################################################
######################## C: LOG10 GGCORRPLOT ######################
# Only going to do Pearson because the number of genes is so high CLT applies and parametric tests can be used

# Samples are 
# originalTHP1_1e6_4_DualrRNA
# originalTHP1_1e6_5_DualrRNA
# originalTHP1_1e6_6_DualrRNA

# Select the samples of interest
Pool.C_Log10 <- my_tpm_Log10 %>%
  select(originalTHP1_1e6_4_DualrRNA, originalTHP1_1e6_5_DualrRNA, originalTHP1_1e6_6_DualrRNA)

# Make the correlation
my_cor_pearson_Pool.C_Log10 <- cor(Pool.C_Log10, method = "pearson")
min(my_cor_pearson_Pool.C_Log10) # 0.7788043

# Plot pearson
Pool.C_PearsonLog10 <- my_cor_pearson_Pool.C_Log10 %>% 
  ggcorrplot(hc.order = FALSE, 
             lab = TRUE, lab_size = 5,
             type = c("full")) + 
  scale_fill_gradient2(limit = c(0.77,1), low = "blue", high =  "red", mid = "white", midpoint = 0.88) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "C Pool: THP1 Pearson Correlation Log10 transformed", 
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra, All half library prepped, DualrRNADep \nSamples 4+5 were pooled", 
       fill = "Correlation")
Pool.C_PearsonLog10
ggsave(Pool.C_PearsonLog10,
       file = "Pool.C_PearsonLog10.pdf",
       path = "Pooling_Figures",
       width = 7, height = 6, units = "in")

# C Scatter: originalTHP1_1e6_4_DualrRNA vs originalTHP1_1e6_5_DualrRNA
Sample1 <- "originalTHP1_1e6_4_DualrRNA"
Sample2 <- "originalTHP1_1e6_5_DualrRNA"
ScatterCorr <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = paste0(Sample1, " vs ", Sample2),
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, ".pdf"),
       path = "Pooling_Figures",
       width = 7, height = 5, units = "in")

# C Scatter: originalTHP1_1e6_4_DualrRNA vs originalTHP1_1e6_6_DualrRNA
Sample1 <- "originalTHP1_1e6_4_DualrRNA"
Sample2 <- "originalTHP1_1e6_6_DualrRNA"
ScatterCorr <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = paste0(Sample1, " vs ", Sample2),
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, ".pdf"),
       path = "Pooling_Figures",
       width = 7, height = 5, units = "in")

# C Scatter: originalTHP1_1e6_4_DualrRNA vs originalTHP1_1e6_6_DualrRNA
Sample1 <- "originalTHP1_1e6_5_DualrRNA"
Sample2 <- "originalTHP1_1e6_6_DualrRNA"
ScatterCorr <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = paste0(Sample1, " vs ", Sample2),
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, ".pdf"),
       path = "Pooling_Figures",
       width = 7, height = 5, units = "in")




