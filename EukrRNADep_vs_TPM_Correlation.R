# Compare eukaryotic rRNA depletion to TPM with correlation plot
# E. Lamont
# 12/11/24

source("Import_data.R") # to get my_tpm
my_tpm$Gene <- rownames(my_tpm)

# Log10 transform the data
my_tpm_Log10 <- my_tpm %>% 
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

###########################################################
########## THP1 HALF LIBPREP CORRELATION SAMPLES ##########

# SampleIDs are: originalTHP1_1e6_4_DualrRNA, originalTHP1_1e6_5_DualrRNA, originalTHP1_1e6_6_DualrRNA, orginalTHP1_1e6_1_MtbrRNA, orginalTHP1_1e6_2_MtbrRNA, orginalTHP1_1e6_3_MtbrRNA
# all are the same, and are half library prepped
# Some have been pooled and some not, but I don't think that makes much of a difference


###########################################################
################ GGCORRPLOT THP1 HALF PREP ################
# Only going to do Pearson because the number of genes is so high CLT applies and parametric tests can be used

my_tpm_THP1Subset <- my_tpm %>% select(originalTHP1_1e6_4_DualrRNA, originalTHP1_1e6_5_DualrRNA, originalTHP1_1e6_6_DualrRNA, orginalTHP1_1e6_1_MtbrRNA, orginalTHP1_1e6_2_MtbrRNA, orginalTHP1_1e6_3_MtbrRNA)
my_cor_pearson_THP1Subset <- cor(my_tpm_THP1Subset, method = "pearson")

min(my_cor_pearson_THP1Subset) # 0.9518504

# Plot pearson
pearson_plot_THP1Subset <- my_cor_pearson_THP1Subset %>% 
  ggcorrplot(hc.order = FALSE, 
             lab = TRUE, lab_size = 1.5,
             type = c("full")) + 
  scale_fill_gradient2(limit = c(0.95,1), low = "blue", high =  "red", mid = "white", midpoint = 0.97) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "THP1 Pearson Correlation", 
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra, All half library prepped", 
       fill = "Correlation")
pearson_plot_THP1Subset

ggsave(pearson_plot_THP1Subset,
       file = "Pearson_CorrelationPlot_THP1Subset.pdf",
       path = "EukrRNADepletion_Figures",
       width = 8, height = 8, units = "in")



###########################################################
############# LOG10 GGCORRPLOT THP1 HALF PREP #############
# Only going to do Pearson because the number of genes is so high CLT applies and parametric tests can be used

my_tpm_THP1Subset_Log10 <- my_tpm %>% 
  select(originalTHP1_1e6_4_DualrRNA, originalTHP1_1e6_5_DualrRNA, originalTHP1_1e6_6_DualrRNA, orginalTHP1_1e6_1_MtbrRNA, orginalTHP1_1e6_2_MtbrRNA, orginalTHP1_1e6_3_MtbrRNA) %>%
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values
my_cor_pearson_THP1Subset_Log10 <- cor(my_tpm_THP1Subset_Log10, method = "pearson")

min(my_cor_pearson_THP1Subset_Log10) # 0.7788043

# Plot pearson
pearson_plot_THP1Subset_Log10 <- my_cor_pearson_THP1Subset_Log10 %>% 
  ggcorrplot(hc.order = FALSE, 
             lab = TRUE, lab_size = 4,
             type = c("full")) + 
  scale_fill_gradient2(limit = c(0.75,1), low = "blue", high =  "red", mid = "white", midpoint = 0.875) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "THP1 Pearson Correlation Log10 transformed", 
       subtitle = "originalTHP1 spiked with 1e6 cells H37Ra, All half library prepped", 
       fill = "Correlation")
pearson_plot_THP1Subset_Log10

ggsave(pearson_plot_THP1Subset_Log10,
       file = "Pearson_CorrelationPlot_THP1_rRNADep_Log10.pdf",
       path = "EukrRNADepletion_Figures",
       width = 8, height = 8, units = "in")


###########################################################
########## SCATTER THP1 Mtb vs Dual Depletion #############

# Just look at the not-pooled ones as a comparison
# originalTHP1_1e6_6_DualrRNA vs orginalTHP1_1e6_3_MtbrRNA

originalTHP1_1e6_DualDepvsMtbDep_Log10 <- my_tpm_Log10 %>% 
  ggplot(aes(x = originalTHP1_1e6_6_DualrRNA, y = orginalTHP1_1e6_3_MtbrRNA)) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = "orginalTHP1_1e6_3_MtbrRNA vs originalTHP1_1e6_6_DualrRNA; unpooled, Half library prep",
       x = "Log10(TPM) (Mtb+human rRNA depletion)", y = "Log10(TPM) (Mtb rRNA depletion only)") + 
  
  # scale_y_continuous(limits = c(0,20000), breaks = seq(0, 20000, 5000)) +
  # scale_x_continuous(limits = c(0,20000), breaks = seq(0, 20000, 5000)) +
  
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
originalTHP1_1e6_DualDepvsMtbDep_Log10
ggplotly(originalTHP1_1e6_DualDepvsMtbDep_Log10)

ggsave(originalTHP1_1e6_DualDepvsMtbDep_Log10,
       file = "originalTHP1_1e6_DualDepvsMtbDep_Log10.pdf",
       path = "EukrRNADepletion_Figures",
       width = 6, height = 4, units = "in")



###########################################################
############ SCATTER THP1 Mtb depletion only ##############

# Just look at the not-pooled ones as a comparison
# orginalTHP1_1e6_2_MtbrRNA vs orginalTHP1_1e6_3_MtbrRNA

originalTHP1_1e6_MtbDep_Log10 <- my_tpm_Log10 %>% 
  ggplot(aes(x = orginalTHP1_1e6_2_MtbrRNA, y = orginalTHP1_1e6_3_MtbrRNA)) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "THP1 with 1e6 cells H37Ra Pearson correlation",
       subtitle = "orginalTHP1_1e6_3_MtbrRNA vs orginalTHP1_1e6_2_MtbrRNA; Half library prep",
       x = "orginalTHP1_1e6_2_MtbrRNA Log10(TPM) (Mtb rRNA depletion only)", y = "orginalTHP1_1e6_3_MtbrRNA Log10(TPM) (Mtb rRNA depletion only)") + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
originalTHP1_1e6_MtbDep_Log10
ggplotly(originalTHP1_1e6_MtbDep_Log10)

ggsave(originalTHP1_1e6_MtbDep_Log10,
       file = "originalTHP1_1e6_MtbDep_Log10.pdf",
       path = "EukrRNADepletion_Figures",
       width = 7, height = 5, units = "in")










