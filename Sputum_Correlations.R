# Compare all sputum correlations
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
################# LOG10 GGCORRPLOT SPUTUM #################
# Only going to do Pearson because the number of genes is so high CLT applies and parametric tests can be used

my_tpm_AllSputum_Log10 <- my_tpm %>% 
  select(contains("S_")) %>%
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values
my_cor_pearson_AllSputum_Log10 <- cor(my_tpm_AllSputum_Log10, method = "pearson")

min(my_cor_pearson_AllSputum_Log10) # 0.09304044

# Plot pearson
pearson_plot_AllSputum_Log10 <- my_cor_pearson_AllSputum_Log10 %>% 
  ggcorrplot(hc.order = FALSE, 
             lab = TRUE, lab_size = 2,
             type = c("full")) + 
  scale_fill_gradient2(limit = c(0.09,1), low = "blue", high =  "red", mid = "white", midpoint = 0.55) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Sputum Pearson Correlation Log10 transformed", 
       subtitle = "All sputum sequenced", 
       fill = "Correlation")
pearson_plot_AllSputum_Log10

ggsave(pearson_plot_AllSputum_Log10,
       file = "Pearson_plot_AllSputum_Log100.pdf",
       path = "Sputum_Correlations_Figures",
       width = 8, height = 8, units = "in")


###########################################################
######################## S_503557 #########################
# Week 0

# S_503557_DualrRNA
# S_503557 - This sample was library prepped for the september run and recaptured here (MtbrRNA)

S_503557_CorrelationLog10 <- my_tpm_Log10 %>% 
  ggplot(aes(x = S_503557, y = S_503557_DualrRNA)) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "S_503557 (Week 0) Pearson correlation Log10 Transformmed",
       subtitle = "S_503557 (MtbrRNA) vs S_503557_DualrRNA",
       x = "S_503557 Log10(TPM)", y = "S_503557_DualrRNA Log10(TPM)") + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
S_503557_CorrelationLog10
# ggplotly(S_503557_CorrelationLog10)

ggsave(S_503557_CorrelationLog10,
       file = "S_503557_CorrelationLog10.pdf",
       path = "Sputum_Correlations_Figures",
       width = 7, height = 5, units = "in")


###########################################################
######################## S_575533 #########################
# Week 2

# S_575533_DualrRNA
# S_575533_MtbrRNA - This sample was library prepped for the september run and recaptured here (MtbrRNA)

S_575533_CorrelationLog10 <- my_tpm_Log10 %>% 
  ggplot(aes(x = S_575533_MtbrRNA, y = S_575533_DualrRNA)) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = "S_575533 (Week 2) Pearson correlation Log10 Transformmed",
       subtitle = "S_575533_MtbrRNA vs S_575533_DualrRNA",
       x = "S_575533_MtbrRNA Log10(TPM)", y = "S_575533_DualrRNA Log10(TPM)") + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
S_575533_CorrelationLog10
# ggplotly(S_503557_CorrelationLog10)

ggsave(S_575533_CorrelationLog10,
       file = "S_575533_CorrelationLog10.pdf",
       path = "Sputum_Correlations_Figures",
       width = 7, height = 5, units = "in")




