# Compare TPMs of samples sequenced on the Sept and Nov runs (Scatterplot)
# E. Lamont
# 12/11/24

source("Import_data.R") # to get my_tpm and Sept_tpm

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default


my_tpm$Gene <- rownames(my_tpm)
Sept_tpm$Gene <- rownames(Sept_tpm)
my_tpm_NotScaled$Gene <- rownames(my_tpm_NotScaled)
Sept_tpm_NotScaled$Gene <- rownames(Sept_tpm_NotScaled)

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=14), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )


###########################################################
##################### NEW MERGED DF #######################
# Make a new dataframe with just the sputum samples that were run twice
# These are the exact same library prepped samples that were re-captured and run (no difference in rRNA depletion!)
# Also look at THP1 samples that were run twice!

# S_575533 (W2), S_687338 (W4), S_503557 (W0), S_250754 (W0)
# THP1_1e6_3, THP1_1e6_3_Probe_3D_25

my_tpm_RunSubset <- my_tpm %>%
  select(S_575533_MtbrRNA, S_687338_MtbrRNA, S_503557, S_250754, THP1_1e6_3, Gene)
Sept_tpm_RunSubset <- Sept_tpm %>% 
  select(S_575533_Probe_3A, S_687338_Probe_4A_100, S_503557_Probe_3D_10, S_250754_Probe_4A_50, THP1_1e6_3_Probe_3D_25, Gene)

# Merge the dataframes
multiRun_tpm <- merge(my_tpm_RunSubset, Sept_tpm_RunSubset, by = "Gene", all = T)

# Log10 transform the data
multiRun_tpm_Log10 <- multiRun_tpm %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values


###########################################################
######### THP1_1e6_3 vs THP1_1e6_3_Probe_3D_25 ############

Sample1 <- "THP1_1e6_3"
Sample2 <- "THP1_1e6_3_Probe_3D_25"
ScatterCorr <- multiRun_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Log10 transformed, Pearson correlation",
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, "_v2.pdf"),
       path = "CompareAcrossRuns_Figures",
       width = 6, height = 4, units = "in")


###########################################################
######################## S_503557 #########################
# S_503557 (Nov) and S_503557_Probe_3D_10 (Sept)

Sample1 <- "S_503557"
Sample2 <- "S_503557_Probe_3D_10"
ScatterCorr <- multiRun_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Log10 transformed, Pearson correlation",
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, ".pdf"),
       path = "CompareAcrossRuns_Figures",
       width = 7, height = 5, units = "in")



###########################################################
######################## S_575533 #########################
# S_575533_MtbrRNA (Nov) and S_575533_Probe_3A (Sept)

Sample1 <- "S_575533_MtbrRNA"
Sample2 <- "S_575533_Probe_3A"
ScatterCorr <- multiRun_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Log10 transformed, Pearson correlation",
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, ".pdf"),
       path = "CompareAcrossRuns_Figures",
       width = 7, height = 5, units = "in")


###########################################################
######################## S_250754 #########################
# S_250754 (Nov) and S_250754_Probe_4A_50 (Sept)

Sample1 <- "S_250754"
Sample2 <- "S_250754_Probe_4A_50"
ScatterCorr <- multiRun_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Log10 transformed, Pearson correlation",
       x = paste0(Sample1, " Log10(TPM+1)"), y = paste0(Sample2, " Log10(TPM+1)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, "_v2.pdf"),
       path = "CompareAcrossRuns_Figures",
       width = 6, height = 4, units = "in")

###########################################################
####################### NOT SCALED ########################
# Compare the same samples that are not scaled

# Make a new merged DF
my_tpm_NotScaled_RunSubset <- my_tpm_NotScaled %>%
  select(S_575533_MtbrRNA, S_687338_MtbrRNA, S_503557, S_250754, THP1_1e6_3, Gene)
Sept_tpm_NotScaled_RunSubset <- Sept_tpm_NotScaled %>% 
  select(S_575533_Probe_3A, S_687338_Probe_4A_100, S_503557_Probe_3D_10, S_250754_Probe_4A_50, THP1_1e6_3_Probe_3D_25, Gene)

# Merge the dataframes
multiRun_tpm_NotScaled <- merge(my_tpm_NotScaled_RunSubset, Sept_tpm_NotScaled_RunSubset, by = "Gene", all = T)

# Log10 transform the data
multiRun_tpm_NotScaled_Log10 <- multiRun_tpm_NotScaled %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# THP1_1e6_3 vs THP1_1e6_3_Probe_3D_25 
Sample1 <- "THP1_1e6_3"
Sample2 <- "THP1_1e6_3_Probe_3D_25"
ScatterCorr <- multiRun_tpm_NotScaled_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Not Scaled! Log10 transformed, Pearson correlation",
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, "_NotScaled.pdf"),
       path = "CompareAcrossRuns_Figures",
       width = 7, height = 5, units = "in")

# S_503557 (Nov) and S_503557_Probe_3D_10 (Sept)
Sample1 <- "S_503557"
Sample2 <- "S_503557_Probe_3D_10"
ScatterCorr <- multiRun_tpm_NotScaled_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Not Scaled! Log10 transformed, Pearson correlation",
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, "_NotScaled.pdf"),
       path = "CompareAcrossRuns_Figures",
       width = 7, height = 5, units = "in")


# S_575533_MtbrRNA (Nov) and S_575533_Probe_3A (Sept)
Sample1 <- "S_575533_MtbrRNA"
Sample2 <- "S_575533_Probe_3A"
ScatterCorr <- multiRun_tpm_NotScaled_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Not Scaled! Log10 transformed, Pearson correlation",
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, "_NotScaled.pdf"),
       path = "CompareAcrossRuns_Figures",
       width = 7, height = 5, units = "in")

