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

poster_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=20), 
        axis.text.x = element_text(angle = 0, size=20, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
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
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "SCALED Log10 transformed, Pearson correlation",
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, "_v3.pdf"),
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
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") +
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
##################### SCALED AVERAGES #####################

ProbeTest4_THP1_tpm <- my_tpm %>% select(THP1_1e6_1, THP1_1e6_2, THP1_1e6_3, THP1_1e6_4, THP1_1e6_5, Gene)
ProbeTest3_THP1_tpm <- Sept_tpm %>% select(THP1_1e6_1_Probe_3D_100, THP1_1e6_2_Probe_3D_50, THP1_1e6_3_Probe_3D_25, THP1_1e6_4_Probe_3D_10, THP1_1e6_5_Probe_1, Gene)

THP1_SCALED_Combined <- inner_join(ProbeTest4_THP1_tpm, ProbeTest3_THP1_tpm, by = "Gene")

THP1_SCALED_Combined_Log10 <- THP1_SCALED_Combined %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Add average columns
THP1_SCALED_Combined_Log10 <- THP1_SCALED_Combined_Log10 %>% mutate(
  ProbeTest3_Averages = rowMeans(select(., c(THP1_1e6_1_Probe_3D_100, THP1_1e6_2_Probe_3D_50, THP1_1e6_3_Probe_3D_25, THP1_1e6_4_Probe_3D_10, THP1_1e6_5_Probe_1)), na.rm = TRUE),
  ProbeTest4_Averages = rowMeans(select(., c(THP1_1e6_1, THP1_1e6_2, THP1_1e6_3, THP1_1e6_4, THP1_1e6_5)), na.rm = TRUE),
)

# Compare the averages! 
Sample1 <- "ProbeTest4_Averages" # Captured
Sample2 <- "ProbeTest3_Averages" # Not Captured
ScatterCorr <- THP1_SCALED_Combined_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("THP1 ProbeTest 3 vs 4: Not scaled Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; 5 samples: THP1 1e6 Ra spiked ",
       x = paste0("SCALED Log10(TPM+1) ProbeTest4 samples averaged"), y = paste0("SCALED Log10(TPM+1) ProbeTest3 averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, "_SCALED.pdf"),
       path = "CompareAcrossRuns_Figures",
       width = 7, height = 5, units = "in")


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
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Not Scaled! Log10 transformed, Pearson correlation",
       x = paste0(Sample1, " Log10(TPM+1)"), y = paste0(Sample2, " Log10(TPM+1)")) + 
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

# S_250754 (Nov) and S_250754_Probe_4A_50 (Sept)
Sample1 <- "S_250754"
Sample2 <- "S_250754_Probe_4A_50"
ScatterCorr <- multiRun_tpm_NotScaled_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Not Scaled! Log10 transformed, Pearson correlation",
       x = paste0(Sample1, " Log10(TPM+1)"), y = paste0(Sample2, " Log10(TPM+1)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, "_NotScaled.pdf"),
       path = "CompareAcrossRuns_Figures",
       width = 7, height = 5, units = "in")
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, "_NotScaled.png"),
       path = "CompareAcrossRuns_Figures",
       width = 7, height = 5, units = "in")


###########################################################
################# NOT SCALED THP1 AVERAGES ################

ProbeTest4_THP1_NotScaled_tpm <- my_tpm_NotScaled %>% select(THP1_1e6_1, THP1_1e6_2, THP1_1e6_3, THP1_1e6_4, THP1_1e6_5, Gene)
ProbeTest3_THP1_NotScaled_tpm <- Sept_tpm_NotScaled %>% select(THP1_1e6_1_Probe_3D_100, THP1_1e6_2_Probe_3D_50, THP1_1e6_3_Probe_3D_25, THP1_1e6_4_Probe_3D_10, THP1_1e6_5_Probe_1, Gene)

THP1_Combined <- inner_join(ProbeTest4_THP1_NotScaled_tpm, ProbeTest3_THP1_NotScaled_tpm, by = "Gene")

THP1_Combined_Log10 <- THP1_Combined %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Add average columns
THP1_Combined_Log10 <- THP1_Combined_Log10 %>% mutate(
    ProbeTest3_Averages = rowMeans(select(., c(THP1_1e6_1_Probe_3D_100, THP1_1e6_2_Probe_3D_50, THP1_1e6_3_Probe_3D_25, THP1_1e6_4_Probe_3D_10, THP1_1e6_5_Probe_1)), na.rm = TRUE),
    ProbeTest4_Averages = rowMeans(select(., c(THP1_1e6_1, THP1_1e6_2, THP1_1e6_3, THP1_1e6_4, THP1_1e6_5)), na.rm = TRUE),
  )

# Compare the averages! 
Sample1 <- "ProbeTest4_Averages" # Captured
Sample2 <- "ProbeTest3_Averages" # Not Captured
ScatterCorr <- THP1_Combined_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("THP1 ProbeTest 3 vs 4: Not scaled Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; 5 samples: THP1 1e6 Ra spiked ",
       x = paste0("Log10(TPM+1) ProbeTest4 samples averaged"), y = paste0("Log10(TPM+1) ProbeTest3 averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, "_NotScaled.pdf"),
       path = "CompareAcrossRuns_Figures",
       width = 7, height = 5, units = "in")

# For poster
ScatterCorr_poster <- THP1_Combined_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.55, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = NULL,
       subtitle = NULL,
       x = paste0("Captured with probe set A Log10(TPM+1)"), y = paste0("Captured with probe set B \nLog10(TPM+1)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  poster_plot_themes
ScatterCorr_poster
ggsave(ScatterCorr_poster,
       file = paste0("ScatterCorr_", Sample1, "_vs_", Sample2, "_NotScaled.pdf"),
       path = "Poster_Figures",
       width = 7.5, height = 4.5, units = "in")



# THP1_1e6_1 vs THP1_1e6_1_Probe_3D_100
Sample1 <- "THP1_1e6_5" # Captured
Sample2 <- "THP1_1e6_5_Probe_1" # Not Captured
ScatterCorr <- THP1_Combined_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") +
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

# THP1_1e6_2 vs THP1_1e6_2_Probe_3D_50
Sample1 <- "THP1_1e6_2" # Captured
Sample2 <- "THP1_1e6_2_Probe_3D_50" # Not Captured
ScatterCorr <- THP1_Combined_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") +
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

# THP1_1e6_2 vs THP1_1e6_2_Probe_3D_50
Sample1 <- "THP1_1e6_2" # Captured
Sample2 <- "THP1_1e6_2_Probe_3D_50" # Not Captured
ScatterCorr <- THP1_Combined_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") +
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


###########################################################
#################### WEIRD GENE CLUSTER ###################

# There is a weird gene cluster of genes that are more expressed in ProbeTest3 compared to ProbeTest4
# Pulled the names out by hand from the averages correlation plot

# Weird extra gene cluster:
WeirdGeneList <- c("MT3275.1","MT2726", "Rv2106", "Rv2649", "Rv3475", "Rv2355", "Rv3185", "Rv0796", "Rv3187", "Rv2279", "Rv1764", "Rv3326", "Rv2479c", "Rv1756c", "Rv3380c", "Rv2167c", "Rv1369c", "Rv2814c", "Rv2105", "Rv2648", "Rv3186", "Rv3474", "Rv2354", "Rv3184", "Rv0795", "Rv2278", "Rv1763", "Rv3325", "Rv1370c", "Rv3381c", "Rv2168c", "Rv1757c", "Rv2480c", "Rv2815c")

# Import the scaling information
ScalingInfo <- read.csv("Pulldown.Enrichment.FinalScaleFactor.csv")

# See what the scale factor is for the weird genes
WeirdGenes_ScalingInfo <- ScalingInfo %>% filter(GENE_ID %in% WeirdGeneList)

# Add a column to the ScalingInfo with the weird genes
ScalingInfo <- ScalingInfo %>% 
  mutate(WeirdGene = ifelse(GENE_ID %in% WeirdGeneList, "Yes", "No"))

# Okay... These are all ones with really low scale factors! But why am I seeing this when comparing two captured samples
# These are also all transposases
# I see this comparing all ProbeTest3 to ProbeTest4. Not about probe concentration (which I would have seen within ProbeTest3) but something to do with the different DNA probes?
## So maybe there is an effect of probe prep? DNA for ProbeTest4 probes has less of these transposase genes?
