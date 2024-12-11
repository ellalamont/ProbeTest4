# Compare TPMs of samples sequenced on the Sept and Nov runs
# E. Lamont
# 12/11/24

source("Import_data.R") # to get my_tpm and Sept_tpm
my_tpm$Gene <- rownames(my_tpm)
Sept_tpm$Gene <- rownames(Sept_tpm)

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
##################### NEW MERGED DF #######################
# Make a new dataframe with just the sputum samples that were run twice
# These are the exact same library prepped samples that were re-captured and run (no difference in rRNA depletion!)

# S_575533 (W2), S_687338 (W4), S_503557 (W0), S_250754 (W0)

my_tpm_RunSubset <- my_tpm %>%
  select(S_575533_MtbrRNA, S_687338_MtbrRNA, S_503557, S_250754)
Sept_tpm_RunSubset <- Sept_tpm %>% 
  select(S_575533_Probe_3A, S_687338_Probe_4A_100, S_503557_Probe_3D_10, S_250754_Probe_4A_50)

# Merge the dataframes
multiRun_tpm <- merge(my_tpm_RunSubset, Sept_tpm_RunSubset, by = "row.names", all = T)

###########################################################
######################## S_503557 #########################
# Week 0




