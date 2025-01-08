# All the Sputum samples from the Sept and Nov sequencing runs 
# Why are some lower read counts than others in the same week? Compare metadata
# E. Lamont 
# 1/8/25

# Compare total_RNA_ng, mRNA_ng, ct, ttd, Hyb_Time, Hyb_Restarted (maybe), 
# Using N_Genomic or P_Genomic
# Facet by Week and just look at W0 and W2!

source("Import_data.R") # to get AllSputum_pipeSummary

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
################### N_GENOMIC vs WEEK #####################

