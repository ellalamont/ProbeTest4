# Do corrlation plot of the differential expression between different samples
# E. Lamont
# 1/20/25

# Like in this paper: https://www.nature.com/articles/s41598-019-55633-6
# Lastly, to confirm the ability of PatH-Cap to preserve differences in gene expression between two conditions, we applied PatH-Cap to RNA-seq libraries made from axenic antibiotic-treated bacterial cultures and untreated controls

# This actually won't work because I need the same think to compare each to. Will be able to do this maybe when I have the broth results!

source("Import_data.R") # to get list_dfs_2


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
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )









