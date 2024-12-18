# Read depth for all the Sputum samples from the Sept and Nov sequencing runs 
# E. Lamont 
# 12/17/24

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

WeekvsReads_sputum1 <- AllSputum_pipeSummary %>% 
  
  ggplot(aes(x = Week, y = N_Genomic)) + 
  geom_point(aes(fill = Sputum_Number, shape = Sputum_Number), shape = 21, size = 6, alpha = 0.8, stroke = 0.8) + 
  # scale_fill_manual("Sample type", values = c(`Marmoset` = "#CAB2D6", `Sputum` = "#0072B2", `Saliva` = "#009E73", `THP1` = "#FF7F00")) +  
  # scale_shape_manual(values=c(`0` = 15, `2` = 0, `4` = 3)) + 
  
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  # geom_text(aes(label = Probe_ng), size= 1.5, nudge_x = 0.07) + 
  
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  
  # facet_grid(~ Sample_Type, scales = "free", space = "free") + 
  
  # scale_y_continuous(limits = c(0,5000000), breaks = seq(0, 5000000, 1000000)) +
  
  labs(title = "Sputum: Week vs number reads aligned to Mtb",
       subtitle = "Label is number of reads aligned to Mtb", 
       x = "Weeks after start of antibiotics", 
       y = "# reads aligning to Mtb genome") + 
  
  my_plot_themes

WeekvsReads_sputum1







