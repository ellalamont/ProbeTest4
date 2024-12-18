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
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = c14) +  
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,5500000), breaks = seq(0, 5500000, 1000000)) +
  labs(title = "Sputum: Week vs number reads aligned to Mtb",
       subtitle = "Label is number of reads aligned to Mtb", 
       x = "Weeks after start of antibiotics", 
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes
WeekvsReads_sputum1
ggplotly(WeekvsReads_sputum1)
ggsave(WeekvsReads_sputum1,
       file = "WeekvsReads_sputum1.pdf",
       path = "AllSputum_Figures",
       width = 8, height = 5, units = "in")


###########################################################
################### P_GENOMIC vs WEEK #####################

WeekvsPercent_sputum1 <- AllSputum_pipeSummary %>% 
  ggplot(aes(x = Week, y = P_Genomic)) + 
  geom_point(aes(fill = Sputum_Number, shape = SeqRun), size = 6, alpha = 0.7, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = c14) +  
  guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend to show fill colors
  scale_shape_manual(values=c(`Sept` = 22, `Nov` = 21)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,75), breaks = seq(0, 75, 10)) +
  labs(title = "Sputum: Week vs percent reads aligned to Mtb",
       subtitle = "Label is number of reads aligned to Mtb", 
       x = "Weeks after start of antibiotics", 
       y = "% reads aligning to Mtb genome") + 
  my_plot_themes
WeekvsPercent_sputum1
ggplotly(WeekvsPercent_sputum1)
ggsave(WeekvsPercent_sputum1,
       file = "WeekvsPercent_sputum1.pdf",
       path = "AllSputum_Figures",
       width = 8, height = 5, units = "in")
