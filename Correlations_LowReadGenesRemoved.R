# Correlations with low read depth genes removed
# E. Lamont
# 12/17/24

# TESTING ONLY!!! Not appropriate for analysis!!

# Are the hight DEG the ones that have few or no reads aligning to some of the samples?


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

# Samples to try with
# orginalTHP1_1e6_1_DualrRNA - The sample was fully library prepped for the september sequencing AND the hydridization was good 16hr
# originalTHP1_1e6_6_DualrRNA - This sample was Half library prep done by JA in Nov, AND the hybridization was stopped then restarted


###############################################################################
############################### REMOVE GENES ##################################

test_tpm <- my_tpm %>% select(orginalTHP1_1e6_1_DualrRNA, originalTHP1_1e6_6_DualrRNA, Gene)
# test_tpm_GenesRemoved <- test_tpm %>%
#   filter(!(orginalTHP1_1e6_1_DualrRNA == 0 & originalTHP1_1e6_6_DualrRNA < 10)) & !# (orginalTHP1_1e6_1_DualrRNA < 10 & originalTHP1_1e6_6_DualrRNA == 0)

# Remove Genes with tpm 0 and < 10 for the two columns
test_tpm_GenesRemoved <- test_tpm %>%
  filter(orginalTHP1_1e6_1_DualrRNA >= 10) %>%
  filter(originalTHP1_1e6_6_DualrRNA >= 10)
  # filter(!(orginalTHP1_1e6_1_DualrRNA == 0 | originalTHP1_1e6_6_DualrRNA == 0)) # %>%
  # filter(!(orginalTHP1_1e6_1_DualrRNA < 10 & originalTHP1_1e6_6_DualrRNA == 0))

# Log10 transform
test_tpm_GenesRemoved_Log10 <- test_tpm_GenesRemoved %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values


###############################################################################
######################### MAKE CORRELATION SCATTER ############################

Sample1 <- "orginalTHP1_1e6_1_DualrRNA"
Sample2 <- "originalTHP1_1e6_6_DualrRNA"
ScatterCorr <- test_tpm_GenesRemoved_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Genes removed if TPM < 10 in either sample; Pearson correlation",
       x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("GENES.REMOVED.ScatterCorr_", Sample1, "_vs_", Sample2, ".pdf"),
#        path = "GenesRemoved_Figures",
#        width = 7, height = 5, units = "in")


###############################################################################
########################### make_volcano_function #############################

make_volcano_function <- function(my_df, graph_title) {
  
  ## Make a volcano plot using output from Bob's pipeline
  
  my_volcano <- my_df %>%
    ggplot(aes(x = LOG2FOLD, y = -log10(AVG_PVALUE), col = DE, label = DE_labels, text = GENE_NAME, label2 = GENE_ID)) + # text is for plotly, could be GENE_ID
    geom_point() + 
    labs(title = graph_title) + 
    geom_vline(xintercept = c(-1,1), col = "grey", linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
    
    # Need it this way so the colors aren't messed up by not having significant up or down
    # scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) + 
    scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00")) +
    
    geom_text_repel(max.overlaps = 10, size = 3) # Can do geom_text_repel or geom_label_rebel
  
  # Determine the max and min axes values for labeling 
  plot_build <- ggplot_build(my_volcano)
  y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
  x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
  x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
  
  # Add the gene number annotations
  text_up <- my_df %>% filter(DE == "significant up") %>% nrow()
  text_down <- my_df %>% filter(DE == "significant down") %>% nrow()
  my_volcano_annotated <- my_volcano +
    annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
    annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
  
  final_volcano <- my_volcano_annotated + my_plot_themes
  
}

###############################################################################
############################# MAKE VOLCANO PLOT ###############################

# orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA is volcano plot name

# Remove the genes that are not in test_tpm_GenesRemoved 
test <- make_volcano_function(list_dfs_2[["orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA"]] %>% 
                                filter(GENE_ID %in% test_tpm_GenesRemoved$Gene),
                              "orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA")
test
