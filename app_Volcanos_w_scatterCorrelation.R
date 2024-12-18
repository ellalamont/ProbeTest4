library(shiny)
library(gmodels)

source("Import_data.R") # for list_dfs_2 and my_tpm

# Generate the tpm file needed for the scatter plots
my_tpm$Gene <- rownames(my_tpm)
# Log10 transform the data
my_tpm_Log10 <- my_tpm %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=12), 
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(angle = 0, size=12, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12), 
        plot.subtitle = element_text(size=12), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank())

# Function to clean the sample names
clean_sample <- function(sample) {
  # Remove the "W#." prefix and add an underscore between "S" and the number
  gsub("W[0-9]+\\.", "", sample) %>%           # Remove the "W#."
    gsub("S([0-9]+)", "S_\\1", .)              # Add underscore between "S" and the number
}



# Define UI ----
ui <- fluidPage(
  titlePanel("EL All DEG from ProbeTest 4 12/13/24"),
  
  # First row
  fluidRow(
    column(width = 7,
           # h3("Choose the condition"),
           selectInput("my_comparison",
                       label = "Choose comparisons",
                       choices = df_names,
                       width = "100%"),
           textInput("my_GeneID", 
                     label = "Gene ID",
                     value = "Rv")
    )
  ),
  
  # Second row
  fluidRow(
    column(width = 6,
           plotlyOutput("volcano_plot",
                        width = "100%", height = "600px")),
    column(width = 6,
           plotlyOutput("ScatterCorrelation_plot",
                        width = "100%", height = "600px")),
    )
)


# Define server logic ----
server <- function(input, output) {
  
  # Volcano Plot
  output$volcano_plot <- renderPlotly({
    
    single_gene <- list_dfs_2[[input$my_comparison]] %>% 
      filter(GENE_ID == input$my_GeneID)
    
    my_volcano <- list_dfs_2[[input$my_comparison]] %>%
      ggplot(aes(x = LOG2FOLD, y = -log10(AVG_PVALUE), col = DE, label = DE_labels, text = GENE_ID)) + 
      geom_point() + 
      
      # Add a differently colored point
      geom_point(data = single_gene, color = "yellow", aes(col = DE, label = DE_labels, text = GENE_ID)) + 
      
      labs(title = input$my_comparison) + 
      geom_vline(xintercept = c(-1,1), col = "grey", linetype = "dashed") + 
      geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
      scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
    # geom_label_repel(max.overlaps = 10) # Can do geom_text_repel or geom_label_rebel
    
    # Determine the max and min axes values for labeling 
    plot_build <- ggplot_build(my_volcano)
    y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
    x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
    x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
    
    # Add the gene number annotations
    text_up <- list_dfs_2[[input$my_comparison]] %>% filter(DE == "significant up") %>% nrow()
    text_down <- list_dfs_2[[input$my_comparison]] %>% filter(DE == "significant down") %>% nrow()
    my_volcano_annotated <- my_volcano +
      annotate("text", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
      annotate("text", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
    
    final_plot <- my_volcano_annotated + my_plot_themes 
    final_plot
  })
  
  # Scatter Correlation Plot
  output$ScatterCorrelation_plot <- renderPlotly({
    
    # Filter for the selected gene in the scatter data
    single_gene <- my_tpm_Log10 %>% 
      filter(Gene == input$my_GeneID)
    
    input_string <- input$my_comparison
    samples <- strsplit(input_string, "_ComparedTo_")[[1]]
    Sample1 <- clean_sample(samples[1])
    Sample2 <- clean_sample(samples[2])
    
    ScatterCorr <- my_tpm_Log10 %>% 
      ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
      geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
      
      # Add a differently colored point
      geom_point(data = single_gene, color = "yellow", aes(text = Gene)) + 
      
      labs(title = paste0(Sample1, " vs ", Sample2),
           subtitle = "Pearson correlation",
           x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
      stat_cor(method="pearson") + # add a correlation to the plot
      my_plot_themes
    
    
    
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)