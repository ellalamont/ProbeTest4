# Make a correlation scatter plot for ALL comparisons
# E. Lamont 
# 12/15/24

# DON't RERUN!! Will make many excess graphs

source("Import_data.R") # to get my_tpm and Sept_tpm
my_tpm$Gene <- rownames(my_tpm)
Sept_tpm$Gene <- rownames(Sept_tpm)

# Log10 transform the data
my_tpm_Log10 <- my_tpm %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Remove the high/low transmission strains
my_tpm_Log10 <- my_tpm_Log10 %>% select(-contains("_Thp"))

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
############### LOOP THROUGH ALL NOV SAMPLES ##############

my_path <- "All_Correlation_Scatter_Figures"

for (i in 1:(length(colnames(my_tpm_Log10)) -1 -1)) { 
  for (j in (i + 1):(length(colnames(my_tpm_Log10))-1)) {
    if (i != j) { # Avoid comparing the same samples or repeating comparisons
    # Access the samples
      Sample1 <- colnames(my_tpm_Log10[i])
      Sample2 <- colnames(my_tpm_Log10[j])
      # cat("Comparing:", Sample1, "with", Sample2, "\n")
      filename <- paste0(Sample1, "_ComparedTo_", Sample2, ".pdf")
      
      ScatterCorr <- my_tpm_Log10 %>% 
        ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
        geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
        labs(title = paste0(Sample1, " vs ", Sample2),
             subtitle = "Pearson correlation",
             x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) + 
        stat_cor(method="pearson") + # add a correlation to the plot
        my_plot_themes
      
      ggsave(ScatterCorr,
             file = filename,
             path = my_path,
             width = 7, height = 5, units = "in")
    }
  }
}
# This works but it obviously making a bunch that I don't need, will go through and delete, so don't rerun!!


###########################################################
############# TRY CHOOSING VIA VOLCANO NAMES ##############

input_string <- "newTHP1_1e5_1_DualrRNA_ComparedTo_newTHP1_1e5_2_DualrRNA"
input_string <- "W0.S354850_DualrRNA_ComparedTo_W0.S354851_DualrRNA"

# Split the string using "_ComparedTo_" as the delimiter
samples <- strsplit(input_string, "_ComparedTo_")[[1]]

# Function to clean the sample names
clean_sample <- function(sample) {
  # Remove the "W#." prefix and add an underscore between "S" and the number
  gsub("W[0-9]+\\.", "", sample) %>%           # Remove the "W#."
    gsub("S([0-9]+)", "S_\\1", .)              # Add underscore between "S" and the number
}

# Assign the two parts to separate variables
sample1 <- clean_sample(samples[1])
sample2 <- clean_sample(samples[2])















