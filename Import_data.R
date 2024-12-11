# Probe Prep testing 4
# Probes made by Jessica, run included THP1 test samples and sputum
# See Nov2024_RNALibraryPrep
# 12/10/24


################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(plotly)
library(ggprism) # for add_pvalue()
library(rstatix) # for adjust_pvalue
library(ggpmisc) # https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
library(ggrepel)
library(pheatmap)
library(ggplotify) # To convert pheatmaps to ggplots
library(corrplot)
library(ggcorrplot)
library(ggfortify) # To make pca plots with plotly
library(scales) # To add commas to the y axis... scale_y_continuous(label=comma, )



cbPalette_1 <- c("#999999", "#E69F00") # Gold and Grey
cbPalette_1.5 <- c("#E69F00", "#999999") # Gold and Grey
cbPalette_2 <- c( "#0072B2", "#999999") # Blue and Grey
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 <-  c("#bfbfbf", "#56B4E9")
cbPalette3 <-  c("#bfbfbf", "#E69F00")
cbPalette4 <- c("#56B4E9", "#009E73", "#F0E442")
cbPalette5 <- c("#009E73", "#FF7F00")
c25 <- c(
  "dodgerblue2", "#E31A1C", "green4",
  "#6A3D9A","#FF7F00","black", "gold1",
  "skyblue2", "#FB9A99","palegreen2","#CAB2D6",
  "#FDBF6F","gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown"
)
c12 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "palegreen2", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4") 
c16 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black","gold1", "#FB9A99", "#CAB2D6", "palegreen2", "gray70", "maroon", "orchid1", "blue1", "darkturquoise", "darkorange4") 
c11 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "palegreen2", "gold1", "maroon", "orchid1", "darkturquoise", "darkorange4", "gray70")
c7 <- c("gray70", "#E0F7FA", "#B2EBF2", "#81D4FA", "#03A9F4","#0288D1", "#01579B")
c8 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "palegreen2", "gray70", "maroon", "black")
Orange_Gradient <- c("#FFD4A3", "#FFA64D", "#FF7F00", "#E66900", "#993D00")

# Stop scientific notation
options(scipen = 999) 
# options(scipen = 0) # To revert back to default

###########################################################
############### IMPORT PIPELINE SUMMARY DATA ##############

my_pipeSummary <- read.csv("Pipeline.Summary.Details.csv")
# This has been edited to include more metadata!

my_pipeSummary$Hyb_Time <- as.character(my_pipeSummary$Hyb_Time)
ordered_Hyb_Time <- c("4", "16")
my_pipeSummary$Hyb_Time <- factor(my_pipeSummary$Hyb_Time, levels = ordered_Hyb_Time)

my_pipeSummary$Week <- as.character(my_pipeSummary$Week)
ordered_Week <- c("0", "2", "4")
my_pipeSummary$Week <- factor(my_pipeSummary$Week, levels = ordered_Week)

my_pipeSummary$EukrRNADep <- as.character(my_pipeSummary$EukrRNADep)
ordered_EukrRNADep <- c("MtbrRNA", "DualrRNA")
my_pipeSummary$EukrRNADep <- factor(my_pipeSummary$EukrRNADep, levels = ordered_EukrRNADep)

my_pipeSummary <- my_pipeSummary[-nrow(my_pipeSummary), ] # Remove the undetermined, which is the last row

my_pipeSummary$SampleID <- gsub(x = my_pipeSummary$SampleID, pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

###########################################################
############ IMPORT AND PROCESS ALL TPM VALUES ############

my_tpm <- read.csv("Mtb.Expression.Gene.Data.SCALED.TPM.csv")

my_tpm <- my_tpm[,-ncol(my_tpm)] # remove the last column which is the Undetermined

# Adjust the names so they are slightly shorter
# names(my_tpm) <- gsub(x = names(my_tpm), pattern = "_2_2_ng_mL", replacement = "")
names(my_tpm) <- gsub(x = names(my_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# Grab the metadata I added to my_pipeSummary
my_metadata <- my_pipeSummary %>% select(2, 14:26)

# Adjust the metadata names so they are the same
my_metadata$SampleID <- sub(x = my_metadata$SampleID, pattern = "_S.*", replacement = "")

# add rownames to the tpm and metadata dataframes
rownames(my_tpm) <- my_tpm[,1] # add the rownames
my_tpm <- my_tpm[,-1] # Remove the old column of rownames
rownames(my_metadata) <- my_metadata[,1] # add the rownames
# my_metadata <- my_metadata[,-1] # Remove the old column of rownames



###########################################################
################ IMPORT SEPT_Seq TPM VALUES ###############

Sept_tpm <- read.csv("/Users/elamont/Documents/RProjects/Sputum/ProbeTest3/Mtb.Expression.Gene.Data.SCALED.TPM.csv")

Sept_tpm <- Sept_tpm[,-ncol(Sept_tpm)] # remove the last column which is the Undetermined

# Adjust the names so they are slightly shorter
# names(my_tpm) <- gsub(x = names(my_tpm), pattern = "_2_2_ng_mL", replacement = "")
names(Sept_tpm) <- gsub(x = names(Sept_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)




# BELOW NOT DONE AS OF 12/10/24


###########################################################
############## IMPORT BOBS DE DATA (SPUTUM) ############### 

W0.S250754_ComparedTo_INDIGO.Broth <- read.delim("JOINED_BobAverages/MTb.MetaResults.S250754_vs_broth/S_250754_Probe_4A_50_S24.MTb.Meta.JOINED.txt")
W0.S250754_ComparedTo_W2.S349941 <- read.delim("JOINED_BobAverages/MTb.MetaResults.S250754_vs_S349941/sputum_250754_W0.MTb.Meta.JOINED.txt")
W0.S250754_ComparedTo_W0.S503557 <- read.delim("JOINED_BobAverages/MTb.MetaResults.S250754_vs_S503557/sputum_250754_W0.MTb.Meta.JOINED.txt")
W0.S250754_ComparedTo_W4.S687338 <- read.delim("JOINED_BobAverages/MTb.MetaResults.S250754_vs_S687338/S_250754_Probe_4A_50_S24.MTb.Meta.JOINED.txt")
W2.S349941_ComparedTo_W2.S575533 <- read.delim("JOINED_BobAverages/MTb.MetaResults.S349941_vs_S575533/S_349941_Probe_3D_25_S21.MTb.Meta.JOINED.txt")

###########################################################
################ MAKE A LIST OF ALL DFs ###################

list_dfs <- list(W0.S250754_ComparedTo_INDIGO.Broth,
                 W0.S250754_ComparedTo_W2.S349941, 
                 W0.S250754_ComparedTo_W0.S503557, 
                 W0.S250754_ComparedTo_W4.S687338, 
                 W2.S349941_ComparedTo_W2.S575533)

# Make a list of all the names
df_names <- c("W0.S250754_ComparedTo_INDIGO.Broth",
              "W0.S250754_ComparedTo_W2.S349941",
              "W0.S250754_ComparedTo_W0.S503557", 
              "W0.S250754_ComparedTo_W4.S687338", 
              "W2.S349941_ComparedTo_W2.S575533")

# Give the df list the correct df names
names(list_dfs) <- df_names

###########################################################
############### ADD COLUMNS OF DE VALUES ##################

# Make a new list to hold dataframes with extra columns
list_dfs_2 <- list()

ordered_DE <- c("significant down", "not significant", "significant up")

# Add extra DE columns to each dataframe
for (i in 1:length(list_dfs)) {
  
  current_df <- list_dfs[[i]]
  current_df_name <- df_names[i]
  
  # Make the column pointing out which ones are differential expressed
  current_df$DE <- ifelse(current_df$LOG2FOLD < -1 & current_df$AVG_PVALUE < 0.05, "significant down",
                          ifelse(current_df$LOG2FOLD > 1 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE <- factor(current_df$DE, levels = ordered_DE)
  
  # Make the column with DE gene names for plotting on graph
  current_df$DE_labels <- ifelse(current_df$DE != "not significant", current_df$GENE_NAME, NA)
  
  list_dfs_2[[current_df_name]] <- current_df
  
}
  
# View(list_dfs_2[[1]])



###########################################################
################### IMPORT MARK's DATA ####################
# 11/26/24

# Import the tpm
mark_tpm <- read.csv("DataFromMark/Data/Updated.Mtb.Expression.Gene.Data.TPM.csv")
# add rownames to the tpm and metadata dataframes
rownames(mark_tpm) <- mark_tpm[,1] # add the rownames
mark_tpm <- mark_tpm[,-1] # Remove the old column of rownames

# Import the metatdata
mark_metadata <- read.csv("DataFromMark/Data/Annotation_including_new_timecourse_batch.csv")

# I think I just need the GroupStrain the Mimic_timepoint and Batch
# GroupStrain should be the same as Sample_Type
colnames(mark_metadata)[colnames(mark_metadata) == "GroupStrain"] ="Sample_Type"
mark_metadata_2 <- mark_metadata %>% select(SampleID, Sample_Type, Batch)

# Merge the TPMs
combined_TPM <- merge(my_tpm, mark_tpm, by = "row.names", all = T)
rownames(combined_TPM) <- combined_TPM[,1] # add the rownames
combined_TPM <- combined_TPM[,-1] # Remove the old column of rownames

# Merge the metadatas
combined_metadata <- merge(my_metadata, mark_metadata_2, all = T)


