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
c14 <- c("dodgerblue2", "#E31A1C", "green4", "#FF7F00", "black","gold1", "#CAB2D6", "palegreen2", "gray70", "maroon", "orchid1", "blue1", "darkturquoise", "darkorange4") 
c3 <- c("#56B4E9", "#E66900", "#009E73")

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

my_pipeSummary <- my_pipeSummary %>% mutate(Sputum_Number = str_extract(SampleID, "S_[0-9]+"))


###########################################################
############ IMPORT SEPT PIPELINE SUMMARY DATA ############

Sept_pipeSummary <- read.csv("ProbeTest3_Pipeline.Summary.Details.csv")

Sept_pipeSummary$Week <- as.character(Sept_pipeSummary$Week)
ordered_Week <- c("0", "2", "4")
Sept_pipeSummary$Week <- factor(Sept_pipeSummary$Week, levels = ordered_Week)

Sept_pipeSummary <- Sept_pipeSummary[-nrow(Sept_pipeSummary), ] # Remove the undetermined, which is the last row

Sept_pipeSummary$SampleID <- gsub(x = Sept_pipeSummary$SampleID, pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)


# Combine the sputum samples only 
AllSputum_pipeSummary <- merge(my_pipeSummary %>% filter(Sample_Type == "Sputum") %>% 
                                 mutate(SeqRun = "Nov"), 
                               Sept_pipeSummary %>% filter(Sample_Type == "Sputum") %>% 
                                 mutate(SeqRun = "Sept"),
                               all = T) %>%
  mutate(Sputum_Number = str_extract(SampleID, "S_[0-9]+")) # Regular expression (regex)

# Add a column with Hex colors, Did this manually to make sure same samples were same color
AllSputum_pipeSummary$Colours <- c("dodgerblue2", "#E31A1C", "green4", "#FF7F00", "black",
                                   "darkorange4", "dodgerblue2", "gold1",  "#CAB2D6", "palegreen2", 
                                   "gray70", "#FF7F00", "#FF7F00", "maroon", "orchid1", 
                                   "black", "black", "blue1", "darkturquoise", "darkorange4", 
                                   "darkorange4") 


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
############## IMPORT NOT SCALED TPM VALUES ###############

my_tpm_NotScaled <- read.csv("Mtb.Expression.Gene.Data.TPM.csv")

my_tpm_NotScaled <- my_tpm_NotScaled[,-ncol(my_tpm_NotScaled)] # remove the last column which is the Undetermined

# Adjust the names so they are slightly shorter
names(my_tpm_NotScaled) <- gsub(x = names(my_tpm_NotScaled), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# add rownames to the tpm and metadata dataframes
rownames(my_tpm_NotScaled) <- my_tpm_NotScaled[,1] # add the rownames
my_tpm_NotScaled <- my_tpm_NotScaled[,-1] # Remove the old column of rownames

###########################################################
################ IMPORT SEPT_Seq TPM VALUES ###############

# Sept_tpm <- read.csv("/Users/elamont/Documents/RProjects/Sputum/ProbeTest3/Mtb.Expression.Gene.Data.SCALED.TPM.csv")
# Sept_tpm <- read.csv("/Users/snork-maiden/Documents/Micro_grad_school/Sherman_Lab/R_projects/Sputum/ProbeTest3/Mtb.Expression.Gene.Data.SCALED.TPM.csv")

# A little complicated because I am working across two computers
possible_paths <- c(
  "/Users/elamont/Documents/RProjects/Sputum/ProbeTest3/Mtb.Expression.Gene.Data.SCALED.TPM.csv",
  "/Users/snork-maiden/Documents/Micro_grad_school/Sherman_Lab/R_projects/Sputum/ProbeTest3/Mtb.Expression.Gene.Data.SCALED.TPM.csv"
)
# Find the first valid path
file_path <- possible_paths[file.exists(possible_paths)][1]
if (!is.null(file_path)) {
  Sept_tpm <- read.csv(file_path)
} else {
  stop("File not found in any of the expected locations.")
}

Sept_tpm <- Sept_tpm[,-ncol(Sept_tpm)] # remove the last column which is the Undetermined

# Adjust the names so they are slightly shorter
# names(my_tpm) <- gsub(x = names(my_tpm), pattern = "_2_2_ng_mL", replacement = "")
names(Sept_tpm) <- gsub(x = names(Sept_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# add rownames to the tpm and metadata dataframes
rownames(Sept_tpm) <- Sept_tpm[,1] # add the rownames
Sept_tpm <- Sept_tpm[,-1] # Remove the old column of rownames


# Combine the sputum samples only 
AllSputum_tpm <- merge(my_tpm %>% select(starts_with("S_")),
                       Sept_tpm %>% select(starts_with("S_")),
                       by = "row.names")
rownames(AllSputum_tpm) <- AllSputum_tpm$Row.names
AllSputum_tpm <- AllSputum_tpm %>% select(-Row.names)


# Grab AllSputum metadata
AllSputum_metadata <- AllSputum_pipeSummary # %>% select(2, 14:30)
# Adjust the metadata names so they are the same
AllSputum_metadata$SampleID <- sub(x = AllSputum_metadata$SampleID, pattern = "_S.*", replacement = "")
rownames(AllSputum_metadata) <- AllSputum_metadata[,1] # add the rownames


###########################################################
############ IMPORT SEPT NOT SCALED TPM VALUES ############

Sept_tpm_NotScaled <- read.csv("ProbeTest3_Mtb.Expression.Gene.Data.TPM.csv")

Sept_tpm_NotScaled <- Sept_tpm_NotScaled[,-ncol(Sept_tpm_NotScaled)] # remove the last column which is the Undetermined

# Adjust the names so they are slightly shorter
names(Sept_tpm_NotScaled) <- gsub(x = names(Sept_tpm_NotScaled), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# add rownames to the tpm and metadata dataframes
rownames(Sept_tpm_NotScaled) <- Sept_tpm_NotScaled[,1] # add the rownames
Sept_tpm_NotScaled <- Sept_tpm_NotScaled[,-1] # Remove the old column of rownames





###########################################################
################## IMPORT BOBS DE DATA #################### 

# Naming these by hand which is a pain but not sure how else to do it to make sure I am keeping track of which one is being compared to which

newTHP1_1e5_1_DualrRNA_ComparedTo_newTHP1_1e5_2_DualrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.newTHP1_1e5_1_DualrRNA_vs_newTHP1_1e5_2_DualrRNA/newTHP1_1e5_1_DualrRNA_S13.MTb.Meta.JOINED.txt")
newTHP1_1e6_1_DualrRNA_ComparedTo_newTHP1_1e6_2_DualrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.newTHP1_1e6_1_DualrRNA_vs_newTHP1_1e6_2_DualrRNA/newTHP1_1e6_1_DualrRNA_S15.MTb.Meta.JOINED.txt")

orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_2_DualrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.newTHP1_1e6_1_DualrRNA_vs_newTHP1_1e6_2_DualrRNA/newTHP1_1e6_1_DualrRNA_S15.MTb.Meta.JOINED.txt")
orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_3_DualrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.orginalTHP1_1e6_1_DualrRNA_vs_originalTHP1_1e6_3_DualrRNA/orginalTHP1_1e6_1_DualrRNA_S22.MTb.Meta.JOINED.txt")
orginalTHP1_1e6_2_DualrRNA_ComparedTo_originalTHP1_1e6_3_DualrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.originalTHP1_1e6_2_DualrRNA_vs_originalTHP1_1e6_3_DualrRNA/originalTHP1_1e6_2_DualrRNA_S23.MTb.Meta.JOINED.txt")
orginalTHP1_1e6_4_DualrRNA_ComparedTo_originalTHP1_1e6_5_DualrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.originalTHP1_1e6_4_DualrRNA_vs_originalTHP1_1e6_5_DualrRNA/originalTHP1_1e6_4_DualrRNA_S25.MTb.Meta.JOINED.txt")
orginalTHP1_1e6_4_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.originalTHP1_1e6_6_DualrRNA_vs_originalTHP1_1e6_4_DualrRNA/originalTHP1_1e6_4_DualrRNA_S25.MTb.Meta.JOINED.txt")
orginalTHP1_1e6_5_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.originalTHP1_1e6_6_DualrRNA_vs_originalTHP1_1e6_5_DualrRNA/originalTHP1_1e6_5_DualrRNA_S26.MTb.Meta.JOINED.txt")
orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.orginalTHP1_1e6_1_DualrRNA_vs_originalTHP1_1e6_6_DualrRNA/orginalTHP1_1e6_1_DualrRNA_S22.MTb.Meta.JOINED.txt")

orginalTHP1_1e6_1_MtbrRNA_ComparedTo_orginalTHP1_1e6_2_MtbrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.orginalTHP1_1e6_1_MtbrRNA_vs_orginalTHP1_1e6_2_MtbrRNA/orginalTHP1_1e6_1_MtbrRNA_S28.MTb.Meta.JOINED.txt")
orginalTHP1_1e6_1_MtbrRNA_ComparedTo_orginalTHP1_1e6_3_MtbrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.orginalTHP1_1e6_1_MtbrRNA_vs_orginalTHP1_1e6_3_MtbrRNA/orginalTHP1_1e6_1_MtbrRNA_S28.MTb.Meta.JOINED.txt")
orginalTHP1_1e6_2_MtbrRNA_ComparedTo_orginalTHP1_1e6_3_MtbrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.orginalTHP1_1e6_2_MtbrRNA_vs_orginalTHP1_1e6_3_MtbrRNA/orginalTHP1_1e6_2_MtbrRNA_S29.MTb.Meta.JOINED.txt")

W0.S354850_DualrRNA_ComparedTo_W0.S354851_DualrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.S_354850_DualrRNA_vs_S_354851_DualrRNA/S_354850_DualrRNA_S4.MTb.Meta.JOINED.txt")
W0.S503557_DualrRNA_ComparedTo_W0.S503557 <- read.delim("JOINED_BobAverages/MTb.MetaResults.S_503557_vs_S_503557_DualrRNA/S_503557_DualrRNA_S1.MTb.Meta.JOINED.txt")
W2.S575533_DualrRNA_ComparedTo_W2.S575533_MtbrRNA <- read.delim("JOINED_BobAverages/MTb.MetaResults.S_575533_DualrRNA_vs_S_575533_MtbrRNA/S_575533_DualrRNA_S2.MTb.Meta.JOINED.txt")

THP1_1e6_1_ComparedTo_THP1_1e6_2 <- read.delim("JOINED_BobAverages/MTb.MetaResults.THP1_1e6_1_vs_THP1_1e6_2/THP1_1e6_1_S41.MTb.Meta.JOINED.txt")
THP1_1e6_1_ComparedTo_THP1_1e6_3 <- read.delim("JOINED_BobAverages/MTb.MetaResults.THP1_1e6_1_vs_THP1_1e6_3/THP1_1e6_1_S41.MTb.Meta.JOINED.txt")
THP1_1e6_1_ComparedTo_THP1_1e6_4 <- read.delim("JOINED_BobAverages/MTb.MetaResults.THP1_1e6_1_vs_THP1_1e6_4/THP1_1e6_1_S41.MTb.Meta.JOINED.txt")
THP1_1e6_1_ComparedTo_THP1_1e6_5 <- read.delim("JOINED_BobAverages/MTb.MetaResults.THP1_1e6_1_vs_THP1_1e6_5/THP1_1e6_1_S41.MTb.Meta.JOINED.txt")
THP1_1e6_2_ComparedTo_THP1_1e6_3 <- read.delim("JOINED_BobAverages/MTb.MetaResults.THP1_1e6_2_vs_THP1_1e6_3/THP1_1e6_2_S42.MTb.Meta.JOINED.txt")
THP1_1e6_2_ComparedTo_THP1_1e6_4 <- read.delim("JOINED_BobAverages/MTb.MetaResults.THP1_1e6_2_vs_THP1_1e6_4/THP1_1e6_2_S42.MTb.Meta.JOINED.txt")
THP1_1e6_2_ComparedTo_THP1_1e6_5 <- read.delim("JOINED_BobAverages/MTb.MetaResults.THP1_1e6_2_vs_THP1_1e6_5/THP1_1e6_2_S42.MTb.Meta.JOINED.txt")
THP1_1e6_3_ComparedTo_THP1_1e6_4 <- read.delim("JOINED_BobAverages/MTb.MetaResults.THP1_1e6_3_vs_THP1_1e6_4/THP1_1e6_3_S43.MTb.Meta.JOINED.txt")
THP1_1e6_3_ComparedTo_THP1_1e6_5 <- read.delim("JOINED_BobAverages/MTb.MetaResults.THP1_1e6_3_vs_THP1_1e6_5/THP1_1e6_3_S43.MTb.Meta.JOINED.txt")
THP1_1e6_4_ComparedTo_THP1_1e6_5 <- read.delim("JOINED_BobAverages/MTb.MetaResults.THP1_1e6_4_vs_THP1_1e6_5/THP1_1e6_4_S44.MTb.Meta.JOINED.txt")

orginalTHP1_1e6_3_MtbrRNA_ComparedTo_THP1_1e6_5 <- read.delim("JOINED_BobAverages/MTb.MetaResults.orginalTHP1_1e6_3_MtbrRNA_vs_THP1_1e6_5/orginalTHP1_1e6_3_MtbrRNA_S30.MTb.Meta.JOINED.txt")
orginalTHP1_1e6_1_MtbrRNA_ComparedTo_THP1_1e6_5 <- read.delim("JOINED_BobAverages/MTb.MetaResults.orginalTHP1_1e6_1_MtbrRNA_vs_THP1_1e6_5/orginalTHP1_1e6_1_MtbrRNA_S28.MTb.Meta.JOINED.txt")
orginalTHP1_1e6_2_MtbrRNA_ComparedTo_THP1_1e6_5 <- read.delim("JOINED_BobAverages/MTb.MetaResults.orginalTHP1_1e6_2_MtbrRNA_vs_THP1_1e6_5/orginalTHP1_1e6_2_MtbrRNA_S29.MTb.Meta.JOINED.txt")

# Grouped sputum comparison
W0_DualrRNADep_ComparedTo_W2_DualrRNADep <- read.delim("JOINED_BobAverages/MTb.MetaResults.W0_vs_W2_DualrRNADep/W0.MTb.Meta.JOINED.txt")

# Compare across This Nov run and the Sept run
W0.S250754_ComparedTo_W0.S250754_Probe_4A_50 <- read.delim("JOINED_BobAverages/DE_Across_Runs/MTb.MetaResults.S_250754_vs_S_250754_Probe_4A_50/S_250754_S47.MTb.Meta.JOINED.txt")
W0.S503557_ComparedTo_W0.S503557_Probe_3D_10 <- read.delim("JOINED_BobAverages/DE_Across_Runs/MTb.MetaResults.S_503557_vs_S_503557_Probe_3D_10/S_503557_S46.MTb.Meta.JOINED.txt")
W2.S575533_ComparedTo_W0.S575533_Probe_3A <- read.delim("JOINED_BobAverages/DE_Across_Runs/MTb.MetaResults.S_575533_MtbrRNA_vs_S_575533_Probe_3A/S_575533_MtbrRNA_S39.MTb.Meta.JOINED.txt")
THP1_1e6_5_ComparedTo_THP1_1e6_3_Probe_3D_25 <- read.delim("JOINED_BobAverages/DE_Across_Runs/MTb.MetaResults.THP1_1e6_5_vs_THP1_1e6_3_Probe_3D_25/THP1_1e6_5_S45.MTb.Meta.JOINED.txt")
THP1_1e6_3_ComparedTo_THP1_1e6_3_Probe_3D_25 <- read.delim("JOINED_BobAverages/MTb.MetaResults.THP1_1e6_3_vs_THP1_1e6_3_Probe_3D_25/THP1_1e6_3_S43.MTb.Meta.JOINED.txt")
# Grouped sputum comparison - All but unique patients
W0_ComparedTo_W2_AllUnique <- read.delim("JOINED_BobAverages/DE_Across_Runs/MTb.MetaResults.W0_vs_W2_AllUnique/W0.MTb.Meta.JOINED.txt")



###########################################################
################ MAKE A LIST OF ALL DFs ###################

list_dfs <- list(newTHP1_1e5_1_DualrRNA_ComparedTo_newTHP1_1e5_2_DualrRNA,
                 newTHP1_1e6_1_DualrRNA_ComparedTo_newTHP1_1e6_2_DualrRNA, #2
                 
                 orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_2_DualrRNA, 
                 orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_3_DualrRNA, 
                 orginalTHP1_1e6_2_DualrRNA_ComparedTo_originalTHP1_1e6_3_DualrRNA,
                 orginalTHP1_1e6_4_DualrRNA_ComparedTo_originalTHP1_1e6_5_DualrRNA,
                 orginalTHP1_1e6_4_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA,
                 orginalTHP1_1e6_5_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA, #8
                 orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA,
                 
                 orginalTHP1_1e6_1_MtbrRNA_ComparedTo_orginalTHP1_1e6_2_MtbrRNA,
                 orginalTHP1_1e6_1_MtbrRNA_ComparedTo_orginalTHP1_1e6_3_MtbrRNA,
                 orginalTHP1_1e6_2_MtbrRNA_ComparedTo_orginalTHP1_1e6_3_MtbrRNA,
                 
                 W0.S354850_DualrRNA_ComparedTo_W0.S354851_DualrRNA,
                 W0.S503557_DualrRNA_ComparedTo_W0.S503557,
                 W2.S575533_DualrRNA_ComparedTo_W2.S575533_MtbrRNA, #15
                 
                 THP1_1e6_1_ComparedTo_THP1_1e6_2,
                 THP1_1e6_1_ComparedTo_THP1_1e6_3,
                 THP1_1e6_1_ComparedTo_THP1_1e6_4,
                 THP1_1e6_1_ComparedTo_THP1_1e6_5,
                 THP1_1e6_2_ComparedTo_THP1_1e6_3,
                 THP1_1e6_2_ComparedTo_THP1_1e6_4,
                 THP1_1e6_2_ComparedTo_THP1_1e6_5,
                 THP1_1e6_3_ComparedTo_THP1_1e6_4,
                 THP1_1e6_3_ComparedTo_THP1_1e6_5,  
                 THP1_1e6_4_ComparedTo_THP1_1e6_5, 
                 
                 orginalTHP1_1e6_3_MtbrRNA_ComparedTo_THP1_1e6_5,
                 orginalTHP1_1e6_1_MtbrRNA_ComparedTo_THP1_1e6_5,
                 orginalTHP1_1e6_2_MtbrRNA_ComparedTo_THP1_1e6_5, 
                 
                 W0_DualrRNADep_ComparedTo_W2_DualrRNADep,
                 
                 W0.S250754_ComparedTo_W0.S250754_Probe_4A_50,
                 W0.S503557_ComparedTo_W0.S503557_Probe_3D_10,
                 W2.S575533_ComparedTo_W0.S575533_Probe_3A,
                 THP1_1e6_5_ComparedTo_THP1_1e6_3_Probe_3D_25,
                 THP1_1e6_3_ComparedTo_THP1_1e6_3_Probe_3D_25,
                 
                 W0_ComparedTo_W2_AllUnique) 

# Make a list of all the names
df_names <- c("newTHP1_1e5_1_DualrRNA_ComparedTo_newTHP1_1e5_2_DualrRNA",
              "newTHP1_1e6_1_DualrRNA_ComparedTo_newTHP1_1e6_2_DualrRNA", #2
              
              "orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_2_DualrRNA", 
              "orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_3_DualrRNA", 
              "orginalTHP1_1e6_2_DualrRNA_ComparedTo_originalTHP1_1e6_3_DualrRNA",
              "orginalTHP1_1e6_4_DualrRNA_ComparedTo_originalTHP1_1e6_5_DualrRNA",
              "orginalTHP1_1e6_4_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA",
              "orginalTHP1_1e6_5_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA", #8
              "orginalTHP1_1e6_1_DualrRNA_ComparedTo_originalTHP1_1e6_6_DualrRNA",
              
              "orginalTHP1_1e6_1_MtbrRNA_ComparedTo_orginalTHP1_1e6_2_MtbrRNA",
              "orginalTHP1_1e6_1_MtbrRNA_ComparedTo_orginalTHP1_1e6_3_MtbrRNA",
              "orginalTHP1_1e6_2_MtbrRNA_ComparedTo_orginalTHP1_1e6_3_MtbrRNA",
              
              "W0.S354850_DualrRNA_ComparedTo_W0.S354851_DualrRNA",
              "W0.S503557_DualrRNA_ComparedTo_W0.S503557",
              "W2.S575533_DualrRNA_ComparedTo_W2.S575533_MtbrRNA", #15
              
              "THP1_1e6_1_ComparedTo_THP1_1e6_2",
              "THP1_1e6_1_ComparedTo_THP1_1e6_3",
              "THP1_1e6_1_ComparedTo_THP1_1e6_4",
              "THP1_1e6_1_ComparedTo_THP1_1e6_5",
              "THP1_1e6_2_ComparedTo_THP1_1e6_3",
              "THP1_1e6_2_ComparedTo_THP1_1e6_4",
              "THP1_1e6_2_ComparedTo_THP1_1e6_5",
              "THP1_1e6_3_ComparedTo_THP1_1e6_4",
              "THP1_1e6_3_ComparedTo_THP1_1e6_5", #23
              "THP1_1e6_4_ComparedTo_THP1_1e6_5", #24
              
              "orginalTHP1_1e6_3_MtbrRNA_ComparedTo_THP1_1e6_5", # 25
              "orginalTHP1_1e6_1_MtbrRNA_ComparedTo_THP1_1e6_5",
              "orginalTHP1_1e6_2_MtbrRNA_ComparedTo_THP1_1e6_5",
              
              "W0_DualrRNADep_ComparedTo_W2_DualrRNADep",
              
              "W0.S250754_ComparedTo_W0.S250754_Probe_4A_50",
              "W0.S503557_ComparedTo_W0.S503557_Probe_3D_10",
              "W2.S575533_ComparedTo_W0.S575533_Probe_3A",
              "THP1_1e6_5_ComparedTo_THP1_1e6_3_Probe_3D_25", #31
              "THP1_1e6_3_ComparedTo_THP1_1e6_3_Probe_3D_25",
              
              "W0_ComparedTo_W2_AllUnique"
  
)

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




