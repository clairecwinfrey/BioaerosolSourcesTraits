# (Robust) ANCOMBC-2 for foliar surface versus bioaerosols
# Updated September and October 2025

# Description: Using un-rarefied foliar surface and bioaerosol bacterial and fungal datasets, this script performs ANCOMBC2 analyses to look for 
# differentially abundant taxa in foliar surface versus bioaerosol samples. The criteria for considering a differentially
# abundant taxon, whether a structural or "main" ANCOM zero, is that it must occur in at least 10% of samples and comprise at least 0.01% of reads 
# from the respective group. AND CRUCIALLY, ONLY GETS THOSE TAXA THAT MEET "ROBUST" THRESHOLD

# Results (10% samples and 0.05% of reads):
# 16S:
# — bioaerosol: 94 differentially abundant ASVs
# — foliar: 126 differentially abundant ASVs

# ITS: 
# — bioaerosol: 152 differentially abundant ASVs
# — foliar: 243 differentially abundant ASVs

##################
# I. SET-UP 
##################

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager", lib="/data/winfreyc/R/x86_64-pc-linux-gnu-library/4.3")
# BiocManager::install(version = "3.18", lib="/data/winfreyc/R/x86_64-pc-linux-gnu-library/4.3")
# 
# BiocManager::install("ANCOMBC") #downloaded 09-28-2025 on my computer
# Downloaded packages are here: /var/folders/8_/n60bcy2n3cj0q6bs08chjq280000gq/T//Rtmp5FCVvO/downloaded_packages
# BiocManager::install("microbiome") #downloaded 09-28-2025 on my computer

# 1. LIBRARIES
library(vegan)
library(phyloseq); packageVersion("phyloseq") #‘1.52.0’
library(ANCOMBC); packageVersion("ANCOMBC") #‘2.10.1’
library("microbiome")
library(tidyr)
library(tibble)
library(dplyr)
library(lmerTest)
library(stringr)

# 2. LOAD PHYLOSEQ OBJECTS (made in scripts on my computer called full16S_EDArarefied_Part2_Sept.R and equivalent fungal script)
##### RAREFIED #####
# Read in 16S data
# I6S_ALL_phyloseq <-  readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_noAirSingsDoubs_ps_sept25.rds")
# # Read in ITS data
# ITS_ALL_phyloseq <- readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/sallITSr_noAirSingsDoubs_ps.rds") #yes, this is correct version, even with "sallITS" typo!
# tax_table(ITS_ALL_phyloseq)

##### NON-RAREFIED #####
# Load unrarefied data for all samples (air singletons and doubletons removed and samples lower than 5.5K reads also removed)
I6Sall_5.5Kfiltered.ps <-  readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airI6S_5.5Kfiltered_phyloseq.rds")
# Load unrarefied data for all samples (air singletons and doubletons removed and samples lower than 8.5K reads alos removed)
ITSall_8.5Kfiltered.ps <-  readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airITS_8.5Kfiltered_phyloseq.rds")

########################################################################################################################
#####################################
# II BACTERIA:
# ANCOM FOR AIR VERSUS LEAF SURFACES/PHYLLOSPHERE
####################################
########## 1. FORMAT!! ##########
# subset no soil
unique(sample_data(I6Sall_5.5Kfiltered.ps)$sampleType)
sort(colSums(as.data.frame(as.matrix(otu_table(I6Sall_5.5Kfiltered.ps))))) #shows that all samples all vary (i.e. not rarefied)
min(colSums(as.data.frame(as.matrix(otu_table(I6Sall_5.5Kfiltered.ps))))) #5517
max(colSums(as.data.frame(as.matrix(otu_table(I6Sall_5.5Kfiltered.ps))))) #117945
I6Sall_5.5Kfiltered.psAIR_LEAF <- subset_samples(I6Sall_5.5Kfiltered.ps, sampleType != "soil")
unique(sample_data(I6Sall_5.5Kfiltered.psAIR_LEAF)$sampleType)
# remove ASVs that do not have taxa present
I6Sall_5.5Kfiltered.psAIR_LEAFZeros <- which(rowSums(otu_table(I6Sall_5.5Kfiltered.psAIR_LEAF))==0)
unique(rowSums(otu_table(I6Sall_5.5Kfiltered.psAIR_LEAF)[I6Sall_5.5Kfiltered.psAIR_LEAFZeros,])) #all zeros
I6Sall_5.5Kfiltered.psAIR_LEAFZerosNames <- names(I6Sall_5.5Kfiltered.psAIR_LEAFZeros)
# Confirm that this worked!
unique(colSums(otu_table(I6Sall_5.5Kfiltered.psAIR_LEAF)[rownames(otu_table(I6Sall_5.5Kfiltered.psAIR_LEAF)) %in% I6Sall_5.5Kfiltered.psAIR_LEAFZerosNames,]))
# Get what is unique in I6Sall_5.5Kfiltered.psAIR_LEAF that is not in the zero ones!
I6Sall_5.5Kfiltered.psAIR_LEAFASVsToKeep <- setdiff(rownames(otu_table(I6Sall_5.5Kfiltered.psAIR_LEAF)), I6Sall_5.5Kfiltered.psAIR_LEAFZerosNames) 
length(I6Sall_5.5Kfiltered.psAIR_LEAFASVsToKeep) #10,540
I6Sall_5.5Kfiltered.psAIR_LEAF.ps <- prune_taxa(taxa= I6Sall_5.5Kfiltered.psAIR_LEAFASVsToKeep, x=I6Sall_5.5Kfiltered.psAIR_LEAF) #9012 taxa and 169 samples 
which(rowSums(otu_table(I6Sall_5.5Kfiltered.psAIR_LEAF.ps)) == 0) #none, as expected

# Add "Species" at ASV level and make a new phyloseq
I6Sall_5.5Kfiltered.psAIR_LEAF_taxSpecies <- as.data.frame(as.matrix(tax_table(I6Sall_5.5Kfiltered.psAIR_LEAF.ps)))
class(I6Sall_5.5Kfiltered.psAIR_LEAF_taxSpecies) #great, it's not a phyloseq object!
head(I6Sall_5.5Kfiltered.psAIR_LEAF_taxSpecies)
I6S_airLeaf_ASVnames <- rownames(I6Sall_5.5Kfiltered.psAIR_LEAF_taxSpecies)
I6Sall_5.5Kfiltered.psAIR_LEAF_taxSpecies$Species <- rownames(I6Sall_5.5Kfiltered.psAIR_LEAF_taxSpecies)
head(I6Sall_5.5Kfiltered.psAIR_LEAF_taxSpecies)
# Turn this new tax table into a tax table for phyloseq
I6Sall_5.5Kfiltered.psAIR_LEAF_taxForPS <- phyloseq::tax_table(I6Sall_5.5Kfiltered.psAIR_LEAF_taxSpecies)
colnames(I6Sall_5.5Kfiltered.psAIR_LEAF_taxForPS) <- c("Kingdom","Phylum","Class", "Order","Family","Genus", "Species")
rownames(I6Sall_5.5Kfiltered.psAIR_LEAF_taxForPS) <- I6S_airLeaf_ASVnames #re-set rownames
#View(I6Sall_5.5Kfiltered.psAIR_LEAF_taxForPS)
head(I6Sall_5.5Kfiltered.psAIR_LEAF_taxForPS)
class(I6Sall_5.5Kfiltered.psAIR_LEAF_taxForPS)

# Finally, re-make phyloseq object
I6S_leafAir_notR_ANCOM.ps <- merge_phyloseq(I6Sall_5.5Kfiltered.psAIR_LEAF.ps, I6Sall_5.5Kfiltered.psAIR_LEAF_taxForPS)
dim(tax_table(I6S_leafAir_notR_ANCOM.ps)) #10540     7
tax_table(I6S_leafAir_notR_ANCOM.ps) #looks correct!
I6S_leafAir_ASVsTab <- as.data.frame(as.matrix(otu_table(I6S_leafAir_notR_ANCOM.ps)))
which(rowSums(I6S_leafAir_ASVsTab)<1) #none of these ASVs are NOT found in these groups! Great
dim(I6S_leafAir_ASVsTab) #10540  ASVs
unique(sample_data(I6S_leafAir_notR_ANCOM.ps)$sampleType) #"air"          "phyllosphere"
# For later, make and save an ASV table with only these ASVs
I6SairPhylloOnly.ps_ASVs <- t(as.data.frame(as.matrix(otu_table(I6S_leafAir_notR_ANCOM.ps))))
# Make sure that ASVs are columns
if (ncol(I6SairPhylloOnly.ps_ASVs)> nrow(I6SairPhylloOnly.ps_ASVs)){
  print("Good! ASVs are columns!")
} else {
  I6SairPhylloOnly.ps_ASVs <- t(I6SairPhylloOnly.ps_ASVs)
}
I6SairPhylloOnly_tax <- as.data.frame(as.matrix(tax_table(I6S_leafAir_notR_ANCOM.ps)))
I6SairPhylloOnly_tax <- I6SairPhylloOnly_tax %>% 
  tibble::rownames_to_column(var="ASV_name")
# Get names of air and phyllo samples to use later
I6SphylloSampNames <- rownames(sample_data(I6S_leafAir_notR_ANCOM.ps))[which(sample_data(I6S_leafAir_notR_ANCOM.ps)$sampleType == "phyllosphere")]
I6SairSampNames <- rownames(sample_data(I6S_leafAir_notR_ANCOM.ps))[which(sample_data(I6S_leafAir_notR_ANCOM.ps)$sampleType == "air")]

# Add in taxonomy
I6S_leafAir_notR_taxTab <- as.data.frame(as.matrix(tax_table(I6S_leafAir_notR_ANCOM.ps))) 
I6S_leafAir_notR_taxTab <- I6S_leafAir_notR_taxTab %>% 
  tibble::rownames_to_column(var="taxon")

# GET READ NUMBERS TO USE LATER
readNumberI6SairFoliar <- rowSums(I6SairPhylloOnly.ps_ASVs) #get read numbers from whole dataset
readNumberI6SairFoliar_df <- as.data.frame(matrix(nrow=length(readNumberI6SairFoliar), ncol=2))
colnames(readNumberI6SairFoliar_df) <- c("sampleName", "TotalReads")
readNumberI6SairFoliar_df[,1] <- names(readNumberI6SairFoliar)
readNumberI6SairFoliar_df[,2] <- unname(readNumberI6SairFoliar)

# September 30, 2025
# saveRDS(I6S_leafAir_notR_ANCOM.ps, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_leafAir_notR_ANCOM_phyloseq.rds")

########## 2. RUN ANCOM!!! ##########
# BACTERIA
I6S_leafAir_notR_ANCOM.ps <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_leafAir_notR_ANCOM_phyloseq.rds")
# Following lines make sure that the factor is done correctly
I6SsampDat_Check <- as.data.frame(sample_data(I6S_leafAir_notR_ANCOM.ps))
sort(colnames(I6SsampDat_Check))
I6SsampDat_Check$sampleType <- stats::relevel(factor(I6SsampDat_Check$sampleType), ref = "air") 
levels(I6SsampDat_Check$sampleType)[1] 
sample_data(I6S_leafAir_notR_ANCOM.ps)$sampleType <- I6SsampDat_Check$sampleType
sample_data(I6S_leafAir_notR_ANCOM.ps)$sampleType #looks good! Air is first level and so will be intercept reference (i.e., - for log fold change..)

# # Depending on version, inspect where LFCs live:
# names(res)                     # map the structure
# colnames(res$res$lfc)          # or: colnames(res$beta) / colnames(res$res$beta)

# RUN ANCOM!
set.seed(19)
# Grayed out since took a long time to run!
# I6S_leafAir_notR_ANCOM_ASVlvl <- ancombc2(data = I6S_leafAir_notR_ANCOM.ps, group = "sampleType", tax_level = "Species", fix_formula = "sampleType", struc_zero = TRUE, neg_lb = FALSE) #bacteria
# Saved/run October 1, 2025
# saveRDS(I6S_leafAir_notR_ANCOM_ASVlvl, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_leafAir_notR_ANCOM_ASVlvl.rds")
I6S_leafAir_notR_ANCOM_ASVlvl <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_leafAir_notR_ANCOM_ASVlvl.rds")
# View(I6S_leafAir_notR_ANCOM_ASVlvl$res)

########## 3. SET UP/INVESTIGATE MAIN (NON STRUCTURAL ZEROS OUTPUT) ##########
# TRY ONLY THOSE NOT SENSITIVE TO PSEUDOCOUNTS
head(I6S_leafAir_notR_ANCOM_ASVlvl$res)
length(which(I6S_leafAir_notR_ANCOM_ASVlvl$res$diff_robust_sampleTypephyllosphere == TRUE)) #91 are not affected by addition of pseudocounts

# Get a subset. Note that lfc are positive are phyllosphere relative to air
I6S_ANCOM_airLeaf_results_robust <- I6S_leafAir_notR_ANCOM_ASVlvl$res[which(I6S_leafAir_notR_ANCOM_ASVlvl$res$diff_robust_sampleTypephyllosphere == TRUE),]
# View(I6S_ANCOM_airLeaf_results_robust)

I6S_ANCOM_airLeaf_results_robust <- left_join(I6S_ANCOM_airLeaf_results_robust, I6S_leafAir_notR_taxTab, by= "taxon")

# Look at phyllosphere enriched taxa
length(which(I6S_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere > 0)) #76
unique(I6S_ANCOM_airLeaf_results_robust$Family[which(I6S_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere > 0)]) #greater than zero are phyllosphere
table(I6S_ANCOM_airLeaf_results_robust$Class[which(I6S_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere > 0)]) #only 1 Bacilli

# Look at bioaerosol-enriched taxa
length(which(I6S_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere < 0)) #15
unique(I6S_ANCOM_airLeaf_results_robust$Family[which(I6S_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere < 0)]) #less than zero are bioaerosol
table(I6S_ANCOM_airLeaf_results_robust$Class[which(I6S_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere < 0)]) #3 Bacilli!

# Saved October 9, 2025:
# saveRDS(I6S_ANCOM_airLeaf_results_robust, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_ANCOM_airLeaf_results_robust_10-9-25.rds")

# Make dataframe smaller
colnames(I6S_ANCOM_airLeaf_results_robust)
# STOP #
colsToKeepANCOM_robust <- c("taxon", "lfc_sampleTypephyllosphere", "W_sampleTypephyllosphere", "q_sampleTypephyllosphere", "diff_robust_sampleTypephyllosphere", "Phylum", "Class",
                     "Order","Family","Genus","Species")
I6S_leafAir_ANCOM_small_robust <- I6S_ANCOM_airLeaf_results_robust[,colnames(I6S_ANCOM_airLeaf_results_robust) %in% colsToKeepANCOM_robust]
colnames(I6S_leafAir_ANCOM_small_robust)[colnames(I6S_leafAir_ANCOM_small_robust) %in% "taxon"] <- "ASV_name" #change to ASV name for merging
head(I6S_leafAir_ANCOM_small_robust)
# Add in category for bioaerosol or foliar surface
I6S_leafAir_ANCOM_small_robust$ANCOMcat[which(I6S_leafAir_ANCOM_small_robust$lfc_sampleTypephyllosphere > 0)] <- "foliar surface" #if lfc greater than 0, then foliar surface
I6S_leafAir_ANCOM_small_robust$ANCOMcat[which(I6S_leafAir_ANCOM_small_robust$lfc_sampleTypephyllosphere < 0)] <- "bioaerosol" #if lfc less than 0, then bioaerosol
# Double check
sort(I6S_leafAir_ANCOM_small_robust$lfc_sampleTypephyllosphere[which(I6S_leafAir_ANCOM_small_robust$ANCOMcat == "bioaerosol")]) #great, all are negative!
sort(I6S_leafAir_ANCOM_small_robust$lfc_sampleTypephyllosphere[which(I6S_leafAir_ANCOM_small_robust$ANCOMcat == "foliar surface")]) #great, all are positive!

foliarI6SANCOMASVsName_robust <- I6S_leafAir_ANCOM_small_robust$ASV_name[I6S_leafAir_ANCOM_small_robust$ANCOMcat == "foliar surface"]
bioaerosolI6SANCOMASVsName_robust <- I6S_leafAir_ANCOM_small_robust$ASV_name[I6S_leafAir_ANCOM_small_robust$ANCOMcat == "bioaerosol"]

# HOW DO THESE COMPARE TO ORIGINAL ANCOM?? Captured 40/49 OG phyllo, but only 8/37 OG air
OG_ANCOM_I6S <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/importedFromServer/I6S_MEAN_airLeafANCOM.Rdata")
head(OG_ANCOM_I6S)
OG_I6SairASVs <- OG_ANCOM_I6S$ASV_name[which(OG_ANCOM_I6S$ANCOMcat == "air")] #37
OG_I6SphylloASVs <- OG_ANCOM_I6S$ASV_name[which(OG_ANCOM_I6S$ANCOMcat == "phyllo")] #49
# Foliar
length(which(foliarI6SANCOMASVsName_robust %in% OG_I6SphylloASVs ==TRUE)) #40 are from last time!
length(intersect(OG_I6SphylloASVs,foliarI6SANCOMASVsName_robust)) #40 are from last time!
# Bioaerosols
length(which(bioaerosolI6SANCOMASVsName_robust %in% OG_I6SairASVs ==TRUE)) #7 are from last time!
length(intersect(OG_I6SairASVs,bioaerosolI6SANCOMASVsName_robust)) #7 are from last time!

########## 3. INVESTIGATE STRUCTURAL ZEROS ##########
# UNDERSTAND THIS DATA FRAME. Overview of the findings below. All $res identified taxa are also in this dataframe, but 
# they are not structural zeros in either group. Structural zeros for phyllo are those detected only in air but not phyllo, and vice versa. However, these
# can vary WIDELY with their abundance, so important to not take all of these as diff. abundant.
I6S_airLeafANCOMtabzero <- I6S_leafAir_notR_ANCOM_ASVlvl$zero_ind
# View(I6S_airLeafANCOMtabzero)
# 1. Is there any overlap between these taxa and those in the main results? YES, ALL MAIN RESULTS PLUS MANY MORE
length(intersect(I6S_ANCOM_airLeaf_results_robust$taxon,I6S_airLeafANCOMtabzero$taxon)) #91 overlap
unique(I6S_ANCOM_airLeaf_results_robust$taxon %in% I6S_airLeafANCOMtabzero$taxon) #all of the original ones are in the structural zeros
length(I6S_airLeafANCOMtabzero$taxon) #I6S_airLeafANCOMtabzero 
# 2. What do the main results taxa look like in the zero results? ALL FALSE, IN OTHER WORDS, THEY ARE NOT STRUCTURAL ZEROS, MAKES SENSE
unique(I6S_airLeafANCOMtabzero$`structural_zero (sampleType = air)`[I6S_airLeafANCOMtabzero$taxon %in% I6S_ANCOM_airLeaf_results_robust$taxon]) #all FALSE
unique(I6S_airLeafANCOMtabzero$`structural_zero (sampleType = phyllosphere)`[I6S_airLeafANCOMtabzero$taxon %in% I6S_ANCOM_airLeaf_results_robust$taxon]) #all FALSE
# 3. TRUE structural_zero (sampleType = air)`, what are these like? All are F for phyllo, makes sense
I6S_airLeafTab0_air <- I6S_airLeafANCOMtabzero[which(I6S_airLeafANCOMtabzero$`structural_zero (sampleType = air)` == TRUE),]
unique(I6S_airLeafTab0_air$`structural_zero (sampleType = phyllosphere)`) #FALSE
# 4. TRUE structural_zero (sampleType = phyllosphere)`, what are these like? All are F for air, makes sense
I6S_airLeafTab0_foliar <- I6S_airLeafANCOMtabzero[which(I6S_airLeafANCOMtabzero$`structural_zero (sampleType = phyllosphere)` == TRUE),]
unique(I6S_airLeafTab0_foliar$`structural_zero (sampleType = air)`) #all FALSE, makes sense
# 5. Are all of these air structural zeros not found in the air at all?
I6Sair0ASV_names_robust <- unique(I6S_airLeafTab0_air$taxon)
# Look for these in air samples -- correct, they are never found
unique(colSums(I6SairPhylloOnly.ps_ASVs[rownames(I6SairPhylloOnly.ps_ASVs) %in% I6SairSampNames,colnames(I6SairPhylloOnly.ps_ASVs) %in% I6Sair0ASV_names_robust]))
# 6. Are all of these air structural zeros found in the phyllosphere, as expected?
# Look for these in phyllo samples -- correct, are found many times! But, in some cases, these are in fact really rare, which is an issue because I don't 
# really want to use these as "differentially abundant" because they make no sense to analyze in a trait-based way.
sort(unique(colSums(I6SairPhylloOnly.ps_ASVs[rownames(I6SairPhylloOnly.ps_ASVs) %in% I6SphylloSampNames,colnames(I6SairPhylloOnly.ps_ASVs) %in% I6Sair0ASV_names_robust])))
# 7. Are all of these phyllo structural zeros not found in the phyllo at all?
I6Sphyllo0ASV_names <- unique(I6S_airLeafTab0_foliar$taxon)
# Look for these in phyllo samples -- correct, they are never found
unique(colSums(I6SairPhylloOnly.ps_ASVs[rownames(I6SairPhylloOnly.ps_ASVs) %in% I6SphylloSampNames,colnames(I6SairPhylloOnly.ps_ASVs) %in% I6Sphyllo0ASV_names]))
# 8. Are all of these phyllo structural zeros found in the air, as expected?
# Look for these in phyllo samples -- correct, are found many times! But, in some cases, these are in fact really rare, which is an issue because I don't 
# really want to use these as "differentially abundant" because they make no sense to analyze in a trait-based way. But then some are super abundant! Yikes, this is the issue!
sort(unique(colSums(I6SairPhylloOnly.ps_ASVs[rownames(I6SairPhylloOnly.ps_ASVs) %in% I6SairSampNames,colnames(I6SairPhylloOnly.ps_ASVs) %in% I6Sphyllo0ASV_names])))

# Given all of the above, get a zero data frame that just has the structural zeros
# get an index for all phyllo or air struct zeros
szIndexI6S <- c(which(I6S_airLeafANCOMtabzero$`structural_zero (sampleType = air)`== TRUE), which(I6S_airLeafANCOMtabzero$`structural_zero (sampleType = phyllosphere)`== TRUE))
# Get this dataframe and then clean it up some
I6S_AFS_zeros <- I6S_airLeafANCOMtabzero[szIndexI6S,]
# Get air and phyllosphere structural zeros
airZeros16S <- I6S_AFS_zeros$taxon[which(I6S_AFS_zeros$`structural_zero (sampleType = air)` == TRUE)]
length(airZeros16S) #5967 in phyllo but not air
phylloZeros16S <- I6S_AFS_zeros$taxon[which(I6S_AFS_zeros$`structural_zero (sampleType = phyllosphere)` == TRUE)]
length(phylloZeros16S) #3473 in air but not phyllo
nrow(I6S_AFS_zeros) == length(unique(I6S_AFS_zeros$taxon))
# Make a cleaner df for this
I6S_AFS_zeros_df <- as.data.frame(matrix(ncol=2, nrow=length(unique(I6S_AFS_zeros$taxon))))
colnames(I6S_AFS_zeros_df) <- c("ASV_name", "struct_ZeroType")
I6S_AFS_zeros_df$ASV_name <- I6S_AFS_zeros$taxon #add in names
I6S_AFS_zeros_df <- I6S_AFS_zeros_df %>%
  mutate( #if ASV name in air zeros, struct tye is bioaerosol, is in phyllo, phyllo
    struct_ZeroType = case_when(
      ASV_name %in% airZeros16S ~ "bioaerosol",
      ASV_name %in% phylloZeros16S ~ "foliar surfaces",
      TRUE ~ NA_character_
    )
  )
unique(sort(I6S_AFS_zeros_df$ASV_name[which(I6S_AFS_zeros_df$struct_ZeroType == "bioaerosol")]) == sort(airZeros16S))  #double check

# HOW DO THESE COMPARE TO ORIGINAL ANCOM ANALYSIS? 6/49 phyllo were in these structrual, 27/37 air!
OG_I6SairASVs <- OG_ANCOM_I6S$ASV_name[which(OG_ANCOM_I6S$ANCOMcat == "air")]; length(OG_I6SairASVs) #37
OG_I6SphylloASVs <- OG_ANCOM_I6S$ASV_name[which(OG_ANCOM_I6S$ANCOMcat == "phyllo")]; length(OG_I6SphylloASVs) #49
# Foliar
length(which(airZeros16S %in% OG_I6SphylloASVs ==TRUE)) #6 are from last time!
length(intersect(OG_I6SphylloASVs,airZeros16S)) #6, just another way of checking this
# Bioaerosols
length(which(phylloZeros16S %in% OG_I6SairASVs ==TRUE)) #27 are from last time!
length(intersect(OG_I6SairASVs,phylloZeros16S)) #still 27, so another way of checking this

########## 4. INVESTIGATE HOW ANCOM MODEL AND ACTUAL DATA COMPARE ##########
##### i. JUST MAIN $RES RESULTS #####
# 1. Create ASV table with only ANCOM ASVs in air and leaf samples. Can use this to compare how ANCOM model and actual data compare
# Double check all are robust and they are! (0 should be what I get below)
which(I6S_leafAir_ANCOM_small_robust$diff_robust_sampleTypephyllosphere == FALSE)
I6S_leafAir_ANCOM_ASVs_robust <- I6S_leafAir_ANCOM_small_robust$ASV_name
length(I6S_leafAir_ANCOM_ASVs_robust) #91 different ASVs, as above
# Pull out only ASV table for significantly different ASVs (ASVs are columns!)
I6S_airLeaf_sigASVtab_robust <- I6SairPhylloOnly.ps_ASVs[,colnames(I6SairPhylloOnly.ps_ASVs) %in% I6S_leafAir_ANCOM_ASVs_robust]
# View(I6S_airLeaf_sigASVtab_robust)
class(head(I6S_airLeaf_sigASVtab_robust))
# But which samples are air versus leaves?
#View(I6S_airLeaf_sigASVtab_robust)
I6S_airLeaf_sigASVtab_robust <- as.data.frame(I6S_airLeaf_sigASVtab_robust) %>% 
  tibble::rownames_to_column(var= "sampleName")
dim(I6S_airLeaf_sigASVtab_robust) #rows are samples, columns are ASVs
# View(I6S_airLeaf_sigASVtab_robust)

# Add in air or leaf, depending on what sample type it is
I6S_airLeaf_sigASVtab_robust <- I6S_airLeaf_sigASVtab_robust %>%
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "air",
      str_detect(sampleName, "phyllo") ~ "phyllo",
      TRUE ~ NA_character_
    )
  )
head(I6S_airLeaf_sigASVtab_robust)
cbind(I6S_airLeaf_sigASVtab_robust$sampleName, I6S_airLeaf_sigASVtab_robust$sampleType) #quick check to make sure code above worked... and it did!
# One last check before calculating...
all.equal(I6S_leafAir_ANCOM_ASVs_robust, colnames(I6S_airLeaf_sigASVtab_robust)[2:(ncol(I6S_airLeaf_sigASVtab_robust)-1)]) #okay these are still the significant ones

# 2. Add in number of reads per sample
#Add to working dataframe
head(I6S_airLeaf_sigASVtab_robust)
I6S_airLeaf_sigASVs_robust <- merge(I6S_airLeaf_sigASVtab_robust, readNumberI6SairFoliar_df, by = "sampleName")
# View(I6S_airLeaf_sigASVs_robust)

# 3. MAKE A "LONGER" DATAFRAME AND ADD IN TAXONOMY
head(I6S_airLeaf_sigASVs_robust) 
# First, pivot longer so that ASVs are column 1 and Sample is second column 
colnames(I6S_airLeaf_sigASVs_robust)
I6S_airLeaf_sigASVs_long_robust <- I6S_airLeaf_sigASVs_robust %>% 
  pivot_longer(cols = ASV_4:ASV_5328, names_to= "ASV_name", values_to = "ASV_abundance")
head(I6S_airLeaf_sigASVs_long_robust)
# View(I6S_airLeaf_sigASVs_long_robust)
dim(I6S_airLeaf_sigASVs_long_robust) #12922      5
unique(I6S_airLeaf_sigASVs_long_robust$sampleType) #great, only air and phyllo!
colnames(I6S_airLeaf_sigASVs_long_robust)
# Result is "long" object with the abundance of each ANCOM ASV in each air and phyllo sample, plus taxonomy!
I6S_airLeaf_ANCOM_long_robust <- left_join(I6S_airLeaf_sigASVs_long_robust, I6SairPhylloOnly_tax, by = "ASV_name")
head(I6S_airLeaf_ANCOM_long_robust)
# View(I6S_airLeaf_ANCOM_long_robust)
# Add in ANCOM cat of each 
colnames(I6S_leafAir_ANCOM_small_robust)
I6S_airLeaf_ANCOM_long_robust <- merge(I6S_airLeaf_ANCOM_long_robust, I6S_leafAir_ANCOM_small_robust[,colnames(I6S_leafAir_ANCOM_small_robust) %in% c("ASV_name", "ANCOMcat")], by ="ASV_name")

# 4. Get mean prop of each ASV in bioaerosol and foliar
I6S_airLeaf_ANCOM_long_robust <- I6S_airLeaf_ANCOM_long_robust %>% 
  mutate(propASVsamp = ASV_abundance/TotalReads) #get each ASVs proportion in each sample. Double checked that this worked!
# View(I6S_airLeaf_ANCOM_long_robust)
I6S_airLeaf_ANCOM_long2_robust <- I6S_airLeaf_ANCOM_long_robust %>% 
  group_by(ASV_name, ANCOMcat, sampleType) %>% 
  summarize(meanASVPctSampType = mean(propASVsamp)*100)
# View(I6S_airLeaf_ANCOM_long2_robust)

# Make this wider for direct comparison
colnames(I6S_airLeaf_ANCOM_long2_robust)
I6S_airLeaf_ANCOM_wide_robust <- I6S_airLeaf_ANCOM_long2_robust %>% 
  pivot_wider(names_from= sampleType, values_from = meanASVPctSampType)
head(I6S_airLeaf_ANCOM_wide_robust)
# View(I6S_airLeaf_ANCOM_wide_robust)

# When is foliar greater than air and vice versa?
I6S_airLeaf_ANCOM_wide2_robust <- I6S_airLeaf_ANCOM_wide_robust %>%
  mutate(whichGreater = case_when(
    air  > phyllo ~ "air",
    phyllo > air  ~ "foliar",
    air == phyllo ~ "ties",
    TRUE          ~ NA_character_     # any NA comparisons land here
  ))
# View(I6S_airLeaf_ANCOM_wide2_robust)
I6SbioFoliarBiggerIndex <- intersect(which(I6S_airLeaf_ANCOM_wide2_robust$ANCOMcat == "bioaerosol"), which(I6S_airLeaf_ANCOM_wide2_robust$whichGreater == "foliar"))
length(I6SbioFoliarBiggerIndex)/251*100 # 2% of the time, foliar more abundant in a biaoerosol indicator
I6SfoliarAirBiggerIndex <- intersect(which(I6S_airLeaf_ANCOM_wide2_robust$ANCOMcat == "foliar surface"), which(I6S_airLeaf_ANCOM_wide2_robust$whichGreater == "air"))
length(I6SfoliarAirBiggerIndex)/170*100 #3% of the time, bioaerosol more abundant in foliar indicator

# Add this to plot
I6S_airLeaf_ANCOM_wide2_robust$discrepencyANCOM <- NA
I6S_airLeaf_ANCOM_wide2_robust$discrepencyANCOM[I6SbioFoliarBiggerIndex] <- "foliar more abundant in bioaerosol indicator!"
I6S_airLeaf_ANCOM_wide2_robust$discrepencyANCOM[I6SfoliarAirBiggerIndex] <- "bioaerosol more abundant in foliar indicator!"

### MAIN RESULTS CHECKING PLOTS ####
head(I6S_airLeaf_ANCOM_long_robust)
# View(I6S_airLeaf_ANCOM_long_robust)
# GET THE PROPORTION OF READS IN EACH SAMPLE THAT COME FROM BIOAEROSOL AND FOLIAR ANCOM TAXA
I6S_airLeaf_ANCOM_longBySamp_robust <- I6S_airLeaf_ANCOM_long_robust %>% 
  group_by(sampleName, ANCOMcat) %>% 
  summarise(propSampANCOMcat = sum(propASVsamp)) #propASVsamp is proportion
# View(I6S_airLeaf_ANCOM_longBySamp_robust)
colnames(I6S_airLeaf_ANCOM_longBySamp_robust)
# Add back in sample type
# Add in air or leaf, depending on what sample type it is
I6S_airLeaf_ANCOM_longBySamp_robust <- I6S_airLeaf_ANCOM_longBySamp_robust %>%
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "bioaerosol",
      str_detect(sampleName, "phyllo") ~ "foliar surface",
      TRUE ~ NA_character_
    )
  )
I6S_airLeaf_ANCOM_longBySamp_robust
# Double check
# air_16S_105 and phyllo_16S_43 proportion of reads from two types of ANCOM taxa
ASVs_air105andPhyllo43_I6S <- I6SairPhylloOnly.ps_ASVs[which(rownames(I6SairPhylloOnly.ps_ASVs) %in% c("air_16S_105", "phyllo_16S_43")),]
totalReadsair105andPhyllo43_I6S <- rowSums(ASVs_air105andPhyllo43_I6S)
# Bioaerosol ASVs
bioaerosolANCOMASVs_air105andPhyllo43_I6S <- ASVs_air105andPhyllo43_I6S[,colnames(ASVs_air105andPhyllo43_I6S) %in% bioaerosolI6SANCOMASVsName_robust]
# Foliar surface ASVs
foliarANCOMASVs_air105andPhyllo43_I6S <- ASVs_air105andPhyllo43_I6S[,colnames(ASVs_air105andPhyllo43_I6S) %in% foliarI6SANCOMASVsName_robust]
propBioaerosol_air105andPhyllo43_I6S <- rowSums(bioaerosolANCOMASVs_air105andPhyllo43_I6S)/totalReadsair105andPhyllo43_I6S
propBioaerosol_air105andPhyllo43_I6S
propFoliar_air105andPhyllo43_I6S <- rowSums(foliarANCOMASVs_air105andPhyllo43_I6S)/totalReadsair105andPhyllo43_I6S
propFoliar_air105andPhyllo43_I6S
# Get indices
I6Sair105bioIndex <- intersect(which(I6S_airLeaf_ANCOM_longBySamp_robust$sampleName == "air_16S_105"),which(I6S_airLeaf_ANCOM_longBySamp_robust$ANCOMcat == "bioaerosol"))
I6Sair105FoliarIndex <- intersect(which(I6S_airLeaf_ANCOM_longBySamp_robust$sampleName == "air_16S_105"),which(I6S_airLeaf_ANCOM_longBySamp_robust$ANCOMcat == "foliar surface"))
I6Sfoliar43bioIndex <- intersect(which(I6S_airLeaf_ANCOM_longBySamp_robust$sampleName == "phyllo_16S_43"),which(I6S_airLeaf_ANCOM_longBySamp_robust$ANCOMcat == "bioaerosol"))
I6Sfoliar43FoliarIndex <- intersect(which(I6S_airLeaf_ANCOM_longBySamp_robust$sampleName == "phyllo_16S_43"),which(I6S_airLeaf_ANCOM_longBySamp_robust$ANCOMcat == "foliar surface"))

# These are correct, so calculation was correct
all.equal(I6S_airLeaf_ANCOM_longBySamp_robust$propSampANCOMcat[I6Sair105bioIndex],unname(propBioaerosol_air105andPhyllo43_I6S[2]))
all.equal(I6S_airLeaf_ANCOM_longBySamp_robust$propSampANCOMcat[I6Sair105FoliarIndex],unname(propFoliar_air105andPhyllo43_I6S[2]))
all.equal(I6S_airLeaf_ANCOM_longBySamp_robust$propSampANCOMcat[I6Sfoliar43bioIndex],unname(propBioaerosol_air105andPhyllo43_I6S[1]))
all.equal(I6S_airLeaf_ANCOM_longBySamp_robust$propSampANCOMcat[I6Sfoliar43FoliarIndex],unname(propFoliar_air105andPhyllo43_I6S[1]))

# 1. For each sample, what proportion of reads are bioaerosol versus foliar surface taxa?
colnames(I6S_airLeaf_ANCOM_longBySamp_robust)

propANCOMByReadsPlot_I6S_robust <- ggplot(I6S_airLeaf_ANCOM_longBySamp_robust, aes(x = sampleType, y = propSampANCOMcat, fill=sampleType)) +
  geom_boxplot(outlier.shape = NA) +           
  scale_fill_manual(values=c("cornflowerblue", "forestgreen")) +
  geom_jitter(width = 0.15, size = 1.8, height =0) + #each sample point
  facet_wrap(~ ANCOMcat, nrow = 1)   +          # ANCOM categoey by sample
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(name= "Proportion bioaerosol or \nfoliar surface ANCOM taxa")
propANCOMByReadsPlot_I6S_robust
# saveRDS(propANCOMByReadsPlot_I6S_robust, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/propANCOMByReadsPlot_I6S_robust_I6S_10-9-25.rds")

##### ii. JUST STRUCTURAL ZEROS #####
I6S_AFS_zeros_df
# 1. Create ASV table with only these struct zero ANCOM ASVs in air and leaf samples. Can use this to compare how ANCOM model and actual data compare
I6S_leafAir_zero_ASVs <- unique(I6S_AFS_zeros_df$ASV_name)
length(I6S_leafAir_zero_ASVs) #9440 different ASVs
# Pull out only ASV table for significantly different ASVs (ASVs are columns!)
I6S_AFS_zerosASVtab <- I6SairPhylloOnly.ps_ASVs[,colnames(I6SairPhylloOnly.ps_ASVs) %in% I6S_leafAir_zero_ASVs]
setdiff(colnames(I6S_AFS_zerosASVtab), I6S_leafAir_zero_ASVs ) #these are the same, as expected
# View(I6S_AFS_zerosASVtab)
head(I6S_AFS_zerosASVtab)
# But which samples are air versus leaves?
#View(I6S_AFS_zerosASVtab)
I6S_AFS_zerosASVtab <- as.data.frame(I6S_AFS_zerosASVtab) %>% 
  rownames_to_column(var= "sampleName") %>% 
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "air",
      str_detect(sampleName, "phyllo") ~ "phyllo",
      TRUE ~ NA_character_
    )
  )
cbind(I6S_AFS_zerosASVtab$sampleName, I6S_AFS_zerosASVtab$sampleType) #quick check to make sure code above worked... and it did!
# One last check before calculating...
all.equal(sort(I6S_leafAir_zero_ASVs), sort(colnames(I6S_AFS_zerosASVtab)[2:(ncol(I6S_AFS_zerosASVtab)-1)])) #okay these are still the significant ones

# 2. Add in number of reads per sample (readNumberI6SairFoliar_df made higher above)
#Add to working dataframe
head(I6S_AFS_zerosASVtab)
I6S_AFS_zeroASVs <- merge(I6S_AFS_zerosASVtab, readNumberI6SairFoliar_df, by = "sampleName")
# View(I6S_AFS_zeroASVs)

# 3. MAKE A "LONGER" DATAFRAME AND ADD IN TAXONOMY
head(I6S_AFS_zeroASVs) 
# i. First, pivot longer so that ASVs are column 1 and Sample is second column 
head(colnames(I6S_AFS_zeroASVs)); tail(colnames(I6S_AFS_zeroASVs)) #get first and last one to put in below
I6S_AFS_zeroASVs_long <- I6S_AFS_zeroASVs %>% 
  pivot_longer(cols = ASV_10:ASV_39616, names_to= "ASV_name", values_to = "ASV_abundance")
head(I6S_AFS_zeroASVs_long)
# View(I6S_AFS_zeroASVs_long)
dim(I6S_AFS_zeroASVs_long) #1,340,480       5
# ii. Below adds in taxonomy with result of abundance of each ANCOM ASV in each air and phyllo sample, plus taxonomy!
I6S_AFSzeros_long <- left_join(I6S_AFS_zeroASVs_long, I6SairPhylloOnly_tax, by = "ASV_name")
head(I6S_AFSzeros_long)
# View(I6S_AFSzeros_long)
# iii. Add in which sample the structural zero came from
colnames(I6S_AFS_zeros_df)
I6S_AFSzeros_long <- merge(I6S_AFSzeros_long, I6S_AFS_zeros_df, by ="ASV_name")
colnames(I6S_AFSzeros_long)

# 4. Get mean prop of each ASV in bioaerosol and foliar
# i. Get each ASVs prop in each sample
I6S_AFSzeros_long <- I6S_AFSzeros_long %>% 
  mutate(propASVsamp = ASV_abundance/TotalReads) #get each ASVs proportion in each sample. Double checked that this worked!
head(I6S_AFSzeros_long)
# ii. Get each ASVs mean percent of both sample types, keeping struct. zero info too
I6S_AFSzeros_long2 <- I6S_AFSzeros_long %>% 
  group_by(ASV_name, struct_ZeroType, sampleType) %>% 
  summarize(meanASVPctSampType = mean(propASVsamp)*100)
# View(I6S_AFSzeros_long2) #these match those I checked in Table S5!

# Make this wider for direct comparison
colnames(I6S_AFSzeros_long2)
I6S_AFS_zeros_wide <- I6S_AFSzeros_long2 %>% 
  pivot_wider(names_from= sampleType, values_from = meanASVPctSampType)
head(I6S_AFS_zeros_wide)
# View(I6S_AFS_zeros_wide)

# When is foliar greater than air and vice versa?
I6S_AFS_zeros_wide2 <- I6S_AFS_zeros_wide %>%
  mutate(whichGreater = case_when(
    air  > phyllo ~ "air",
    phyllo > air  ~ "foliar",
    air == phyllo ~ "ties",
    TRUE          ~ NA_character_     # any NA comparisons land here
  ))
# View(I6S_AFS_zeros_wide2)
I6SbioFoliarBiggerIndex <- intersect(which(I6S_AFS_zeros_wide2$struct_ZeroType == "bioaerosol"), which(I6S_AFS_zeros_wide2$whichGreater == "foliar"))
length(I6SbioFoliarBiggerIndex)/length(airZeros16S ) # As expected, this is always true because foliar always air str. zeros
I6SfoliarAirBiggerIndex <- intersect(which(I6S_AFS_zeros_wide2$struct_ZeroType == "foliar surfaces"), which(I6S_AFS_zeros_wide2$whichGreater == "air"))
length(I6SfoliarAirBiggerIndex)/length(phylloZeros16S) #100%, as expected and confirmed above

### STRUCTURAL ZEROS CHECKING PLOTS ####
head(I6S_AFSzeros_long)
# View(I6S_AFSzeros_long)
# GET THE PROPORTION OF READS IN EACH SAMPLE THAT COME FROM BIOAEROSOL AND FOLIAR ANCOM TAXA
I6S_AFSzeros_longBySamp <- I6S_AFSzeros_long %>% 
  group_by(sampleName, struct_ZeroType) %>% 
  summarise(propSampStrZeroType = sum(propASVsamp)) #propASVsamp is proportion
# View(I6S_AFSzeros_longBySamp)
colnames(I6S_AFSzeros_longBySamp)
# Add back in sample type
# Add in air or leaf, depending on what sample type it is
I6S_AFSzeros_longBySamp <- I6S_AFSzeros_longBySamp %>%
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "bioaerosol",
      str_detect(sampleName, "phyllo") ~ "foliar surface",
      TRUE ~ NA_character_
    )
  )
I6S_AFSzeros_longBySamp

# 1. For each sample, what proportion of reads are bioaerosol versus foliar surface taxa?
# Make copy to edit
I6S_AFSzeros_longBySampPlot_df <- I6S_AFSzeros_longBySamp
# make this a factor to re-order facets
I6S_AFSzeros_longBySampPlot_df$struct_ZeroType <- factor(I6S_AFSzeros_longBySampPlot_df$struct_ZeroType, levels = c("foliar surfaces", "bioaerosol"))
propZeroASVsByReadsPlot_I6S <- ggplot(I6S_AFSzeros_longBySampPlot_df, aes(x = sampleType, y = propSampStrZeroType, fill=sampleType)) +
  geom_boxplot(outlier.shape = NA) +           
  scale_fill_manual(values=c("cornflowerblue", "forestgreen")) +
  geom_jitter(width = 0.15, size = 1.8, height =0) + #each sample point
  facet_wrap(~ struct_ZeroType, nrow = 1, labeller = labeller(struct_ZeroType =  c("bioaerosol" = "taxa not found in bioaerosols", "foliar surfaces" = "taxa not found on foliar surfaces")))   +          #structural zero type by category
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(name= "Proportion bioaerosol or \nfoliar surface from str. zero Taxa")
propZeroASVsByReadsPlot_I6S
# saveRDS(propZeroASVsByReadsPlot_I6S, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/propZeroASVsByReadsPlot_I6S_10-6-25.rds")

########## 5. CHANGE THRESHOLDS FOR INCLUDING TAXA, AFTER ANCOMBC ##########
###### 1. PULL OUT ONLY ASV TABLE WITH ALL ANCOM ASVS ######
head(I6S_AFS_zeros_df) #has structural zeros for bioaerosols and foliar surfaces
bioaerosolI6SANCOMASVsName_robust #bioaerosol 16S ANCOM$res ASV names
foliarI6SANCOMASVsName_robust #foliarI6S ANCOM$res names
# i. Concatenate all of these ASVs together in order to get all of the ANCOM-identifed ASVs (struc zero and main)
allI6S_AFS_ANCOM_ASVs_robust <- c(I6S_AFS_zeros_df$ASV_name, bioaerosolI6SANCOMASVsName_robust, foliarI6SANCOMASVsName_robust)
length(allI6S_AFS_ANCOM_ASVs_robust) == length(unique(allI6S_AFS_ANCOM_ASVs_robust)) #no repeats, as expected
# ii. Get an ASV table with only these ANCOM ASVs
I6S_BAFS_ancom_sigASVtab_robust <- I6SairPhylloOnly.ps_ASVs[,colnames(I6SairPhylloOnly.ps_ASVs) %in% allI6S_AFS_ANCOM_ASVs_robust]
# iii. Transpose so that ASVs are now rows to add in additional ASV info
I6S_BAFS_ancom_sigASVtab_t_robust <- as.data.frame(t(I6S_BAFS_ancom_sigASVtab_robust))
if (nrow(I6S_BAFS_ancom_sigASVtab_t_robust) < ncol(I6S_BAFS_ancom_sigASVtab_t_robust)) {
  print("STOP! ASVs need to be rows!")
}
# iv. Add in sum in bioaerosol and phyllosphere samples
# Sum of each ASV in foliar surface
I6S_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceTotal <- rowSums(I6S_BAFS_ancom_sigASVtab_t_robust[,which(colnames(I6S_BAFS_ancom_sigASVtab_t_robust)%in% I6SphylloSampNames)])
# Occupancy in foliar surface
I6S_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceOcc <- rowSums(I6S_BAFS_ancom_sigASVtab_t_robust[, which(colnames(I6S_BAFS_ancom_sigASVtab_t_robust)%in% I6SphylloSampNames)] > 0, na.rm = TRUE)
# Total read counts in bioaerosols
I6S_BAFS_ancom_sigASVtab_t_robust$bioaerosolTotal <- rowSums(I6S_BAFS_ancom_sigASVtab_t_robust[,which(colnames(I6S_BAFS_ancom_sigASVtab_t_robust)%in% I6SairSampNames)])
# Occupancy in bioaerosol
I6S_BAFS_ancom_sigASVtab_t_robust$bioaerosolOcc <- rowSums(I6S_BAFS_ancom_sigASVtab_t_robust[, which(colnames(I6S_BAFS_ancom_sigASVtab_t_robust)%in% I6SairSampNames)] > 0, na.rm = TRUE)

# v. Add in which kind of ANCOM taxon it is
I6S_BAFS_ancom_sigASVtab_t_robust <- I6S_BAFS_ancom_sigASVtab_t_robust %>% 
  tibble::rownames_to_column(var= "ASV_name") %>% #make rownames a column
  mutate(ANCOMcat = case_when(
    ASV_name %in% c(phylloZeros16S, bioaerosolI6SANCOMASVsName_robust) ~ "bioaerosol", #bioaerosol if not in phyllo or $res bioaerosol
    ASV_name %in% c(airZeros16S, foliarI6SANCOMASVsName_robust) ~ "foliar surface", #foliar if not in air or $res bioaerosol
  )) %>% 
  mutate(resType = case_when(
    ASV_name %in% c(bioaerosolI6SANCOMASVsName_robust, foliarI6SANCOMASVsName_robust) ~ "mainNotStructZero",
    ASV_name %in% c(phylloZeros16S, airZeros16S) ~ "structZero"
  )
  )
# vi. Perform a few checks:
# 1. Get indices
I6S_BAFtotal_ASV10_index_robust <- which(I6S_BAFS_ancom_sigASVtab_t_robust$ASV_name == "ASV_10")
I6S_BAFtotal_ASV231_index_robust <- which(I6S_BAFS_ancom_sigASVtab_t_robust$ASV_name == "ASV_231")
# 2. Is total number working?
I6S_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceTotal[I6S_BAFtotal_ASV10_index_robust] == sum(I6S_BAFS_ancom_sigASVtab_t_robust[I6S_BAFtotal_ASV10_index_robust, colnames(I6S_BAFS_ancom_sigASVtab_t_robust)%in% I6SphylloSampNames])
I6S_BAFS_ancom_sigASVtab_t_robust$bioaerosolTotal[I6S_BAFtotal_ASV231_index_robust] == sum(I6S_BAFS_ancom_sigASVtab_t_robust[I6S_BAFtotal_ASV231_index_robust, colnames(I6S_BAFS_ancom_sigASVtab_t_robust)%in% I6SairSampNames])
# 3. Is occupancy working?
I6S_BAFS_ancom_sigASVtab_t_robust$bioaerosolOcc[I6S_BAFtotal_ASV10_index_robust] == length(which(I6S_BAFS_ancom_sigASVtab_t_robust[I6S_BAFtotal_ASV10_index_robust, colnames(I6S_BAFS_ancom_sigASVtab_t_robust)%in% I6SairSampNames]>0))
# 4. Are they labeled correctly?
setdiff(I6S_BAFS_ancom_sigASVtab_t_robust$ASV_name[which(I6S_BAFS_ancom_sigASVtab_t_robust$resType == "structZero")], c(phylloZeros16S, airZeros16S)) #no differences!
# 5. Are ANCOM categories correct?
setdiff(I6S_BAFS_ancom_sigASVtab_t_robust$ASV_name[which(I6S_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "foliar surface")], c(airZeros16S, foliarI6SANCOMASVsName_robust)) #no differences!
# View(I6S_BAFS_ancom_sigASVtab_t_robust)

# vii. Add in taxonomy!
I6S_BAFS_ancom_sigASVtab_t_robust <- left_join(I6S_BAFS_ancom_sigASVtab_t_robust, I6SairPhylloOnly_tax[,1:7], by = "ASV_name") #don't have last column because "species" is just ASV name for this analysis

# viii. Add in the proportion of reads in each category that an ASV contributes
# Use this, made above, to get total number of reads in folair surface and air samples
totalReadsType <- readNumberI6SairFoliar_df %>% 
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "air",
      str_detect(sampleName, "phyllo") ~ "phyllo",
      TRUE ~ NA_character_
    )
  ) %>% 
  group_by(sampleType) %>% 
  summarize(nReadsType = sum(TotalReads))
totalReadsType
# sampleType nReadsType
# <chr>           <dbl>
#   1 air           1577986
# 2 phyllo        1562809

# If I did 0.05%, or 5 out of every 10,000 reads it's roughly the same as 800 times at least
0.05/100*1577986 #788.993
0.05/100*1562809 #781.4045

colnames(I6S_BAFS_ancom_sigASVtab_t_robust)
I6S_BAFS_ancom_sigASVtab_t_robust <- I6S_BAFS_ancom_sigASVtab_t_robust %>% 
  mutate(pctBioaerosol = bioaerosolTotal/1577986*100) %>% 
  mutate(pctFoliar = foliarSurfaceTotal/1562809*100)

# DECISION AND FILTERING TIME:
# 1. To be considered diff. abundant, an ASV must be in at about 10% of respective group's samples (8 bioaerosol or 6 foliar surface) AND
# 2. Must be present at a certain threshold (see ways I've tried this below)
##### TEST BIAEROSOL ####
round(length(I6SairSampNames)/10) #8 samples 
airASVs_8occ50Abund_index_robust <- intersect(intersect(which(I6S_BAFS_ancom_sigASVtab_t_robust$bioaerosolOcc >=8), which(I6S_BAFS_ancom_sigASVtab_t_robust$bioaerosolTotal >=50)), which(I6S_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "bioaerosol"))
length(airASVs_8occ50Abund_index_robust ) #211
airASVs_8occ75Abund_index_robust <- intersect(intersect(which(I6S_BAFS_ancom_sigASVtab_t_robust$bioaerosolOcc >=8), which(I6S_BAFS_ancom_sigASVtab_t_robust$bioaerosolTotal >=75)), which(I6S_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "bioaerosol"))
length(airASVs_8occ75Abund_index_robust) #202
# Found in at least 8 samples and .01% of reads
airASVs_8occ.01pct_index_robust <- intersect(intersect(which(I6S_BAFS_ancom_sigASVtab_t_robust$bioaerosolOcc >=8), which(I6S_BAFS_ancom_sigASVtab_t_robust$pctBioaerosol >=0.01)), which(I6S_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "bioaerosol"))
length(airASVs_8occ.01pct_index_robust) #184
# View(I6S_BAFS_ancom_sigASVtab_t_robust[airASVs_8occ.01pct_index_robust,])
bioaerosolI6S_0.01pct_ASVsNames <- I6S_BAFS_ancom_sigASVtab_t_robust$ASV_name[airASVs_8occ.01pct_index_robust]
# Found in at least 8 samples and .05% of reads
airASVs_8occ.05pct_index_robust <- intersect(intersect(which(I6S_BAFS_ancom_sigASVtab_t_robust$bioaerosolOcc >=8), which(I6S_BAFS_ancom_sigASVtab_t_robust$pctBioaerosol >=0.05)), which(I6S_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "bioaerosol"))
length(airASVs_8occ.05pct_index_robust) #94
# View(I6S_BAFS_ancom_sigASVtab_t_robust[airASVs_8occ.05pct_index_robust,])
bioaerosolI6S_0.05pct_ASVsNames <- I6S_BAFS_ancom_sigASVtab_t_robust$ASV_name[airASVs_8occ.05pct_index_robust]
length(bioaerosolI6S_0.05pct_ASVsNames) #94

##### TEST FOLIAR SURFACE ####
round(length(I6SphylloSampNames)/10) #6 samples
foliarASVs_6occ50Abund_index_robust <- intersect(intersect(which(I6S_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceOcc >=6), which(I6S_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceTotal >=50)), which(I6S_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "foliar surface"))
length(foliarASVs_6occ50Abund_index_robust) #1855
foliarASVs_6occ100Abund_index_robust <- intersect(intersect(which(I6S_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceOcc >=6), which(I6S_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceTotal >=100)), which(I6S_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "foliar surface"))
length(foliarASVs_6occ100Abund_index_robust) #952
# Found in at least 6 samples and at least 0.01% of reads
foliarASVs_6occ.01pct_index_robust <- intersect(intersect(which(I6S_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceOcc >=6), which(I6S_BAFS_ancom_sigASVtab_t_robust$pctFoliar >=0.01)), which(I6S_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "foliar surface"))
length(foliarASVs_6occ.01pct_index_robust) #580
# View(I6S_BAFS_ancom_sigASVtab_t_robust[foliarASVs_6occ.01pct_index_robust,])
foliarI6S_0.01pct_ASVsNames <- I6S_BAFS_ancom_sigASVtab_t_robust$ASV_name[foliarASVs_6occ.01pct_index_robust]
# Found in at least 6 samples and at least 0.05% of reads
foliarASVs_6occ.05pct_index_robust <- intersect(intersect(which(I6S_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceOcc >=6), which(I6S_BAFS_ancom_sigASVtab_t_robust$pctFoliar >=0.05)), which(I6S_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "foliar surface"))
length(foliarASVs_6occ.05pct_index_robust) #126
# View(I6S_BAFS_ancom_sigASVtab_t_robust[foliarASVs_6occ.05pct_index_robust,])
foliarI6S_0.05pct_ASVsNames <- I6S_BAFS_ancom_sigASVtab_t_robust$ASV_name[foliarASVs_6occ.05pct_index_robust]
length(foliarI6S_0.05pct_ASVsNames) #126

###### PLOT AND EXPLORE NEW RESULTS! #####
# GO with 0.05%
# 1. Make a new dataframe that just has these new, differentially abundant taxa by concatenating indices from above
I6S_ANCOMall_.05pct_df <- I6S_BAFS_ancom_sigASVtab_t_robust[c(airASVs_8occ.05pct_index_robust, foliarASVs_6occ.05pct_index_robust),]
setdiff(I6S_ANCOMall_.05pct_df$ASV_name, c(bioaerosolI6S_0.05pct_ASVsNames, foliarI6S_0.05pct_ASVsNames)) #these are the same!
# View(I6S_ANCOMall_.05pct_df)
colnames(I6S_ANCOMall_.05pct_df)
# A few checks!
# all look good! (Minimums about than 0.05)
sort(I6S_ANCOMall_.05pct_df$pctBioaerosol[which(I6S_ANCOMall_.05pct_df$ANCOMcat == "bioaerosol")]) #minimum looks good!
sort(I6S_ANCOMall_.05pct_df$pctFoliar[which(I6S_ANCOMall_.05pct_df$ANCOMcat == "foliar surface")]) #minimum looks good!

length(I6S_ANCOMall_.05pct_df$ASV_name) #220
# Unique families and orders
table(I6S_ANCOMall_.05pct_df$ANCOMcat, I6S_ANCOMall_.05pct_df$Order)

# Saved October 11, 2025:
# saveRDS(I6S_ANCOMall_.05pct_df, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_ANCOMall_.05pct_df.rds")

# 2. Get proportion of each sample that is foliar surface or bioaerosols ANCOM taxa
# View(I6S_ANCOMall_.05pct_df)
colnames(I6S_ANCOMall_.05pct_df)
# i. Get each ASV's abundance in each sample
I6S_ANCOMall_.05pct_df_long <- I6S_ANCOMall_.05pct_df %>% 
  pivot_longer(cols = phyllo_16S_1:air_16S_97, names_to = "sampleName", values_to = "ASV_abundance")
# ii. Grouping by sampName and ANCOMtype, get sum of ANCOM type in each sample
head(I6S_ANCOMall_.05pct_df_long)
I6S_nANCOMcatSamp <- I6S_ANCOMall_.05pct_df_long %>% 
  select(c(ASV_name, ANCOMcat, sampleName, ASV_abundance)) %>% 
  group_by(ANCOMcat, sampleName) %>% 
  summarize(nANCOMsamp = sum(ASV_abundance))
head(I6S_nANCOMcatSamp)
# Check for a random sample to ensure it's working. Get index where samp is 89 and ANCOMcat is bioaerosol. Then pull out ASV_abundance and sum it. True so this is working!
sum(I6S_ANCOMall_.05pct_df_long$ASV_abundance[intersect(which(I6S_ANCOMall_.05pct_df_long$sampleName == "air_16S_89"), which(I6S_ANCOMall_.05pct_df_long$ANCOMcat == "bioaerosol"))]) == 
  I6S_nANCOMcatSamp$nANCOMsamp[intersect(which(I6S_nANCOMcatSamp$ANCOMcat == "bioaerosol"), which(I6S_nANCOMcatSamp$sampleName == "air_16S_89"))]
# iii. Add in with read number data and then divide by this for proportion
head(readNumberI6SairFoliar_df) #made above, is total number of reads per sample
I6S_nANCOMcatSamp2 <- left_join(I6S_nANCOMcatSamp, readNumberI6SairFoliar_df, by = "sampleName")
head(I6S_nANCOMcatSamp2)
# Add in prop and sample Type
I6S_nANCOMcatSamp2 <- I6S_nANCOMcatSamp2 %>% 
  mutate(propANCOMcatSamp = nANCOMsamp/TotalReads) %>% 
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "air",
      str_detect(sampleName, "phyllo") ~ "phyllo",
      TRUE ~ NA_character_
    )
  )

# 3. Plot the proportion of eaach sample that is foliar surface or bioaerosols!
head(I6S_nANCOMcatSamp2)
I6S_propANCOM_filt_Reads_plot <- ggplot(I6S_nANCOMcatSamp2, aes(x = sampleType, y = propANCOMcatSamp, fill=sampleType)) +
  geom_boxplot(outlier.shape = NA) +           
  scale_fill_manual(values=c("cornflowerblue", "forestgreen")) +
  geom_jitter(width = 0.15, size = 1.8, height =0) + #each sample point
  facet_wrap(~ ANCOMcat, nrow = 1)   +          # ANCOM category by sample
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(name= "Proportion bioaerosol or \nfoliar surface ANCOM taxa") +
  ggtitle("16S after filtering ASVs: proportion of each sample's reads by ANCOM type")
I6S_propANCOM_filt_Reads_plot

# saveRDS(I6S_propANCOM_filt_Reads_plot, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_propANCOM_filt_Reads_plot_10-11-25.rds")

# 2. Get proportion of various ANCOM ASVs in foliar surface or bioaerosol reads
# This dataframe, made above, has the ASV abundance of each ANCOM-ASV in each sample 
head(I6S_ANCOMall_.05pct_df_long)
# i. Add in number of reads per sample
head(readNumberI6SairFoliar_df) #made above, is total number of reads per sample
I6S_ANCOMall_.05pct_df_long_2 <- left_join(I6S_ANCOMall_.05pct_df_long, readNumberI6SairFoliar_df, by = "sampleName")
head(I6S_ANCOMall_.05pct_df_long_2)
# ii. Add a column to get proportion of each ASV in each sample
# Add in prop and sample Type
I6S_ANCOMall_.05pct_df_long_2 <- I6S_ANCOMall_.05pct_df_long_2 %>% 
  mutate(propASVPerSamp = ASV_abundance/TotalReads) %>% 
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "air",
      str_detect(sampleName, "phyllo") ~ "phyllo",
      TRUE ~ NA_character_
    )
  )
head(I6S_ANCOMall_.05pct_df_long_2)
### NOTE: CODE BELOW IS CORRECT BUT IS COMMENTED OUT BECAUSE ADDS SO MUCH TO THE 
# MEMORY. ALSO, I LOOKED OVER THESE PLOTS AND ALL SEEMS GOOD!

# iii. Make a list of dataframes, where each dataframe in the list has the data for a particular species
I6S_ANCOMASVs_list <- I6S_ANCOMall_.05pct_df_long_2 %>%
  group_by(ASV_name) %>%
  group_split(.keep = TRUE)
# View(I6S_ANCOMASVs_list[[1]])
# Name each list element by the ASV_name
names(I6S_ANCOMASVs_list) <- character(length(I6S_ANCOMASVs_list))
for (i in seq_along(I6S_ANCOMASVs_list)) {
  names(I6S_ANCOMASVs_list)[i] <- paste0(
    unique(I6S_ANCOMASVs_list[[i]]$ASV_name), "_",
    unique(I6S_ANCOMASVs_list[[i]]$ANCOMcat), "ANCOM"
  )
}

# # 1. Here is a plot just to get the syntax correct
# head(I6S_ANCOMASVs_list[[1]])
# I6S_propASVsamp_1 <- ggplot(I6S_ANCOMASVs_list[[1]], aes(x = sampleType, y = propASVPerSamp, fill=sampleType)) +
#   geom_boxplot(outlier.shape = NA) +           
#   scale_fill_manual(values=c("cornflowerblue", "forestgreen")) +
#   geom_jitter(width = 0.15, size = 1.8, height =0) + #each sample point
#   #facet_wrap(~ ANCOMcat, nrow = 1)   +          # ANCOM category by sample
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   scale_y_continuous(name= "Proportion in each bioaerosol or \nfoliar surface sample") +
#   ggtitle(paste0("16S_", names(I6S_ANCOMASVs_list)[1], "_", unique(I6S_ANCOMASVs_list[[1]]$Class)))
# I6S_propASVsamp_1
# 
# # 2. Make a separate list of plots for biaoerosol and foliar surface ANCOM taxa, to make
# # visualization easier
# First get separate lists of data for bioaerosol and foliar surface taxa
I6SbioaerosolANCOMlist <- I6S_ANCOMASVs_list[grepl(x= names(I6S_ANCOMASVs_list), pattern = "bioaerosol")]
I6SFoliarANCOMlist <- I6S_ANCOMASVs_list[grepl(x= names(I6S_ANCOMASVs_list), pattern = "foliar surface")]

# ### NOTE: PLOT BELOW IS CORRECT BUT IS COMMENTED OUT BECAUSE ADDS SO MUCH TO THE
# # MEMORY. ALSO, I LOOKED OVER THESE PLOTS AND ALL SEEMS GOOD!
# # 3. Bioaerosol plots
# I6S_ANCOM_BA_Plots <- vector(mode = "list", length = length(I6SbioaerosolANCOMlist))
# for (i in seq_along(I6SbioaerosolANCOMlist)) {
#   ANCOM16S_BA_df <- I6SbioaerosolANCOMlist[[i]]
#   ANCOM16S_BAttl <- paste0(
#     "16S_", names(I6SbioaerosolANCOMlist)[i], "_",
#     paste(unique(ANCOM16S_BA_df$Class), collapse = ", ")
#   )
# 
#   I6S_ANCOM_BA_Plots[[i]] <- try(
#     ggplot(ANCOM16S_BA_df, aes(x = sampleType, y = propASVPerSamp, fill = sampleType)) +
#       geom_boxplot(outlier.shape = NA) +
#       scale_fill_manual(values = c("cornflowerblue", "forestgreen")) +
#       geom_jitter(width = 0.15, size = 1.8, height = 0) +
#       theme_bw() +
#       theme(panel.grid = element_blank()) +
#       scale_y_continuous(name = "Proportion in each bioaerosol or \nfoliar surface sample") +
#       ggtitle(ANCOM16S_BAttl),
#     silent = TRUE
#   )
# }
# 
# # NOTE: PLOT BELOW IS CORRECT BUT IS COMMENTED OUT BECAUSE ADDS SO MUCH TO THE
# #MEMORY. ALSO, I LOOKED OVER THESE PLOTS AND ALL SEEMS GOOD!
# # 3. Foliar surface plots. Break these into 5 for memory/plotting reasons
# I6S_ANCOM_FS_Plots <- vector(mode = "list", length = length(I6SFoliarANCOMlist))
# for (i in seq_along(I6SFoliarANCOMlist)) {
#   ANCOM16S_FS_df <- I6SFoliarANCOMlist[[i]]
#   ANCOM16S_BAttl <- paste0(
#     "16S_", names(I6SFoliarANCOMlist)[i], "_",
#     paste(unique(ANCOM16S_FS_df$Class), collapse = ", ")
#   )
# 
#   I6S_ANCOM_FS_Plots[[i]] <- try(
#     ggplot(ANCOM16S_FS_df, aes(x = sampleType, y = propASVPerSamp, fill = sampleType)) +
#       geom_boxplot(outlier.shape = NA) +
#       scale_fill_manual(values = c("cornflowerblue", "forestgreen")) +
#       geom_jitter(width = 0.15, size = 1.8, height = 0) +
#       theme_bw() +
#       theme(panel.grid = element_blank()) +
#       scale_y_continuous(name = "Proportion in each bioaerosol or \nfoliar surface sample") +
#       ggtitle(ANCOM16S_BAttl),
#     silent = TRUE
#   )
# }


################################################################################################
# FUNGI
# subset no soil
unique(sample_data(ITSall_8.5Kfiltered.ps)$sampleType)
sort(colSums(as.data.frame(as.matrix(otu_table(ITSall_8.5Kfiltered.ps))))) #shows that all samples all vary (i.e. not rarefied)
min(colSums(as.data.frame(as.matrix(otu_table(ITSall_8.5Kfiltered.ps))))) #8846
max(colSums(as.data.frame(as.matrix(otu_table(ITSall_8.5Kfiltered.ps))))) #94248
ITSall_8.5Kfiltered.psAIR_LEAF <- subset_samples(ITSall_8.5Kfiltered.ps, sampleType != "soil")
unique(sample_data(ITSall_8.5Kfiltered.psAIR_LEAF)$sampleType)
# remove ASVs that do not have taxa present
ITSall_8.5Kfiltered.psAIR_LEAFZeros <- which(rowSums(otu_table(ITSall_8.5Kfiltered.psAIR_LEAF))==0)
unique(rowSums(otu_table(ITSall_8.5Kfiltered.psAIR_LEAF)[ITSall_8.5Kfiltered.psAIR_LEAFZeros,])) #all zeros
ITSall_8.5Kfiltered.psAIR_LEAFZerosNames <- names(ITSall_8.5Kfiltered.psAIR_LEAFZeros)
# Get what is unique in ITSall_8.5Kfiltered.psAIR_LEAF that is not in the zero ones!
ITSall_8.5Kfiltered.psAIR_LEAFASVsToKeep <- setdiff(rownames(otu_table(ITSall_8.5Kfiltered.psAIR_LEAF)), ITSall_8.5Kfiltered.psAIR_LEAFZerosNames) 
length(ITSall_8.5Kfiltered.psAIR_LEAFASVsToKeep) #9,012
ITSall_8.5Kfiltered.psAIR_LEAF.ps <- prune_taxa(taxa= ITSall_8.5Kfiltered.psAIR_LEAFASVsToKeep, x=ITSall_8.5Kfiltered.psAIR_LEAF) #9012 taxa and 169 samples 
which(rowSums(otu_table(ITSall_8.5Kfiltered.psAIR_LEAF.ps)) == 0) #none, as expected

# Add "Species" at ASV level and make a new phyloseq
ITSall_8.5Kfiltered.psAIR_LEAF_taxSpecies <- as.data.frame(as.matrix(tax_table(ITSall_8.5Kfiltered.psAIR_LEAF.ps)))
class(ITSall_8.5Kfiltered.psAIR_LEAF_taxSpecies) #great, it's not a phyloseq object!
head(ITSall_8.5Kfiltered.psAIR_LEAF_taxSpecies)
ITS_airLeaf_ASVnames <- rownames(ITSall_8.5Kfiltered.psAIR_LEAF_taxSpecies)
ITSall_8.5Kfiltered.psAIR_LEAF_taxSpecies$Species <- rownames(ITSall_8.5Kfiltered.psAIR_LEAF_taxSpecies)
head(ITSall_8.5Kfiltered.psAIR_LEAF_taxSpecies)
# Turn this new tax table into a tax table for phyloseq
ITSall_8.5Kfiltered.psAIR_LEAF_taxForPS <- phyloseq::tax_table(ITSall_8.5Kfiltered.psAIR_LEAF_taxSpecies)
colnames(ITSall_8.5Kfiltered.psAIR_LEAF_taxForPS) <- c("Kingdom","Phylum","Class", "Order","Family","Genus", "Species")
rownames(ITSall_8.5Kfiltered.psAIR_LEAF_taxForPS) <- ITS_airLeaf_ASVnames #re-set rownames
# View(ITSall_8.5Kfiltered.psAIR_LEAF_taxForPS)
head(ITSall_8.5Kfiltered.psAIR_LEAF_taxForPS) #great species are the kind they need to be!
class(ITSall_8.5Kfiltered.psAIR_LEAF_taxForPS)

# Finally, re-make phyloseq object, Have to force a bit below to keep the Species as ASV_X...
ITS_leafAir_notR_ANCOM.ps <- merge_phyloseq(otu_table(ITSall_8.5Kfiltered.psAIR_LEAF.ps), sample_data(ITSall_8.5Kfiltered.psAIR_LEAF.ps), ITSall_8.5Kfiltered.psAIR_LEAF_taxForPS)
head(tax_table(ITS_leafAir_notR_ANCOM.ps)) #9012    7
dim(otu_table(ITS_leafAir_notR_ANCOM.ps))
ITS_leafAir_ASVsTab <- as.data.frame(as.matrix(otu_table(ITS_leafAir_notR_ANCOM.ps)))
which(rowSums(ITS_leafAir_ASVsTab)<1) #none of these ASVs are NOT found in these groups! Great
dim(ITS_leafAir_ASVsTab) #9012  ASVs
unique(sample_data(ITS_leafAir_notR_ANCOM.ps)$sampleType) #"air"          "phyllosphere"

ITSairPhylloOnly.ps_ASVs <- t(as.data.frame(as.matrix(otu_table(ITS_leafAir_notR_ANCOM.ps))))
# Make sure that ASVs are columns
if (ncol(ITSairPhylloOnly.ps_ASVs)> nrow(ITSairPhylloOnly.ps_ASVs)){
  print("Good! ASVs are columns!")
} else {
  ITSairPhylloOnly.ps_ASVs <- t(ITSairPhylloOnly.ps_ASVs)
}

# Get names of air and phyllo samples to use later
ITSphylloSampNames <- rownames(sample_data(ITS_leafAir_notR_ANCOM.ps))[which(sample_data(ITS_leafAir_notR_ANCOM.ps)$sampleType == "phyllosphere")]
ITSairSampNames <- rownames(sample_data(ITS_leafAir_notR_ANCOM.ps))[which(sample_data(ITS_leafAir_notR_ANCOM.ps)$sampleType == "air")]

# taxonomy to use later with original species
ITSall_8.5Kfiltered_taxTab <- as.data.frame(as.matrix(tax_table(ITSall_8.5Kfiltered.ps))) 
head(ITSall_8.5Kfiltered_taxTab)
ITSall_8.5Kfiltered_taxTab <- ITSall_8.5Kfiltered_taxTab %>% 
  tibble::rownames_to_column(var="taxon")

# GET READ NUMBERS TO USE LATER
readNumberITSairFoliar <- rowSums(ITSairPhylloOnly.ps_ASVs) #get read numbers from whole dataset
readNumberITSairFoliar_df <- as.data.frame(matrix(nrow=length(readNumberITSairFoliar), ncol=2))
colnames(readNumberITSairFoliar_df) <- c("sampleName", "TotalReads")
readNumberITSairFoliar_df[,1] <- names(readNumberITSairFoliar)
readNumberITSairFoliar_df[,2] <- unname(readNumberITSairFoliar)

# September 30, 2025
#saveRDS(ITS_leafAir_notR_ANCOM.ps, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_leafAir_notR_ANCOM_phyloseq.rds")
# Read in data (made earlier in this script)
ITS_leafAir_notR_ANCOM.ps <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_leafAir_notR_ANCOM_phyloseq.rds")
head(tax_table(ITS_leafAir_notR_ANCOM.ps)) #double check that Species is ASV_X
# Following lines make sure that the factor is done correctly
ITSsampDat_Check <- as.data.frame(sample_data(ITS_leafAir_notR_ANCOM.ps))
sort(colnames(ITSsampDat_Check))
ITSsampDat_Check$sampleType <- stats::relevel(factor(ITSsampDat_Check$sampleType), ref = "air") 
levels(ITSsampDat_Check$sampleType)[1] 
sample_data(ITS_leafAir_notR_ANCOM.ps)$sampleType <- ITSsampDat_Check$sampleType
sample_data(ITS_leafAir_notR_ANCOM.ps)$sampleType #looks good! Air is first level and so will be intercept reference

# RUN ANCOM!
set.seed(19) #Grayed out since takes a long time to run!
# ITS_leafAir_notR_ANCOM_ASVlvl <- ancombc2(data = ITS_leafAir_notR_ANCOM.ps, group = "sampleType", tax_level = "Species", fix_formula = "sampleType", struc_zero = TRUE, neg_lb = FALSE) #fungi
# Saved/run October 1, 2025
# saveRDS(ITS_leafAir_notR_ANCOM_ASVlvl, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_leafAir_notR_ANCOM_ASVlvl.rds")

ITS_leafAir_notR_ANCOM_ASVlvl <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_leafAir_notR_ANCOM_ASVlvl.rds")
colnames(ITS_leafAir_notR_ANCOM_ASVlvl$res)

########## 3. SET UP/INVESTIGATE MAIN (NON STRUCTURAL ZEROS OUTPUT) ##########
# USE ROBUST
ITS_ANCOM_airLeaf_results_robust <- ITS_leafAir_notR_ANCOM_ASVlvl$res[which(ITS_leafAir_notR_ANCOM_ASVlvl$res$diff_robust_sampleTypephyllosphere == TRUE),]
# Good, shows that working as expected
length(which(ITS_ANCOM_airLeaf_results_robust$q_sampleTypephyllosphere[which(ITS_ANCOM_airLeaf_results_robust$diff_sampleTypephyllosphere == TRUE)]< 0.05)) #372 #still has high adj. p-value
length(which(ITS_ANCOM_airLeaf_results_robust$diff_sampleTypephyllosphere == TRUE)) #372 (incorporates less robust results)

# Add in taxonomy
# Get taxa (with original species)
ITS_ANCOM_airLeaf_results_robust <- left_join(ITS_ANCOM_airLeaf_results_robust, ITSall_8.5Kfiltered_taxTab, by= "taxon")
# View(ITS_ANCOM_airLeaf_results_robust)

# Look at phyllosphere enriched taxa
length(which(ITS_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere > 0)) #203
unique(ITS_ANCOM_airLeaf_results_robust$Order[which(ITS_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere > 0)]) #greater than zero are phyllosphere. 26 unique
sort(table(ITS_ANCOM_airLeaf_results_robust$Order[which(ITS_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere > 0)])) #no polypores. Capnodiales and Dothideales, Tremalles and Pelosporales are the largest, as before!

# Look at bioaerosol-enriched taxa
length(which(ITS_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere < 0)) #169
unique(ITS_ANCOM_airLeaf_results_robust$Order[which(ITS_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere < 0)]) #less than zero are bioaerosol. 14 unique orders
sort(table(ITS_ANCOM_airLeaf_results_robust$Order[which(ITS_ANCOM_airLeaf_results_robust$lfc_sampleTypephyllosphere < 0)])) #Polyporales are 91, then Hymenochaetales (26), Russulales (11)

# Saved October 9, 2025:
# saveRDS(ITS_ANCOM_airLeaf_results_robust, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_ANCOM_airLeaf_results_robust_10-9-25.rds")

# Make dataframe smaller
colsToKeepANCOM_robust <- c("taxon", "lfc_sampleTypephyllosphere", "W_sampleTypephyllosphere", "q_sampleTypephyllosphere", "diff_robust_sampleTypephyllosphere", "Phylum", "Class",
                     "Order","Family","Genus","Species")
ITS_leafAir_ANCOM_small_robust <- ITS_ANCOM_airLeaf_results_robust[,colnames(ITS_ANCOM_airLeaf_results_robust) %in% colsToKeepANCOM_robust]
colnames(ITS_leafAir_ANCOM_small_robust)[colnames(ITS_leafAir_ANCOM_small_robust) %in% "taxon"] <- "ASV_name" #change to ASV name for merging
head(ITS_leafAir_ANCOM_small_robust)
# Add in category for bioaerosol or foliar surface
ITS_leafAir_ANCOM_small_robust$ANCOMcat[which(ITS_leafAir_ANCOM_small_robust$lfc_sampleTypephyllosphere > 0)] <- "foliar surface" #if lfc greater than 0, then foliar surface
ITS_leafAir_ANCOM_small_robust$ANCOMcat[which(ITS_leafAir_ANCOM_small_robust$lfc_sampleTypephyllosphere < 0)] <- "bioaerosol" #if lfc less than 0, then bioaerosol
# Double check
sort(ITS_leafAir_ANCOM_small_robust$lfc_sampleTypephyllosphere[which(ITS_leafAir_ANCOM_small_robust$ANCOMcat == "bioaerosol")]) #great, all are negative!
sort(ITS_leafAir_ANCOM_small_robust$lfc_sampleTypephyllosphere[which(ITS_leafAir_ANCOM_small_robust$ANCOMcat == "foliar surface")]) #great, all are positive!

foliarITSANCOMASVsName_robust <- ITS_leafAir_ANCOM_small_robust$ASV_name[ITS_leafAir_ANCOM_small_robust$ANCOMcat == "foliar surface"]
length(foliarITSANCOMASVsName_robust) #203
bioaerosolITSANCOMASVsName_robust <- ITS_leafAir_ANCOM_small_robust$ASV_name[ITS_leafAir_ANCOM_small_robust$ANCOMcat == "bioaerosol"]
length(bioaerosolITSANCOMASVsName_robust) #169
# Saved October 11, 2025
# saveRDS(ITS_leafAir_ANCOM_small_robust, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_leafAir_ANCOM_small_robust.rds")

# HOW DO THESE COMPARE TO ORIGINAL ANCOM?? Captured 102/237 OG phyllo, but only 121/150 OG air
OG_ANCOM_ITS <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/importedFromServer/ITS_MEAN_airLeafANCOM.Rdata")
head(OG_ANCOM_ITS)
OG_ITSairASVs <- OG_ANCOM_ITS$ASV_name[which(OG_ANCOM_ITS$ANCOMcat == "air")]; length(OG_ITSairASVs) #150
OG_ITSphylloASVs <- OG_ANCOM_ITS$ASV_name[which(OG_ANCOM_ITS$ANCOMcat == "phyllo")]; length(OG_ITSphylloASVs) #237
# Foliar
length(which(foliarITSANCOMASVsName_robust %in% OG_ITSphylloASVs ==TRUE)) #102 are from last time!
length(intersect(OG_ITSphylloASVs,foliarITSANCOMASVsName_robust)) #102are from last time!
# Bioaerosols
length(which(bioaerosolITSANCOMASVsName_robust %in% OG_ITSairASVs ==TRUE)) #114 are from last time!
length(intersect(OG_ITSairASVs,bioaerosolITSANCOMASVsName_robust)) #114 are from last time!

########## 3. INVESTIGATE STRUCTURAL ZEROS ##########
# UNDERSTAND THIS DATA FRAME. Overview of the findings below. All $res identified taxa are also in this dataframe, but 
# they are not structural zeros in either group. Structural zeros for phyllo are those detected only in air but not phyllo, and vice versa. However, these
# can vary WIDELY with their abundance, so important to not take all of these as diff. abundant.
ITS_airLeafANCOMtabzero <- ITS_leafAir_notR_ANCOM_ASVlvl$zero_ind
# View(ITS_airLeafANCOMtabzero)
# 1. Is there any overlap between these taxa and those in the main results? YES, ALL MAIN RESULTS PLUS MANY MORE
length(intersect(ITS_ANCOM_airLeaf_results_robust$taxon,ITS_airLeafANCOMtabzero$taxon)) #372
unique(ITS_ANCOM_airLeaf_results_robust$taxon %in% ITS_airLeafANCOMtabzero$taxon) #all of the original ones are in the structural zeros
ITS_airLeafANCOMtabzero$taxon[ITS_airLeafANCOMtabzero$taxon %in% ITS_ANCOM_airLeaf_results_robust$taxon] #there are many extra that are not found in main results
# 2. What do the main results taxa look like in the zero results? ALL FALSE, IN OTHER WORDS, THEY ARE NOT STRUCTURAL ZEROS, MAKES SENSE
unique(ITS_airLeafANCOMtabzero$`structural_zero (sampleType = air)`[ITS_airLeafANCOMtabzero$taxon %in% ITS_ANCOM_airLeaf_results_robust$taxon]) #all FALSE
unique(ITS_airLeafANCOMtabzero$`structural_zero (sampleType = phyllosphere)`[ITS_airLeafANCOMtabzero$taxon %in% ITS_ANCOM_airLeaf_results_robust$taxon]) #all FALSE
# 3. TRUE structural_zero (sampleType = air)`, what are these like? All are F for phyllo, makes sense
ITS_airLeafTab0_air <- ITS_airLeafANCOMtabzero[which(ITS_airLeafANCOMtabzero$`structural_zero (sampleType = air)` == TRUE),]
unique(ITS_airLeafTab0_air$`structural_zero (sampleType = phyllosphere)`) #FALSE
# 4. TRUE structural_zero (sampleType = phyllosphere)`, what are these like? All are F for air, makes sense
ITS_airLeafTab0_foliar <- ITS_airLeafANCOMtabzero[which(ITS_airLeafANCOMtabzero$`structural_zero (sampleType = phyllosphere)` == TRUE),]
unique(ITS_airLeafTab0_foliar$`structural_zero (sampleType = air)`) #all FALSE, makes sense
# 5. Are all of these air structural zeros not found in the air at all?
ITSair0ASV_names <- unique(ITS_airLeafTab0_air$taxon)
# Look for these in air samples -- correct, they are never found
unique(colSums(ITSairPhylloOnly.ps_ASVs[rownames(ITSairPhylloOnly.ps_ASVs) %in% ITSairSampNames,colnames(ITSairPhylloOnly.ps_ASVs) %in% ITSair0ASV_names]))
# 6. Are all of these air structural zeros found in the phyllosphere, as expected?
# Look for these in phyllo samples -- correct, are found many times! But, in some cases, these are in fact really rare, which is an issue because I don't 
# really want to use these as "differentially abundant" because they make no sense to analyze in a trait-based way.
sort(unique(colSums(ITSairPhylloOnly.ps_ASVs[rownames(ITSairPhylloOnly.ps_ASVs) %in% ITSphylloSampNames,colnames(ITSairPhylloOnly.ps_ASVs) %in% ITSair0ASV_names])))
# 7. Are all of these phyllo structural zeros not found in the phyllo at all?
ITSphyllo0ASV_names <- unique(ITS_airLeafTab0_foliar$taxon)
# Look for these in phyllo samples -- correct, they are never found
unique(colSums(ITSairPhylloOnly.ps_ASVs[rownames(ITSairPhylloOnly.ps_ASVs) %in% ITSphylloSampNames,colnames(ITSairPhylloOnly.ps_ASVs) %in% ITSphyllo0ASV_names]))
# 8. Are all of these phyllo structural zeros found in the air, as expected?
# Look for these in phyllo samples -- correct, are found many times! But, in some cases, these are in fact really rare, which is an issue because I don't 
# really want to use these as "differentially abundant" because they make no sense to analyze in a trait-based way. But then some are super abundant! Yikes, this is the issue!
sort(unique(colSums(ITSairPhylloOnly.ps_ASVs[rownames(ITSairPhylloOnly.ps_ASVs) %in% ITSairSampNames,colnames(ITSairPhylloOnly.ps_ASVs) %in% ITSphyllo0ASV_names])))

# Given all of the above, get a zero data frame that just has the structural zeros
# get an index for all phyllo or air struct zeros
szIndexITS <- c(which(ITS_airLeafANCOMtabzero$`structural_zero (sampleType = air)`== TRUE), which(ITS_airLeafANCOMtabzero$`structural_zero (sampleType = phyllosphere)`== TRUE))
# Get this dataframe and then clean it up some
ITS_AFS_zeros <- ITS_airLeafANCOMtabzero[szIndexITS,]
# Get air and phyllosphere structural zeros
airZerosITS <- ITS_AFS_zeros$taxon[which(ITS_AFS_zeros$`structural_zero (sampleType = air)` == TRUE)]
length(airZerosITS) #3400 in phyllo but not air
phylloZerosITS <- ITS_AFS_zeros$taxon[which(ITS_AFS_zeros$`structural_zero (sampleType = phyllosphere)` == TRUE)]
length(phylloZerosITS) #4018 in air but not phyllo
nrow(ITS_AFS_zeros) == length(unique(ITS_AFS_zeros$taxon))
# Make a cleaner df for this
ITS_AFS_zeros_df <- as.data.frame(matrix(ncol=2, nrow=length(unique(ITS_AFS_zeros$taxon))))
colnames(ITS_AFS_zeros_df) <- c("ASV_name", "struct_ZeroType")
ITS_AFS_zeros_df$ASV_name <- ITS_AFS_zeros$taxon #add in names
ITS_AFS_zeros_df <- ITS_AFS_zeros_df %>%
  mutate( #if ASV name in air zeros, struct tye is bioaerosol, is in phyllo, phyllo
    struct_ZeroType = case_when(
      ASV_name %in% airZerosITS ~ "bioaerosol",
      ASV_name %in% phylloZerosITS ~ "foliar surfaces",
      TRUE ~ NA_character_
    )
  )
unique(sort(ITS_AFS_zeros_df$ASV_name[which(ITS_AFS_zeros_df$struct_ZeroType == "bioaerosol")]) == sort(airZerosITS))  #double check

# HOW DO THESE COMPARE TO ORIGINAL ANCOM ANALYSIS? 123/237 phyllo were in these structrual, 29/150 air!
# Foliar
length(which(airZerosITS %in% OG_ITSphylloASVs ==TRUE)) #123 are from last time!
length(intersect(OG_ITSphylloASVs,airZerosITS)) #123, just another way of checking this
# Bioaerosols
length(which(phylloZerosITS %in% OG_ITSairASVs ==TRUE)) #29 are from last time!
length(intersect(OG_ITSairASVs,phylloZerosITS)) #still 29, so another way of checking this

########## 4. INVESTIGATE HOW ANCOM MODEL AND ACTUAL DATA COMPARE ##########
##### i. JUST MAIN $RES RESULTS #####
# 1. Create ASV table with only ANCOM ASVs in air and leaf samples. Can use this to compare how ANCOM model and actual data compare
head(ITS_leafAir_ANCOM_small_robust)
which(ITS_leafAir_ANCOM_small_robust$diff_robust_sampleTypephyllosphere ==FALSE) #good, all of these are true
ITS_leafAir_ANCOM_ASVs_robust <- ITS_leafAir_ANCOM_small_robust$ASV_name
length(ITS_leafAir_ANCOM_ASVs_robust) #372 different ASVs
# Pull out only ASV table for significantly different ASVs (ASVs are columns!)
ITS_airLeaf_sigASVtab_robust <- ITSairPhylloOnly.ps_ASVs[,colnames(ITSairPhylloOnly.ps_ASVs) %in% ITS_leafAir_ANCOM_ASVs_robust]
# View(ITS_airLeaf_sigASVtab_robust)
head(ITS_airLeaf_sigASVtab_robust)
# But which samples are air versus leaves?
#View(ITS_airLeaf_sigASVtab_robust)
ITS_airLeaf_sigASVtab_robust <- as.data.frame(ITS_airLeaf_sigASVtab_robust) %>% 
  rownames_to_column(var= "sampleName")
dim(ITS_airLeaf_sigASVtab_robust) #rows are samples, columns are ASVs
# View(ITS_airLeaf_sigASVtab_robust)

# Add in air or leaf, depending on what sample type it is
ITS_airLeaf_sigASVtab_robust <- ITS_airLeaf_sigASVtab_robust %>%
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "air",
      str_detect(sampleName, "phyllo") ~ "phyllo",
      TRUE ~ NA_character_
    )
  )
head(ITS_airLeaf_sigASVtab_robust)
cbind(ITS_airLeaf_sigASVtab_robust$sampleName, ITS_airLeaf_sigASVtab_robust$sampleType) #quick check to make sure code above worked... and it did!
# One last check before calculating...
all.equal(ITS_leafAir_ANCOM_ASVs_robust, colnames(ITS_airLeaf_sigASVtab_robust)[2:(ncol(ITS_airLeaf_sigASVtab_robust)-1)]) #okay these are still the significant ones

# 2. Add in number of reads per sample
#Add to working dataframe
head(ITS_airLeaf_sigASVtab_robust)
ITS_airLeaf_sigASVs_robust <- merge(ITS_airLeaf_sigASVtab_robust, readNumberITSairFoliar_df, by = "sampleName")
# View(ITS_airLeaf_sigASVs_robust)

# 3. MAKE A "LONGER" DATAFRAME AND ADD IN TAXONOMY
head(ITS_airLeaf_sigASVs_robust) 
# First, pivot longer so that ASVs are column 1 and Sample is second column 
colnames(ITS_airLeaf_sigASVs_robust)
ITS_airLeaf_sigASVs_robust_long <- ITS_airLeaf_sigASVs_robust %>% 
  pivot_longer(cols = ASV_1:ASV_7792, names_to= "ASV_name", values_to = "ASV_abundance")
head(ITS_airLeaf_sigASVs_robust_long)
# View(ITS_airLeaf_sigASVs_robust_long)
dim(ITS_airLeaf_sigASVs_robust_long) #62868     5
unique(ITS_airLeaf_sigASVs_robust_long$sampleType) #great, only air and phyllo!
colnames(ITS_airLeaf_sigASVs_robust_long)
# Result is "long" object with the abundance of each ANCOM ASV in each air and phyllo sample, plus taxonomy!
# Re-name taxon as ASV_name
ITSall_8.5Kfiltered_taxTab2 <- ITSall_8.5Kfiltered_taxTab
colnames(ITSall_8.5Kfiltered_taxTab2)[colnames(ITSall_8.5Kfiltered_taxTab2) %in% "taxon"] <- "ASV_name"
head(ITSall_8.5Kfiltered_taxTab2)
ITS_airLeaf_ANCOM_long_robust <- left_join(ITS_airLeaf_sigASVs_robust_long, ITSall_8.5Kfiltered_taxTab2, by = "ASV_name")
head(ITS_airLeaf_ANCOM_long_robust)
# View(ITS_airLeaf_ANCOM_long_robust)
# Add in ANCOM cat of each 
colnames(ITS_leafAir_ANCOM_small_robust)
ITS_airLeaf_ANCOM_long_robust <- merge(ITS_airLeaf_ANCOM_long_robust, ITS_leafAir_ANCOM_small_robust[,colnames(ITS_leafAir_ANCOM_small_robust) %in% c("ASV_name", "ANCOMcat")], by ="ASV_name")

# 4. Get mean prop of each ASV in bioaerosol and foliar
ITS_airLeaf_ANCOM_long_robust <- ITS_airLeaf_ANCOM_long_robust %>% 
  mutate(propASVsamp = ASV_abundance/TotalReads) #get each ASVs proportion in each sample. Double checked that this worked!
# View(ITS_airLeaf_ANCOM_long_robust)
# USE THIS DATAFRAME FOR LATER PLOTS
ITS_airLeaf_ANCOM_long_robust2 <- ITS_airLeaf_ANCOM_long_robust %>% 
  group_by(ASV_name, ANCOMcat, sampleType) %>% 
  summarize(meanASVPctSampType = mean(propASVsamp)*100)
# View(ITS_airLeaf_ANCOM_long_robust2)

# Make this wider for direct comparison
colnames(ITS_airLeaf_ANCOM_long_robust2)
ITS_airLeaf_ANCOM_wide_robust <- ITS_airLeaf_ANCOM_long_robust2 %>% 
  pivot_wider(names_from= sampleType, values_from = meanASVPctSampType)
head(ITS_airLeaf_ANCOM_wide_robust)
# View(ITS_airLeaf_ANCOM_wide_robust)

# When is foliar greater than air and vice versa?
ITS_airLeaf_ANCOM_wide_robust2 <- ITS_airLeaf_ANCOM_wide_robust %>%
  mutate(whichGreater = case_when(
    air  > phyllo ~ "air",
    phyllo > air  ~ "foliar",
    air == phyllo ~ "ties",
    TRUE          ~ NA_character_     # any NA comparisons land here
  ))
# View(ITS_airLeaf_ANCOM_wide_robust2)
ITSbioFoliarBiggerIndex <- intersect(which(ITS_airLeaf_ANCOM_wide_robust2$ANCOMcat == "bioaerosol"), which(ITS_airLeaf_ANCOM_wide_robust2$whichGreater == "foliar"))
length(ITSbioFoliarBiggerIndex)/length(which(ITS_airLeaf_ANCOM_wide_robust2$ANCOMcat == "bioaerosol"))*100 # 4.7% of the time, foliar more abundant in a biaoerosol indicator
ITSfoliarAirBiggerIndex <- intersect(which(ITS_airLeaf_ANCOM_wide_robust2$ANCOMcat == "foliar surface"), which(ITS_airLeaf_ANCOM_wide_robust2$whichGreater == "air"))
length(ITSfoliarAirBiggerIndex)/length(which(ITS_airLeaf_ANCOM_wide_robust2$ANCOMcat == "foliar surface"))*100 #never is bioaerosol more abundant in foliar indicator

# Add this to dataframe
ITS_airLeaf_ANCOM_wide_robust2$discrepencyANCOM <- NA
ITS_airLeaf_ANCOM_wide_robust2$discrepencyANCOM[ITSbioFoliarBiggerIndex] <- "foliar more abundant in bioaerosol indicator!"
ITS_airLeaf_ANCOM_wide_robust2$discrepencyANCOM[ITSfoliarAirBiggerIndex] <- "bioaerosol more abundant in foliar indicator!"

### MAIN RESULTS CHECKING PLOTS ####
head(ITS_airLeaf_ANCOM_long_robust)
# View(ITS_airLeaf_ANCOM_long_robust)
# GET THE PROPORTION OF READS IN EACH SAMPLE THAT COME FROM BIOAEROSOL AND FOLIAR ANCOM TAXA
ITS_airLeaf_ANCOM_long_robustBySamp <- ITS_airLeaf_ANCOM_long_robust %>% 
  group_by(sampleName, ANCOMcat) %>% 
  summarise(propSampANCOMcat = sum(propASVsamp)) #propASVsamp is proportion
# View(ITS_airLeaf_ANCOM_long_robustBySamp)
colnames(ITS_airLeaf_ANCOM_long_robustBySamp)
# Add back in sample type
# Add in air or leaf, depending on what sample type it is
ITS_airLeaf_ANCOM_long_robustBySamp <- ITS_airLeaf_ANCOM_long_robustBySamp %>%
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "bioaerosol",
      str_detect(sampleName, "phyllo") ~ "foliar surface",
      TRUE ~ NA_character_
    )
  )
ITS_airLeaf_ANCOM_long_robustBySamp
# Double check
# air_ITS_105 and phyllo_ITS_43 proportion of reads from two types of ANCOM taxa
ASVs_air105andPhyllo43_ITS <- ITSairPhylloOnly.ps_ASVs[which(rownames(ITSairPhylloOnly.ps_ASVs) %in% c("air_ITS_105", "phyllo_ITS_43")),]
totalReadsair105andPhyllo43_ITS <- rowSums(ASVs_air105andPhyllo43_ITS)
# Bioaerosol ASVs
bioaerosolANCOMASVs_air105andPhyllo43_ITS <- ASVs_air105andPhyllo43_ITS[,colnames(ASVs_air105andPhyllo43_ITS) %in% bioaerosolITSANCOMASVsName_robust]
# Foliar surface ASVs
foliarANCOMASVs_air105andPhyllo43_ITS <- ASVs_air105andPhyllo43_ITS[,colnames(ASVs_air105andPhyllo43_ITS) %in% foliarITSANCOMASVsName_robust]
propBioaerosol_air105andPhyllo43_ITS <- rowSums(bioaerosolANCOMASVs_air105andPhyllo43_ITS)/totalReadsair105andPhyllo43_ITS
propBioaerosol_air105andPhyllo43_ITS
propFoliar_air105andPhyllo43_ITS <- rowSums(foliarANCOMASVs_air105andPhyllo43_ITS)/totalReadsair105andPhyllo43_ITS
propFoliar_air105andPhyllo43_ITS
# Get indices
ITSair105bioIndex <- intersect(which(ITS_airLeaf_ANCOM_long_robustBySamp$sampleName == "air_ITS_105"),which(ITS_airLeaf_ANCOM_long_robustBySamp$ANCOMcat == "bioaerosol"))
ITSair105FoliarIndex <- intersect(which(ITS_airLeaf_ANCOM_long_robustBySamp$sampleName == "air_ITS_105"),which(ITS_airLeaf_ANCOM_long_robustBySamp$ANCOMcat == "foliar surface"))
ITSfoliar43bioIndex <- intersect(which(ITS_airLeaf_ANCOM_long_robustBySamp$sampleName == "phyllo_ITS_43"),which(ITS_airLeaf_ANCOM_long_robustBySamp$ANCOMcat == "bioaerosol"))
ITSfoliar43FoliarIndex <- intersect(which(ITS_airLeaf_ANCOM_long_robustBySamp$sampleName == "phyllo_ITS_43"),which(ITS_airLeaf_ANCOM_long_robustBySamp$ANCOMcat == "foliar surface"))

# These are correct, so calculation was correct
all.equal(ITS_airLeaf_ANCOM_long_robustBySamp$propSampANCOMcat[ITSair105bioIndex],unname(propBioaerosol_air105andPhyllo43_ITS[2]))
all.equal(ITS_airLeaf_ANCOM_long_robustBySamp$propSampANCOMcat[ITSair105FoliarIndex],unname(propFoliar_air105andPhyllo43_ITS[2]))
all.equal(ITS_airLeaf_ANCOM_long_robustBySamp$propSampANCOMcat[ITSfoliar43bioIndex],unname(propBioaerosol_air105andPhyllo43_ITS[1]))
all.equal(ITS_airLeaf_ANCOM_long_robustBySamp$propSampANCOMcat[ITSfoliar43FoliarIndex],unname(propFoliar_air105andPhyllo43_ITS[1]))

# 1. For each sample, what proportion of reads are bioaerosol versus foliar surface taxa?
colnames(ITS_airLeaf_ANCOM_long_robustBySamp)

propANCOMByReadsPlot_ITS_robust <- ggplot(ITS_airLeaf_ANCOM_long_robustBySamp, aes(x = sampleType, y = propSampANCOMcat, fill=sampleType)) +
  geom_boxplot(outlier.shape = NA) +           
  scale_fill_manual(values=c("cornflowerblue", "forestgreen")) +
  geom_jitter(width = 0.15, size = 1.8, height =0) + #each sample point
  facet_wrap(~ ANCOMcat, nrow = 1)   +          # ANCOM category by sample
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(name= "Proportion bioaerosol or \nfoliar surface ANCOM taxa")
propANCOMByReadsPlot_ITS_robust
# saveRDS(propANCOMByReadsPlot_ITS_robust, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/propANCOMByReadsPlot_ITS_robust_10-8-25.rds")

##### ii. JUST STRUCTURAL ZEROS #####
ITS_AFS_zeros_df
# 1. Create ASV table with only these struct zero ANCOM ASVs in air and leaf samples. Can use this to compare how ANCOM model and actual data compare
ITS_leafAir_zero_ASVs <- unique(ITS_AFS_zeros_df$ASV_name)
length(ITS_leafAir_zero_ASVs) #7418 different ASVs
# Pull out only ASV table for significantly different ASVs (ASVs are columns!)
ITS_AFS_zerosASVtab <- ITSairPhylloOnly.ps_ASVs[,colnames(ITSairPhylloOnly.ps_ASVs) %in% ITS_leafAir_zero_ASVs]
setdiff(colnames(ITS_AFS_zerosASVtab), ITS_leafAir_zero_ASVs ) #these are the same, as expected
# View(ITS_AFS_zerosASVtab)
head(ITS_AFS_zerosASVtab)
# But which samples are air versus leaves?
#View(ITS_AFS_zerosASVtab)
ITS_AFS_zerosASVtab <- as.data.frame(ITS_AFS_zerosASVtab) %>% 
  rownames_to_column(var= "sampleName") %>% 
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "air",
      str_detect(sampleName, "phyllo") ~ "phyllo",
      TRUE ~ NA_character_
    )
  )
cbind(ITS_AFS_zerosASVtab$sampleName, ITS_AFS_zerosASVtab$sampleType) #quick check to make sure code above worked... and it did!
# One last check before calculating...
all.equal(sort(ITS_leafAir_zero_ASVs), sort(colnames(ITS_AFS_zerosASVtab)[2:(ncol(ITS_AFS_zerosASVtab)-1)])) #okay these are still the significant ones

# 2. Add in number of reads per sample (readNumberITSairFoliar_df made higher above)
#Add to working dataframe
head(ITS_AFS_zerosASVtab)
ITS_AFS_zeroASVs <- merge(ITS_AFS_zerosASVtab, readNumberITSairFoliar_df, by = "sampleName")
# View(ITS_AFS_zeroASVs)

# 3. MAKE A "LONGER" DATAFRAME AND ADD IN TAXONOMY
head(ITS_AFS_zeroASVs) 
# i. First, pivot longer so that ASVs are column 1 and Sample is second column 
head(colnames(ITS_AFS_zeroASVs)); tail(colnames(ITS_AFS_zeroASVs)) #get first and last one to put in below
ITS_AFS_zeroASVs_long <- ITS_AFS_zeroASVs %>% 
  pivot_longer(cols = ASV_10:ASV_16858, names_to= "ASV_name", values_to = "ASV_abundance")
head(ITS_AFS_zeroASVs_long)
# View(ITS_AFS_zeroASVs_long)
dim(ITS_AFS_zeroASVs_long) #1253642       5
# ii. Below adds in taxonomy with result of abundance of each ANCOM ASV in each air and phyllo sample, plus taxonomy!
ITS_AFSzeros_long <- left_join(ITS_AFS_zeroASVs_long, ITSall_8.5Kfiltered_taxTab2, by = "ASV_name")
head(ITS_AFSzeros_long)
# View(ITS_AFSzeros_long)
# iii. Add in which sample the structural zero came from
colnames(ITS_AFS_zeros_df)
ITS_AFSzeros_long <- merge(ITS_AFSzeros_long, ITS_AFS_zeros_df, by ="ASV_name")
colnames(ITS_AFSzeros_long)

# 4. Get mean prop of each ASV in bioaerosol and foliar
# i. Get each ASVs prop in each sample
ITS_AFSzeros_long <- ITS_AFSzeros_long %>% 
  mutate(propASVsamp = ASV_abundance/TotalReads) #get each ASVs proportion in each sample. Double checked that this worked!
head(ITS_AFSzeros_long)
# ii. Get each ASVs mean percent of both sample types, keeping struct. zero info too
# USE THIS ONE TO PLOT OUT EACH ASV
ITS_AFSzeros_long2 <- ITS_AFSzeros_long %>% 
  group_by(ASV_name, struct_ZeroType, sampleType) %>% 
  summarize(meanASVPctSampType = mean(propASVsamp)*100)
# View(ITS_AFSzeros_long2) #these match those I checked in Table S6, e.g. ASV_56!

# Make this wider for direct comparison
colnames(ITS_AFSzeros_long2)
ITS_AFS_zeros_wide <- ITS_AFSzeros_long2 %>% 
  pivot_wider(names_from= sampleType, values_from = meanASVPctSampType)
head(ITS_AFS_zeros_wide)
# View(ITS_AFS_zeros_wide)

# When is foliar greater than air and vice versa?
ITS_AFS_zeros_wide2 <- ITS_AFS_zeros_wide %>%
  mutate(whichGreater = case_when(
    air  > phyllo ~ "air",
    phyllo > air  ~ "foliar",
    air == phyllo ~ "ties",
    TRUE          ~ NA_character_     # any NA comparisons land here
  ))
# View(ITS_AFS_zeros_wide2)
ITSbioFoliarBiggerIndex <- intersect(which(ITS_AFS_zeros_wide2$struct_ZeroType == "bioaerosol"), which(ITS_AFS_zeros_wide2$whichGreater == "foliar"))
length(ITSbioFoliarBiggerIndex)/length(airZerosITS ) # As expected, this is always true because foliar always air str. zeros
ITSfoliarAirBiggerIndex <- intersect(which(ITS_AFS_zeros_wide2$struct_ZeroType == "foliar surfaces"), which(ITS_AFS_zeros_wide2$whichGreater == "air"))
length(ITSfoliarAirBiggerIndex)/length(phylloZerosITS) #100%, as expected and confirmed above

### STRUCTURAL ZEROS CHECKING PLOTS ####
head(ITS_AFSzeros_long)
# View(ITS_AFSzeros_long)
# GET THE PROPORTION OF READS IN EACH SAMPLE THAT COME FROM BIOAEROSOL AND FOLIAR ANCOM TAXA
ITS_AFSzeros_longBySamp <- ITS_AFSzeros_long %>% 
  group_by(sampleName, struct_ZeroType) %>% 
  summarise(propSampStrZeroType = sum(propASVsamp)) #propASVsamp is proportion
# View(ITS_AFSzeros_longBySamp)
colnames(ITS_AFSzeros_longBySamp)
# Add back in sample type
# Add in air or leaf, depending on what sample type it is
ITS_AFSzeros_longBySamp <- ITS_AFSzeros_longBySamp %>%
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "bioaerosol",
      str_detect(sampleName, "phyllo") ~ "foliar surface",
      TRUE ~ NA_character_
    )
  )
ITS_AFSzeros_longBySamp

#Double check
# air_ITS_105 and phyllo_ITS_43 proportion of reads from two types of ANCOM taxa
ASVs_air105andPhyllo43_ITS <- ITSairPhylloOnly.ps_ASVs[which(rownames(ITSairPhylloOnly.ps_ASVs) %in% c("air_ITS_105", "phyllo_ITS_43")),]
rownames(ASVs_air105andPhyllo43_ITS) #First is phyllo 43, second is air 105
totalReadsair105andPhyllo43_ITS <- rowSums(ASVs_air105andPhyllo43_ITS)
# First is phyllo 43, second is air 105
# Bioaerosol zero ASVs
BA_ASVs_air105andPhyllo43_ITS <- ASVs_air105andPhyllo43_ITS[,colnames(ASVs_air105andPhyllo43_ITS) %in% airZerosITS]
rownames(BA_ASVs_air105andPhyllo43_ITS) # First is phyllo 43, second is air 105
# Foliar surface zero ASVs
FS_ASVs_air105andPhyllo43_ITS <- ASVs_air105andPhyllo43_ITS[,colnames(ASVs_air105andPhyllo43_ITS) %in% phylloZerosITS]
propBioaerosol_air105andPhyllo43_ITS <- rowSums(BA_ASVs_air105andPhyllo43_ITS)/totalReadsair105andPhyllo43_ITS
propBioaerosol_air105andPhyllo43_ITS # First is phyllo 43, second is air 105
propFoliar_air105andPhyllo43_ITS <- rowSums(FS_ASVs_air105andPhyllo43_ITS)/totalReadsair105andPhyllo43_ITS
propFoliar_air105andPhyllo43_ITS # First is phyllo 43, second is air 105
# Get indices
# i. bioaerosol structural zeros, i.e. only in foliar for air_ITS_105
ITSair105bioIndex <- intersect(which(ITS_AFSzeros_longBySamp$sampleName == "air_ITS_105"),which(ITS_AFSzeros_longBySamp$struct_ZeroType == "bioaerosol"))
# ii. foliar structural zeros, i.e. only in bioaerosol, for air_ITS_105
ITSair105FoliarIndex <- intersect(which(ITS_AFSzeros_longBySamp$sampleName == "air_ITS_105"),which(ITS_AFSzeros_longBySamp$struct_ZeroType == "foliar surfaces"))
# iii. Bioaerosol structural zeros, i.e. only in foliar, for phyllo_ITS_43
ITSfoliar43bioIndex <- intersect(which(ITS_AFSzeros_longBySamp$sampleName == "phyllo_ITS_43"),which(ITS_AFSzeros_longBySamp$struct_ZeroType == "bioaerosol"))
# iv. foliar structural zeros, i.e. only in bioaerosol, for phyllo_ITS_43
ITSfoliar43FoliarIndex <- intersect(which(ITS_AFSzeros_longBySamp$sampleName == "phyllo_ITS_43"),which(ITS_AFSzeros_longBySamp$struct_ZeroType == "foliar surfaces"))

# # These are correct, so calculation was correct
# i. proportion bioaerosol struc zeros for air 105
all.equal(ITS_AFSzeros_longBySamp$propSampStrZeroType[ITSair105bioIndex],unname(propBioaerosol_air105andPhyllo43_ITS[2]))
# Bioaerosol type structural zero, 
all.equal(ITS_AFSzeros_longBySamp$propSampStrZeroType[ITSair105FoliarIndex],unname(propFoliar_air105andPhyllo43_ITS)[2])
all.equal(ITS_AFSzeros_longBySamp$propSampStrZeroType[ITSfoliar43bioIndex],unname(propBioaerosol_air105andPhyllo43_ITS[1]))
all.equal(ITS_AFSzeros_longBySamp$propSampStrZeroType[ITSfoliar43FoliarIndex],unname(propFoliar_air105andPhyllo43_ITS[1]))

# 1. For each sample, what proportion of reads are bioaerosol versus foliar surface taxa?
# Make copy to edit
ITS_AFSzeros_longBySampPlot_df <- ITS_AFSzeros_longBySamp
# make this a factor to re-order facets
ITS_AFSzeros_longBySampPlot_df$struct_ZeroType <- factor(ITS_AFSzeros_longBySampPlot_df$struct_ZeroType, levels = c("foliar surfaces", "bioaerosol"))
propZeroASVsByReadsPlot_ITS <- ggplot(ITS_AFSzeros_longBySampPlot_df, aes(x = sampleType, y = propSampStrZeroType, fill=sampleType)) +
  geom_boxplot(outlier.shape = NA) +           
  scale_fill_manual(values=c("cornflowerblue", "forestgreen")) +
  geom_jitter(width = 0.15, size = 1.8, height =0) + #each sample point
  facet_wrap(~ struct_ZeroType, nrow = 1, labeller = labeller(struct_ZeroType =  c("bioaerosol" = "taxa not found in bioaerosols", "foliar surfaces" = "taxa not found on foliar surfaces")))   +          #structural zero type by category
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(name= "Proportion bioaerosol or \nfoliar surface from str. zero Taxa")
propZeroASVsByReadsPlot_ITS
# saveRDS(propZeroASVsByReadsPlot_ITS, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/propZeroASVsByReadsPlot_ITS_10-6-25.rds")

########## 5. CHANGE THRESHOLDS FOR INCLUDING TAXA, AFTER ANCOMBC ##########
###### 1. PULL OUT ONLY ASV TABLE WITH ALL ANCOM ASVS ######
head(ITS_AFS_zeros_df) #has structural zeros for bioaerosols and foliar surfaces
bioaerosolITSANCOMASVsName_robust #bioaerosol 16S ANCOM$res ASV names
foliarITSANCOMASVsName_robust #foliarITS ANCOM$res names
# i. Concatenate all of these ASVs together in order to get all of the ANCOM-identifed ASVs (struc zero and main)
allITS_AFS_ANCOM_ASVs_robust <- c(ITS_AFS_zeros_df$ASV_name, bioaerosolITSANCOMASVsName_robust, foliarITSANCOMASVsName_robust)
length(allITS_AFS_ANCOM_ASVs_robust) == length(unique(allITS_AFS_ANCOM_ASVs_robust)) #no repeats, as expected
# ii. Get an ASV table with only these ANCOM ASVs
ITS_BAFS_ancom_sigASVtab_robust <- ITSairPhylloOnly.ps_ASVs[,colnames(ITSairPhylloOnly.ps_ASVs) %in% allITS_AFS_ANCOM_ASVs_robust]
# iii. Transpose so that ASVs are now rows to add in additional ASV info
ITS_BAFS_ancom_sigASVtab_t_robust <- as.data.frame(t(ITS_BAFS_ancom_sigASVtab_robust))
if (nrow(ITS_BAFS_ancom_sigASVtab_t_robust) < ncol(ITS_BAFS_ancom_sigASVtab_t_robust)) {
  print("STOP! ASVs need to be rows!")
}
# iv. Add in sum in bioaerosol and phyllosphere samples
# Sum of each ASV in foliar surface
ITS_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceTotal <- rowSums(ITS_BAFS_ancom_sigASVtab_t_robust[,which(colnames(ITS_BAFS_ancom_sigASVtab_t_robust)%in% ITSphylloSampNames)])
# Occupancy in foliar surface
ITS_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceOcc <- rowSums(ITS_BAFS_ancom_sigASVtab_t_robust[, which(colnames(ITS_BAFS_ancom_sigASVtab_t_robust)%in% ITSphylloSampNames)] > 0, na.rm = TRUE)
# Total read counts in bioaerosols
ITS_BAFS_ancom_sigASVtab_t_robust$bioaerosolTotal <- rowSums(ITS_BAFS_ancom_sigASVtab_t_robust[,which(colnames(ITS_BAFS_ancom_sigASVtab_t_robust)%in% ITSairSampNames)])
# Occupancy in bioaerosol
ITS_BAFS_ancom_sigASVtab_t_robust$bioaerosolOcc <- rowSums(ITS_BAFS_ancom_sigASVtab_t_robust[, which(colnames(ITS_BAFS_ancom_sigASVtab_t_robust)%in% ITSairSampNames)] > 0, na.rm = TRUE)

# v. Add in which kind of ANCOM taxon it is
ITS_BAFS_ancom_sigASVtab_t_robust <- ITS_BAFS_ancom_sigASVtab_t_robust %>% 
  tibble::rownames_to_column(var= "ASV_name") %>% #make rownames a column
  mutate(ANCOMcat = case_when(
    ASV_name %in% c(phylloZerosITS, bioaerosolITSANCOMASVsName_robust) ~ "bioaerosol", #bioaerosol if not in phyllo or $res bioaerosol
    ASV_name %in% c(airZerosITS, foliarITSANCOMASVsName_robust) ~ "foliar surface", #foliar if not in air or $res bioaerosol
  )) %>% 
  mutate(resType = case_when(
    ASV_name %in% c(bioaerosolITSANCOMASVsName_robust, foliarITSANCOMASVsName_robust) ~ "mainNotStructZero",
    ASV_name %in% c(phylloZerosITS, airZerosITS) ~ "structZero"
  )
  )
# vi. Perform a few checks:
# 1. Get indices
ITS_BAFtotal_ASV10_index_robust <- which(ITS_BAFS_ancom_sigASVtab_t_robust$ASV_name == "ASV_10")
ITS_BAFtotal_ASV230_index_robust <- which(ITS_BAFS_ancom_sigASVtab_t_robust$ASV_name == "ASV_230")
# 2. Is total number working?
ITS_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceTotal[ITS_BAFtotal_ASV10_index_robust] == sum(ITS_BAFS_ancom_sigASVtab_t_robust[ITS_BAFtotal_ASV10_index_robust, colnames(ITS_BAFS_ancom_sigASVtab_t_robust)%in% ITSphylloSampNames])
ITS_BAFS_ancom_sigASVtab_t_robust$bioaerosolTotal[ITS_BAFtotal_ASV230_index_robust] == sum(ITS_BAFS_ancom_sigASVtab_t_robust[ITS_BAFtotal_ASV230_index_robust, colnames(ITS_BAFS_ancom_sigASVtab_t_robust)%in% ITSairSampNames])
# 3. Is occupancy working?
ITS_BAFS_ancom_sigASVtab_t_robust$bioaerosolOcc[ITS_BAFtotal_ASV10_index_robust] == length(which(ITS_BAFS_ancom_sigASVtab_t_robust[ITS_BAFtotal_ASV10_index_robust, colnames(ITS_BAFS_ancom_sigASVtab_t_robust)%in% ITSairSampNames]>0))
# 4. Are they labeled correctly?
setdiff(ITS_BAFS_ancom_sigASVtab_t_robust$ASV_name[which(ITS_BAFS_ancom_sigASVtab_t_robust$resType == "structZero")], c(phylloZerosITS, airZerosITS)) #no differences!
# 5. Are ANCOM categories correct?
setdiff(ITS_BAFS_ancom_sigASVtab_t_robust$ASV_name[which(ITS_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "foliar surface")], c(airZerosITS, foliarITSANCOMASVsName_robust)) #no differences!
# View(ITS_BAFS_ancom_sigASVtab_t_robust)

# vii. Add in taxonomy!
head(ITSall_8.5Kfiltered_taxTab2)
ITS_BAFS_ancom_sigASVtab_t_robust <- left_join(ITS_BAFS_ancom_sigASVtab_t_robust, ITSall_8.5Kfiltered_taxTab2, by = "ASV_name") 
head(ITS_BAFS_ancom_sigASVtab_t_robust) #has species, not ASV name, as I want!

# viii. Add in the proportion of reads in each category that an ASV contributes
# Use this, made above, to get total number of reads in folair surface and air samples
totalReadsType <- readNumberITSairFoliar_df %>% 
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "air",
      str_detect(sampleName, "phyllo") ~ "phyllo",
      TRUE ~ NA_character_
    )
  ) %>% 
  group_by(sampleType) %>% 
  summarize(nReadsType = sum(TotalReads))
totalReadsType 
# sampleType nReadsType
# <chr>           <dbl>
#   1 air           4007444
# 2 phyllo        2423240
totalReadsType$nReadsType[totalReadsType$sampleType %in% "air"]

colnames(ITS_BAFS_ancom_sigASVtab_t_robust)
ITS_BAFS_ancom_sigASVtab_t_robust <- ITS_BAFS_ancom_sigASVtab_t_robust %>% 
  mutate(pctBioaerosol = bioaerosolTotal/totalReadsType$nReadsType[totalReadsType$sampleType %in% "air"]*100) %>% 
  mutate(pctFoliar = foliarSurfaceTotal/totalReadsType$nReadsType[totalReadsType$sampleType %in% "phyllo"]*100)
head(ITS_BAFS_ancom_sigASVtab_t_robust)
# DECISION AND FILTERING TIME:
# 1. To be considered diff. abundant, an ASV must be in at about 10% of respective group's samples (11 bioaerosol or 6 foliar surface) AND
# 2. Must be present at a certain threshold (see ways I've tried this below)

##### TEST BIAEROSOL ####
round(length(ITSairSampNames)/10) #11 samples 
# Found in at least 11 samples and .01% of reads
ITSairASVs_11occ.01pct_index_robust <- intersect(intersect(which(ITS_BAFS_ancom_sigASVtab_t_robust$bioaerosolOcc >=11), which(ITS_BAFS_ancom_sigASVtab_t_robust$pctBioaerosol >=0.01)), which(ITS_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "bioaerosol"))
length(ITSairASVs_11occ.01pct_index_robust) #517
ITSairASVs_11occ.05pct_index_robust <- intersect(intersect(which(ITS_BAFS_ancom_sigASVtab_t_robust$bioaerosolOcc >=11), which(ITS_BAFS_ancom_sigASVtab_t_robust$pctBioaerosol >=0.05)), which(ITS_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "bioaerosol"))
length(ITSairASVs_11occ.05pct_index_robust) #152
bioaerosolITS_0.01pct_ASVsNames <- ITS_BAFS_ancom_sigASVtab_t_robust$ASV_name[ITSairASVs_11occ.01pct_index_robust]
# View(ITS_BAFS_ancom_sigASVtab_t_robust[ITSairASVs_11occ.01pct_index_robust,])
# View(ITS_BAFS_ancom_sigASVtab_t_robust[ITSairASVs_11occ.05pct_index_robust,])
bioaerosolITS_0.05pct_ASVsNames <- ITS_BAFS_ancom_sigASVtab_t_robust$ASV_name[ITSairASVs_11occ.05pct_index_robust]
length(bioaerosolITS_0.05pct_ASVsNames) #152
##### TEST FOLIAR SURFACE ####
round(length(ITSphylloSampNames)/10) #6 samples
# Found in at least 6 samples and at least 0.01% of reads
ITS_foliarASVs_6occ.01pct_index_robust <- intersect(intersect(which(ITS_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceOcc >=6), which(ITS_BAFS_ancom_sigASVtab_t_robust$pctFoliar >=0.01)), which(ITS_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "foliar surface"))
length(ITS_foliarASVs_6occ.01pct_index_robust) #735
ITS_foliarASVs_6occ.05pct_index_robust <- intersect(intersect(which(ITS_BAFS_ancom_sigASVtab_t_robust$foliarSurfaceOcc >=6), which(ITS_BAFS_ancom_sigASVtab_t_robust$pctFoliar >=0.05)), which(ITS_BAFS_ancom_sigASVtab_t_robust$ANCOMcat == "foliar surface"))
length(ITS_foliarASVs_6occ.05pct_index_robust) #243
foliarITS_0.05pct_ASVsNames <- ITS_BAFS_ancom_sigASVtab_t_robust$ASV_name[ITS_foliarASVs_6occ.05pct_index_robust]
length(foliarITS_0.05pct_ASVsNames) #243
# View(ITS_BAFS_ancom_sigASVtab_t_robust[ITS_foliarASVs_6occ.05pct_index_robust,])
# View(ITS_BAFS_ancom_sigASVtab_t_robust)

###### PLOT AND EXPLORE NEW RESULTS! #####
# 1. Make a new dataframe that just has these new, differentially abundant taxa by concatenating indices from above
ITS_ANCOMall_.05pct_df <- ITS_BAFS_ancom_sigASVtab_t_robust[c(ITSairASVs_11occ.05pct_index_robust, ITS_foliarASVs_6occ.05pct_index_robust),]
setdiff(ITS_ANCOMall_.05pct_df$ASV_name, c(bioaerosolITS_0.05pct_ASVsNames, foliarITS_0.05pct_ASVsNames)) #these are the same!
# View(ITS_ANCOMall_.05pct_df)
colnames(ITS_ANCOMall_.05pct_df)
# A few checks!
# all look good!
sort(unique(ITS_ANCOMall_.05pct_df$pctBioaerosol[which(ITS_ANCOMall_.05pct_df$ANCOMcat == "bioaerosol")])) #minimum looks good!
sort(unique(ITS_ANCOMall_.05pct_df$pctFoliar[which(ITS_ANCOMall_.05pct_df$ANCOMcat == "foliar surface")])) #minimum looks good!

# View(ITS_ANCOMall_.05pct_df)
# Saved October 11, 2025:
# saveRDS(ITS_ANCOMall_.05pct_df, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_ANCOMall_.05pct_df.rds")

# 2. Get proportion of each sample that is foliar surface or bioaerosols ANCOM taxa
# View(ITS_ANCOMall_.05pct_df)
colnames(ITS_ANCOMall_.05pct_df)
# i. Get each ASV's abundance in each sample
ITS_ANCOMall_.05pct_df_long <- ITS_ANCOMall_.05pct_df %>% 
  pivot_longer(cols = phyllo_ITS_1:air_ITS_99, names_to = "sampleName", values_to = "ASV_abundance")
# ii. Grouping by sampName and ANCOMtype, get sum of ANCOM type in each sample
head(ITS_ANCOMall_.05pct_df_long)
ITS_nANCOMcatSamp <- ITS_ANCOMall_.05pct_df_long %>% 
  select(c(ASV_name, ANCOMcat, sampleName, ASV_abundance)) %>% 
  group_by(ANCOMcat, sampleName) %>% 
  summarize(nANCOMsamp = sum(ASV_abundance))
head(ITS_nANCOMcatSamp)
# Check for a random sample to ensure it's working. Get index where samp is 112 and ANCOMcat is bioaerosol. Then pull out ASV_abundance and sum it. True so this is working!
sum(ITS_ANCOMall_.05pct_df_long$ASV_abundance[intersect(which(ITS_ANCOMall_.05pct_df_long$sampleName == "air_ITS_112"), which(ITS_ANCOMall_.05pct_df_long$ANCOMcat == "bioaerosol"))]) == 
  ITS_nANCOMcatSamp$nANCOMsamp[intersect(which(ITS_nANCOMcatSamp$ANCOMcat == "bioaerosol"), which(ITS_nANCOMcatSamp$sampleName == "air_ITS_112"))]
# iii. Add in with read number data and then divide by this for proportion
head(readNumberITSairFoliar_df) #made above, is total number of reads per sample
ITS_nANCOMcatSamp2 <- left_join(ITS_nANCOMcatSamp, readNumberITSairFoliar_df, by = "sampleName")
head(ITS_nANCOMcatSamp2)
# Add in prop and sample Type
ITS_nANCOMcatSamp2 <- ITS_nANCOMcatSamp2 %>% 
  mutate(propANCOMcatSamp = nANCOMsamp/TotalReads) %>% 
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "air",
      str_detect(sampleName, "phyllo") ~ "phyllo",
      TRUE ~ NA_character_
    )
  )

# 3. Plot the proportion of eaach sample that is foliar surface or bioaerosols!
head(ITS_nANCOMcatSamp2)
ITS_propANCOM_filt_Reads_plot <- ggplot(ITS_nANCOMcatSamp2, aes(x = sampleType, y = propANCOMcatSamp, fill=sampleType)) +
  geom_boxplot(outlier.shape = NA) +           
  scale_fill_manual(values=c("cornflowerblue", "forestgreen")) +
  geom_jitter(width = 0.15, size = 1.8, height =0) + #each sample point
  facet_wrap(~ ANCOMcat, nrow = 1)   +          # ANCOM category by sample
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(name= "Proportion bioaerosol or \nfoliar surface ANCOM taxa") +
  ggtitle("ITS after filtering ASVs: proportion of each sample's reads by ANCOM type")
ITS_propANCOM_filt_Reads_plot

# saveRDS(ITS_propANCOM_filt_Reads_plot, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_propANCOM_filt_Reads_plot_10-11-25.rds")

head(ITS_ANCOMall_.05pct_df_long)
# i. Add in number of reads per sample
head(readNumberITSairFoliar_df) #made above, is total number of reads per sample
ITS_ANCOMall_.05pct_df_long_2 <- left_join(ITS_ANCOMall_.05pct_df_long, readNumberITSairFoliar_df, by = "sampleName")
head(ITS_ANCOMall_.05pct_df_long_2)
# ii. Add a column to get proportion of each ASV in each sample
# Add in prop and sample Type
ITS_ANCOMall_.05pct_df_long_2 <- ITS_ANCOMall_.05pct_df_long_2 %>% 
  mutate(propASVPerSamp = ASV_abundance/TotalReads) %>% 
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "air",
      str_detect(sampleName, "phyllo") ~ "phyllo",
      TRUE ~ NA_character_
    )
  )
head(ITS_ANCOMall_.05pct_df_long_2)

# 2. Get proportion of various ANCOM ASVs in foliar surface or bioaerosol reads
# This dataframe, made above, has the ASV abundance of each ANCOM-ASV in each sample 

### NOTE: CODE BELOW IS CORRECT BUT IS COMMENTED OUT BECAUSE ADDS SO MUCH TO THE
# MEMORY. ALSO, I LOOKED OVER THESE PLOTS AND ALL SEEMS GOOD!

# iii. Make a list of dataframes. First, divide these into types of ANCOM, then orders
# 1. Sort into bioaerosols and foliar surface
ITS_ANCOMASVs_ANCOMcatlist <- ITS_ANCOMall_.05pct_df_long_2 %>%
  group_by(ANCOMcat) %>%
  group_split(.keep = TRUE)
length(ITS_ANCOMASVs_ANCOMcatlist)
# View(ITS_ANCOMASVs_ANCOMcatlist[[1]])

# Name each list element by ANCOMtype
names(ITS_ANCOMASVs_ANCOMcatlist) <- character(length(ITS_ANCOMASVs_ANCOMcatlist)) #kind of pre-allocate this
for (i in seq_along(ITS_ANCOMASVs_ANCOMcatlist)) {
  names(ITS_ANCOMASVs_ANCOMcatlist)[i] <- paste0(unique(ITS_ANCOMASVs_ANCOMcatlist[[i]]$ANCOMcat), "ANCOM")
}

# 2. Divide these two lists by order and name them according to order
# Bioaerosols
ITS_ANCOMorders_BA_list <- ITS_ANCOMASVs_ANCOMcatlist$bioaerosolANCOM %>%
  group_by(Order) %>%
  group_split(.keep = TRUE)
names(ITS_ANCOMorders_BA_list) <- character(length(ITS_ANCOMorders_BA_list)) #kind of pre-allocate this
for (i in seq_along(ITS_ANCOMorders_BA_list)) {
  names(ITS_ANCOMorders_BA_list)[i] <- paste0(unique(ITS_ANCOMorders_BA_list[[i]]$Order), "BA_ANCOM")
}
# Foliar surfaces
ITS_ANCOMorders_FS_list <- ITS_ANCOMASVs_ANCOMcatlist$`foliar surfaceANCOM` %>%
  group_by(Order) %>%
  group_split(.keep = TRUE)
names(ITS_ANCOMorders_FS_list) <- character(length(ITS_ANCOMorders_FS_list)) #kind of pre-allocate this
for (i in seq_along(ITS_ANCOMorders_FS_list)) {
  names(ITS_ANCOMorders_FS_list)[i] <- paste0(unique(ITS_ANCOMorders_FS_list[[i]]$Order), "FS_ANCOM")
}

# 3. Make lists to be nested inside orders lists
# NOW! First [[]] is Order, second is specific ASV
# BIOAEROSOLS
ITS_BA_orders_nestedASVsList <- lapply(ITS_ANCOMorders_BA_list, function(df) {
  df %>% tidyr::nest(.by = ASV_name) %>%
    { setNames(.$data, .$ASV_name) }
})
names(ITS_BA_orders_nestedASVsList)[1] #e.g. AgaricalesBA_ANCOM
names(ITS_BA_orders_nestedASVsList[[1]]) #shows ASVs within AgaricalesBA_ANCOM
# FOLIAR SURFACES
ITS_FS_orders_nestedASVsList <- lapply(ITS_ANCOMorders_FS_list, function(df) {
  df %>% tidyr::nest(.by = ASV_name) %>%
    { setNames(.$data, .$ASV_name) }
})
names(ITS_FS_orders_nestedASVsList)[2] #e.g. CapnodialesFS_ANCOM
names(ITS_FS_orders_nestedASVsList[[2]]) #shows ASVs within CapnodialesFS_ANCOM
length(ITS_FS_orders_nestedASVsList[[2]]) #49

# Just look at a few orders within each to check (don't do all of them since there are so many!)
# Agaricales bioaerosol ANCOMs
Agaricales_BA_PlotsList <- vector(mode = "list", length = length(ITS_BA_orders_nestedASVsList[[1]]))
for (i in seq_along(ITS_BA_orders_nestedASVsList[[1]])) {
  df <- ITS_BA_orders_nestedASVsList[[1]][[i]]
  ttl <- paste0("ITS_BA_",names(ITS_BA_orders_nestedASVsList[[1]])[i], "_",
                paste(unique(df$Family), collapse = ", "))

  Agaricales_BA_PlotsList[[i]] <- try(
    ggplot(df, aes(x = sampleType, y = propASVPerSamp, fill = sampleType)) +
      geom_boxplot(outlier.shape = NA) +
      scale_fill_manual(values = c("cornflowerblue", "forestgreen")) +
      geom_jitter(width = 0.15, size = 1.8, height = 0) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      scale_y_continuous(name = "Proportion in each bioaerosol or \nfoliar surface sample") +
      ggtitle(ttl),
    silent = TRUE
  )
}
# Grayed out since a lot of plots
# Agaricales_BA_PlotsList

# Polyporales B.A. indicators
names(ITS_BA_orders_nestedASVsList)[9]
Polyporales_BA_PlotsList <- vector(mode = "list", length = length(ITS_BA_orders_nestedASVsList[[9]]))
for (i in seq_along(ITS_BA_orders_nestedASVsList[[9]])) {
  df <- ITS_BA_orders_nestedASVsList[[9]][[i]]
  ttl <- paste0("ITS_BA_",names(ITS_BA_orders_nestedASVsList[[9]])[i], "_",
                paste(unique(df$Family), collapse = ", "))

  Polyporales_BA_PlotsList[[i]] <- try(
    ggplot(df, aes(x = sampleType, y = propASVPerSamp, fill = sampleType)) +
      geom_boxplot(outlier.shape = NA) +
      scale_fill_manual(values = c("cornflowerblue", "forestgreen")) +
      geom_jitter(width = 0.15, size = 1.8, height = 0) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      scale_y_continuous(name = "Proportion in each bioaerosol or \nfoliar surface sample") +
      ggtitle(ttl),
    silent = TRUE
  )
}
# Grayed out to not run automatically since so many!!
# Polyporales_BA_PlotsList
length(Polyporales_BA_PlotsList) #84!!!
# Capnodiales F.S. indicators
names(ITS_FS_orders_nestedASVsList)[2]
Capnodiales_FS_PlotsList <- vector(mode = "list", length = length(ITS_FS_orders_nestedASVsList[[2]]))
for (i in seq_along(ITS_FS_orders_nestedASVsList[[2]])) {
  df <- ITS_FS_orders_nestedASVsList[[2]][[i]]
  ttl <- paste0("ITS_FS_",names(ITS_FS_orders_nestedASVsList[[2]])[i], "_",
                paste(unique(df$Family), collapse = ", "))

  Capnodiales_FS_PlotsList[[i]] <- try(
    ggplot(df, aes(x = sampleType, y = propASVPerSamp, fill = sampleType)) +
      geom_boxplot(outlier.shape = NA) +
      scale_fill_manual(values = c("cornflowerblue", "forestgreen")) +
      geom_jitter(width = 0.15, size = 1.8, height = 0) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      scale_y_continuous(name = "Proportion in each bioaerosol or \nfoliar surface sample") +
      ggtitle(ttl),
    silent = TRUE
  )
}
length(Capnodiales_FS_PlotsList) #49
# Grayed out to not automatically plot since there's so many!
