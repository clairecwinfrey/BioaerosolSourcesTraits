# fullITS_EDA_part1_Sept.R
# Air, soil, and phyllosphere--ITS Exploratory Data Analysis and Data Clean Up
# (re-started) started September 16, 2025

# SCRIPT ORDER: bioinformatics_final_ITS.R --> THIS SCRIPT --> fullITS_EDArarefied_Part2.R

# Description from original fullITS_EDA_part1.R

# This script processed raw data (after dada2/bioinformatics in script called "bioinformatics_final_ITS.R", 
# which was run on the Fierer Lab "microbe" server), to generate a final, cleaned and rarefied dataset to be used in all down-
# stream analyses on fungal data. First, I removed anything which wasn't classified as Kingdom Fungi and which wasn't at least
# classified to phylum level (i.e. phylum = NA). Second, I removed air samples not sampled at least 20 hours. Third, I explored
# ASVs that could be contaminants in the air data (i.e., those that showed up in at least 19/39 of the AIR controls/blanks), but
# I did NOT remove any ASVs from the dataset because none matched this criterion Fourth, I rarefied all (i.e., air, soil, and 
# phyllosphere) at 8,500 reads.

# The end result of this script is:

# This script's bacterial counterpart is called full16S_EDA_part1_Sept.R

#######################################
#             SCRIPT SET UP
#######################################

# FUNCTIONS CREATED IN THIS SCRIPT:
# 1. Function to get the taxonomy table out of phyloseq:
# (Inspired by similar function at https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# 2. Function to get ASV table out of phyloseq so that we can #View it better
# (Inspired by similar function at https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  return(as.data.frame(ASVTable))
}

#######################################################################################
# Set working directory
list.files()
# Read in libraries
library("phyloseq"); packageVersion("phyloseq") #‘1.46.0’
library("tidyverse"); packageVersion("tidyverse") #‘2.0.0’
library("vegan"); packageVersion("vegan") #2.7.1’
library("gridExtra"); packageVersion("gridExtra")  #‘2.3’  
# library("xlsx"); packageVersion("xlsx")  #‘2.3’  
library(Polychrome); packageVersion("Polychrome")  #‘1.5.1’

##################################################################################
# I. SET-UP, DATA CLEANING, RAREFACTION, AND FIRST TAXONOMIC & ORDINATION PLOTS
##################################################################################
allITS_seqtab_wTax_mctoolsr <- read.table("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/bioinformaticsOutput/correct18April2023/allITS_seqtab_wTax_mctoolsr.txt", header=TRUE)
str(allITS_seqtab_wTax_mctoolsr)
head(allITS_seqtab_wTax_mctoolsr)
##View(allITS_seqtab_wTax_mctoolsr)
# First column is X.ASV_ID. Change to ASV_ID
colnames(allITS_seqtab_wTax_mctoolsr)[1] <- "ASV_ID"
# All of the sample names that start with "X" are soils. Edit their names here:
colnames(allITS_seqtab_wTax_mctoolsr) <- gsub(colnames(allITS_seqtab_wTax_mctoolsr), pattern = "X", replacement = "soil_ITS_") 
head(allITS_seqtab_wTax_mctoolsr)
rownames(allITS_seqtab_wTax_mctoolsr) #numbers

allITS_seqtab <- read.table("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/bioinformaticsOutput/correct18April2023/allITS_seqtab_final.txt", header=T)
dim(allITS_seqtab) #16866   413
# View(allITS_seqtab)

#### LOAD FILES IN FOR PHYLOSEQ #####
# 1. OTU table: ASVs/OTUs as rows, samples as columns. This is "seqtab"
allITS_otu_mat <- allITS_seqtab
allITS_seqtab$X #X is the ASV names
allITS_otu_mat <- allITS_otu_mat %>%
  tibble::column_to_rownames("X") #make it so that "X"(which is column with ASV names is the rownames columnn!)
colnames(allITS_otu_mat) <- gsub(colnames(allITS_otu_mat), pattern = "X", replacement = "soil_ITS_") 
colnames(allITS_otu_mat) #colnames are sample names
rownames(allITS_otu_mat) #ASV names
##View(allITS_otu_mat)

# 2. OTU taxonomy
# This is basically needs to be the same as the first and last columns of allITS_seqtab_wTax_mctoolsr
# So none of the guts that are essentially the ASV table
dim(allITS_seqtab_wTax_mctoolsr) #16866   415
colnames(allITS_seqtab_wTax_mctoolsr) #samples 
head(allITS_seqtab_wTax_mctoolsr)
head(allITS_seqtab_wTax_mctoolsr[,1]) #THE ASV names!
head(allITS_seqtab_wTax_mctoolsr[,ncol(allITS_seqtab_wTax_mctoolsr)]) #"s__thailandicum" NA                "s__furfuracea"   NA                NA                NA 
head(allITS_seqtab_wTax_mctoolsr[,(ncol(allITS_seqtab_wTax_mctoolsr)-1)]) #the last and second to last columns are the taxonomy 
# The last column  sometimes has species info
# Bind together ASV IS names and taxonomy  
step1 <- cbind(allITS_seqtab_wTax_mctoolsr[,1], allITS_seqtab_wTax_mctoolsr[,(ncol(allITS_seqtab_wTax_mctoolsr)-1)]) 
head(step1)
# Bind step 1 with species identification, if available 
step2 <- cbind(step1, allITS_seqtab_wTax_mctoolsr[,ncol(allITS_seqtab_wTax_mctoolsr)])
head(step2) # columns are now "#ASV_ID", "taxonomy", and "Species"
colnames(step2) <- c("ASV_ID", "taxonomy","Species") #make ASV_ID and taxonomy the column names
head(step2)
dim(step2) #16866, 3

# Make new columns based on the ";" separator
require(tidyr)
allITS_tax_sep <- separate(as.data.frame(step2), col = taxonomy, into= c("Kingdom", "Phylum", "Class", 
                                                                         "Order", "Family", "Genus", "ASV_ID"),
                           sep = ";")
head(allITS_tax_sep)
colnames(allITS_tax_sep) #"Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "ASV_ID"  "Species"
# Continuing tidying up taxonomy table
# #View(allITS_tax_sep)
rownames(allITS_tax_sep)
allITS_tax_sep <- allITS_tax_sep %>% #make ASV_ID column the rownames to make it compatible with phyloseq
  tibble::column_to_rownames("ASV_ID")
head(allITS_tax_sep)

# 3. Sample metadata--
all_ITS_Metadata <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/aeroAllMetadat_Apr20_2023.csv", row.names=1) #made in ITS_metadataAndMappingFileSetUp.R
head(all_ITS_Metadata)
# View(all_ITS_Metadata)
colnames(all_ITS_Metadata) #first column in sampleNumber, EU, TransectMeter, DateSetOut
rownames(all_ITS_Metadata) #numbers 
airSamplesALL_df <- all_ITS_Metadata %>%
  tibble::column_to_rownames("sampleNumber")
head(airSamplesALL_df) #sample names are NOW rows
# View(airSamplesALL_df) 

# Which samples differ betweeen OTU table and metadata table? Clean these up!
setdiff(sort(colnames(allITS_otu_mat)), sort(rownames(airSamplesALL_df))) # differences below
#  "air_ITS_PCR_NTC_2_H12" "ITS_NTC_P1"            "ITS_NTC_P2"            "ITS_NTC_P3" 
air_ITS_PCR_NTC_2_H12index <- which(colnames(allITS_otu_mat) == "air_ITS_PCR_NTC_2_H12")
sum(allITS_otu_mat[,air_ITS_PCR_NTC_2_H12index]) #since this only has only read, it can safely be dropped!
ITS_NTC_P1index <- which(colnames(allITS_otu_mat) == "ITS_NTC_P1")
sum(allITS_otu_mat[,ITS_NTC_P1index]) #since this also only has only read, it can safely be dropped!
ITS_NTC_P2index <- which(colnames(allITS_otu_mat) == "ITS_NTC_P2")
sum(allITS_otu_mat[,ITS_NTC_P2index]) #since this also only has two reads, it can safely be dropped!
ITS_NTC_P3index <- which(colnames(allITS_otu_mat) == "ITS_NTC_P3")
sum(allITS_otu_mat[,ITS_NTC_P3index]) #since this also only has 9 reads, it can safely be dropped!
allITS_otu_mat <- allITS_otu_mat[,-c(air_ITS_PCR_NTC_2_H12index, ITS_NTC_P1index, ITS_NTC_P2index, ITS_NTC_P3index)]
setdiff(sort(colnames(allITS_otu_mat)), sort(rownames(airSamplesALL_df))) #Cool, names are now the same
setdiff(sort(rownames(airSamplesALL_df)), sort(colnames(allITS_otu_mat))) #as expected, "air_ITS_21" and "air_ITS_PCR_NTC_1" are in metadata but not in otu_table,m
# so these need to be dropped from metadata
index21andNTC1 <- c(which(rownames(airSamplesALL_df) == "air_ITS_21"), which(rownames(airSamplesALL_df) == "air_ITS_PCR_NTC_1")) #140, 137
airSamples_df <- airSamplesALL_df[-c(index21andNTC1),] #remove air_ITS_21" and "air_ITS_PCR_NTC_1" 
# airSamples_df <- airSamples_df[,c(1:9, 11:13)] #drop some extra columns
# View(airSamples_df)
setdiff(sort(rownames(airSamples_df)), sort(colnames(allITS_otu_mat))) #now the same!

# CLEAN UP allITS_tax_sep by removing pesky "k__", "p__", etc. throughout
remove_substrings <- function(x) { #make a lil function
  gsub("k__|p__|c__|o__|f__|g__|s__", "", x)
}
# Apply it to each column of allITS_tax_sep
allITS_tax_sepFixed <- as.data.frame(lapply(allITS_tax_sep, remove_substrings))
head(allITS_tax_sepFixed) #looks good!
rownames(allITS_tax_sepFixed) <- rownames(allITS_tax_sep) #add the rownames back!

# Transform otu table and taxonomy tables into matrices
allITS_otu_mat <- as.matrix(allITS_otu_mat)
allITS_tax_mat <- as.matrix(allITS_tax_sepFixed)
# View(allITS_tax_mat)

# Transform to phyloseq objects
allITS_OTU <- otu_table(allITS_otu_mat, taxa_are_rows = TRUE)
allITS_TAX <- tax_table(allITS_tax_mat)
allITS_samples <- sample_data(airSamplesALL_df)
allITS_raw.ps <- phyloseq(allITS_OTU, allITS_TAX, allITS_samples)
allITS_raw.ps #16866 taxa and 408 samples
colnames(otu_table(allITS_raw.ps)) #sample names as expected!
sample_names(allITS_raw.ps) 
rank_names(allITS_raw.ps)
sample_variables(allITS_raw.ps)

# saved April 23, 2023
# save(allITS_raw.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allITS_raw_phyloseq") 
# saved on server June 16, 2024 (with cleaned up taxonomy table!)
# save(allITS_raw.ps, file ="~/SRS_aeromicrobiome_2022/RobjectsSaved/allITS_raw_phyloseq")

# saved again on own computer Sept. 16, 2025 for double-chcecking purposes: 
# save(allITS_raw.ps, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/allITS_raw_phyloseq")

##################################################################
# DATA CLEANING PART I
# Check for anything non-fungal and ASVs not classified to at least phylum level
# and remove if present-
##################################################################
# Get taxonomy table out of phyloseq:
allITS_taxTable <- taxtable_outta_ps(allITS_raw.ps)
allITS_ASVs <-  ASVs_outta_ps(allITS_raw.ps)
allITS_ASVs$ASVabundance <- rowSums(allITS_ASVs)
# View(allITS_taxTable) #raw has 16,866 ASVs
# View(allITS_ASVs)

# Removals of NAs and non-fungi
allITS_tax_noNAs <- allITS_taxTable %>% 
  filter(Kingdom != "NA") %>% #Remove ASVs where kingdom is unknown ("NA" in first column) 
  filter(Phylum != "NA") #Remove ASVs where phylum is unknown ("NA" in first column) 
# #View(allITS_tax_noNAs) 
# Are any of these remaining a phylum other than fungi?
kingdomFung <- allITS_tax_noNAs %>% 
  filter(Kingdom == "Fungi")
dim(kingdomFung)[1] == dim(allITS_tax_noNAs)[1] #no, there were no kingdoms other than fungi

# Filtered out a total of 2,946 ASVS, for a resulting 13,920 ASVs.
dim(allITS_taxTable)[1] - dim(allITS_tax_noNAs)[1]

# We now have to trim the OTU table, because it still has the taxa in it that were filtered out
# above. It currently has 16,866 ASVs, whereas allITS_tax_noNAs has 13,920 ASVs. We'll want it this 
# way so that we re-make the phyloseq object and then rarefy correctly. 
dim(otu_table(allITS_raw.ps)) #16,866 ASVs

allITS_ASVsToKeep <- rownames(allITS_tax_noNAs) #grab names of ASVs to keep (that were filtered out of tax
# table above)

allITS_noNAs.ps <- prune_taxa(allITS_ASVsToKeep, allITS_raw.ps) #new phyloseq object with just the stuff we want!
# Now 13920 taxa 
# Check to make sure new phyloseq object looks as expected
unique(rownames(allITS_tax_noNAs) == rownames(otu_table(allITS_noNAs.ps))) #yep! 

# TAKING A LOOK AT THE DATA
# How are number of reads distributed across the samples?
colnames(otu_table(allITS_noNAs.ps))
allITS_seqspersample <- colSums(otu_table(allITS_noNAs.ps))
allITS_SeqNumberPlot <-barplot(allITS_seqspersample, main="All ITS: Total Sequences Per Sample", ylab= "Number of Sequences", xlab= "Sample Number")

##################################################################
# DATA CLEANING PART II
# (Add more metadata (controls versus samples) and remove samples
# not sampled at least 20 hours)
##################################################################
allITS_noNAs.ps
sampDat <- as.data.frame(as.matrix(sample_data(allITS_noNAs.ps)))
unique(sampDat$sampleType)

# Make a new column to distinguish all controls from all samples
sampDat$isControl <- rep(NA, nrow(sampDat))
for (i in 1:nrow(sampDat)){
  if (grepl(pattern= "ontrol", sampDat$sampleType[i])) {
    sampDat$isControl[i] <- "control" 
  } else {
    sampDat$isControl[i] <- "sample"
  }
}

# Now make the others controls too
for (i in 1:nrow(sampDat)){
  if (grepl(pattern= "kitome", sampDat$sampleType[i])) {
    sampDat$isControl[i] <- "control" 
  } 
  if (grepl(pattern= "Blank", sampDat$sampleType[i])) {
    sampDat$isControl[i] <- "control"
  }
}

cbind(sampDat$sampleType, sampDat$isControl) #controls look correctly defined

# Make allITS_noNAs.ps have this new metadata -- allITS_noNAs2
allITS_noNAs2.ps <- phyloseq(tax_table(allITS_noNAs.ps), otu_table(allITS_noNAs.ps), sample_data(sampDat))
allITS_noNAs2.ps
# 13920 taxa and 408 samples

# Saved Sept 5, 2023
#save(allITS_noNAs2.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allITS_noNAs2phyloseq") 
# saved on server June 17, 2024 
# save(allITS_noNAs2.ps, file ="~/SRS_aeromicrobiome_2022/RobjectsSaved/allITS_noNAs2phyloseq") #saved on server June 17, 2024

# Saved again on own computer Sept. 16, 2025 for double-chcecking purposes: 
# save(allITS_noNAs2.ps, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/allITS_noNAs2_phyloseq")

#### Drop air samples that ran less than 20 hours (air)
colnames(sample_data(allITS_noNAs2.ps))
airSamps <- which(sample_data(allITS_noNAs2.ps)[,21] == "air") 
# View(as.data.frame(as.matrix(sample_data(allITS_noNAs2.ps)[airSamps,])))
greater20h <- which(sample_data(allITS_noNAs2.ps)[,15] >= 20.00) #column 15 is sampledRunTime
# Double check by viewing that the time was no less than 20 hours. Looks good below!
# View(as.data.frame(as.matrix(sample_data(allITS_noNAs2.ps)[intersect(airSamps, greater20h),])))
dim(as.data.frame(as.matrix(sample_data(allITS_noNAs2.ps)[intersect(airSamps, greater20h),]))) 
# This shows that 110 out of the 126 air samples (not counting controls) sampled for at least 20 hours.
# In fact, the lowest in this group (air_ITS_19) sampled for 23.410 hours!

# Find which air samples to drop by finding those air samples that did not sample at least 20 hours
badAir <- setdiff(airSamps, greater20h)
# View(as.data.frame(as.matrix(sample_data(allITS_noNAs2.ps)[badAir,])))  #double check-- looks good!
rownames(sample_data(allITS_noNAs2.ps)[-badAir,])
# View(as.data.frame(as.matrix(sample_data(allITS_noNAs2.ps)[-badAir,])))  #double check-- looks good!
goodNames <- rownames(sample_data(allITS_noNAs2.ps)[-badAir,])
allITS_20hrs.ps <- prune_samples(goodNames, allITS_noNAs2.ps)
# View(as.data.frame(as.matrix(sample_data(allITS_20hrs.ps)))) #looks good!

# What are the names of all of the air samples (including ALL controls) that were sampled at least 20 hours?
ITS_airAllIndex <- grep(pattern= "air", x= rownames(as.data.frame(as.matrix(sample_data(allITS_20hrs.ps))))) #get all of the rownames that contain "air"
rownames(as.data.frame(as.matrix(sample_data(allITS_20hrs.ps))))[ITS_airAllIndex] #yes, this pulled out everything with air in the title
length(ITS_airAllIndex) #147 (controls/blanks and samples )

# saved this September 30, 2025
# saveRDS(allITS_20hrs.ps, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allITS_20hrs_phyloseq.rds")

##################################################################
# DATA CLEANING PART III-- more investigation of controls and 
# (BEFORE contaminant removal)
##################################################################
# Get sample names that have "air" in them 
airITS_OTUs20h <- as.data.frame(as.matrix(otu_table(allITS_20hrs.ps)[,ITS_airAllIndex]))
airITS_OTUs20h <- airITS_OTUs20h[which(rowSums(airITS_OTUs20h) > 0),] #keep only those rows (i.e., ASVs) for which there is at least one copy of the ASV found in the air
# View(airITS_OTUs20h)
# Double check
sort(unique(rowSums(airITS_OTUs20h))) #none of these were zeros, as desired!
# Get taxa table 
tax_allITS_20hrs <- as.data.frame(as.matrix(tax_table(allITS_20hrs.ps)))
# merge ASV taxonomy onto the end of the trimmed air only OTU table
ITS_airTaxOTUs <- merge(airITS_OTUs20h, tax_allITS_20hrs, all=FALSE, by=0)
# written and saved July 21, 2023
# write.csv(ITS_airTaxOTUs, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_airTaxOTUs.csv") #on personal computer ofc

colnames(ITS_airTaxOTUs) #columns 2 through 148 are samples/controls
# Add in read count for each sample and make this the last row
airITS_TaxOTUsSums <- rbind(ITS_airTaxOTUs,c(NA,colSums(ITS_airTaxOTUs[,2:148]), rep(NA, 7))) #2:148 is just the columns representing samples. NA and rep(NA, 6) just so rownames and taxonomy cols have empty row below it and R doesnt get upset
# View(airITS_TaxOTUsSums)

# Get names and indices of different types of air samples
washBufferControlNames <- rownames(sample_data(allITS_20hrs.ps))[which(sample_data(allITS_20hrs.ps)$sampleType == "washBufferControl")]
washBufferControlIndex <- which(colnames(airITS_TaxOTUsSums) %in% washBufferControlNames == TRUE)
fieldControlNames <- rownames(sample_data(allITS_20hrs.ps))[which(sample_data(allITS_20hrs.ps)$sampleType == "fieldControl")]
fieldControlIndex <- which(colnames(airITS_TaxOTUsSums) %in% fieldControlNames == TRUE)
kitomeAirControlNames <- rownames(sample_data(allITS_20hrs.ps))[which(sample_data(allITS_20hrs.ps)$sampleType == "kitome")]
kitomeIndex <- which(colnames(airITS_TaxOTUsSums) %in% kitomeAirControlNames == TRUE)
PCRAirControlNames <- rownames(sample_data(allITS_20hrs.ps))[which(sample_data(allITS_20hrs.ps)$sampleType == "PCRcontrol")]
PCR_Index <- which(colnames(airITS_TaxOTUsSums) %in% PCRAirControlNames == TRUE)
airNames <- rownames(sample_data(allITS_20hrs.ps))[which(sample_data(allITS_20hrs.ps)$sampleType == "air")]
airSampsIndex <- which(colnames(airITS_TaxOTUsSums) %in% airNames == TRUE)
# Make a big index based on all of the control samples 
controlsIndex <- c(fieldControlIndex, washBufferControlIndex, kitomeIndex, PCR_Index)
# get names of all af these to use later
ITS_airControlAllNames <- c(washBufferControlNames, fieldControlNames, kitomeAirControlNames, PCRAirControlNames) #37 of 'em

# Re-order for easier comparison so that air samples, then control samples, then taxonomy.
colnames(airITS_TaxOTUsSums)
ITS_airTaxOTUsOrdered <- airITS_TaxOTUsSums[, c(1, airSampsIndex, controlsIndex, 149:155)] #149:155 are taxonomy categories
# View(ITS_airTaxOTUsOrdered)
ITS_airTaxOTUsOrderedReNamed <- ITS_airTaxOTUsOrdered #make another one to rename based on control type

# Rename for easier comparison
colnames(ITS_airTaxOTUsOrderedReNamed)[which(colnames(ITS_airTaxOTUsOrderedReNamed) %in% washBufferControlNames == TRUE)] <- "bufferCtrl"
colnames(ITS_airTaxOTUsOrderedReNamed)[which(colnames(ITS_airTaxOTUsOrderedReNamed) %in% fieldControlNames == TRUE)] <- "fieldCtrl"
colnames(ITS_airTaxOTUsOrderedReNamed)[which(colnames(ITS_airTaxOTUsOrderedReNamed) %in% kitomeAirControlNames == TRUE)] <- "kitome"
colnames(ITS_airTaxOTUsOrderedReNamed)[which(colnames(ITS_airTaxOTUsOrderedReNamed) %in% PCRAirControlNames == TRUE)] <- "PCRcontrol"
# View(ITS_airTaxOTUsOrderedReNamed)

# Just air + associated controls after this:
airITS_Index_all <- which(sample_data(allITS_20hrs.ps)$sampleType== "air")
fieldControlITS_Index_all <- which(sample_data(allITS_20hrs.ps)$sampleType== "fieldControl")
washBufferControlITS_Index_all <- which(sample_data(allITS_20hrs.ps)$sampleType== "washBufferControl")
kitomeAirITS_Index_all <- which(sample_data(allITS_20hrs.ps)$sampleType== "kitome")
PCRcontrolAirITS_Index_all <- which(sample_data(allITS_20hrs.ps)$sampleType== "PCRcontrol") 
sample_data(allITS_20hrs.ps)[PCRcontrolAirITS_Index_all,] #double check, is just air, i.e., no phyllosphere-associated
# combine all of these 
airSampAndCtrlITS_Index_all <- c(airITS_Index_all, fieldControlITS_Index_all, washBufferControlITS_Index_all, kitomeAirITS_Index_all, PCRcontrolAirITS_Index_all)
ITS_airCtrlSampNames <- rownames(sample_data(allITS_20hrs.ps))[airSampAndCtrlITS_Index_all]
# Make a new phyloseq object
ITS_airCtrlSamps_20hrs.ps  <-   prune_samples(samples= ITS_airCtrlSampNames, x=allITS_20hrs.ps)
setdiff(ITS_airCtrlSampNames, rownames(sample_data(ITS_airCtrlSamps_20hrs.ps))) #looks right
setdiff(rownames(sample_data(ITS_airCtrlSamps_20hrs.ps)), ITS_airCtrlSampNames) #looks right
# Remove those taxa without reads
ITS_airCtrlSamps_20hrs.ps <- prune_taxa(taxa_sums(ITS_airCtrlSamps_20hrs.ps) > 0, ITS_airCtrlSamps_20hrs.ps)
ITS_airCtrlSamps_20hrs.ps # 5763 taxa and 147 samples
unique(sample_data(ITS_airCtrlSamps_20hrs.ps)$sampleType) #looks correct!
ITS_airCtrlSamps_20hrsASVs <- as.data.frame(as.matrix(otu_table(ITS_airCtrlSamps_20hrs.ps)))
# View(ITS_airCtrlSamps_20hrsASVs)

# TOTAL READ COUNTS IN SAMPLES VERSUS CONTROLS
sum(ITS_airCtrlSamps_20hrsASVs) #4,069,592 reads across ALL air samples and controls
sum(as.matrix(ITS_airCtrlSamps_20hrsASVs[sapply(ITS_airCtrlSamps_20hrsASVs, is.numeric)])) #matches above
# Get ASV tables for only samples and only blanks
sampsAirITS_20hrs <- ITS_airCtrlSamps_20hrsASVs[,which(sample_data(ITS_airCtrlSamps_20hrs.ps)$isControl != "control")]
CtrlAirITS_20hrs <- ITS_airCtrlSamps_20hrsASVs[,which(sample_data(ITS_airCtrlSamps_20hrs.ps)$isControl == "control")]
# Correct samples present
sort(colnames(sampsAirITS_20hrs)) == sort(rownames(sample_data(ITS_airCtrlSamps_20hrs.ps))[which(sample_data(ITS_airCtrlSamps_20hrs.ps)$sampleType == "air")])
sum(as.matrix(sampsAirITS_20hrs[sapply(sampsAirITS_20hrs, is.numeric)])) #4,007,616 just from air. So
sum(ITS_airCtrlSamps_20hrsASVs) - sum(as.matrix(sampsAirITS_20hrs[sapply(sampsAirITS_20hrs, is.numeric)])) #61,976 from controls
(sum(ITS_airCtrlSamps_20hrsASVs) - sum(as.matrix(sampsAirITS_20hrs[sapply(sampsAirITS_20hrs, is.numeric)])))/
  sum(ITS_airCtrlSamps_20hrsASVs)*100 # blanks = 1.522905 = 1.52% of all reads
61976/4069592*100 #gets same answer as above

# MEDIAN READ COUNT ACROSS CONTROLS (before contaminant removal)
median(colSums(CtrlAirITS_20hrs)) #582
sort(colSums(CtrlAirITS_20hrs)) #uh-oh,  air_ITS_114 and air_ITS_37 are high.  Both are, unsurprisingly, fieldControls :(
# Double check these
sample_data(allITS_20hrs.ps)$sampleType[rownames(sample_data(allITS_20hrs.ps)) %in% c("air_ITS_114", "air_ITS_37")] #yep
colSums(otu_table(allITS_20hrs.ps)[,colnames(otu_table(allITS_20hrs.ps)) %in% c("air_ITS_114", "air_ITS_37")])
# air_ITS_114  air_ITS_37 
# 10381       22163 
sample_data(allITS_20hrs.ps)[rownames(sample_data(allITS_20hrs.ps)) %in% c("air_ITS_114", "air_ITS_37"),] 

median(colSums(sampsAirITS_20hrs)) #35,188
# Controls, a deeper look
ITS_controlsTypeSum <- as.data.frame(matrix(ncol=2,nrow=ncol(CtrlAirITS_20hrs)))
colnames(ITS_controlsTypeSum) <- c("sampleName", "numberReads")
ITS_controlsTypeSum[,1] <- colnames(CtrlAirITS_20hrs)
ITS_controlsTypeSum[,2] <- as.numeric(colSums(CtrlAirITS_20hrs))
ITS_controlsTypeSum <- ITS_controlsTypeSum %>% 
  mutate(
    ctrlType = case_when(
      sampleName %in% washBufferControlNames ~ "washBuffer",
      sampleName %in% fieldControlNames ~ "fieldControl",
      sampleName %in% kitomeAirControlNames ~ "kitome",
      sampleName %in% PCRAirControlNames ~ "PCRAirControl",
      TRUE ~ NA_character_
    )) %>% 
      arrange(desc(numberReads))
ITS_controlsTypeSum <- ITS_controlsTypeSum[,c(1,3,2)] #re-order columns
# Get averages per type and add to dataframe 
ITS_controlsTypeSum <- ITS_controlsTypeSum %>% 
  group_by(ctrlType) %>% 
  mutate(meanType = mean(numberReads)) %>% 
  mutate(medianType =median(numberReads))
# View(ITS_controlsTypeSum)
# saved this September 16, 2025
# write.csv(ITS_controlsTypeSum, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_controlsTypeSum.csv")

# Below written and saved July 21, 2023
# write.csv(ITS_airTaxOTUsOrderedReNamed, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_airTaxOTUsOrderedDecontam.csv")

##################################################################
# DATA CLEANING PART IV
# (dealing with air contamination)
# (For bioaerosols, same approach as in fullITS_EDA_part1.R)
# Note: no ASVs were removed in this step, since no prevalent, probable 'contaminant'
# ASVs were found!
##################################################################
# REMOVE ASVS THAT OCCUR IN AT LEAST HALF OF CONTROLS (find which ASVs (i.e. rows) are present in at lest half of controls?)
# First, get a subset with just the controls
colnames(ITS_airTaxOTUsOrderedReNamed)
ITS_airTaxOTUsOrderedReNamed[nrow(ITS_airTaxOTUsOrderedReNamed),]
ITS_justControlsSubset <- ITS_airTaxOTUsOrdered[1:nrow(ITS_airTaxOTUsOrdered)-1,c(1, 112:148)] %>% # columns 112:148 are the controls. 
  #Remove last row since it is read count in each sample
  column_to_rownames(var="Row.names")
# View(ITS_justControlsSubset) #37 columns, as expected (one column for each control/blank)

sort(colSums(ITS_justControlsSubset)) #as explored earlier, there are a few blanks with many reads, but all but 3 have less than 3K
# check out samples with stupid high reads... As expected, it's a field control (code below shows this) so the filter probably came
# to me contaminated. Or, 37 actually had a cracked PEtri dish, but sample didn't seem to be exposed
sample_data(allITS_20hrs.ps)[which(rownames(sample_data(allITS_20hrs.ps)) %in% "air_ITS_37"),] #field control (22163 reads! D: )
sample_data(allITS_20hrs.ps)[which(rownames(sample_data(allITS_20hrs.ps)) %in% "air_ITS_114"),] #field control (10381 reads! D: )

# This for loop counts the number of times a given ASV shows up in the controls and places this number in 
ITS_removalDataFrame <- as.data.frame(matrix(data=NA,nrow= nrow(ITS_justControlsSubset),ncol=3))
colnames(ITS_removalDataFrame) <- c("ASV_name", "controlInstances", "removalDecision")
ITS_removalDataFrame$ASV_name <- rownames(ITS_justControlsSubset)

for (i in 1:nrow(ITS_justControlsSubset)){
  ITS_removalDataFrame$controlInstances[i] <- sum(ITS_justControlsSubset[i,] > 0)
  if (ITS_removalDataFrame$controlInstances[i] >= 19){
    ITS_removalDataFrame$removalDecision[i] <- "remove ASV"
  } else {
    ITS_removalDataFrame$removalDecision[i] <- "retain ASV"
  }
}

# Some checks to show that the above for loop worked 
ITS_justControlsSubset[1,]
ITS_removalDataFrame[1,]
ITS_justControlsSubset[2,]
ITS_removalDataFrame[2,]
ITS_justControlsSubset[10,]
ITS_removalDataFrame[10,]
ITS_justControlsSubset[1334,]
ITS_removalDataFrame[1334,]

# Pull out only the ASVs in the ITS_removalDataFrame to remove
ITS_ASVstoRemoveNames <- ITS_removalDataFrame$ASV_name[which(ITS_removalDataFrame$removalDecision == "remove ASV")]
length(ITS_ASVstoRemoveNames) #0
# THERE WERE NO TAXA AMONG FUNGI THAT APPEARED TO BE CONTAMINANTS BASED ON CRITERIA, SO NONE WILL BE REMOVED.

##################################################################
# DATA CLEANING PART V
# (dealing with phyllosphere-associated sample contamination)
# New to this script
##################################################################
#### 2. FOLIAR CONTROLS #####
# Pull out just foliar surface samples -- note that 2 PCR NTCs for phyllo failed (as shown above)
sort(unique(sample_data(allITS_20hrs.ps)$sampleType)) 
sample_data(allITS_20hrs.ps)[which(sample_data(allITS_20hrs.ps)$sampleType == "PCR_NegControl"),] #PCR controls did make it thru sequencing :(
# 1. GET DATA SUBSETS
## i. all phyllo and phyllo controls
phylloITS_20hrs.ps <- subset_samples(allITS_20hrs.ps, sampleType %in% c("phyllosphere", "phyllo_NegFieldControl", "phyllo_ExtBlank", "PCR_NegControl"))
unique(sample_data(phylloITS_20hrs.ps)$sampleType)
## ii. Just samples
samps_phylloITS_20hrs.ps <- subset_samples(phylloITS_20hrs.ps, isControl == "sample")
sort(colSums(otu_table(samps_phylloITS_20hrs.ps ))) #reads in each
median(colSums(otu_table(samps_phylloITS_20hrs.ps))) #42,182 added to manuscript (since no contaminants)

## iii. Just controls
Ctrl_phylloITS_20hrs.ps <- subset_samples(phylloITS_20hrs.ps, isControl == "control")
sort(colSums(otu_table(Ctrl_phylloITS_20hrs.ps))) #reads in each
# phyllo_ITS_PCR_NTC_1        phyllo_ITS_65 phyllo_ITS_PCR_NTC_2         phyllo_ITS_6        phyllo_ITS_66        phyllo_ITS_52 
#               2                    6                    6                  237                  290                  347 
# phyllo_ITS_28        phyllo_ITS_27        phyllo_ITS_51 
# 428                  503                  527 
median(colSums(otu_table(Ctrl_phylloITS_20hrs.ps)) ) #290, added to manuscript (since no contaminants)

# View(as.data.frame(as.matrix(otu_table(Ctrl_phylloITS_20hrs.ps))))
length(which(rowSums(otu_table(Ctrl_phylloITS_20hrs.ps))>0)) #37 ASVs were found in any abundance in these controls
length(which(rowSums(otu_table(Ctrl_phylloITS_20hrs.ps))>5)) #16 ASVs were found more than 5 times in these controls
length(which(rowSums(otu_table(Ctrl_phylloITS_20hrs.ps))>100)) #7 ASVs were found more than 100 times in these controls

# 1. IDENTIFY ANY ASVS THAT WERE DETECTED IN AT LEAST 4 CONTROLS (4/ORIGINAL 9)
ctrlsPhyllo_df <- as.data.frame(as.matrix(otu_table(Ctrl_phylloITS_20hrs.ps)))
head(ctrlsPhyllo_df)
# Subject phyllo to same procedure
for (i in 1:nrow(ctrlsPhyllo_df)) {
  ctrlsPhyllo_df$controlInstances[i] <- sum(ctrlsPhyllo_df[i, 1:9] > 0) #columns 1:7 because these are the samples
  if (ctrlsPhyllo_df$controlInstances[i] >= 4) {
    ctrlsPhyllo_df$removalDecision[i] <- "remove ASV"
  } else {
    ctrlsPhyllo_df$removalDecision[i] <- "retain ASV"
  }
}
# View(ctrlsPhyllo_df)
# Add number of times each ASV is found
colnames(ctrlsPhyllo_df)
ctrlsPhyllo_df$ASVsumCtrls <- rowSums(ctrlsPhyllo_df[,1:7])
# View(ctrlsPhyllo_df)
sort(rownames(ctrlsPhyllo_df)[which(ctrlsPhyllo_df$removalDecision == "remove ASV")]) 
# none!

##################################################################
# RAREFYING (EXPLORING, THEN RAREFYING AT 8.5K READS TO RETAIN (air) 112 SAMPLES 
##################################################################
# EXPLORE READS ACROSS SAMPLES
colnames(otu_table(allITS_20hrs.ps)) #show sample names (across all sample types)
allITS_seqspersample <- sort(colSums(otu_table(allITS_20hrs.ps)))
allITS_seqspersample
airSeqsIndex <- grep(x= names(allITS_seqspersample), pattern= "air")
length(airSeqsIndex) == length(ITS_airAllIndex) #147, as expected 
allITS_SeqNumberPlot <-barplot(allITS_seqspersample, main="ITS: Total Sequences Per Sample", ylab= "Number of Sequences", xlab= "Sample Number")
airITS_seqsPerSample <- allITS_seqspersample[airSeqsIndex]

# NOTE: the number in these indices are indices NOT number of seqs
airDrop5000 <- which(airITS_seqsPerSample < 5000) # drops 35 samples/blanks
length(airDrop5000) #35
airDrop6000 <- which(airITS_seqsPerSample < 6000) # drops 35 samples/blanks
length(airDrop6000) #35 samples
airDrop8000 <- which(airITS_seqsPerSample < 8000) # drops 35 samples/blanks
length(airDrop8000) #35 samples
airDrop8500 <- which(airITS_seqsPerSample < 8500) # drops 35 samples/blanks
length(airDrop8000) #35 samples
airDrop9000 <- which(airITS_seqsPerSample < 9000) # drops 37 samples/blanks
length(airDrop9000) #37 samples

ITS_airControlAllNames #made above near where I was making airTaxOTUsOrdered
# WHICH SAMPLES WOULD BE DROPPED IF USED THE FOLLOWING CRITERIA??
# 8.5K
length(which(names(airDrop8500) %in% ITS_airControlAllNames == TRUE)) #35 of these 35 dropped are controls!

# 9K
length(which(names(airDrop9000) %in% ITS_airControlAllNames == TRUE)) #35 of these 36 dropped are controls,
# so I'll rarefy at 8,500 reads so that I keep that sample/non-control that would be dropped if I rarefied at
# 9,000 reads

names(airDrop8500)
ITS_samplesToDrop <- names(airDrop8500) #names of the samples to drop-- i.e. rarefying at 8.5K reads
length(names(airITS_seqsPerSample)) - length(ITS_samplesToDrop) #this will leave us with 112 total samples/blanks

# Get sample data of the samples to be dropped - shows that no air samples are dropped, but 14 field controls, 
# 10 wash buffer controls, 10 kitomes, and 1 PCR control will be dropped following rarefying
ITS_droppedMetaData <- sample_data(allITS_20hrs.ps)[which(rownames(sample_data(allITS_20hrs.ps)) %in% ITS_samplesToDrop == TRUE),]
unique(ITS_droppedMetaData$sampleType)
length(which(ITS_droppedMetaData$sampleType == "air")) #0 are air samples
length(which(ITS_droppedMetaData$sampleType == "fieldControl")) #14 are field controls
length(which(ITS_droppedMetaData$sampleType == "washBufferControl")) #10 are wash buffer controls
length(which(ITS_droppedMetaData$sampleType == "kitome")) #10 are kitome
length(which(ITS_droppedMetaData$sampleType == "PCRcontrol")) #1 are PCRcontrol

# Samples that will be retained after rarefying!
samplesToKeepNames <- names(which(allITS_seqspersample >= 8500)) #get all of these ACROSS sample types
samplesToKeepNames 
length(grep(pattern= "air", x=samplesToKeepNames)) #shows 112, and since only 2 field blanks remained, then 110 are samples

set.seed(19)
allITS_rarefied.ps <- rarefy_even_depth(allITS_20hrs.ps, sample.size = 8500, replace = FALSE, trimOTUs = TRUE)
# 776 ASVs were removed because they are no longer present in any sample after random subsampling
# In addition, 64 samples (across all types, i.e., air, soil, phyllosphere, and controls).
# Now, there are 328 samples (across air, phyllo, and soil)

table((sample_data(allITS_rarefied.ps)$sampleType))
# air        fieldControl        phyllosphere                soil       soil_ExtBlank soilExtControlWater 
# 110                   2                  59                 155                   1                   1 

# June 17, 2024 -- saved on server
# save(allITS_rarefied.ps, file ="~/SRS_aeromicrobiome_2022/RobjectsSaved/allITS_rarefiedFinal_phyloseq") #saved on server June 17, 2024

# saved July 23, 2023 on own computer
# NOTE: Below is the most recent version on the computer. There is also an original version called "allITSrarefied_phyloseq"
# that is in same folder and is not the version to use.
# save(allITS_rarefied.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allITS_rarefied_phyloseq_July23") 

# saved once again on own computer to double-check September 16, 2025
# save(allITS_rarefied.ps, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/allITS_rarefied_phyloseq")
