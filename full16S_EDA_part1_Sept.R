# 16S full EDA rarefied - PART 1
# (re-started) JULY 12, 2023

# SCRIPT ORDER: bioinformatics_16S.R --> THIS SCRIPT --> full16S_EDArarefied_Part2.R

# This script investigates contamination more, and begins to deviate from the original full16S_EDA_part1July12.R at the section labeled DATA CLEANING PART III

# Description from original full16S_EDA_part1July12.R
# This script explores the bioinformatics output of the air 16S data (following dada2/bioinformatics on microbe server in script
# called bioinformatics_16S.R). Here, I first removed anything that: i) wasn't classified as Kingdom Bacteria, ii) wasn't at least classified to phylum level
# (i.e. phylum = NA), iii) removed chloroplasts and mitochondria. Second, I removed air samples not sampled at least 20 
# hours. Third, I removed ASVs determined to be AIR-associated contaminants (i.e., those that showed up in at least 19/39
# of the AIR controls/blanks). Fourth, I rarefied all (i.e., air, soil, and phyllosphere) at 5,500 reads.

# Main file saved in this script is the final, rarefied dataset: 


# This script's fungal counterpart is called fullITS_EDA_part1_Sept.R 

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

# 2. Function to get ASV table out of phyloseq so that we can view it better
# (Inspired by similar function at https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  return(as.data.frame(ASVTable))
}

#######################################################################################
# 
options(scipen=999) #remove pesky scientific notation in plots

# Set working directory
# setwd("~/SRS_aeromicrobiome_2022") (was for server)
list.files()
# Read in libraries
library("phyloseq"); packageVersion("phyloseq") #'1.52.0’
library("tidyverse"); packageVersion("tidyverse") #‘2.0.0’
library("vegan"); packageVersion("vegan") #‘2.7.1’
library("gridExtra"); packageVersion("gridExtra")  #‘2.3’  
library(Polychrome); packageVersion("Polychrome")  #‘1.5.4’

##################################################################################
# I. SET-UP, DATA CLEANING, RAREFACTION, AND FIRST TAXONOMIC & ORDINATION PLOTS
##################################################################################
all16S_seqtab_wTax_mctoolsr <- read.table("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/bioinformaticsOutput/all16S_seqtab_wTax_mctoolsr.txt", header=TRUE)
str(all16S_seqtab_wTax_mctoolsr)
head(all16S_seqtab_wTax_mctoolsr)
#View(all16S_seqtab_wTax_mctoolsr)
colnames(all16S_seqtab_wTax_mctoolsr)
rownames(all16S_seqtab_wTax_mctoolsr)

all16S_seqtab <- read.table("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/bioinformaticsOutput/all16S_seqtab_final.txt", header=T)
dim(all16S_seqtab) #40582 ASVs  418 samples 
# View(all16S_seqtab)

#### LOAD FILES IN FOR PHYLOSEQ #####
# 1. OTU table: ASVs/OTUs as rows, samples as columns. This is "seqtab"
all16S_otu_mat <- all16S_seqtab
all16S_seqtab$X #X is the ASV names
all16S_otu_mat <- all16S_otu_mat %>%
  tibble::column_to_rownames("X") #make it so that "X"(which is column with ASV names is the rownames columnn!)
# View(all16S_otu_mat)
head(all16S_otu_mat)
# "X" samples are soils, re-format them below
colnames(all16S_otu_mat) <- gsub(colnames(all16S_otu_mat), pattern = "X", replacement = "soil_16S_") 
colnames(all16S_otu_mat) #colnames are sample names
rownames(all16S_otu_mat) #ASV names
# View(all16S_otu_mat)

# 2. OTU TAXONOMY
# This is basically needs to be the same as the first and last columns of seqtab_wTax_mctoolsr
# So none of the guts that are essentially the ASV table
dim(all16S_seqtab_wTax_mctoolsr)
colnames(all16S_seqtab_wTax_mctoolsr) #samples (but first colname is "X.ASV_ID")
head(all16S_seqtab_wTax_mctoolsr)
head(all16S_seqtab_wTax_mctoolsr[,1]) #THE ASV names!
head(all16S_seqtab_wTax_mctoolsr[,ncol(all16S_seqtab_wTax_mctoolsr)])
# The last column  sonetimes has species info
step1 <- cbind(all16S_seqtab_wTax_mctoolsr[,1], all16S_seqtab_wTax_mctoolsr[,(ncol(all16S_seqtab_wTax_mctoolsr))]) #bind together ASV ID names and taxonomy  
head(step1)
colnames(step1) <- c("ASV_ID", "taxonomy") #make ASV_ID and taxonomy the column names
head(step1)

# Make new columns based on the ";" separator
require(tidyr)
all16S_tax_sep <- separate(as.data.frame(step1), col = taxonomy, into= c("Kingdom", "Phylum", "Class", 
                                                                         "Order", "Family", "Genus", "ASV_ID"),
                           sep = ";")
head(all16S_tax_sep)
colnames(all16S_tax_sep) #"Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "ASV_ID" 
# Continuing tidying up taxonomy table
# View(all16S_tax_sep)
rownames(all16S_tax_sep)
all16S_tax_sep <- all16S_tax_sep %>% #make ASV_ID column the rownames to make it compatible with phyloseq
  tibble::column_to_rownames("ASV_ID")
head(all16S_tax_sep)

# 3. SAMPLE METADATA--
allMetadata <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/aeroAllMetadat_Apr20_2023.csv", row.names=1) #made in ITS_metadataAndMappingFileSetUp.R
# View(allMetadata)
colnames(allMetadata)
rownames(allMetadata)
airSamplesALL_df <- allMetadata %>%
  tibble::column_to_rownames("sampleNumber")
# View(airSamplesALL_df)
# Replace "ITS" with "16S"-- sample names and metadata are otherwise the same (processed was ITS first)
rownames(airSamplesALL_df) <- gsub(x=rownames(airSamplesALL_df), pattern="_ITS_", replacement = "_16S_")

# WHAT ARE THE ORGINAL METADATA NUMBERS (BEFORE ACCOUNTING FOR SAMPLES LOST, BELOW)?
# View(allMetadata)
colnames(allMetadata)
# Wee lil for loop to print out type of sample and number
for (i in 1:length(unique(allMetadata$sampleType))){
  print(cbind(unique(allMetadata$sampleType)[i],length(which(allMetadata$sampleType== unique(allMetadata$sampleType)[i]))))
}
# Total number of samples not counting soil controls:
126 + 10 + 2 + 11 + 16 + 59 + 4 +2 + 3 #233

# WHICH SAMPLES ARE IN THE OTU TABLE BUT NOT IN THE METADATA?
setdiff(sort(colnames(all16S_otu_mat)), sort(rownames(airSamplesALL_df))) # differences below
# "air_16S_PCR_NTC_1_H12" "ExtControlWater_11"    "ExtControlWater_16"    "ExtControlWater_18"    "ExtControlWater_3"     "ExtControlWater_5"    
# "ExtControlWater_6"     "ExtControlWater_7"     "ExtControlWater_8"     "PCR_NTC"

# These are in all16S_otu_mat but not in airSamplesALL_df
colnames(all16S_otu_mat)[which(colnames(all16S_otu_mat) %in% rownames(airSamplesALL_df) == FALSE)]
# [1] "air_16S_PCR_NTC_1_H12" "ExtControlWater_11"    "ExtControlWater_16"    "ExtControlWater_18"    "ExtControlWater_3"    
# [6] "ExtControlWater_5"     "ExtControlWater_6"     "ExtControlWater_7"     "ExtControlWater_8"     "PCR_NTC" 

# These are in airSamplesALL_df but not in all16S_otu_mat
rownames(airSamplesALL_df)[which(rownames(airSamplesALL_df) %in% colnames(all16S_otu_mat) == FALSE)]

# INVESTIGATE THESE CONTROLS (which are currently not in the metadata)
air_16S_PCR_NTC_2_H12index <- which(colnames(all16S_otu_mat) == "air_16S_PCR_NTC_1_H12")
sum(all16S_otu_mat[,air_16S_PCR_NTC_2_H12index]) # 2,249, so need to keep and investigate further 
ExtControlWater_11index <- which(colnames(all16S_otu_mat) == "ExtControlWater_11")
sum(all16S_otu_mat[,ExtControlWater_11index]) #276 reads
ExtControlWater_16index <- which(colnames(all16S_otu_mat) == "ExtControlWater_16")
sum(all16S_otu_mat[,ExtControlWater_16index]) #972 reads
ExtControlWater_18index <- which(colnames(all16S_otu_mat) == "ExtControlWater_18")
sum(all16S_otu_mat[,ExtControlWater_18index]) #431 reads 
ExtControlWater_3index <- which(colnames(all16S_otu_mat) == "ExtControlWater_3")
sum(all16S_otu_mat[,ExtControlWater_3index]) #519 reads 
ExtControlWater_5index <- which(colnames(all16S_otu_mat) == "ExtControlWater_5")
sum(all16S_otu_mat[,ExtControlWater_5index]) #305 reads 
ExtControlWater_6index <- which(colnames(all16S_otu_mat) == "ExtControlWater_6")
sum(all16S_otu_mat[,ExtControlWater_6index]) #165 reads 
ExtControlWater_7index <- which(colnames(all16S_otu_mat) == "ExtControlWater_7")
sum(all16S_otu_mat[,ExtControlWater_7index]) #332 reads 
ExtControlWater_8index <- which(colnames(all16S_otu_mat) == "ExtControlWater_8")
sum(all16S_otu_mat[,ExtControlWater_8index]) #77 reads 
PCR_NTCindex <- which(colnames(all16S_otu_mat) == "PCR_NTC")
sum(all16S_otu_mat[,PCR_NTCindex]) #335 reads 

# ALL of these ExtControlWaters and PCR_NTC are soil and hve a low number of reads. Further, since they were dropped from my
# earlier soil edge effects project, they are not on NCBI. So, I will drop them here.
colsToDropOTUmat <- which(colnames(all16S_otu_mat) %in% rownames(airSamplesALL_df) == FALSE)[2:10] #2:10 to not drop air_16S_PCR_NTC_2_H12index
all16S_otu_matTrimmed <- all16S_otu_mat[,-colsToDropOTUmat]
# Air grepping
length(grep(x= colnames(all16S_otu_matTrimmed), pattern= "air")) #165 air samples (and controls) still at this point

# CHECK OUT AND FURTHER EVALUATE DIFFERENCES
# Now PCR control in OTU table but not in metadata
setdiff(sort(colnames(all16S_otu_matTrimmed)), sort(rownames(airSamplesALL_df))) # "air_16S_PCR_NTC_1_H12"
setdiff(sort(rownames(airSamplesALL_df)), sort(colnames(all16S_otu_matTrimmed)))
# [1] "air_16S_PCR_NTC_1"    "phyllo_16S_PCR_NTC_1" "phyllo_16S_PCR_NTC_2"
# All of these samples did not make it through bioinformatics pipeline, so can be dropped from metadata. However,
# since "air_16S_PCR_NTC_1_H12" replaces "air_16S_PCR_NTC_1", I'll keep it in the metadata but just re-name it
# "air_16S_PCR_NTC_1_H12"
# Get which samples are in airSamplesALL_df that are not in all16S_otu_matTrimmed
dfToDrop <- which(rownames(airSamplesALL_df) %in% colnames(all16S_otu_matTrimmed) == FALSE)
rownames(airSamplesALL_df)[dfToDrop] #"air_16S_PCR_NTC_1_H12" "phyllo_16S_PCR_NTC_1"  "phyllo_16S_PCR_NTC_2" 
rownames(airSamplesALL_df)[dfToDrop[1]] <- "air_16S_PCR_NTC_1_H12" #rename sample
airSamplesALL_df <- airSamplesALL_df[-dfToDrop[2:3],] 
rownames(airSamplesALL_df)
# The next two lines confirm that these are the same
setdiff(sort(colnames(all16S_otu_matTrimmed)), sort(rownames(airSamplesALL_df)))
setdiff(sort(rownames(airSamplesALL_df)), sort(colnames(all16S_otu_matTrimmed)))

# Transform otu table and taxonomy tables into matrices
all16S_otu_mat <- as.matrix(all16S_otu_matTrimmed)
all16S_tax_mat <- as.matrix(all16S_tax_sep)
# View(all16S_tax_mat)

# Transform to phyloseq objects
all16S_OTU <- otu_table(all16S_otu_matTrimmed, taxa_are_rows = TRUE)
all16S_TAX <- tax_table(all16S_tax_mat)
all16S_samples <- sample_data(airSamplesALL_df)
all16S_raw.ps <- phyloseq(all16S_OTU, all16S_TAX, all16S_samples)
all16S_raw.ps #40582 taxa and 408 samples
colnames(otu_table(all16S_raw.ps))
sample_variables(all16S_raw.ps)

# SAVED IN ORIGINAL full16S_EDA_part1July12.R, AS DOCUMENTED BELOW
# saved July 3, 2023 (own computer) and June 10, 2024 (on server)
# save(all16S_raw.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16S_raw_phyloseq") #computer
# save(all16S_raw.ps, file ="~/SRS_aeromicrobiome_2022/RobjectsSaved/all16S_raw_phyloseq") #saved on server June 10, 2024

##################################################################
# DATA CLEANING PART I
# (remove anything not bacteria those not classfied
# to at least phylum level, and chloroplasts and mitochondria)
##################################################################
# SOME SET UP
# Get taxonomy table out of phyloseq:
all16S_taxTable <- taxtable_outta_ps(all16S_raw.ps)
all16S_ASVs <-  ASVs_outta_ps(all16S_raw.ps)
all16S_ASVs$ASVabundance <- rowSums(all16S_ASVs)
# View(all16S_taxTable) #raw has 40,582 ASVs
# View(all16S_ASVs)

# REMOVE EVERYTHING THAT ISN'T BACTERIA AND THAT ISN'T CLASSIFIED TO AT LEAST PHYLUM:
all16S_taxTrimmed1 <- all16S_taxTable %>% 
  filter(Kingdom == "Bacteria") %>% #Remove ASVs where kingdom is not Bacteria
  filter(Phylum != "NA") #Remove ASVs where phylum is unknown ("NA" in first column) 
unique(all16S_taxTrimmed1$Kingdom) #all bacteria
# View(all16S_taxTrimmed1) 
length(which(is.na(all16S_taxTable$Phylum))) # 0 because removed above
length(which(all16S_taxTable$Kingdom == "Archaea")) #only removed  182 Archaea

# INVESTIGATE ARCHAEA (before removing any reads or samples)
# What was the average # of archaea per air sample?
head(all16S_taxTable)
archaeaASVnames <- rownames(all16S_taxTable)[which(all16S_taxTable$Kingdom == "Archaea")] #get Archaea names
airNamesTrimmingStep <- rownames(as.matrix(sample_data(all16S_raw.ps)))[as.data.frame(as.matrix(sample_data(all16S_raw.ps)))$sampleType == "air"]
head(all16S_ASVs)
# Get numbers of archaea in each air sample by looking at only archaea rows and only air columns
sort(colSums(all16S_ASVs[rownames(all16S_ASVs) %in% archaeaASVnames,colnames(all16S_ASVs) %in% airNamesTrimmingStep]))
# this shows that these samples below have a lot (more than 200 archaea). Which ASVs?
# Below shows that it doesn't tend to be just one ASV across all of these
# View(all16S_ASVs[rownames(all16S_ASVs) %in% archaeaASVnames,colnames(all16S_ASVs) %in% c("air_16S_153", "air_16S_2", "air_16S_60", "air_16S_39", "air_16S_27", "air_16S_152", "air_16S_137")])

# What percentage were archaea of 16S rRNA reads in bioaerosol samples?
colSums(all16S_ASVs[rownames(all16S_ASVs) %in% archaeaASVnames,colnames(all16S_ASVs) %in% airNamesTrimmingStep])
# Below gets ASV table for just archaea ASVs in air. Get sum of all reads across this smaller table, 
# divide it by sum of all reads across full 16S dataset, and then multiply by 100 for "% of 16S rRNA gene reads in 
# raw bioaerosol samples" reported in the manuscript:
# Sum of all Archaea reads across air samples
archAirReadsSum <- sum(all16S_ASVs[rownames(all16S_ASVs) %in% archaeaASVnames,colnames(all16S_ASVs) %in% airNamesTrimmingStep])
# shows that this summing worked as expected
archAirReadsSum == sum(unname(colSums(all16S_ASVs[rownames(all16S_ASVs) %in% archaeaASVnames,colnames(all16S_ASVs) %in% airNamesTrimmingStep])))
# Sum of all 16S rRNA reads across air samples
allAirReadsSum <- sum(all16S_ASVs[,colnames(all16S_ASVs) %in% airNamesTrimmingStep]) #keep all rows (all ASVs) and only air samp columns
unique(rownames(all16S_ASVs[,colnames(all16S_ASVs) %in% airNamesTrimmingStep]) == rownames(all16S_ASVs)) #check all ASVs included
# Divide # archaea reads by # of all reads and get %
archAirReadsSum/allAirReadsSum *100 #shows less than 0.1%!

# REMOVE CHLOROPLASTS AND MITOCHONDRIA
all16S_taxTrimmed2 <- all16S_taxTrimmed1 %>% 
  filter(Family != "Mitochondria") %>% #remove mitochondria
  filter(Order != "Chloroplast") #remove chloroplasts

# Filtered out a total of 6,808 ASVS, for a resulting 33,774 ASVs.
dim(all16S_taxTrimmed2)[1]
dim(all16S_taxTable)[1] - dim(all16S_taxTrimmed2)[1]

# TRIM PHYLOSEQ OBJECT (AND MAKE A NEW PHYLOSEQ OBJECT)
# We now have to trim the OTU table, because it still has the taxa in it that were filtered out
# above. 

all16S_ASVsToKeep1 <- rownames(all16S_taxTrimmed2) #grab names of ASVs to keep (that were filtered out of tax
# table above)
length(all16S_ASVsToKeep1) == length(all16S_ASVsToKeep1) #check shows this is true 

# Make a new phyloseq object with only these taxa. Note that contamination air taxa haven't been dealt with yet!
all16S_trimmed1.ps <- prune_taxa(all16S_ASVsToKeep1, all16S_raw.ps) #new phyloseq object with just the stuff we want!
all16S_trimmed1.ps # Now 33774 taxa and 408 samples
# Check to make sure new phyloseq object looks as expected
unique(rownames(all16S_taxTrimmed2) == rownames(otu_table(all16S_trimmed1.ps))) #yep! 
length(which(sample_data(all16S_trimmed1.ps)$sampleType=="air")) #126 air (samples, not blanks) still

##################################################################
# DATA CLEANING PART II
# (Add more metadata (controls versus samples) and remove samples
# not sampled at least 20 hours)
##################################################################
# Get sample data
sampDat <- as.data.frame(as.matrix(sample_data(all16S_trimmed1.ps)))
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

cbind(sampDat$sampleType, sampDat$isControl) #quick little spot check shows that controls look correctly defined

# Make all16S_trimmed1.ps have this new metadata -- all16S_noNAs2
all16S_trimmed1.ps <- phyloseq(tax_table(all16S_trimmed1.ps), otu_table(all16S_trimmed1.ps), sample_data(sampDat))
all16S_trimmed1.ps
# 33774 taxa and 408 samples

#### 1. Drop air samples that ran less than 20 hours (air)
colnames(sample_data(all16S_trimmed1.ps))
airSamps <- which(sample_data(all16S_trimmed1.ps)[,21] == "air") #column 21 is sampleType
length(airSamps) #126
# #View(as.data.frame(as.matrix(sample_data(all16S_trimmed1.ps)[airSamps,])))
colnames(sample_data(all16S_trimmed1.ps))[15] #"SampledRunTime"
greater20h <- which(sample_data(all16S_trimmed1.ps)[,15] >= 20.00)
length(greater20h)
##View(as.data.frame(as.matrix(sample_data(all16S_trimmed1.ps)[intersect(airSamps, greater20h),])))
dim(as.data.frame(as.matrix(sample_data(all16S_trimmed1.ps)[intersect(airSamps, greater20h),]))) 
as.numeric(as.data.frame(as.matrix(sample_data(all16S_trimmed1.ps)[intersect(airSamps, greater20h),]))$SampledRunTime)
# This shows that 110 out of the 126 air samples (not counting controls) sampled for at least 20 hours.
sort(as.numeric(as.data.frame(as.matrix(sample_data(all16S_trimmed1.ps)[intersect(airSamps, greater20h),]))$SampledRunTime))
# In fact, the lowest in this group (air_16S_19) sampled for 23.410 hours!

# Find which air samples to drop by finding those air samples that did not sample at least 20 hours
badAir <- setdiff(airSamps, greater20h)
length(badAir) #16 bad air samples (as expected!)
badAirNames <- sort(rownames(sample_data(all16S_trimmed1.ps)[badAir,])) #gives names of "bad samples"
badAirNamesNumbersOnly <- gsub(x=rownames(sample_data(all16S_trimmed1.ps)[badAir,]), pattern = "air_16S_", replacement = "")
sort(badAirNamesNumbersOnly)
# Double check that these do, in fact, have less than 20 hours 
sample_data(all16S_trimmed1.ps)$SampledRunTime[badAir] #yep, this check worked
sample_data(all16S_trimmed1.ps)$SampledRunTime[rownames(sample_data(all16S_trimmed1.ps)) %in% badAirNames] #this check also worked!
# Finally visualize one more time:
sample_data(all16S_trimmed1.ps)[rownames(sample_data(all16S_trimmed1.ps)) %in% badAirNames,] #this check also worked!

# SAVED IN ORIGINAL full16S_EDA_part1July12.R, AS DOCUMENTED BELOW
# saved Aug 4, 2023 in order to drop these in chloroplast exploration (see chloroplast_16SASVs_EDA.R)
# saveRDS(badAirNames, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/RObjectsSaved/badAirNames")

# View(as.data.frame(as.matrix(sample_data(all16S_trimmed1.ps)[badAir,])))  #double check-- looks good!
rownames(sample_data(all16S_trimmed1.ps)[-badAir,]) 
length(rownames(sample_data(all16S_trimmed1.ps)[-badAir,])) == (408-16)
# View(as.data.frame(as.matrix(sample_data(all16S_trimmed1.ps)[-badAir,])))  #double check-- looks good!
goodNames <- rownames(sample_data(all16S_trimmed1.ps)[-badAir,])
length(goodNames) #as expected, this is 392 samples
all16S_20hrs.ps <- prune_samples(goodNames, all16S_trimmed1.ps)
# View(as.data.frame(as.matrix(sample_data(all16S_20hrs.ps)))) #looks good!
length(which(sample_data(all16S_20hrs.ps)$sampleType == "air")) #110 samples still

# Saved Dec. 17, 2025 (so that I could scrape the metadata to make NCBI thingy)
# saveRDS(all16S_20hrs.ps, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allI6S_20hrs_phyloseq.rds")

# What are the names of all of the air samples (including ALL controls) that were sampled at least 20 hours?
airAllIndex <- grep(pattern= "air", x= rownames(as.data.frame(as.matrix(sample_data(all16S_20hrs.ps))))) #get all of the rownames that contain "air"
onlyGoodAirNames <- rownames(as.data.frame(as.matrix(sample_data(all16S_20hrs.ps))))[airAllIndex] #yes, this pulled out everything with air in the title
onlyGoodAirNames
length(airAllIndex) #149 air samples and air-related controls

# Get sample names that have "air" in them 
allAirSampsIndex <- grep(x= colnames(otu_table(all16S_20hrs.ps)), pattern= "air")
airOTUs20h <- as.data.frame(as.matrix(otu_table(all16S_20hrs.ps)[,allAirSampsIndex]))
airOTUs20h <- airOTUs20h[which(rowSums(airOTUs20h) > 0),] #keep only those rows for which there is at least one count of the air ASV
# View(airOTUs20h)
# Get taxa table 
taxall16S_20hrs <- as.data.frame(as.matrix(tax_table(all16S_20hrs.ps)))
# merge ASV taxonomy onto the end of the trimmed air only OTU table
airTaxOTUs <- merge(airOTUs20h, taxall16S_20hrs, all=FALSE, by=0)
# SAVED IN ORIGINAL full16S_EDA_part1July12.R, AS DOCUMENTED BELOW
# written and saved July 12, 2023 on own computer to examine in Excel
# write.csv(airTaxOTUs, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airTaxOTUsDecontam_July12.csv")

colnames(airTaxOTUs)
# Add in read count for each sample and make this the last row
airTaxOTUsSums <- rbind(airTaxOTUs,c(NA,colSums(airTaxOTUs[,2:150]), rep(NA, 6))) #2:150 is just the columns representing samples. rep(NA, 6) just so taxonomy has empty row below it and R doesnt get upset
# View(airTaxOTUsSums)

##################################################################
# DATA CLEANING PART III-- more investigation of controls and 
# archaea that were filtered out (BEFORE contaminant removal)
##################################################################
#### 1. BIOAEROSOLS: A DEEPER LOOK AT CONTROLS VERSUS SAMPLES BEFORE CONTAMINANT REMOVAL #####
# Get names and indices of different types of air samples
washBufferControlNames <- rownames(sample_data(all16S_20hrs.ps))[which(sample_data(all16S_20hrs.ps)$sampleType == "washBufferControl")]
washBufferControlIndex <- which(colnames(airTaxOTUsSums) %in% washBufferControlNames == TRUE)
fieldControlNames <- rownames(sample_data(all16S_20hrs.ps))[which(sample_data(all16S_20hrs.ps)$sampleType == "fieldControl")]
fieldControlIndex <- which(colnames(airTaxOTUsSums) %in% fieldControlNames == TRUE)
kitomeAirControlNames <- rownames(sample_data(all16S_20hrs.ps))[which(sample_data(all16S_20hrs.ps)$sampleType == "kitome")]
kitomeIndex <- which(colnames(airTaxOTUsSums) %in% kitomeAirControlNames == TRUE)
PCRAirControlNames <- rownames(sample_data(all16S_20hrs.ps))[which(sample_data(all16S_20hrs.ps)$sampleType == "PCRcontrol")]
PCR_Index <- which(colnames(airTaxOTUsSums) %in% PCRAirControlNames == TRUE)
airNames <- rownames(sample_data(all16S_20hrs.ps))[which(sample_data(all16S_20hrs.ps)$sampleType == "air")]
airSampsIndex <- which(colnames(airTaxOTUsSums) %in% airNames == TRUE)
# Make a big index based on all of the control samples 
controlsIndex <- c(fieldControlIndex, washBufferControlIndex, kitomeIndex, PCR_Index)
length(controlsIndex)
# get names of all af these to use later
airControlAllNames <- c(washBufferControlNames, fieldControlNames, kitomeAirControlNames, PCRAirControlNames)

# Just air + associated controls after this:
airIndex_all <- which(sample_data(all16S_20hrs.ps)$sampleType== "air")
fieldControlIndex_all <- which(sample_data(all16S_20hrs.ps)$sampleType== "fieldControl")
washBufferControlIndex_all <- which(sample_data(all16S_20hrs.ps)$sampleType== "washBufferControl")
kitomeAirIndex_all <- which(sample_data(all16S_20hrs.ps)$sampleType== "kitome")
PCRcontrolAirIndex_all <- which(sample_data(all16S_20hrs.ps)$sampleType== "PCRcontrol") 
sample_data(all16S_20hrs.ps)[PCRcontrolAirIndex_all,] #double check, is just air
# combine all of these 
airSampAndCtrlIndex_all <- c(airIndex_all, fieldControlIndex_all, washBufferControlIndex_all, kitomeAirIndex_all, PCRcontrolAirIndex_all)
airCtrlSampNames <- rownames(sample_data(all16S_20hrs.ps))[airSampAndCtrlIndex_all]
# Make a new phyloseq object
I6S_airCtrlSamps_20hrs.ps  <-   prune_samples(samples= airCtrlSampNames, x=all16S_20hrs.ps)
setdiff(airCtrlSampNames, rownames(sample_data(I6S_airCtrlSamps_20hrs.ps))) #looks right
setdiff(rownames(sample_data(I6S_airCtrlSamps_20hrs.ps)), airCtrlSampNames) #looks right
# Remove those taxa without reads
I6S_airCtrlSamps_20hrs.ps <- prune_taxa(taxa_sums(I6S_airCtrlSamps_20hrs.ps) > 0, I6S_airCtrlSamps_20hrs.ps)
I6S_airCtrlSamps_20hrs.ps # 5278 taxa and 149 samples
I6S_airCtrlSamps_20hrsASVs <- as.data.frame(as.matrix(otu_table(I6S_airCtrlSamps_20hrs.ps)))
# View(I6S_airCtrlSamps_20hrsASVs)

# TOTAL READ COUNTS IN SAMPLES VERSUS CONTROLS
sum(I6S_airCtrlSamps_20hrsASVs) #1,898,574 reads across ALL air samples and controls
sum(as.matrix(I6S_airCtrlSamps_20hrsASVs[sapply(I6S_airCtrlSamps_20hrsASVs, is.numeric)]))
sampsAir16S_20hrs <- I6S_airCtrlSamps_20hrsASVs[,which(sample_data(I6S_airCtrlSamps_20hrs.ps)$isControl != "control")]
CtrlAir16S_20hrs <- I6S_airCtrlSamps_20hrsASVs[,which(sample_data(I6S_airCtrlSamps_20hrs.ps)$isControl == "control")]
# Correct samples present
sort(colnames(sampsAir16S_20hrs)) == sort(rownames(sample_data(I6S_airCtrlSamps_20hrs.ps))[which(sample_data(I6S_airCtrlSamps_20hrs.ps)$sampleType == "air")])
sum(as.matrix(sampsAir16S_20hrs[sapply(sampsAir16S_20hrs, is.numeric)])) #1,732,090 just from air. So
sum(I6S_airCtrlSamps_20hrsASVs) - sum(as.matrix(sampsAir16S_20hrs[sapply(sampsAir16S_20hrs, is.numeric)])) #166,484 from controls
166484/1898574*100 # blanks = 8.77% of all reads
39/ncol(I6S_airCtrlSamps_20hrsASVs)*100 # 26.1745%, or more than a quarter of samples were controls

# MEDIAN READ COUNT ACROSS CONTROLS (before contaminant removal)
median(colSums(CtrlAir16S_20hrs)) #3036
sort(colSums(CtrlAir16S_20hrs))
median(colSums(sampsAir16S_20hrs)) #11317
# Controls, a deeper look
controlsTypeSum <- as.data.frame(matrix(ncol=2,nrow=39))
colnames(controlsTypeSum) <- c("sampleName", "numberReads")
controlsTypeSum[,1] <- colnames(CtrlAir16S_20hrs)
controlsTypeSum[,2] <- colSums(CtrlAir16S_20hrs)
controlsTypeSum <- controlsTypeSum %>% 
  mutate(
    ctrlType = case_when(
      sampleName %in% washBufferControlNames ~ "washBuffer",
      sampleName %in% fieldControlNames ~ "fieldControl",
      sampleName %in% kitomeAirControlNames ~ "kitome",
      sampleName %in% PCRAirControlNames ~ "PCRAirControl",
      TRUE ~ NA_character_
    ) 
  )
controlsTypeSum <- controlsTypeSum[,c(1,3,2)] #re-order columns
# Get means per type and add to dataframe 
controlsTypeSum <- controlsTypeSum %>% 
  group_by(ctrlType) %>% 
  mutate(meanType = mean(numberReads)) %>% 
  mutate(medianType =median(numberReads))
# View(controlsTypeSum)
# saved this September 10, 2025
# write.csv(controlsTypeSum, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/controlsTypeSum.csv")

#### 2. PERCENTAGES OF ARCHAEA ASVS ACROSS ALL RETAINED AIR SAMPLES (I.E., ONLY THOSE SAMPLED TO 20 HOURS) #####
all16S_20hrsSampDat <- as.data.frame(as.matrix(sample_data(all16S_20hrs.ps))) #get retained samples
goodAirOnlySampNames <- rownames(all16S_20hrsSampDat)[all16S_20hrsSampDat$sampleType == "air"]
all16S_20hrsSampDat[rownames(all16S_20hrsSampDat) %in% goodAirOnlySampNames,21] #air, as expected
# Calculate using the dataframe made before removing archaea (i.e. upstream of removing chloroplasts, mitos, less than 20 hour samples)
mean(colSums(all16S_ASVs[rownames(all16S_ASVs) %in% archaeaASVnames,colnames(all16S_ASVs) %in% goodAirOnlySampNames])/
       colSums(all16S_ASVs[,colnames(all16S_ASVs) %in% goodAirOnlySampNames]) * 100)
# 0.06900057% of reads in raw air samples were archaea
names(which(rowSums(all16S_ASVs[rownames(all16S_ASVs) %in% archaeaASVnames,colnames(all16S_ASVs) %in% goodAirOnlySampNames])>0))
# these are the only archaeal ASVs that were present in raw air samples: "ASV_391"   "ASV_809"   "ASV_1479"  "ASV_1751"  "ASV_2853"  "ASV_3433"  "ASV_3469"  "ASV_7299"  "ASV_10799" "ASV_12186" "ASV_13838"
# "ASV_14605" "ASV_16131" "ASV_20562" "ASV_26511" "ASV_28084"
which(rowSums(all16S_ASVs[rownames(all16S_ASVs) %in% archaeaASVnames,colnames(all16S_ASVs) %in% goodAirOnlySampNames])>0)

##################################################################
# DATA CLEANING PART IV
# (dealing with air contamination)
# (For bioaerosols, same approach as in full16S_EDA_part1July12.R)
##################################################################
#### 1. BIOAEROSOL CONTAMINANTS #####
# Re-order for easier comparison
colnames(airTaxOTUsSums) #cols 151:156 is taxonomy
airTaxOTUsOrdered <- airTaxOTUsSums[, c(1, airSampsIndex, controlsIndex, 151:156)]
#View(airTaxOTUsOrdered)
airTaxOTUsOrderedReNamed <- airTaxOTUsOrdered #make another one to rename based on control type

# Rename for easier comparison
colnames(airTaxOTUsOrderedReNamed)[which(colnames(airTaxOTUsOrderedReNamed) %in% washBufferControlNames == TRUE)] <- "bufferCtrl"
colnames(airTaxOTUsOrderedReNamed)[which(colnames(airTaxOTUsOrderedReNamed) %in% fieldControlNames == TRUE)] <- "fieldCtrl"
colnames(airTaxOTUsOrderedReNamed)[which(colnames(airTaxOTUsOrderedReNamed) %in% kitomeAirControlNames == TRUE)] <- "kitome"
colnames(airTaxOTUsOrderedReNamed)[which(colnames(airTaxOTUsOrderedReNamed) %in% PCRAirControlNames == TRUE)] <- "PCRcontrol"
# View(airTaxOTUsOrderedReNamed)
airTaxOTUsOrderedReNamed_clean <- airTaxOTUsOrderedReNamed
colnames(airTaxOTUsOrderedReNamed_clean)[1] <- "ASV_name"
colnames(airTaxOTUsOrderedReNamed_clean[,2:150])
airTaxOTUsOrderedReNamed_clean$ASVsum <- rowSums(airTaxOTUsOrderedReNamed_clean[,2:150])
dim(airTaxOTUsOrderedReNamed_clean)
# Note that last row with no is the total number in each sample, i.e. row 5279
tail(rownames(airTaxOTUsOrderedReNamed_clean))
# Sept 10, 2025
# write.csv(airTaxOTUsOrderedReNamed_clean, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airTaxOTUsOrderedReNamed_clean.csv")

### LOOK AT JUST ASV_10, Staphylococcus,ACROSS CONTROLS ###
colnames(airTaxOTUsOrderedReNamed)
# Confirm that it is Staphylococcus
colnames(airTaxOTUsOrderedReNamed)[111]; colnames(airTaxOTUsOrderedReNamed)[112]; colnames(airTaxOTUsOrderedReNamed)[150]; colnames(airTaxOTUsOrderedReNamed)[151];  #112 is first control, 150 is last control
airTaxOTUsOrderedReNamed[which(airTaxOTUsOrderedReNamed[,1] %in% "ASV_10"),112:ncol(airTaxOTUsOrderedReNamed)]
# Sort the prevalence in controls
staphControls <- airTaxOTUsOrderedReNamed[which(airTaxOTUsOrderedReNamed[,1] %in% "ASV_10"),112:150] # ASV_10 is Staph
staphControls
rownames(staphControls) <- NULL
length(sort(as.matrix(staphControls))) #39 controls, as expected
length(which(as.matrix(staphControls) >0)) #13,so present in 13 samples
staphControls[which(as.matrix(staphControls) >0)] #pull out the controls that contained the ASV. Shows found min of 1 time, max of 11
length(staphControls[which(as.matrix(staphControls) >0)]) #present in 13, matches above
39-13 #not in 26
# How many reads did these staph-containing controls have overall (before removal of contaminants)?
airOTUsJustControls <- airOTUs20h[, colnames(airOTUs20h) %in% airControlAllNames]
sort(colnames(airOTUsJustControls)) == sort(airControlAllNames) #correct subsetting, shows that there are only controls
staphIndexSamps <- which(airOTUsJustControls[which(rownames(airOTUsJustControls) %in% "ASV_10"==TRUE),]>0)
# Pull out the ASV table for only those samples with ASV_10
staphyControlsOnlyASVtab <- airOTUsJustControls[,staphIndexSamps]
# View(staphyControlsOnlyASVtab)
staphyControlsOnlyASVtab$ASVabund <- rowSums(staphyControlsOnlyASVtab)
# Only retain ASVs with any reads
length(which(staphyControlsOnlyASVtab$ASVabund > 0)) #527 ASVs in these control samples with any reads
staphyControlsOnlyASVtab <- staphyControlsOnlyASVtab[which(staphyControlsOnlyASVtab$ASVabund > 0), ] 
# View(staphyControlsOnlyASVtab)
staphyControlsOnly_colSums <- colSums(staphyControlsOnlyASVtab[,1:(ncol(staphyControlsOnlyASVtab)-1)])
staphyControlsOnly_colSums_df <- as.data.frame(matrix(ncol=2, nrow= length(staphyControlsOnly_colSums)))
colnames(staphyControlsOnly_colSums_df) <- c("sampName", "numberReads")
staphyControlsOnly_colSums_df[,1] <- names(staphyControlsOnly_colSums)
staphyControlsOnly_colSums_df[,2] <- unname(staphyControlsOnly_colSums)
# Add in sampleType
staphyControlsOnly_colSums_df <- staphyControlsOnly_colSums_df %>%
  mutate(
    sampType = case_when(
      sampName %in% washBufferControlNames ~ "washBuffer",
      sampName %in% fieldControlNames ~ "fieldControl",
      sampName %in% kitomeAirControlNames ~ "kitome",
      sampName %in% PCRAirControlNames ~ "PCRAirControl",
      TRUE ~ NA_character_
    ) 
  ) %>% 
  arrange(numberReads)
staphyControlsOnly_colSums_df

### REMOVE ASVS THAT OCCUR IN AT LEAST 19/39 BIOAEROSOL CONTROLS ### 
# REMOVE ASVS THAT OCCUR IN AT LEAST HALF OF CONTROLS (find which ASVs (i.e. rows) are present in at least half of controls?)
# First, get a subset with just the controls
colnames(airTaxOTUsOrderedReNamed)
justControlsSubset <- airTaxOTUsOrdered[1:nrow(airTaxOTUsOrdered)-1,c(1, 112:150)] %>% # columns 112:150 are the controls
  column_to_rownames(var="Row.names")
# View(justControlsSubset) #39 columns, as expected

# i. BRIEF DETOUR: which ASVs had any reads at all in the controls?
presASVsAirCtrl <- justControlsSubset[which(rowSums(justControlsSubset) > 0),] #804!
presASVsAirCtrl$ASVsum <- rowSums(presASVsAirCtrl)
presASVsAirCtrl <- merge(presASVsAirCtrl, taxall16S_20hrs, by=0)
nrow(justControlsSubset[which(rowSums(justControlsSubset) > 1),]) #604
nrow(justControlsSubset[which(rowSums(justControlsSubset) > 10),]) #414
nrow(justControlsSubset[which(rowSums(justControlsSubset) > 10000),]) #1

# ii.FOR LOOP TO COUNT the number of times a given ASV shows up in the controls and places this number in 
removalDataFrame <- as.data.frame(matrix(data=NA,nrow= nrow(justControlsSubset),ncol=3)) #5,278 rows = number of air ASVs
colnames(removalDataFrame) <- c("ASV_name", "controlInstances", "removalDecision")
removalDataFrame$ASV_name <- rownames(justControlsSubset)

for (i in 1:nrow(justControlsSubset)){ #looping over all 5,278 rows of ASVs
  removalDataFrame$controlInstances[i] <- sum(justControlsSubset[i,] > 0) #count number of times ASV appears in a column
  if (removalDataFrame$controlInstances[i] >= 19){ #19/39 = 49% of ASVs
    removalDataFrame$removalDecision[i] <- "remove ASV"
  } else {
    removalDataFrame$removalDecision[i] <- "retain ASV"
  }
}

# Some checks to show that the above for loop worked 
justControlsSubset[1,]
removalDataFrame[1,]
justControlsSubset[2,]
removalDataFrame[2,]
justControlsSubset[10,]
removalDataFrame[10,]
justControlsSubset[1334,]
removalDataFrame[1334,]

# Pull out only the ASVs in the removalDataFrame to remove
ASVstoRemoveNames <- removalDataFrame$ASV_name[which(removalDataFrame$removalDecision == "remove ASV")]
length(ASVstoRemoveNames) #10 
# ASVstoRemoveNames

# iii. MAKE NEW PHYLOSEQ OBJECT, MADE BY REMOVING CONTAMINANTS FROM WHOLE BIG DATASET
ASVsToKeepAll <- setdiff(rownames(otu_table(all16S_20hrs.ps)), ASVstoRemoveNames) #get all of the ASV names OTHER THAN THE REMOVE ASV
length(ASVsToKeepAll) #33,764 taxa (which is 10 less than earlier, 33,774 taxa in all16S_20hrs.ps)
all16S_20hrs_decontam.ps <- prune_taxa(ASVsToKeepAll, all16S_20hrs.ps)

# iv. FUTHER INVESTIGATION OF CONTAMINANTS:
# HOW MANY READS DO THE CONTROLS HAVE IF THESE CONTAMINANTS ARE REMOVED?
allControlsNames <- rownames(sample_data(all16S_20hrs_decontam.ps))[which(sample_data(all16S_20hrs_decontam.ps)$isControl == "control")]

# Air-associated controls: looking at how many reads post decontamination
decontamControlIndex <- which(colnames(otu_table(all16S_20hrs_decontam.ps))%in% airControlAllNames == TRUE)
decontamControls_ASVtab <- otu_table(all16S_20hrs_decontam.ps)[,decontamControlIndex] #get ASV table for just air controls
decontamControls_ASVtab_2 <- as.data.frame(as.matrix(decontamControls_ASVtab)) #one that isn't a phyloseq object
decontamControls_ASVtab_2$ASVsummedReads <- rowSums(decontamControls_ASVtab_2) #number of times each ASV appears in controls
# View(decontamControls_ASVtab_2)
setdiff(colnames(decontamControls_ASVtab), airControlAllNames) #shows that I subsetted all control samples correctly!
sort(colSums(decontamControls_ASVtab))
# Make a dataframe showing number of reads within each control sample post decontamination
readsControlsDeContam <- as.data.frame(matrix(nrow=length(colSums(decontamControls_ASVtab)), ncol=2))
colnames(readsControlsDeContam) <- c("sampName", "readNum_POST")
readsControlsDeContam[,1] <- names(colSums(decontamControls_ASVtab))
readsControlsDeContam[,2] <- unname(colSums(decontamControls_ASVtab))

# Show just controls with ASV_10 Staph in them after removing contamination
staphyControls_readsBeforeAfter_decontam <- left_join(staphyControlsOnly_colSums_df, readsControlsDeContam, by= "sampName")
staphyControls_readsBeforeAfter_decontam

# Get the before and after of just air controls before and after decontamination
airControlAllNames
readsControlsDeContam
airOTUs20h #this was made earlier in this script and has ASV table for each air sample pre-contam filtering (controls and samples) 
airCtrlpreDecontam <- airOTUs20h[,colnames(airOTUs20h) %in% airControlAllNames]
setdiff(colnames(airCtrlpreDecontam),airControlAllNames) #these are the same, as they should be
# Make a dataframe showing number of reads within each control sample post decontamination
airControlsPreDecontam <- as.data.frame(matrix(nrow=length(colSums(airCtrlpreDecontam)), ncol=2))
colnames(airControlsPreDecontam) <- c("sampName", "readNum_PRE")
airControlsPreDecontam[,1] <- names(colSums(airCtrlpreDecontam))
airControlsPreDecontam[,2] <- unname(colSums(airCtrlpreDecontam))
# Merge together pre and post decontamination 
airCtrlPrePostDecontam <- merge(airControlsPreDecontam, readsControlsDeContam, by= "sampName")
# Add in sample type information
airCtrlPrePostDecontam <- airCtrlPrePostDecontam %>%
  mutate(
    sampType = case_when(
      sampName %in% washBufferControlNames ~ "washBuffer",
      sampName %in% fieldControlNames ~ "fieldControl",
      sampName %in% kitomeAirControlNames ~ "kitome",
      sampName %in% PCRAirControlNames ~ "PCRAirControl",
      TRUE ~ NA_character_
    ) 
  ) %>% 
  arrange(desc(readNum_POST)) %>% 
  mutate(
    percReadsPOST = (readNum_POST/readNum_PRE)*100
  )
airCtrlPrePostDecontam
airCtrlPrePostDecontam %>% arrange(desc(percReadsPOST))
mean(airCtrlPrePostDecontam$percReadsPOST) #74.84696... this is not good

# Get the before and after of all (samples and controls included)
# Make a dataframe showing number of reads within each sample PRE decontamination
all16S_20hrs_decontamASVs <- as.data.frame(as.matrix(otu_table(all16S_20hrs_decontam.ps)))
airNames <- colnames(all16S_20hrs_decontamASVs)[grep(x=colnames(all16S_20hrs_decontamASVs), pattern = "air")] #looks correct
# Pull out only air samples (all samples, control and samples)
AIR16S_decontamASVs <- all16S_20hrs_decontamASVs[,grep(x=colnames(all16S_20hrs_decontamASVs), pattern = "air")] #pull out all rows but only air samples
# Make a dataframe for storing post-decontamination air
all_AIR_POST_Decontam <- as.data.frame(matrix(nrow=length(colSums(AIR16S_decontamASVs)), ncol=2))
colnames(all_AIR_POST_Decontam) <- c("sampName", "readNum_POST")
all_AIR_POST_Decontam[,1] <- names(colSums(AIR16S_decontamASVs))
all_AIR_POST_Decontam[,2] <- unname(colSums(AIR16S_decontamASVs))

# Make a dataframe showing number of reads within each sample PRE decontamination
airOTUs20h #this was made earlier in this script and has ASV table for each air sample pre-contam filtering (controls and samples) 
setdiff(colnames(airOTUs20h), airNames) #same as above

# Make a dataframe for storing PRE-decontamination air
all_AIR_PRE_Decontam <- as.data.frame(matrix(nrow=length(colSums(airOTUs20h)), ncol=2))
colnames(all_AIR_PRE_Decontam) <- c("sampName", "readNum_PRE")
all_AIR_PRE_Decontam[,1] <- names(colSums(airOTUs20h))
all_AIR_PRE_Decontam[,2] <- unname(colSums(airOTUs20h))

# Merge together pre and post decontamination 
all_AIR_PrePostDecontam <- merge(all_AIR_PRE_Decontam, all_AIR_POST_Decontam, by= "sampName")
# Add in sample type information
all_AIR_PrePostDecontam <- all_AIR_PrePostDecontam %>%
  mutate(
    sampType = case_when(
      sampName %in% washBufferControlNames ~ "washBuffer",
      sampName %in% fieldControlNames ~ "fieldControl",
      sampName %in% kitomeAirControlNames ~ "kitome",
      sampName %in% PCRAirControlNames ~ "PCRAirControl",
      TRUE ~ "sample"
    ) 
  ) %>% 
  arrange(desc(readNum_POST)) %>% 
  mutate(
    percReadsPOST = (readNum_POST/readNum_PRE)*100
  )
all_AIR_PrePostDecontam

all_AIR_PrePostDecontam %>% 
  group_by(sampType) %>% 
  summarise(medianPre = median(readNum_PRE))
# sampType      medianPre
# <chr>             <dbl>
#   1 PCRAirControl     1520.
# 2 fieldControl      3955 
# 3 kitome            2272.
# 4 sample           11317 
# 5 washBuffer        2813 

# GET SOME NUMBERS FOR PAPER:
colnames(all_AIR_PrePostDecontam)
# 1. Median % of reads in bioaerosol controls and samples that were present after decontamination
100-median(all_AIR_PrePostDecontam[which(all_AIR_PrePostDecontam$sampType != "sample"),5]) #get the median percReadsPOST for controls
# 21.23319 = 21.2% of control
100-median(all_AIR_PrePostDecontam[which(all_AIR_PrePostDecontam$sampType == "sample"),5]) #get the median percReadsPOST for samples
# 2.456829 = 2.5% of bioaerosol samples
# 2. Median # of reads per sample type
median(all_AIR_PrePostDecontam$readNum_POST[which(all_AIR_PrePostDecontam$sampType != "sample")]) #controls are 2,077
median(all_AIR_PrePostDecontam$readNum_POST[which(all_AIR_PrePostDecontam$sampType == "sample")]) #samples are 10,583

# Written Sept 8, 2025 -- original for the original rule to remove an ASV if in more than half of controls
# write.csv(all_AIR_PrePostDecontam, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/all_AIR_PrePostDecontam_original.csv")

all_AIR_PrePostDecontam %>% 
  group_by(sampType) %>% 
  summarize(medianPre = median(readNum_PRE)) 
# sampType      medianPre
# <chr>             <dbl>
#   1 PCRAirControl     1520.
# 2 fieldControl      3955 
# 3 kitome            2272.
# 4 sample           11317 
# 5 washBuffer        2813 

all_AIR_PrePostDecontam %>% 
  group_by(sampType) %>% 
  summarize(medianPost = median(readNum_POST))
# sampType      medianPost
# <chr>              <dbl>
#   1 PCRAirControl       690.
# 2 fieldControl       3456.
# 3 kitome             1648 
# 4 sample            10583 
# 5 washBuffer         1899 

## HOW ABOUT overlap b/w contaminant ASVs and most common in controls
colnames(presASVsAirCtrl)[colnames(presASVsAirCtrl) %in% "Row.names"] <- "ASV_name"
colnames(presASVsAirCtrl)
# Add how many occurrences of each ASV
for (i in 1:nrow(presASVsAirCtrl)){
  presASVsAirCtrl$nSampOcc[i] <- sum(presASVsAirCtrl[i,2:40] > 0)
}
# View(presASVsAirCtrl)
# Re-order:
colnames(presASVsAirCtrl)
presASVsAirCtrl <- presASVsAirCtrl[,c(1:41,48,42:47)]

# View(presASVsAirCtrl)
# Add sample type
presASVsAirCtrl
sampTypesDf <- as.data.frame(matrix(ncol=48, nrow=1)) #48 to match presASVsAirCtrl
colnames(sampTypesDf) <- c("ASV_name", washBufferControlNames, fieldControlNames, kitomeAirControlNames, PCRAirControlNames, 
                           "ASVsum","nSampOcc", "Kingdom", "Phylum", "Class", "Order","Family","Genus")
sampTypesDf[,] <- c("ctrlType", rep("washBuffer", length(washBufferControlNames)), rep("fieldCtrl", length(fieldControlNames)),
                    rep("kitome", length(kitomeAirControlNames)), rep("PCRAirCtrl", length(PCRAirControlNames)), rep(NA, 8)) #NA 8 times for extra colums on end of presASVsAirCtrl

# Re-order presASVsAirCtrl to match sampTypesDf
colnames(presASVsAirCtrl)
presASVsAirCtrl <- presASVsAirCtrl[ , c("ASV_name",washBufferControlNames, fieldControlNames, kitomeAirControlNames, PCRAirControlNames,
                                        "ASVsum","nSampOcc", "Kingdom", "Phylum", "Class", "Order","Family","Genus")]
colnames(presASVsAirCtrl) == colnames(sampTypesDf) #since all true can bind below
head(presASVsAirCtrl)
# Arrange by ASVsum BEFORE binding below, when it will no longer be able to be sorted
presASVsAirCtrl <- presASVsAirCtrl %>% 
  arrange(desc(ASVsum))
# Bind together
ASVsAirCtrl_df <- rbind(sampTypesDf, presASVsAirCtrl)
# View(ASVsAirCtrl_df)
colnames(ASVsAirCtrl_df)
ASVsAirCtrl_df[2:nrow(ASVsAirCtrl_df),c(2:42)] <- as.numeric(ASVsAirCtrl_df[2:nrow(ASVsAirCtrl_df),c(2:42)])
# This is no good b/c then can't order numerically.
# September 10, 2025
# write.csv(ASVsAirCtrl_df, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ASVsAirCtrl_df.csv")

presASVsAirCtrl_25 <- presASVsAirCtrl %>% 
  arrange(desc(ASVsum)) %>% 
  slice_head(n= 25)
# Half of top 10 most abundant weren't even all removed-- almost definitely due to cross contamination from real samples at PCR step
setdiff(presASVsAirCtrl_25$ASV_name[1:10], ASVstoRemoveNames) #"ASV_90"  "ASV_485" "ASV_269" "ASV_301" "ASV_638" "ASV_69" 

# SAVED IN ORIGINAL full16S_EDA_part1July12.R, AS DOCUMENTED BELOW
# written and saved July 5, 2023- airTaxOTUsOrderedDecontam.csv
# Below written and saved July 12, 2023 on personal computer 
# write.csv(airTaxOTUsOrderedReNamed, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airTaxOTUsOrderedDecontam_July12.csv")


# ## THIS LOOKS AT HOW OFTEN THE MOST COMMON AIR ASVS APPEAR IN CONTROLS, AND BECAUSE IT SORT OF BORROWS FROM A LATER SCRIPT,
# # I WILL GRAY IT OUT FOR NOW, EVEN THOUGH IT IS VERY USEFUL!
# top100bactTax_TableS4 <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/manuscript/I6S_topAirTable_final_Dec18.csv", header=TRUE)
# # View(top100bactTax_TableS4)
# top100bactTax_wCtrls <- left_join(top100bactTax_TableS4, presASVsAirCtrl, by="ASV_name")
# colnames(top100bactTax_wCtrls)
# top100bactTax_wCtrls2 <- top100bactTax_wCtrls[,c(2:13,53,54)]
# # View(top100bactTax_wCtrls2)
# length(which(is.na(top100bactTax_wCtrls2$ASVsum))) #44 are not found at all in the controls
# length(which(top100bactTax_wCtrls2$ASVsum < 5)) #18, so 18+44 = 62/100
# length(which(top100bactTax_wCtrls2$ASVsum < 80)) #18, so 30+44 = 74/100 are found on average less than 2 times per control
# length(which(top100bactTax_wCtrls2$ASVsum < 400)) # 37+44 = 81/100. So this is anything less than 10 times on average.
# 
# # Clean this up
# colnames(top100bactTax_wCtrls2)[colnames(top100bactTax_wCtrls2) %in% "ASVsum"] <- "ctrlASVsum"
# colnames(top100bactTax_wCtrls2)[colnames(top100bactTax_wCtrls2) %in% "nSampOcc"] <- "nCtrlSampOcc"
# colnames(top100bactTax_wCtrls2)[colnames(top100bactTax_wCtrls2) %in% "Phylum.x"] <- "Phylum"
# colnames(top100bactTax_wCtrls2)[colnames(top100bactTax_wCtrls2) %in% "Class.x"] <- "Class"
# colnames(top100bactTax_wCtrls2)[colnames(top100bactTax_wCtrls2) %in% "Order.x"] <- "Order"
# colnames(top100bactTax_wCtrls2)[colnames(top100bactTax_wCtrls2) %in% "Family.x"] <- "Family"
# colnames(top100bactTax_wCtrls2)[colnames(top100bactTax_wCtrls2) %in% "Genus.x"] <- "Genus"
# # replace NAs with 0s
# top100bactTax_wCtrls2$ctrlASVsum[which(is.na(top100bactTax_wCtrls2$ctrlASVsum))] <- 0
# top100bactTax_wCtrls2$nCtrlSampOcc[which(is.na(top100bactTax_wCtrls2$nCtrlSampOcc))] <- 0
# # View(top100bactTax_wCtrls2)
# # Save to send to discuss with Noah (saved Sept 10, 2025)
# # write.csv(top100bactTax_wCtrls2, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/top100bactTax_wCtrls2.csv")

##################################################################
# DATA CLEANING PART V
# (dealing with phyllosphere-associated sample contamination)
# New to this script
##################################################################
#### 2. FOLIAR CONTROLS #####
# Pull out just foliar surface samples -- note that 2 PCR NTCs for phyllo failed (as shown above)
unique(sample_data(all16S_20hrs.ps)$sampleType) # sampleType == "PCRcontrol" are only air since phyllo PCR controls failed
# 1. GET DATA SUBSETS
## i. all phyllo and phyllo controls
phyllo16S_20hrs.ps <- subset_samples(all16S_20hrs.ps, sampleType %in% c("phyllosphere", "phyllo_NegFieldControl", "phyllo_ExtBlank"))
## ii. Just samples
samps_phyllo16S_20hrs.ps <- subset_samples(phyllo16S_20hrs.ps, isControl == "sample")
sort(colSums(otu_table(samps_phyllo16S_20hrs.ps ))) #reads in each
median(colSums(otu_table(samps_phyllo16S_20hrs.ps )) ) #24,177

## iii. Just controls
Ctrl_phyllo16S_20hrs.ps <- subset_samples(phyllo16S_20hrs.ps, isControl == "control")
sort(colSums(otu_table(Ctrl_phyllo16S_20hrs.ps))) #reads in each
# phyllo_16S_66 phyllo_16S_52  phyllo_16S_6 phyllo_16S_28 phyllo_16S_27 phyllo_16S_51 phyllo_16S_65 
# 270           685           733          1333          6919          7924         10006 
median(colSums(otu_table(Ctrl_phyllo16S_20hrs.ps)) ) #1,333, added to manuscript

# View(as.data.frame(as.matrix(otu_table(Ctrl_phyllo16S_20hrs.ps))))
length(which(rowSums(otu_table(Ctrl_phyllo16S_20hrs.ps))>0)) #157 ASVs were found in any abundance in these controls
length(which(rowSums(otu_table(Ctrl_phyllo16S_20hrs.ps))>5)) #74 ASVs were found more than 5 times in these controls
length(which(rowSums(otu_table(Ctrl_phyllo16S_20hrs.ps))>100)) #40 ASVs were found more than 100 times in these controls

# 1. IDENTIFY ANY ASVS THAT WERE DETECTED IN AT LEAST 4 CONTROLS (4/ORIGINAL 9, WHICH IS NOW 7)
ctrlsPhyllo_df <- as.data.frame(as.matrix(otu_table(Ctrl_phyllo16S_20hrs.ps)))
head(ctrlsPhyllo_df)
# Subject phyllo to same procedure
for (i in 1:nrow(ctrlsPhyllo_df)) {
  ctrlsPhyllo_df$controlInstances[i] <- sum(ctrlsPhyllo_df[i, 1:7] > 0) #columns 1:7 because these are the samples
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
# Above shows "ASV_11"  "ASV_110" "ASV_132" "ASV_79"  "ASV_88". Overlaps with air removal taxa (79, 88, 110), so unique to foliar are "ASV_11" and  "ASV_132" 
ASVsToRemovePhyllo <- sort(rownames(ctrlsPhyllo_df)[which(ctrlsPhyllo_df$removalDecision == "remove ASV")]) #these are the bad phyllosphere taxa

# Add in taxonomy to View(ctrlsPhyllo_df)
phylloTax <- as.data.frame(as.matrix(tax_table(phyllo16S_20hrs.ps)))
head(phylloTax); head(ctrlsPhyllo_df)
ctrlsPhylloWtax_df <- merge(ctrlsPhyllo_df, phylloTax, by=0)
ctrlsPhylloWtax_df[which(ctrlsPhylloWtax_df$removalDecision== "remove ASV"),"Genus"]
# [1] "Methylobacterium-Methylorubrum"             "Burkholderia-Caballeronia-Paraburkholderia" "Staphylococcus"                            
# [4] "Geobacillus"                                "Tetragenococcus" 
# ASV 178 is Staphylococcus, not the same as above (ASV_10)

# BEFORE REMOVAL, WHAT PERCENTAGE OF TOTAL READS IN SAMPS AND CONTROLS DID THESE CONTAMINANTS COMPRISE?
# Samples
badPhylloTaxCounts_samps <- colSums(otu_table(samps_phyllo16S_20hrs.ps)[rownames(otu_table(samps_phyllo16S_20hrs.ps)) %in% ASVsToRemovePhyllo,])
allPhylloTaxCounts_samps <- colSums(otu_table(samps_phyllo16S_20hrs.ps))
names(badPhylloTaxCounts_samps) == names(allPhylloTaxCounts_samps) #okay so good to divide!
median((badPhylloTaxCounts_samps/allPhylloTaxCounts_samps)*100) #1.902227, 1.9% added to manuscript

# Controls
badPhylloTaxCounts_ctrls <- colSums(otu_table(Ctrl_phyllo16S_20hrs.ps)[rownames(otu_table(Ctrl_phyllo16S_20hrs.ps)) %in% ASVsToRemovePhyllo,])
allPhylloTaxCounts_ctrls <- colSums(otu_table(Ctrl_phyllo16S_20hrs.ps))
names(badPhylloTaxCounts_ctrls) == names(allPhylloTaxCounts_ctrls) #okay so good to divide!
median((badPhylloTaxCounts_ctrls/allPhylloTaxCounts_ctrls)*100) #6.396162, 6.4% added to manuscript

##################################################################
# DATA CLEANING PART VI
# New phyloseq objects
##################################################################
# MAKE NEW PHYLOSEQ OBJECT, MADE BY REMOVING AIR ASSOCIATED CONTAMINANTS FROM AIR AND PHYLLO FROM PHYLLO
ASVstoRemoveNames #these are the 'bad' air taxa
ASVsToRemovePhyllo #bad phyllosphere

# 1. Remove from bioaerosols
I6S_airCtrlSamps_20hrs.ps #just air + all air-associated controls, also made from all16S_20hrs.ps
min(rowSums(otu_table(I6S_airCtrlSamps_20hrs.ps))) #shows all ASVs have at least some reads, as we want
# Get air-associated ASVs to keep by finding the ASVs that differ b/w all and removal ASVs
goodAirASVs <- setdiff(rownames(otu_table(I6S_airCtrlSamps_20hrs.ps)), ASVstoRemoveNames)
# Double check
length(goodAirASVs) == length(rownames(otu_table(I6S_airCtrlSamps_20hrs.ps))) - 10 #this is true
unique(goodAirASVs %in% ASVstoRemoveNames) #great these are not the same
# New phyloseq object for just air that removes them
I6S_decontam_justAir.ps <- prune_taxa(taxa= goodAirASVs, x= I6S_airCtrlSamps_20hrs.ps) #new one does have 10 fewer, as expected

# 2. Remove from foliar surfaces
phyllo16S_20hrs.ps #just phyllo samples + control dataset, made from all16S_20hrs.ps
# Remove those taxa without reads
phyllo16S_20hrs.ps <- prune_taxa(taxa_sums(phyllo16S_20hrs.ps) > 0, phyllo16S_20hrs.ps) #keep only those with at least 1 read
min(rowSums(otu_table(phyllo16S_20hrs.ps))) #shows all ASVs have at least some reads, as we want
# Get foliar surface-associated ASVs to keep by finding the ASVs that differ b/w all and removal ASVs
goodFoliarASVs <- setdiff(rownames(otu_table(phyllo16S_20hrs.ps)), ASVsToRemovePhyllo)
# Double check
length(goodFoliarASVs) == length(rownames(otu_table(phyllo16S_20hrs.ps))) - 5 #this is true
unique(goodFoliarASVs %in% ASVsToRemovePhyllo) #great these are not the same
# New phyloseq object for just foliar surface that removes them
I6S_decontam_justFoliar.ps <- prune_taxa(taxa = goodFoliarASVs, phyllo16S_20hrs.ps) #7091 taxa, 5 less than 7096 in original
#. Number of reads in foliar surface controls after this:
allPhylloCtrlNames <- rownames(sample_data(I6S_decontam_justFoliar.ps))[sample_data(I6S_decontam_justFoliar.ps)$isControl == "control"]
median(colSums(otu_table(I6S_decontam_justFoliar.ps)[,colnames(otu_table(I6S_decontam_justFoliar.ps)) %in% allPhylloCtrlNames])) #1,083 added to manuscript
#. Number of reads in foliar surface SAMPLES after this:
allPhylloSampNames <- rownames(sample_data(I6S_decontam_justFoliar.ps))[sample_data(I6S_decontam_justFoliar.ps)$isControl == "sample"]
median(colSums(otu_table(I6S_decontam_justFoliar.ps)[,colnames(otu_table(I6S_decontam_justFoliar.ps)) %in% allPhylloSampNames])) #23,587 added to paper

# 3. Bring back in soil taxa (will drop soil controls, as we did for initial paper, since these were so rare compared to soil samples)
which(sample_data(all16S_20hrs.ps)$sampleType == "soil")
soil_20hrs.ps <- subset_samples(all16S_20hrs.ps, sampleType == "soil")
soil_20hrs.ps <- prune_taxa(taxa_sums(soil_20hrs.ps) > 0, soil_20hrs.ps) #keep only those with at least 1 read. This dropped 33774-25285 = 8489 ASVs

# Make new, consensus phyloseq object
I6Sall_d_foliarAir.ps <- merge_phyloseq(I6S_decontam_justAir.ps, I6S_decontam_justFoliar.ps, soil_20hrs.ps)
airCheckPS_ASVs <- otu_table(I6Sall_d_foliarAir.ps)[,sample_data(I6Sall_d_foliarAir.ps)$sampleType == "air"]
rowSums(airCheckPS_ASVs[which(rownames(airCheckPS_ASVs) %in% ASVstoRemoveNames),]) #great, these are not present in air samples!
foliarCheckPS_ASVs <- otu_table(I6Sall_d_foliarAir.ps)[,sample_data(I6Sall_d_foliarAir.ps)$sampleType == "phyllosphere"]
rowSums(foliarCheckPS_ASVs[which(rownames(foliarCheckPS_ASVs) %in% ASVsToRemovePhyllo),]) #great, these are not present in foliar samples!

# Saved September 15, 2025
# save(I6Sall_d_foliarAir.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6Sall_d_foliarAir_ps_Sept15") 

##################################################################
# RAREFYING (EXPLORING, THEN RAREFYING AT 5.5K READS TO RETAIN (air) 84 SAMPLES AND 8 CONTROLS)
##################################################################
### ORIGINAL, WITH NO PHYLLOSPHERE CONTROLS REMOVED ### 
#i.  EXPLORE READS ACROSS SAMPLES
colnames(otu_table(all16S_20hrs_decontam.ps))
all16S_seqspersample <- sort(colSums(otu_table(all16S_20hrs_decontam.ps)))
all16S_seqspersample
airSeqsIndex <- grep(x= names(all16S_seqspersample), pattern= "air")
length(airSeqsIndex) == length(airAllIndex) #149. As expected 
all16S_SeqNumberPlot <-barplot(all16S_seqspersample, main="16S: Total Sequences Per Sample", ylab= "Number of Sequences", xlab= "Sample Number")
air16S_seqsPerSample <- all16S_seqspersample[airSeqsIndex]

# NOTE: the number in these indices are indices NOT number of seqs
airDrop5000 <- which(air16S_seqsPerSample < 5000) # drops 54 samples/blanks
length(airDrop5000) #54 sample
airDrop5500 <- which(air16S_seqsPerSample < 5500) # drops 57 samples/blanks
length(airDrop5500) #57 sample
airDrop6000 <- which(air16S_seqsPerSample < 6000) # drops 60 samples/blanks
length(airDrop6000) #60 samples
airDrop8000 <- which(air16S_seqsPerSample < 8000) # drops 75 samples/blanks
length(airDrop8000) #75 samples

airControlAllNames #made above near where I was making airTaxOTUsOrdered
# WHICH SAMPLES WOULD BE DROPPED IF USED THE FOLLOWING CRITERIA??
# 8K
length(which(names(airDrop8000) %in% airControlAllNames == TRUE)) #33 of these 75 dropped are controls. 
length(airDrop8000)-33 #this would drop 42 samples!! 

# 6K
length(which(names(airDrop6000) %in% airControlAllNames == TRUE)) #32 of these 60 dropped are controls. 
length(airDrop6000)-32 #this would drop 28 samples!! 
110 -28 #would only have 82 samples left :(

# 5K
length(which(names(airDrop5000) %in% airControlAllNames == TRUE)) #31 of these 54 dropped are controls. 
length(airDrop5000)-31 #this would drop 23 samples!! 

# 5.5K
length(which(names(airDrop5500) %in% airControlAllNames == TRUE)) #31 of these 54 dropped are controls. 
length(airDrop5500)-31 #this would drop 26 samples!! 
110 -26 #84 samples left 

names(airDrop5500)
samplesToDrop <- names(airDrop5500) #names of the samples to drop-- i.e. rarefying at 5.5K reads
length(names(air16S_seqsPerSample)) - length(samplesToDrop) #this will leave us with 92 total samples/blanks,
# including 84 samples and 8 controls

# saved Aug 4, 2023 (on personal computer, not server) in order to drop these in chloroplast exploration (see chloroplast_16SASVs_EDA.R)
#saveRDS(samplesToDrop, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/RObjectsSaved/samplesToDrop550016Sreads")

# Get sample data of the samples to be dropped 
droppedMetaData <- sample_data(all16S_20hrs_decontam.ps)[which(rownames(sample_data(all16S_20hrs_decontam.ps)) %in% samplesToDrop == TRUE),]
length(which(droppedMetaData$sampleType == "air")) #26 are air samples
length(which(droppedMetaData$sampleType == "fieldControl")) #10 are field controls
length(which(droppedMetaData$sampleType == "washBufferControl")) #9 are wash buffer controls
length(which(droppedMetaData$sampleType == "kitome")) #10 are kitome
length(which(droppedMetaData$sampleType == "PCRcontrol")) #2 are PCRcontrol (both)

# Samples that will be retained after rarefying!
samplesToKeepNames <- names(which(all16S_seqspersample >= 5500)) #get all of these ACROSS sample types
samplesToKeepNames #

set.seed(19)
all16S_rarefied.ps <- rarefy_even_depth(all16S_20hrs_decontam.ps, sample.size = 5500, replace = FALSE, trimOTUs = TRUE)
# 3,468 ASVs were removed because they are no longer present in any sample after random subsampling. 77 samples removed (across all types)
length(which(sample_data(all16S_rarefied.ps)$sampleType== "air")) #84 air samples (not controls)
length(which(sample_data(all16S_rarefied.ps)$sampleType== "fieldControl")) #6 field controls
length(which(sample_data(all16S_rarefied.ps)$sampleType== "washBufferControl")) #2 field controls
length(which(sample_data(all16S_rarefied.ps)$sampleType== "phyllo_NegFieldControl")) #3 phyllosphere field controls

# saved July 12, 2023 on personal computer, saved June 10, 2024 on server
# save(all16S_rarefied.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16S_rarefied_phyloseq_July12") 
# save(all16S_rarefied.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/all16S_rarefied_phyloseq_July12") 
# Saved again Sept. 15, 2025 for double checking that it's the same.
# save(all16S_rarefied.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/all16S_rarefied_phyloseq_Sept15_2025") 

### SECOND WAY WITH PHYLLOSPHERE CONTROLS REMOVED ### 
# Will remove forward with this version
set.seed(19)
I6S_dcAP_rarefied.ps <- rarefy_even_depth(I6Sall_d_foliarAir.ps, sample.size = 5500, replace = FALSE, trimOTUs = TRUE)
# 3161OTUs were removed because they are no longer present in any sample after random subsampling

table((sample_data(I6S_dcAP_rarefied.ps)$sampleType))

# air           fieldControl phyllo_NegFieldControl           phyllosphere 
# 84                      6                      3                     58 
# soil      washBufferControl 
# 157                      2 


# Saved September 15, 2025
# save(I6S_dcAP_rarefied.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_dcAP_rarefied_ps_Sept15") 

##################################################################
# CONTAMINANTS IN BIOAEROSOLS AND FOLIAR SURFACES (TABLES)
##################################################################
#### DEEPER LOOK AT ASVS THAT WERE CONTAMINANTS IN BIOAEROSOLS. ####
# ALSO, MAKE CLEANED UP TABLE OF THESE FOR SUPPLEMENT (TABLE S2)
# Add a column with the mean and range or SD of the percent of total reads in the samples that these ASVs comprise
tableS2_I6S_air_contams_1 <- as.data.frame(tax_table(all16S_20hrs.ps)[which(rownames(tax_table(all16S_20hrs.ps)) %in% ASVstoRemoveNames == TRUE),])

# Get only bioaerosol samples that were retained after rarefying and then look at the original amounts that they had of these ASVs
keptAirNames <- rownames(sample_data(I6S_dcAP_rarefied.ps))[sample_data(I6S_dcAP_rarefied.ps)$sampleType== "air"] #84 air samples (not controls)
all16S_20hrsASVs <- as.data.frame(as.matrix(otu_table(all16S_20hrs.ps)))
which(rownames(all16S_20hrsASVs) %in% ASVstoRemoveNames==TRUE)
which(colnames(all16S_20hrsASVs) %in% keptAirNames == TRUE)
airI6S_pre_ASVs <- all16S_20hrsASVs[, which(colnames(all16S_20hrsASVs) %in% keptAirNames == TRUE)]
airI6S_pre_totalReads <- colSums(airI6S_pre_ASVs) #these are the number of reads before rarefying and contaminant removal
# Now for each of these 10 ASVs, want abundance in each air sample
airI6S_contam_ASVs <- all16S_20hrsASVs[which(rownames(all16S_20hrsASVs) %in% ASVstoRemoveNames==TRUE), which(colnames(all16S_20hrsASVs) %in% keptAirNames == TRUE)]

# FOR LOOP TO GET CALCULATIONS
# preallocate dataframe with contam ASVs as rows and air sample names as columns
contamInAir <- as.data.frame(matrix(nrow=length(ASVstoRemoveNames), ncol=length(keptAirNames)))
rownames(contamInAir) <- ASVstoRemoveNames
colnames(contamInAir) <- keptAirNames

# For loop divides number of reads that a sample [i] has for contam ASV [j] and divides this by total reads sample [i] has
# and then places this in the pre-allocated dataframe in the row of the contam ASV [j] and the column for sample [i]
for (i in seq_along(keptAirNames)){ #loop over all of the samples 84 samples in keptAirNames/airI6S_contam_ASVs
  for (j in seq_along(ASVstoRemoveNames)){ # loop over all of the contaminant ASVs in ASVstoRemoveNames/airI6S_contam_ASVs
    contamInAir[rownames(contamInAir) %in% ASVstoRemoveNames[j], colnames(contamInAir) %in% keptAirNames[i]] <-
      airI6S_contam_ASVs[rownames(airI6S_contam_ASVs) %in% ASVstoRemoveNames[j],colnames(airI6S_contam_ASVs) %in% keptAirNames[i]]/airI6S_pre_totalReads[names(airI6S_pre_totalReads) %in% keptAirNames[i]]
  } 
} 
contamInAir

# DOUBLE CHECK FOR LOOP
# De-construct nested for loop above to make sure calculations working as expected 
airI6S_contam_ASVs[rownames(airI6S_contam_ASVs) %in% ASVstoRemoveNames[1],] #pull out all ASV_110
airI6S_contam_ASVs[,colnames(airI6S_contam_ASVs) %in% keptAirNames[1]] #pull out all of the first air sample (this is all ASV count for col 1, "air_16S_1")
# Put together-- yes this did it correctly! ASV_110 for air sample 1 is 0 in airI6S_contam_ASVs
airI6S_contam_ASVs[rownames(airI6S_contam_ASVs) %in% ASVstoRemoveNames[1],colnames(airI6S_contam_ASVs) %in% keptAirNames[1]] #pull out all ASV_110 for "air_16S_1"
# Get colSums for this sample
airI6S_pre_totalReads[names(airI6S_pre_totalReads) %in% keptAirNames[1]] #this is working correctly
airI6S_contam_ASVs[rownames(airI6S_contam_ASVs) %in% ASVstoRemoveNames[1],colnames(airI6S_contam_ASVs) %in% keptAirNames[1]]/airI6S_pre_totalReads[names(airI6S_pre_totalReads) %in% keptAirNames[1]]
# this is zero, as we wanted

# Now just check calculations --looking good!
contamInAir[rownames(contamInAir) %in% ASVstoRemoveNames[1], colnames(contamInAir) %in% keptAirNames[1]] == airI6S_contam_ASVs[rownames(airI6S_contam_ASVs) %in% ASVstoRemoveNames[1],colnames(airI6S_contam_ASVs) %in% keptAirNames[1]]/airI6S_pre_totalReads[names(airI6S_pre_totalReads) %in% keptAirNames[1]]
contamInAir[rownames(contamInAir) %in% ASVstoRemoveNames[5], colnames(contamInAir) %in% keptAirNames[10]] == airI6S_contam_ASVs[rownames(airI6S_contam_ASVs) %in% ASVstoRemoveNames[5],colnames(airI6S_contam_ASVs) %in% keptAirNames[10]]/airI6S_pre_totalReads[names(airI6S_pre_totalReads) %in% keptAirNames[10]]
# Double check the above even more manually for ASV_2292 and air_16S_116
ASVstoRemoveNames[5]; keptAirNames[10] 
airI6S_pre_totalReads[which(names(airI6S_pre_totalReads)== keptAirNames[10])] #this is air_16S_116
airI6S_contam_ASVs[which(rownames(airI6S_contam_ASVs)== ASVstoRemoveNames[5]), which(colnames(airI6S_contam_ASVs)== keptAirNames[10])]/
  airI6S_pre_totalReads[which(names(airI6S_pre_totalReads)== keptAirNames[10])] == contamInAir[rownames(contamInAir) %in% ASVstoRemoveNames[5], colnames(contamInAir) %in% keptAirNames[10]]

# FINISH MAKING TABLE
# What is median %abundance in each of these samples for each contaminant ASV
# 1. Make it longer!
contamInAir_longer <- contamInAir %>% 
  rownames_to_column(var="ASV_name") %>% 
  pivot_longer(cols = starts_with("air_"),
               names_to = "sampName",
               values_to = "propInSamp")

# 2. Group by ASV and then get median percentage
contamInAir_longer <- contamInAir_longer %>% 
  group_by(ASV_name) %>% 
  mutate(medianPercentage = median(propInSamp)*100) %>% #does not matter if median taken before or after *100
  mutate(minPercentage = min(propInSamp)*100) %>% 
  mutate(maxPercentage = max(propInSamp)*100) %>% 
  mutate(sdPercentage = sd(propInSamp)*100)

contamInAir[which(rownames(contamInAir) %in% "ASV_110"),] # a lot of zeros indeed!

# Finally, clean this up to make it the table
colnames(contamInAir_longer)
# Make ASV name a column and then below, remove Kingdom since it's not necessary
tableS2_I6S_air_contams_1 <- tableS2_I6S_air_contams_1 %>% 
  rownames_to_column(var="ASV_name")
tableS2_I6S_air_contams_2 <- left_join(tableS2_I6S_air_contams_1[,c(1,3:ncol(tableS2_I6S_air_contams_1))] ,contamInAir_longer[,c(1,4,7)], by = "ASV_name")
tableS2_I6S_air_contams_final <- tableS2_I6S_air_contams_2 %>% 
  distinct() #remove duplicate rows which were created by pivoting longer above
# write the .csv to finish this
# Written September 15, 2025
# write.csv(tableS2_I6S_air_contams_final, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/manuscript/tableS2_I6S_air_contams_final_Sept15.csv")

#### DEEPER LOOK AT ASVS THAT WERE CONTAMINANTS IN FOLIAR SURFACES. ####
# ALSO, MAKE CLEANED UP TABLE OF THESE FOR SUPPLEMENT (TABLE S3, I think)
# Add a column with the mean and range or SD of the percent of total reads in the samples that these ASVs comprise
tableS2_I6S_foliar_contams_1 <- as.data.frame(tax_table(all16S_20hrs.ps)[which(rownames(tax_table(all16S_20hrs.ps)) %in% ASVsToRemovePhyllo == TRUE),])

# Get only bioaerosol samples that were retained after rarefying and then look at the original amounts that they had of these ASVs
keptPhylloNames <- rownames(sample_data(I6S_dcAP_rarefied.ps))[sample_data(I6S_dcAP_rarefied.ps)$sampleType== "phyllosphere"] #58 phyllo samples (not controls)
which(rownames(all16S_20hrsASVs) %in% ASVsToRemovePhyllo==TRUE)
which(colnames(all16S_20hrsASVs) %in% keptPhylloNames == TRUE)
phylloI6S_pre_ASVs <- all16S_20hrsASVs[, which(colnames(all16S_20hrsASVs) %in% keptPhylloNames == TRUE)]
phylloI6S_pre_totalReads <- colSums(phylloI6S_pre_ASVs) #these are the number of reads before rarefying and contaminant removal
# Now for each of these 5 ASVs, want abundance in each phyllo sample
phylloI6S_contam_ASVs <- all16S_20hrsASVs[which(rownames(all16S_20hrsASVs) %in% ASVsToRemovePhyllo==TRUE), which(colnames(all16S_20hrsASVs) %in% keptPhylloNames == TRUE)]
dim(phylloI6S_contam_ASVs) #5, 58, as expected

# FOR LOOP TO GET CALCULATIONS
# preallocate dataframe with contam ASVs as rows and phyllo sample names as columns
contamInphyllo <- as.data.frame(matrix(nrow=length(ASVsToRemovePhyllo), ncol=length(keptPhylloNames)))
rownames(contamInphyllo) <- ASVsToRemovePhyllo
colnames(contamInphyllo) <- keptPhylloNames

# For loop divides number of reads that a sample [i] has for contam ASV [j] and divides this by total reads sample [i] has
# and then places this in the pre-allocated dataframe in the row of the contam ASV [j] and the column for sample [i]
for (i in seq_along(keptPhylloNames)){ #loop over all of the samples 84 samples in keptPhylloNames/phylloI6S_contam_ASVs
  for (j in seq_along(ASVsToRemovePhyllo)){ # loop over all of the contaminant ASVs in ASVsToRemovePhyllo/phylloI6S_contam_ASVs
    contamInphyllo[rownames(contamInphyllo) %in% ASVsToRemovePhyllo[j], colnames(contamInphyllo) %in% keptPhylloNames[i]] <-
      phylloI6S_contam_ASVs[rownames(phylloI6S_contam_ASVs) %in% ASVsToRemovePhyllo[j],colnames(phylloI6S_contam_ASVs) %in% keptPhylloNames[i]]/phylloI6S_pre_totalReads[names(phylloI6S_pre_totalReads) %in% keptPhylloNames[i]]
  } 
} 
contamInphyllo

# DOUBLE CHECK FOR LOOP
# De-construct nested for loop above to make sure calculations working as expected 
phylloI6S_contam_ASVs[rownames(phylloI6S_contam_ASVs) %in% ASVsToRemovePhyllo[1],] #pull out all ASV_11
phylloI6S_contam_ASVs[,colnames(phylloI6S_contam_ASVs) %in% keptPhylloNames[1]] #pull out all of the first phyllo sample (this is all ASV count for col 1, "phyllo_16S_1")
# Put together-- yes this did it correctly! ASV_11 for phyllo sample 1 is 385 in phylloI6S_contam_ASVs
phylloI6S_contam_ASVs[rownames(phylloI6S_contam_ASVs) %in% ASVsToRemovePhyllo[1],colnames(phylloI6S_contam_ASVs) %in% keptPhylloNames[1]] #pull out all ASV_11 for "phyllo_16S_1"
# Get colSums for this sample
phylloI6S_pre_totalReads[names(phylloI6S_pre_totalReads) %in% keptPhylloNames[1]] #this is working correctly
phylloI6S_contam_ASVs[rownames(phylloI6S_contam_ASVs) %in% ASVsToRemovePhyllo[1],colnames(phylloI6S_contam_ASVs) %in% keptPhylloNames[1]]/phylloI6S_pre_totalReads[names(phylloI6S_pre_totalReads) %in% keptPhylloNames[1]]
# this is 0.01863504 , as we wanted (same as contamInphyllo[1,1])

# Now just check calculations --looking good!
contamInphyllo[rownames(contamInphyllo) %in% ASVsToRemovePhyllo[1], colnames(contamInphyllo) %in% keptPhylloNames[1]] == phylloI6S_contam_ASVs[rownames(phylloI6S_contam_ASVs) %in% ASVsToRemovePhyllo[1],colnames(phylloI6S_contam_ASVs) %in% keptPhylloNames[1]]/phylloI6S_pre_totalReads[names(phylloI6S_pre_totalReads) %in% keptPhylloNames[1]]
contamInphyllo[rownames(contamInphyllo) %in% ASVsToRemovePhyllo[5], colnames(contamInphyllo) %in% keptPhylloNames[10]] == phylloI6S_contam_ASVs[rownames(phylloI6S_contam_ASVs) %in% ASVsToRemovePhyllo[5],colnames(phylloI6S_contam_ASVs) %in% keptPhylloNames[10]]/phylloI6S_pre_totalReads[names(phylloI6S_pre_totalReads) %in% keptPhylloNames[10]]
# Double check the above even more manually for ASV_2292 and phyllo_16S_116
ASVsToRemovePhyllo[5]; keptPhylloNames[10] 
phylloI6S_pre_totalReads[which(names(phylloI6S_pre_totalReads)== keptPhylloNames[10])] #this is phyllo_16S_19 
phylloI6S_contam_ASVs[which(rownames(phylloI6S_contam_ASVs)== ASVsToRemovePhyllo[5]), which(colnames(phylloI6S_contam_ASVs)== keptPhylloNames[10])]/
  phylloI6S_pre_totalReads[which(names(phylloI6S_pre_totalReads)== keptPhylloNames[10])] == contamInphyllo[rownames(contamInphyllo) %in% ASVsToRemovePhyllo[5], colnames(contamInphyllo) %in% keptPhylloNames[10]]

# FINISH MAKING TABLE
# What is median %abundance in each of these samples for each contaminant ASV
# 1. Make it longer!
contamInPhyllo_longer <- contamInphyllo %>% 
  rownames_to_column(var="ASV_name") %>% 
  pivot_longer(cols = starts_with("phyllo_"),
               names_to = "sampName",
               values_to = "propInSamp")

# 2. Group by ASV and then get median percentage
contamInPhyllo_longer <- contamInPhyllo_longer %>% 
  group_by(ASV_name) %>% 
  mutate(medianPercentage = median(propInSamp)*100) %>% #does not matter if median taken before or after *100
  mutate(minPercentage = min(propInSamp)*100) %>% 
  mutate(maxPercentage = max(propInSamp)*100) %>% 
  mutate(sdPercentage = sd(propInSamp)*100)

# Finally, clean this up to make it the table
colnames(contamInPhyllo_longer)
# Make ASV name a column and then below, remove Kingdom since it's not necessary
tableS2_I6S_foliar_contams_1 <- tableS2_I6S_foliar_contams_1 %>% 
  rownames_to_column(var="ASV_name")
tableS2_I6S_phyllo_contams_2 <- left_join(tableS2_I6S_foliar_contams_1[,c(1,3:ncol(tableS2_I6S_foliar_contams_1))] ,contamInPhyllo_longer[,c(1,4,7)], by = "ASV_name")
tableS2_I6S_phyllo_contams_final <- tableS2_I6S_phyllo_contams_2 %>% 
  distinct() #remove duplicate rows which were created by pivoting longer above
# write the .csv to finish this
# Written September 15, 2025
# write.csv(tableS2_I6S_phyllo_contams_final, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/manuscript/tableS2_I6S_phyllo_contams_final_Sept15.csv")

