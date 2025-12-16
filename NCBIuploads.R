# NCBI uploads 

library("phyloseq")
library("tidyverse")
# Read in metadata made in weather_qPCRTable.R and clean it up
final_I6SITS_table <- read.csv(file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/supplement_formatted/final_I6SITS_table_Nov21_2025.csv")
head(final_I6SITS_table)
final_I6SITS_table$X <- NULL

# For the NCBI SRA upload get sample names with ITS or 16S (for relevant samples) appended in front:
## ITS ##
ITSonlyMetaFinal <- final_I6SITS_table
colnames(ITSonlyMetaFinal)
ITSsampleNames <- paste("ITS", final_I6SITS_table$sampleName, sep = "_")
ITSonlyMetaFinal <- cbind(ITSonlyMetaFinal, ITSsampleNames)
head(ITSonlyMetaFinal)
# Add in column for "*env_broad_scale" column based on "HabitatAir"
ITSonlyMetaFinal <- ITSonlyMetaFinal %>% 
  mutate(
    env_broad_scale = case_when(
      HabitatAir %in% "savanna" ~ "ENVO:01000187",
      HabitatAir %in% "forest" ~ "ENVO:01000222|ENVO:00000117"
    )
  )
head(ITSonlyMetaFinal)

# Dec. 14, 2025
# write.csv(x = ITSonlyMetaFinal, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/ITSonlyMetaFinal.csv")

## 16S ##
# Check that no samples that failed for qPCR are actually used samples
length(which(is.na(final_I6SITS_table$I6Scopies)== FALSE)) #84 had qPCR data so we are good!
# Remove those samples without 16S data, i.e., keep only those that do NOT have NA in I6Scopies column
I6SonlyMetaFinal <- final_I6SITS_table[which(is.na(final_I6SITS_table$I6Scopies)== FALSE), ]
dim(I6SonlyMetaFinal) #84
I6SsampleNames <- paste("16S", I6SonlyMetaFinal$sampleName, sep = "_")
I6SonlyMetaFinal <- cbind(I6SonlyMetaFinal, I6SsampleNames)
# Add in column for "*env_broad_scale" column based on "HabitatAir"
I6SonlyMetaFinal <- I6SonlyMetaFinal %>% 
  mutate(
    env_broad_scale = case_when(
      HabitatAir %in% "savanna" ~ "ENVO:01000187",
      HabitatAir %in% "forest" ~ "ENVO:01000222|ENVO:00000117"
    )
  )
head(I6SonlyMetaFinal)
# View(I6SonlyMetaFinal)

# Dec. 14, 2025
# write.csv(x = I6SonlyMetaFinal, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/I6SonlyMetaFinal.csv")

################################
# PHYLLOSPHERE -- see supplementalFigures_Oct16_2025 for the track changes versions to make sure that this metadata is up to date and w/o typos
# in the metadata file
# 1. Load in data file (made May 16-17, 2025, see notes there) and make sure it's correct
# Even though it says 16S, this is mostly the same for both samples
phylloMeta1 <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/phyllo16SMetaDataFinal.csv")
head(phylloMeta1)
sort(phylloMeta1$GenusSpeciesCode)
# i. Check that "DICVIL" is "DICVOL" as it should be
which(phylloMeta1$GenusSpeciesCode == "DICVOL") #great this was already updated!
phylloMeta1$SampleLongerID[which(phylloMeta1$GenusSpeciesCode == "DICVOL")] #and this was fixed too, excellent!
# ii. Make sure that ASTPAT is Sericocarpus	asteroides to match convention of the rest we used.
# Great these both looked as expected!
phylloMeta1$Genus[which(phylloMeta1$GenusSpeciesCode == "ASTPAT")] #"Sericocarpus"
phylloMeta1$Species[which(phylloMeta1$GenusSpeciesCode == "ASTPAT")] #"asteroides"
# iii. There are some typos in the species names
# a. Vitis rotundafolia should be "rotundifolia"
phylloMeta1$Species[which(phylloMeta1$Species == "rotundafolia")] <- rep("rotundifolia", times = length(which(phylloMeta1$Species == "rotundafolia")))
phylloMeta1$Species[which(phylloMeta1$Species == "rotundifolia")] #looks good!
# b. Centrosema virginanum should be Centrosema virginianum
phylloMeta1$Species[which(phylloMeta1$Species == "virginanum")] <- rep("virginianum", times = length(which(phylloMeta1$Species == "virginanum")))
phylloMeta1$Species[which(phylloMeta1$Species == "virginianum")] #looks good!
# c. While is the correct species, I have to use an older version so that NCBI doesn't get mad:
# See https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=3042424 and https://fsus.ncbg.unc.edu/show-taxon-detail.php?taxonid=1955
# So Dichanthelium should be Panicum for villosissimum
phylloMeta1[which(phylloMeta1$Genus == "Dichanthelium"),] #there is species = villosissimum and aciculare, but
# only villosissimum gave the error
phylloMeta1$Genus[which(phylloMeta1$Species == "villosissimum")] <- "Panicum"
phylloMeta1$Genus[which(phylloMeta1$Species == "villosissimum")] #Great it is now Panicum!

# 2. Add in other columns needed for the NCBI metadata that are different from row to row,
# Make a genus and species column for "host" as the plant it was sampled from
# First get rid of any weird spaces in the genus and species
phylloMeta1$Genus <- gsub(x = phylloMeta1$Genus, pattern = " ", replacement ="")
phylloMeta1$Species <- gsub(x = phylloMeta1$Species, pattern = " ", replacement ="")
# Add in a new host column by binding these together
phylloMeta1$Host <- paste(phylloMeta1$Genus, phylloMeta1$Species, sep = " ")

# 3. Make separate for ITS and append names
phylloMeta_ITS <- phylloMeta1
# Add in a sampleName column that will be sample_name for NCBI datasheet
phylloMeta_ITS$sampleName <- paste("ITS", phylloMeta_ITS$SampleLongerID, sep = "_")
# Get rid of "EU" because it's too long
phylloMeta_ITS$sampleName <- gsub(x=phylloMeta_ITS$sampleName, pattern = "EU_", replacement = "")
phylloMeta_ITS$sampleName # looks good!
# December 16, 2025
# write.csv(phylloMeta_ITS, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/phylloMeta_ITS.csv")

# 4. Make separate for 16S. Remove sample that didn't work and append 16S to new sample names
phyllo16Smeta1 <- phylloMeta_ITS
dim(phyllo16Smeta1) #This is 59 long, but EUPCOM sampled in EU 52 returned too few16S rRNA reads to be considered in bacterial analyses
head(phyllo16Smeta1)
# Remove EUPCOM
phyllo16Smeta1$EU[which(phyllo16Smeta1$GenusSpeciesCode == "EUPCOM")] #only 1 and it's EU52. So drop this one!
phyllo16Smeta2 <- phyllo16Smeta1[-which(phyllo16Smeta1$GenusSpeciesCode == "EUPCOM"),]
dim(phyllo16Smeta2) #58 long now.
phyllo16Smeta2$EU[which(phyllo16Smeta2$GenusSpeciesCode == "EUPCOM")] #no longer there!
# Add in sampleName by replacing ITS with 16S
phyllo16Smeta2$sampleName <- gsub(x= phyllo16Smeta2$sampleName, pattern = "ITS_", replacement = "16S_")
phylloMeta_16S <- phyllo16Smeta2
# December 16, 2025
# write.csv(phylloMeta_16S, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/phylloMeta_16S.csv")

