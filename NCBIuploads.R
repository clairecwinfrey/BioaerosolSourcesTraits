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
# d. Andropogon unknown needs to be Andropogon sp.
phylloMeta1$Species[which(phylloMeta1$Genus == "Andropogon")] <- rep("sp.", time= length(which(phylloMeta1$Genus == "Andropogon")))

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

######################
# CONTROLS
# ITS
allITS_20hrs.ps <- readRDS(file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allITS_20hrs_phyloseq.rds")
ITS_ctrls.ps <- subset_samples(physeq= allITS_20hrs.ps, isControl == "control")
ITSctrlMeta <- as.data.frame(as.matrix(sample_data(ITS_ctrls.ps)))
ITSctrlMeta
colnames(ITSctrlMeta)
# Get a smaller dataframe
ITSctrlMeta <- ITSctrlMeta[,colnames(ITSctrlMeta) %in% c("EU", "DateSetOut", "sampleType", "SampleDescription", "DateDNAextracted",
                             "Sample.ID", "SampleLongerID", "isControl" )]
# View(ITSctrlMeta)
# Drop the soil-associated controls (since these were uploaded with soil data)
ITSctrlMeta <- ITSctrlMeta[-which(ITSctrlMeta$sampleType %in% c("soil_ExtBlank", "soilExtControlWater")),]
unique(ITSctrlMeta$sampleType) # [1] "fieldControl"           "washBufferControl"      "kitome"                 "PCRcontrol"            
# [5] "phyllo_NegFieldControl" "phyllo_ExtBlank"        "PCR_NegControl"  
# All of these are good! PCR control is actually a PCR negative control, but for air-associated samples

# Continue to clean up
# a. Make rownames into column
ITSctrlMeta <- ITSctrlMeta %>% 
  rownames_to_column("var" = "sampleName")
ITSctrlMeta$sampleName #all of these look fine!
# b. Make all of the field controls descriptions formatted the same way
ITSctrlMeta$SampleDescription[which(ITSctrlMeta$sampleName == "air_ITS_101")] <- "field control 7/8/22"
ITSctrlMeta$SampleDescription[which(ITSctrlMeta$sampleType %in% "fieldControl")] #almost looks good!
# Remove open parentheses
ITSctrlMeta$SampleDescription[which(ITSctrlMeta$sampleType %in% "fieldControl")] <- gsub(x= ITSctrlMeta$SampleDescription[which(ITSctrlMeta$sampleType %in% "fieldControl")],
                                                                                         pattern = "\\(", replacement = "")
# Remove close parentheses
ITSctrlMeta$SampleDescription[which(ITSctrlMeta$sampleType %in% "fieldControl")] <- gsub(x= ITSctrlMeta$SampleDescription[which(ITSctrlMeta$sampleType %in% "fieldControl")],
                                                                                         pattern = "\\)", replacement = "")
# Get rid of "DateSetOut" and "DateDNAextracted" 
ITSctrlMeta <- ITSctrlMeta[,-which(colnames(ITSctrlMeta) %in% c("DateSetOut", "DateDNAextracted", "Sample.ID", "isControl"))]
# View(ITSctrlMeta)

# Find whichever rows have something in SampleLongerID and put this in the SampleDescription (which are all blank for these controls)
ITSsampleLongIDindex <- which(!is.na(ITSctrlMeta$SampleLongerID))
ITSctrlMeta$SampleLongerID[ITSsampleLongIDindex] #as expected
ITSctrlMeta$SampleDescription[ITSsampleLongIDindex] #shows nothing there, as expected
ITSctrlMeta$SampleDescription[ITSsampleLongIDindex] <- ITSctrlMeta$SampleLongerID[ITSsampleLongIDindex]
# View(ITSctrlMeta) #confirms this looks good so drop that column
ITSctrlMeta <- ITSctrlMeta[,-which(colnames(ITSctrlMeta) %in% c("SampleLongerID"))]

# Add in description for air_ITS_PCR_NTC_2, currently the only one without a description
ITSctrlMeta$SampleDescription[which(ITSctrlMeta$sampleName == "air_ITS_PCR_NTC_2")] <- "air PCR negative control 2"

# Make this described like the others to be "blank" instead of control
ITSctrlMeta$SampleDescription[which(ITSctrlMeta$sampleName == "air_ITS_102")]
ITSctrlMeta$SampleDescription[which(ITSctrlMeta$sampleName == "air_ITS_102")] <- "wash blank #7"

# Make those without EU info have "not applicable" for EU
ITSctrlMeta$EU[which(is.na(ITSctrlMeta$EU))] <- "not applicable"

# To make the terminology more like the paper, make some changes
ITSctrlMeta$sampleType[which(ITSctrlMeta$sampleType == "kitome")] <- "air_ExtBlank" #change "kitome" to "extraction blank"
ITSctrlMeta$sampleType[which(ITSctrlMeta$sampleType == "PCRcontrol")] <- "air_PCR_NegControl"
ITSctrlMeta$sampleType[which(ITSctrlMeta$sampleType == "PCR_NegControl")] <- "phyllo_PCR_NegControl"
ITSctrlMeta$sampleType[which(ITSctrlMeta$sampleType == "fieldControl")] <- "air_fieldControl"
# View(ITSctrlMeta)

# Make all kitome in description "air DNA extraction blank"
ITSctrlMeta$SampleDescription <- gsub(x=ITSctrlMeta$SampleDescription, pattern = "kitome", 
                                      replacement= "air DNA extraction blank")

# Make all kitome in description "air DNA extraction blank"
ITSctrlMeta$SampleDescription <- gsub(x=ITSctrlMeta$SampleDescription, pattern = "NegPhylloCtrl", 
                                      replacement= "phylloNegFieldCtrl")

# Make all ExtBlank_ in description "phyllo DNA extraction blank"
ITSctrlMeta$SampleDescription <- gsub(x=ITSctrlMeta$SampleDescription, pattern = "ExtBlank_", 
                                      replacement= "phyllo DNA extraction blank")

# Remove all # in sampleDescription
ITSctrlMeta$SampleDescription <- gsub(x=ITSctrlMeta$SampleDescription, pattern = "\\#", 
                                      replacement= "")

# Remove all blank space in sampleDescription, make "_"
ITSctrlMeta$SampleDescription <- gsub(x=ITSctrlMeta$SampleDescription, pattern = " ", 
                                      replacement= "_")

# Make any double underscores one
ITSctrlMeta$SampleDescription <- gsub(x=ITSctrlMeta$SampleDescription, pattern = "__", 
                                      replacement= "_")
# View(ITSctrlMeta)

# Export
# Dec 17, 2025
# write.csv(ITSctrlMeta, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/ITSctrlMeta.csv")

##########################
# I6S
allI6S_20hrs.ps <- readRDS(file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allI6S_20hrs_phyloseq.rds")
I6S_ctrls.ps <- subset_samples(physeq= allI6S_20hrs.ps, isControl == "control")
I6SctrlMeta <- as.data.frame(as.matrix(sample_data(I6S_ctrls.ps)))
I6SctrlMeta
colnames(I6SctrlMeta)
# Get a smaller dataframe
I6SctrlMeta <- I6SctrlMeta[,colnames(I6SctrlMeta) %in% c("EU", "DateSetOut", "sampleType", "SampleDescription", "DateDNAextracted",
                                                         "Sample.ID", "SampleLongerID", "isControl" )]
# View(I6SctrlMeta)
# Drop the soil-associated controls (since these were uploaded with soil data)
I6SctrlMeta <- I6SctrlMeta[-which(I6SctrlMeta$sampleType %in% c("soil_ExtBlank", "soilExtControlWater")),]
unique(I6SctrlMeta$sampleType) # [1] "fieldControl"           "washBufferControl"      "kitome"                 "PCRcontrol"            
# [5] "phyllo_NegFieldControl" "phyllo_ExtBlank"  
# All of these are good! PCR control is actually a PCR negative control, but for air-associated samples

# Continue to clean up
# a. Make rownames into column
I6SctrlMeta <- I6SctrlMeta %>% 
  rownames_to_column("var" = "sampleName")
I6SctrlMeta$sampleName #all of these look fine!
# b. Make all of the field controls descriptions formatted the same way
I6SctrlMeta$SampleDescription[which(I6SctrlMeta$sampleName == "air_16S_101")] <- "field control 7/8/22"
I6SctrlMeta$SampleDescription[which(I6SctrlMeta$sampleType %in% "fieldControl")] #almost looks good!
# Remove open parentheses
I6SctrlMeta$SampleDescription[which(I6SctrlMeta$sampleType %in% "fieldControl")] <- gsub(x= I6SctrlMeta$SampleDescription[which(I6SctrlMeta$sampleType %in% "fieldControl")],
                                                                                         pattern = "\\(", replacement = "")
# Remove close parentheses
I6SctrlMeta$SampleDescription[which(I6SctrlMeta$sampleType %in% "fieldControl")] <- gsub(x= I6SctrlMeta$SampleDescription[which(I6SctrlMeta$sampleType %in% "fieldControl")],
                                                                                         pattern = "\\)", replacement = "")
# Get rid of "DateSetOut" and "DateDNAextracted" 
I6SctrlMeta <- I6SctrlMeta[,-which(colnames(I6SctrlMeta) %in% c("DateSetOut", "DateDNAextracted", "Sample.ID", "isControl"))]
# View(I6SctrlMeta)

# Find whichever rows have something in SampleLongerID and put this in the SampleDescription (which are all blank for these controls)
I6SsampleLongIDindex <- which(!is.na(I6SctrlMeta$SampleLongerID))
I6SctrlMeta$SampleLongerID[I6SsampleLongIDindex] #as expected
I6SctrlMeta$SampleDescription[I6SsampleLongIDindex] #shows nothing there, as expected
I6SctrlMeta$SampleDescription[I6SsampleLongIDindex] <- I6SctrlMeta$SampleLongerID[I6SsampleLongIDindex]
# View(I6SctrlMeta) #confirms this looks good so drop that column
I6SctrlMeta <- I6SctrlMeta[,-which(colnames(I6SctrlMeta) %in% c("SampleLongerID"))]

# Add in description for air_16S_PCR_NTC_1_H12 and air_16S_PCR_NTC_2, currently the only ones without a description
I6SctrlMeta$SampleDescription[which(I6SctrlMeta$sampleName == "air_16S_PCR_NTC_1_H12")] <- "air PCR negative control 1"
I6SctrlMeta$SampleDescription[which(I6SctrlMeta$sampleName == "air_16S_PCR_NTC_2")] <- "air PCR negative control 2"

# Make this described like the others to be "blank" instead of control
I6SctrlMeta$SampleDescription[which(I6SctrlMeta$sampleName == "air_16S_102")]
I6SctrlMeta$SampleDescription[which(I6SctrlMeta$sampleName == "air_16S_102")] <- "wash blank #7"

# Make those without EU info have "not applicable" for EU
I6SctrlMeta$EU[which(is.na(I6SctrlMeta$EU))] <- "not applicable"

# To make the terminology more like the paper, make some changes
I6SctrlMeta$sampleType[which(I6SctrlMeta$sampleType == "kitome")] <- "air_ExtBlank" #change "kitome" to "extraction blank"
I6SctrlMeta$sampleType[which(I6SctrlMeta$sampleType == "fieldControl")] <- "air_fieldControl"
# View(I6SctrlMeta)

# Make all kitome in description "air DNA extraction blank"
I6SctrlMeta$SampleDescription <- gsub(x=I6SctrlMeta$SampleDescription, pattern = "kitome", 
                                      replacement= "air DNA extraction blank")

# Make all kitome in description "air DNA extraction blank"
I6SctrlMeta$SampleDescription <- gsub(x=I6SctrlMeta$SampleDescription, pattern = "NegPhylloCtrl", 
                                      replacement= "phylloNegFieldCtrl")

# Make all ExtBlank_ in description "phyllo DNA extraction blank"
I6SctrlMeta$SampleDescription <- gsub(x=I6SctrlMeta$SampleDescription, pattern = "ExtBlank_", 
                                      replacement= "phyllo DNA extraction blank")

# Remove all # in sampleDescription
I6SctrlMeta$SampleDescription <- gsub(x=I6SctrlMeta$SampleDescription, pattern = "\\#", 
                                      replacement= "")

# Remove all blank space in sampleDescription, make "_"
I6SctrlMeta$SampleDescription <- gsub(x=I6SctrlMeta$SampleDescription, pattern = " ", 
                                      replacement= "_")

# Make any double underscores one
I6SctrlMeta$SampleDescription <- gsub(x=I6SctrlMeta$SampleDescription, pattern = "__", 
                                      replacement= "_")
# View(I6SctrlMeta)

# Export
# Dec 17, 2025
# write.csv(I6SctrlMeta, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/I6SctrlMeta.csv")

