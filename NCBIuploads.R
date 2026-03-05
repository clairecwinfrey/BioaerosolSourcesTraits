# NCBI uploads 
# This script was created in order to get metadata in shape for the substantial
# amount of files required for NCBI BioProject, BioSamples, and finally the 
# Sequence Read Archive.
# The link for these on NCBI's website:

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

###################################################
# SRA METADATA FILES
# (For filling out I6S_SRA_metadata_acc.xlsx and ITS_SRA_metadata_acc.xlsx,
# both of which are here: ~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets)

# LOAD BioSample attributes files for each sample type:
list.files("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/downloadedBioSampleAttributes")
ITSairBioSamps <-read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/downloadedBioSampleAttributes/ITS_bioAerosolsBioSampleAttributes.csv")
ITSphylloBioSamps <-read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/downloadedBioSampleAttributes/ITS_PhylloBioSampleAttributes.csv")
ITSCtrlsBioSamps <-read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/downloadedBioSampleAttributes/ITS_ControlsBioSampleAttributes.csv")
I6SairBioSamps <-read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/downloadedBioSampleAttributes/I6S_bioAerosolsBioSampleAttributes.csv")
I6SphylloBioSamps <-read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/downloadedBioSampleAttributes/I6S_PhylloBioSampleAttributes.csv")
I6SCtrlsBioSamps <-read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/NCBIuploadsAndDataSheets/downloadedBioSampleAttributes/ControlsBioSampleAttributes.csv")

# Read in fastq.gz names. Note that these include PCR controls too, except in
# the case of phyllo 16S since these did not PCR controls did not make it through sequencing.
# Control file names will have to be dealt with separately.
bioaersolITSfileNames <- read.delim("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/ITSdata/airITS_demultiplexed/bioaerosol_ITS_fastq_files.txt", header = FALSE) #ITS bioaerosols
bioaersol16SfileNames <- read.delim("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/I6Sdata/air16S_demultiplexed/bioaerosol_16S_fastq_files.txt", header = FALSE) #16S bioaerosols
phylloITSfileNames <- read.delim("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/ITSdata/phylloITS_demultiplexed/phyllo_ITS_fastq_files.txt", header = FALSE) #ITS phyllo
phyllo16SfileNames <- read.delim("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/I6Sdata/phyllo16S_demultiplexed/phyllo_16S_fastq_files.txt", header = FALSE) #16S phyllo

########### I. ITS bioaerosols ###########
head(ITSairBioSamps)
sort(colnames(ITSairBioSamps))

# Make library ID column
ITSairBioSamps$library_ID <- paste0(ITSairBioSamps$sample_name, "_", "bioaerosol")
# Make "Title" column
ITSairBioSamps$title <- paste0("ITS region rRNA of bioaerosol sample:", ITSairBioSamps$sample_name)

# Get sample name and whether its R1 or R2 
bioaersolITSfileNames_clean <- bioaersolITSfileNames %>%
  mutate(
    sample_name = str_extract(V1, "ITS_\\d+"),
    read = str_extract(V1, "^R[12]")
  ) 
# Before pivoting wider, I need to:
# i. drop Unassigned reads
bioaersolITSfileNames_clean <- bioaersolITSfileNames_clean[-which(grepl(x=bioaersolITSfileNames_clean$V1, pattern="Unassigned_reads")),]
unique(grepl(x=bioaersolITSfileNames_clean$V1, pattern = "Unassigned")) #removed!
# ii. assign sample names to those with NAs
bioaersolITSfileNames_clean[which(is.na(bioaersolITSfileNames_clean$sample_name)),]
bioaersolITSfileNames_clean$sample_name[which(is.na(bioaersolITSfileNames_clean$sample_name))] <-
  c("air_ITS_PCR_NTC_1", "air_ITS_PCR_NTC_2_H12", "air_ITS_PCR_NTC_2", "air_ITS_PCR_NTC_1", 
    "air_ITS_PCR_NTC_2_H12", "air_ITS_PCR_NTC_2")

# Now pivot wider to make this look like what we want!
bioaersolITSfileNames_clean <- bioaersolITSfileNames_clean %>%
select(sample_name, read, V1) %>%
pivot_wider(
  names_from = read,
  values_from = V1
)
#View(bioaersolITSfileNames_clean )

# Finally merge with ITSairBioSamps
colnames(ITSairBioSamps)
colsOfInterest <- c("accession", "sample_name", "title", "bioproject_accession", "experimentalunit", "library_ID", "title")
ITSairBioSamps[,colnames(ITSairBioSamps) %in% colsOfInterest]
ITSairBioSamps_df <- left_join(ITSairBioSamps[,colnames(ITSairBioSamps) %in% colsOfInterest], bioaersolITSfileNames_clean, by = "sample_name")
# View(ITSairBioSamps_df)

########### II. 16S bioaerosols ###########
head(I6SairBioSamps)
sort(colnames(I6SairBioSamps))

# Make library ID column
I6SairBioSamps$library_ID <- paste0(I6SairBioSamps$sample_name, "_", "bioaerosol")
# Make "Title" column
I6SairBioSamps$title <- paste0("16S rRNA gene of bioaerosol sample:", I6SairBioSamps$sample_name)

# Get sample name and whether its R1 or R2 
bioaersol16SfileNames_clean <- bioaersol16SfileNames %>%
  mutate(
    sample_name = str_extract(V1, "16S_\\d+"),
    read = str_extract(V1, "^R[12]")
  ) 
# Before pivoting wider, I need to:
# i. drop Unassigned reads
bioaersol16SfileNames_clean <- bioaersol16SfileNames_clean[-which(grepl(x=bioaersol16SfileNames_clean$V1, pattern="Unassigned_reads")),]
unique(grepl(x=bioaersol16SfileNames_clean$V1, pattern = "Unassigned")) #removed!
# ii. assign sample names to those with NAs
bioaersol16SfileNames_clean[which(is.na(bioaersol16SfileNames_clean$sample_name)),]
bioaersol16SfileNames_clean$sample_name[which(is.na(bioaersol16SfileNames_clean$sample_name))] <-
  c("air_16S_PCR_NTC_1_H12", "air_16S_PCR_NTC_1", "air_16S_PCR_NTC_2_H11", "R1_air_16S_PCR_NTC_2", 
    "air_16S_PCR_NTC_1_H12", "air_16S_PCR_NTC_1", "air_16S_PCR_NTC_2_H11", "air_16S_PCR_NTC_2")

# Now pivot wider to make this look like what we want!
bioaersol16SfileNames_clean <- bioaersol16SfileNames_clean %>%
  select(sample_name, read, V1) %>%
  pivot_wider(
    names_from = read,
    values_from = V1
  )
#View(bioaersol16SfileNames_clean)

# Finally merge with I6SairBioSamps, by left joining so that only those in the metadata online get put in there!
colnames(I6SairBioSamps)
colsOfInterest <- c("accession", "sample_name", "title", "bioproject_accession", "experimentalunit", "library_ID", "title")
I6SairBioSamps[,colnames(I6SairBioSamps) %in% colsOfInterest]
I6SairBioSamps_df <- left_join(I6SairBioSamps[,colnames(I6SairBioSamps) %in% colsOfInterest], bioaersol16SfileNames_clean, by = "sample_name")
# View(I6SairBioSamps_df)

########### III. ITS leaf surfaces ###########
head(ITSphylloBioSamps) #the sample names here (so from NCBI's perspective) are
# based on EU and the type of plant
head(phylloITSfileNames) #this is based on the sample names in the DNA plate
# So, I'll have to bring in a different metadata file to find the link between
# the sequence file name and the current sample name I have in NCBI. 

# Made for ITS but the keys are the same for 16S and ITS (and ITS had one more sample than 16S)
all_ITS_Metadata <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/aeroAllMetadat_Apr20_2023.csv", row.names=1) #made in ITS_metadataAndMappingFileSetUp.R
sort(colnames(all_ITS_Metadata))
# Pull out just the relevant columns and the phyllosphere
phylloOnlyMetaDat <- all_ITS_Metadata %>% 
  select(sampleNumber, sampleType, SampleLongerID, PlantSpecies) %>% 
  filter(sampleType == "phyllosphere")
# sampleNumber goes with the phylloITSfileNames numbers 
head(phylloOnlyMetaDat)
# Rename sampleNumber in the metadata sample_name for later merging
phylloOnlyMetaDat$sampleNumber
colnames(phylloOnlyMetaDat)[which(colnames(phylloOnlyMetaDat) %in% "sampleNumber"== TRUE)] <- "sample_name"

# Make library ID column
ITSphylloBioSamps$library_ID <- paste0(ITSphylloBioSamps$sample_name, "_", "foliarSurface")
# Make "Title" column
ITSphylloBioSamps$title <- paste0("ITS region rRNA of foliar surface sample:", ITSphylloBioSamps$sample_name)

# Get sample name and whether its R1 or R2 
phylloITSfileNames_clean <- phylloITSfileNames %>%
  mutate(
    sample_name = str_extract(V1, "phyllo_ITS_\\d+"),
    read = str_extract(V1, "^R[12]")
  ) 
head(phylloITSfileNames_clean) #looks good!

# Before pivoting wider, I need to:
# i. Drop Unassigned reads
phylloITSfileNames_clean <- phylloITSfileNames_clean[-which(grepl(x=phylloITSfileNames_clean$V1, pattern="Unassigned_reads")),]
unique(grepl(x=phylloITSfileNames_clean$V1, pattern = "Unassigned")) #removed!
# ii. Assign sample names to those with NAs
phylloITSfileNames_clean[which(is.na(phylloITSfileNames_clean$sample_name)),]
phylloITSfileNames_clean$sample_name[which(is.na(phylloITSfileNames_clean$sample_name))] <-
  c("phyllo_ITS_PCR_NTC_1", "phyllo_ITS_PCR_NTC_2", "phyllo_ITS_PCR_NTC_1", "phyllo_ITS_PCR_NTC_2")
# View(phylloITSfileNames_clean)

# Now pivot wider to make this look like what we want!
phylloITSfileNames_clean <- phylloITSfileNames_clean %>%
  select(sample_name, read, V1) %>%
  pivot_wider(
    names_from = read,
    values_from = V1
  )
# View(phylloITSfileNames_clean)

# Finally merge to get metadata.
# 1. Is anything mising?
setdiff(phylloOnlyMetaDat$sample_name, phylloITSfileNames_clean$sample_name) #no samples in metadata not in file names
setdiff(phylloITSfileNames_clean$sample_name, phylloOnlyMetaDat$sample_name) #But there are samples in files not in metadata.
# # My best guess is that these are wash buffer and extraction blank controls. And I think I remember PCR NTCs failing
# [1] "phyllo_ITS_27"        "phyllo_ITS_28"        "phyllo_ITS_51"        "phyllo_ITS_52"        "phyllo_ITS_6"        
# [6] "phyllo_ITS_65"        "phyllo_ITS_66"        "phyllo_ITS_PCR_NTC_1" "phyllo_ITS_PCR_NTC_2"
# 2. Look at what we need for biosamples.
colsOfInterest <- c("accession", "sample_name", "title", "bioproject_accession", "experimentalunit", "library_ID", "title")
nrow(ITSphylloBioSamps[,colnames(ITSphylloBioSamps) %in% colsOfInterest]) #59 samples
length(phylloOnlyMetaDat$sample_name) #also 59, so I am guessing that we are fine above
# 2. Left join the phylloITSfileNames_clean and phylloOnlyMetaDat
fileNameMeta1 <- left_join(phylloITSfileNames_clean, phylloOnlyMetaDat, by= "sample_name")

# After the left_join, “sample_name” in ITSphylloBioSamps is almost SampleLongerID in the
# original metadata, except the EU in SampleLongerID should be replaced with “ITS”
head(fileNameMeta1)
fileNameMeta1$SampleLongerID
ITSphylloBioSamps$sample_name

# For merging (not left_joining!), first change the "sample_name" from fileNameMeta1 to "anotherSampName", change SampleLongerID
# colname to sample_name, and then change "EU_" in newly named "sample_name" to "ITS"
colnames(fileNameMeta1)[which(colnames(fileNameMeta1) %in% "sample_name"==TRUE)] <- "anotherSampName"
colnames(fileNameMeta1)[which(colnames(fileNameMeta1) %in% "SampleLongerID"==TRUE)] <- "sample_name"
fileNameMeta1$sample_name <- gsub(x=fileNameMeta1$sample_name, pattern = "EU_", replacement = "ITS_")

# Now, merge fileNameMeta1 with ITSphylloBioSamps[,colnames(ITSphylloBioSamps) %in% colsOfInterest]
ITSphylloBioSamps_df <- merge(ITSphylloBioSamps[,colnames(ITSphylloBioSamps) %in% colsOfInterest], fileNameMeta1, by = "sample_name")

########### IV. 16S leaf surfaces ###########
head(I6SphylloBioSamps) #the sample names here (so from NCBI's perspective) are
# based on EU and the type of plant
head(phyllo16SfileNames) #this is based on the sample names in the DNA plate
# So, I'll have to bring in a different metadata file to find the link between
# the sequence file name and the current sample name I have in NCBI. 

# Made for ITS but the keys are the same for 16S and ITS (and ITS had one more sample than 16S)
phylloOnlyMetaDat #edited above for ITS

# Make library ID column
I6SphylloBioSamps$library_ID <- paste0(I6SphylloBioSamps$sample_name, "_", "foliarSurface")
# Make "Title" column
I6SphylloBioSamps$title <- paste0("16S region rRNA of foliar surface sample:", I6SphylloBioSamps$sample_name)

# Get sample name and whether its R1 or R2 
phyllo16SfileNames_clean <- phyllo16SfileNames %>%
  mutate(
    sample_name = str_extract(V1, "phyllo_16S_\\d+"),
    read = str_extract(V1, "^R[12]")
  ) 
head(phyllo16SfileNames_clean) #looks good!

# Before pivoting wider, I need to:
# i. Drop Unassigned reads
phyllo16SfileNames_clean <- phyllo16SfileNames_clean[-which(grepl(x=phyllo16SfileNames_clean$V1, pattern="Unassigned_reads")),]
unique(grepl(x=phyllo16SfileNames_clean$V1, pattern = "Unassigned")) #removed!
# ii. Assign sample names to those with NAs
phyllo16SfileNames_clean[which(is.na(phyllo16SfileNames_clean$sample_name)),] #no NAs
# View(phyllo16SfileNames_clean)

# Now pivot wider to make this look like what we want!
phyllo16SfileNames_clean <- phyllo16SfileNames_clean %>%
  select(sample_name, read, V1) %>%
  pivot_wider(
    names_from = read,
    values_from = V1
  )
# View(phyllo16SfileNames_clean)

# Finally merge to get metadata.
# 1. Change ITS to be 16S in the metadata file
phylloOnlyMetaDat_16S <- phylloOnlyMetaDat #make a copy
phylloOnlyMetaDat_16S$sample_name <- gsub(x=phylloOnlyMetaDat_16S$sample_name , pattern = "_ITS_", replacement = "_16S_")
# 2. Is anything mising?
setdiff(phylloOnlyMetaDat_16S$sample_name, phyllo16SfileNames_clean$sample_name) #no samples in metadata not in file names
setdiff(phyllo16SfileNames_clean$sample_name, phylloOnlyMetaDat_16S$sample_name) #But there are samples in files not in metadata.
# My best guess is that these are wash buffer and extraction blank controls
# "phyllo_16S_27" "phyllo_16S_28" "phyllo_16S_51" "phyllo_16S_52" "phyllo_16S_6"  "phyllo_16S_65" "phyllo_16S_66"
# 3. Look at what we need for biosamples.
nrow(I6SphylloBioSamps[,colnames(I6SphylloBioSamps) %in% colsOfInterest]) #58 samples
length(phylloOnlyMetaDat_16S$sample_name) #59, so this includes a sample that worked for ITS but not 16S 
# (which was EUPCOM sampled in EU 52), but this will be dropped by left_joining next step
# 3. Left join the phyllo16SfileNames_clean and phylloOnlyMetaDat_16S
fileNameMeta16S_1 <- left_join(phyllo16SfileNames_clean, phylloOnlyMetaDat_16S, by= "sample_name")

# After the left_join, “sample_name” in I6SphylloBioSamps is almost SampleLongerID in the
# original metadata, except the EU in SampleLongerID should be replaced with “ITS”
head(fileNameMeta16S_1)
fileNameMeta16S_1$SampleLongerID
I6SphylloBioSamps$sample_name

# For merging (not left_joining!), first change the "sample_name" from fileNameMeta16S_1 to "anotherSampName", change SampleLongerID
# colname to sample_name, and then change "EU_" in newly named "sample_name" to "ITS"
colnames(fileNameMeta16S_1)[which(colnames(fileNameMeta16S_1) %in% "sample_name"==TRUE)] <- "anotherSampName"
colnames(fileNameMeta16S_1)[which(colnames(fileNameMeta16S_1) %in% "SampleLongerID"==TRUE)] <- "sample_name"
fileNameMeta16S_1$sample_name <- gsub(x=fileNameMeta16S_1$sample_name, pattern = "EU_", replacement = "16S_")

# Now, merge fileNameMeta16S_1 with I6SphylloBioSamps[,colnames(I6SphylloBioSamps) %in% colsOfInterest]
I6SphylloBioSamps_df <- merge(I6SphylloBioSamps[,colnames(I6SphylloBioSamps) %in% colsOfInterest], fileNameMeta16S_1, by = "sample_name")
# View(I6SphylloBioSamps_df) #Looks good!

########### V. ITS CONTROLS! ########### 
# View(ITSCtrlsBioSamps)
# These samples are spread across bioaersolITSfileNames and phylloITSfileNames,
# and it's also possible that not all of them are even in these, but hopefully 
# they are!
ITSCtrlsBioSamps$sample_name 
head(bioaersolITSfileNames)
head(phylloITSfileNames)

# Make library ID column
ITSCtrlsBioSamps$library_ID <- paste0(ITSCtrlsBioSamps$sample_name, "_", "blank") #may want to add in bioaerosol or phyllo control later
# Make "Title" column
ITSCtrlsBioSamps$title <- paste0("ITS region rRNA of blank:", ITSCtrlsBioSamps$sample_name)

# Look at previously made dataframes that contain these blank sequence file locations
head(bioaersolITSfileNames_clean)
head(phylloITSfileNames_clean)
intersect(bioaersolITSfileNames_clean$sample_name, phylloITSfileNames_clean$sample_name) #none shared so can rbind
# Rbind these for left_joining later
allITSfileNames <- rbind(bioaersolITSfileNames_clean, phylloITSfileNames_clean)
#View(allITSfileNames)

# Merge with ITSCtrlsBioSamps
sort(colnames(ITSCtrlsBioSamps))
ITSCtrlsBioSamps[,colnames(ITSCtrlsBioSamps) %in% colsOfInterest]
ITSCtrlsBioSamps_df <- left_join(ITSCtrlsBioSamps[,colnames(ITSCtrlsBioSamps) %in% colsOfInterest], allITSfileNames, by = "sample_name")
# View(ITSCtrlsBioSamps_df)

# Get indices and then add in details about what kind of control it is
# 1. Those that have EU info and are NOT "not applicable" are field blanks 
fieldBlanksITS_index <- which(ITSCtrlsBioSamps_df$experimentalunit != "not applicable")
# 2. Get only the phyllo sample index
phylloITS_index <- which(grepl(x=ITSCtrlsBioSamps_df$sample_name, pattern= "phyllo")==TRUE)
# 3. Get only the bioaerosol sample index
airITS_index <- which(grepl(x=ITSCtrlsBioSamps_df$sample_name, pattern= "air")==TRUE)
# 4. Get the indices for the PCR controls
PCRITS_index <- which(grepl(x=ITSCtrlsBioSamps_df$sample_name, pattern= "PCR_NTC")==TRUE)
# 5. Intersection of those with EU info and air are air field controls
airFieldCtrlIndex <- intersect(fieldBlanksITS_index, airITS_index); length(airFieldCtrlIndex) #16 
ITSCtrlsBioSamps_df[airFieldCtrlIndex,]
ITSCtrlsBioSamps_df$title[airFieldCtrlIndex] <- gsub(x=ITSCtrlsBioSamps_df$title[airFieldCtrlIndex], pattern = "of blank", replacement = "of bioaerosol field blank")
ITSCtrlsBioSamps_df$title[airFieldCtrlIndex]
# 6. Intersection of PCR and phyllo are phyllo PCR controls
phylloPCRCtrlIndex <- intersect(PCRITS_index, phylloITS_index) #2
ITSCtrlsBioSamps_df[phylloPCRCtrlIndex ,]
ITSCtrlsBioSamps_df$title[phylloPCRCtrlIndex] <- gsub(x=ITSCtrlsBioSamps_df$title[phylloPCRCtrlIndex], pattern = "of blank", replacement = "of foliar surface PCR no-template control")
ITSCtrlsBioSamps_df$title[phylloPCRCtrlIndex]
# 7. Intersection of PCR and air are air PCR controls
airPCRCtrlIndex <- intersect(PCRITS_index, airITS_index) #just one
ITSCtrlsBioSamps_df[airPCRCtrlIndex ,]
ITSCtrlsBioSamps_df$title[airPCRCtrlIndex] <- gsub(x=ITSCtrlsBioSamps_df$title[airPCRCtrlIndex], pattern = "of blank", replacement = "of bioaerosol PCR no-template control")
ITSCtrlsBioSamps_df$title[airPCRCtrlIndex]

# For the extraction blanks and wash buffer blanks, will have to look at longer metadata:
# Pull out just the relevant columns for controls
sort(colnames(all_ITS_Metadata))
unique(all_ITS_Metadata$sampleType)
ctrlsOnlyMetaDat <- all_ITS_Metadata %>% 
  select(sampleNumber, sampleType, SampleLongerID, SampleDescription) %>% 
  filter(sampleType %in% c("kitome", "PCRcontrol", "washBufferControl", "fieldControl", "phyllo_ExtBlank",
                           "PCR_NegControl" , "phyllo_NegFieldControl"))
# View(ctrlsOnlyMetaDat)
# 8. Extraction blanks for biaoerosol samples, which I called "DNA extraction blanks" in the paper
# Note that there are 10 here
ctrlsOnlyMetaDat[which(ctrlsOnlyMetaDat$sampleType == "kitome"),]
airKitomeNames <- ctrlsOnlyMetaDat$sampleNumber[which(ctrlsOnlyMetaDat$sampleType == "kitome")]
airKitomeIndex <- which(ITSCtrlsBioSamps_df$sample_name %in% airKitomeNames==TRUE) #these are the kitome samples
ITSCtrlsBioSamps_df[airKitomeIndex,]
ITSCtrlsBioSamps_df$title[airKitomeIndex] <- gsub(x=ITSCtrlsBioSamps_df$title[airKitomeIndex], pattern = "of blank", replacement = "of bioaerosol DNA extraction blank")

# 9. Wash buffer blanks for biaoerosol samples, which I called "wash buffer blanks" in the paper
# Note that there are 10 here
unique(ctrlsOnlyMetaDat$sampleType)
ctrlsOnlyMetaDat[which(ctrlsOnlyMetaDat$sampleType == "washBufferControl"),]
airWBBNames_ITS <- ctrlsOnlyMetaDat$sampleNumber[which(ctrlsOnlyMetaDat$sampleType == "washBufferControl")]
airwashBufferControlIndex_ITS <- which(ITSCtrlsBioSamps_df$sample_name %in% airWBBNames_ITS==TRUE) #these are the washBufferControl samples
ITSCtrlsBioSamps_df[airwashBufferControlIndex_ITS,]
ITSCtrlsBioSamps_df$title[airwashBufferControlIndex_ITS] <- gsub(x=ITSCtrlsBioSamps_df$title[airwashBufferControlIndex_ITS], pattern = "of blank", replacement = "of bioaerosol wash buffer blank")

# 10. Last remaining phyllosphere samples?? Which are these?

