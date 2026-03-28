# ITS full EDA rarefied - PART 2

# Script description: This script follows fullITS_EDA_part1_Sept.R (which cleans data and gets rarefied dataset).

# *** MAY NEED TO UPDATE THE REST OF THIS:
# Specifically, First: got convenient subsets of the data (samples only, controls only, air only, phyllosphere only, 
# soil only, singletons and doubletons removed from air, etc. ). Second: NMDS ordinations and visualizations of 
# different subsets of the data. Third, PERMANOVAs showuing that samples differ based on sample Type. Fourth,
# (preliminary, i.e. not the final version used in MS) permanovas to look for difference between forest and savanna (in several 
# different ways).Fourth: plots and calculations of top fungal classes and orders. 

# Original version of this script was called: full16S_EDArarefied_Part2.R. Originally run on server
# Differences between that script and this one include new treatment of phyllopshere (removing contaminants) and no more in-depth 
# dispersion or wind effects analyses.

#######################################
#             SCRIPT SET UP
#######################################
# Read in libraries
library("phyloseq") #‘1.32.0’
library("tidyverse") #‘1.3.1’
library("vegan") #vegan 2.5-7
library("gridExtra")    #allows you to make multiple plots on the same page with ggplot
library("Polychrome") #pretty colors for data viz

load(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/allITS_rarefied_phyloseq") #Made in fullITS_EDA_part1_Sept.R
allITS_rarefied.ps #13,144 taxa and 328 samples (same as before)

length(grep(x= colnames(otu_table(allITS_rarefied.ps)), pattern = "air")) #112 air samples/blanks, as expected

# Load unrarefied data (called allITS_20hrs.ps)
allITS_20hrs.ps <- readRDS(file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allITS_20hrs_phyloseq.rds")

######
# NON-RAREFIED: 
####
# 1. REMOVE SAMPLES WITH FEWER THAN 8,500 READS (TO MATCH RAREFIED VERSION)
# Identify the samples to keep
samplesToKeep_ITS <- sample_sums(allITS_20hrs.ps) >= 8500 #sample_sums() returns the total counts for each sample
# Prune out these samples
ITSall_8.5K.ps <- prune_samples(samplesToKeep_ITS, allITS_20hrs.ps)
min(sample_sums(ITSall_8.5K.ps)) #8846 is new minimum
# Remove taxa without reads
ITSall_8.5K.ps  <- prune_taxa(taxa_sums(ITSall_8.5K.ps) > 0, ITSall_8.5K.ps )
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 13778 taxa and 328 samples ]
# sample_data() Sample Data:       [ 328 samples by 40 sample variables ]
# tax_table()   Taxonomy Table:    [ 13778 taxa by 7 taxonomic ranks ]

# Saved September 30, 2025
# saveRDS(ITSall_8.5K.ps, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITSall_8.5K_phyloseq.rds")

# 2. AIR SUBSETTED DATASET: REMOVE BIOAEROSOL SINGELTONS AND DOUBLETONS 
airITS_8.5K.ps <- subset_samples(ITSall_8.5K.ps, sampleType == "air")
sample_data(airITS_8.5K.ps)$sampleType
# Remove taxa without reads
airITS_8.5K.ps <- prune_taxa(taxa_sums(airITS_8.5K.ps) > 0, airITS_8.5K.ps)
# Find singletons and doubletons (do not do < because we want to keep those that don't appear because they are in soil or phyllosphere!)
ITS_NotRareSingletons <- which(rowSums(otu_table(airITS_8.5K.ps)) == 1) 
ITS_NotRareDoubletons <- which(rowSums(otu_table(airITS_8.5K.ps)) == 2) 
ITS_notRareSingDoubs <- c(names(ITS_NotRareSingletons), names(ITS_NotRareDoubletons))
length(ITS_notRareSingDoubs) #106 taxa
unique(rowSums(otu_table(airITS_8.5K.ps)[rownames(otu_table(airITS_8.5K.ps)) %in% ITS_notRareSingDoubs,])) #confirms that all are one or two!
# Prune out these taxa
goodNotRareAirASVsITS <- setdiff(rownames(otu_table(airITS_8.5K.ps)),ITS_notRareSingDoubs)
airITS_8.5K.ps <- prune_taxa(x= airITS_8.5K.ps, taxa = goodNotRareAirASVsITS)

# Add in a column that has habitat type for air
ITS_airOnly_notR_meta <- as.data.frame(as.matrix(sample_data(airITS_8.5K.ps)))
airSavIndex <- intersect(which(ITS_airOnly_notR_meta$sampleType == "air"), which(ITS_airOnly_notR_meta$Habitat == "savanna"))
airForIndex <- intersect(which(ITS_airOnly_notR_meta$sampleType == "air"), which(ITS_airOnly_notR_meta$Habitat == "forest"))
ITS_airOnly_notR_meta$HabitatAir <- rep(NA, nrow(ITS_airOnly_notR_meta))
ITS_airOnly_notR_meta$HabitatAir[airSavIndex] <- "savanna"
ITS_airOnly_notR_meta$HabitatAir[airForIndex] <- "forest"
ITS_airOnly_notR_meta$HabitatAir[which(ITS_airOnly_notR_meta$sampleType == "soil")] <- "soil"
ITS_airOnly_notR_meta$HabitatAir[which(ITS_airOnly_notR_meta$sampleType == "phyllosphere")] <- "phyllosphere"

ITS_airOnly_notR_meta <- sample_data(ITS_airOnly_notR_meta) #make this into "sample data" compatible with phyloseq
airITS_8.5K.ps <- merge_phyloseq(airITS_8.5K.ps, sample_data(ITS_airOnly_notR_meta)) 
# otu_table()   OTU Table:         [ 5612 taxa and 110 samples ]
# sample_data() Sample Data:       [ 110 samples by 41 sample variables ]
# tax_table()   Taxonomy Table:    [ 5612 taxa by 7 taxonomic ranks ]

sort(colnames(sample_data(airITS_8.5K.ps))) #shows HabitatAir was incorporated!

# Saved September 30, 2025
# saveRDS(airITS_8.5K.ps, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airITS_8.5K_phyloseq.rds")

# 3. ALL ITS WITH SINGLETONS AND DOUBLETONS FOUND ABOVE REMOVED FOR BIOAEROSOLS
# (Here remove leftover controls too since they were investigated and found to be fine in last project (soil) or were simly
# unable to get rid of (one field control))
ITSall_8.5K_noAir <- subset_samples(ITSall_8.5K.ps, sampleType!= "air")
unique(sample_data(ITSall_8.5K_noAir)$sampleType)
ITSall_8.5K_noAir <- subset_samples(ITSall_8.5K_noAir, sampleType!= "soil_ExtBlank")
ITSall_8.5K_noAir <- subset_samples(ITSall_8.5K_noAir, sampleType!= "soilExtControlWater")
ITSall_8.5K_noAir <- subset_samples(ITSall_8.5K_noAir, sampleType!= "fieldControl")
unique(sample_data(ITSall_8.5K_noAir)$sampleType) # "phyllosphere" "soil"        
# Make a new phyloseq
ITSall_8.5Kfiltered.ps <- merge_phyloseq(ITSall_8.5K_noAir, airITS_8.5K.ps)
# Add in habitat Air, as soil or phyllsphere, for soil or phyllosphere taxa
sample_data(ITSall_8.5Kfiltered.ps)$HabitatAir[which(sample_data(ITSall_8.5Kfiltered.ps)$sampleType == "phyllosphere")] <- "phyllosphere"
sample_data(ITSall_8.5Kfiltered.ps)$HabitatAir[which(sample_data(ITSall_8.5Kfiltered.ps)$sampleType == "soil")] <- "soil"
# Double check:
cbind(sample_data(ITSall_8.5Kfiltered.ps)$HabitatAir, sample_data(ITSall_8.5Kfiltered.ps)$sampleType)

# ITSall_8.5Kfiltered.ps
# otu_table()   OTU Table:         [ 13778 taxa and 324 samples ]
# sample_data() Sample Data:       [ 324 samples by 41 sample variables ]
# tax_table()   Taxonomy Table:    [ 13778 taxa by 7 taxonomic ranks ]

# Saved September 30, 2025
# saveRDS(ITSall_8.5Kfiltered.ps, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airITS_8.5Kfiltered_phyloseq.rds")

######
# GET CONVENIENT SUBSETS OF THE DATA 
#### i. all data, no controls 
allITSr_noControls.ps <- subset_samples(allITS_rarefied.ps, isControl == "sample") #new phyloseq object without controls/blanks
# View(as.data.frame((as.matrix(sample_data(allITSr_noControls.ps))))) #looks correct
meta_allITSr_noControls <- as.data.frame(as.matrix(sample_data(allITSr_noControls.ps)))

# Add in a column that has habitat type for air, but not for soil or phyllosphere (for ease of plotting later!)
airSavIndex <- intersect(which(meta_allITSr_noControls$sampleType == "air"), which(meta_allITSr_noControls$Habitat == "savanna"))
airForIndex <- intersect(which(meta_allITSr_noControls$sampleType == "air"), which(meta_allITSr_noControls$Habitat == "forest"))
meta_allITSr_noControls$HabitatAir <- rep(NA, nrow(meta_allITSr_noControls))
meta_allITSr_noControls$HabitatAir[airSavIndex] <- "savanna"
meta_allITSr_noControls$HabitatAir[airForIndex] <- "forest"
meta_allITSr_noControls$HabitatAir[which(meta_allITSr_noControls$sampleType == "soil")] <- "soil"
meta_allITSr_noControls$HabitatAir[which(meta_allITSr_noControls$sampleType == "phyllosphere")] <- "phyllosphere"

meta_allITSr_noControls <- sample_data(meta_allITSr_noControls) #make this into "sample data" compatible with phyloseq
allITSr_noControls.ps <- merge_phyloseq(allITSr_noControls.ps, meta_allITSr_noControls) 
allITSr_noControls.ps #13144 taxa and 324 samples 
# Remove ASVs without reads (i.e., those that only showed up in controls!)
noControlsZeros <- which(rowSums(otu_table(allITSr_noControls.ps))==0) #these are all of the ASVs that DID NOT show up in samples (only in controls)
unique(rowSums(otu_table(allITSr_noControls.ps)[noControlsZeros,])) #all zeros
noControlsZerosNames <- names(noControlsZeros)
noControlsASVsToKeep <- setdiff(rownames(otu_table(allITSr_noControls.ps)), noControlsZerosNames) 
allITSr_noControls.ps <- prune_taxa(taxa= noControlsASVsToKeep, x=allITSr_noControls.ps) 
allITSr_noControls.ps #13129 taxa and 324 samples
which(rowSums(otu_table(allITSr_noControls.ps)) == 0) #none, as expected

# saved Sept. 6, 2023
# save(allITSr_noControls.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allITSr_noControls_phyloseq") 
# re-saved June 17, 2024 on server
# save(allITSr_noControls.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/allITSr_noControls_phyloseq") 

#### ii. all data, only controls
allITSr_onlyControls.ps <- subset_samples(allITS_rarefied.ps, isControl == "control")
# View(as.data.frame((as.matrix(sample_data(allITSr_onlyControls.ps))))) #looks correct, AND shows that only have 2 field controls (air),
# 1 water extraction blank (soil), and 1 (soil) extraction blank.
onlyControlsZeros <- which(rowSums(otu_table(allITSr_onlyControls.ps))==0) #these are all of the ASVs that DID NOT show up in the retained controls
unique(rowSums(otu_table(allITSr_onlyControls.ps)[onlyControlsZeros,])) #all zeros
onlyControlsZerosNames <- names(onlyControlsZeros)
onlyControlsASVsToKeep <- setdiff(rownames(otu_table(allITSr_onlyControls.ps)), onlyControlsZerosNames) 
allITSr_onlyControlsTrimmed.ps <- prune_taxa(taxa= onlyControlsASVsToKeep, x=allITSr_onlyControls.ps) #152 taxa and 4 samples!
which(rowSums(otu_table(allITSr_onlyControlsTrimmed.ps)) == 0) #none, as expected
# View(as.data.frame(tax_table(allITSr_onlyControlsTrimmed.ps)))

# saved Sept. 6, 2023
# save(allITSr_onlyControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allITSr_onlyControlsTrimmed_phyloseq") 
# Re-saved on server June 17, 2024
# save(allITSr_onlyControlsTrimmed.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/allITSr_onlyControlsTrimmed_phyloseq") 

#### iii. only air, no controls
airITSr_noControls.ps <- subset_samples(allITSr_noControls.ps, sampleType == "air")
ASVsairITSr_noControls <- as.data.frame(as.matrix(otu_table(airITSr_noControls.ps)))
zeroASVs <- rownames(ASVsairITSr_noControls[which(rowSums(ASVsairITSr_noControls)==0),]) #ASVs with no appearances in (air) data!
zeroIndexCheck2 <- which(rownames(ASVsairITSr_noControls) %in% zeroASVs == TRUE) #double check that index is correct
unique(rowSums(ASVsairITSr_noControls[zeroIndexCheck2,])) #all zeros. So these names are correct

# Get ASVs that are not = to zero
ASVsNotZeroInAirNames <- setdiff(rownames(ASVsairITSr_noControls), zeroASVs) 
length(ASVsNotZeroInAirNames) # 5,440 ASVs that aren't zero!
# Double check, again before subsetting data
sort(rowSums(ASVsairITSr_noControls[which(rownames(ASVsairITSr_noControls) %in% ASVsNotZeroInAirNames == TRUE),])) #great, none of these are zero

# Using ASVsNotZeroInAirNames, drop ASVs not found in air
airITSr_noControlsTrimmed.ps <- prune_taxa(taxa= ASVsNotZeroInAirNames, x=airITSr_noControls.ps) 
airITSr_noControlsTrimmed.ps #this now has 5,440 taxa across 110 samples!

meta_airITSr_noControlsTrimmed <- as.data.frame(as.matrix(sample_data(airITSr_noControlsTrimmed.ps))) 
meta_airITSr_noControlsTrimmed$sampleType #just air!

## Add in a column to the metadata that takes into account the date set out (with first date set out as day 1)
DateSetOuts <- sort(unique(meta_airITSr_noControlsTrimmed$DateSetOut)) #this is the actual date the UPAS samplers were set out.
daysOut <- vector(length=length(DateSetOuts)) #this refers to the chronological order of the days, counting
# both days sampling and days off. In other words, June 16 is "1" because start of sampling, and the start of
# the second sampling period, June 18, is "3" since it occurred on the 3rd calendar day of the sampling period.
daysOut_df <- as.data.frame(cbind(DateSetOuts, daysOut))
daysOut_df$daysOut <- c(16, 1, 3, 17, 5, 7, 8, 9, 11, 12, 14, 15, 19, 20, 21, 23) #organized in Excel sheet called
# samplingDateOrganization

# Make nested for loops with an if statement to add this to metadata
for (i in 1:nrow(meta_airITSr_noControlsTrimmed)){ #loop over rows in metadata, 110 (equal to number of samples remaining)
  for (j in 1:nrow(daysOut_df)){ #equal to 16, as in the total # of days sampled
    # if actual date of each sample (so outer loop, 110 samples) is equal to the jth DateSetOut (n=16), then
    # daysOut in the metadata gets the daysOut (i.e., chronological order) associated with the DateSetOut
    if (meta_airITSr_noControlsTrimmed$DateSetOut[i] == daysOut_df$DateSetOuts[j]) {
      meta_airITSr_noControlsTrimmed$daysOut[i] <- daysOut_df$daysOut[j]
    }
  }
}
# Check that for loop worked-- looks good against Excel sheet referenced above
cbind(unique(meta_airITSr_noControlsTrimmed$DateSetOut), unique(meta_airITSr_noControlsTrimmed$daysOut))
cbind(meta_airITSr_noControlsTrimmed$DateSetOut, meta_airITSr_noControlsTrimmed$daysOut,meta_airITSr_noControlsTrimmed$EU)

airITSr_noControlsTrimmed_meta <- sample_data(meta_airITSr_noControlsTrimmed) #make this metadata formatted for phyloseq
airITSr_noControlsTrimmed.ps <- merge_phyloseq(airITSr_noControlsTrimmed.ps, airITSr_noControlsTrimmed_meta)
# sample_data(airITSr_noControlsTrimmed.ps)
unique(sample_data(airITSr_noControlsTrimmed.ps)$sampleType)
unique(sample_data(airITSr_noControlsTrimmed.ps)$isControl)
min(rowSums(otu_table(airITSr_noControlsTrimmed.ps))) #great, lowest number of times an ASV shows up is 1
# sample_data(airITSr_noControlsTrimmed.ps)

# saved September 6, 2023
# save(airITSr_noControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airITSr_noControls_trimmed_phyloseq") 
# re-saved June 18, 2024 on server
# save(airITSr_noControlsTrimmed.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/airITSr_noControls_trimmed_phyloseq") 

#### iv. only air, no controls, no singletons or doubletons 
# How many taxa does each air sample have? Not that many but seems accurate!
airITS_richObserved <- estimate_richness(airITSr_noControlsTrimmed.ps, split = TRUE, measures = "Observed")
mean(airITS_richObserved$Observed) #670.8636
max(airITS_richObserved$Observed) #1242
min(airITS_richObserved$Observed) #106

airSingletons <- which(rowSums(otu_table(airITSr_noControlsTrimmed.ps))== 1) 
length(airSingletons) #400
airSingletonsNames <- names(airSingletons) 
airDoubletons <- which(rowSums(otu_table(airITSr_noControlsTrimmed.ps))== 2)
length(airDoubletons) #367
airDoubletonsNames <- names(airDoubletons)
airSingDoubNames <- c(airDoubletonsNames, airSingletonsNames)
airSingDoubNames

# As is unsurprising, most ASVs are rare
ASVoccurrences <- sort(rowSums(as.data.frame(as.matrix(otu_table(airITSr_noControlsTrimmed.ps)))), decreasing = TRUE)
length(ASVoccurrences)
ASVoccurrencesPlot <-barplot(ASVoccurrences, main="ITS: number of reads per air ASV ", ylab= "number of reads")

# Remove singletons and doubletons from air and make a new dataset
airSingDoubIndex <- which(rownames(otu_table(airITSr_noControlsTrimmed.ps)) %in% airSingDoubNames)
unique(rowSums(otu_table(airITSr_noControlsTrimmed.ps)[airSingDoubIndex,])) #these are correct since these are all 1 or 2

airNoSingOrDoubsNames <- setdiff(rownames(otu_table(airITSr_noControlsTrimmed.ps)), airSingDoubNames)
length(airNoSingOrDoubsNames) #4,673 are not singletons or doubletons 

airITS_noSingDoubs.ps <- prune_taxa(airNoSingOrDoubsNames, airITSr_noControlsTrimmed.ps) #now only 4,673 taxa
airITS_noSingDoubs.ps

airITS_richObserved_noSingDoubs <- estimate_richness(airITS_noSingDoubs.ps, split = TRUE, measures = "Observed")
mean(airITS_richObserved_noSingDoubs$Observed) #662.8091 (as it be expected, the richness in each sample went down
# a little with the removal of singletons/doubletons, but not much (was 670.8636))

# saved September 6, 2023
# save(airITS_noSingDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airITS_noSingDoubs_phyloseq") 
sort(colnames(sample_data(airITS_noSingDoubs.ps)))
load(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airITS_noSingDoubs_phyloseq")
# saved Feb 14, 2024
#saveRDS(airITS_noSingDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airOnly_ITS_noAirSingsDoubs_phyloseq.rds")

# re-saved June 18, 2024 on server (and as an rds!)
# saveRDS(airITS_noSingDoubs.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/airOnly_ITS_noAirSingsDoubs_phyloseq.rds") 

# re-saved September 17, 2025 on computer to double check (and as an rds!)
# saveRDS(airITS_noSingDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/airOnly_ITS_noAirSingsDoubs_ps.rds") 

#### vi. Phyllosphere only dataset
phylloITSr_noControls.ps <- subset_samples(allITSr_noControls.ps, sampleType == "phyllosphere")
phylloZeros <- which(rowSums(otu_table(phylloITSr_noControls.ps))==0)
unique(rowSums(otu_table(phylloITSr_noControls.ps)[phylloZeros,])) #all zeros
phylloZerosNames <- names(phylloZeros)
phylloASVsToKeep <- setdiff(rownames(otu_table(phylloITSr_noControls.ps)), phylloZerosNames) 
length(phylloASVsToKeep) #4,703
phylloITSr_noControlsTrimmed.ps <- prune_taxa(taxa= phylloASVsToKeep, x=phylloITSr_noControls.ps) #4,703 taxa and 59 samples
which(rowSums(otu_table(phylloITSr_noControlsTrimmed.ps)) == 0) #none, as expected

# saved September 6, 2023
# save(phylloITSr_noControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/phylloITSr_noControlsTrimmed_phyloseq") 
# re-saved June 18, 2024 on server (and as an rds!)
# saveRDS(phylloITSr_noControlsTrimmed.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/phylloITSr_noControlsTrimmed_phyloseq.rds") 

# re-saved September 17, 2025 on computer to double check (and as an rds!)
# saveRDS(phylloITSr_noControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/phylloITSr_noControlsTrimmed_ps.rds") 

#### vii. soils only dataset 
soilITSr_noControls.ps <- subset_samples(allITSr_noControls.ps, sampleType == "soil")
soilZeros <- which(rowSums(otu_table(soilITSr_noControls.ps))==0)
unique(rowSums(otu_table(soilITSr_noControls.ps)[soilZeros,])) #all zeros
soilZerosNames <- names(soilZeros)
soilASVsToKeep <- setdiff(rownames(otu_table(soilITSr_noControls.ps)), soilZerosNames) 
length(soilASVsToKeep) #5,612 #WOW! This is higher than phyllosphere, but not by a crazy amount!
soilITSr_noControlsTrimmed.ps <- prune_taxa(taxa= soilASVsToKeep, x=soilITSr_noControls.ps) #5,612 taxa and 155 samples
which(rowSums(otu_table(soilITSr_noControlsTrimmed.ps)) == 0) #none, as expected

#### viii. only samples (air w/ sings and doubs removed, all soil, all phyllo (i.e. sings and doubs not removed for soil
# and phyllosphere))
allITSr_noAirSingsDoubs.ps <- merge_phyloseq(airITS_noSingDoubs.ps, soilITSr_noControlsTrimmed.ps, phylloITSr_noControlsTrimmed.ps)

# saved September 6, 2023
# save(soilITSr_noControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/soilITSr_noControlsTrimmed_phyloseq") 
# saved September 25, 2023
# re-saved June 18, 2024 on server (and as an rds!)
# saveRDS(soilITSr_noControlsTrimmed.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/soilITSr_noControlsTrimmed_phyloseq.rds")

# save(allITSr_noAirSingsDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allITSr_noAirSingsDoubs.ps") 
# saved February 20, 2024 (server needs .rds data file)
#saveRDS(allITSr_noAirSingsDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allITSr_noAirSingsDoubs.rds")
# Re-saved June 18, 2024 on server (and as an rds!)
# saveRDS(allITSr_noAirSingsDoubs.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/allITSr_noAirSingsDoubs.rds") 

# re-saved September 17, 2025 on computer to double check (and as an rds!)
# saveRDS(soilITSr_noControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/soilITSr_noControlsTrimmed_ps.rds") 
# saveRDS(allITSr_noAirSingsDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/sallITSr_noAirSingsDoubs_ps.rds") 


##############################################################
# FIRST VISUALIZATIONS AND STATS OF ALL DATA
##############################################################
############## ORDINATIONS ##############
#### 1. "Big" ordination of all remaining samples (including blanks)
# Shows that bioaerosol field controls are closest to air, but are on their own.
set.seed(21)
allITSrarefied.ord <- ordinate(allITS_rarefied.ps, "NMDS", "bray")
ITS_rarefiedOrd1 <- plot_ordination(allITS_rarefied.ps, allITSrarefied.ord, type="samples", color="sampleType") 
ITS_rarefiedOrd1 + geom_polygon(aes(fill=sampleType), alpha = 1/5) + geom_point(size=5) + ggtitle("samples") + theme_bw()

#### 2. PERMANOVAS
# Are samples different based on sample type (i.e., air, soil, versus phyllosphere)? Need to do Bonferroni corrections each
# time 
# i. SET UP FOR PERMANOVAS (Get subsets of data, B-C distances of each, and metadata for each type.)
# Soil versus foliar surfaces
ITSr_SFS.ps <- subset_samples(allITSr_noAirSingsDoubs.ps, sampleType != "air") #only soil and f.s.
unique(sample_data(ITSr_SFS.ps)$sampleType) #"soil", "phyllosphere"
ITSr_SFS_ASVs <- t(as.data.frame(as.matrix(otu_table(ITSr_SFS.ps)))) #get ASV table
if (ncol(ITSr_SFS_ASVs) < nrow(ITSr_SFS_ASVs)) {
  stop("samples need to be rows for vegdist!!") #samples must be rows!!
} else {
  ITSr_SFS_BCs <- vegdist(ITSr_SFS_ASVs, method = "bray") # get B-C distances
}
ITSr_SFS_meta <- as.data.frame(as.matrix(sample_data(ITSr_SFS.ps)))  #get metadata

# Bioaerosols versus soil
ITSr_BAS.ps <- subset_samples(allITSr_noAirSingsDoubs.ps, sampleType != "phyllosphere") #only soil and f.s.
unique(sample_data(ITSr_BAS.ps)$sampleType) #"air"  "soil"
ITSr_BAS_ASVs <- t(as.data.frame(as.matrix(otu_table(ITSr_BAS.ps)))) #get ASV table
if (ncol(ITSr_BAS_ASVs) < nrow(ITSr_BAS_ASVs)) {
  stop("samples need to be rows for vegdist!!") #samples must be rows!!
} else {
  ITSr_BAS_BCs <- vegdist(ITSr_BAS_ASVs, method = "bray") # get B-C distances
}
ITSr_BAS_meta <- as.data.frame(as.matrix(sample_data(ITSr_BAS.ps)))  #get metadata

# Bioaerosols versus phyllosphere
ITSr_BAFS.ps <- subset_samples(allITSr_noAirSingsDoubs.ps, sampleType != "soil") #only air and f.s.
unique(sample_data(ITSr_BAFS.ps)$sampleType) #"air", "phyllosphere"
ITSr_BAFS_ASVs <- t(as.data.frame(as.matrix(otu_table(ITSr_BAFS.ps)))) #get ASV table
if (ncol(ITSr_BAFS_ASVs) < nrow(ITSr_BAFS_ASVs)) {
  stop("samples need to be rows for vegdist!!") #samples must be rows!!
} else {
  ITSr_BAFS_BCs <- vegdist(ITSr_BAFS_ASVs, method = "bray") # get B-C distances
}
ITSr_BAFS_meta <- as.data.frame(as.matrix(sample_data(ITSr_BAFS.ps)))  #get metadata

# ii. PERFORM PERMANOVAS WITH BONFERRONI CORRECTIONS
# Soil versus foliar surfaces
set.seed(1121)
ITSr_SFS_PERMANOVA <- adonis2(ITSr_SFS_BCs ~ sampleType, data=ITSr_SFS_meta, permutations= 9999)
ITSr_SFS__Ps_BonF <- p.adjust(ITSr_SFS_PERMANOVA$`Pr(>F)`, method = "bonferroni", n=3)
ITSr_SFS__Ps_BonF #3e-04

# Bioaerosol verus soil
set.seed(1121)
ITSr_BAS_PERMANOVA <- adonis2(ITSr_BAS_BCs ~ sampleType, data=ITSr_BAS_meta, permutations= 9999)
ITSr_BAS_Ps_BonF <- p.adjust(ITSr_BAS_PERMANOVA$`Pr(>F)`, method = "bonferroni", n=3)
ITSr_BAS_Ps_BonF #3e-04

# Bioaerosol verus foliar
set.seed(1121)
ITSr_BAFS_PERMANOVA <- adonis2(ITSr_BAFS_BCs ~ sampleType, data=ITSr_BAFS_meta, permutations= 9999)
ITSr_BAFS_Ps_BonF <- p.adjust(ITSr_BAFS_PERMANOVA$`Pr(>F)`, method = "bonferroni", n=3)
ITSr_BAFS_Ps_BonF #3e-04

#### 3. Do samples (non-control) cluster together based on EU?  -- (all samples). There isn't a super clear pattern.
allITSr_EUs_ordPlot <- plot_ordination(allITSr_noAirSingsDoubs.ps, allITSr_onlySamps.ord, type="samples", color="sampleType", shape = "EU") 
# quartz()
allITSr_EUs_ordPlot + theme_bw() + geom_point(size=3.5)

# Bracketed by EU 
allITSr_EUs_ordPlot2 <- plot_ordination(allITSr_noAirSingsDoubs.ps, allITSr_onlySamps.ord, type="samples", color="sampleType", shape="Habitat") 
# quartz()
allITSr_EUs_ordPlot2 + facet_wrap(~EU, 2) + geom_point(size=3.5) + theme_bw() +                                                                # Change font size
  theme(strip.text.x = element_text(size = 12))

#### 4. Ordination of only non-controls- this shows that air, soil, and the phyllosphere are *very* distinct 
# (i.e. no overlap) which is overall expected and good. Furthermore, it suggests that air samples are closer
# to the phyllosphere samples than they are to the soil!
unique(as.data.frame((as.matrix(sample_data(allITSr_noAirSingsDoubs.ps))))$sampleType) #shows no controls
allITSr_onlySamps.ord <- ordinate(allITSr_noAirSingsDoubs.ps, "NMDS", "bray") # 0.09624302 
allITSr_onlySamps_ordPlot <- plot_ordination(allITSr_noAirSingsDoubs.ps, allITSr_onlySamps.ord, type="samples", color="sampleType") 
allITSr_onlySamps_ordPlot + geom_polygon(aes(fill=sampleType)) + geom_point(size=5) + ggtitle("samples")
allITSr_onlySamps_ordPlot + theme_bw() + geom_polygon(aes(fill=sampleType, alpha = 1/5)) + geom_point(size=5) + ggtitle("NMDS based on sample origin") + theme_bw()

# New version for manuscript (Sept. 17, 2025), added to figureBeautifying.R
ITS_bySampTypeOrd <- allITSr_onlySamps_ordPlot + 
  scale_color_manual(values =c("cornflowerblue", "forestgreen","chocolate"), 
                     labels = c("bioaerosol\nn = 110 (ITS)\nn = 84 (16S)\n", "foliar surfaces\nn=59 (ITS)\nn = 58 (16S)\n", "soil\nn = 155 (ITS)\nn = 157 (16S)"),
                     name = NULL) + #remove legend title name
  geom_point(size=2) +
  theme_bw() +
  theme(panel.grid = element_blank())
  
ITS_bySampTypeOrd
# Saved October 25, 2025
# saveRDS(ITS_bySampTypeOrd, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_bySampTypeOrd_10-25-2025.rds")

#### 5. Just air- Do samples cluster together based on EU?  -- (just air samples)
# Plot shows that there does not seem to be a pattern based on EU. 
# Get better colors
colors4 <- glasbey.colors(6) #do 17 so that I can avoid white, which is G1 nd G5, which looks black
swatch(colors4)
colsForEUs <- unname(colors4)[c(2:4,6)]

airITSr_onlySamps.ord <- ordinate(airITS_noSingDoubs.ps, "NMDS", "bray")
airITSr_EUs_ordPlot <- plot_ordination(airITS_noSingDoubs.ps, airITSr_onlySamps.ord, type="samples", color="EU", shape= "Habitat") + 
  theme_bw() + geom_point(size=4) + theme(axis.text.x=element_text(size=13)) + theme(axis.text.y=element_text(size=13)) + scale_color_manual(values = colsForEUs)

# quartz()
airITSr_EUs_ordPlot + geom_point(size=3.5) + theme_bw()

#### 6. Do air samples seem to differ based on savanna versus forest (all days, all EUs)??
# NOPE!
airITSr_EUs_ordPlot2 <- plot_ordination(airITS_noSingDoubs.ps, airITSr_onlySamps.ord, type="samples", color= "Habitat", shape= "Habitat") + geom_point(size=4) + theme_bw()
airITSr_EUs_ordPlot2

#### 6. Exploring time -- (just air samples)
# a. separate out by facets, depending on time
unique(sample_data(airITS_noSingDoubs.ps)$DateSetOut) #16 unique days out (as expected!!!)
## With day set out as color ##
# Get good, sequential colors for different orders
colors16 <- glasbey.colors(17) #do 17 so that I cn avoid white, which is G1
swatch(colors16)
colorsForDaysOut <- unname(colors16)[2:17] #remove first color, which looks white

airITSr_Ord_byDateSetOut_ordPlot <- plot_ordination(airITS_noSingDoubs.ps, airITSr_onlySamps.ord, type="samples", color= "DateSetOut", shape= "Habitat") + 
  scale_color_manual(values = colorsForDaysOut) + theme_bw() + geom_point(size=4) + theme(axis.text.x=element_text(size=13)) + theme(axis.text.y=element_text(size=13))

# quartz()
grid.arrange(airITSr_EUs_ordPlot, airITSr_Ord_byDateSetOut_ordPlot, ncol=2)

## Faceted by date ##
# This lists the daysOut by EU, with the four dates for EU 52, EU 53S, EU54S, 
# and EU 8, respectively.
daysOutByEU <- c("16-Jun-2022", "18-Jun-2022", "26-Jun-2022", "1-Jul-2022", "22-Jun-2022", "27-Jun-2022", "2-Jul-2022",
                 "6-Jul-2022", "24-Jun-2022", "30-Jun-2022", "5-Jul-2022", "8-Jul-2022", "20-Jun-2022", "23-Jun-2022",
                 "29-Jun-2022", "4-Jul-2022")
airITSr_Ord_byDateSetOut_ordPlot_faceted <- airITSr_EUs_ordPlot2 + facet_wrap(~factor(DateSetOut, c("16-Jun-2022", "18-Jun-2022", "26-Jun-2022", "1-Jul-2022", "22-Jun-2022", "27-Jun-2022", "2-Jul-2022",
                                                                                                    "6-Jul-2022", "24-Jun-2022", "30-Jun-2022", "5-Jul-2022", "8-Jul-2022", "20-Jun-2022", "23-Jun-2022",
                                                                                                    "29-Jun-2022", "4-Jul-2022"))) + geom_point(size=3.5) + theme_bw() +                                                                # Change font size
  theme(strip.text.x = element_text(size = 12)) +
  theme(legend.key.size = unit(2.3, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size

# quartz()
airITSr_Ord_byDateSetOut_ordPlot_faceted

# b. Do samples cluster together based on time?
# maybe put wave thingys on top of plot, like in code from Julian
# see help here: https://chrischizinski.github.io/rstats/ordisurf/
# quartz()
airITSr_EUs_ordPlot + geom_point(size=3.5) + theme_bw()

# c. As another way of doing this (based on Noah's suggestions), split all of the days into separate phyloseq objects 
# and make ordinations for these. I added a tryCatch so that everything would run (even if some days failed).
# IMPORTANT: these are made in order of all of EU 52's dates, then all of EU 8's dates, etc.
ITS_byDate_ps_list <- vector("list", length=length(daysOutByEU)) # pre-allocate list; not that daysOutByEU in NOT in complete chronological order;
# i.e., it has the chronological order of the 4 EU 52 dates, then the four from EU 8, etc.
names(ITS_byDate_ps_list)
ITS_byDateOrds_list <- vector("list", length=length(daysOutByEU)) # pre-allocate list 
ITS_byDatePlot_list <- vector("list", length=length(daysOutByEU)) # pre-allocate list
for (j in 1:length(ITS_byDate_ps_list)) {
  tryCatch({
    ITS_byDate_ps_list[[j]] <- subset_samples(airITS_noSingDoubs.ps, DateSetOut == daysOutByEU[j])
    names(ITS_byDate_ps_list)[[j]] <- daysOutByEU[j]
    ITS_byDateOrds_list[[j]] <- ordinate(ITS_byDate_ps_list[[j]], "NMDS", "bray")
    ITS_byDatePlot_list[[j]] <- plot_ordination(ITS_byDate_ps_list[[j]], ITS_byDateOrds_list[[j]], type="samples", color="Habitat") +
      scale_color_manual(values = c("darkgreen", "goldenrod")) + theme(axis.title.y = element_blank()) + theme_bw() + geom_point(size=4)
  }, error=function(e){})
}
names(ITS_byDate_ps_list) #this confirms that they are daysOutByEU
# These are all the numbers of samples in ITS_byDate_ps_list. They equal 84, which is the number of total samples
# 5 + 7 + 5 + 2 + 6 + 4 + 5 + 4 + 4 + 7 + 6 + 5 + 6 + 6 + 6 + 6

# Fix these figures by removing axis titles and legend (couldn't do it above for some reason!)
ITS_byDatePlot_list[[1]] <- ITS_byDatePlot_list[[1]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[2]] <- ITS_byDatePlot_list[[2]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[3]] <- ITS_byDatePlot_list[[3]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[4]] <- ITS_byDatePlot_list[[4]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[5]] <- ITS_byDatePlot_list[[5]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[6]] <- ITS_byDatePlot_list[[6]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[7]] <- ITS_byDatePlot_list[[7]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[8]] <- ITS_byDatePlot_list[[8]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[9]] <- ITS_byDatePlot_list[[9]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[10]] <- ITS_byDatePlot_list[[10]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[11]] <- ITS_byDatePlot_list[[11]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[12]] <- ITS_byDatePlot_list[[12]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[13]] <- ITS_byDatePlot_list[[13]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[14]] <- ITS_byDatePlot_list[[14]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[15]] <- ITS_byDatePlot_list[[15]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
ITS_byDatePlot_list[[16]] <- ITS_byDatePlot_list[[16]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")

# Plot them all!

# quartz()
grid.arrange(ITS_byDatePlot_list[[1]], ITS_byDatePlot_list[[2]], ITS_byDatePlot_list[[3]], ITS_byDatePlot_list[[4]],
             ITS_byDatePlot_list[[5]], ITS_byDatePlot_list[[6]], ITS_byDatePlot_list[[7]], ITS_byDatePlot_list[[8]],
             ITS_byDatePlot_list[[9]], ITS_byDatePlot_list[[10]], ITS_byDatePlot_list[[11]], ITS_byDatePlot_list[[12]],
             ITS_byDatePlot_list[[13]], ITS_byDatePlot_list[[14]], ITS_byDatePlot_list[[15]], ITS_byDatePlot_list[[16]], ncol=4)

# Test stats
# First pre-allocate some stuff:
ITS_metaDat_list <- vector("list", length=length(daysOutByEU)) # pre-allocate list 
ITS_permanResults_list <- vector("list", length=length(daysOutByEU)) # pre-allocate list 
# Then run for loop to actually run statistics
set.seed(19)
for (k in 1:length(ITS_byDate_ps_list)){
  ITS_metaDat_list[[k]] <- as(sample_data(ITS_byDate_ps_list[[k]]), "data.frame")
  names(ITS_metaDat_list)[[k]] <- names(ITS_byDate_ps_list)[[k]]
  ITS_permanResults_list[[k]] <- adonis2(phyloseq::distance(ITS_byDate_ps_list[[k]], method="bray") ~ HabitatAir,
                                         data = ITS_metaDat_list[[k]])
}
# Make an object to hold all of the Fs and Ps
ITS_permanResultsByDate <- rbind(unlist(ITS_permanResults_list[[1]][,4:5])[c(1,4)], unlist(ITS_permanResults_list[[2]][,4:5])[c(1,4)],
                                 unlist(ITS_permanResults_list[[3]][,4:5])[c(1,4)], unlist(ITS_permanResults_list[[4]][,4:5])[c(1,4)],
                                 unlist(ITS_permanResults_list[[5]][,4:5])[c(1,4)], unlist(ITS_permanResults_list[[6]][,4:5])[c(1,4)],
                                 unlist(ITS_permanResults_list[[7]][,4:5])[c(1,4)], unlist(ITS_permanResults_list[[8]][,4:5])[c(1,4)],
                                 unlist(ITS_permanResults_list[[9]][,4:5])[c(1,4)], unlist(ITS_permanResults_list[[10]][,4:5])[c(1,4)],
                                 unlist(ITS_permanResults_list[[11]][,4:5])[c(1,4)], unlist(ITS_permanResults_list[[12]][,4:5])[c(1,4)],
                                 unlist(ITS_permanResults_list[[13]][,4:5])[c(1,4)], unlist(ITS_permanResults_list[[14]][,4:5])[c(1,4)],
                                 unlist(ITS_permanResults_list[[15]][,4:5])[c(1,4)], unlist(ITS_permanResults_list[[16]][,4:5])[c(1,4)])
rownames(ITS_permanResultsByDate) <- daysOutByEU

ITS_permanResultsByDate 
# F1 Pr(>F)1
# 16-Jun-2022 1.1142514   0.200
# 18-Jun-2022 0.8880555   0.813
# 26-Jun-2022 1.5305107   0.022 #SIGNIFICANT
# 1-Jul-2022  1.2822814   0.032 #SIGNIFICANT
# 22-Jun-2022 0.7097026   1.000
# 27-Jun-2022 1.0770781   0.500
# 2-Jul-2022  0.9563617   0.589
# 6-Jul-2022  0.9732587   0.461
# 24-Jun-2022 0.9722853   0.536
# 30-Jun-2022 1.1883763   0.033 #SIGNIFICANT
# 5-Jul-2022  1.0351936   0.400
# 8-Jul-2022  0.8786437   0.861
# 20-Jun-2022 0.9672045   0.571
# 23-Jun-2022 1.6385249   0.100
# 29-Jun-2022 1.1583629   0.198
# 4-Jul-2022  1.3245284   0.059 #MARGINALLY SIGNIFICANT

#####################
# TOP TAXA IN AIR (BY HABITAT), SOIL, AND PHYLLOSPHERE SAMPLES
#####################
# For air, doing this for forest and savanna. For soil and phyllosphere, keeping them separate.
# So four panels
sample_data(allITSr_noAirSingsDoubs.ps)

##### CLASS #####
allITSr_noAirSingsDoubs.ps.class.glom <-  tax_glom(allITSr_noAirSingsDoubs.ps, taxrank = "Class") 
dim(tax_table(allITSr_noAirSingsDoubs.ps.class.glom)) # 58 orders

# Transform sample counts on just glommed samples (unlike what we did at first)
allITSr_noAirSingsDoubs.ps.class.0 <- transform_sample_counts(allITSr_noAirSingsDoubs.ps.class.glom , function(x) x / sum(x) )
rownames(otu_table(allITSr_noAirSingsDoubs.ps.class.0 )) #ASVs are just a representative from each class

# Merge samples so that we only have combined abundances for forest (air), savanna (air), phyllosphere, and soil
allITSr_relabun.class.1 <- merge_samples(allITSr_noAirSingsDoubs.ps.class.0 , group = "HabitatAir")
sample_data(allITSr_relabun.class.1) #sample categories are correct

# Convert to proportions again b/c total abundance of each site will equal number of species that were merged
allITSr_relabun.class.2 <- transform_sample_counts(allITSr_relabun.class.1, function(x) x / sum(x))
sample_data(allITSr_relabun.class.2)

# Now get only taxa that comprise a certain proportion of the abundance (the rest will be grouped together as one color)
allITSr_relabun.class.df <-psmelt(allITSr_relabun.class.2)
dim(allITSr_relabun.class.df) 
colnames(allITSr_relabun.class.df) #metadata variables 
allITSr_relabun.class.top99 <- allITSr_relabun.class.df 
allITSr_relabun.class.top95 <- allITSr_relabun.class.df 

allITSr_relabun.class.top99$Class[allITSr_relabun.class.top99$Abundance < 0.01] <- "< 1% abund."
allITSr_relabun.class.top95$Class[allITSr_relabun.class.top95$Abundance < 0.05] <- "< 5% abund."

allITSr_top_99p_class <- unique(allITSr_relabun.class.top99$Class)
allITSr_top_99p_class 

# PLOT!
allITSr_class.plot <- ggplot(data=allITSr_relabun.class.top99, aes(x=Sample, y=Abundance, fill=Class)) + theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold"))
# quartz()()
allITSr_class.plot + geom_bar(aes(), stat="identity", position="fill") +
  #scale_fill_manual(values = c("#a50026", "#fee090",  "#74add1", "#313695","#5e4fa2", "#f46d43", "#fdae61", "#4575b4", "#FFFF00", "#EE82EE","#8B008B", "#9ACD32", "#CD5C5C", "grey"), 
  #name= "Family", breaks= c("D_4__Enterococcaceae", "D_4__Moraxellaceae", "D_4__Porphyromonadaceae", "D_4__Enterobacteriaceae", "D_4__Planococcaceae", "D_4__Flavobacteriaceae", "D_4__Micrococcaceae", "D_4__Streptococcaceae", "D_4__Comamonadaceae", "D_4__Nocardiaceae", "D_4__Fusobacteriaceae", "D_4__Family XI", "D_4__Rickettsiales Incertae Sedis", "< 1% abund."), 
  #labels =c("Enterococcaceae", "Moraxellaceae", "Porphyromonadaceae", "Enterobacteriaceae", "Planococcaceae", "Flavobacteriaceae", "Micrococcaceae", "Streptococcaceae", "Comamonadaceae", "Nocardiaceae", "Fusobacteriaceae", "Family XI", "Rickettsiales Incertae Sedis", "< 1% abund.")) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5)) + theme(legend.text = element_text(colour="black", size = 16))  + theme(legend.title = element_blank())

top_classes <- allITSr_relabun.class.df %>%
  group_by(Sample, Class) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) %>%
  #View()
  
  ###################################
##### ORDER #####
require("magrittr")
allITSr_noAirSingsDoubs.ps.order.glom <-  phyloseq::tax_glom(allITSr_noAirSingsDoubs.ps, taxrank = "Order") 
length(unique(tax_table(allITSr_noAirSingsDoubs.ps.order.glom)[,4])) #264

# Transform sample counts on just glommed samples (unlike what we did at first)
allITSr_noAirSingsDoubs.ps.order.0 <- transform_sample_counts(allITSr_noAirSingsDoubs.ps.order.glom , function(x) x / sum(x) )
#ASVs in this above are just a representative from each class

# Merge samples so that we only have combined abundances for site and different kinds of controls
allITSr_relabun.order.1 <- merge_samples(allITSr_noAirSingsDoubs.ps.order.0 , group = "HabitatAir") 
sample_data(allITSr_relabun.order.1) #sample categories are correct

# Convert to proportions again b/c total abundance of each site will equal number of species that were merged
allITSr_relabun.order.2 <- transform_sample_counts(allITSr_relabun.order.1, function(x) x / sum(x))
sample_data(allITSr_relabun.order.2)

# Now get only taxa that comprise a certain proportion of the abundance (the rest will be grouped together as one color)
allITSr_relabun.order.df <-psmelt(allITSr_relabun.order.2)
colnames(allITSr_relabun.order.df) #metadata variables 
allITSr_relabun.order.top99 <- allITSr_relabun.order.df 
allITSr_relabun.order.top98 <- allITSr_relabun.order.df 
allITSr_relabun.order.top95 <- allITSr_relabun.order.df 

allITSr_relabun.order.top99$Order[allITSr_relabun.order.top99$Abundance < 0.01] <- "< 1% abund."
allITSr_relabun.order.top98$Order[allITSr_relabun.order.top98$Abundance < 0.02] <- "< 2% abund."
allITSr_relabun.order.top95$Order[allITSr_relabun.order.top95$Abundance < 0.05] <- "< 5% abund."

allITSr_top_99p_order <- unique(allITSr_relabun.order.top99$Order)
length(allITSr_top_99p_order) #32
allITSr_top_98p_order <- unique(allITSr_relabun.order.top98$Order)
length(allITSr_top_98p_order) #27
allITSr_top_95p_order <- unique(allITSr_relabun.order.top95$Order)
length(allITSr_top_95p_order) #14

# Tidy up the dataframe a little bit before plotting
allITSr_relabun.order.top99$Order <- gsub(pattern= "Microbotryomycetes_ord_Incertae_sedis", x=allITSr_relabun.order.top99$Order, replacement = "Microbotryomycetes, inc. sed.")
allITSr_relabun.order.top99$Order <- gsub(pattern= "o__", x=allITSr_relabun.order.top99$Order, replacement = "") #remove silly "
unique(allITSr_relabun.order.top99$Order) #fixed
# drop taxa unclassified at order level
indexNAs <- which(allITSr_relabun.order.top99$Order == "NA")
allITSr_relabun.order.top99 <- allITSr_relabun.order.top99[-indexNAs,]
unique(allITSr_relabun.order.top99$Order) #NO NAs

# Get good, sequential colors for different orders
length(unique(allITSr_relabun.order.top99$Order)) #31 unique
colors31 <- glasbey.colors(32)
swatch(colors31)
colorsForOrders <- unname(colors31)[2:32] #remove first color, which looks white

# PLOT!
fungiOrder.plot <- ggplot(data=allITSr_relabun.order.top99, aes(x=Sample, y=Abundance, fill=Order)) + theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold"))
positions <- c("forest", "savanna", "phyllosphere", "soil") #for re-rodering 
# quartz()
fungiOrder.plot + geom_bar(aes(), stat="identity", position="fill") + theme_bw() +
  scale_fill_manual(values = colorsForOrders) + ylab("Relative abundance") + theme(axis.title.y = element_text(size = 14)) + xlab("") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=8)) + theme(legend.text = element_text(colour="black", size = 10))  + theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(size=14)) + theme(axis.text.y =element_text(size=14)) + scale_x_discrete(limits = positions)

allITSr_top_orders <- allITSr_relabun.order.df %>%
  group_by(Sample, Order) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

allITSr_top_orders
