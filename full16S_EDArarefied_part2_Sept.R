# full16S_EDArarefied_part2_Sept.R
# 16S full EDA rarefied - PART 2

# Script description: This script follows full16S_EDA_part1_Sept.R (which cleans data and gets rarefied dataset).
# Specifically, First: got convenient subsets of the data (samples only, controls only, air only, phyllosphere only, 
# soil only, singletons and doubletons removed from air, etc. ). Second: NMDS ordinations and visualizations of 
# different subsets of the data. Third, PERMANOVAs showuing that samples differ based on sample Type. Fourth,
# (preliminary, i.e. not the final version used in MS) permanovas to look for difference between forest and savanna (in several 
# different ways).Fourth: plots and calculations of top bacterial classes.

# Original version of this script was called: full16S_EDArarefied_Part2.R. Originally run on server
# Differences between that script and this one include new treatment of phyllopshere (removing contaminants) and no more in-depth 
# dispersion or wind effects analyses. Also, no more visualizing top bacterial orders (since not in manuscript)
# ADDITIONAL NOTE: This script gets certain important visualizations for the BROADN symposium, July 13, 2023.

# Orignal versions of objects saved in previous versions of this script: 1) all16Sr_noAirSingsDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_noAirSingsDoubs.ps_phyloseq"

#######################################
#             SCRIPT SET UP
#######################################
# Read in libraries
library("phyloseq") 
library("tidyverse") 
library("vegan") 
library("gridExtra")  
library("Polychrome") 
library("Matrix")
library("reshape2")

# Load rarefied data:
load(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_dcAP_rarefied_ps_Sept15") #called I6S_dcAP_rarefied.ps. Made in full16S_EDA_part1_Sept.R
I6S_dcAP_rarefied.ps #name of object just loaded

length(grep(x= colnames(otu_table(I6S_dcAP_rarefied.ps)), pattern = "air")) #92 air samples/blanks, as expected. = 84 samples + 8 blanks

# Load unrarefied data (called I6Sall_d_foliarAir.ps)
load(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6Sall_d_foliarAir_ps_Sept15") 
I6Sall_d_foliarAir.ps

######
# NON-RAREFIED: 
####
# 1. REMOVE SAMPLES WITH FEWER THAN 5,500 READS (TO MATCH RAREFIED VERSION)
# Identify the samples to keep
samplesToKeep_16S <- sample_sums(I6Sall_d_foliarAir.ps) >= 5500 #sample_sums() returns the total counts for each sample
# Prune out these samples
I6Sall_d_5.5K.ps <- prune_samples(samplesToKeep_16S, I6Sall_d_foliarAir.ps)
min(sample_sums(I6Sall_d_5.5K.ps)) #5517 is new minimum

# 2. AIR SUBSETTED DATASET: REMOVE BIOAEROSOL SINGELTONS AND DOUBLETONS 
I6Sall_d_5.5K_justAir.ps <- subset_samples(I6Sall_d_5.5K.ps, sampleType == "air")
sample_data(I6Sall_d_5.5K_justAir.ps)$sampleType
# Remove taxa without reads
I6Sall_d_5.5K_justAir.ps <- prune_taxa(taxa_sums(I6Sall_d_5.5K_justAir.ps) > 0, I6Sall_d_5.5K_justAir.ps)
# Find singletons and doubletons (do not do < because we want to keep those that don't appear because they are in soil or phyllosphere!)
I6S_NotRareSingletons <- which(rowSums(otu_table(I6Sall_d_5.5K_justAir.ps)) == 1) 
I6S_NotRareDoubletons <- which(rowSums(otu_table(I6Sall_d_5.5K_justAir.ps)) == 2) 
notRareSingDoubs <- c(names(I6S_NotRareSingletons), names(I6S_NotRareDoubletons))
unique(rowSums(otu_table(I6Sall_d_5.5K_justAir.ps)[rownames(otu_table(I6Sall_d_5.5K_justAir.ps)) %in% notRareSingDoubs,])) #confirms that all are one or two!
# Prune out these taxa
goodNotRareAirASVs16S <- setdiff(rownames(otu_table(I6Sall_d_5.5K_justAir.ps)),notRareSingDoubs)
I6Sall_d_5.5K_justAir.ps <- prune_taxa(x= I6Sall_d_5.5K_justAir.ps, taxa = goodNotRareAirASVs16S)

# Add in a column that has habitat type for air
I6S_airOnly_notR_meta <- as.data.frame(as.matrix(sample_data(I6Sall_d_5.5K_justAir.ps)))
airSavIndex <- intersect(which(I6S_airOnly_notR_meta$sampleType == "air"), which(I6S_airOnly_notR_meta$Habitat == "savanna"))
airForIndex <- intersect(which(I6S_airOnly_notR_meta$sampleType == "air"), which(I6S_airOnly_notR_meta$Habitat == "forest"))
I6S_airOnly_notR_meta$HabitatAir <- rep(NA, nrow(I6S_airOnly_notR_meta))
I6S_airOnly_notR_meta$HabitatAir[airSavIndex] <- "savanna"
I6S_airOnly_notR_meta$HabitatAir[airForIndex] <- "forest"
I6S_airOnly_notR_meta$HabitatAir[which(I6S_airOnly_notR_meta$sampleType == "soil")] <- "soil"
I6S_airOnly_notR_meta$HabitatAir[which(I6S_airOnly_notR_meta$sampleType == "phyllosphere")] <- "phyllosphere"

I6S_airOnly_notR_meta <- sample_data(I6S_airOnly_notR_meta) #make this into "sample data" compatible with phyloseq
I6Sall_d_5.5K_justAir.ps <- merge_phyloseq(I6Sall_d_5.5K_justAir.ps, sample_data(I6S_airOnly_notR_meta)) #299 samples
sort(colnames(sample_data(I6Sall_d_5.5K_justAir.ps))) #shows HabitatAir was incorporated!

# Saved September 30, 2025
# saveRDS(I6Sall_d_5.5K_justAir.ps, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6Sall_d_5.5K_justAir_phyloseq.rds")
# otu_table()   OTU Table:         [ 4573 taxa and 84 samples ]
# sample_data() Sample Data:       [ 84 samples by 41 sample variables ]
# tax_table()   Taxonomy Table:    [ 4573 taxa by 6 taxonomic ranks ]

# 3. ALL I6S WITH SINGLETONS AND DOUBLETONS FOUND ABOVE REMOVED FOR BIOAEROSOLS
# (Here remove leftover controls too )
I6Sall_5.5K_noAir <- subset_samples(I6Sall_d_5.5K.ps, sampleType!= "air")
unique(sample_data(I6Sall_5.5K_noAir)$sampleType)
I6Sall_5.5K_noAir <- subset_samples(I6Sall_5.5K_noAir, sampleType!= "washBufferControl")
I6Sall_5.5K_noAir <- subset_samples(I6Sall_5.5K_noAir, sampleType!= "fieldControl")
I6Sall_5.5K_noAir <- subset_samples(I6Sall_5.5K_noAir, sampleType!= "phyllo_NegFieldControl")
unique(sample_data(I6Sall_5.5K_noAir)$sampleType) # "phyllosphere" "soil"        
# Make a new phyloseq
I6Sall_5.5Kfiltered.ps <- merge_phyloseq(I6Sall_5.5K_noAir, I6Sall_d_5.5K_justAir.ps)
# Add in habitat Air, as soil or phyllsphere, for soil or phyllosphere taxa
sample_data(I6Sall_5.5Kfiltered.ps)$HabitatAir[which(sample_data(I6Sall_5.5Kfiltered.ps)$sampleType == "phyllosphere")] <- "phyllosphere"
sample_data(I6Sall_5.5Kfiltered.ps)$HabitatAir[which(sample_data(I6Sall_5.5Kfiltered.ps)$sampleType == "soil")] <- "soil"
# Double check:
cbind(sample_data(I6Sall_5.5Kfiltered.ps)$HabitatAir, sample_data(I6Sall_5.5Kfiltered.ps)$sampleType)
# Remove taxa without reads
I6Sall_5.5Kfiltered.ps <- prune_taxa(taxa_sums(I6Sall_5.5Kfiltered.ps) > 0, I6Sall_5.5Kfiltered.ps)

I6Sall_5.5Kfiltered.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 32744 taxa and 299 samples ]
# sample_data() Sample Data:       [ 299 samples by 41 sample variables ]
# tax_table()   Taxonomy Table:    [ 32744 taxa by 6 taxonomic ranks ]

# Saved September 30, 2025
# saveRDS(I6Sall_5.5Kfiltered.ps, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airI6S_5.5Kfiltered_phyloseq.rds")

######
# RAREFIED: GET CONVENIENT SUBSETS OF THE DATA 
#### i. all data, no controls 
all16Sr_noControls.ps <- subset_samples(I6S_dcAP_rarefied.ps, isControl == "sample") #new phyloseq object without controls/blanks
# View(as.data.frame((as.matrix(sample_data(all16Sr_noControls.ps))))) #looks correct
meta_r_noControls <- as.data.frame(as.matrix(sample_data(all16Sr_noControls.ps)))

# Add in a column that has habitat type for air, but not for soil or phyllosphere (for ease of plotting later!)
airSavIndex <- intersect(which(meta_r_noControls$sampleType == "air"), which(meta_r_noControls$Habitat == "savanna"))
airForIndex <- intersect(which(meta_r_noControls$sampleType == "air"), which(meta_r_noControls$Habitat == "forest"))
meta_r_noControls$HabitatAir <- rep(NA, nrow(meta_r_noControls))
meta_r_noControls$HabitatAir[airSavIndex] <- "savanna"
meta_r_noControls$HabitatAir[airForIndex] <- "forest"
meta_r_noControls$HabitatAir[which(meta_r_noControls$sampleType == "soil")] <- "soil"
meta_r_noControls$HabitatAir[which(meta_r_noControls$sampleType == "phyllosphere")] <- "phyllosphere"

meta_all16Sr_noControls <- sample_data(meta_r_noControls) #make this into "sample data" compatible with phyloseq
all16Sr_noControls.ps <- merge_phyloseq(all16Sr_noControls.ps, sample_data(meta_r_noControls)) #299 samples
sort(colnames(sample_data(all16Sr_noControls.ps))) #shows HabitatAir was incorporated!

# Remove ASVs without reads (i.e., those that only showed up in controls!)
noControlsZeros <- which(rowSums(otu_table(all16Sr_noControls.ps))==0) #these are all of the ASVs that DID NOT show up in samples (only in controls)
length(noControlsZeros) #129 ASVs were found only in controls.
unique(rowSums(otu_table(all16Sr_noControls.ps)[noControlsZeros,])) #all zeros
noControlsZerosNames <- names(noControlsZeros)
noControlsASVsToKeep <- setdiff(rownames(otu_table(all16Sr_noControls.ps)), noControlsZerosNames) 
all16Sr_noControls.ps <- prune_taxa(taxa= noControlsASVsToKeep, x=all16Sr_noControls.ps) #remove these ASVs
all16Sr_noControls.ps #30091 taxa and 299 samples
which(rowSums(otu_table(all16Sr_noControls.ps)) == 0) #none, as expected

# saved Dec. 4, 2023 on personal computer
# save(all16Sr_noControls.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_noControls_phyloseq") 
# re-saved June 18, 2024 on server (note that this version has ASVs without reads in non-control samples removed)
# save(all16Sr_noControls.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/all16Sr_noControls_phyloseq") 

# re-saved September 17, 2025 on personal comp (note that this version has ASVs without reads in non-control samples removed)
# save(all16Sr_noControls.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_noControls_ps_sept25") 

#### ii. all data, only controls
all16Sr_onlyControls.ps <- subset_samples(I6S_dcAP_rarefied.ps, isControl == "control")
# View(as.data.frame((as.matrix(sample_data(all16Sr_onlyControls.ps))))) #looks correct, AND shows that only have 6 field controls (air),
# 2 wash buffer controls (air), 3 (phyllo) negative field control
onlyControlsZeros <- which(rowSums(otu_table(all16Sr_onlyControls.ps))==0)
unique(rowSums(otu_table(all16Sr_onlyControls.ps)[onlyControlsZeros,])) #all zeros
onlyControlsZerosNames <- names(onlyControlsZeros)
onlyControlsASVsToKeep <- setdiff(rownames(otu_table(all16Sr_onlyControls.ps)), onlyControlsZerosNames) 
all16Sr_onlyControlsTrimmed.ps <- prune_taxa(taxa= onlyControlsASVsToKeep, x=all16Sr_onlyControls.ps) #501 taxa and 11 samples
which(rowSums(otu_table(all16Sr_onlyControlsTrimmed.ps)) == 0) #none, as expected

# saved Dec. 4, 2023
# save(all16Sr_onlyControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_onlyControlsTrimmed_phyloseq") 
# re-saved June 11, 2024 on server
# save(all16Sr_onlyControlsTrimmed.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/all16Sr_onlyControlsTrimmed_phyloseq") 
# re-saved September 17, 2025 on personal comp (note that this version has ASVs without reads in non-control samples removed)
# save(all16Sr_onlyControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_onlyControlsTrimmed_ps_sept25") 

#### iii. only air, no controls
air16Sr_noControls.ps <- subset_samples(all16Sr_noControls.ps, sampleType == "air")
ASVsair16Sr_noControls <- as.data.frame(as.matrix(otu_table(air16Sr_noControls.ps)))
zeroASVs <- rownames(ASVsair16Sr_noControls[which(rowSums(ASVsair16Sr_noControls)==0),]) #ASVs with no appearances in data (but are in soil and/or phyllosphere!
zeroIndexCheck2 <- which(rownames(ASVsair16Sr_noControls) %in% zeroASVs == TRUE) #double check that index is correct
unique(rowSums(ASVsair16Sr_noControls[zeroIndexCheck2,])) #all zeros. So these names are correct
sample_data(air16Sr_noControls.ps)$sampleType #only air
sample_data(air16Sr_noControls.ps)$isControl #no controls!

# Get ASVs that are not = to zero
ASVsNotZeroInAirNames <- setdiff(rownames(ASVsair16Sr_noControls), zeroASVs) 
length(ASVsNotZeroInAirNames) # 4,641 ASVs that aren't zero!
# Double check, again before subsetting data
min(rowSums(ASVsair16Sr_noControls[which(rownames(ASVsair16Sr_noControls) %in% ASVsNotZeroInAirNames == TRUE),])) #great, none of these are zero (min is 1)

# Using ASVsNotZeroInAirNames, drop ASVs not found in air
air16Sr_noControlsTrimmed.ps <- prune_taxa(taxa= ASVsNotZeroInAirNames, x=air16Sr_noControls.ps) #this now has 4641 taxa!
dim(otu_table(air16Sr_noControlsTrimmed.ps)) #4641, 84
length(ASVsNotZeroInAirNames) == dim(otu_table(air16Sr_noControlsTrimmed.ps))[1] #TRUE so the correct ASVs were filtered out.
meta_air16Sr_noControlsTrimmed <- as.data.frame(as.matrix(sample_data(air16Sr_noControlsTrimmed.ps))) 
unique(meta_air16Sr_noControlsTrimmed$sampleType) #just air!

## Add in a column to the metadata that takes into account the date set out (with first date set out as day 1)
DateSetOuts <- sort(unique(meta_air16Sr_noControlsTrimmed$DateSetOut)) #this is the actual date the UPAS samplers were set out.
daysOut <- vector(length=length(DateSetOuts)) #this refers to the chronological order of the days, counting
# both days sampling and days off. In other words, June 16 is "1" because start of sampling, and the start of
# the second sampling period, June 18, is "3" since it occurred on the 3rd calendar day of the sampling period.
daysOut_df <- as.data.frame(cbind(DateSetOuts, daysOut))
daysOut_df$daysOut <- c(16, 1, 3, 17, 5, 7, 8, 9, 11, 12, 14, 15, 19, 20, 21, 23) #organized in Excel sheet called
# samplingDateOrganization
daysOut_df #double checked, this looks correct

# Make nested for loops with an if statement to add this to metadata
for (i in 1:nrow(meta_air16Sr_noControlsTrimmed)){ #loop over rows in metadata, 84 (equal to number of samples remaining)
  for (j in 1:nrow(daysOut_df)){ #equal to 16, as in the total # of days sampled
    # if actual date of each sample (so outer loop, 84 samples) is equal to the jth DateSetOut (n=16), then
    # daysOut in the metadata gets the daysOut (i.e., chronological order) associated with the DateSetOut
    if (meta_air16Sr_noControlsTrimmed$DateSetOut[i] == daysOut_df$DateSetOuts[j]) {
      meta_air16Sr_noControlsTrimmed$daysOut[i] <- daysOut_df$daysOut[j]
    }
  }
}
# Check that for loop worked-- looks good against Excel sheet referenced above
cbind(unique(meta_air16Sr_noControlsTrimmed$DateSetOut), unique(meta_air16Sr_noControlsTrimmed$daysOut))
cbind(meta_air16Sr_noControlsTrimmed$DateSetOut, meta_air16Sr_noControlsTrimmed$daysOut, meta_air16Sr_noControlsTrimmed$EU)

air16Sr_noControlsTrimmed_meta <- sample_data(meta_air16Sr_noControlsTrimmed) #make this metadata formatted for phyloseq
air16Sr_noControlsTrimmed.ps <- merge_phyloseq(air16Sr_noControlsTrimmed.ps, air16Sr_noControlsTrimmed_meta)
# sample_data(air16Sr_noControlsTrimmed.ps)
unique(sample_data(air16Sr_noControlsTrimmed.ps)$sampleType) #air
unique(sample_data(air16Sr_noControlsTrimmed.ps)$isControl) #sample
min(rowSums(otu_table(air16Sr_noControlsTrimmed.ps))) #great, lowest number of times an ASV shows up is 1

# saved Dec. 4, 2023
# save(air16Sr_noControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/air16Sr_noControls_trimmed_phyloseq") 
# re-saved June 11, 2024 on server
# save(air16Sr_noControlsTrimmed.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/air16Sr_noControls_trimmed_phyloseq") 

# re-saved September 17, 2025 on personal comp (note that this version has ASVs without reads in non-control samples removed)
# save(air16Sr_noControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/air16Sr_noControlsTrimmed.ps_sept25") 

#### iv. only air, no controls, no singletons or doubletons 
# How many taxa does each air sample have? Not many :(
air16S_richObserved <- estimate_richness(air16Sr_noControlsTrimmed.ps, split = TRUE, measures = "Observed")
mean(air16S_richObserved$Observed) #143.2381
max(air16S_richObserved$Observed) #534
min(air16S_richObserved$Observed) #54

airSingletons <- which(rowSums(otu_table(air16Sr_noControlsTrimmed.ps))== 1) 
length(airSingletons) #345
airSingletonsNames <- names(airSingletons) 
airDoubletons <- which(rowSums(otu_table(air16Sr_noControlsTrimmed.ps))== 2) 
length(airDoubletons) #191
airDoubletonsNames <- names(airDoubletons)
airSingDoubNames <- c(airDoubletonsNames, airSingletonsNames)
air3times <- which(rowSums(otu_table(air16Sr_noControlsTrimmed.ps))== 3)
length(air3times) #149
air4times <- which(rowSums(otu_table(air16Sr_noControlsTrimmed.ps))== 4)
length(air4times) #127
air5times <- which(rowSums(otu_table(air16Sr_noControlsTrimmed.ps))== 5)
length(air5times) #141

# As is unsurprising, most ASVs are rare
ASVoccurrences <- sort(rowSums(as.data.frame(as.matrix(otu_table(air16Sr_noControlsTrimmed.ps)))), decreasing = TRUE)
length(ASVoccurrences) #4641, as it should be
ASVoccurrencesPlot <-barplot(ASVoccurrences, main="16S: number of reads per air ASV ", ylab= "number of reads")

# Remove singletons and doubletons from air and make a new dataset
airSingDoubIndex <- which(rownames(otu_table(air16Sr_noControlsTrimmed.ps)) %in% airSingDoubNames)
unique(rowSums(otu_table(air16Sr_noControlsTrimmed.ps)[airSingDoubIndex,])) #these are correct since these are all 1 or 2
# Get subset that aren't singletons or doubeltons
airNoSingOrDoubsNames <- setdiff(rownames(otu_table(air16Sr_noControlsTrimmed.ps)), airSingDoubNames)

air16S_noSingDoubs.ps <- prune_taxa(airNoSingOrDoubsNames, air16Sr_noControlsTrimmed.ps) #now only 4105 taxa

air16S_richObserved_noSingDoubs <- estimate_richness(air16S_noSingDoubs.ps, split = TRUE, measures = "Observed")
mean(air16S_richObserved_noSingDoubs$Observed) #136.3929
max(air16S_richObserved_noSingDoubs$Observed) #481
min(air16S_richObserved_noSingDoubs$Observed) #51

# re-saved Feb 14, 2024 because HabitatAir stuff above was not incorporated!
# save(air16S_noSingDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/air16S_noSingDoubs_phyloseq") 

# saved Feb 14, 2024
#saveRDS(air16S_noSingDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airOnly_16Sr_noAirSingsDoubs_phyloseq.rds")

# re-saved June 11, 2024 on server (and as an rds!)
# saveRDS(air16S_noSingDoubs.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/airOnly_16Sr_noAirSingsDoubs_phyloseq.rds") 

# re-saved September 17, 2025 on personal comp (note that this version has ASVs without reads in non-control samples removed)
# saveRDS(air16S_noSingDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/air16S_noSingDoubs_ps_sept25.rds") 

#### vi. Phyllosphere only dataset
phyllo16Sr_noControls.ps <- subset_samples(all16Sr_noControls.ps, sampleType == "phyllosphere")
unique(sample_data(phyllo16Sr_noControls.ps)$sampleType)
phylloZeros <- which(rowSums(otu_table(phyllo16Sr_noControls.ps))==0)
unique(rowSums(otu_table(phyllo16Sr_noControls.ps)[phylloZeros,])) #all zeros, as expected
phylloZerosNames <- names(phylloZeros)
phylloASVsToKeep <- setdiff(rownames(otu_table(phyllo16Sr_noControls.ps)), phylloZerosNames) #get ASV names that are in all names but are not zeros! 
length(phylloASVsToKeep) #6,520
phyllo16Sr_noControlsTrimmed.ps <- prune_taxa(taxa= phylloASVsToKeep, x=phyllo16Sr_noControls.ps) #6520 taxa and 58 samples
which(rowSums(otu_table(phyllo16Sr_noControlsTrimmed.ps)) == 0) #none, as expected

# saved Dec. 4, 2023, 2023
# save(phyllo16Sr_noControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/phyllo16Sr_noControlsTrimmed_phyloseq") 
# re-saved June 11, 2024 on server (and as an rds!)
# saveRDS(phyllo16Sr_noControlsTrimmed.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/phyllo16Sr_noControlsTrimmed_phyloseq.rds") 

# re-saved September 17, 2025 on personal comp (note that this version has ASVs without reads in non-control samples removed)
# saveRDS(phyllo16Sr_noControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/phyllo16Sr_noControlsTrimmed_ps_sept25.rds") 

#### vii. soils only dataset 
soil16Sr_noControls.ps <- subset_samples(all16Sr_noControls.ps, sampleType == "soil")
unique(sample_data(soil16Sr_noControls.ps)$sampleType)
soilZeros <- which(rowSums(otu_table(soil16Sr_noControls.ps))==0)
unique(rowSums(otu_table(soil16Sr_noControls.ps)[soilZeros,])) #all zeros
soilZerosNames <- names(soilZeros)
soilASVsToKeep <- setdiff(rownames(otu_table(soil16Sr_noControls.ps)), soilZerosNames) 
length(soilASVsToKeep) #22,868
soil16Sr_noControlsTrimmed.ps <- prune_taxa(taxa= soilASVsToKeep, x=soil16Sr_noControls.ps) #22845 taxa and 157 samples 
soil16Sr_noControlsTrimmed.ps #22868  taxa, 157 samples
which(rowSums(otu_table(soil16Sr_noControlsTrimmed.ps)) == 0) #none, as expected

#### viii. only samples (air w/ sings and doubs removed, all soil, all phyllo (i.e. sings and doubs NOT removed for soil
# and phyllosphere))
all16Sr_noAirSingsDoubs.ps <- merge_phyloseq(air16S_noSingDoubs.ps, soil16Sr_noControlsTrimmed.ps, phyllo16Sr_noControlsTrimmed.ps) #29640 taxa and 299 samples 
sort(colnames(sample_data(all16Sr_noAirSingsDoubs.ps))) #HabitatAir was retained!

# saved Dec. 4, 2023 (on personal computer)
# save(soil16Sr_noControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/soil16Sr_noControlsTrimmed_phyloseq") 
# re-saved June 11, 2024 on server (and as an rds!)
# saveRDS(soil16Sr_noControlsTrimmed.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/soil16Sr_noControlsTrimmed_phyloseq.rds")

# saved Dec. 4, 2023 (on personal computer)
# save(all16Sr_noAirSingsDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_noAirSingsDoubs.ps") 

# saved Feb. 20, 2024 (needs to be rds for ANCOM on server!)
# saveRDS(all16Sr_noAirSingsDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_noAirSingsDoubs.rds") 

# re-saved June 11, 2024 on server (and as an rds!)
# saveRDS(all16Sr_noAirSingsDoubs.ps, file="~/SRS_aeromicrobiome_2022/RobjectsSaved/all16Sr_noAirSingsDoubs.rds") 

# Re-saved with newest version September 17, 2025
# saveRDS(all16Sr_noAirSingsDoubs.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_noAirSingsDoubs_ps_sept25.rds") 
# saveRDS(soil16Sr_noControlsTrimmed.ps, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/soil16Sr_noControlsTrimmed_ps_sept25.rds") 

##############################################################
# FIRST VISUALIZATIONS OF ALL DATA
##############################################################
############## ORDINATIONS ##############
# NOTE: Although these are useful, they are NOT on Hellinger-transformed data
#### 1. "Big" ordination of all remaining samples
# Overall, this shows that bioaerosol-associated blanks seem undistinguishable from bioaerosol samples,
# whereas phyllo neg. field controls are slightly closest to phyllo, but still congregate by themselves.
# This shows cross-contamination with bioaerosol controls, unsurprising for low biomasss
unique(as.data.frame((as.matrix(sample_data(all16Sr_noAirSingsDoubs.ps))))$sampleType)
set.seed(10)
all16Srarefied.ord <- ordinate(I6S_dcAP_rarefied.ps, "NMDS", "bray")

rarefiedOrd1 <- plot_ordination(I6S_dcAP_rarefied.ps, all16Srarefied.ord, type="samples", color="sampleType") 

# quartz()
rarefiedOrd1 + geom_polygon(aes(fill=sampleType), alpha = 1/5) + geom_point(size=5) + ggtitle("samples") + theme_bw()

#### 2. PERMANOVAS
# Are samples different based on sample type (i.e., air, soil, versus phyllosphere)? Need to do Bonferroni corrections each
# time 
# i. SET UP FOR PERMANOVAS (Get subsets of data, B-C distances of each, and metadata for each type.)
# Soil versus foliar surfaces
I6Sr_SFS.ps <- subset_samples(all16Sr_noAirSingsDoubs.ps, sampleType != "air") #only soil and f.s.
unique(sample_data(I6Sr_SFS.ps)$sampleType) #"soil", "phyllosphere"
I6Sr_SFS_ASVs <- t(as.data.frame(as.matrix(otu_table(I6Sr_SFS.ps)))) #get ASV table
if (ncol(I6Sr_SFS_ASVs) < nrow(I6Sr_SFS_ASVs)) {
  stop("samples need to be rows for vegdist!!") #samples must be rows!!
} else {
  I6Sr_SFS_BCs <- vegdist(I6Sr_SFS_ASVs, method = "bray") # get B-C distances
}
I6Sr_SFS_meta <- as.data.frame(as.matrix(sample_data(I6Sr_SFS.ps)))  #get metadata

# Bioaerosols versus soil
I6Sr_BAS.ps <- subset_samples(all16Sr_noAirSingsDoubs.ps, sampleType != "phyllosphere") #only soil and f.s.
unique(sample_data(I6Sr_BAS.ps)$sampleType) #"air"  "soil"
I6Sr_BAS_ASVs <- t(as.data.frame(as.matrix(otu_table(I6Sr_BAS.ps)))) #get ASV table
if (ncol(I6Sr_BAS_ASVs) < nrow(I6Sr_BAS_ASVs)) {
  stop("samples need to be rows for vegdist!!") #samples must be rows!!
} else {
  I6Sr_BAS_BCs <- vegdist(I6Sr_BAS_ASVs, method = "bray") # get B-C distances
}
I6Sr_BAS_meta <- as.data.frame(as.matrix(sample_data(I6Sr_BAS.ps)))  #get metadata

# Bioaerosols versus phyllosphere
I6Sr_BAFS.ps <- subset_samples(all16Sr_noAirSingsDoubs.ps, sampleType != "soil") #only air and f.s.
unique(sample_data(I6Sr_BAFS.ps)$sampleType) #"air", "phyllosphere"
I6Sr_BAFS_ASVs <- t(as.data.frame(as.matrix(otu_table(I6Sr_BAFS.ps)))) #get ASV table
if (ncol(I6Sr_BAFS_ASVs) < nrow(I6Sr_BAFS_ASVs)) {
  stop("samples need to be rows for vegdist!!") #samples must be rows!!
} else {
  I6Sr_BAFS_BCs <- vegdist(I6Sr_BAFS_ASVs, method = "bray") # get B-C distances
}
I6Sr_BAFS_meta <- as.data.frame(as.matrix(sample_data(I6Sr_BAFS.ps)))  #get metadata

# ii. PERFORM PERMANOVAS WITH BONFERRONI CORRECTIONS
# Soil versus foliar surfaces
set.seed(1121)
I6Sr_SFS_PERMANOVA <- adonis2(I6Sr_SFS_BCs ~ sampleType, data=I6Sr_SFS_meta, permutations= 9999)
I6Sr_SFS__Ps_BonF <- p.adjust(I6Sr_SFS_PERMANOVA$`Pr(>F)`, method = "bonferroni", n=3)
I6Sr_SFS__Ps_BonF #3e-04

# Bioaerosol verus soil
set.seed(1121)
I6Sr_BAS_PERMANOVA <- adonis2(I6Sr_BAS_BCs ~ sampleType, data=I6Sr_BAS_meta, permutations= 9999)
I6Sr_BAS_Ps_BonF <- p.adjust(I6Sr_BAS_PERMANOVA$`Pr(>F)`, method = "bonferroni", n=3)
I6Sr_BAS_Ps_BonF #3e-04

# Bioaerosol verus foliar
set.seed(1121)
I6Sr_BAFS_PERMANOVA <- adonis2(I6Sr_BAFS_BCs ~ sampleType, data=I6Sr_BAFS_meta, permutations= 9999)
I6Sr_BAFS_Ps_BonF <- p.adjust(I6Sr_BAFS_PERMANOVA$`Pr(>F)`, method = "bonferroni", n=3)
I6Sr_BAFS_Ps_BonF #3e-04

#### 3. Ordination of only non-controls- this shows that air, soil, and the phyllosphere are *very* distinct 
# (i.e. no overlap) which is overall expected and good. Furthermore, it suggests that air samples are closer
# to the phyllosphere samples than they are to the soil!
unique(as.data.frame((as.matrix(sample_data(all16Sr_noAirSingsDoubs.ps))))$sampleType)
set.seed(12)
all16Sr_onlySamps.ord <- ordinate(all16Sr_noAirSingsDoubs.ps, "NMDS", "bray") #stress is 0.1148366 
all16Sr_onlySamps_ordPlot <- plot_ordination(all16Sr_noAirSingsDoubs.ps, all16Sr_onlySamps.ord, type="samples", color="sampleType") 
all16Sr_onlySamps_ordPlot + geom_polygon(aes(fill=sampleType)) + geom_point(size=5) + ggtitle("samples")
# quartz()

# New version for manuscript (October 25, 2025), added to figureBeautifying.R
I6S_bySampTypeOrd <- all16Sr_onlySamps_ordPlot +
  scale_color_manual(values =c("cornflowerblue", "forestgreen","chocolate"), 
                    labels = c("bioaerosol\nn = 110 (ITS)\nn = 84 (16S)\n", "foliar surfaces\nn=59 (ITS)\nn = 58 (16S)\n", "soil\nn = 155 (ITS)\nn = 157 (16S)"),
                    name = NULL) + #remove legend title name
  geom_point(size=2) +
  theme_bw() +
  theme(panel.grid = element_blank())
I6S_bySampTypeOrd

# Saved September 17, 2025
# saveRDS(I6S_bySampTypeOrd, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_bySampTypeOrd_10-25-2025.rds")

#### 4. Do samples (non-control) cluster together based on EU?  -- (all samples). There isn't a super clear pattern.
all16Sr_EUs_ordPlot <- plot_ordination(all16Sr_noAirSingsDoubs.ps, all16Sr_onlySamps.ord, type="samples", color="sampleType", shape = "EU") 
# quartz()
all16Sr_EUs_ordPlot + theme_bw() + geom_point(size=3.5)

# Bracketed by EU 
all16Sr_EUs_ordPlot2 <- plot_ordination(all16Sr_noAirSingsDoubs.ps, all16Sr_onlySamps.ord, type="samples", color="sampleType", shape="Habitat") 
# quartz()
all16Sr_EUs_ordPlot2 + facet_wrap(~EU, 2) + geom_point(size=3.5) + theme_bw() +                                                                # Change font size
  theme(strip.text.x = element_text(size = 12))

#### 5. Just air- Do samples cluster together based on EU?  -- (just air samples)
# Plot shows that there does not seem to be a pattern based on EU. 
# Get better colors
colors4 <- glasbey.colors(6) #do 17 so that I can avoid white, which is G1 nd G5, which looks black
swatch(colors4)
colsForEUs <- unname(colors4)[c(2:4,6)]

air16Sr_onlySamps.ord <- ordinate(air16S_noSingDoubs.ps, "NMDS", "bray")
air16Sr_EUs_ordPlot <- plot_ordination(air16S_noSingDoubs.ps, air16Sr_onlySamps.ord, type="samples", color="EU", shape= "Habitat") + 
  theme_bw() + geom_point(size=4) + theme(axis.text.x=element_text(size=13)) + theme(axis.text.y=element_text(size=13)) + scale_color_manual(values = colsForEUs)

# quartz()
air16Sr_EUs_ordPlot + geom_point(size=3.5) + theme_bw()

#### 6. Do air samples seem to differ based on savanna versus forest (all days, all EUs)??
# NOPE! not really 
air16Sr_EUs_ordPlot2 <- plot_ordination(air16S_noSingDoubs.ps, air16Sr_onlySamps.ord, type="samples", color= "Habitat", shape= "Habitat") + geom_point(size=4)
air16Sr_EUs_ordPlot2

#### 7. Exploring time -- (just air samples)
# a. separate out by facets, depending on time
unique(sample_data(air16S_noSingDoubs.ps)$DateSetOut) #16 unique days out (as expected!!!)
## With day set out as color ##
# Get good, sequential colors for different orders
colors16 <- glasbey.colors(17) #do 17 so that I cn avoid white, which is G1
swatch(colors16)
colorsForDaysOut <- unname(colors16)[2:17] #remove first color, which looks white

air16Sr_Ord_byDateSetOut_ordPlot <- plot_ordination(air16S_noSingDoubs.ps, air16Sr_onlySamps.ord, type="samples", color= "DateSetOut", shape= "Habitat") + 
  scale_color_manual(values = colorsForDaysOut) + theme_bw() + geom_point(size=4) + theme(axis.text.x=element_text(size=13)) + theme(axis.text.y=element_text(size=13))

# quartz()
grid.arrange(air16Sr_EUs_ordPlot, air16Sr_Ord_byDateSetOut_ordPlot, ncol=2)

## Faceted by date ##
# This lists the daysOut by EU, with the four dates for EU 52, EU 53S, EU54S, 
# and EU 8, respectively.
daysOutByEU <- c("16-Jun-2022", "18-Jun-2022", "26-Jun-2022", "1-Jul-2022", "22-Jun-2022", "27-Jun-2022", "2-Jul-2022",
                 "6-Jul-2022", "24-Jun-2022", "30-Jun-2022", "5-Jul-2022", "8-Jul-2022", "20-Jun-2022", "23-Jun-2022",
                 "29-Jun-2022", "4-Jul-2022")
air16Sr_Ord_byDateSetOut_ordPlot_faceted <- air16Sr_EUs_ordPlot2 + facet_wrap(~factor(DateSetOut, c("16-Jun-2022", "18-Jun-2022", "26-Jun-2022", "1-Jul-2022", "22-Jun-2022", "27-Jun-2022", "2-Jul-2022",
                                                                                                    "6-Jul-2022", "24-Jun-2022", "30-Jun-2022", "5-Jul-2022", "8-Jul-2022", "20-Jun-2022", "23-Jun-2022",
                                                                                                    "29-Jun-2022", "4-Jul-2022"))) + geom_point(size=3.5) + theme_bw() +                                                                # Change font size
  theme(strip.text.x = element_text(size = 12)) +
  theme(legend.key.size = unit(2.3, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size

# quartz()
air16Sr_Ord_byDateSetOut_ordPlot_faceted

# b. Do samples cluster together based on time?
# maybe put wave thingys on top of plot, like in code from Julian
# see help here: https://chrischizinski.github.io/rstats/ordisurf/
# quartz()
air16Sr_EUs_ordPlot + geom_point(size=3.5) + theme_bw()

# c. As another way of doing this (based on Noah's suggestions), split all of the days into separate phyloseq objects 
# and make ordinations for these. I added a tryCatch so that everything would run (even if some days failed)
# NOTE, I COULD ONLY GET 15 ORDINATIONS, PROBABLY BECAUSE 1-Jul-2022` has only two remaining samples!
# THESE are ordered by EU in rows and days 1-4 in EU as columns
byDate_ps_list <- vector("list", length=length(daysOutByEU)) # pre-allocate list 
names(byDate_ps_list) <- daysOutByEU #name them based on day set out, i.e., how they were allocated below
byDateOrds_list <- vector("list", length=length(daysOutByEU)) # pre-allocate list 
names(byDateOrds_list) <- daysOutByEU #name them based on day set out, i.e., how they were allocated below
byDatePlot_list <- vector("list", length=length(daysOutByEU)) # pre-allocate list
names(byDatePlot_list) <- daysOutByEU #name them based on day set out, i.e., how they were allocated below
for (j in 1:length(byDate_ps_list)) { #loop over the 16 days
  tryCatch({ #add tryCatch so that the code keeps going even when looking at broken samples
    # Make a list of phyloseq objects that are divided by daysOutByEU
    byDate_ps_list[[j]] <- subset_samples(air16S_noSingDoubs.ps, DateSetOut == daysOutByEU[j])
    names(byDate_ps_list)[[j]] <- daysOutByEU[j]
    byDateOrds_list[[j]] <- ordinate(byDate_ps_list[[j]], "NMDS", "bray")
    byDatePlot_list[[j]] <- plot_ordination(byDate_ps_list[[j]], byDateOrds_list[[j]], type="samples", color="Habitat") +
      scale_color_manual(values = c("darkgreen", "goldenrod")) + theme(axis.title.y = element_blank()) + theme_bw() + geom_point(size=4) +
      ggtitle(paste0("plot for", names(byDate_ps_list)[[j]]))
  }, error=function(e){})
}
names(byDatePlot_list)
sampsByDate <- matrix(ncol=2, nrow=length(byDate_ps_list))
# These are all the numbers of samples in byDate_ps_list. They equal 84, which is the number of total samples
# 5 + 7 + 5 + 2 + 6 + 4 + 5 + 4 + 4 + 7 + 6 + 5 + 6 + 6 + 6 + 6

# Fix these figures by removing axis titles and legend (couldn't do it above for some reason!)
# Note: one day did not work above, because it had only 2 samples, so that's why there's only 15 below
byDatePlot_list$`16-Jun-2022` <- byDatePlot_list$`16-Jun-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`18-Jun-2022` <- byDatePlot_list$`18-Jun-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`26-Jun-2022` <- byDatePlot_list$`26-Jun-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`1-Jul-2022` <- byDatePlot_list$`1-Jul-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`22-Jun-2022` <- byDatePlot_list$`22-Jun-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`27-Jun-2022` <- byDatePlot_list$`27-Jun-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`2-Jul-2022` <- byDatePlot_list$`2-Jul-2022`+ theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`6-Jul-2022` <- byDatePlot_list$`6-Jul-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`24-Jun-2022` <- byDatePlot_list$`24-Jun-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`30-Jun-2022` <- byDatePlot_list$`30-Jun-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`5-Jul-2022` <- byDatePlot_list$`5-Jul-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`8-Jul-2022` <- byDatePlot_list$`8-Jul-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`20-Jun-2022` <- byDatePlot_list$`20-Jun-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`23-Jun-2022` <- byDatePlot_list$`23-Jun-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`29-Jun-2022` <- byDatePlot_list$`29-Jun-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")
byDatePlot_list$`4-Jul-2022` <- byDatePlot_list$`4-Jul-2022` + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())  + theme(legend.position = "none")

# Plot them all!

# quartz()
do.call(grid.arrange, c(byDatePlot_list, ncol = 4, nrow = 4))

# These are the ones that has enough stress to plot
# quartz()
grid.arrange(byDatePlot_list[[1]], byDatePlot_list[[2]], byDatePlot_list[[3]],
             byDatePlot_list[[5]], byDatePlot_list[[6]], byDatePlot_list[[7]], byDatePlot_list[[8]],
             byDatePlot_list[[9]], ncol=4)

# Test stats
# First pre-allocate some stuff:
metaDat_list <- vector("list", length=length(daysOutByEU)) # pre-allocate list 
names(metaDat_list) <- daysOutByEU
permanResults_list <- vector("list", length=length(daysOutByEU)) # pre-allocate list 
names(permanResults_list) <- daysOutByEU
# Then run for loop to actually run statistics
set.seed(19)
for (k in 1:length(byDate_ps_list)){
  metaDat_list[[k]] <- as(sample_data(byDate_ps_list[[k]]), "data.frame")
  permanResults_list[[k]] <- adonis2(phyloseq::distance(byDate_ps_list[[k]], method="bray") ~ HabitatAir,
                                     data = metaDat_list[[k]])
}

# Make an object to hold all of the Fs and Ps
permanResultsByDate <- do.call(rbind, lapply(1:16, function(i) unlist(permanResults_list[[i]][,4:5])[c(1,4)]))
rownames(permanResultsByDate) <- daysOutByEU

permanResultsByDate 
# F1   Pr(>F)1
# 16-Jun-2022 0.9773623 0.6000000
# 18-Jun-2022 1.1715641 0.0650000 #MARGINAL
# 26-Jun-2022 1.0118131 0.4000000
# 1-Jul-2022         NA        NA
# 22-Jun-2022 0.9489477 0.7333333
# 27-Jun-2022 0.9551194 0.6666667
# 2-Jul-2022  1.0050147 0.4000000
# 6-Jul-2022  0.9992704 1.0000000
# 24-Jun-2022 0.9223469 1.0000000
# 30-Jun-2022 0.9908992 0.5520000
# 5-Jul-2022  0.9695110 0.8000000
# 8-Jul-2022  1.1130889 0.1000000
# 20-Jun-2022 1.0260497 0.2666667
# 23-Jun-2022 1.0463631 0.1000000
# 29-Jun-2022 1.3092004 0.1000000
# 4-Jul-2022  0.9973838 0.5333333

#####################
# TOP TAXA IN AIR (BY HABITAT), SOIL, AND PHYLLOSPHERE SAMPLES
### *** HERE NEED TO UPDATE WITH air16S_noSingDoubs.ps
#####################
# For air, doing this for forest and savanna. For soil and phyllosphere, keeping them separate.
# So four panels
sample_data(all16Sr_noAirSingsDoubs.ps)

##### CLASS #####
all16Sr_noAirSingsDoubs.ps.class.glom <-  tax_glom(all16Sr_noAirSingsDoubs.ps, taxrank = "Class") 
dim(tax_table(all16Sr_noAirSingsDoubs.ps.class.glom)) # 125 orders

# Transform sample counts on just glommed samples (unlike what we did at first)
all16Sr_noAirSingsDoubs.ps.class.0 <- transform_sample_counts(all16Sr_noAirSingsDoubs.ps.class.glom , function(x) x / sum(x) )
rownames(otu_table(all16Sr_noAirSingsDoubs.ps.class.0 )) #ASVs are just a representative from each class

# Merge samples so that we only have combined abundances for forest (air), savanna (air), phyllosphere, and soil
all16Sr_relabun.class.1 <- merge_samples(all16Sr_noAirSingsDoubs.ps.class.0 , group = "HabitatAir")
sample_data(all16Sr_relabun.class.1) #sample categories are correct

# Convert to proportions again b/c total abundance of each site will equal number of species that were merged
all16Sr_relabun.class.2 <- transform_sample_counts(all16Sr_relabun.class.1, function(x) x / sum(x))
sample_data(all16Sr_relabun.class.2)

# Now get only taxa that comprise a certain proportion of the abundance (the rest will be grouped together as one color)
all16Sr_relabun.class.df <-psmelt(all16Sr_relabun.class.2)
dim(all16Sr_relabun.class.df) 
colnames(all16Sr_relabun.class.df) #metadata variables 
all16Sr_relabun.class.top99 <- all16Sr_relabun.class.df 
all16Sr_relabun.class.top95 <- all16Sr_relabun.class.df 

all16Sr_relabun.class.top99$Class[all16Sr_relabun.class.top99$Abundance < 0.01] <- "< 1% abund."
all16Sr_relabun.class.top95$Class[all16Sr_relabun.class.top95$Abundance < 0.05] <- "< 5% abund."

all16Sr_top_99p_class <- unique(all16Sr_relabun.class.top99$Class)
all16Sr_top_99p_class 

# PLOT!
all16Sr_class.plot <- ggplot(data=all16Sr_relabun.class.top99, aes(x=Sample, y=Abundance, fill=Class)) + theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold"))
#quartz()
all16Sr_class.plot + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5)) + theme(legend.text = element_text(colour="black", size = 16))  + theme(legend.title = element_blank())

top_classes <- all16Sr_relabun.class.df %>%
  group_by(Sample, Class) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 
# View(top_classes)
  