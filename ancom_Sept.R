# ANCOM for SRS aeromicrobiome-- bioaerosols in patch versus matrix
# January 2024, updated September 2025

# Description: Using un-rarefied bioaerosol bacterial and fungal datasets, this script performs ANCOMBC2 analyses to look for 
# differentially abundant taxa in the open patch versus the forested matrix habitat. Crucially, the tests look for habitat 
# differences, while also having nested random effects of date sampled within EU. The criteria for considering a differentially
# abundant taxon, whether a structural or "main" ANCOM zero, is that it must occur in at least 10% of samples and comprise at least 
# 0.05% of reads from the respective group.

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
library(dplyr)
library(lmerTest)

# 2. LOAD PHYLOSEQ OBJECTS (made in scripts on my computer called full16S_EDArarefied_Part2_Sept.R and equivalent fungal script)
# ##### RAREFIED #####
# # Read in 16S data
# I6S_ALL_phyloseq <-  readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_noAirSingsDoubs_ps_sept25.rds")
# # Read in ITS data
# ITS_ALL_phyloseq <- readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/sallITSr_noAirSingsDoubs_ps.rds") #yes, this is correct version, even with "sallITS" typo!
# tax_table(ITS_ALL_phyloseq)
# 
# # Phyloseq objects of bioaerosols only
# # air16S_noSingDoubs.ps
# load(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/air16S_noSingDoubs_phyloseq")
# # airITS_noSingDoubs.ps
# load(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airITS_noSingDoubs_phyloseq")
# head(tax_table(airITS_noSingDoubs.ps)) #somehow this order got turned around, but can fix later

##### NON-RAREFIED #####
# NON-RAREFIED DATA AIR DATA 
# (From full16S_EDArarefied_part2_Sept.R, has singletons and doubletons and samples with < 5500 reads removed)
I6Sall_d_5.5K_justAir.ps <- readRDS(file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6Sall_d_5.5K_justAir_phyloseq.rds")
# (From fullITS_EDArarefied_part2_Sept.R, has singletons and doubletons and samples with < 8500 reads removed)
airITS_8.5K.ps <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airITS_8.5K_phyloseq.rds")

# Load unrarefied data for all samples (air singletons and doubletons removed and samples lower than 5.5K reads also removed)
I6Sall_5.5Kfiltered.ps <-  readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airI6S_5.5Kfiltered_phyloseq.rds")
# Load unrarefied data for all samples (air singletons and doubletons removed and samples lower than 8.5K reads alos removed)
ITSall_8.5Kfiltered.ps <-  readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airITS_8.5Kfiltered_phyloseq.rds")

#####################################
# II. BACTERIA: ANCOM FOR OPEN PATCH VERSUS FORESTED MATRIX
####################################
# Note about ANCOM statistics: W-statistic for each taxon, representing the number of pairwise comparisons 
# where the null hypothesis (no difference in relative abundance) is rejected. A higher W-statistic indicates
# a higher likelihood that the taxon is differentially abundant.
# One way of wording: "Individual ASVs that were in differential abundance between cases and controls in both matched 
# and unmatched disease cohorts were identified using the R implementation of ANCOM11 V2.0. In brief, W scores were calculated
# for each ASV, representing the number of times the null hypothesis was rejected based on log ratio abundances with the taxon
# of interest and each taxon within the dataset. W scores that indicated rejection of the null hypothesis for ≥90% of log ratios
# were designated as having an ANCOM threshold of 0.9, while ≥80% corresponds to a 0.8 ANCOM threshold, and so on."

# Based on discussions here and the fact that ANCOMBC2 accounts for differences in sequencing depth
########## 1. FORMAT!! ##########
# i. BACTERIA
# Right now, there is no "ASV level" in the taxonomy table, but this is what I want to do the analysis on. 
# So, I will need to add in a species level to the taxa table 
colnames(tax_table(I6Sall_d_5.5K_justAir.ps))
head(tax_table(I6Sall_d_5.5K_justAir.ps))

# Make a copy to edit
I6S_airOnly_notR_taxSpecies <- as.data.frame(as.matrix(tax_table(I6Sall_d_5.5K_justAir.ps)))
class(I6S_airOnly_notR_taxSpecies) #great, it's not a phyloseq object!
head(I6S_airOnly_notR_taxSpecies)
I6S_notR_ASVnames <- rownames(I6S_airOnly_notR_taxSpecies)
I6S_airOnly_notR_taxSpecies$Species <- rownames(I6S_airOnly_notR_taxSpecies)
head(I6S_airOnly_notR_taxSpecies)
# Turn this new tax table into a tax table for phyloseq
I6S_airOnly_notR_taxSpeciesForPS <- phyloseq::tax_table(I6S_airOnly_notR_taxSpecies)
colnames(I6S_airOnly_notR_taxSpeciesForPS) <- c("Kingdom","Phylum","Class", "Order","Family","Genus", "Species")
rownames(I6S_airOnly_notR_taxSpeciesForPS) <- I6S_notR_ASVnames #re-set rownames
# View(I6S_airOnly_notR_taxSpeciesForPS)
head(I6S_airOnly_notR_taxSpeciesForPS)
class(I6S_airOnly_notR_taxSpeciesForPS)

# Finally, re-make phyloseq object
I6S_airOnly_notR_ANCOM.ps <- merge_phyloseq(I6Sall_d_5.5K_justAir.ps, I6S_airOnly_notR_taxSpeciesForPS)
tax_table(I6S_airOnly_notR_ANCOM.ps)
otu_table(I6S_airOnly_notR_ANCOM.ps)
# ASV table to use later
I6Sair_notR_ASVsTab <- as.data.frame(as.matrix(otu_table(I6S_airOnly_notR_ANCOM.ps)))
which(rowSums(I6Sair_notR_ASVsTab)<1) #none of these ASVs are NOT found in the air!
dim(I6Sair_notR_ASVsTab) #4,573 ASVs, 84 samples
# Get the tax table to use later
I6SairOnly_notR_taxaTab <- as.data.frame(as.matrix(tax_table(I6S_airOnly_notR_ANCOM.ps)))
I6SairOnly_notR_taxaTab <- tibble::rownames_to_column(I6SairOnly_notR_taxaTab, var= "ASV_name")
# Get sample names to use later
I6SsavAirNames <- rownames(sample_data(I6S_airOnly_notR_ANCOM.ps))[which(sample_data(I6S_airOnly_notR_ANCOM.ps)$HabitatAir == "savanna")]
sample_data(I6S_airOnly_notR_ANCOM.ps)$HabitatAir[rownames(sample_data(I6S_airOnly_notR_ANCOM.ps)) %in% I6SsavAirNames] #double check looks good
I6SForAirNames <- rownames(sample_data(I6S_airOnly_notR_ANCOM.ps))[which(sample_data(I6S_airOnly_notR_ANCOM.ps)$HabitatAir == "forest")]

########## 2. USE ANCOM! ##########
# 1. BACTERIA: differences between habitats
# For reference, here is glmmTMB(GenomeEquiv_rounded ~ HabitatAir + (1|EU/samplingDay), data=allI6S_qPCR, family="nbinom2")
# Need nested random effects of sampling day (DateSetOut is equivalent, but it was called samplingDay for qCPR test because of data re-formatting) within EU.  
sort(colnames(sample_data(I6S_airOnly_notR_ANCOM.ps)))
sample_data(I6S_airOnly_notR_ANCOM.ps)$DateSetOut
set.seed(120)
# Grayed out since takes a while, but code works great!
# I6S_habitatAir_nonR_ANCOM <- ancombc2(data = I6S_airOnly_notR_ANCOM.ps, group = "HabitatAir", tax_level = "Species", fix_formula = "HabitatAir", rand_formula = "(1|EU/DateSetOut)", struc_zero = TRUE, neg_lb = FALSE) #bacteria
class(I6S_habitatAir_nonR_ANCOM)

# Saved/run September 30, 2025
# Added _2 for setting seed
# saveRDS(I6S_habitatAir_nonR_ANCOM, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_habitatAir_nonR_ANCOM_Sept2025_2.rds")

# Look at results:
# Load previously saved results if need be
I6S_habitatAir_nonR_ANCOM <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_habitatAir_nonR_ANCOM_Sept2025_2.rds")
# Shows that no taxa were differentially abundant but that all passed sensitivity analysis (so pseudo count addition which ANOCMBC2 uses is okay)
unique(I6S_habitatAir_nonR_ANCOM$res$diff_HabitatAirsavanna) #no differentially abundant taxa
unique(I6S_habitatAir_nonR_ANCOM$res$diff_robust_HabitatAirsavanna) #no differentially abundant taxa
unique(I6S_habitatAir_nonR_ANCOM$res$`passed_ss_(Intercept)`) #all passed sensitivity analysis
unique(I6S_habitatAir_nonR_ANCOM$res$passed_ss_HabitatAirsavanna) #all passed sensitivity analysis

###### 1. PULL OUT ONLY ASV TABLE WITH ALL ANCOM ASVS ######
structZeroForestNames <- I6S_habitatAir_nonR_ANCOM$zero_ind$taxon[which(I6S_habitatAir_nonR_ANCOM$zero_ind$`structural_zero (HabitatAir = forest)` ==TRUE)]
length(structZeroForestNames) #1126, that's a lot
structZeroSavNames <- I6S_habitatAir_nonR_ANCOM$zero_ind$taxon[which(I6S_habitatAir_nonR_ANCOM$zero_ind$`structural_zero (HabitatAir = savanna)`==TRUE)]
length(structZeroSavNames) #1518

# i. Concatenate all of these ASVs together in order to get all of the ANCOM-identifed ASVs (only struc zero since no main)
allI6S_airHab_ANCOM_ASVs <- c(structZeroForestNames, structZeroSavNames)
length(allI6S_airHab_ANCOM_ASVs) == length(unique(allI6S_airHab_ANCOM_ASVs)) #no repeats, as expected
# ii. Get an ASV table with only these ANCOM ASVs
I6S_airHab_ancom_sigASVtab <- I6Sair_notR_ASVsTab[rownames(I6Sair_notR_ASVsTab) %in% allI6S_airHab_ANCOM_ASVs,]
# iii. Transpose so that ASVs are now rows to add in additional ASV info
I6S_airHab_ancom_sigASVtab_t <- as.data.frame(t(I6S_airHab_ancom_sigASVtab))
if (nrow(I6S_airHab_ancom_sigASVtab_t) < ncol(I6S_airHab_ancom_sigASVtab_t)) {
  I6S_airHab_ancom_sigASVtab_t <- t(I6S_airHab_ancom_sigASVtab_t) 
} else {
  print("yay! ASVs are rows")
}
I6S_airHab_ancom_sigASVtab_t <- as.data.frame(I6S_airHab_ancom_sigASVtab_t); class(I6S_airHab_ancom_sigASVtab_t)
setdiff(rownames(I6S_airHab_ancom_sigASVtab_t), allI6S_airHab_ANCOM_ASVs) #looks good

# iv. Add in sum in savanna and matrix samples
# Sum of each ASV in air savanna 
I6S_airHab_ancom_sigASVtab_t$savAirTotal <- rowSums(I6S_airHab_ancom_sigASVtab_t[,which(colnames(I6S_airHab_ancom_sigASVtab_t)%in% I6SsavAirNames)])
# Occupancy in air savanna 
I6S_airHab_ancom_sigASVtab_t$savAirOcc <- rowSums(I6S_airHab_ancom_sigASVtab_t[, which(colnames(I6S_airHab_ancom_sigASVtab_t)%in% I6SsavAirNames)] > 0, na.rm = TRUE)
# Total read counts in forest matrix air
I6S_airHab_ancom_sigASVtab_t$forestAirTotal <- rowSums(I6S_airHab_ancom_sigASVtab_t[,which(colnames(I6S_airHab_ancom_sigASVtab_t) %in% I6SForAirNames)])
# Occupancy in forest matrix air
I6S_airHab_ancom_sigASVtab_t$forestAirOcc <- rowSums(I6S_airHab_ancom_sigASVtab_t[, which(colnames(I6S_airHab_ancom_sigASVtab_t)%in% I6SForAirNames)] > 0, na.rm = TRUE)

# v. Add in which kind of ANCOM taxon it is
I6S_airHab_ancom_sigASVtab_t <- I6S_airHab_ancom_sigASVtab_t %>% 
  tibble::rownames_to_column(var= "ASV_name") %>% #make rownames a column
  mutate(ANCOMcat = case_when(
    ASV_name %in% structZeroForestNames ~ "open patch", #open patch if not in forest
    ASV_name %in% structZeroSavNames ~ "forested matrix", #forested matrix if not in savanna
  )) %>% 
  mutate(resType = case_when(
    #ASV_name %in% c(bioaerosolI6SANCOMASVsName, foliarI6SANCOMASVsName) ~ "mainNotStructZero", no main from $res so skip this!
    ASV_name %in% c(structZeroForestNames, structZeroSavNames) ~ "structZero"
  )
  )
# vi. Perform a few checks:
# 1. Get indices
I6S_HabAirTotal_ASV53_index <- which(I6S_airHab_ancom_sigASVtab_t$ASV_name == "ASV_53")
I6S_HabAirTotal_ASV385_index <- which(I6S_airHab_ancom_sigASVtab_t$ASV_name == "ASV_385")
# 2. Is total number working?
colnames(I6S_airHab_ancom_sigASVtab_t)
I6S_airHab_ancom_sigASVtab_t$savAirTotal[I6S_HabAirTotal_ASV53_index] == sum(I6S_airHab_ancom_sigASVtab_t[I6S_HabAirTotal_ASV53_index, colnames(I6S_airHab_ancom_sigASVtab_t)%in% I6SsavAirNames])
I6S_airHab_ancom_sigASVtab_t$forestAirTotal[I6S_HabAirTotal_ASV385_index] == sum(I6S_airHab_ancom_sigASVtab_t[I6S_HabAirTotal_ASV385_index, colnames(I6S_airHab_ancom_sigASVtab_t)%in% I6SForAirNames])
# 3. Is occupancy working?
I6S_airHab_ancom_sigASVtab_t$forestAirOcc[I6S_HabAirTotal_ASV53_index] == length(which(I6S_airHab_ancom_sigASVtab_t[I6S_HabAirTotal_ASV53_index, colnames(I6S_airHab_ancom_sigASVtab_t)%in% I6SForAirNames]>0))
# 4. Are they labeled correctly?
setdiff(I6S_airHab_ancom_sigASVtab_t$ASV_name[which(I6S_airHab_ancom_sigASVtab_t$resType == "structZero")], c(structZeroForestNames, structZeroSavNames)) #no differences!
# 5. Are ANCOM categories correct?
setdiff(I6S_airHab_ancom_sigASVtab_t$ASV_name[which(I6S_airHab_ancom_sigASVtab_t$ANCOMcat == "open patch")], structZeroForestNames) #no differences!
# View(I6S_airHab_ancom_sigASVtab_t)

# vii. Add in taxonomy!
head(I6SairOnly_notR_taxaTab)
I6S_airHab_ancom_sigASVtab_t <- left_join(I6S_airHab_ancom_sigASVtab_t, I6SairOnly_notR_taxaTab[,1:7], by = "ASV_name") #don't have last column because "species" is just ASV name for this analysis

# viii. Add in the proportion of reads in each category that an ASV contributes
# Use this, made above, to get total number of reads in savanna air and forested matrix samples
head(I6Sair_notR_ASVsTab)
I6StotReadsSavAir <- sum(unname(colSums(I6Sair_notR_ASVsTab[,colnames(I6Sair_notR_ASVsTab) %in% I6SsavAirNames]))) # 692923
length(I6SsavAirNames) #38 samples
I6StotReadsForAir <- sum(unname(colSums(I6Sair_notR_ASVsTab[,colnames(I6Sair_notR_ASVsTab) %in% I6SForAirNames]))) # 885063
length(I6SForAirNames) #46 samples

# Look at what getting 0.0005 or 0.05% of reads would be about (or 5 out of every 10,000 reads). Do it this way since forest had more reads and more samples
round(0.0005*I6StotReadsSavAir) # needs to occur roughly 346 times in savanna 
round(0.0005*I6StotReadsForAir) # needs to occur roughly 443 times in forest

# Add in percent that each ASV comprises
I6S_airHab_ancom_sigASVtab_t <- I6S_airHab_ancom_sigASVtab_t %>% 
  mutate(pctSavAir = savAirTotal/I6StotReadsSavAir*100) %>% 
  mutate(pctForAir = forestAirTotal/I6StotReadsForAir*100)

# DECISION AND FILTERING TIME:
# 1. To be considered diff. abundant, an ASV must be in at about 10% of respective group's samples (4 savanna or 5 forest) AND
# 2. Must be present at a certain threshold (0.05% of reads)
colnames(I6S_airHab_ancom_sigASVtab_t)
##### AIR SAVANNA ####
round(length(I6SsavAirNames)/10) #4 samples 
# Found in at least 4 samples and .05% of reads
savAir_I6S4occ.05pct_Index <- intersect(intersect(which(I6S_airHab_ancom_sigASVtab_t$savAirOcc >=4), which(I6S_airHab_ancom_sigASVtab_t$pctSavAir >=0.05)), which(I6S_airHab_ancom_sigASVtab_t$ANCOMcat == "open patch"))
length(savAir_I6S4occ.05pct_Index) #6
# View(I6S_airHab_ancom_sigASVtab_t[savAir_I6S4occ.05pct_Index,])
I6SsavAirOnlyASVNames <- I6S_airHab_ancom_sigASVtab_t$ASV_name[savAir_I6S4occ.05pct_Index]
# Which samples did these occur in?
sampColsI6S_airHab <- c(2:85) #these are the columns with data, which is important for thing below
I6SsavAirOnlySamplesList <- vector("list", length = length(I6SsavAirOnlyASVNames))
names(I6SsavAirOnlySamplesList) <- I6SsavAirOnlyASVNames
for (i in 1:length(I6SsavAirOnlySamplesList)){
  I6SsavAirOnlySamplesList[[i]] <- colnames(I6S_airHab_ancom_sigASVtab_t[,sampColsI6S_airHab])[which(I6S_airHab_ancom_sigASVtab_t[I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SsavAirOnlySamplesList)[i],sampColsI6S_airHab]>0)]
}

# Testing for loop:
I6S_airHab_ancom_sigASVtab_t[I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SsavAirOnlySamplesList)[1],]
which(I6S_airHab_ancom_sigASVtab_t[I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SsavAirOnlySamplesList)[1],sampColsI6S_airHab]>0)
# Confirmed that these are the samples that names(I6SsavAirOnlySamplesList)[1], "ASV_161", occurs in
colnames(I6S_airHab_ancom_sigASVtab_t[,sampColsI6S_airHab])[which(I6S_airHab_ancom_sigASVtab_t[I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SsavAirOnlySamplesList)[1],sampColsI6S_airHab]>0)]
I6SsavAirOnlySamplesList[[1]] # "air_16S_108" "air_16S_118" "air_16S_132" "air_16S_86" 

# Shows that 3/4 savanna air samples with fecal taxon ([Ruminococcus] torques group) are from same day very close together (EU 52 near center of patch on June 18)
sampsI6S_ASV_1565 <- as.vector(I6SsavAirOnlySamplesList$ASV_1565)
sample_data(I6S_airOnly_notR_ANCOM.ps)$DateSetOut[rownames(sample_data(I6S_airOnly_notR_ANCOM.ps)) %in% sampsI6S_ASV_1565] #"30-Jun-2022" "18-Jun-2022" "18-Jun-2022" "18-Jun-2022"
sample_data(I6S_airOnly_notR_ANCOM.ps)$TransectMeter[rownames(sample_data(I6S_airOnly_notR_ANCOM.ps)) %in% sampsI6S_ASV_1565] #"R0" "B0" "L0" "T0" 
sample_data(I6S_airOnly_notR_ANCOM.ps)$EU[rownames(sample_data(I6S_airOnly_notR_ANCOM.ps)) %in% sampsI6S_ASV_1565] #"EU_54S" "EU_52"  "EU_52"  "EU_52" 
I6S_airHab_ancom_sigASVtab_t[I6S_airHab_ancom_sigASVtab_t$ASV_name %in% "ASV_1565", colnames(I6S_airHab_ancom_sigASVtab_t) %in% sampsI6S_ASV_1565]
# LOADED in 118
#     air_16S_118 air_16S_150 air_16S_39 air_16S_86
#     1308           2          2          1
sum(I6S_airHab_ancom_sigASVtab_t[,colnames(I6S_airHab_ancom_sigASVtab_t) %in% "air_16S_118"]) #1308/2665
1308/sum(I6S_airHab_ancom_sigASVtab_t[,colnames(I6S_airHab_ancom_sigASVtab_t) %in% "air_16S_118"])*100 # = 49.1% in this one sample!!

# Not really a pattern with Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium 
sampsI6S_ASV_2271 <- as.vector(I6SsavAirOnlySamplesList$ASV_2271)
sample_data(I6S_airOnly_notR_ANCOM.ps)$DateSetOut[rownames(sample_data(I6S_airOnly_notR_ANCOM.ps)) %in% sampsI6S_ASV_2271] #""29-Jun-2022" "22-Jun-2022" "18-Jun-2022" "18-Jun-2022"
sample_data(I6S_airOnly_notR_ANCOM.ps)$TransectMeter[rownames(sample_data(I6S_airOnly_notR_ANCOM.ps)) %in% sampsI6S_ASV_2271] #"B0" "L0" "L0" "T0"
sample_data(I6S_airOnly_notR_ANCOM.ps)$EU[rownames(sample_data(I6S_airOnly_notR_ANCOM.ps)) %in% sampsI6S_ASV_2271] #"EU_8"   "EU_53S" "EU_52"  "EU_52" 
I6S_airHab_ancom_sigASVtab_t[I6S_airHab_ancom_sigASVtab_t$ASV_name %in% "ASV_2271", colnames(I6S_airHab_ancom_sigASVtab_t) %in% sampsI6S_ASV_2271]
#       air_16S_154 air_16S_33 air_16S_39 air_16S_86
#           6        383          2          1

# Other four ASVs
# "ASV_161": Verruco in Candidatus Xiphinematobacter
I6S_airHab_ancom_sigASVtab_t[which(I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SsavAirOnlySamplesList)[1]), which(colnames(I6S_airHab_ancom_sigASVtab_t) %in% I6SsavAirOnlySamplesList[[1]])]
#   air_16S_108 air_16S_118 air_16S_132 air_16S_86
#         292           1         168          1

# "ASV_244" Acidobacteriota in Pyrinomonadaceae, RB41 genus
I6S_airHab_ancom_sigASVtab_t[which(I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SsavAirOnlySamplesList)[2]), which(colnames(I6S_airHab_ancom_sigASVtab_t) %in% I6SsavAirOnlySamplesList[[2]])]
#   air_16S_1 air_16S_150 air_16S_2 air_16S_38 air_16S_39
#         51         114         1        126         65

# "ASV_2396" Actino, Streptomycetaceae, Streptacidiphilus
I6S_airHab_ancom_sigASVtab_t[which(I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SsavAirOnlySamplesList)[5]), which(colnames(I6S_airHab_ancom_sigASVtab_t) %in% I6SsavAirOnlySamplesList[[5]])]
# # air_16S_108 air_16S_132 air_16S_150 air_16S_2 air_16S_38 air_16S_46 air_16S_65
#         103           1         287         3          2         87          8

# ASV_3495: Bacteroidota,  Hymenobacteraceae, Hymenobacter
I6S_airHab_ancom_sigASVtab_t[which(I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SsavAirOnlySamplesList)[6]), which(colnames(I6S_airHab_ancom_sigASVtab_t) %in% I6SsavAirOnlySamplesList[[6]])]
# air_16S_130 air_16S_132 air_16S_143 air_16S_150 air_16S_65
#           1           1         527           1          1

# Summary of these: ASV_1565 [Ruminococcus] torques group. 3/4 samples from one day, close to each other. However, very rare except nearly 50% of reads
# in one sample, presumably a pig pooped nearby.
# Maybe of interest: "ASV_244" Acidobacteriota in Pyrinomonadaceae, RB41 genus and "ASV_2396" Actino, Streptomycetaceae, Streptacidiphilus. Found in good amounts in
# 3-4 samples
# Likely not of interest: "ASV_161": Verruco in Candidatus Xiphinematobacter found 2 times in at least 168 samples, but then only 1 in other two
# Not of interest: ASV_2271 (Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium) and ASV_3495: Bacteroidota,  Hymenobacteraceae, Hymenobacter. In one of samples, found at least 383 times, rest no more than 6. Not very interesting

##### FOREST AIR ####
round(length(I6SForAirNames)/10) #5 samples
# Found in at least 5 samples and at least 0.05% of reads
forAir_I6S5occ.05pct_Index <- intersect(intersect(which(I6S_airHab_ancom_sigASVtab_t$forestAirOcc >=5), which(I6S_airHab_ancom_sigASVtab_t$pctForAir >=0.05)), which(I6S_airHab_ancom_sigASVtab_t$ANCOMcat == "forested matrix"))
length(forAir_I6S5occ.05pct_Index) #3
# View(I6S_airHab_ancom_sigASVtab_t[forAir_I6S5occ.05pct_Index,])
# 3 ASVs: 
# ASV_1088: Proteobacteria, Alphaproteobacteria, Caulobacterales, Caulobacteraceae. 
# ASV_1319: Proteobacteria (all else NA), Actinobacteriota, Actinobacteria, Frankiales
# ASV_1599: Proteobacteria, all the rest unclassified
I6SforAirOnlyASVNames <- I6S_airHab_ancom_sigASVtab_t$ASV_name[forAir_I6S5occ.05pct_Index]
# Which samples did these occur in?
sampColsI6S_airHab <- c(2:85) #these are the columns with data, which is important for thing below
I6SforAirOnlySamplesList <- vector("list", length = length(I6SforAirOnlyASVNames))
names(I6SforAirOnlySamplesList) <- I6SforAirOnlyASVNames
for (i in 1:length(I6SforAirOnlySamplesList)){
  I6SforAirOnlySamplesList[[i]] <- colnames(I6S_airHab_ancom_sigASVtab_t[,sampColsI6S_airHab])[which(I6S_airHab_ancom_sigASVtab_t[I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SforAirOnlySamplesList)[i],sampColsI6S_airHab]>0)]
}

# Testing for loop:
I6S_airHab_ancom_sigASVtab_t[I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SforAirOnlySamplesList)[1],]
which(I6S_airHab_ancom_sigASVtab_t[I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SforAirOnlySamplesList)[1],sampColsI6S_airHab]>0)
# Confirmed that these are the samples that names(I6SforAirOnlySamplesList)[1], "ASV_1088", occurs in
colnames(I6S_airHab_ancom_sigASVtab_t[,sampColsI6S_airHab])[which(I6S_airHab_ancom_sigASVtab_t[I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SforAirOnlySamplesList)[1],sampColsI6S_airHab]>0)]
I6SforAirOnlySamplesList[[1]]

# ASV_1088
I6S_airHab_ancom_sigASVtab_t[which(I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SforAirOnlySamplesList)[1]), which(colnames(I6S_airHab_ancom_sigASVtab_t) %in% I6SforAirOnlySamplesList[[1]])]
# air_16S_103 air_16S_105 air_16S_28 air_16S_72 air_16S_88
#       387           5          1         99        152
# ASV_1319
I6S_airHab_ancom_sigASVtab_t[which(I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SforAirOnlySamplesList)[2]), which(colnames(I6S_airHab_ancom_sigASVtab_t) %in% I6SforAirOnlySamplesList[[2]])]
# air_16S_149 air_16S_36 air_16S_49 air_16S_60 air_16S_72
#        1576          1          1          1          1
# ASV_1599
I6S_airHab_ancom_sigASVtab_t[which(I6S_airHab_ancom_sigASVtab_t$ASV_name %in% names(I6SforAirOnlySamplesList)[3]), which(colnames(I6S_airHab_ancom_sigASVtab_t) %in% I6SforAirOnlySamplesList[[3]])]
# air_16S_127 air_16S_141 air_16S_145 air_16S_54 air_16S_90
#         1273           1           1          1          5

# Summary of "forested matrix only": ASV_1319 (Frankiales) and ASV_1599 (unknown Proteobacterium) occur at least 1200 in one sample,
# with the other 4 samples having it 5 times or fewer. The only one of interest, ASV_1088 (Caulobacteraceae), found in 5 samples. 3 of 
# these samples at least 99 times, two 5 or fewer reads. 

#####################################
# III. FUNGI: ANCOM FOR OPEN PATCH VERSUS FORESTED MATRIX
####################################

# Right now, there is a Species level, but some are NAs which is no good in the taxonomy table.
# So I will change all species to ASV name, as with bacteria and then check for species information again later
colnames(tax_table(airITS_8.5K.ps)) #in correct order
head(tax_table(airITS_8.5K.ps))
# Get a tax table (WITH SPECIES) to use later
airITS_8.5K_taxTab <- as.data.frame(as.matrix(tax_table(airITS_8.5K.ps)))
airITS_8.5K_taxTab  <- airITS_8.5K_taxTab %>% 
  tibble::rownames_to_column(var = "ASV_name")
head(airITS_8.5K_taxTab)

# Make a copy to edit
ITS_airOnly_notR_taxSpecies <- as.data.frame(as.matrix(tax_table(airITS_8.5K.ps)))
class(ITS_airOnly_notR_taxSpecies) #great, it's not a phyloseq object!
head(ITS_airOnly_notR_taxSpecies)
ITS_notR_ASVnames <- rownames(ITS_airOnly_notR_taxSpecies)
ITS_airOnly_notR_taxSpecies$Species <- rownames(ITS_airOnly_notR_taxSpecies)
head(ITS_airOnly_notR_taxSpecies)
# Turn this new tax table into a tax table for phyloseq
ITS_airOnly_notR_taxSpeciesForPS <- phyloseq::tax_table(ITS_airOnly_notR_taxSpecies)
colnames(ITS_airOnly_notR_taxSpeciesForPS) <- c("Kingdom","Phylum","Class", "Order","Family","Genus", "Species")
rownames(ITS_airOnly_notR_taxSpeciesForPS) <- ITS_notR_ASVnames #re-set rownames
# View(ITS_airOnly_notR_taxSpeciesForPS)
head(ITS_airOnly_notR_taxSpeciesForPS)
class(ITS_airOnly_notR_taxSpeciesForPS)

# Finally, re-make phyloseq object. Need to brute force this a bit since species made above arren't taking with jsut the merge.
class(ITS_airOnly_notR_taxSpeciesForPS)
ITS_airOnly_notR_ANCOM.ps <- merge_phyloseq(otu_table(airITS_8.5K.ps), sample_data(airITS_8.5K.ps), ITS_airOnly_notR_taxSpeciesForPS)
tax_table(ITS_airOnly_notR_ANCOM.ps) #now looks correct!
ITSair_notR_ASVsTab <- as.data.frame(as.matrix(otu_table(ITS_airOnly_notR_ANCOM.ps)))
which(rowSums(ITSair_notR_ASVsTab)<1) #none of these ASVs are NOT found in the air!
dim(ITSair_notR_ASVsTab) #5612 ASVs,  110 samples
# Get sample names to use later
ITSsavAirNames <- rownames(sample_data(ITS_airOnly_notR_ANCOM.ps))[which(sample_data(ITS_airOnly_notR_ANCOM.ps)$HabitatAir == "savanna")]
sample_data(ITS_airOnly_notR_ANCOM.ps)$HabitatAir[rownames(sample_data(ITS_airOnly_notR_ANCOM.ps)) %in% ITSsavAirNames] #double check looks good
ITSForAirNames <- rownames(sample_data(ITS_airOnly_notR_ANCOM.ps))[which(sample_data(ITS_airOnly_notR_ANCOM.ps)$HabitatAir == "forest")]

# 2. FUNGI: differences between habitats
# Need nested random effects of sampling day (DateSetOut is equivalent) within EU.  
sort(colnames(sample_data(ITS_airOnly_notR_ANCOM.ps)))
head(tax_table(ITS_airOnly_notR_ANCOM.ps))
sample_data(ITS_airOnly_notR_ANCOM.ps)$DateSetOut
set.seed(121)
# Grayed out since takes a while, but code works great!
#ITS_habitatAir_nonR_ANCOM <- ancombc2(data = ITS_airOnly_notR_ANCOM.ps, group = "HabitatAir", tax_level = "Species", fix_formula = "HabitatAir", rand_formula = "(1|EU/DateSetOut)", struc_zero = TRUE, neg_lb = FALSE) #FUNGI

# Saved/run October 8, 2025
# (Re-run because Species was not ASV_name, which I need for below)
# saveRDS(ITS_habitatAir_nonR_ANCOM, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_habitatAir_nonR_ANCOM_Oct8-2025.rds")
class(ITS_habitatAir_nonR_ANCOM)

# Look at results:
# Load if need be
ITS_habitatAir_nonR_ANCOM <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_habitatAir_nonR_ANCOM_Oct8-2025.rds")
# Shows that no taxa were differentially abundant but that all passed sensitivity analysis (so pseudo count addition which ANOCMBC2 uses is okay)
unique(ITS_habitatAir_nonR_ANCOM$res$diff_HabitatAirsavanna) #no differentially abundant taxa
unique(ITS_habitatAir_nonR_ANCOM$res$`passed_ss_(Intercept)`) #all passed sensitivity analysis
unique(ITS_habitatAir_nonR_ANCOM$res$passed_ss_HabitatAirsavanna) #all passed sensitivity analysis
# View(ITS_habitatAir_nonR_ANCOM$zero_ind)

###### 1. PULL OUT ONLY ASV TABLE WITH ALL ANCOM ASVS ######
ITS_structZeroForestNames <- ITS_habitatAir_nonR_ANCOM$zero_ind$taxon[which(ITS_habitatAir_nonR_ANCOM$zero_ind$`structural_zero (HabitatAir = forest)` ==TRUE)]
length(ITS_structZeroForestNames) #606
ITS_structZeroSavNames <- ITS_habitatAir_nonR_ANCOM$zero_ind$taxon[which(ITS_habitatAir_nonR_ANCOM$zero_ind$`structural_zero (HabitatAir = savanna)`==TRUE)]
length(ITS_structZeroSavNames) #669

# i. Concatenate all of these ASVs together in order to get all of the ANCOM-identifed ASVs (only struc zero since no main)
allITS_airHab_ANCOM_ASVs <- c(ITS_structZeroForestNames, ITS_structZeroSavNames)
length(allITS_airHab_ANCOM_ASVs) == length(unique(allITS_airHab_ANCOM_ASVs)) #no repeats, as expected
# ii. Get an ASV table with only these ANCOM ASVs
ITS_airHab_ancom_sigASVtab <- ITSair_notR_ASVsTab[rownames(ITSair_notR_ASVsTab) %in% allITS_airHab_ANCOM_ASVs,]
# iii. Transpose so that ASVs are now rows to add in additional ASV info
ITS_airHab_ancom_sigASVtab_t <- as.data.frame(t(ITS_airHab_ancom_sigASVtab))
if (nrow(ITS_airHab_ancom_sigASVtab_t) < ncol(ITS_airHab_ancom_sigASVtab_t)) {
  ITS_airHab_ancom_sigASVtab_t <- t(ITS_airHab_ancom_sigASVtab_t) 
} else {
  print("yay! ASVs are rows")
}
ITS_airHab_ancom_sigASVtab_t <- as.data.frame(ITS_airHab_ancom_sigASVtab_t); class(ITS_airHab_ancom_sigASVtab_t)
setdiff(rownames(ITS_airHab_ancom_sigASVtab_t), allITS_airHab_ANCOM_ASVs) #looks good

# iv. Add in sum in savanna and matrix samples
# Sum of each ASV in air savanna 
ITS_airHab_ancom_sigASVtab_t$savAirTotal <- rowSums(ITS_airHab_ancom_sigASVtab_t[,which(colnames(ITS_airHab_ancom_sigASVtab_t)%in% ITSsavAirNames)])
# Occupancy in air savanna 
ITS_airHab_ancom_sigASVtab_t$savAirOcc <- rowSums(ITS_airHab_ancom_sigASVtab_t[, which(colnames(ITS_airHab_ancom_sigASVtab_t)%in% ITSsavAirNames)] > 0, na.rm = TRUE)
# Total read counts in forest matrix air
ITS_airHab_ancom_sigASVtab_t$forestAirTotal <- rowSums(ITS_airHab_ancom_sigASVtab_t[,which(colnames(ITS_airHab_ancom_sigASVtab_t) %in% ITSForAirNames)])
# Occupancy in forest matrix air
ITS_airHab_ancom_sigASVtab_t$forestAirOcc <- rowSums(ITS_airHab_ancom_sigASVtab_t[, which(colnames(ITS_airHab_ancom_sigASVtab_t)%in% ITSForAirNames)] > 0, na.rm = TRUE)

# v. Add in which kind of ANCOM taxon it is
ITS_airHab_ancom_sigASVtab_t <- ITS_airHab_ancom_sigASVtab_t %>% 
  tibble::rownames_to_column(var= "ASV_name") %>% #make rownames a column
  mutate(ANCOMcat = case_when(
    ASV_name %in% ITS_structZeroForestNames ~ "open patch", #open patch if not in forest
    ASV_name %in% ITS_structZeroSavNames ~ "forested matrix", #forested matrix if not in savanna
  )) %>% 
  mutate(resType = case_when(
    #ASV_name %in% c(bioaerosolITSANCOMASVsName, foliarITSANCOMASVsName) ~ "mainNotStructZero", no main from $res so skip this!
    ASV_name %in% c(ITS_structZeroForestNames, ITS_structZeroSavNames) ~ "structZero"
  )
  )
#View(ITS_airHab_ancom_sigASVtab_t)
# vi. Perform a few checks:
# 1. Get indices for a few ASVs in ITS_airHab_ancom_sigASVtab_t
ITS_HabAirTotal_ASV_50_index <- which(ITS_airHab_ancom_sigASVtab_t$ASV_name == "ASV_50")
ITS_HabAirTotal_ASV_921_index <- which(ITS_airHab_ancom_sigASVtab_t$ASV_name == "ASV_921")
# 2. Is total number working?
colnames(ITS_airHab_ancom_sigASVtab_t)
ITS_airHab_ancom_sigASVtab_t$savAirTotal[ITS_HabAirTotal_ASV_50_index] == sum(ITS_airHab_ancom_sigASVtab_t[ITS_HabAirTotal_ASV_50_index, colnames(ITS_airHab_ancom_sigASVtab_t)%in% ITSsavAirNames])
ITS_airHab_ancom_sigASVtab_t$forestAirTotal[ITS_HabAirTotal_ASV_921_index] == sum(ITS_airHab_ancom_sigASVtab_t[ITS_HabAirTotal_ASV_921_index, colnames(ITS_airHab_ancom_sigASVtab_t)%in% ITSForAirNames])
# 3. Is occupancy working?
ITS_airHab_ancom_sigASVtab_t$forestAirOcc[ITS_HabAirTotal_ASV_50_index] == length(which(ITS_airHab_ancom_sigASVtab_t[ITS_HabAirTotal_ASV_50_index, colnames(ITS_airHab_ancom_sigASVtab_t)%in% ITSForAirNames]>0))
# 4. Are they labeled correctly?
setdiff(ITS_airHab_ancom_sigASVtab_t$ASV_name[which(ITS_airHab_ancom_sigASVtab_t$resType == "structZero")], c(ITS_structZeroForestNames, ITS_structZeroSavNames)) #no differences!
# 5. Are ANCOM categories correct?
setdiff(ITS_airHab_ancom_sigASVtab_t$ASV_name[which(ITS_airHab_ancom_sigASVtab_t$ANCOMcat == "open patch")], ITS_structZeroForestNames) #no differences!
# View(ITS_airHab_ancom_sigASVtab_t)

# vii. Add in taxonomy!
head(airITS_8.5K_taxTab) #USE THIS ONE BECAUSE HAS TAXONOMIC SPECIES NAME NOT "ASV_X"
ITS_airHab_ancom_sigASVtab_t <- left_join(ITS_airHab_ancom_sigASVtab_t, airITS_8.5K_taxTab, by = "ASV_name")

# viii. Add in the proportion of reads in each category that an ASV contributes
# Use this, made above, to get total number of reads in savanna air and forested matrix samples
head(ITSair_notR_ASVsTab)
ITStotReadsSavAir <- sum(unname(colSums(ITSair_notR_ASVsTab[,colnames(ITSair_notR_ASVsTab) %in% ITSsavAirNames]))); ITStotReadsSavAir  #1917987 reads
length(ITSsavAirNames) #51 samples
ITStotReadsForAir <- sum(unname(colSums(ITSair_notR_ASVsTab[,colnames(ITSair_notR_ASVsTab) %in% ITSForAirNames]))); ITStotReadsForAir #2089457
length(ITSForAirNames) #59 samples

# Look at what getting 0.0005 or 0.05% of reads would be about. Do it this way since forest had more reads and more samples
round(0.0005*ITStotReadsSavAir) # needs to occur roughly 959 times in savanna 
round(0.0005*ITStotReadsForAir) # needs to occur roughly 1045 times in forest

# Add in percent that each ASV comprises
ITS_airHab_ancom_sigASVtab_t <- ITS_airHab_ancom_sigASVtab_t %>% 
  mutate(pctSavAir = savAirTotal/ITStotReadsSavAir*100) %>% 
  mutate(pctForAir = forestAirTotal/ITStotReadsForAir *100)

# DECISION AND FILTERING TIME:
# 1. To be considered diff. abundant, an ASV must be in at about 10% of respective group's samples (5 savanna or 6 forest) AND
# 2. Must comprise at least 0.05% of reads in respective group
colnames(ITS_airHab_ancom_sigASVtab_t)
##### AIR SAVANNA ####
round(length(ITSsavAirNames)/10) #5 samples 
# Found in at least 5 samples and .05% of reads
savAir_ITS5occ.05pct_Index <- intersect(intersect(which(ITS_airHab_ancom_sigASVtab_t$savAirOcc >=5), which(ITS_airHab_ancom_sigASVtab_t$pctSavAir >=0.05)), which(ITS_airHab_ancom_sigASVtab_t$ANCOMcat == "open patch"))
length(savAir_ITS5occ.05pct_Index) #0 none detected

##### FOREST AIR ####
round(length(ITSForAirNames)/10) #6 samples
# Found in at least 5 samples and at least 0.01% of reads
forAir_ITS6occ.05pct_Index <- intersect(intersect(which(ITS_airHab_ancom_sigASVtab_t$forestAirOcc >=6), which(ITS_airHab_ancom_sigASVtab_t$pctForAir >=0.05)), which(ITS_airHab_ancom_sigASVtab_t$ANCOMcat == "forested matrix"))
length(forAir_ITS6occ.05pct_Index) #0
