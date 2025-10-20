# sporeSizeSept
# October 16, 2025

# based on version from August 26, 2025
# Just fungal spore size data and figures to edit and review 

#######################################
# I.  SCRIPT SET UP
#######################################
##### 1. LOAD LIBRARIES #####
library(tidyverse)
library(grDevices)
library(dplyr)

# LOAD data
ITS_ANCOMall_.05pct_df <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_ANCOMall_.05pct_df.rds")

########################################################################
#   VI. FUNGAL SPORE SIZE
########################################################################
####### 1. SET UP PART 1 AND INVESTIGATION INTO DATA AVAILABLE #######
head(ITS_ANCOMall_.05pct_df)
# First, make a new column of genus and species for ease later on.
ITS_ANCOMall_.05pct_dfgspMerge <- ITS_ANCOMall_.05pct_df %>% #merge genus and species
  mutate(genusSpecies = paste(Genus, Species, sep = " "))
nrow(ITS_ANCOMall_.05pct_dfgspMerge) - nrow(distinct(ITS_ANCOMall_.05pct_dfgspMerge)) #no duplicate rows

# Which genus and species represent multiple distinct ASVs?
dupGenSpecies <- ITS_ANCOMall_.05pct_df %>%
  count(Genus, Species, sort = TRUE) %>%
  filter(n > 1)
dupGenSpecies #51 genus species have multiple matches!
# View(dupGenSpecies)
# Check a few-- seems to be working. It is (but note for some Genus NA species I can't check as below)
length(which(ITS_ANCOMall_.05pct_df$Species == "carniolicum")) == dupGenSpecies[which(dupGenSpecies$Species == "carniolicum"),3]
length(which(ITS_ANCOMall_.05pct_df$Species == "nigra")) == dupGenSpecies[which(dupGenSpecies$Species == "nigra"),3]

allFungiSpores <-readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/scripts/AllFungi.RDS") #from https://github.com/aguilart/Symbiotic-status-and-fungal-spore-size/blob/main/AllFungi.RDS
head(allFungiSpores)
tail(allFungiSpores)
# Note that there are some taxa with the same volume in the original dataset. For example,
# Aspergillus savannensis, Cordyceps erotyli, Mixtoconidium insidens, Roselliniella euparmeliicola, Trichoderma hausknechtii,
# Trichoderma saturnisporum, Endococcus NA all have a spore volume of 0.22089323 in the original dataset. Same area too: 0.4417865...
# There are other instances of this, too.
# Strange decimal places are most likely due to multiplying by pi to get volume and area.
# View(allFungiSpores) 
# Are there some TOTAL duplicates here? NO! This means that some of the other metadata that's removed later 
# distinguishes rows.
nrow(dplyr::distinct(allFungiSpores)) == nrow(allFungiSpores) #TRUE no duplicates

# Assign spore types for cleaning up data.
unique(allFungiSpores$SporeName)
sexualSpores <- c("Ascospores", "Basidiospores", "Teliospores", "Zygospores")
asexualSpores <- c("Azygospores", "Conidia", "Chlamydospores", "bulbil", "Sporangiospores")
neitherSexualOrAsexual <- c("papulospore")
unique(allFungiSpores$SporeType)
multinucleate <- c("Multinucleate sexual spores", "Multinucleate asexual spores")
uninucleate <- c("Meiospores", "Mitospores")

colnames(allFungiSpores)
# Make a smaller dataset for A.T. et al. dataset
# THIS CAUSES THERE TO BE DUPLICATES, SO THIS SHOULD MEAN THAT ONE OF THE DROPPED COLUMNS DISTINGUISHED ENTRIES.
# STILL, I ONLY CARE ABOUT MATCHING TO GENUS/SPECIES, so I'll just drop all of the duplicates
allFungiSporesTrimmed <- allFungiSpores[,c("SporeName","family","genus", "species","infraspecies", "Specific_sporeName", "codigo",
                                           "names_to_use", "SporeArea","SporeVolume","SporeType")]

# Are there weird instances of where taxonomy differs between species and names_to_use?
table(allFungiSporesTrimmed$species == allFungiSporesTrimmed$names_to_use) #11 false instances

# These seem to be due to having different subspecies or varieties, so can ignore and just used species
# View(allFungiSporesTrimmed[which((allFungiSporesTrimmed$species == allFungiSporesTrimmed$names_to_use) == FALSE),])

# In almost all cases, species and "names_to_use" are the same. Where are they different?
nrow(dplyr::distinct(allFungiSporesTrimmed)) #39216
# Re-name a few of these columns to match the ANCOM dataframe
colnames(allFungiSporesTrimmed)
colnames(allFungiSporesTrimmed)[2] <- "Family"
colnames(allFungiSporesTrimmed)[3] <- "Genus"
colnames(allFungiSporesTrimmed)[4] <- "Species"
head(allFungiSporesTrimmed)

# Why are there some duplicates?
nrow(allFungiSporesTrimmed) - nrow(dplyr::distinct(allFungiSporesTrimmed)) 
# Investigate a duplicate. 
# View(allFungiSpores[which(allFungiSpores$codigo =="Dothistroma septosporum_Conidia_microconidia"),])

# allFungiSporesTrimmed has Genus and Species name in the "species" column. Need to fix for proper merging
allFungiSporesTrimmed$Species
# Separate this out and then check to make sure that there are no discrepancies between "Genus" and this new "Genus2"
allFungiSporesTrimmed <- allFungiSporesTrimmed %>% 
  separate(col = Species, into= c("Genus2", "Species"),
           sep = " ")
# Investigate where these don't match. Just looks like a weird typo, perhaps prompted by a un-recognized character?
allFungiSporesTrimmed[which(allFungiSporesTrimmed$Genus != allFungiSporesTrimmed$Genus2),]
# Simply remove this weird row and Genus2 (which is not needed since it's all the same!)
colnames(allFungiSporesTrimmed) #Genus2 is column 4
allFungiSporesTrimmed <- allFungiSporesTrimmed[-which(allFungiSporesTrimmed$Genus != allFungiSporesTrimmed$Genus2),-4]
allFungiSporesTrimmed[which(allFungiSporesTrimmed$Genus == "Candelariella"),] #Great this shows that the weird one was removed
# ARE THERE DUPLICATE ROWS? -- YES, BUT THIS SHOULD BE FINE, AS LONG AS MEAN IS TAKEN FOR MULTIPLE SPECIES MATCHES
nrow(allFungiSporesTrimmed) - nrow(distinct(allFungiSporesTrimmed)) 

# INVESTIGATE A FEW QUESTIONS, TO DETERMINE IF USING SEXUAL VERSUS ASSEXUAL FORMS, OR MULTINUCLEATE VS. UNINUCLEATE MATTERS
colnames(allFungiSporesTrimmed)
# 1) Explore interactions among SporeType and SporeName
table(allFungiSporesTrimmed$SporeName, allFungiSporesTrimmed$SporeType)
# This shows that all SporeNames (e.g., Conidia, Ascospores, etc.) are only ONE Type of Spore Type
MeisporeSporeNames <- c("Ascospores", "Basidiospores")
MitosporesSporeNames <- c("Chlamydospores", "Conidia", "Sporangiospores")
MultinucleateAsexSporeNames <- "Azygospores"
MultinucleateSexSporeNames <- c("Teliospores", "Zygospores")

# 2) Is it the case that within a given Genus or Species, there are multiple "SporeType" or "SporeName"?
# RESULTS: YES! Essentially, within both species and genus, there are many different sporeTypes and sporeNames.
# View(table(allFungiSporesTrimmed$SporeName, allFungiSporesTrimmed$Genus))
# Find genera with multiple distinct "SporeTypes"
generaMultipleSporeTypes <- allFungiSporesTrimmed %>%
  group_by(Genus) %>%
  summarize(distinct_sporeTypes = n_distinct(SporeType)) %>%
  filter(distinct_sporeTypes > 1)
unique(generaMultipleSporeTypes$Genus) #many, many genera have different sporeTypes (i.e., meispores, mitospores, etc.), which
# ecologically of course makes sense.

# Find genera with multiple distinct "SporeNames"
generaMultipleSporeNames <- allFungiSporesTrimmed %>%
  group_by(Genus) %>%
  summarize(distinct_sporeNames = n_distinct(SporeName)) %>%
  filter(distinct_sporeNames > 1)
unique(generaMultipleSporeNames$Genus) #many, many genera have different sporeNames (e.g., Conidia, Ascospores, etc.)

# Find Species with multiple distinct "SporeTypes"
SpeciesMultipleSporeTypes <- allFungiSporesTrimmed %>%
  group_by(Species) %>%
  summarize(distinct_sporeTypes = n_distinct(SporeType)) %>%
  filter(distinct_sporeTypes > 1)
unique(SpeciesMultipleSporeTypes$Species) #many, many Species have different sporeTypes (i.e., meispores, mitospores, etc.)

# Find Species with multiple distinct "SporeNames"
SpeciesMultipleSporeNames <- allFungiSporesTrimmed %>%
  group_by(Species) %>%
  summarize(distinct_sporeNames = n_distinct(SporeName)) %>%
  filter(distinct_sporeNames > 1)
unique(SpeciesMultipleSporeNames$Species) #many, many Species have different sporeNames (e.g., Conidia, Ascospores, etc.)

# 3) Which SporeTypes and SporeNames have the most data?
sort(table(allFungiSporesTrimmed$SporeName)) #Conidia, followed by Basidiospores, then Ascospores, then Chlamydospores.
sort(table(allFungiSporesTrimmed$SporeType)) # Mitospores and Meiospores (by FAR!), then Multinucleate sexual spores then Multinucleate asexual spores

##### SUMMARY OF DECISION TO FOCUS ON UNINUCLEATE AND SPLIT SEXUAL AND ASEXUAL SPORES #####
# So, all in all, these investigations show that I should likely focus on  Mitospores and Meiospores, which, given that they are 
# classified differently than multinucleate, that they are likely/implied to be "uniNucleate". Given that:
# MeisporeSpores are Ascospores and Basidiospores, and that Mitospores are "Chlamydospores", "Conidia", and "Sporangiospores", these are
# the groups I'll go with. Concerning sexual or asexual, MeisporeSporeNames are sexual,
# and MitosporesSporeNames are asexual. 
MeisporeSporeNames
MitosporesSporeNames

####### 2. SET UP PART 2-- CLEAN UP DATABASE DATAFRAME #######
# So, to tackle this, I will remove all mulitnucleate spores and then divide allFungiSporesTrimmed into sexual and asexual. 
# i. KEEP ONLY "Mitospores" and Meiospores" SPORE TYPE to remove multinucleate spores
unique(allFungiSporesTrimmed$SporeType)
allFungiSporesUniNuc <- allFungiSporesTrimmed[which(allFungiSporesTrimmed$SporeType %in% c("Mitospores", "Meiospores")),]
unique(allFungiSporesUniNuc$SporeType) #"Mitospores" "Meiospores"
nrow(allFungiSporesUniNuc) - nrow(distinct(allFungiSporesUniNuc)) #there are duplicate rows

# ii. DIVIDE allFungiSporesTrimmedUniNuc INTO SEXUAL AND ASEXUAL
# Sexual 
unique(allFungiSporesTrimmed$SporeName)
allFungiSporesUniNucSexual <- allFungiSporesUniNuc[allFungiSporesUniNuc$SporeName %in% MeisporeSporeNames,]
nrow(allFungiSporesUniNucSexual) - nrow(distinct(allFungiSporesUniNucSexual)) #many are repeated, but should be okay with matching steps below
unique(allFungiSporesUniNucSexual$SporeName) # "Ascospores"    "Basidiospores"
# Asexual
allFungiSporesUniNucAsex <- allFungiSporesUniNuc[allFungiSporesUniNuc$SporeName %in% MitosporesSporeNames,]
unique(allFungiSporesUniNucAsex$SporeName) #"Conidia","Chlamydospores","Sporangiospores"

# iii. Check how many duplicates that there are for genus and species
# Which genus and species represent multiple distinct ASVs?-- many many duplicates of this
dupGenSpeciesSex <- allFungiSporesUniNucSexual %>%
  count(Genus, Species, sort = TRUE) %>%
  filter(n > 1)

dupGenSpeciesAsex <- allFungiSporesUniNucAsex %>%
  count(Genus, Species, sort = TRUE) %>%
  filter(n > 1)

# iv. Make a new column with genus and species for ease later on. 
# Sexual 
colnames(allFungiSporesUniNucSexual)
allFungiSporesUniNucSexualgspMerge <- allFungiSporesUniNucSexual %>% #sexual spores
  mutate(genusSpecies = paste(Genus, Species, sep = " "))
colnames(allFungiSporesUniNucSexualgspMerge)
# Asexual 
colnames(allFungiSporesUniNucAsex)
allFungiSporesUniNucAsexualgspMerge <- allFungiSporesUniNucAsex %>%
  mutate(genusSpecies = paste(Genus, Species, sep = " "))
head(allFungiSporesUniNucAsexualgspMerge)

####### 3.SEXUAL SPORES: EXACT GENUS AND SPECIES MATCH #######
# 1. This specifically gets ONLY those ANCOM taxa with an exact match in database, allowing for mutliple
# matches in cases where there are multiple rows of same genus species in either dataframe.
sexSporeExactMatches <- ITS_ANCOMall_.05pct_dfgspMerge %>%
  inner_join(allFungiSporesUniNucSexualgspMerge, by = "genusSpecies", relationship = "many-to-many") %>%
  mutate(match = "exact")
# View(sexSporeExactMatches)
dim(sexSporeExactMatches) #867 198, reflecting many matches
# Remove those with NAs in the species column, since we want only species matches
sexSporeExactMatches <- sexSporeExactMatches[-which(is.na(sexSporeExactMatches$Species.x)==TRUE),]
# View(sexSporeExactMatches)
length(which(is.na(sexSporeExactMatches$Species.x)==TRUE)) #all removed

# 2. INVESTIGATE to make sure it all makes sense. do all of these columns with .x and .y have the same values? (x is ANCOM, y is database)
colnames(sexSporeExactMatches)
colnames(sexSporeExactMatches) %in% c("Family.x", "Family.y")
# i. FAMILY - some variation, but this is okay because of taxonomic fluctuations. Genus and species always are the same
colnames(sexSporeExactMatches)
sexSporeExactMatches[which((sexSporeExactMatches$Family.x == sexSporeExactMatches$Family.y) == FALSE),colnames(sexSporeExactMatches) %in% c("Family.x", "Family.y")]
sexSporeExactMatches[which((sexSporeExactMatches$Family.x == sexSporeExactMatches$Family.y) == FALSE),colnames(sexSporeExactMatches) %in% c("Genus.x", "Genus.y")]
# View(sexSporeExactMatches[which((sexSporeExactMatches$Family.x == sexSporeExactMatches$Family.y) == FALSE),])
# ii. GENUS - none, as expected for these exact matches
unique((sexSporeExactMatches$Genus.x == sexSporeExactMatches$Genus.y)) #TRUE
sexSporeExactMatches[which((sexSporeExactMatches$Genus.x == sexSporeExactMatches$Genus.y) == FALSE),] #no mismatches
# iii. SPECIES - none, as expected for these exact matches
unique(sexSporeExactMatches$Species.x == sexSporeExactMatches$Species.y) #TRUE
sexSporeExactMatches[which((sexSporeExactMatches$Species.x == sexSporeExactMatches$Species.y) == FALSE),]

# 3. CLEAN UP
# So remove un-needed columns for cleaner merging below (all .y all from database, not OG UNITE with ANCOM):
colnames(sexSporeExactMatches) %in% c("Family.y",  "Genus.y", "Species.y", "infraspecies") #list those I want to drop
sexSporeExactMatches_cleaned <- sexSporeExactMatches[,-which(colnames(sexSporeExactMatches) %in% c("Family.y",  "Genus.y", "Species.y", "infraspecies")==TRUE)]
colnames(sexSporeExactMatches_cleaned) #looks good!
# Change names to remove ".x"
colnames(sexSporeExactMatches_cleaned)[colnames(sexSporeExactMatches_cleaned) %in% c("Family.x",  "Genus.x", "Species.x")]
colnames(sexSporeExactMatches_cleaned)[colnames(sexSporeExactMatches_cleaned) %in% c("Family.x",  "Genus.x", "Species.x")] <- c("Family", "Genus", "Species")
sort(colnames(sexSporeExactMatches_cleaned))
# View(sexSporeExactMatches_cleaned)
# I expect more rows than unique ASVs, reflecting instances where one ANCOM ASV had mulitple matches
# 608 matches compared to 68 unique ASVs
length(sexSporeExactMatches_cleaned$ASV_name) > length(unique(sexSporeExactMatches_cleaned$ASV_name))
# REMOVE ROWS with the same exact info
sexSporeExactMatches_cleaned <- distinct(sexSporeExactMatches_cleaned) #remove duplicate rows
dim(sexSporeExactMatches_cleaned)[1] #147 now!
# View(sexSporeExactMatches_cleaned) #this shows that some spores have multiple entries for volume, which I'll deal with later.
length(sexSporeExactMatches_cleaned$ASV_name) #now 147
length(unique(sexSporeExactMatches_cleaned$ASV_name)) #68 ASVs
# This shows that same spore volumes are being assigned to different ASVs, but upon looking, these are different
# ASVs with the same taxonomic info
length(sexSporeExactMatches_cleaned$SporeVolume) == length(unique(sexSporeExactMatches_cleaned$SporeVolume)) 
# For example, Bjerkandera adusta "smoky polypore" corresponds to different ASVs AND has different entries in database
# View(sexSporeExactMatches_cleaned[which(sexSporeExactMatches_cleaned$genusSpecies == "Bjerkandera adusta"),])
# ASV_447 and ASV_410
# And 4 different spore volumes, likely because of different sources in the database, which is fine.
unique(sexSporeExactMatches_cleaned$SporeVolume[which(sexSporeExactMatches_cleaned$genusSpecies == "Bjerkandera adusta")])
mean(unique(sexSporeExactMatches_cleaned$SporeVolume[which(sexSporeExactMatches_cleaned$genusSpecies == "Bjerkandera adusta")]))
median(unique(sexSporeExactMatches_cleaned$SporeVolume[which(sexSporeExactMatches_cleaned$genusSpecies == "Bjerkandera adusta")]))

# 4. SET UP GENUS MATCHING IN NEXT STEP:
# Remove exact matches from ITS_ANCOMall_.05pct_dfgspMerge to match with genus (SO USE THIS IN NEXT STEPS)
ITS_ANCOMall_.05pct_dfgspMerge_genusOnly <- ITS_ANCOMall_.05pct_dfgspMerge[which(ITS_ANCOMall_.05pct_dfgspMerge$ASV_name %in% sexSporeExactMatches_cleaned$ASV_name == FALSE),]
unique(ITS_ANCOMall_.05pct_dfgspMerge_genusOnly$ASV_name %in% sexSporeExactMatches_cleaned$ASV_name) # FALSE so, these were removed
# Check if no duplicate rows after doing this-- CONFIRMED!
nrow(ITS_ANCOMall_.05pct_dfgspMerge_genusOnly) - nrow(distinct(ITS_ANCOMall_.05pct_dfgspMerge_genusOnly))  #0

####### 4.SEXUAL SPORES: GENUS LEVEL MATCH #######
# 1. GENUS LEVEL MATCH and then clean up this dataframe using DF THAT DOES NOT HAVE SPECIES MATCHES
nrow(ITS_ANCOMall_.05pct_dfgspMerge_genusOnly) - nrow(distinct(ITS_ANCOMall_.05pct_dfgspMerge_genusOnly)) # NO DUPLICATES HERE
# With the join below, only taxa that they both share, ancom dataset and spore size dataset
sexSporeGenusMatches <- ITS_ANCOMall_.05pct_dfgspMerge_genusOnly %>%
  inner_join(allFungiSporesUniNucSexualgspMerge, by = "Genus",relationship = "many-to-many") %>%
  mutate(match = "genusMatch")
# View(sexSporeGenusMatches)
colnames(sexSporeGenusMatches)

# 2. INVESTIGATE to make sure it all makes sense. ANd this is fine given that genus names have to be unique within family by naming convention
colnames(sexSporeGenusMatches)
# i. FAMILY - some do vary but since species and genus are the same, not going to worry about in-flux taxonomic assignments
sexSporeGenusMatches[which((sexSporeGenusMatches$Family.x == sexSporeGenusMatches$Family.y) == FALSE),colnames(sexSporeGenusMatches) %in% c("Family.x", "Family.y")]
# View(sexSporeGenusMatches[which((sexSporeGenusMatches$Family.x == sexSporeGenusMatches$Family.y) == FALSE),])
# 49 distinct ASVs where this happens examples (double check that these aren't totally different ). Get ASV name, family x and family y and genus
famGenMatchDifferences <- distinct(sexSporeGenusMatches[which((sexSporeGenusMatches$Family.x == sexSporeGenusMatches$Family.y) == FALSE),colnames(sexSporeGenusMatches) %in% c("ASV_name", "Family.x", "Family.y", "Genus")])
# View(famGenMatchDifferences)
# Check to make sure that these 20 genera are the same (i.e., that different families don't have v different
# genera names...)
unique(famGenMatchDifferences$Genus) #19 unique genera. 

# ii. SPECIES - expect all different (unless both unassigned )
sexSporeGenusMatches[which((sexSporeGenusMatches$Species.x == sexSporeGenusMatches$Species.y) == FALSE),]
sexSporeGenusMatches[which((sexSporeGenusMatches$Species.x == sexSporeGenusMatches$Species.y) == TRUE),] #in no cases is this is same, as expected!
# Only one genus column so these are all the same
# So remove un-needed columns for cleaner merging below:
colnames(sexSporeGenusMatches) %in% c("Family.y", "Species.y", "infraspecies") #list those I want to drop
sexSporeGenusMatches_cleaned <- sexSporeGenusMatches[,-which(colnames(sexSporeGenusMatches) %in% c("Family.y", "Species.y", "infraspecies")==TRUE)]
colnames(sexSporeGenusMatches_cleaned) #looks good!
# Change names to remove ".x"
colnames(sexSporeGenusMatches_cleaned)[colnames(sexSporeGenusMatches_cleaned) %in% c("Family.x", "Species.x")]
colnames(sexSporeGenusMatches_cleaned)[colnames(sexSporeGenusMatches_cleaned) %in% c("Family.x", "Species.x")] <- c("Family", "Species")
# Keep ANCOM genusSpecies.x as genusSpecies and rename the other one as "genusSpeciesDatabase"
colnames(sexSporeGenusMatches_cleaned)[colnames(sexSporeGenusMatches_cleaned) %in% c("genusSpecies.x", "genusSpecies.y")] <- c("genusSpecies", "genusSpeciesDatabase")
sort(colnames(sexSporeGenusMatches_cleaned))
# View(sexSporeGenusMatches_cleaned)
# Are there any duplicate rows?
nrow(sexSporeGenusMatches_cleaned) # 5760
nrow(distinct(sexSporeGenusMatches_cleaned)) #2453
# Remove totally duplicate rows
sexSporeGenusMatches_cleaned <- distinct(sexSporeGenusMatches_cleaned)
nrow(sexSporeGenusMatches_cleaned) #2453

# 3. SET UP NO MATCH AVAILABLE  Assign still unmatched rows as "noMatchAvailable". ITS_ANCOMall_.05pct_dfgspMerge_genusOnly has all without species match,sexSporeGenusMatches_cleaned is just genus.
# So, this line gets those ASVs for which there are no genus or species matches
ITS_ANCOMall_.05pct_dfgspMerge_noMatches <- ITS_ANCOMall_.05pct_dfgspMerge_genusOnly[which(ITS_ANCOMall_.05pct_dfgspMerge_genusOnly$ASV_name %in% sexSporeGenusMatches_cleaned$ASV_name == FALSE),]
unique(ITS_ANCOMall_.05pct_dfgspMerge_noMatches$ASV_name %in%sexSporeGenusMatches_cleaned$ASV_name) #FALSE great, these were removed
unique(ITS_ANCOMall_.05pct_dfgspMerge_noMatches$ASV_name %in%sexSporeExactMatches_cleaned$ASV_name) #FALSE great, these were removed
ITS_ANCOMall_.05pct_dfgspMerge_noMatches <- ITS_ANCOMall_.05pct_dfgspMerge_noMatches %>% 
  mutate(
    size = "noDatabaseMatch",
    match = "noMatchAvailable"
  )
sort(colnames(ITS_ANCOMall_.05pct_dfgspMerge_noMatches))
# NO duplicates here, as expected (below returns TRUE)
nrow(ITS_ANCOMall_.05pct_dfgspMerge_noMatches) == nrow(distinct(ITS_ANCOMall_.05pct_dfgspMerge_noMatches))

# 4. CHECK NO OVERLAP IN SPECIES MATCHES, GENUS MATCHES, AND NON MATCHED
# Check no overlap in ASVs to insure that ASVs aren't being double counted
intersect(ITS_ANCOMall_.05pct_dfgspMerge_noMatches$ASV_name, sexSporeGenusMatches_cleaned$ASV_name) #no overlap
intersect(sexSporeGenusMatches_cleaned$ASV_name, sexSporeExactMatches_cleaned$ASV_name) #no overlap 
colnames(ITS_ANCOMall_.05pct_dfgspMerge_noMatches); colnames(sexSporeGenusMatches_cleaned); colnames(sexSporeExactMatches_cleaned)
dim(ITS_ANCOMall_.05pct_dfgspMerge_noMatches); dim(sexSporeGenusMatches); dim(sexSporeExactMatches)
# View(ITS_ANCOMall_.05pct_dfgspMerge_noMatches)

# 5. CHECK SPECIFIC SPORE TYPES
table(sexSporeGenusMatches_cleaned$Specific_sporeName, sexSporeGenusMatches_cleaned$ANCOMcat)
# CURRENT: Ascos only in foliar, basidio on in bioaerosol. Matches ANCOM taxonomy!
#                 bioaerosol foliar surface
# ascospores             0            472
# basidiospores       1940             41

####### 4.SUMMARIZING SEXUAL SPORES BY FULL (sexSporeExactMatches_cleaned) AND GENUS (sexSporeGenusMatches_cleaned) MATCHES #######
# (this is where diverges most substantially from ANCOM_traits_I6SandITS.R)
# 1. INFER Mean, median, range, etc. per ASV for sexSporeExactMatches_cleaned
sort(colnames(sexSporeExactMatches_cleaned))
sexualExact_sum <- sexSporeExactMatches_cleaned %>%
  group_by(ASV_name, Phylum, Class, Order, Family, Genus, Species, ANCOMcat, bioaerosolOcc, foliarSurfaceOcc, pctBioaerosol, pctFoliar) %>%
  # group_by all of of these to retain rows, and is fine since each ASV matches in this respect
  # should retain ASVs that have the same species designation, e.g. ASV_447 and ASV_410, both B. adjusta with 4 separate sizes)
  summarize(
    numberOfMatches = n_distinct(SporeVolume), #this can be that it is in sexSporeExactMatches_cleaned, sexSporeGenusMatches_cleaned, ITS_ANCOMall_.05pct_dfgspMerge_noMatches. So all should have this
    meanVol = mean(SporeVolume, na.rm = TRUE),
    medianVol = median(SporeVolume, na.rm = TRUE), 
    minVol = min(SporeVolume, na.rm = TRUE),
    maxVol = max(SporeVolume, na.rm = TRUE),
    rangeVol = max(SporeVolume, na.rm = TRUE) - min(SporeVolume, na.rm = TRUE),
    .groups = "drop"
  )
# View(sexualExact_sum)
length(unique(sexSporeExactMatches_cleaned$ASV_name)) == nrow(sexualExact_sum); nrow(sexualExact_sum)  #68, matches before as expected!
sort(unique(sexualExact_sum$ASV_name)) == sort(unique(sexSporeExactMatches_cleaned$ASV_name)) #these are the same!
which(is.nan(sexualExact_sum$meanVol)) #all have amounts

# 2. INFER Mean, median, range, etc. per ASV for GENUS LEVEL MATCHES (sexSporeGenusMatches_cleaned)
sort(colnames(sexSporeGenusMatches_cleaned))
sexualGenLevel_sum <- sexSporeGenusMatches_cleaned %>%
  group_by(ASV_name, Phylum, Class, Order, Family, Genus, ANCOMcat, bioaerosolOcc, foliarSurfaceOcc, pctBioaerosol, pctFoliar) %>%
  # Do not group by Species
  # should retain ASVs that have the same species designation, e.g. ASV_447 and ASV_410, both B. adjusta with 4 separate sizes)
  summarize(
    numberOfMatches = n_distinct(SporeVolume), #this can be that it is in sexSporeGenusMatches_cleaned, sexSporeGenusMatches_cleaned, ITS_ANCOMall_.05pct_dfgspMerge_noMatches. So all should have this
    meanVol = mean(SporeVolume, na.rm = TRUE),
    medianVol = median(SporeVolume, na.rm = TRUE), 
    minVol = min(SporeVolume, na.rm = TRUE),
    maxVol = max(SporeVolume, na.rm = TRUE),
    rangeVol = max(SporeVolume, na.rm = TRUE) - min(SporeVolume, na.rm = TRUE),
    .groups = "drop"
  )
# View(sexualGenLevel_sum)
length(unique(sexSporeGenusMatches_cleaned$ASV_name)) == nrow(sexualGenLevel_sum) #124, matches before as expected!
unique(sort(unique(unique(sexualGenLevel_sum$ASV_name)) == sort(unique(sexSporeGenusMatches_cleaned$ASV_name)))) #these are the same!
which(is.nan(sexualGenLevel_sum$meanVol)) #all have amounts

# 3. INVESTIGATE Double check that these are working as expected -- yay! 
# i. (these show that each genera has same mean volume, even though species (UNITE database) differs and is still shown in df)
# Genus == "Hyphodontia"
median(sexSporeGenusMatches_cleaned$SporeVolume[which(sexSporeGenusMatches_cleaned$Genus == "Hyphodontia")]) == sexualGenLevel_sum$medianVol[which(sexualGenLevel_sum$Genus == "Hyphodontia")]
sexualGenLevel_sum$numberOfMatches[which(sexualGenLevel_sum$Genus == "Hyphodontia")] #47 matches wow
# Genus == "Trametes"
median(sexSporeGenusMatches_cleaned$SporeVolume[which(sexSporeGenusMatches_cleaned$Genus == "Trametes")]) == sexualGenLevel_sum$medianVol[which(sexualGenLevel_sum$Genus == "Trametes")]
sexualGenLevel_sum$numberOfMatches[which(sexualGenLevel_sum$Genus == "Trametes")] #65 matches wow
# ii. Does each ASV have same mean spore volume across the multiple matches per ASV?
sexualGenLevel_sumCHECK <- sexualGenLevel_sum %>% 
  group_by(ASV_name) %>% 
  summarize(nVolsPerASV = n_distinct(medianVol))
unique(sexualGenLevel_sumCHECK$nVolsPerASV) #all 1, confirms that matching as expected.
# iii. Number of matches expected?
nrow(sexualGenLevel_sum) == nrow(distinct(sexualGenLevel_sum)); nrow(sexualGenLevel_sum) #yes! 124

# 4. KEEP ONLY MATCHES (exact then genus matches) WHERE MATCH VARIATION RELATIVE TO THE MEDIAN IS SMALL
# In edge effects paper, "retained the ASV only if mean calculated genome size was less than 
# 15% of the range between the smallest and the largest matched genomes.
# ASV only if the matched genome sizes were similar: ((max−min)/median < 0.15)

# i. Exact (species-level) matches
head(sexualExact_sum)
colnames(sexualExact_sum)
sexExact_trim <-  sexualExact_sum %>% 
  mutate(rangeDivMed = rangeVol/medianVol) %>% #keep only spores where matches are close together, relative to the median
  mutate(keep = rangeVol < 0.15 * medianVol) 
table(sexExact_trim$keep, sexExact_trim$ANCOMcat) #CURRENT: 25 bioaerosol and 6 foliar. 
sexExact_trim2 <- sexExact_trim %>% #Filter out ones not to keep
  filter(keep == TRUE)
dim(sexExact_trim2)[1] #31, as expected (25 bioaerosol + 6 foliar)

# ii. Genus-level matches
head(sexualGenLevel_sum)
colnames(sexualGenLevel_sum)
sexGen_trim <-  sexualGenLevel_sum %>% 
  mutate(rangeDivMed = rangeVol/medianVol) %>% #keep only spores where matches are close together, relative to the median
  mutate(keep = rangeVol < 0.15 * medianVol) 
# View(sexGen_trim)
table(sexGen_trim$keep, sexGen_trim$ANCOMcat) #CURRENT: total 15. Only 0 bioaerosol and 15 foliar. 
sexGen_trim2 <- sexGen_trim %>% #filter out ones not to keep
  filter(keep == TRUE)
dim(sexGen_trim2)[1] #15, as expected

####### 5. JOIN GENUS AND SPECIES MATCHES INTO ONE DATAFRAME TO USE LATER ON! ######
# 1. CHECK Next 2 lines show only diff is sexualExact_sum has Species column
setdiff(colnames(sexualExact_sum), colnames(sexualGenLevel_sum))
setdiff(colnames(sexualGenLevel_sum), colnames(sexualExact_sum))
# 2. ADD DETAILS BEFORE MERGING
sexualGenLevel_sum$Species <- NA
sexualGenLevel_sum$matchType <- "genusMatch"
sexualExact_sum$matchType <- "speciesMatch"
# 3. ROW BIND!
sexSporeMatches <- rbind(sexualExact_sum,sexualGenLevel_sum)
# Saved October 16, 2025
# saveRDS(sexSporeMatches, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/sexSporeMatches_October16.RData")
# 4. EXPLORE AND THEN GET MEDIANS FOR MANUSCRIPT
colnames(sexSporeMatches) %in% "medianVol"
median(sexSporeMatches$medianVol[sexSporeMatches$ANCOMcat == "foliar surface"]) #175.9292 (175.9 reported in MS)
median(sexSporeMatches$medianVol[sexSporeMatches$ANCOMcat == "bioaerosol"]) #21.49005 (21.5 reported in MS)

####### 6. RESTRICTED BY VARIATION: JOIN GENUS AND SPECIES MATCHES INTO ONE DATAFRAME TO USE LATER ON! ######
# 1. CHECK Next 2 lines show only diff is sexExact_trim2 has Species column
setdiff(colnames(sexExact_trim2), colnames(sexGen_trim2))
setdiff(colnames(sexGen_trim2), colnames(sexExact_trim2))
# 2. ADD DETAILS BEFORE MERGING
sexGen_trim2$Species <- NA
sexGen_trim2$matchType <- "genusMatch"
sexExact_trim2$matchType <- "speciesMatch"
# 3. ROW BIND!
sexSporeMatches_trim <- rbind(sexExact_trim2,sexGen_trim2)
# View(sexSporeMatches_trim)
# saveRDS(sexSporeMatches_trim, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/sexSporeMatches_trim_Oct16.RData")
# 4. EXPLORE AND THEN GET MEDIANS
colnames(sexSporeMatches_trim) %in% "medianVol"
median(sexSporeMatches_trim$medianVol[sexSporeMatches_trim$ANCOMcat == "foliar surface"]) #154.8543
median(sexSporeMatches_trim$medianVol[sexSporeMatches_trim$ANCOMcat == "bioaerosol"]) #22.38385 
table(sexSporeMatches_trim$ANCOMcat)
# Number of matches in each group
# bioaerosol foliar surface 
# 25             21 

####### 7. VIOLIN PLOT ######
# 1. Genus and species matches
head(sexSporeMatches)
unique(sexSporeMatches$Order) #18 unique
log(max(sexSporeMatches$medianVol)) #max is 7.007117, informative for ylim parameter
# View(sexSporeMatches)
dim(sexSporeMatches)[1] #192
sexSpore_ViolinPlot <- ggplot(data=sexSporeMatches, aes(x=ANCOMcat, y=log(medianVol))) + 
  geom_violin(linewidth = 1.25) + # (aes(color= ANCOMcat), linewidth =1.5) this would change outline
  geom_jitter(aes(color = ANCOMcat), size=2, alpha=0.6, height = 0, width = 0.35) +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.text.y = element_text(size = 10, color= "black"),
    axis.title.y = element_text(size = 10), color= "black",
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  )+
  scale_color_manual(values = c("bioaerosol" = "#004CFF", "foliar surface" = "darkgreen")) + 
  labs(y = "Log of Median Spore Volume") +
  ggtitle("Fungal Uninucleate Sexual Spores") +
  ylim(-1, 10)  #Change y-axis so that there is enough room on top to add mean spore sizes
sexSpore_ViolinPlot 
# Saved October 16, 2025
# saveRDS(sexSpore_ViolinPlot, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/sexSpore_ViolinPlot_Oct16.rds")

# 2. Just species match
dim(sexualExact_sum) #68 rows 
sexSpore_JustExact_ViolinPlot <- ggplot(data=sexualExact_sum, aes(x=ANCOMcat, y=log(medianVol))) + 
  geom_violin(linewidth = 1.25) + # (aes(color= ANCOMcat), linewidth =1.5) this would change outline
  geom_jitter(aes(color = ANCOMcat), size=4, alpha=0.6, height = 0, width = 0.35) +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position="none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16)) +
  scale_color_manual(values = c("bioaerosol" = "#004CFF", "foliar surface" = "darkgreen")) + 
  labs(y = "Log of Median Spore Volume") +
  ggtitle("Fungal Uninucleate Sexual Spores") +
  ylim(-1, 10)  #Change y-axis so that there is enough room on top to add mean spore sizes
sexSpore_JustExact_ViolinPlot # not saved since this is not the verson for the script

# 3. With new, restricted dataset (keep only those with max-min vol divided by median < 0.15)
sexSporeMatches_trim
dim(sexSporeMatches_trim) #46 ASVs
sexSpore_TrimViolPlot <- ggplot(data=sexSporeMatches_trim, aes(x=ANCOMcat, y=log(medianVol))) + 
  geom_violin(linewidth = 1.25) + # (aes(color= ANCOMcat), linewidth =1.5) this would change outline
  geom_jitter(aes(color = ANCOMcat), size=4, alpha=0.6, height = 0, width = 0.35) +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position="none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16)) +
  scale_color_manual(values = c("bioaerosol" = "#004CFF", "foliar surface" = "darkgreen")) + 
  #scale_x_discrete(labels = c("air" = "air", "phyllo" = "leaf surfaces")) +  
  labs(y = "Log of Median Spore Volume") +
  ggtitle("Fungal Uninucleate Sexual Spores") +
  ylim(-1, 10)  #Change y-axis so that there is enough room on top to add mean spore sizes
sexSpore_TrimViolPlot # not saved since not in manuscript

####### 8. STATS ######
which(is.na(sexSporeMatches$medianVol)) #no NAs in data
# Check to see if volume is normally distributed with Shapiro-Wilk test
shapiro.test(x = sexSporeMatches$medianVol) 
# Shapiro-Wilk normality test
# data:  sexSporeMatches$medianVol
# W = 0.55522, p-value < 2.2e-16
# (reject null that is normal)
# Does log transforming these make them normal (this also matches the likely plots...). NO!
# Q-Q plot - looks pretty good with plots.
qqnorm(log(sexSporeMatches$medianVol))
qqline(log(sexSporeMatches$medianVol))
shapiro.test(x = log(sexSporeMatches$medianVol)) #W = 0.9594, p-value = 2.513e-05 these data are still NOT normal (reject null that is normal)
# Given that these are non-normal, perform non-parametric test (Wilcoxon rank sum test/Mann-Whitney test)
# Test, do bioaerosols have smaller spores?
set.seed(19)
sexSporeWilcoxon <- wilcox.test(medianVol ~ ANCOMcat, data = sexSporeMatches, alternative = "less") #less because biaoerosol first alphabeticaly 
sexSporeWilcoxon # Wilcoxon rank sum test with continuity correction
# data:  medianVol by ANCOMcat
# W = 818, p-value < 2.2e-16
# alternative hypothesis: true location shift is less than 0 SHOWS THAT PHYLLOSPHERE DOES HAVE BIGGER SPORES THAN AIR
# (AMONG SEXUAL SPORES ANYWAY) # added to manuscript October 5, 2025

# 2. How many of the original foliar and air ANCOM taxa could I get info for? (added to MS Aug. 28. 2025)
sexSporeMatches #matches with database
ITS_ANCOMall_.05pct_df #original ANCOM

table(sexSporeMatches$ANCOMcat) #116 air and 76 phyllo
table(ITS_ANCOMall_.05pct_df$ANCOMcat) #152 air and 243 phyllo
116/152*100 #air is 76.31579, or 76.3%, reported in MS
76/243*100 #foliar is 31.27572, or 31.3%, reported in MS Oct. 5
table(sexSporeMatches$Order, sexSporeMatches$ANCOMcat)

# 2. Just exact, species match
colnames(sexualExact_sum)
sexualExact_sum$medianVol
which(is.na(sexualExact_sum$medianVol)) #no NAs in data
# Check to see if volume is normally distributed with Shapiro-Wilk test
shapiro.test(x = sexualExact_sum$medianVol)# W = 0.55175, p-value = 4.297e-13 these data are NOT normal (reject null that is normal)
# Does log transforming these make them normal (this also matches the likely plots...). NO!
# Q-Q plot - not quite
qqnorm(log(sexualExact_sum$medianVol))
qqline(log(sexualExact_sum$medianVol))
shapiro.test(x = log(sexualExact_sum$medianVol)) #W = 0.95292, p-value = 0.01195 these data are not normal (reject null that is normal)

# Non-parametric test (Wilcoxon rank sum test/Mann-Whitney test)
sexualSporeWilcoxon_exact <- wilcox.test(medianVol ~ ANCOMcat, data = sexualExact_sum, alternative = "less")
sexualSporeWilcoxon_exact #W = 81, p-value = 0.0004329 Shows that bioaersol is smaller than foliar surface

# 3. With new, more restrictive matches
colnames(sexSporeMatches_trim)
sexSporeMatches_trim$medianVol
which(is.na(sexSporeMatches_trim$medianVol)) #no NAs in data
# Check to see if volume is normally distributed with Shapiro-Wilk test
shapiro.test(x = sexSporeMatches_trim$medianVol) #W = 0.73862, p-value = 1.076e-07 these data are NOT normal (reject null that is normal)
# Does log transforming these make them normal (this also matches the likely plots...). NO!
# Q-Q plot - not quite
qqnorm(log(sexSporeMatches_trim$medianVol))
qqline(log(sexSporeMatches_trim$medianVol))
shapiro.test(x = log(sexSporeMatches_trim$medianVol)) #W = 0.93979, p-value = 0.01922 no

# Perform non-parametric test (Wilcoxon rank sum test/Mann-Whitney test)
sexualSporeWilcoxon_trim <- wilcox.test(medianVol ~ ANCOMcat, data = sexSporeMatches_trim, alternative = "less")
sexualSporeWilcoxon_trim #W = 142, p-value = 0.004016 Shows that there is difference among groups, air is 
# smaller

####### 9. INFO BY TAXONOMIC GROUP ######
# i. In full spore volume datasets 
unique(sexSporeMatches$Order) #18 now
ordersANCOM <- as.data.frame.matrix(table(sexSporeMatches$Order, sexSporeMatches$ANCOMcat))
ordersANCOM <- rownames_to_column(ordersANCOM)
colnames(ordersANCOM) <- c("Order", "bioaerosol_all", "foliar_all")
ordersANCOM
# Get median volume for each order
orderMed <- sexSporeMatches %>% 
  group_by(Order) %>% 
  summarize(orderMed = median(medianVol))
ordersANCOM <- merge(ordersANCOM, orderMed, by="Order")
ordersANCOM

# ii. In subsetted spore dataset
unique(sexSporeMatches_trim$Order) #11
trimOrdersANCOM <- as.data.frame.matrix(table(sexSporeMatches_trim$Order, sexSporeMatches_trim$ANCOMcat))
trimOrdersANCOM <- rownames_to_column(trimOrdersANCOM)
colnames(trimOrdersANCOM) <- c("Order", "bioaerosol_trim", "foliar_trim")
trimOrdersANCOM
trimOrderMed <- sexSporeMatches_trim %>% 
  group_by(Order) %>% 
  summarize(orderMed_trim = median(medianVol))
trimOrdersANCOM <- merge(trimOrdersANCOM, trimOrderMed, by="Order")
trimOrdersANCOM

# iii. Combine dataframes together
orderSporeSize_all <- left_join(ordersANCOM, trimOrdersANCOM, by = "Order")
# View(orderSporeSize_all) #note that NAs in trimmed part should be zeros
# write.csv(orderSporeSize_all, file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/orderSporeSize_all_October16.csv")

########## ASEXUAL SPORES ##########
# NOTE SEPT. 3, 2025: NOT USED IN MANUSCRIPT
####### 3.ASEXUAL SPORES: EXACT GENUS AND SPECIES MATCH #######
# 1. This specifically gets ONLY those ANCOM taxa with an exact match in database, allowing for mutliple
# matches in cases where there are multiple rows of same genus species in either dataframe.
asexSporeExactMatches <- ITS_ANCOMall_.05pct_dfgspMerge %>%
  inner_join(allFungiSporesUniNucAsexualgspMerge, by = "genusSpecies", relationship = "many-to-many") %>%
  mutate(match = "exact")
# View(asexSporeExactMatches)
dim(asexSporeExactMatches) #843  25, reflecting many matches
# Remove those with NAs in the species column, since we want only species matches (genusSpecies above has some without species)
asexSporeExactMatches <- asexSporeExactMatches[-which(is.na(asexSporeExactMatches$Species.x)==TRUE),]
# View(asexSporeExactMatches)
length(which(is.na(asexSporeExactMatches$Species.x)==TRUE)) #all removed

# 2. INVESTIGATE to make sure it all makes sense. do all of these columns with .x and .y have the same values? (x is ANCOM, y is database)
colnames(asexSporeExactMatches)
# i. FAMILY - some variation, but this is okay because of taxonomic fluctuations
asexSporeExactMatches[which((asexSporeExactMatches$Family.x == asexSporeExactMatches$Family.y) == FALSE),c(6,15)]
# View(asexSporeExactMatches[which((asexSporeExactMatches$Family.x == asexSporeExactMatches$Family.y) == FALSE),])
# ii. GENUS - none, as expected for these exact matches
unique((asexSporeExactMatches$Genus.x == asexSporeExactMatches$Genus.y)) #TRUE
asexSporeExactMatches[which((asexSporeExactMatches$Genus.x == asexSporeExactMatches$Genus.y) == FALSE),] #no mismatches
# iii. SPECIES - none, as expected for these exact matches
unique(asexSporeExactMatches$Species.x == asexSporeExactMatches$Species.y) #TRUE
asexSporeExactMatches[which((asexSporeExactMatches$Species.x == asexSporeExactMatches$Species.y) == FALSE),]

# 3. CLEAN UP
# So remove un-needed columns for cleaner merging below (all .y all from database, not OG UNITE with ANCOM):
colnames(asexSporeExactMatches) %in% c("Family.y",  "Genus.y", "Species.y", "infraspecies") #list those I want to drop
asexSporeExactMatches_cleaned <- asexSporeExactMatches[,-which(colnames(asexSporeExactMatches) %in% c("Family.y",  "Genus.y", "Species.y", "infraspecies")==TRUE)]
colnames(asexSporeExactMatches_cleaned) #looks good!
# Change names to remove ".x"
colnames(asexSporeExactMatches_cleaned)[colnames(asexSporeExactMatches_cleaned) %in% c("Family.x",  "Genus.x", "Species.x")] <- c("Family", "Genus", "Species")
sort(colnames(asexSporeExactMatches_cleaned))
# View(asexSporeExactMatches_cleaned)
# I expect more rows than unique ASVs, reflecting instances where one ANCOM ASV had mulitple matches
# 210 matches compared to 43 unique ASVs
length(asexSporeExactMatches_cleaned$ASV_name) > length(unique(asexSporeExactMatches_cleaned$ASV_name)); length(asexSporeExactMatches_cleaned$ASV_name); length(unique(asexSporeExactMatches_cleaned$ASV_name))
# REMOVE ROWS with the same exact info
asexSporeExactMatches_cleaned <- distinct(asexSporeExactMatches_cleaned) #remove duplicate rows
dim(asexSporeExactMatches_cleaned) #67 now!
# View(asexSporeExactMatches_cleaned) #this shows that some spores have multiple entries for volume, which I'll deal with later.
length(asexSporeExactMatches_cleaned$ASV_name) #now 67
length(unique(asexSporeExactMatches_cleaned$ASV_name)) #43 ASVs (same as above)
# This shows that same spore volumes are being assigned to different ASVs, but upon looking, these are different
# ASVs with the same taxonomic info
length(asexSporeExactMatches_cleaned$SporeVolume) == length(unique(asexSporeExactMatches_cleaned$SporeVolume)) 

# 4. SET UP GENUS MATCHING IN NEXT STEP:
# Remove exact matches from ITS_ANCOMall_.05pct_dfgspMerge to match with genus (SO USE THIS IN NEXT STEPS)
# (keeps only those ASVs which were not represented in exact matches above)
airLeafANCOM_genusOnlyASex <- ITS_ANCOMall_.05pct_dfgspMerge[which(ITS_ANCOMall_.05pct_dfgspMerge$ASV_name %in% asexSporeExactMatches_cleaned$ASV_name == FALSE),]
unique(airLeafANCOM_genusOnlyASex$ASV_name %in% asexSporeExactMatches_cleaned$ASV_name) #FALSE, so these were removed
# Check if no duplicate rows after doing this-- CONFIRMED!
nrow(airLeafANCOM_genusOnlyASex) - nrow(distinct(airLeafANCOM_genusOnlyASex)) 

# 5. CHECK SPECIFIC SPORE TYPES -- this may explain part of issue: Chlamydospores are resistant,
# usually just sit in soil, with wiki saying no dispersal mechanism. 
table(asexSporeExactMatches_cleaned$Specific_sporeName, asexSporeExactMatches_cleaned$ANCOMcat)
#               bioaerosol foliar surface
# arthroconidia            6              0
# asteroconidia            0              1
# ballistoconidia          0              3
# blastoconidia            0              1
# chlamydospores          23              2
# conidia                 18             59
# ramoconidia              0              6
table(asexSporeExactMatches_cleaned$SporeName, asexSporeExactMatches_cleaned$ANCOMcat)
#                 bioaerosol foliar surface
# Chlamydospores         23              2
# Conidia                24             70
####### 4.ASEXUAL SPORES: GENUS LEVEL MATCH #######
# 1. GENUS LEVEL MATCH and then clean up this dataframe using DF THAT DOES NOT HAVE SPECIES MATCHES
nrow(airLeafANCOM_genusOnlyASex) - nrow(distinct(airLeafANCOM_genusOnlyASex)) # NO DUPLICATES HERE
# With the join below, only taxa that they both share, ancom dataset and spore size dataset
asexSporeGenusMatches <- airLeafANCOM_genusOnlyASex %>%
  inner_join(allFungiSporesUniNucAsexualgspMerge, by = "Genus",relationship = "many-to-many") %>%
  mutate(match = "genusMatch")
# View(asexSporeGenusMatches)
table(asexSporeGenusMatches$SporeName, asexSporeGenusMatches$ANCOMcat)
                  # bioaerosol foliar surface
# Chlamydospores        1181            389
# Conidia               2156           8722
# Sporangiospores          0             20

# 2. INVESTIGATE to make sure it all makes sense.
colnames(asexSporeGenusMatches)
# i. FAMILY - some do vary but since species and genus are the same, not going to worry about in-flux taxonomic assignments
unique(asexSporeGenusMatches[which((asexSporeGenusMatches$Family.x == asexSporeGenusMatches$Family.y) == FALSE),c(9,15)])
# View(asexSporeGenusMatches[which((asexSporeGenusMatches$Family.x == asexSporeGenusMatches$Family.y) == FALSE),])
# Only 46 distinct ASVs where this happens examples (double check that these aren't totally different ). Get ASV name, family x and family y and genus
famGenMatchDifferences_Asex <- distinct(asexSporeGenusMatches[which((asexSporeGenusMatches$Family.x == asexSporeGenusMatches$Family.y) == FALSE),c(1,9,10, 15)])
# View(famGenMatchDifferences_Asex)
# Check to make sure that these 25 genera are the same (i.e., that different families don't have v different
# genera names...)
unique(famGenMatchDifferences_Asex$Genus) #44 unique genera. 

# ii. SPECIES - expect all different (unless both unassigned )
asexSporeGenusMatches[which((asexSporeGenusMatches$Species.x == asexSporeGenusMatches$Species.y) == FALSE),]
asexSporeGenusMatches[which((asexSporeGenusMatches$Species.x == asexSporeGenusMatches$Species.y) == TRUE),] #in no cases is this is same, as expected!
# Only one genus column so these are all the same
# So remove un-needed columns for cleaner merging below:
colnames(asexSporeGenusMatches) %in% c("Family.y", "Species.y", "infraspecies") #list those I want to drop
asexSporeGenusMatches_cleaned <- asexSporeGenusMatches[,-which(colnames(asexSporeGenusMatches) %in% c("Family.y", "Species.y", "infraspecies")==TRUE)]
colnames(asexSporeGenusMatches_cleaned) #looks good!
# Change names to remove ".x"
colnames(asexSporeGenusMatches_cleaned)[colnames(asexSporeGenusMatches_cleaned) %in% c("Family.x", "Species.x")]
colnames(asexSporeGenusMatches_cleaned)[colnames(asexSporeGenusMatches_cleaned) %in% c("Family.x", "Species.x")] <- c("Family", "Species")
# Keep ANCOM genusSpecies.x as genusSpecies and rename the other one as "genusSpeciesDatabase"
colnames(asexSporeGenusMatches_cleaned)[colnames(asexSporeGenusMatches_cleaned) %in% c("genusSpecies.x", "genusSpecies.y")] <- c("genusSpecies", "genusSpeciesDatabase")
sort(colnames(asexSporeGenusMatches_cleaned))
# View(asexSporeGenusMatches_cleaned)
# Are there any duplicate rows?
nrow(asexSporeGenusMatches_cleaned) #12468
nrow(distinct(asexSporeGenusMatches_cleaned)) #7147
# Remove totally duplicate rows
asexSporeGenusMatches_cleaned <- distinct(asexSporeGenusMatches_cleaned)
nrow(asexSporeGenusMatches_cleaned) #7147, as above

# 3. SET UP NO MATCH AVAILABLE  Assign still unmatched rows as "noMatchAvailable". airLeafANCOM_genusOnlyASex has all without species match,asexSporeGenusMatches_cleaned is just genus.
# So, this line gets those ASVs for which there are no genus or species matches
airLeafANCOM_noMatchesASex <- airLeafANCOM_genusOnlyASex[which(airLeafANCOM_genusOnlyASex$ASV_name %in% asexSporeGenusMatches_cleaned$ASV_name == FALSE),]
unique(airLeafANCOM_noMatchesASex$ASV_name %in%asexSporeGenusMatches_cleaned$ASV_name) #FALSE great, these were removed
unique(airLeafANCOM_noMatchesASex$ASV_name %in%asexSporeExactMatches_cleaned$ASV_name) #FALSE great, these were removed
airLeafANCOM_noMatchesASex <- airLeafANCOM_noMatchesASex %>% 
  mutate(
    size = "noDatabaseMatch",
    match = "noMatchAvailable"
  )
sort(colnames(airLeafANCOM_noMatchesASex))
# NO duplicates here, as expected
nrow(airLeafANCOM_noMatchesASex) == nrow(distinct(airLeafANCOM_noMatchesASex))

# 4. CHECK NO OVERLAP IN SPECIES MATCHES, GENUS MATCHES, AND NON MATCHED
# Check no overlap in ASVs to insure that ASVs aren't being double counted
intersect(airLeafANCOM_noMatchesASex$ASV_name, asexSporeGenusMatches_cleaned$ASV_name) #no overlap
intersect(asexSporeGenusMatches_cleaned$ASV_name, asexSporeExactMatches_cleaned$ASV_name) #no overlap 
colnames(airLeafANCOM_noMatchesASex); colnames(asexSporeGenusMatches_cleaned); colnames(asexSporeExactMatches_cleaned)
dim(airLeafANCOM_noMatchesASex); dim(asexSporeGenusMatches); dim(asexSporeExactMatches)
# View(airLeafANCOM_noMatchesASex)

# 5. CHECK SPECIFIC SPORE TYPES -- this may explain part of issue: Most in air are chlamydospores,
# which are under fundamentally different selection pressures (resistant structures, long term in soil)
table(asexSporeGenusMatches$Specific_sporeName, asexSporeGenusMatches$ANCOMcat)

####### 4.SUMMARIZING ASEXUAL SPORES BY FULL (asexSporeExactMatches_cleaned) AND GENUS (asexSporeGenusMatches_cleaned) MATCHES #######
# (this is where diverges most substantially from ANCOM_traits_I6SandITS.R)
# 1. INFER Mean, median, range, etc. per ASV for asexSporeExactMatches_cleaned
sort(colnames(asexSporeExactMatches_cleaned))
asexualExact_sum <- asexSporeExactMatches_cleaned %>%
  group_by(ASV_name, Phylum, Class, Order, Family, Genus, Species, lfc_sampleTypephyllosphere, ANCOMcat) %>%
  # group_by all of of these to retain rows, and is fine since each ASV matches in this respect
  # should retain ASVs that have the same species designation, e.g. ASV_447 and ASV_410, both B. adjusta with 4 separate sizes)
  summarize(
    numberOfMatches = n_distinct(SporeVolume), #this can be that it is in asexSporeExactMatches_cleaned, asexSporeGenusMatches_cleaned, airLeafANCOM_noMatchesASex. So all should have this
    meanVol = mean(SporeVolume, na.rm = TRUE),
    medianVol = median(SporeVolume, na.rm = TRUE), 
    minVol = min(SporeVolume, na.rm = TRUE),
    maxVol = max(SporeVolume, na.rm = TRUE),
    rangeVol = max(SporeVolume, na.rm = TRUE) - min(SporeVolume, na.rm = TRUE),
    .groups = "drop"
  )
# View(asexualExact_sum)
length(unique(asexSporeExactMatches_cleaned$ASV_name)) == nrow(asexualExact_sum); nrow(asexualExact_sum)  #43, matches before as expected!
sort(unique(asexualExact_sum$ASV_name)) == sort(unique(asexSporeExactMatches_cleaned$ASV_name)) #these are the same!
which(is.nan(asexualExact_sum$meanVol)) #all have amounts


# 2. INFER Mean, median, range, etc. per ASV for GENUS LEVEL MATCHES (asexSporeGenusMatches_cleaned)
sort(colnames(asexSporeGenusMatches_cleaned))
asexualGenLevel_sum <- asexSporeGenusMatches_cleaned %>%
  group_by(ASV_name, Phylum, Class, Order, Family, Genus, lfc_sampleTypephyllosphere, ANCOMcat) %>%
  # Do not group by Species
  # should retain ASVs that have the same species designation, e.g. ASV_447 and ASV_410, both B. adjusta with 4 separate sizes)
  summarize(
    numberOfMatches = n_distinct(SporeVolume), #this can be that it is in asexSporeGenusMatches_cleaned, asexSporeGenusMatches_cleaned, airLeafANCOM_noMatchesASex. So all should have this
    meanVol = mean(SporeVolume, na.rm = TRUE),
    medianVol = median(SporeVolume, na.rm = TRUE), 
    minVol = min(SporeVolume, na.rm = TRUE),
    maxVol = max(SporeVolume, na.rm = TRUE),
    rangeVol = max(SporeVolume, na.rm = TRUE) - min(SporeVolume, na.rm = TRUE),
    .groups = "drop"
  )
# View(asexualGenLevel_sum)
length(unique(asexSporeGenusMatches_cleaned$ASV_name)) == nrow(asexualGenLevel_sum) #254, matches before as expected!
unique(sort(unique(asexualGenLevel_sum$ASV_name)) == sort(unique(asexSporeGenusMatches_cleaned$ASV_name))) #these are the same!
which(is.nan(asexualGenLevel_sum$meanVol)) #all have amounts

# 3. INVESTIGATE Double check that these are working as expected -- yay! 
# i. (these show that each genera has same mean volume, even though species (UNITE database) differs and is still shown in df)
# Genus == "Didymella"
median(asexSporeGenusMatches_cleaned$SporeVolume[which(asexSporeGenusMatches_cleaned$Genus == "Didymella")]) == asexualGenLevel_sum$medianVol[which(asexualGenLevel_sum$Genus == "Didymella")]
asexualGenLevel_sum$numberOfMatches[which(asexualGenLevel_sum$Genus == "Didymella")] #49 matches wow
# Genus == "Pseudocercospora"
median(asexSporeGenusMatches_cleaned$SporeVolume[which(asexSporeGenusMatches_cleaned$Genus == "Pseudocercospora")]) == asexualGenLevel_sum$medianVol[which(asexualGenLevel_sum$Genus == "Pseudocercospora")]
asexualGenLevel_sum$numberOfMatches[which(asexualGenLevel_sum$Genus == "Pseudocercospora")] #578 matches wow
# ii. Does each ASV have same mean spore volume across the multiple matches per ASV?
asexualGenLevel_sumCHECK <- asexualGenLevel_sum %>% 
  group_by(ASV_name) %>% 
  summarize(nVolsPerASV = n_distinct(medianVol))
unique(asexualGenLevel_sumCHECK$nVolsPerASV) #all 1, confirms that matching as expected.
# iii. Number of matches expected?
nrow(asexualGenLevel_sum) == nrow(distinct(asexualGenLevel_sum)); nrow(asexualGenLevel_sum) #yes! 254

# 4. KEEP ONLY MATCHES (exact then genus matches) WHERE MATCH VARIATION RELATIVE TO THE MEDIAN IS SMALL
# In edge effects paper, "retained the ASV only if mean calculated genome size was less than 
# 15% of the range between the smallest and the largest matched genomes.
# ASV only if the matched genome sizes were similar: ((max−min)/median < 0.15)
# i. Exact (species-level) matches
head(asexualExact_sum)
colnames(asexualExact_sum)
asexExact_trim <-  asexualExact_sum %>% 
  mutate(rangeDivMed = rangeVol/medianVol) %>% #keep only spores where matches are close together, relative to the median
  mutate(keep = rangeVol < 0.15 * medianVol) 
table(asexExact_trim$keep, asexExact_trim$ANCOMcat) #keep only 38
# NEW
# bioaerosol foliar surface
# FALSE         10             19
# TRUE          19             19

# OLD
# air phyllo
# FALSE   4      9
# TRUE   13     17
asexExact_trim2 <- asexExact_trim %>% 
  filter(keep == TRUE)
table(asexExact_trim2$ANCOMcat) #matches above

# ii. Genus-level matches
head(asexualGenLevel_sum)
colnames(asexualGenLevel_sum)
asexGen_trim <-  asexualGenLevel_sum %>% 
  mutate(rangeDivMed = rangeVol/medianVol) %>% #keep only spores where matches are close together, relative to the median
  mutate(keep = rangeVol < 0.15 * medianVol) 
# View(asexGen_trim)
table(asexGen_trim$keep, asexGen_trim$ANCOMcat) #keep only 32 
# NEW:
# bioaerosol foliar surface
# FALSE        103             87
# TRUE          54             10
# OLD:    air phyllo
# FALSE  33     70
# TRUE   23      9
asexGen_trim2 <- asexGen_trim %>% #keep only those that aren't super variable
  filter(keep == TRUE)
table(asexGen_trim2$ANCOMcat) #54             10  same as above

####### 5. JOIN GENUS AND SPECIES MATCHES INTO ONE DATAFRAME TO USE LATER ON! ######
# 1. CHECK Next 2 lines show only diff is asexualExact_sum has Species column
setdiff(colnames(asexualExact_sum), colnames(asexualGenLevel_sum))
setdiff(colnames(asexualGenLevel_sum), colnames(asexualExact_sum))
# 2. ADD DETAILS BEFORE MERGING
asexualGenLevel_sum$Species <- NA
asexualGenLevel_sum$matchType <- "genusMatch"
asexualExact_sum$matchType <- "speciesMatch"
# 3. ROW BIND!
asexSporeMatches <- rbind(asexualExact_sum,asexualGenLevel_sum)

# 4. EXPLORE AND THEN GET MEDIANS FOR MANUSCRIPT
median(asexSporeMatches$medianVol[asexSporeMatches$ANCOMcat == "foliar surface"]) #89.79719 (_____ reported in MS)
median(asexSporeMatches$medianVol[asexSporeMatches$ANCOMcat == "bioaerosol"]) #331.8307 (_____ reported in MS)

####### 6. RESTRICTED BY VARIATION: JOIN GENUS AND SPECIES MATCHES INTO ONE DATAFRAME TO USE LATER ON! ######
# 1. CHECK Next 2 lines show only diff is asexExact_trim2 has Species column
setdiff(colnames(asexExact_trim2), colnames(asexGen_trim2))
setdiff(colnames(asexGen_trim2), colnames(asexExact_trim2))
# 2. ADD DETAILS BEFORE MERGING
asexGen_trim2$Species <- NA
asexGen_trim2$matchType <- "genusMatch"
asexExact_trim2$matchType <- "speciesMatch"
# 3. ROW BIND!
asexSporeMatches_trim <- rbind(asexExact_trim2,asexGen_trim2)
# saveRDS(asexSporeMatches_trim, file = "RobjectsSaved/asexSporeMatches_trim_Sept2.RData")
# 4. EXPLORE AND THEN GET MEDIANS-- LIKELY NOT GOING IN MANUSCRIPT
median(asexSporeMatches_trim$medianVol[asexSporeMatches_trim$ANCOMcat == "phyllo"]) #247.7277 (247.7)
median(asexSporeMatches_trim$medianVol[asexSporeMatches_trim$ANCOMcat == "air"]) #1059.6 (1059.6)
table(asexSporeMatches_trim$ANCOMcat)

####### 7. PLOTS! ######
# 1. Genus and species matches
colnames(asexSporeMatches)
asexSpore_ViolinPlot <- ggplot(data=asexSporeMatches, aes(x=ANCOMcat, y=log(medianVol))) +
  geom_violin(linewidth = 1.25) + # (aes(color= ANCOMcat), linewidth =1.5) thsi would change outline
  geom_jitter(aes(color = ANCOMcat), size=4, alpha=0.6, height = 0, width = 0.35) +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    legend.position="none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16)) +
  scale_color_manual(values = c("air" = "#004CFF", "phyllo" = "darkgreen")) +
  scale_x_discrete(labels = c("air" = "air", "phyllo" = "leaf surfaces")) +
  labs(y = "Log of Median Spore Volume") +
  ggtitle("Fungal Uninucleate Asexual Spores") +
  ylim(-1, 18.5)  #Change y-axis so that there is enough room on top to add mean spore sizes
asexSpore_ViolinPlot

# 2. Just species
asexSpore_ViolinPlot_exactMatch <- ggplot(data=asexualExact_sum, aes(x=ANCOMcat, y=log(medianVol))) +
  geom_violin(linewidth = 1.25) + # (aes(color= ANCOMcat), linewidth =1.5) thsi would change outline
  geom_jitter(aes(color = ANCOMcat), size=4, alpha=0.6, height = 0, width = 0.35) +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    legend.position="none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16)) +
  scale_color_manual(values = c("air" = "#004CFF", "phyllo" = "darkgreen")) +
  scale_x_discrete(labels = c("air" = "air", "phyllo" = "leaf surfaces")) +
  labs(y = "Log of Median Spore Volume") +
  ggtitle("Fungal Uninucleate Asexual Spores") +
  ylim(-1, 18.5)  #Change y-axis so that there is enough room on top to add mean spore sizes
asexSpore_ViolinPlot_exactMatch

# 3. MEAN VOL (as in, as it was originally performed)
log(max(asexSporeMatches$meanVol)); log(min(asexSporeMatches$meanVol)) #check out what means and mins to plot
asexSpore_VioPlot_MEAN <- ggplot(data=asexSporeMatches, aes(x=ANCOMcat, y=log(meanVol))) +
  geom_violin(linewidth = 1.25) + # (aes(color= ANCOMcat), linewidth =1.5) thsi would change outline
  geom_jitter(aes(color = ANCOMcat), size=4, alpha=0.6, height = 0, width = 0.35) +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    legend.position="none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16)) +
  scale_color_manual(values = c("air" = "#004CFF", "phyllo" = "darkgreen")) +
  scale_x_discrete(labels = c("air" = "air", "phyllo" = "leaf surfaces")) +
  labs(y = "Log of Mean Spore Volume") +
  ggtitle("Fungal Uninucleate Asexual Spores-MEAN") +
  ylim(-1, 18.5)  #Change y-axis so that there is enough room on top to add mean spore sizes
asexSpore_VioPlot_MEAN

# 4. WITH MORE RESTRICTED MATCHES
log(max(asexSporeMatches_trim$meanVol)); log(min(asexSporeMatches_trim$meanVol)) #check out what means and mins to plot
asexSpore_VioPlot_restricted <- ggplot(data=asexSporeMatches_trim, aes(x=ANCOMcat, y=log(meanVol))) +
  geom_violin(linewidth = 1.25) + # (aes(color= ANCOMcat), linewidth =1.5) thsi would change outline
  geom_jitter(aes(color = ANCOMcat), size=4, alpha=0.6, height = 0, width = 0.35) +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    legend.position="none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16)) +
  scale_color_manual(values = c("air" = "#004CFF", "phyllo" = "darkgreen")) +
  scale_x_discrete(labels = c("air" = "air", "phyllo" = "leaf surfaces")) +
  labs(y = "Log of Mean Spore Volume") +
  ggtitle("Fungal Uninucleate Asexual Spores-restricted") +
  ylim(-1, 18.5)  #Change y-axis so that there is enough room on top to add mean spore sizes
asexSpore_VioPlot_restricted

####### 8. STATS! ######
# 1. Genus and species match
colnames(asexSporeMatches)
asexSporeMatches$medianVol
which(is.na(asexSporeMatches$medianVol)) #no NAs in data
# Check to see if volume is normally distributed with Shapiro-Wilk test
qqnorm(asexSporeMatches$medianVol)
shapiro.test(x = asexSporeMatches$medianVol) # W = 0.14204, p-value < 2.2e-16 these data are NOT normal (reject null that is normal)
# Does log transforming these make them normal (this also matches the likely plots...). NO!
# Q-Q plot - not quite, there are still those weird, GIANT phyllosphere ones...
qqnorm(log(asexSporeMatches$medianVol))
qqline(log(asexSporeMatches$medianVol))
shapiro.test(x = log(asexSporeMatches$medianVol)) # W = 0.8361, p-value = 7.192e-13 these data are NOT normal (reject null that is normal)
# Given that these are non-normal, perform non-parametric test (Wilcoxon rank sum test/Mann-Whitney test)
asexualSporeWilcoxon <- wilcox.test(medianVol ~ ANCOMcat, data = asexSporeMatches)
asexualSporeWilcoxon # W = 4945, p-value = 0.001003 

# 2. Just exact, species match
colnames(asexualExact_sum)
asexualExact_sum$medianVol
which(is.na(asexualExact_sum$medianVol)) #no NAs in data
# Check to see if volume is normally distributed with Shapiro-Wilk test
shapiro.test(x = asexualExact_sum$medianVol) # W = 0.47578, p-value = 3.405e-11 these data are NOT normal (reject null that is normal)
# Does log transforming these make them normal (this also matches the likely plots...). NO!
# Q-Q plot - not quite
qqnorm(log(asexualExact_sum$medianVol))
qqline(log(asexualExact_sum$medianVol))
shapiro.test(x = log(asexualExact_sum$medianVol)) #W = 0.81514, p-value = 7.741e-06 these data are NOT normal (reject null that is normal)
# Given that these are non-normal, perform non-parametric test (Wilcoxon rank sum test/Mann-Whitney test)
asexualSporeWilcoxon_exact <- wilcox.test(medianVol ~ ANCOMcat, data = asexualExact_sum)
asexualSporeWilcoxon_exact # W = 220, p-value = 0.9901 Shows that there is no difference among groups

# 3. More restricted genus species match
colnames(asexSporeMatches_trim)
asexSporeMatches_trim$medianVol
which(is.na(asexSporeMatches_trim$medianVol)) #no NAs in data
# Check to see if volume is normally distributed with Shapiro-Wilk test
shapiro.test(x = asexSporeMatches_trim$medianVol) # W = 0.39888, p-value = 1.646e-14 these data are NOT normal (reject null that is normal)
# Does log transforming these make them normal (this also matches the likely plots...). NO!
# Q-Q plot - not quite
qqnorm(log(asexSporeMatches_trim$medianVol))
qqline(log(asexSporeMatches_trim$medianVol))
shapiro.test(x = log(asexSporeMatches_trim$medianVol)) #W = 0.8824, p-value = 2.423e-05 these data are NOT normal (reject null that is normal)
# Given that these are non-normal, perform non-parametric test (Wilcoxon rank sum test/Mann-Whitney test)
asexualSporeWilcoxon_restricted <- wilcox.test(medianVol ~ ANCOMcat, data = asexSporeMatches_trim)
asexualSporeWilcoxon_restricted # W = 513, p-value = 0.523 Shows that there is no difference among groups

