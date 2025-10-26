# I6S_souceTracking_Sept.R
# This version started September 22, 2025
# Earlier version of this script was I6S_souceTracking_final.R.

# Script description: Aims to explore where bacterial reads in the air are coming from (i.e., from soil or phyllosphere sources)

# 1. Indicator analyses (indcator species analysis to test for different taxa between soil and phyllosphere):
## a. Plots: plots showing proportions of soil and phyllo indicator species in air, soil, and phyllosphere samples
### i. Two paneled plot, with the top plot showing soil indicator taxa across air (matrix and patch combined),
# phyllosphere, and soil samples, and bottom phyllosphere indicator taxa across air, phyllosphere, and soil samples
### ii. Proportions of only air samples, but separated out by forest/savanna
### iii. NEED TO ADD: PROPORTION IN EACH HABITAT BASED ON 
## b. Analyses (based on initial indicator species analysis results):
#### i. Separate Wilcoxon signed-rank tests to look for differences in the proportion of foliar surface or soil 
# indicator taxa in forested matrix or open patch samples


#######################################
#             SCRIPT SET UP
#######################################
# LOAD REQUIRED PACKAGES:
library("phyloseq") #‘1.32.0’
library("tidyverse") #‘1.3.1’
library("vegan") #vegan 2.5-7
library("gridExtra")    #allows you to make multiple plots on the same page with ggplot, 
library("indicspecies"); packageVersion("indicspecies") #‘1.8.0’
library("Polychrome") #for color palette
library("ggVennDiagram")
library("magrittr")
library("dplyr")
library("ggbeeswarm")
library("ggplot2"); packageVersion("ggplot2") #‘4.0.0’

# HANDY FUNCTION TO GET ASV TABLE OUT OF PHYLOSEQ:
# Function to get ASV table out of phyloseq so that we can #View it better
# (Inspired by similar function at https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  return(as.data.frame(ASVTable))
}

# LOAD DATA (made/saved in full16S_EDArarefied_Part2_Sept.R)
all16Sr_noAirSingsDoubs.ps <- readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_noAirSingsDoubs_ps_sept25.rds")
all16Sr_noAirSingsDoubs.ps
sample_data(all16Sr_noAirSingsDoubs.ps)
# Confirms that this has no bioaerosol singletons or doubletons
air16SSampsColsIndex <- grepl(x=colnames(otu_table(all16Sr_noAirSingsDoubs.ps)), pattern= "air_")
unique(sort(rowSums(otu_table(all16Sr_noAirSingsDoubs.ps)[,air16SSampsColsIndex]))) #(no values of 1 or 2)

length(which(sample_data(all16Sr_noAirSingsDoubs.ps)$sampleType == "soil")) #157 samples
length(which(sample_data(all16Sr_noAirSingsDoubs.ps)$sampleType == "air")) #84 samples
length(which(sample_data(all16Sr_noAirSingsDoubs.ps)$sampleType == "phyllosphere")) #58 samples
plantDat_df <- sample_data(all16Sr_noAirSingsDoubs.ps)[which(sample_data(all16Sr_noAirSingsDoubs.ps)$sampleType == "phyllosphere"),]

table(plantDat_df$EU)
table(plantDat_df$EU, plantDat_df$PlantSpecies)


############################################
# I. SOIL AND PHYLLOSPHERE ONLY INDICATOR SPECIES ANALYSIS 
############################################
# For this, I perform an indicator species analysis on only soil and phyllosphere samples

##### 1. PERFORM INDICATOR SPECIES ANALYSIS #####
# This will perform an indicator species analysis to determine which are soil and phyllosphere specialists

# Make phyloseq object with only phyllosphere and soil samples
phylloSoilonly_I6S.ps <- subset_samples(all16Sr_noAirSingsDoubs.ps, sampleType!= "air")

# Get the metadata for these samples
phylloSoilonly_I6S_meta <- as.data.frame(as.matrix(sample_data(phylloSoilonly_I6S.ps)))
unique(phylloSoilonly_I6S_meta$HabitatAir) #only phyllosphere and soil
# Get a string that has the soil versus phyllosphere designation
soilPhyllos <- phylloSoilonly_I6S_meta$HabitatAir

# Get ASV and taxonomy tables out of phyloseq:
phylloSoilonly_I6S_ASVs <- t(ASVs_outta_ps(phylloSoilonly_I6S.ps))
# phylloSoilonly_I6S_ASVs needs to be columns as ASVs and rows as samples
if (nrow(phylloSoilonly_I6S_ASVs)> ncol(phylloSoilonly_I6S_ASVs)) {
  phylloSoilonly_I6S_ASVs <- t(phylloSoilonly_I6S_ASVs)
} else {
  print("As needed for multipatt, ASVs are columns!")
}
phylloSoilonly_I6S_ASVs 

phylloSoilonly_I6S_tax <- as.data.frame(phyloseq::tax_table(phylloSoilonly_I6S.ps), stringsAsFactors = F)
head(phylloSoilonly_I6S_tax) #looks good!

# Before running below. Because this is true, then the ASV table has same order of samples as soilPhyllos, which was built from phylloSoilonly_I6S_meta
unique(rownames(phylloSoilonly_I6S_meta) == rownames(phylloSoilonly_I6S_ASVs))

# Perform indicator analysis --commented out below because if it runs it takes forever
set.seed(9)
# I6S_phylloSoilIndic <- multipatt(x=phylloSoilonly_I6S_ASVs, cluster=soilPhyllos, func="r.g", duleg = TRUE, control=how(nperm = 9999)) #duleg=true means that combinations
# are not considered
summary(I6S_phylloSoilIndic) #shows all of the taxon IDs for ASVs associated with each group, as well as
# the  r.g. value (under stat, this is the correlation index that takes into account the
# variation within and among groups), and the p-value from a permutation test.

# re-run and saved on Aug. 5, 2024 (on server)
# save(I6S_phylloSoilIndic, file = "RobjectsSaved/I6S_phylloSoilIndic_Aug5_2024")
# re-run and re-saved September 22, 2025 on own computer
# save(I6S_phylloSoilIndic, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/I6S_phylloSoilIndic_Sept_22_2025")

# HERE, I JUST LOAD IT IN CASE THE CODE ABOVE HASN'T BEEN RUN (WHICH DOES TAKE A WHILE)
load(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/I6S_phylloSoilIndic_Sept_22_2025")
summary(I6S_phylloSoilIndic)

##### 2. CHECK OUT AND RE-FORMAT RESULTS ##### 
# View(I6S_phylloSoilIndic$sign)
I6S_phylloSoilIndic_Results <- I6S_phylloSoilIndic$sign
I6S_phylloSoilIndic_.05 <- I6S_phylloSoilIndic_Results %>% 
  filter(p.value < 0.05) #filter out ASVs that had a sig rating of 0.05 or greater
dim(I6S_phylloSoilIndic_.05) #10031 ASVs left
# View(I6S_phylloSoilIndic_.05)

# ASSIGN SOIL OR PHYLLOSPHERE BASED ON INDEX
I6S_phylloSoilIndic_.05$group <- NA #add in another column, group
for (j in 1:nrow(I6S_phylloSoilIndic_.05)){
  if (I6S_phylloSoilIndic_.05$index[j] == "1") {
    I6S_phylloSoilIndic_.05$group[j] <- "phyllosphere" 
  } else {
    I6S_phylloSoilIndic_.05$group[j] <- "soil" 
  }
}

# GET INDICES AND NAMES OF PHYLLOSPHERE AND SOIL SPECIALISTS
I6S_soilSpecialistsIndex <- which(I6S_phylloSoilIndic_.05$group == "soil")
length(I6S_soilSpecialistsIndex) #7411 
I6S_soilSpecialists <- rownames(I6S_phylloSoilIndic_.05)[I6S_soilSpecialistsIndex]
I6S_phyllospecialistsIndex <- which(I6S_phylloSoilIndic_.05$group == "phyllosphere") 
length(I6S_phyllospecialistsIndex) #2620
I6S_phylloSpecialists <- rownames(I6S_phylloSoilIndic_.05)[I6S_phyllospecialistsIndex]

# Make two plots that have the % of indicator taxa in each group (soil, air, and phyllo), where
# first plot is % soil and second plot is % phyllo for each of these groups.
phylloSoilonly_I6S_tax <- as.data.frame(phyloseq::tax_table(phylloSoilonly_I6S.ps), stringsAsFactors = F)
head(phylloSoilonly_I6S_tax)

ASVs_all16Sr_noAirSingsDoubs.ps <- ASVs_outta_ps(all16Sr_noAirSingsDoubs.ps)
# Get how to pull the specialists out of the big dataframe
# Soil
I6S_soilSpecialistsBigIndex <- which(rownames(ASVs_all16Sr_noAirSingsDoubs.ps) %in% I6S_soilSpecialists == TRUE)
# Next line shows that by using "BigIndex", I get these same taxa. This proof shows that the for loop later will work.
setdiff((rownames(ASVs_all16Sr_noAirSingsDoubs.ps[I6S_soilSpecialistsBigIndex,])), I6S_soilSpecialists) #these are the same, as expected! 
# Phyllo
I6S_phylloSpecialistsBigIndex <- which(rownames(ASVs_all16Sr_noAirSingsDoubs.ps) %in% I6S_phylloSpecialists == TRUE)
# Next line shows that by using "BigIndex", I get these same taxa. This proof shows that the for loop later will work.
setdiff((rownames(ASVs_all16Sr_noAirSingsDoubs.ps[I6S_phylloSpecialistsBigIndex,])), I6S_phylloSpecialists) #these are the same, as expected! 

# Pull out metadata from phyloseq object
I6Srarefied_metadat <- as.data.frame(as.matrix(sample_data(all16Sr_noAirSingsDoubs.ps)))
colnames(I6Srarefied_metadat)
unique(rownames(I6Srarefied_metadat) == colnames(ASVs_all16Sr_noAirSingsDoubs.ps)) #can keep other stuff, since order is same as in ITSrarefied_meta

# PRE-ALLOCATE AND CHECK to store the proportions of specialists in each sample (as many rows as samples)
I6S_specialistsProportion <- as.data.frame(matrix(data=NA, nrow=ncol(ASVs_all16Sr_noAirSingsDoubs.ps), ncol = 5))
colnames(I6S_specialistsProportion) <- c("sampleName", "sampleType", "totalNumbReads", "phylloSpecProp", "soilSpecProp")
I6S_specialistsProportion[,1] <- rownames(I6Srarefied_metadat) #column one gets all sample names
I6S_specialistsProportion[,2] <- I6Srarefied_metadat$sampleType #column 2 gets sampleType
unique(colnames(ASVs_all16Sr_noAirSingsDoubs.ps) == I6S_specialistsProportion[,1]) # shows that can get the number of reads in each sample as below
I6S_specialistsProportion$totalNumbReads <- colSums(ASVs_all16Sr_noAirSingsDoubs.ps)
head(I6S_specialistsProportion)
sort(unique(colSums(ASVs_all16Sr_noAirSingsDoubs.ps))) #this occasionally varies a tiny bit from the rarefied number, 5500, because of removing singletons and doubletons

# FILL IN PROPORTIONS OF EACH TYPE OF SPECIALIST IN EACH SAMPLE
for (k in 1:ncol(ASVs_all16Sr_noAirSingsDoubs.ps)){ #k indexes each sample
  # Add up all of the reads for soil specialists in the k-th column (i.e. kth sample) and divide it by total number of reads in that kth sample.And place it in the kth row for soil specialists 
  I6S_specialistsProportion[k,5] <- sum(ASVs_all16Sr_noAirSingsDoubs.ps[I6S_soilSpecialistsBigIndex,k])/I6S_specialistsProportion$totalNumbReads[k] 
  # Do the same for the phyllosphere specialists
  I6S_specialistsProportion[k,4] <- sum(ASVs_all16Sr_noAirSingsDoubs.ps[I6S_phylloSpecialistsBigIndex,k])/I6S_specialistsProportion$totalNumbReads[k]
}

# View(I6S_specialistsProportion)
# re-saved on server August 6, 2024:
# write_rds(I6S_specialistsProportion, file="RobjectsSaved/I6S_specialistsProportion_Aug6_2024.Rdata")
# saved Oct. 15, 2023
# save(I6S_specialistsProportion, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_specialistsProportion_Oct15")

# MOST RECENT: Saved September 22, 2025 on own computer
# save(I6S_specialistsProportion, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/I6S_specialistsProportion_Sept_22_2025")

##############
# II. PLOTS SHOWING PROPORTIONS OF SOIL AND PHYLLO INDICATOR SPECIES IN AIR, SOIL, AND PHYLLOSPHERE SAMPLES
##############
soilSpecProp_plot_I6S <- ggplot(I6S_specialistsProportion, aes(x=sampleType, y=soilSpecProp, fill=sampleType)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=c("cornflowerblue", "chartreuse4", "chocolate4")) +
  geom_jitter(color="black", size=1, alpha=0.9, height = 0) + #height = 0 means no horizontal jitter which would be very misleading!
  scale_x_discrete(name= "Sample Type", labels=c("air", "phyllosphere", "soil")) +
  scale_y_continuous(name= "Proportion of soil indicator taxa") + # limI6S=c(25, 225) add this argument to change limI6S
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  ggtitle("Soil indicator taxa across samples")  +
  theme(plot.title = element_text(size=20)) + #increase plot title size
  theme(legend.position = "none")  #remove legend 

# ggplot it!
phylloSpecProp_plot_I6S <- ggplot(I6S_specialistsProportion, aes(x=sampleType, y=phylloSpecProp, fill=sampleType)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=c("cornflowerblue", "chartreuse4", "chocolate4")) +
  geom_jitter(color="black", size=1, alpha=0.9, height = 0) +
  scale_x_discrete(name= "Sample Type", labels=c("air", "phyllosphere", "soil")) +
  scale_y_continuous(name= "Proportion of phyllosphere indicator taxa") + # limI6S=c(25, 225) add this argument to change limI6S
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  ggtitle("Phyllosphere indicator taxa across samples")  +
  theme(plot.title = element_text(size=20)) + #increase plot title size
  theme(legend.position = "none")  #remove legend 

# quartz()
grid.arrange(soilSpecProp_plot_I6S, phylloSpecProp_plot_I6S, ncol=1)

# For manuscipt: what is the mean proportion of foliar surface-associated taxa in air samples?
# View(I6S_specialistsProportion)
colnames(I6S_specialistsProportion)
I6S_specialistsProportion %>% 
  filter(sampleType == "air") %>% 
  summarize(meanPhylloAllSample =mean(phylloSpecProp))
# Result of the above is 0.3327747, so 33.3%.

##############
# PLOTS SHOWING PROPORTIONS OF ONLY AIR SAMPLES, BUT SEPARATED OUT BY FOREST/SAVANNA
# Some grayed out and not re-checked (although likely correct) since not in submitted version of manuscript)
##############
#View(I6S_specialistsProportion)
# View(I6Srarefied_metadat) #look at the sample meta data. HabitatAir is what we need to divide the samples up!

# Get savanna air names
I6S_airSavSampNames <- rownames(I6Srarefied_metadat %>%
                                  filter(sampleType == "air" & HabitatAir == "savanna"))

# Get forest air names
I6S_airForestSampNames <- rownames(I6Srarefied_metadat %>%
                                     filter(sampleType == "air" & HabitatAir == "forest"))

intersect(I6S_airSavSampNames, I6S_airForestSampNames) #these are different, as they should be!

# Add in a column that specifies if that air sample is forest or savanna
I6S_specialistsProportion_airOnly <- I6S_specialistsProportion %>%
  filter(sampleType == "air") %>% #pull out only air samples
  mutate(HabitatAir = ifelse(sampleName %in% I6S_airForestSampNames, "forest", "savanna")) #calls it forest if in I6S_airForestSampNames, otherwise, "savanna"

# Double check that this looks good:
savNamesCheck_I6S <- I6S_specialistsProportion_airOnly %>%
  filter(HabitatAir == "savanna") %>%
  select(sampleName)

sort(savNamesCheck_I6S$sampleName) == sort(I6S_airSavSampNames) #looks good!

# # FINALLY, PLOT THEM!
# # SOIL COMPARISON
# airOnlysoilSpecProp_plot_I6S <- ggplot(I6S_specialistsProportion_airOnly, aes(x=HabitatAir, y=soilSpecProp, fill=HabitatAir)) +
#   geom_boxplot() +
#   theme_bw() +
#   scale_fill_manual(values=c("darkgreen", "goldenrod")) +
#   geom_jitter(color="black", size=1, alpha=0.9, height = 0) +
#   scale_x_discrete(name= "Underlying vegetation type", labels=c("forest", "savanna")) +
#   scale_y_continuous(name= "Proportion of soil indicator taxa") + # limI6S=c(25, 225) add this argument to change limI6S
#   theme(axis.text=element_text(size=18),
#         axis.title=element_text(size=18)) +
#   ggtitle("Soil bacterial indicator taxa across samples")  +
#   theme(plot.title = element_text(size=20)) + #increase plot title size
#   theme(legend.position = "none")  #remove legend
#
# # quartz()
# airOnlysoilSpecProp_plot_I6S
#
# # PHYLLOSPHERE COMPARISON
# airOnlyphylloSpecProp_plot_I6S <- ggplot(I6S_specialistsProportion_airOnly, aes(x=HabitatAir, y=phylloSpecProp, fill=HabitatAir)) +
#   geom_boxplot() +
#   theme_bw() +
#   scale_fill_manual(values=c("darkgreen", "goldenrod")) +
#   geom_jitter(color="black", size=1, alpha=0.9, height = 0) +
#   scale_x_discrete(name= "Underlying vegetation type", labels=c("forest", "savanna")) +
#   scale_y_continuous(name= "Proportion of phyllosphere indicator taxa") + # limI6S=c(25, 225) add this argument to change limI6S
#   theme(axis.text=element_text(size=18),
#         axis.title=element_text(size=18)) +
#   ggtitle("Phyllosphere bacterial indicator taxa across samples")  +
#   theme(plot.title = element_text(size=20)) + #increase plot title size
#   theme(legend.position = "none")  #remove legend
#
# # quartz()
# airOnlyphylloSpecProp_plot_I6S
#
# #### MAKE A TWO PANELED PLOT ###
# # Re-shape dataframe to make proportions all in there!
# colnames(I6S_specialistsProportion_airOnly)
# I6S_specialistsProportion_air_long <- I6S_specialistsProportion_airOnly %>%
#   pivot_longer(cols = c(phylloSpecProp, soilSpecProp),
#                names_to = "proportionType",
#                values_to = "proportionOfReads")
# #View(I6S_specialistsProportion_air_long)
#
# # Vector for new facet labels
# facet_labels <- c("phylloSpecProp" = "Leaf surfaces",
#                   "soilSpecProp" = "Soil")
#
# # Create faceted plot
# I6S_airOnly_phylloSoil_2panels <- ggplot(I6S_specialistsProportion_air_long, aes(x=HabitatAir, y=proportionOfReads, fill=HabitatAir)) +
#   geom_boxplot() +
#   facet_wrap(~proportionType, scales = "fixed", labeller = labeller(proportionType = facet_labels)) +
#   theme_bw() +
#   scale_fill_manual(values=c("darkgreen", "goldenrod")) +
#   geom_jitter(color="black", size=1, alpha=0.9, height = 0) + #height = 0 insures no horizontal jitter
#   scale_x_discrete(labels=c("forest", "savanna")) +
#   theme(axis.text=element_text(size=16),
#         axis.title=element_text(size=16)) +
#   ggtitle("Bacterial indicator taxa across air samples")  +
#   labs(y = "Proportion of reads from\nleaf or soil indicator taxa",
#        x = NULL) +
#   theme(plot.title = element_text(size=16)) + #increase plot title size
#   theme(legend.position = "none") + #remove legend
#   ylim(0, 1.00) + #make y limits between 0 and 1
#   theme(strip.text = element_text(size = 16))
#
# I6S_airOnly_phylloSoil_2panels #exported as /RobjectsSaved/I6S_airOnly_phylloSoil_2panels.pdf
# # saveRDS(I6S_airOnly_phylloSoil_2panels, file = "RobjectsSaved/I6S_airOnly_phylloSoil_2panels_Jan8_2024.rds") #saved January 8, 2024


#### SAVE ALL PLOTS ####
# ON SERVER: saved August 6, 2024:
# write_rds(phylloSpecProp_plot_I6S, file = "RobjectsSaved/phylloSpecProp_plot_I6S_Aug6_2024")
# write_rds(soilSpecProp_plot_I6S, file = "RobjectsSaved/soilSpecProp_plot_I6S_Aug6_2024")
# write_rds(airOnlyphylloSpecProp_plot_I6S, file = "RobjectsSaved/airOnlyphylloSpecProp_plot_I6S_Aug6_2024")
# write_rds(airOnlysoilSpecProp_plot_I6S, file = "RobjectsSaved/airOnlysoilSpecProp_plot_I6S_Aug6_2024")

# ON OWN COMPUTER: saved Oct 15, 2023
# save(airOnlyphylloSpecProp_plot_I6S, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airOnlyphylloSpecProp_plot_I6S_Oct15_2023")
# save(airOnlysoilSpecProp_plot_I6S, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airOnlysoilSpecProp_plot_I6S_Oct15_2023")

### PLOT WITH AIR FOREST AND AIR SAVANNA AND FOLIAR SURFACES AND PHYLLOSPHERE AND SOIL ####
# Savanna air names (object made above)
I6S_airSavSampNames
# Forest air names (object made above)
I6S_airForestSampNames
# Get soil names 
I6S_soilSampNames <- rownames(I6Srarefied_metadat %>% 
                                filter(sampleType == "soil"))
# Get foliar surface names 
I6S_foliarSampNames <- rownames(I6Srarefied_metadat %>% 
                                  filter(sampleType == "phyllosphere"))

# Add in a column that specifies what kind of sample it is
colnames(I6S_specialistsProportion)
I6S_specialistsProportion_allTypes <- I6S_specialistsProportion %>% 
  mutate(
    sampleTypeWithAir = case_when(
      sampleName %in% I6S_airSavSampNames ~ "SavAirSample",
      sampleName %in% I6S_airForestSampNames ~ "ForAirSample",
      sampleName %in% I6S_foliarSampNames ~ "foliarSample",
      sampleName %in% I6S_soilSampNames ~ "soilSample",
      TRUE ~ NA_character_
    )
  )
colnames(I6S_specialistsProportion_allTypes)

# Make dataframe "longer" for faceting purposes
colnames(I6S_specialistsProportion_allTypes)
I6S_specialistsProportion_allTypes_long <- I6S_specialistsProportion_allTypes %>%
  pivot_longer(cols = c(phylloSpecProp, soilSpecProp), 
               names_to = "proportionType", 
               values_to = "proportionOfReads")

# Vector for new facet labels
facet_labels <- c("phylloSpecProp" = "Foliar surfaces", 
                  "soilSpecProp" = "Soil")

# Re-order x-axis elements
unique(I6S_specialistsProportion_allTypes_long$sampleTypeWithAir)
I6S_specialistsProportion_allTypes_long$sampleTypeWithAir <- 
  factor(I6S_specialistsProportion_allTypes_long$sampleTypeWithAir, 
         levels= c("ForAirSample", "SavAirSample", "foliarSample", "soilSample"))

# Saved May 20, 2025
# saveRDS(I6S_specialistsProportion_allTypes_long, file="RobjectsSaved/I6S_specialistsProportion_allTypes_long.RData")
# MOST RECENT: Saved September 22, 2025 on own computer
# save(I6S_specialistsProportion_allTypes_long, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/I6S_specialistsProportion_allTypes_long_Sept_22_2025")

# Create faceted plot
I6S_allTypesAir_2panels <- ggplot(I6S_specialistsProportion_allTypes_long, aes(x=sampleTypeWithAir, y=proportionOfReads, fill=sampleTypeWithAir)) + 
  geom_boxplot() + 
  geom_jitter(aes(color = sampleTypeWithAir), size=1, alpha=0.9, height = 0) + #height = 0 insures no horizontal jitter
  facet_wrap(~proportionType, scales = "fixed", labeller = labeller(proportionType = facet_labels)) +
  theme_bw() +
  scale_fill_manual(values=c("cornflowerblue", "cornflowerblue", "chartreuse4", "chocolate4")) +
  scale_color_manual(values = c("black", "black", "black", "black")) + #color for jittered points
  scale_x_discrete(labels=c("bioaerosol\nforested matrix", "bioaerosol\nopen patch", "foliar surface", "soil")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("Bacterial indicator taxa across all samples")  +
  labs(y = "Proportion of reads from\nfoliar surface or soil indicator taxa",
       x = NULL) +
  theme(plot.title = element_text(size=16)) + #increase plot title size
  theme(legend.position = "none") + #remove legend 
  ylim(0, 1.00) + #make y limits between 0 and 1
  theme(strip.text = element_text(size = 16))
I6S_allTypesAir_2panels

# SAVE PLOT (saved September 22, 2025) 
# saveRDS(I6S_allTypesAir_2panels, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_allTypesAir_2panels_I6S_09-22-25")

##############
# STATISTICS

# 1. First, test if air samples collected open patch or forested matrix differed in their proportions of foliar surface 
# indicator taxa or soil indicator taxa.
# Reported APA style as here:
# https://about.illinoisstate.edu/jhkahn/apastats/#:~:text=Chi%2DSquare%20statistics%20are%20reported,35.

# Re-shape dataframe to make proportions all in there!
colnames(I6S_specialistsProportion_airOnly)
I6S_specialistsProportion_air_long <- I6S_specialistsProportion_airOnly %>%
  pivot_longer(cols = c(phylloSpecProp, soilSpecProp),
               names_to = "proportionType",
               values_to = "proportionOfReads")
# View(I6S_specialistsProportion_air_long)

# Vector for new facet labels
facet_labels <- c("phylloSpecProp" = "foliar surfaces",
                  "soilSpecProp" = "Soil")

head(I6S_specialistsProportion_air_long) #this has information on habitat type too
# Get counts of numbers of soil and foliar surface indicators, by 
# multiplying the proportion of each sample's foliar and soil indicator count numbers by the total # of reads
I6S_specialistsProportion_air_long$ReadsByHab <- I6S_specialistsProportion_air_long$proportionOfReads* I6S_specialistsProportion_air_long$totalNumbReads
head(I6S_specialistsProportion_air_long)
# Since these number of samples are so different, will need to do a proportions test
length(which(I6S_specialistsProportion_air_long$HabitatAir== "savanna")) #76
length(which(I6S_specialistsProportion_air_long$HabitatAir== "forest")) #92

# Do a few checks to ensure that I6S_specialistsProportion_air_long is correct
# Try for two samples: "air_16S_103" and "air_16S_66"
head(I6S_specialistsProportion_air_long)
# Make smaller, subsetted ASV table
ASVtab_103_66 <- otu_table(all16Sr_noAirSingsDoubs.ps)[,colnames(otu_table(all16Sr_noAirSingsDoubs.ps)) %in% c("air_16S_103", "air_16S_66")]
# Get sums of reads in each sample
ASVSums_103_66 <- colSums(ASVtab_103_66) 
# Get number of reads from foliar surface and soil indicators
foliar103_ASVs <- ASVtab_103_66[rownames(ASVtab_103_66) %in% I6S_phylloSpecialists,1]
foliar66_ASVs <- ASVtab_103_66[rownames(ASVtab_103_66) %in% I6S_phylloSpecialists,2]
soil103_ASVs <- ASVtab_103_66[rownames(ASVtab_103_66) %in% I6S_soilSpecialists,1]
soil66_ASVs <- ASVtab_103_66[rownames(ASVtab_103_66) %in% I6S_soilSpecialists,2]
# Get proportions
foliarProp_103 <- colSums(foliar103_ASVs)/ASVSums_103_66[1]
foliarProp_66 <- colSums(foliar66_ASVs)/ASVSums_103_66[2]
soilProp_103 <- colSums(soil103_ASVs)/ASVSums_103_66[1]
soilProp_66 <- colSums(soil66_ASVs)/ASVSums_103_66[2]
# Double check these calculations above against those in larger dataframe -- ALL TRUE!
# For air_16S_103
I6S_specialistsProportion_air_long[which(I6S_specialistsProportion_air_long$sampleName %in% "air_16S_103" & I6S_specialistsProportion_air_long$proportionType == "phylloSpecProp"),6] == foliarProp_103
I6S_specialistsProportion_air_long[which(I6S_specialistsProportion_air_long$sampleName %in% "air_16S_103" & I6S_specialistsProportion_air_long$proportionType == "soilSpecProp"),6] == soilProp_103
# For air_16S_66
I6S_specialistsProportion_air_long[which(I6S_specialistsProportion_air_long$sampleName %in% "air_16S_66" & I6S_specialistsProportion_air_long$proportionType == "phylloSpecProp"),6] == foliarProp_66
I6S_specialistsProportion_air_long[which(I6S_specialistsProportion_air_long$sampleName %in% "air_16S_66" & I6S_specialistsProportion_air_long$proportionType == "soilSpecProp"),6] == soilProp_66

# https://www.sthda.com/english/wiki/two-proportions-z-test-in-r
# i. Build a contingency table
I6S_specialistsSummed <- I6S_specialistsProportion_air_long %>% 
  group_by(HabitatAir, proportionType) %>% 
  summarize(numbASVsinHabSourceType = sum(ReadsByHab, na.rm = TRUE), .groups = "drop")
# View(I6S_specialistsSummed)
I6S_specialistsSummed

# Create a contingency table
I6S_specialistsSummed_Table <- as.data.frame(matrix(nrow= 2, ncol=2))
colnames(I6S_specialistsSummed_Table) <- c("foliar", "soil")
rownames(I6S_specialistsSummed_Table) <- c("forest", "savanna")
# Fill in values from the dataframe above
I6S_specialistsSummed_Table[1,1] <- I6S_specialistsSummed[1,3]
I6S_specialistsSummed_Table[1,2] <- I6S_specialistsSummed[2,3]
I6S_specialistsSummed_Table[2,1] <- I6S_specialistsSummed[3,3]
I6S_specialistsSummed_Table[2,2] <- I6S_specialistsSummed[4,3]
I6S_specialistsSummed_Table
I6S_specialistsSummed # compare to make sure these are correct and they are!
I6S_specialistsSummed_Table$totalGroup <- rowSums(I6S_specialistsSummed_Table) #add in a column 
# that specifies total in forest and savanna
I6S_specialistsSummed_Table

# Perform "two-proportions z-test" for foliar surfaces. Prop test below has numbers of foliar in each group,
# then total numbers in each group. Forest is first, so x= goes foliar forest, then foliar savanna, then total forest, then total savanna
I6S_foliarPropTest <- prop.test(x = c(I6S_specialistsSummed_Table[1,1], I6S_specialistsSummed_Table[2,1]), 
                                n = c(I6S_specialistsSummed_Table[1,3], I6S_specialistsSummed_Table[2,3]))
I6S_foliarPropTest # this shows that air samples from matrix have a higher proportion of foliar-associated reads!
# data:  c(I6S_specialistsSummed_Table[1, 1], I6S_specialistsSummed_Table[2, 1]) out of c(I6S_specialistsSummed_Table[1, 3], I6S_specialistsSummed_Table[2, 3])
# X-squared = 1665.3, df = 1, p-value < 0.00000000000000022
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   0.06445602 0.07112170
# sample estimates:
#   prop 1    prop 2 
# 0.8881641 0.8203752 


# Test for soil
I6S_specialistsSummed_Table
# x goes soil forest, then soil savanna. N goes forest total then savanna total
I6S_soilPropTest <- prop.test(x = c(I6S_specialistsSummed_Table[1,2], I6S_specialistsSummed_Table[2,2]), 
                              n = c(I6S_specialistsSummed_Table[1,3], I6S_specialistsSummed_Table[2,3]))
I6S_soilPropTest #shows that soil indicators are higher in the patch than in the matrix.
# data:  c(I6S_specialistsSummed_Table[1, 2], I6S_specialistsSummed_Table[2, 2]) out of c(I6S_specialistsSummed_Table[1, 3], I6S_specialistsSummed_Table[2, 3])
# X-squared = 1665.3, df = 1, p-value < 0.00000000000000022
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.07112170 -0.06445602
# sample estimates:
#   prop 1    prop 2 
# 0.1118359 0.1796248 
sum(I6S_specialistsSummed_Table[,3]) #178,956 entered as the total amount (i.e. N)

# T-tests or non-parametric (Wilcoxon signed-rank) to compare, ignorning habitat type
head(I6S_specialistsProportion_air_long)
str(I6S_specialistsProportion_air_long)
I6S_specialistsProportion_air_long$proportionType <- as.factor(I6S_specialistsProportion_air_long$proportionType)
str(I6S_specialistsProportion_air_long)
# Test to see if I can do a one-way ANOVA
car::leveneTest(ReadsByHab ~ proportionType, data = I6S_specialistsProportion_air_long) # shows that variances are NOT roughly homogenous, since significant
# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value            Pr(>F)    
# group   1  56.572 0.000000000003238 ***
#       166         
# Since assumptions aren't met, try (paired)

# Pivot wider to get two numbers to line up
I6S_specPropAir_wide <- I6S_specialistsProportion_air_long %>% 
  select(-ReadsByHab) %>%  #remove because next step does not "widen" this
  pivot_wider(names_from = proportionType, values_from= proportionOfReads)
I6S_specPropAir_wide
#  Wilcoxon signed-rank test (paired, one-sided: foliar > soil)
# (wilcox.test ignores exact p when there are ties/large n which is ok!)
wilcoxPropReads <- wilcox.test(
  I6S_specPropAir_wide$phylloSpecProp, I6S_specPropAir_wide$soilSpecProp,
  paired = TRUE, #because samples are not independent
  alternative = "greater", #tests greater foliar than soil
  exact = FALSE,    #avoids warnings with ties/large n
  conf.int = TRUE,  #gives a CI for the shift (one-sided CI)
  conf.level = 0.95
)
wilcoxPropReads
# Wilcoxon signed rank test with continuity correction
# data:  I6S_specPropAir_wide$phylloSpecProp and I6S_specPropAir_wide$soilSpecProp
# V = 3484, p-value = 0.00000000000001796
# alternative hypothesis: true location shift is greater than 0
# 95 percent confidence interval:
#   0.2437737       Inf
# sample estimates:
#   (pseudo)median 
# 0.2818454 

##############
# WHAT ARE THE MEAN PROPORTIONS OF LEAF AND SOIL INDICATOR TAXA IN EACH AIR SAMPLE?
##############
# This is the ASV table for all samples
ASVs_all16Sr_noAirSingsDoubs.ps
I6S_phylloSpecialists #names of foliar surface indicator taxa
length(I6S_phylloSpecialists)
I6S_soilSpecialists #names of soil indicator taxa

# GET READS PER SAMPLE
rownames(ASVs_all16Sr_noAirSingsDoubs.ps) #rownames are the ASVs
rownames(ASVs_all16Sr_noAirSingsDoubs.ps)[length(rownames(ASVs_all16Sr_noAirSingsDoubs.ps))] #this is the last ASV
colnames(ASVs_all16Sr_noAirSingsDoubs.ps) #column names are the sample names 
# How many reads are in each sample?
readsPerSample <- colSums(ASVs_all16Sr_noAirSingsDoubs.ps)
readsPerSample #good, all are about 5500, rarefying threshold
readsPerSample2 <- as.data.frame(matrix(ncol=2, nrow=length(readsPerSample)))
colnames(readsPerSample2) <- c("sampleName", "readsPerSample")
readsPerSample2[,1] <- names(readsPerSample)
readsPerSample2[,2] <- unname(readsPerSample)
readsPerSample2

head(rownames(ASVs_all16Sr_noAirSingsDoubs.ps)) #first one is ASV_4
tail(rownames(ASVs_all16Sr_noAirSingsDoubs.ps)) #last one is ASV_ASV_39610
ASVs_all16Sr_noAirSingsDoubs_long <- ASVs_all16Sr_noAirSingsDoubs.ps %>% 
  t() %>% #transpose so that ASV names are columns
  as.data.frame() %>% #re-make a dataframe because t() function un-did this
  rownames_to_column(var= "sampleName") %>% #make the now rows into a new column called sampleName
  pivot_longer(cols = ASV_4:ASV_39610, names_to = "ASV_name", values_to = "readCount")
# View(ASVs_all16Sr_noAirSingsDoubs_long)
# Add in reads per sample
ASVs_all16Sr_noAirSingsDoubs_longer <- left_join(ASVs_all16Sr_noAirSingsDoubs_long, readsPerSample2, by= "sampleName")
# Add in indicator type for each ASV
colnames(ASVs_all16Sr_noAirSingsDoubs_longer)
ASVs_all16Sr_noAirSingsDoubs_longer <- ASVs_all16Sr_noAirSingsDoubs_longer %>%
  mutate(
    indictorType = case_when(
      ASV_name %in% I6S_phylloSpecialists ~ "foliarIndicator",
      ASV_name %in% I6S_soilSpecialists ~ "soilIndicator",
      TRUE ~ NA_character_
    )
  )
# View(ASVs_all16Sr_noAirSingsDoubs_longer) # WOAH! 8,862,360 rows!

# Now get proportions of each ASV
I6S_proportionsEachASV <- ASVs_all16Sr_noAirSingsDoubs_longer #make a copy
head(I6S_proportionsEachASV)
I6S_proportionsEachASV$propASV <- I6S_proportionsEachASV$readCount/I6S_proportionsEachASV$readsPerSample
# View(I6S_proportionsEachASV)

# Add in sample type 
soilSampNames <- rownames(sample_data(all16Sr_noAirSingsDoubs.ps)[which(sample_data(all16Sr_noAirSingsDoubs.ps)$sampleType == "soil")]) 
soilSampNames
airSampNames <- rownames(sample_data(all16Sr_noAirSingsDoubs.ps)[which(sample_data(all16Sr_noAirSingsDoubs.ps)$sampleType == "air")])
airSampNames
phylloSampNames <- rownames(sample_data(all16Sr_noAirSingsDoubs.ps)[which(sample_data(all16Sr_noAirSingsDoubs.ps)$sampleType == "phyllosphere")])
phylloSampNames

I6S_proportionsEachASV <- I6S_proportionsEachASV %>%
  mutate(
    sampleType = case_when(
      sampleName %in% airSampNames ~ "airSample",
      sampleName %in% soilSampNames ~ "soilSample",
      sampleName %in% phylloSampNames ~ "leafSample",
      TRUE ~ NA_character_
    )
  )

# Split into foliar surface and soil indicators
I6S_propsI6S_soilASVs <- I6S_proportionsEachASV[which(I6S_proportionsEachASV$indictorType == "soilIndicator"),]
unique(I6S_propsI6S_soilASVs$indictorType)
I6S_propsFoliarASVs <- I6S_proportionsEachASV[which(I6S_proportionsEachASV$indictorType == "foliarIndicator"),]
unique(I6S_propsFoliarASVs$indictorType)

# Get MEAN ASV abundance for each ASV within each sample, and then sort by this:
# FOLIAR SURFACES (proportions only for foliar surface indicator taxa)
colnames(I6S_propsFoliarASVs)
meanASVPerSample_foliar <- I6S_propsFoliarASVs %>% 
  group_by(ASV_name, sampleType) %>% 
  summarize(meanASV_relAbundPerType = mean(propASV)) %>% 
  arrange(sampleType, meanASV_relAbundPerType)
# Merge with taxonomic information
all16Sr_noAirSingsDoubsTax <- as.data.frame(phyloseq::tax_table(all16Sr_noAirSingsDoubs.ps), stringsAsFactors = F)
colnames(all16Sr_noAirSingsDoubsTax)
all16Sr_noAirSingsDoubsTax <- all16Sr_noAirSingsDoubsTax %>% 
  rownames_to_column(var= "ASV_name")
str(all16Sr_noAirSingsDoubsTax)
meanASVPerSample_foliar <- left_join(meanASVPerSample_foliar, all16Sr_noAirSingsDoubsTax, by = "ASV_name")
# View(meanASVPerSample_foliar)

# Look at only air samples
meanASVPerSample_foliarOnlyAir <- meanASVPerSample_foliar %>% 
  filter(sampleType == "airSample") %>% 
  arrange(desc(meanASV_relAbundPerType))
# View(meanASVPerSample_foliarOnlyAir)

# SOIL (proportions only for soil indicator taxa)
colnames(I6S_propsI6S_soilASVs)
meanASVPerSample_soil <- I6S_propsI6S_soilASVs %>% 
  group_by(ASV_name, sampleType) %>% 
  summarize(meanASV_relAbundPerType = mean(propASV)) %>% 
  arrange(sampleType, meanASV_relAbundPerType)
# Merge with taxonomic information
meanASVPerSample_soil <- left_join(meanASVPerSample_soil, all16Sr_noAirSingsDoubsTax, by = "ASV_name")
# View(meanASVPerSample_soil)
meanASVPerSample_soilOnlyAir <- meanASVPerSample_soil %>% 
  filter(sampleType == "airSample")
# View(meanASVPerSample_soilOnlyAir)

## FINALLY, LOOK AT ALL AT ONCE TO GET TOP AIR ASVs, this will go into supplemental materials
head(I6S_proportionsEachASV)
I6S_proportionsEachASV_justAirMean <- I6S_proportionsEachASV %>% 
  filter(sampleType == "airSample") %>% #consider only air samples
  group_by(ASV_name) %>% #grouping ASVs together
  summarize(meanASV_relAbund_air = mean(propASV)) %>% #get mean propprtion of the ASV (but only across bioaerosol samples)
  arrange(desc(meanASV_relAbund_air)) #arrange so that the one with the highest mean is on top
# View(I6S_proportionsEachASV_justAirMean)

# Get just top 100
I6S_proportionsEachASV_justAir_t100 <- I6S_proportionsEachASV_justAirMean %>% 
  slice_head(n = 100)
# View(I6S_proportionsEachASV_justAir_t100)

# Now merge with taxonomy
I6S_proportionsEachASVtax_Air_100 <- left_join(I6S_proportionsEachASV_justAir_t100, all16Sr_noAirSingsDoubsTax, by = "ASV_name")
# View(I6S_proportionsEachASVtax_Air_100) 

# Add in indicator type
head(I6S_proportionsEachASV[,c(2,5)]) #head shows columns 2 and 5, which are ASV_name and indicator type
I6S_proportionsEachASVtax_Air_100 <- left_join(I6S_proportionsEachASV_justAir_t100, I6S_proportionsEachASV[,c(2,5)], by = "ASV_name")
I6S_proportionsEachASVtax_Air_100 <- I6S_proportionsEachASVtax_Air_100 %>% 
  distinct() #remove duplicate rows that were caused by merge above (was because I6S_proportionsEachASV was in long format)
# View(I6S_proportionsEachASVtax_Air_100)

# Get number of air samples that these 100 taxa occur in
colnames(I6S_proportionsEachASVtax_Air_100)
top100Air_ASVs <- unique(I6S_proportionsEachASVtax_Air_100$ASV_name)
# 1. Subset phyloseq object to be only air samples
unique(sample_data(all16Sr_noAirSingsDoubs.ps)$sampleType)
I6SairOnly.ps <- subset_samples(all16Sr_noAirSingsDoubs.ps, sampleType == "air")
head(otu_table(I6SairOnly.ps))
dim(as.data.frame(as.matrix(otu_table(I6SairOnly.ps)))) #ASVs are rows, as they should be 
# each entry below is TRUE if the value is greater than 0, and FALSE otherwise, effectively counting the occurrences!
occurrencesEachASVair <- rowSums(as.data.frame(as.matrix(otu_table(I6SairOnly.ps))) > 0)
# Edit this to make it mergable with the dataframe above (needs to be a dataframe)
occurrencesEachASVair_2 <- as.data.frame(matrix(ncol=2, nrow= length(occurrencesEachASVair)))
colnames(occurrencesEachASVair_2) <- c("ASV_name", "numAirOccurrences")
occurrencesEachASVair_2$ASV_name <- names(occurrencesEachASVair)
occurrencesEachASVair_2$numAirOccurrences <- unname(occurrencesEachASVair)
# Merge with above dataframe
I6S_proportionsEachASVtax_Air_100_v2 <- left_join(I6S_proportionsEachASVtax_Air_100, occurrencesEachASVair_2, by = "ASV_name")
I6S_proportionsEachASVtax_Air_100_v2$percentOccupancy <- (I6S_proportionsEachASVtax_Air_100_v2$numAirOccurrences/84)*100 #84 samples
# Double check with ASV_78, one of the top bioaerosol ASVs
length(which(otu_table(I6SairOnly.ps)[rownames(otu_table(I6SairOnly.ps)) %in% "ASV_78",]> 0)) #20
I6S_proportionsEachASVtax_Air_100_v2$numAirOccurrences[which(I6S_proportionsEachASVtax_Air_100_v2$ASV_name == "ASV_78")] #20 matches above!
# Double check with ASV_11, one of the top bioaerosol ASVs
length(which(otu_table(I6SairOnly.ps)[rownames(otu_table(I6SairOnly.ps)) %in% "ASV_11",]> 0)) #62 times
I6S_proportionsEachASVtax_Air_100_v2$numAirOccurrences[which(I6S_proportionsEachASVtax_Air_100_v2$ASV_name == "ASV_11")] #62 matches above!

# DOUBLE CHECK I6S_proportionsEachASV_justAirMean MADE A BIT ABOVE 
I6S_proportionsEachASV_justAirMean
# DOUBLE CHECK-- IT WORKED!
airASVsTopAirCheck_df1 <- otu_table(I6SairOnly.ps)[rownames(otu_table(I6SairOnly.ps)) %in% top100Air_ASVs,]
# Look at ASV 11. Divide how much each sample has of ASV 11 by its total read count. 
mean(airASVsTopAirCheck_df1[rownames(airASVsTopAirCheck_df1 ) %in% "ASV_11",]/colSums(otu_table(I6SairOnly.ps))) == I6S_proportionsEachASV_justAirMean[which(I6S_proportionsEachASV_justAirMean$ASV_name == "ASV_11"),2]

# GET NUMBER OF SOIL AND FOLIAR SURFACE SAMPLES THAT THESE TOP AIR TAXA OCCUR IN
# SOIL-- i.e., mean prop ASV within a given ASV for each soil sample
I6S_soilOnly.ps <-  subset_samples(all16Sr_noAirSingsDoubs.ps, sampleType == "soil")
head(otu_table(I6S_soilOnly.ps))
dim(as.data.frame(as.matrix(otu_table(I6S_soilOnly.ps)))) #ASVs are rows, as they should be 
# each entry below is TRUE if the value is greater than 0, and FALSE otherwise, effectively counting the occurrences!
occurrencesEachASV_soil <- rowSums(as.data.frame(as.matrix(otu_table(I6S_soilOnly.ps))) > 0)
# Edit this to make it mergable with the dataframe above (needs to be a dataframe)
occurrencesEachASV_soil_2 <- as.data.frame(matrix(ncol=2, nrow= length(occurrencesEachASV_soil)))
colnames(occurrencesEachASV_soil_2) <- c("ASV_name", "numSoilOccurrences")
occurrencesEachASV_soil_2$ASV_name <- names(occurrencesEachASV_soil)
occurrencesEachASV_soil_2$numSoilOccurrences <- unname(occurrencesEachASV_soil)
occurrencesEachASV_soil_topAir <- occurrencesEachASV_soil_2 %>% 
  filter(ASV_name %in% top100Air_ASVs) #get only those that are top air ASVs
# Double check with ASV_78, one of the top bioaerosol ASVs
which(otu_table(I6S_soilOnly.ps)[rownames(otu_table(I6S_soilOnly.ps)) %in% "ASV_78",]> 0) #this ASV does not appear
occurrencesEachASV_soil_topAir[which(occurrencesEachASV_soil_topAir$ASV_name == "ASV_78"),] #0 matches above!
# Double check with ASV_11, one of the top bioaerosol ASVs
length(which(otu_table(I6S_soilOnly.ps)[rownames(otu_table(I6S_soilOnly.ps)) %in% "ASV_11",]> 0)) #55 times
occurrencesEachASV_soil_topAir[which(occurrencesEachASV_soil_topAir$ASV_name == "ASV_11"),] #55 matches above!
# Add percent occupancy by dividing number of samples in by 157 samples x 100
occurrencesEachASV_soil_topAir$SoilPercentOccupancy <- occurrencesEachASV_soil_topAir$numSoilOccurrences/157*100

# GET NUMBER OF FOLIAR SURFACES AND FOLIAR SURFACE SAMPLES THAT THESE TOP AIR TAXA OCCUR IN
# FOLIAR SURFACES-- i.e., mean prop ASV within a given ASV for each foliar sample
I6S_foliarOnly.ps <-  subset_samples(all16Sr_noAirSingsDoubs.ps, sampleType == "phyllosphere")
head(otu_table(I6S_foliarOnly.ps ))
dim(as.data.frame(as.matrix(otu_table(I6S_foliarOnly.ps)))) #ASVs are rows, as they should be 
# each entry below is TRUE if the value is greater than 0, and FALSE otherwise, effectively counting the occurrences!
occurrencesEachASV_foliar <- rowSums(as.data.frame(as.matrix(otu_table(I6S_foliarOnly.ps))) > 0)
# Edit this to make it mergable with the dataframe above (needs to be a dataframe)
occurrencesEachASV_foliar_2 <- as.data.frame(matrix(ncol=2, nrow= length(occurrencesEachASV_foliar)))
colnames(occurrencesEachASV_foliar_2) <- c("ASV_name", "numFoliarOccurrences")
occurrencesEachASV_foliar_2$ASV_name <- names(occurrencesEachASV_foliar)
occurrencesEachASV_foliar_2$numFoliarOccurrences <- unname(occurrencesEachASV_foliar)
occurrencesEachASV_foliar_topAir <- occurrencesEachASV_foliar_2 %>% 
  filter(ASV_name %in% top100Air_ASVs) #get only those that are top air ASVs
# Double check with ASV_78, one of the top bioaerosol ASVs
length(which(otu_table(I6S_foliarOnly.ps)[rownames(otu_table(I6S_foliarOnly.ps)) %in% "ASV_78",]> 0)) #3
occurrencesEachASV_foliar_topAir[which(occurrencesEachASV_foliar_topAir$ASV_name == "ASV_78"),] #3 matches above!
# Double check with ASV_11, one of the top bioaerosol ASVs
length(which(otu_table(I6S_foliarOnly.ps)[rownames(otu_table(I6S_foliarOnly.ps)) %in% "ASV_11",]> 0)) #0 times
occurrencesEachASV_foliar_topAir[which(occurrencesEachASV_foliar_topAir$ASV_name == "ASV_11"),] #0 matches above!
# Add percent occupancy by dividing number of samples in by 58 foliar surface samples x 100
occurrencesEachASV_foliar_topAir$FoliarPercentOccupancy <- occurrencesEachASV_foliar_topAir$numFoliarOccurrences/58*100

# PROPORTIONS OF TOP BIOAEROSOL ASVS IN SOIL AND FOLIAR SAMPLES
# SOIL-- i.e., mean prop ASV within a given ASV for each soil sample
I6S_proportionsEachASV_justSoil <- I6S_proportionsEachASV %>%
  filter(sampleType == "soilSample") %>% #only soil
  filter(ASV_name %in% top100Air_ASVs) %>%  #only top air ASVs
  group_by(ASV_name) %>% #group by ASV_name
  summarize(meanASV_relAbund_soil = mean(propASV)) %>% #mean prop ASV within a given ASV for each soil sample
  arrange(desc(meanASV_relAbund_soil))
# View(I6S_proportionsEachASV_justSoil)
# DOUBLE CHECK-- IT WORKED!
soilASVsTopAirCheck_df1 <- otu_table(I6S_soilOnly.ps)[rownames(otu_table(I6S_soilOnly.ps)) %in% top100Air_ASVs,]
# Look at ASV 11. Divide how much each sample has of ASV 11 by its total read count. 
mean(soilASVsTopAirCheck_df1[rownames(soilASVsTopAirCheck_df1) %in% "ASV_11",]/colSums(otu_table(I6S_soilOnly.ps))) == I6S_proportionsEachASV_justSoil[which(I6S_proportionsEachASV_justSoil$ASV_name == "ASV_11"),2]

# FOLIAR
I6S_proportionsEachASV_justLeaves <- I6S_proportionsEachASV %>%
  filter(sampleType == "leafSample") %>%
  filter(ASV_name %in% top100Air_ASVs) %>%
  group_by(ASV_name) %>%
  summarize(meanASV_relAbund_foliar = mean(propASV)) %>%
  arrange(desc(meanASV_relAbund_foliar))
# View(I6S_proportionsEachASV_justLeaves)

# MERGE ALL OF THIS TOGETHER AND ADD IN TAXONOMY
I6S_topAirTable1 <- merge(I6S_proportionsEachASVtax_Air_100_v2, I6S_proportionsEachASV_justSoil, by= "ASV_name")
I6S_topAirTable2 <- merge(I6S_topAirTable1, occurrencesEachASV_soil_topAir,by= "ASV_name")
I6S_topAirTable3 <- merge(I6S_topAirTable2, I6S_proportionsEachASV_justLeaves, by= "ASV_name")
I6S_topAirTable4 <- merge(I6S_topAirTable3, occurrencesEachASV_foliar_topAir, by= "ASV_name")
I6S_topAirTable_final <- left_join(I6S_topAirTable4, all16Sr_noAirSingsDoubsTax, by = "ASV_name")
# View(I6S_topAirTable_final)

# Clean up I6S_topAirTable_final
I6S_topAirTable_final <- I6S_topAirTable_final %>% 
  select(-Kingdom) %>% 
  arrange(desc(meanASV_relAbund_air))
# change NA for non-indicators to be "not indicator"
I6S_topAirTable_final$indictorType[which(is.na(I6S_topAirTable_final$indictorType))] <- "not indicator"
# Make all of these proportions percentages
I6S_topAirTable_final$meanASV_relAbund_air <- I6S_topAirTable_final$meanASV_relAbund_air*100
I6S_topAirTable_final$meanASV_relAbund_soil <- I6S_topAirTable_final$meanASV_relAbund_soil*100
I6S_topAirTable_final$meanASV_relAbund_foliar <- I6S_topAirTable_final$meanASV_relAbund_foliar*100
# View(I6S_topAirTable_final)
# Round off all of these to the nearest tenth
I6S_topAirTable_final$meanASV_relAbund_air <- round(I6S_topAirTable_final$meanASV_relAbund_air, 2) 
I6S_topAirTable_final$meanASV_relAbund_soil <- round(I6S_topAirTable_final$meanASV_relAbund_soil, 2) 
I6S_topAirTable_final$meanASV_relAbund_foliar <- round(I6S_topAirTable_final$meanASV_relAbund_foliar, 2) 
I6S_topAirTable_final$percentOccupancy <- round(I6S_topAirTable_final$percentOccupancy, 1)
I6S_topAirTable_final$SoilPercentOccupancy <- round(I6S_topAirTable_final$SoilPercentOccupancy, 1)
I6S_topAirTable_final$FoliarPercentOccupancy <- round(I6S_topAirTable_final$FoliarPercentOccupancy, 1)

# Make sure that indicator makes sense and it does!
colnames(I6S_topAirTable_final)
#View(I6S_topAirTable_final[,c("ASV_name", "indictorType", "meanASV_relAbund_soil" , "SoilPercentOccupancy", 
#                        "meanASV_relAbund_foliar","FoliarPercentOccupancy" )])

colnames(I6S_topAirTable_final)

# Make a streamlined version for the manuscript supplemental x
I6S_topAirTable_forSupp <- I6S_topAirTable_final[,c(1:3,5,8,11:13,15:16)]
I6S_topAirTable_forSupp$indictorType[which(I6S_topAirTable_forSupp$indictorType == "soilIndicator")] <- "soil"
I6S_topAirTable_forSupp$indictorType[which(I6S_topAirTable_forSupp$indictorType == "foliarIndicator")] <- "foliar"
I6S_topAirTable_forSupp$indictorType[which(I6S_topAirTable_forSupp$indictorType == "not indicator")] <- "not"
I6S_topAirTable_forSupp$ASV_name <- gsub(x=I6S_topAirTable_forSupp$ASV_name, pattern= "ASV_", replacement = "")
I6S_topAirTable_forSupp$ASV_name <- as.character(I6S_topAirTable_forSupp$ASV_name)
# View(I6S_topAirTable_forSupp)

# write.csv(I6S_topAirTable_forSupp, file = "Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_topAirTable_forSupp_09-23-2025.csv") #written September 23, 2025

##############
#  A FEW EXTRA CALCULATIONS FOR MANUSCRIPT
##############
### Double-checking this in manuscript:
# A notable bacterial example is the leaf-specializing genus Methylobacterium-Methylorubrum (Knief et al. 2010, Vacher et al. 2016, Zhang et al. 2024), represented in our bioaerosol 
# samples by two ASVs detected in 74% and 57% of samples and which together comprised an average of 1.7% of reads per sample (Supplemental Material: Table S4).
head(I6S_topAirTable_forSupp)
I6S_topAirTable_forSupp[which(I6S_topAirTable_forSupp$Genus == "Methylobacterium-Methylorubrum"),4] # Shows ASVs detected in 74% and 57% of samples
I6S_topAirTable_forSupp[which(I6S_topAirTable_forSupp$Genus == "Methylobacterium-Methylorubrum"),2] # Shows ASVs mean abundance 0.92 0.77
0.92+0.77 #average of 1.7% per sample

# To fill in these sentences:
# Specifically, the fungal reads in bioaerosol samples consisted of a median of ___% foliar surface- and _____% soil-associated taxa, 
# whereas bacterial reads per bioaerosol sample had medians of ___ foliar surface-associated and _____% soil-associated taxa (Fig.4). 
head(ASVs_all16Sr_noAirSingsDoubs_longer) #made earlier
allSampsIndicReads <- ASVs_all16Sr_noAirSingsDoubs_longer #make a quick copy to edit
# Make non-indicators labeled just to make sure calculations don't get messed up
allSampsIndicReads$indictorType[which(is.na(allSampsIndicReads$indictorType))] <- "notIndicator"
which(is.na(allSampsIndicReads$readCount)) #none are NA, but remove it below in code any
PropsPerAirSampPerIndicType <- allSampsIndicReads %>% 
  filter(str_detect(sampleName, "air_16S")) %>%  #get only air samples
  group_by(sampleName, indictorType) %>% 
  summarize(readsInType = sum(readCount), .groups = "drop_last") %>%
  group_by(sampleName) %>%
  mutate(PropSampIndType = readsInType/ sum(readsInType, na.rm = TRUE)) %>%
  ungroup()
# View(PropsPerAirSampPerIndicType)
# DOUBLE CHECK: Test these calculations for a air_16S_106 (has good amount of reads in all categories)
# 1. What is total number of reads?
air_16S_106_ASVtab <- otu_table(I6SairOnly.ps)[,colnames(otu_table(I6SairOnly.ps)) %in% "air_16S_106"] #get ASV table subset
PropsPerAirSampPerIndicType_106 <-  PropsPerAirSampPerIndicType[which(PropsPerAirSampPerIndicType$sampleName == "air_16S_106"),] #get PropsPerAirSampPerIndicType subset
colSums(air_16S_106_ASVtab) == sum(PropsPerAirSampPerIndicType_106$readsInType) # TRUE and is #5497 !
# 2. Check proportions for each type (1 and 3 are foliar and soil specialists)
colSums(air_16S_106_ASVtab[rownames(air_16S_106_ASVtab) %in% I6S_phylloSpecialists,])/colSums(air_16S_106_ASVtab) == PropsPerAirSampPerIndicType_106$PropSampIndType[1]
colSums(air_16S_106_ASVtab[rownames(air_16S_106_ASVtab) %in% I6S_soilSpecialists,])/colSums(air_16S_106_ASVtab) == PropsPerAirSampPerIndicType_106$PropSampIndType[3]

# Medians and means of foliar and soil in air samples
PropsPerAirSampPerIndicType %>% 
  group_by(indictorType) %>% 
  reframe(meanPropIndicType = mean(PropSampIndType))
# indictorType    meanPropIndicType
# <chr>                       <dbl>
#   1 foliarIndicator            0.333 
# 2 notIndicator               0.612 
# 3 soilIndicator              0.0549
# Added to mnanscript: 31.7% foliar and 2.9% soil
PropsPerAirSampPerIndicType %>% 
  group_by(indictorType) %>% 
  reframe(medianPropIndicType = median(PropSampIndType))
# indictorType    medianPropIndicType
# <chr>                         <dbl>
#   1 foliarIndicator              0.317 
# 2 notIndicator                 0.604 
# 3 soilIndicator                0.0290

# ##############
# # SAME ANALYSIS, BUT WITH PRESENCE/ABSENCE 0F INDICATOR TAXA (I.E. FOR EACH SAMPLE, OUT OF ALL OF THE 
# # ASVS THAT OCCURRED, WHAT PROPRPTION WERE SOIL OR PHYLLOSPHERE INDICATORS (e.g. if a sample has 
# # a richness of 100 ASVs and 10 were fungi (abundance doesn't matter), then would have value of .1))
# ##############
# 
I6S_specialistsPresAbs <- as.data.frame(matrix(data=NA, nrow=ncol(ASVs_all16Sr_noAirSingsDoubs.ps), ncol = 7))
colnames(I6S_specialistsPresAbs) <- c("sampleName", "sampleType", "numberPhylloIndic", "numberSoilIndic", "totalASVsGreaterThanZero", "propNumPhyllo", "propNumSoil")
# Add in sample names 
I6S_specialistsPresAbs[,1] <- rownames(I6Srarefied_metadat) #I6Srarefied_metadat was created from ASVs_all16Sr_noAirSingsDoubs.ps earlier in script 
# Add in sample types (i.e., air, soil, or phyllosphere)
I6S_specialistsPresAbs[,2] <- I6Srarefied_metadat$sampleType

# Before starting, confirm once more that these indices are correct (they are!)
unique(rownames(ASVs_all16Sr_noAirSingsDoubs.ps)[I6S_phylloSpecialistsBigIndex] == I6S_phylloSpecialists)
unique(rownames(ASVs_all16Sr_noAirSingsDoubs.ps)[I6S_soilSpecialistsBigIndex] == I6S_soilSpecialists)

# Fill in the PRESENCE/ABSENCE for each
# This ignores abundance, getting for each sample the count of foliar surface or soil indicator taxa that appear in each sample.
# The proportion is then this number divided by the total number of ASVs (as presence/absence) in each sample
for (k in 1:ncol(ASVs_all16Sr_noAirSingsDoubs.ps)){ #k indexes each sample (columns) and there are 299 samples total
  # Next 2 code lines go only on the rows where there are phyllo or soil specialists, respectively, and give each sample a one if it has at least one count of that ASV.
  # In other words, it sums up total presence/absence of phyllosphere and soil specialists in each sample
  I6S_specialistsPresAbs[k,3] <- sum(ASVs_all16Sr_noAirSingsDoubs.ps[I6S_phylloSpecialistsBigIndex,k] > 0)  #add up all of the times that there is a count more than 0
  I6S_specialistsPresAbs[k,4] <- sum(ASVs_all16Sr_noAirSingsDoubs.ps[I6S_soilSpecialistsBigIndex,k] > 0)  #add up all of the times that there is a count more than 0
  I6S_specialistsPresAbs[k,5] <- sum(ASVs_all16Sr_noAirSingsDoubs.ps[,k] > 0) #adds up all of the total number of non-zero ASVs in each sample k
  I6S_specialistsPresAbs[k,6] <- I6S_specialistsPresAbs[k,3]/I6S_specialistsPresAbs[k,5] #divides number of present foliar surface ASVs by number of total present ASVs in sample
  I6S_specialistsPresAbs[k,7] <- I6S_specialistsPresAbs[k,4]/I6S_specialistsPresAbs[k,5] #divides number of present soil ASVs by number of total present ASVs
}

I6S_specialistsPresAbs %>% 
  group_by(sampleType) %>% 
  summarize(medianNumASVs = median(totalASVsGreaterThanZero))
# sampleType   medianNumASVs
# <chr>                <dbl>
#   1 air                   106.
# 2 phyllosphere          592.
# 3 soil                 1974 

I6S_specialistsPresAbs %>% 
  group_by(sampleType) %>% 
  summarize(sdNumASVs = sd(totalASVsGreaterThanZero)) 
# sampleType   sdNumASVs
# <chr>            <dbl>
#   1 air               82.8
# 2 phyllosphere     385. 
# 3 soil             353. 


# Add in a column that specifies if that air sample is forest or savanna
I6S_specialistsPresAbs_airOnly <- I6S_specialistsPresAbs %>%
  filter(sampleType == "air") %>% #pull out only air samples
  mutate(HabitatAir = ifelse(sampleName %in% I6S_airForestSampNames, "forest", "savanna")) #calls it forest if in I6S_airForestSampNames, otherwise, "savanna"
# View(I6S_specialistsPresAbs_airOnly)

# PLOT IT!
# FINALLY, PLOT (AIR ONLY COMPARISONS)
# SOIL COMPARISON
I6S_airOnlyPresAbs_soil_plot <- ggplot(I6S_specialistsPresAbs_airOnly, aes(x=HabitatAir, y=propNumSoil, fill=HabitatAir)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=c("darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=1, alpha=0.9, height = 0) +
  scale_x_discrete(name= "Underlying vegetation type", labels=c("forest", "savanna")) +
  scale_y_continuous(name= "Number of soil indicator taxa") + # limI6S=c(25, 225) add this argument to change limI6S
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  ggtitle("Soil bacterial indicator taxa across air samples")  +
  theme(plot.title = element_text(size=20)) + #increase plot title size
  theme(legend.position = "none")  #remove legend
# quartz()
I6S_airOnlyPresAbs_soil_plot
# 
# PHYLLOSPHERE COMPARISON
I6S_airOnlyPresAbs_phyllo_plot <- ggplot(I6S_specialistsPresAbs_airOnly, aes(x=HabitatAir, y=propNumPhyllo, fill=HabitatAir)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=c("darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=1, alpha=0.9, height = 0) + #height = zero ensures no vertical jitter
  scale_x_discrete(name= "Underlying vegetation type", labels=c("forest", "savanna")) +
  scale_y_continuous(name= "Number of phyllosphere indicator taxa") + # limI6S=c(25, 225) add this argument to change limI6S
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  ggtitle("Phyllosphere bacterial indicator taxa across air samples")  +
  theme(plot.title = element_text(size=20)) + #increase plot title size
  theme(legend.position = "none")  #remove legend

#quartz()
I6S_airOnlyPresAbs_phyllo_plot


### PRESENCE/ABSENCE: PLOT WITH AIR FOREST AND AIR SAVANNA AND FOLIAR SURFACES AND PHYLLOSPHERE AND SOIL ####
# Savanna air names (object made above)
I6S_airSavSampNames
# Forest air names (object made above)
I6S_airForestSampNames
# Get soil names 
I6S_soilSampNames
# Get foliar surface names 
I6S_foliarSampNames 
# Add in a column that specifies what kind of sample it is
colnames(I6S_specialistsProportion)
I6S_specialistsPresAbs_allTypes <- I6S_specialistsPresAbs %>% 
  mutate(
    sampleTypeWithAir = case_when(
      sampleName %in% I6S_airSavSampNames ~ "SavAirSample",
      sampleName %in% I6S_airForestSampNames ~ "ForAirSample",
      sampleName %in% I6S_foliarSampNames ~ "foliarSample",
      sampleName %in% I6S_soilSampNames ~ "soilSample",
      TRUE ~ NA_character_
    )
  )
colnames(I6S_specialistsPresAbs_allTypes)
# View(I6S_specialistsPresAbs_allTypes)

# Make dataframe "longer" for faceting purposes
colnames(I6S_specialistsPresAbs_allTypes)
I6S_specialistsPresAbs_allTypes_long <- I6S_specialistsPresAbs_allTypes %>%
  pivot_longer(cols = c(propNumPhyllo, propNumSoil), 
               names_to = "proportionType", 
               values_to = "proportionOfTaxaByPresAbs")

# Re-order x-axis elements
unique(I6S_specialistsPresAbs_allTypes_long$sampleTypeWithAir)
I6S_specialistsPresAbs_allTypes_long$sampleTypeWithAir <- 
  factor(I6S_specialistsPresAbs_allTypes_long$sampleTypeWithAir, 
         levels= c("ForAirSample", "SavAirSample", "foliarSample", "soilSample"))

# Vector for new facet labels
facet_labels_PresAbs <- c("propNumPhyllo" = "Foliar surfaces",
                          "propNumSoil" = "Soil")

# Create faceted plot
colnames(I6S_specialistsPresAbs_allTypes_long)
I6S_allTypesPresAbs_2panels <- ggplot(I6S_specialistsPresAbs_allTypes_long, aes(x=sampleTypeWithAir, y=proportionOfTaxaByPresAbs, fill=sampleTypeWithAir)) + 
  geom_boxplot() + 
  geom_jitter(aes(color = sampleTypeWithAir), size=1, alpha=0.9, height = 0) + #height = 0 insures no horizontal jitter
  facet_wrap(~proportionType, scales = "fixed", labeller = labeller(proportionType = facet_labels_PresAbs)) +
  theme_bw() +
  scale_fill_manual(values=c("cornflowerblue", "cornflowerblue", "chartreuse4", "chocolate4")) +
  scale_color_manual(values = c("black", "black", "black", "black")) + #color for jittered points
  scale_x_discrete(labels=c("bioaerosol\nmatrix", "bioaerosol\npatch", "foliar surface", "soil")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16)) +
  #ggtitle("Bacterial indicator taxa across all samples")  +
  labs(y = "Proportion of leaf or soil indicator taxa",
       x = NULL) +
  theme(plot.title = element_text(size=16)) + #increase plot title size
  theme(legend.position = "none") + #remove legend 
  ylim(0, 1.00) + #make y limits between 0 and 1
  theme(strip.text = element_text(size = 16))
I6S_allTypesPresAbs_2panels

# SAVE PLOT (saved October 25, 2025) 
# saveRDS(I6S_allTypesPresAbs_2panels, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_allTypesPresAbs_2panels_10-25-25.rds")

### PRESENCE/ABSENCE: STATISTCS
# Paired Wilcoxon signed-rank test
head(I6S_specialistsPresAbs_allTypes)
dim(I6S_specialistsPresAbs_allTypes)
# Get dataframe 
I6S_specialistsPresAbs_air <- I6S_specialistsPresAbs_allTypes %>% 
  filter(sampleType == "air")
head(I6S_specialistsPresAbs_air)

#  Wilcoxon signed-rank test (paired, one-sided: foliar > soil)
# (wilcox.test ignores exact p when there are ties/large n which is ok!)
wilcoxPresAbs <- wilcox.test(
  I6S_specialistsPresAbs_air$propNumPhyllo, I6S_specialistsPresAbs_air$propNumSoil,
  paired = TRUE, #because samples are not independent
  alternative = "greater",
  exact = FALSE,    #avoids warnings with ties/large n
  conf.int = TRUE,  #gives a CI for the shift (one-sided CI)
  conf.level = 0.95
)
wilcoxPresAbs
# data:  I6S_specialistsPresAbs_air$propNumPhyllo and I6S_specialistsPresAbs_air$propNumSoil
# V = 3570, p-value = 0.0000000000000008705
# alternative hypothesis: true location shift is greater than 0
# 95 percent confidence interval:
#   0.2108007       Inf
# sample estimates:
#   (pseudo)median 
# 0.2226462 


# ##############
# MAKE A VENN DIAGRAM
# ##############
# Get phyllosphere only dataset
I6S_phyllo_rarefied_noControls <- subset_samples(all16Sr_noAirSingsDoubs.ps, sampleType == "phyllosphere")
# Get soil only dataset
I6S_soil_rarefied_noControls <- subset_samples(all16Sr_noAirSingsDoubs.ps, sampleType == "soil")

# Get air ASVs
I6S_airASVs <- ASVs_outta_ps(I6SairOnly.ps)
unique(colSums(I6S_airASVs))
I6S_airASVnames <- names(which(rowSums(I6S_airASVs) > 0)) #get all of the ASVs that have at least one sequence
length(I6S_airASVnames) == 2870+230+479+526 #4105
unique(sort(rowSums(I6S_airASVs[rownames(I6S_airASVs) %in% I6S_airASVnames,]))) #none are zero and no singletons and doubletons

# Get foliar surface ASVs
I6S_phylloASVs <- ASVs_outta_ps(I6S_phyllo_rarefied_noControls)
I6S_phylloASVnames <- names(which(rowSums(I6S_phylloASVs) > 0))
length(I6S_phylloASVnames) #6520

# Get soil ASVs
I6S_soilASVs <- ASVs_outta_ps(I6S_soil_rarefied_noControls)
I6S_soilASVnames <- names(which(rowSums(I6S_soilASVs) > 0))
length(I6S_soilASVnames) #22868
unique(sort(rowSums(I6S_soilASVs[rownames(I6S_soilASVs) %in% I6S_soilASVnames,]))) #none are zero
length(I6S_soilASVnames) == sum(20022, 229, 479, 2138) #22868 and this equals what's reported in Venn diagram


# Make a list with an element for each of the categories of ASVs
I6S_ASVsForVenn <- list(
  I6S_air <- I6S_airASVnames,
  I6S_soil <- I6S_soilASVnames,
  I6S_foliarSurface <- I6S_phylloASVnames
)
str(I6S_ASVsForVenn)

# Within blue circles (air only), add in the percentage that each of these number makes up of the total air community (thinking of all air samples merged as one air sample). 

# Using the vignette as a guide, manually recreate and edit Venn Diagram
vignette("fully-customed", package = "ggVennDiagram")

I6S_venn <- Venn(I6S_ASVsForVenn)
I6S_vennData <- process_data(I6S_venn)
I6S_vennPlot <- ggplot() +
  # 1. region count layer
  geom_polygon(aes(X, Y, fill = id, group = id), 
               data = venn_regionedge(I6S_vennData)) +
  # 2. set edge layer
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(I6S_vennData),
            linewidth = 3.5,
            show.legend = FALSE) +
  # 3. set label layer -- remove labels 
  # geom_text(aes(X, Y, label = name), 
  #       data = venn_setlabel(I6S_vennData)) +
  # 4. region label layer-- changing this from geom_label to geom_text removes label outlines
  geom_text(aes(X, Y, label = count), 
            data = venn_regionlabel(I6S_vennData),
            size=5.5) +
  # 5. Change colors so that all of the filled in areas are white
  scale_color_manual(values = c("cornflowerblue","chocolate4","chartreuse4")) +
  scale_fill_manual(values = c("white", "white", "white","white", "white", "white", "white")) +
  coord_equal() +
  theme_void() + 
  theme(legend.position='none') #remove legend

I6S_vennPlot 
# saveRDS(I6S_vennPlot, file = "RobjectsSaved/I6S_vennPlot_Jan8_2024.rds") #saved January 14, 2024
# saveRDS(I6S_vennPlot, file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_vennPlot_09-26-2025.rds") #saved September 26, 2025

# Double check that the calculations that were manually made in the Venn Diagrams are correct and get ASV names for subsets:
# i. Only soil
# Union gets all of the ASVs that are found in the air and phyllo datasets and sees how this differs
# from those that are found only in the soil
I6S_unique_soil <- setdiff(I6S_soilASVnames, union(I6S_airASVnames, I6S_phylloASVnames))
length(I6S_unique_soil) == 20022 #great, this shows that calculations were correct!
# ii. Only foliar surface
I6S_unique_foliar <- setdiff(I6S_phylloASVnames, union(I6S_airASVnames, I6S_soilASVnames))
length(I6S_unique_foliar) == 3375 #great, this shows that calculations were correct!
# iii. Only in bioaerosols
I6S_unique_air <- setdiff(I6S_airASVnames, union(I6S_soilASVnames, I6S_phylloASVnames))
length(I6S_unique_air) #2869, correct, matches Venn calculations
# iv. Overlap between air and soil
length(which(I6S_airASVnames %in% I6S_soilASVnames==TRUE)) == 229 + 479 #overlap between air and soil
airSoilOverlapNames <- I6S_airASVnames[which(I6S_airASVnames %in% I6S_soilASVnames==TRUE)]
# v. Overlap between air and foliar surfaces
length(which(I6S_airASVnames %in% I6S_phylloASVnames==TRUE)) == 528 + 479 #overlap between air and foliar surfaces
airFoliarOverlapNames <- I6S_airASVnames[which(I6S_airASVnames %in% I6S_phylloASVnames==TRUE)]
# vi. Shared among all!
allOverlapNames <- intersect(airFoliarOverlapNames, I6S_soilASVnames)
length(allOverlapNames) == 479
# vii. Overlap between air and soil NO FOLIAR
airSoilOverLapOnlyNames <- setdiff(airSoilOverlapNames, allOverlapNames) #takes the 229+479 shared by air and soil and subtracts the 479
length(airSoilOverLapOnlyNames) == 229
# viii. Overlap between air and foliar (NO SOIL)
airFoliarOverLapOnlyNames <- setdiff(airFoliarOverlapNames, allOverlapNames) #takes the 528 + 479 shared by air and soil and subtracts the 479
length(airFoliarOverLapOnlyNames) == 528

#####################
# SOME ADDITIONAL CALCULATIONS FOR MANUSCRIPT AND VENN DIAGRAMS (TO ADD IN POWERPOINT OR INKSCAPE)
#####################
# In manuscript: what percentage of total foliar surface ASVs overlapped with soil?
((479+2123)/length(I6S_phylloASVnames))*100 #39.90798, compare with 19.17925 in fungal analysis. Added 40% to manuscript

# Set up of air reads
I6S_airASVnames #these are the ASVs present in air 
2869+229+479+528 == length(I6S_airASVnames)
I6S_airASVs #Made above and gets air ASV table
head(I6S_airASVs)
# For manuscript: Fill in text here: Second, there were some taxa present in the bioaerosols that were not found in either the sampled soils or foliar surfaces, but the contribution of these taxa was minimal
# (xxx% of fungal and xxx% of all bacterial ASVs detected in the air samples
length(I6S_unique_air)/length(I6S_airASVnames)*100 #69.89, so 69.9% added to manuscript 

# To add within blue circles: What percentage do all of these ASVs (2869, 229, 479, and 528)
# make up of the total air community? (Thinking of all air samples merged as one air sample)
# BUILD A DATAFRAME TO GET PERCENTAGE OF EACH SUBSET

head(I6S_airASVs) #this is the ASV table with ASVs as rows and samples as columns
colnames(I6S_airASVs)
# i. Get a longer dataframe with ASV name, sample name, and ASV abundance in that sample
I6S_airASVs_df <- I6S_airASVs %>% 
  rownames_to_column(var="ASVname") %>% 
  pivot_longer(air_16S_1:air_16S_97, names_to = "sampleName", values_to = "ASVabundance")
head(I6S_airASVs_df)
# ii. Grouping by each ASV, get the total abundance of that ASV (across only air samples)
unique(I6S_airASVs_df$sampleName)
I6SASVabund_AirSamps <- I6S_airASVs_df %>% 
  group_by(ASVname) %>% 
  summarize(ASVtotalAbund = sum(ASVabundance))
sum(I6SASVabund_AirSamps$ASVtotalAbund) #461273
# iii. Grouping by each ASV, proportion that each ASV makes up of total community is that ASV total 
# abundance divided by the total ASV abundance of all reads in the air samples
I6SASVabund_AirSamps_2 <- I6SASVabund_AirSamps %>% 
  group_by(ASVname) %>% 
  mutate(propASVabund = ASVtotalAbund/461273)
# iv. For the ASVs that are only in air (2869 ASVs), what percentage is this of total reads?
# This code pulls out the proportional abundances of all of the ASVs unique to the air and then 
# adds these together
sumOf_16SuniqueAir <- sum(I6SASVabund_AirSamps_2$propASVabund[I6SASVabund_AirSamps_2$ASVname %in% I6S_unique_air])
sumOf_16SuniqueAir #0.5033006, so 50.3% goes below 2869, i.e., bioaerosol-only taxa
# v. For the ASVs that are shared between soil and air (229 ASVs), what percentage is this of total reads?
sumOf_16SairSoil <- sum(I6SASVabund_AirSamps_2$propASVabund[I6SASVabund_AirSamps_2$ASVname %in% airSoilOverLapOnlyNames])
sumOf_16SairSoil #0.0472345, so 4.7% goes below 229
# vi. For the ASVs that are shared between phyllo and air (528 ASVs), what percentage is this of total reads?
sumOf_16SairPhyllo <- sum(I6SASVabund_AirSamps_2$propASVabund[I6SASVabund_AirSamps_2$ASVname %in% airFoliarOverLapOnlyNames])
sumOf_16SairPhyllo  #0.2594147%, so 25.9% goes below 528, i.e., bioaerosol and foliar surface taxa
# vii. Of the ASVs that are shared among all taxa, what proportion are these among air samples?
sumOf_16Sall <- sum(I6SASVabund_AirSamps_2$propASVabund[I6SASVabund_AirSamps_2$ASVname %in% allOverlapNames])
sumOf_16Sall #0.1900501, so 19.0% goes below 479
# DO THESE ADD UP TO 1, AS THEY SHOULD?
sum(sumOf_16SuniqueAir, sumOf_16SairSoil, sumOf_16SairPhyllo, sumOf_16Sall) #YES!

#### GET PROPORTIONS OF READS OF TAXA THAT ARE ONLY FOUND IN BIOAEROSOL SAMPLES IN BIOAEROSOLS #####
# IS CORRECT BUT GRAYED OUT SINCE NO LONGER USING

# I6S_air_taxa <- as.data.frame(as.matrix(tax_table(I6SairOnly.ps)))
# taxaOnlyInAir <- I6S_air_taxa[rownames(I6S_air_taxa) %in% I6S_unique_air,]
# # View(taxaOnlyInAir)
# taxaOnlyInAir <- taxaOnlyInAir %>% 
#   rownames_to_column(var="ASVname")
# # Confirm that neither are present here-- both are false as they should be 
# unique(taxaOnlyInAir$ASVname %in% I6S_soilASVnames)
# unique(taxaOnlyInAir$ASVname %in% I6S_phylloASVnames)
# # Add in mean relative abundance in air samples
# unique(taxaOnlyInAir$ASVname %in% I6S_airASVs_df$ASVname)
# # I6S_airASVs_df is the abundance of ASVs across all air samples
# I6S_airASVs_df_2 <- I6S_airASVs_df %>%  #first, add in a column of ASV abundance summed by sample (this matches after rarefying!)
#   group_by(sampleName) %>% 
#   mutate(totalASVabund= sum(ASVabundance)) %>% 
#   mutate(relAbundASV = ASVabundance/totalASVabund) #get relative abundance of each ASV
# # View(I6S_airASVs_df_2)
# colnames(I6S_airASVs_df_2)
# # Get mean, min, and max, then remove extra columns that are duplicates (need to do mutate first or summarize would make min and max
# # impossible to get after getting mean)
# allI6S_airASVsRelAbundStats <- I6S_airASVs_df_2 %>% 
#   group_by(ASVname) %>% 
#   mutate(meanRelAbund = mean(relAbundASV)) %>% #mutate all of these below so that can do additional summarizing of them
#   mutate(minRelAbund = min(relAbundASV)) %>% 
#   mutate(maxRelAbund = max(relAbundASV)) %>% 
#   select(-sampleName) %>%  #remove sampleName 
#   select(-ASVabundance) %>% 
#   select(-totalASVabund) %>% 
#   select(-relAbundASV)
# # View(allI6S_airASVsRelAbundStats)
# allI6S_airASVsRelAbundStats <- allI6S_airASVsRelAbundStats %>% 
#   distinct() #remove duplicate rows caused by initial keeping of sample names and associated sample info
# # View(allI6S_airASVsRelAbundStats)
# # MERGE WITH taxaOnlyInAir for INFO ON THESE TAXA
# taxaOnlyInAirFinal <- left_join(taxaOnlyInAir, allI6S_airASVsRelAbundStats, by= "ASVname")
# taxaOnlyInAirFinal <- taxaOnlyInAirFinal %>% 
#   arrange(desc(meanRelAbund))
# # View(taxaOnlyInAirFinal) -- these calculations also match those in Table S4, a good confirmation of results
# # write.csv(x= taxaOnlyInAirFinal, file = "taxaOnlyInAirFinal.csv") #may 23, 2025

#### ANOTHER WAY OF CHECKING PERCENTAGES IN EACH VENN DIAGRAM PART (MATCHES CALCUALTIONS ABOVE SO IS A GOOD CHECK!!)
# ## 2. What percentage of total reads are the 230 ASVs that are shared between soil and air (but not foliar) of the soil?
# # These are the ASVs that are in air and soil
# length(which(I6S_airASVnames %in% I6S_soilASVnames==TRUE)) == 230 + 479 #overlap between air and soil
# # Get the names of these ASVs
# I6S_airSoilSharedNames <- I6S_airASVnames[which(I6S_airASVnames %in% I6S_soilASVnames)]
# # Remove phyllosphere ASVs in this to get only these 169. So these are ones in air and soil, but not phyllosphere
# I6S_airSoilNOphylloNames <- I6S_airSoilSharedNames[which(I6S_airSoilSharedNames %in% I6S_phylloASVnames ==FALSE)]
# length(I6S_airSoilNOphylloNames) == 230 #matches Venn diagram 
# # Get proportion that these ASVs make up of all air ASVs
# sumOf_16S_airSoilnoPhyllo <- sum(I6SASVabund_AirSamps_2$propASVabund[I6SASVabund_AirSamps_2$ASVname %in% I6S_airSoilNOphylloNames])
# sumOf_16S_airSoilnoPhyllo #0.03750057
# 
# ## 3. For the 526 ASVs that are in air and phyllosphere, but not the soil.
# # These are the ASVs that are in air and phyllosphere
# length(which(I6S_airASVnames %in% I6S_phylloASVnames==TRUE)) == 526 + 479 #overlap between air and foliar surfaces
# # Get the names of these ASVs
# I6S_airPhylloSharedNames <- I6S_airASVnames[which(I6S_airASVnames %in% I6S_phylloASVnames)]
# length(I6S_airPhylloSharedNames) ==526 + 479 #matches Venn Diagram
# # Remove soil ASVs in this to get only these 526. So these are ones in air and phyllo, but not soil
# I6S_airPhylloNOsoilNames <- I6S_airPhylloSharedNames[which(I6S_airPhylloSharedNames %in% I6S_soilASVnames ==FALSE)]
# length(I6S_airPhylloNOsoilNames) == 526
# # Get proportion that these ASVs make up of all air ASVs
# sumOf_16S_airPhylloNoSoil <- sum(I6SASVabund_AirSamps_2$propASVabund[I6SASVabund_AirSamps_2$ASVname %in% I6S_airPhylloNOsoilNames])
# sumOf_16S_airPhylloNoSoil #0.2592283
# 
# ## 3. For the 469 ASVs that were found across all sample types
# # Strt with the ASVs that are found in phyllo and air and then get which of these are also in the soil
# I6S_airPhylloAndSoilNames <- I6S_airPhylloSharedNames[which(I6S_airPhylloSharedNames %in% I6S_soilASVnames ==TRUE)]
# length(I6S_airPhylloAndSoilNames) == 479
# # Get proportion that these ASVs make up of all air ASVs
# sumOf_16S_airPhylloAndSoil <- sum(I6SASVabund_AirSamps_2$propASVabund[I6SASVabund_AirSamps_2$ASVname %in% I6S_airPhylloAndSoilNames])
# sumOf_16S_airPhylloAndSoil #0.1982882
# 
# # Check that all of these add up to one!
# sumOf_16S_airPhylloAndSoil + sumOf_16S_airPhylloNoSoil + sumOf_16S_airSoilnoPhyllo + sumOf_16SuniqueAir
# 
