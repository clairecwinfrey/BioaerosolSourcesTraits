# ITS_souceTracking_Sept.R
# 
# This version was started Sept 28, 2025
# 
# Earlier version of this script was ITS_souceTracking_final.R (August 6, 2024). 

# Script description: Aims to explore where fungal reads in the air are coming from (i.e., from soil or phyllosphere sources)

# 1. Indicator analyses (indcator species analysis to test for different taxa between soil and phyllosphere):
## a. Plots: plots showing proportions of soil and phyllo indicator species in air, soil, and phyllosphere samples
### i. Two paneled plot, with the top plot showing soil indicator taxa across air (matrix and patch combined),
# phyllosphere, and soil samples, and bottom phyllosphere indicator taxa across air, phyllosphere, and soil samples
### ii. Proportions of only air samples, but separated out by forest/savanna
### iii. NEED TO ADD: PROPORTION IN EACH HABITAT BASED ON 
## b. Analyses (based on initial indicator species analysis results):
#### i.  Separate Wilcoxon signed-rank tests to look for differences in the proportion of foliar surface or soil 
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

# HANDY FUNCTION TO GET ASV TABLE OUT OF PHYLOSEQ:
# Function to get ASV table out of phyloseq so that we can #View it better
# (Inspired by similar function at https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  return(as.data.frame(ASVTable))
}

# LOAD DATA (made/saved in fullITS_EDArarefied_Part2_Sept.R)
allITSr_noAirSingsDoubs.ps <- readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/sallITSr_noAirSingsDoubs_ps.rds") #yes, this is correct version, even with "sallITS" typo!
allITSr_noAirSingsDoubs.ps
sample_data(allITSr_noAirSingsDoubs.ps)
length(which(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType == "soil")) #155 samples
length(which(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType == "air")) #110 samples
length(which(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType == "phyllosphere")) #59 samples
plantDat_df <- sample_data(allITSr_noAirSingsDoubs.ps)[which(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType == "phyllosphere"),]
airSampsColsIndex <- grepl(x=colnames(otu_table(allITSr_noAirSingsDoubs.ps)), pattern= "air_")
# Confirms that this has no bioaerosol singletons or doubletons
unique(sort(rowSums(otu_table(allITSr_noAirSingsDoubs.ps)[,airSampsColsIndex])))

table(plantDat_df$EU)
# EU_52 EU_53S EU_54S   EU_8 
# 17     12     15     15 
table(plantDat_df$EU, plantDat_df$PlantSpecies)

############################################
# I. SOIL AND PHYLLOSPHERE ONLY INDICATOR SPECIES ANALYSIS 
############################################
# For this, I perform an indicator species analysis on only soil and phyllosphere samples

##### 1. PERFORM INDICATOR SPECIES ANALYSIS #####
# This will perform an indicator species analysis to determine which are soil and phyllosphere specialists

# Make phyloseq object with only phyllosphere and soil samples
phylloSoilonly_ITS.ps <- subset_samples(allITSr_noAirSingsDoubs.ps, sampleType!= "air")

# Get the metadata for these samples
phylloSoilonly_ITS_meta <- as.data.frame(as.matrix(sample_data(phylloSoilonly_ITS.ps)))
unique(phylloSoilonly_ITS_meta$HabitatAir) #only phyllosphere and soil
# Get a string that has the soil versus phyllosphere designation
soilPhyllos_ITS <- phylloSoilonly_ITS_meta$HabitatAir

# Get ASV and taxonomy tables out of phyloseq:
phylloSoilonly_ITS_ASVs <- t(ASVs_outta_ps(phylloSoilonly_ITS.ps))
# phylloSoilonly_ITS_ASVs needs to be columns as ASVs and rows as samples
if (nrow(phylloSoilonly_ITS_ASVs)> ncol(phylloSoilonly_ITS_ASVs)) {
  phylloSoilonly_ITS_ASVs <- t(phylloSoilonly_ITS_ASVs)
} else {
  print("As needed for multipatt, ASVs are columns!")
}
phylloSoilonly_ITS_ASVs 

phylloSoilonly_ITS_tax <- as.data.frame(phyloseq::tax_table(phylloSoilonly_ITS.ps), stringsAsFactors = F)
head(phylloSoilonly_ITS_tax) #looks good!

# Before running below. Because this is true, then the ASV table has same order of samples as soilPhyllos_ITS, which was built from phylloSoilonly_ITS_meta
unique(rownames(phylloSoilonly_ITS_meta) == rownames(phylloSoilonly_ITS_ASVs))

# Perform indicator analysis --commented out below because if it runs it takes forever
set.seed(9)
# ITS_phylloSoilIndic <- multipatt(x=phylloSoilonly_ITS_ASVs, cluster=soilPhyllos_ITS, func="r.g", duleg = TRUE, control=how(nperm = 9999)) #duleg=true means that combinations
# are not considered
summary(ITS_phylloSoilIndic) #shows all of the taxon IDs for ASVs associated with each group, as well as
# the  r.g. value (under stat, this is the correlation index that takes into account the
# variation within and among groups), and the p-value from a permutation test.

# re-run and saved on Aug. 5, 2024 (on server)
# save(ITS_phylloSoilIndic, file = "RobjectsSaved/ITS_phylloSoilIndic_Aug5_2024")
# re-run and re-saved September 28, 2025 on own computer
# save(ITS_phylloSoilIndic, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/ITS_phylloSoilIndic_Sept_28_2025")

# HERE, I JUST LOAD IT IN CASE THE CODE ABOVE HASN'T BEEN RUN (WHICH DOES TAKE A WHILE)
load(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/ITS_phylloSoilIndic_Sept_28_2025")
summary(ITS_phylloSoilIndic)

##### 2. CHECK OUT AND RE-FORMAT RESULTS ##### 
# View(ITS_phylloSoilIndic$sign)
ITS_phylloSoilIndic_Results <- ITS_phylloSoilIndic$sign
ITS_phylloSoilIndic_.05 <- ITS_phylloSoilIndic_Results %>% 
  filter(p.value < 0.05) #filter out ASVs that had a sig rating of 0.05 or greater
dim(ITS_phylloSoilIndic_.05) #3710 ASVs left
# View(ITS_phylloSoilIndic_.05)

# ASSIGN SOIL OR PHYLLOSPHERE BASED ON INDEX
ITS_phylloSoilIndic_.05$group <- NA #add in another column, group
for (j in 1:nrow(ITS_phylloSoilIndic_.05)){
  if (ITS_phylloSoilIndic_.05$index[j] == "1") {
    ITS_phylloSoilIndic_.05$group[j] <- "phyllosphere" 
  } else {
    ITS_phylloSoilIndic_.05$group[j] <- "soil" 
  }
}

# GET INDICES AND NAMES OF PHYLLOSPHERE AND SOIL SPECIALISTS
ITS_soilSpecialistsIndex <- which(ITS_phylloSoilIndic_.05$group == "soil")
length(ITS_soilSpecialistsIndex) #755 
ITS_soilSpecialists <- rownames(ITS_phylloSoilIndic_.05)[ITS_soilSpecialistsIndex]
ITS_phyllospecialistsIndex <- which(ITS_phylloSoilIndic_.05$group == "phyllosphere") 
length(ITS_phyllospecialistsIndex) #2955
ITS_phylloSpecialists <- rownames(ITS_phylloSoilIndic_.05)[ITS_phyllospecialistsIndex]

# Make two plots that have the % of indicator taxa in each group (soil, air, and phyllo), where
# first plot is % soil and second plot is % phyllo for each of these groups.
phylloSoilonly_ITS_tax <- as.data.frame(phyloseq::tax_table(phylloSoilonly_ITS.ps), stringsAsFactors = F)
head(phylloSoilonly_ITS_tax)

ASVs_allITSr_noAirSingsDoubs.ps <- ASVs_outta_ps(allITSr_noAirSingsDoubs.ps)
# Get how to pull the specialists out of the big dataframe
# Soil
ITS_soilSpecialistsBigIndex <- which(rownames(ASVs_allITSr_noAirSingsDoubs.ps) %in% ITS_soilSpecialists == TRUE)
# Next line shows that by using "BigIndex", I get these same taxa. This proof shows that the for loop later will work.
setdiff((rownames(ASVs_allITSr_noAirSingsDoubs.ps[ITS_soilSpecialistsBigIndex,])), ITS_soilSpecialists) #these are the same, as expected! 
# Phyllo
ITS_phylloSpecialistsBigIndex <- which(rownames(ASVs_allITSr_noAirSingsDoubs.ps) %in% ITS_phylloSpecialists == TRUE)
# Next line shows that by using "BigIndex", I get these same taxa. This proof shows that the for loop later will work.
setdiff((rownames(ASVs_allITSr_noAirSingsDoubs.ps[ITS_phylloSpecialistsBigIndex,])), ITS_phylloSpecialists) #these are the same, as expected! 

# Pull out metadata from phyloseq object
ITSrarefied_metadat <- as.data.frame(as.matrix(sample_data(allITSr_noAirSingsDoubs.ps)))
colnames(ITSrarefied_metadat)
unique(rownames(ITSrarefied_metadat) == colnames(ASVs_allITSr_noAirSingsDoubs.ps)) #can keep other stuff, since order is same as in ITSrarefied_meta

# PRE-ALLOCATE AND CHECK to store the proportions of specialists in each sample (as many rows as samples)
ITS_specialistsProportion <- as.data.frame(matrix(data=NA, nrow=ncol(ASVs_allITSr_noAirSingsDoubs.ps), ncol = 5))
colnames(ITS_specialistsProportion) <- c("sampleName", "sampleType", "totalNumbReads", "phylloSpecProp", "soilSpecProp")
ITS_specialistsProportion[,1] <- rownames(ITSrarefied_metadat) #column one gets all sample names
ITS_specialistsProportion[,2] <- ITSrarefied_metadat$sampleType #column 2 gets sampleType
unique(colnames(ASVs_allITSr_noAirSingsDoubs.ps) == ITS_specialistsProportion[,1]) # shows that can get the number of reads in each sample as below
ITS_specialistsProportion$totalNumbReads <- colSums(ASVs_allITSr_noAirSingsDoubs.ps)
head(ITS_specialistsProportion)
sort(unique(colSums(ASVs_allITSr_noAirSingsDoubs.ps))) #this occasionally varies a tiny bit from the rarefied number, 8500, because of removing singletons and doubletons

# FILL IN PROPORTIONS OF EACH TYPE OF SPECIALIST IN EACH SAMPLE
for (k in 1:ncol(ASVs_allITSr_noAirSingsDoubs.ps)){ #k indexes each sample (324 )
  # Add up all of the reads for soil specialists in the k-th column (i.e. kth sample) and divide it by total number of reads in that kth sample.And place it in the kth row for soil specialists 
  ITS_specialistsProportion[k,5] <- sum(ASVs_allITSr_noAirSingsDoubs.ps[ITS_soilSpecialistsBigIndex,k])/ITS_specialistsProportion$totalNumbReads[k] 
  # Do the same for the phyllosphere specialists
  ITS_specialistsProportion[k,4] <- sum(ASVs_allITSr_noAirSingsDoubs.ps[ITS_phylloSpecialistsBigIndex,k])/ITS_specialistsProportion$totalNumbReads[k]
}

# View(ITS_specialistsProportion)
# re-saved on server August 6, 2024:
# write_rds(ITS_specialistsProportion, file="RobjectsSaved/ITS_specialistsProportion_Aug6_2024.Rdata")
# saved Oct. 15, 2023
# save(ITS_specialistsProportion, file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_specialistsProportion_Oct15")

# MOST RECENT: Saved September 28, 2025 on own computer
# save(ITS_specialistsProportion, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/ITS_specialistsProportion_Sept_28_2025")

##############
# II. PLOTS SHOWING PROPORTIONS OF SOIL AND PHYLLO INDICATOR SPECIES IN AIR, SOIL, AND PHYLLOSPHERE SAMPLES
##############
soilSpecProp_plot_ITS <- ggplot(ITS_specialistsProportion, aes(x=sampleType, y=soilSpecProp, fill=sampleType)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=c("cornflowerblue", "chartreuse4", "chocolate4")) +
  geom_jitter(color="black", size=1, alpha=0.9, height = 0) + #height = 0 means no horizontal jitter which would be very misleading!
  scale_x_discrete(name= "Sample Type", labels=c("air", "phyllosphere", "soil")) +
  scale_y_continuous(name= "Proportion of soil indicator taxa") + # limits=c(25, 225) add this argument to change limits
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  ggtitle("Soil indicator taxa across samples")  +
  theme(plot.title = element_text(size=20)) + #increase plot title size
  theme(legend.position = "none")  #remove legend 

# ggplot it!
phylloSpecProp_plot_ITS <- ggplot(ITS_specialistsProportion, aes(x=sampleType, y=phylloSpecProp, fill=sampleType)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=c("cornflowerblue", "chartreuse4", "chocolate4")) +
  geom_jitter(color="black", size=1, alpha=0.9, height = 0) +
  scale_x_discrete(name= "Sample Type", labels=c("air", "phyllosphere", "soil")) +
  scale_y_continuous(name= "Proportion of phyllosphere indicator taxa") + # limITS=c(25, 225) add this argument to change limITS
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  ggtitle("Phyllosphere indicator taxa across samples")  +
  theme(plot.title = element_text(size=20)) + #increase plot title size
  theme(legend.position = "none")  #remove legend 

# quartz()
grid.arrange(soilSpecProp_plot_ITS, phylloSpecProp_plot_ITS, ncol=1)

##############
# PLOTS SHOWING PROPORTIONS OF ONLY AIR SAMPLES, BUT SEPARATED OUT BY FOREST/SAVANNA
##############
#View(ITS_specialistsProportion)
# View(ITSrarefied_metadat) #look at the sample meta data. HabitatAir is what we need to divide the samples up!

# Get savanna air names
ITS_airSavSampNames <- rownames(ITSrarefied_metadat %>%
                                  filter(sampleType == "air" & HabitatAir == "savanna"))

# Get forest air names
ITS_airForestSampNames <- rownames(ITSrarefied_metadat %>%
                                     filter(sampleType == "air" & HabitatAir == "forest"))

intersect(ITS_airSavSampNames, ITS_airForestSampNames) #these are different, as they should be!

# Add in a column that specifies if that air sample is forest or savanna
ITS_specialistsProportion_airOnly <- ITS_specialistsProportion %>%
  filter(sampleType == "air") %>% #pull out only air samples
  mutate(HabitatAir = ifelse(sampleName %in% ITS_airForestSampNames, "forest", "savanna")) #calls it forest if in ITS_airForestSampNames, otherwise, "savanna"

# Double check that this looks good:
savNamesCheck_ITS <- ITS_specialistsProportion_airOnly %>%
  filter(HabitatAir == "savanna") %>%
  select(sampleName)

unique(sort(savNamesCheck_ITS$sampleName) == sort(ITS_airSavSampNames)) #looks good!

# FINALLY, PLOT THEM!
# SOIL COMPARISON
airOnlysoilSpecProp_plot_ITS <- ggplot(ITS_specialistsProportion_airOnly, aes(x=HabitatAir, y=soilSpecProp, fill=HabitatAir)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values=c("darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=1, alpha=0.9, height = 0) +
  scale_x_discrete(name= NULL, labels=c("forested\nmatrix", "open\npatch")) +
  scale_y_continuous(name= "Proportion of reads from\nsoil indicator taxa", limits=c(0.00, 1.00)) + # limITS=c(25, 225) add this argument to change limITS
  theme(axis.text=element_text(size=16, color = "black"),
        axis.title=element_text(size=16, color = "black")) +
  ggtitle("Fungi")  +
  theme(plot.title = element_text(size=16, color = "black")) + #increase plot title size
  theme(legend.position = "none")  #remove legend

# quartz()
airOnlysoilSpecProp_plot_ITS

# Saved Nov. 4, 2025 (made for postdoc talk/exit talk)
# saveRDS(airOnlysoilSpecProp_plot_ITS, file = "~/Desktop/PostCUJobSearch/ResearchTalk/Figures/airOnlysoilSpecProp_plot_ITS.rds")
airOnlysoilSpecProp_plot_ITS <- readRDS("~/Desktop/PostCUJobSearch/ResearchTalk/Figures/airOnlysoilSpecProp_plot_ITS.rds")

# PHYLLOSPHERE COMPARISON
airOnlyphylloSpecProp_plot_ITS <- ggplot(ITS_specialistsProportion_airOnly, aes(x=HabitatAir, y=phylloSpecProp, fill=HabitatAir)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values=c("darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=1, alpha=0.9, height = 0) +
  scale_x_discrete(name= NULL, labels=c("forested\nmatrix", "open\npatch")) +
  scale_y_continuous(name= "Proportion of reads from\nleaf indicator taxa", limits=c(0.00, 1.00)) + # limITS=c(25, 225) add this argument to change limITS
  theme(axis.text=element_text(size=16, color = "black"),
        axis.title=element_text(size=16, color = "black")) +
  ggtitle("Fungi")  +
  theme(plot.title = element_text(size=16, color = "black")) + #increase plot title size
  theme(legend.position = "none")  #remove legend

# quartz()
airOnlyphylloSpecProp_plot_ITS

# Saved Nov. 4, 2025 (made for postdoc talk/exit talk)
# saveRDS(airOnlyphylloSpecProp_plot_ITS, file = "~/Desktop/PostCUJobSearch/ResearchTalk/Figures/airOnlyphylloSpecProp_plot_ITS.rds")
airOnlyphylloSpecProp_plot_ITS <- readRDS(file = "~/Desktop/PostCUJobSearch/ResearchTalk/Figures/airOnlyphylloSpecProp_plot_ITS.rds")
airOnlyphylloSpecProp_plot_ITS

# #### MAKE A TWO PANELED PLOT ###
# # Re-shape dataframe to make proportions all in there!
# colnames(ITS_specialistsProportion_airOnly)
# ITS_specialistsProportion_air_long <- ITS_specialistsProportion_airOnly %>%
#   pivot_longer(cols = c(phylloSpecProp, soilSpecProp),
#                names_to = "proportionType",
#                values_to = "proportionOfReads")
# #View(ITS_specialistsProportion_air_long)
#
# # Vector for new facet labels
# facet_labels <- c("phylloSpecProp" = "Leaf surfaces",
#                   "soilSpecProp" = "Soil")
#
# # Create faceted plot
# ITS_airOnly_phylloSoil_2panels <- ggplot(ITS_specialistsProportion_air_long, aes(x=HabitatAir, y=proportionOfReads, fill=HabitatAir)) +
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
# ITS_airOnly_phylloSoil_2panels #exported as /RobjectsSaved/ITS_airOnly_phylloSoil_2panels.pdf
# saveRDS(ITS_airOnly_phylloSoil_2panels, file = "RobjectsSaved/ITS_airOnly_phylloSoil_2panels_Jan8_2024.rds") #saved January 8, 2024


#### SAVE ALL PLOTS ####
# ON SERVER: saved August 6, 2024:
# write_rds(phylloSpecProp_plot_ITS, file = "RobjectsSaved/phylloSpecProp_plot_ITS_Aug6_2024")
# write_rds(soilSpecProp_plot_ITS, file = "RobjectsSaved/soilSpecProp_plot_ITS_Aug6_2024")
# write_rds(airOnlyphylloSpecProp_plot_ITS, file = "RobjectsSaved/airOnlyphylloSpecProp_plot_ITS_Aug6_2024")
# write_rds(airOnlysoilSpecProp_plot_ITS, file = "RobjectsSaved/airOnlysoilSpecProp_plot_ITS_Aug6_2024")

# ON OWN COMPUTER: saved Oct 15, 2023
# save(airOnlyphylloSpecProp_plot_ITS, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airOnlyphylloSpecProp_plot_ITS_Oct15_2023")
# save(airOnlysoilSpecProp_plot_ITS, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airOnlysoilSpecProp_plot_ITS_Oct15_2023")

### PLOT WITH AIR FOREST AND AIR SAVANNA AND FOLIAR SURFACES AND PHYLLOSPHERE AND SOIL ####
# Savanna air names (object made above)
ITS_airSavSampNames
# Forest air names (object made above)
ITS_airForestSampNames
# Get soil names 
ITS_soilSampNames <- rownames(ITSrarefied_metadat %>% 
                                filter(sampleType == "soil"))
# Get foliar surface names 
ITS_foliarSampNames <- rownames(ITSrarefied_metadat %>% 
                                  filter(sampleType == "phyllosphere"))

# Add in a column that specifies what kind of sample it is
colnames(ITS_specialistsProportion)
ITS_specialistsProportion_allTypes <- ITS_specialistsProportion %>% 
  mutate(
    sampleTypeWithAir = case_when(
      sampleName %in% ITS_airSavSampNames ~ "SavAirSample",
      sampleName %in% ITS_airForestSampNames ~ "ForAirSample",
      sampleName %in% ITS_foliarSampNames ~ "foliarSample",
      sampleName %in% ITS_soilSampNames ~ "soilSample",
      TRUE ~ NA_character_
    )
  )
colnames(ITS_specialistsProportion_allTypes)

# Make dataframe "longer" for faceting purposes
colnames(ITS_specialistsProportion_allTypes)
ITS_specialistsProportion_allTypes_long <- ITS_specialistsProportion_allTypes %>%
  pivot_longer(cols = c(phylloSpecProp, soilSpecProp), 
               names_to = "proportionType", 
               values_to = "proportionOfReads")

# Vector for new facet labels
facet_labels <- c("phylloSpecProp" = "Foliar surfaces", 
                  "soilSpecProp" = "Soil")

# Re-order x-axis elements
unique(ITS_specialistsProportion_allTypes_long$sampleTypeWithAir)
ITS_specialistsProportion_allTypes_long$sampleTypeWithAir <- 
  factor(ITS_specialistsProportion_allTypes_long$sampleTypeWithAir, 
         levels= c("ForAirSample", "SavAirSample", "foliarSample", "soilSample"))

# Saved May 20, 2025
# saveRDS(ITS_specialistsProportion_allTypes_long, file="RobjectsSaved/ITS_specialistsProportion_allTypes_long.RData")
# MOST RECENT: Saved September 22, 2025 on own computer
# save(ITS_specialistsProportion_allTypes_long, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/ITS_specialistsProportion_allTypes_long_Sept_22_2025")

# Create faceted plot
ITS_allTypesAir_2panels <- ggplot(ITS_specialistsProportion_allTypes_long, aes(x=sampleTypeWithAir, y=proportionOfReads, fill=sampleTypeWithAir)) + 
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
ITS_allTypesAir_2panels

# SAVE PLOT (saved Oct 25, 2025) 
# saveRDS(ITS_allTypesAir_2panels, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_allTypesAir_2panels_ITS_09-28-25.rds")

##############
# STATISTICS

# 1. First, test if air samples collected open patch or forested matrix differed in their proportions of foliar surface 
# indicator taxa or soil indicator taxa.
# Reported APA style as here:
# https://about.illinoisstate.edu/jhkahn/apastats/#:~:text=Chi%2DSquare%20statistics%20are%20reported,35.

# Re-shape dataframe to make proportions all in there!
colnames(ITS_specialistsProportion_airOnly)
ITS_specialistsProportion_air_long <- ITS_specialistsProportion_airOnly %>%
  pivot_longer(cols = c(phylloSpecProp, soilSpecProp),
               names_to = "proportionType",
               values_to = "proportionOfReads")
# View(ITS_specialistsProportion_air_long)

# Vector for new facet labels
facet_labels <- c("phylloSpecProp" = "foliar surfaces",
                  "soilSpecProp" = "Soil")

head(ITS_specialistsProportion_air_long) #this has information on habitat type too
# Get counts of numbers of soil and foliar surface indicators, by 
# multiplying the proportion of each sample's foliar and soil indicator count numbers by the total # of reads
ITS_specialistsProportion_air_long$ReadsByHab <- ITS_specialistsProportion_air_long$proportionOfReads* ITS_specialistsProportion_air_long$totalNumbReads
head(ITS_specialistsProportion_air_long)
# Since these number of samples are so different, will need to do a proportions test
length(which(ITS_specialistsProportion_air_long$HabitatAir== "savanna")) #102
length(which(ITS_specialistsProportion_air_long$HabitatAir== "forest")) #118

# Do a few checks to ensure that ITS_specialistsProportion_air_long is correct
# Try for two samples: "air_ITS_103" and "air_ITS_66"
head(ITS_specialistsProportion_air_long)
# Make smaller, subsetted ASV table
ASVtab_103_66 <- otu_table(allITSr_noAirSingsDoubs.ps)[,colnames(otu_table(allITSr_noAirSingsDoubs.ps)) %in% c("air_ITS_103", "air_ITS_66")]
# Get sums of reads in each sample
ASVSums_103_66 <- colSums(ASVtab_103_66) 
# Get number of reads from foliar surface and soil indicators
foliar103_ASVs <- ASVtab_103_66[rownames(ASVtab_103_66) %in% ITS_phylloSpecialists,1]
foliar66_ASVs <- ASVtab_103_66[rownames(ASVtab_103_66) %in% ITS_phylloSpecialists,2]
soil103_ASVs <- ASVtab_103_66[rownames(ASVtab_103_66) %in% ITS_soilSpecialists,1]
soil66_ASVs <- ASVtab_103_66[rownames(ASVtab_103_66) %in% ITS_soilSpecialists,2]
# Get proportions
foliarProp_103 <- colSums(foliar103_ASVs)/ASVSums_103_66[1]
foliarProp_66 <- colSums(foliar66_ASVs)/ASVSums_103_66[2]
soilProp_103 <- colSums(soil103_ASVs)/ASVSums_103_66[1]
soilProp_66 <- colSums(soil66_ASVs)/ASVSums_103_66[2]
# Double check these calculations above against those in larger dataframe -- ALL TRUE!
# For air_ITS_103
ITS_specialistsProportion_air_long[which(ITS_specialistsProportion_air_long$sampleName %in% "air_ITS_103" & ITS_specialistsProportion_air_long$proportionType == "phylloSpecProp"),6] == foliarProp_103
ITS_specialistsProportion_air_long[which(ITS_specialistsProportion_air_long$sampleName %in% "air_ITS_103" & ITS_specialistsProportion_air_long$proportionType == "soilSpecProp"),6] == soilProp_103
# For air_ITS_66
ITS_specialistsProportion_air_long[which(ITS_specialistsProportion_air_long$sampleName %in% "air_ITS_66" & ITS_specialistsProportion_air_long$proportionType == "phylloSpecProp"),6] == foliarProp_66
ITS_specialistsProportion_air_long[which(ITS_specialistsProportion_air_long$sampleName %in% "air_ITS_66" & ITS_specialistsProportion_air_long$proportionType == "soilSpecProp"),6] == soilProp_66

# https://www.sthda.com/english/wiki/two-proportions-z-test-in-r
# i. Build a contingency table
ITS_specialistsSummed <- ITS_specialistsProportion_air_long %>% 
  group_by(HabitatAir, proportionType) %>% 
  summarize(numbASVsinHabSourceType = sum(ReadsByHab, na.rm = TRUE), .groups = "drop")
# View(ITS_specialistsSummed)
ITS_specialistsSummed

# HabitatAir proportionType numbASVsinHabSourceType
# <chr>      <chr>                            <dbl>
#   1 forest     phylloSpecProp                  391927
# 2 forest     soilSpecProp                      1248
# 3 savanna    phylloSpecProp                  339885
# 4 savanna    soilSpecProp                      1321

# Create a contingency table
ITS_specialistsSummed_Table <- as.data.frame(matrix(nrow= 2, ncol=2))
colnames(ITS_specialistsSummed_Table) <- c("foliar", "soil")
rownames(ITS_specialistsSummed_Table) <- c("forest", "savanna")
# Fill in values from the dataframe above
ITS_specialistsSummed_Table[1,1] <- ITS_specialistsSummed[1,3]
ITS_specialistsSummed_Table[1,2] <- ITS_specialistsSummed[2,3]
ITS_specialistsSummed_Table[2,1] <- ITS_specialistsSummed[3,3]
ITS_specialistsSummed_Table[2,2] <- ITS_specialistsSummed[4,3]
ITS_specialistsSummed_Table
ITS_specialistsSummed # compare to make sure these are correct and they are!
ITS_specialistsSummed_Table$totalGroup <- rowSums(ITS_specialistsSummed_Table) #add in a column 
# that specifies total in forest and savanna
ITS_specialistsSummed_Table

# Perform "two-proportions z-test" for foliar surfaces. Prop test below has numbers of foliar in each group,
# then total numbers in each group. Forest is first, so x= goes foliar forest, then foliar savanna, then total forest, then total savanna
ITS_foliarPropTest <- prop.test(x = c(ITS_specialistsSummed_Table[1,1], ITS_specialistsSummed_Table[2,1]), 
                                n = c(ITS_specialistsSummed_Table[1,3], ITS_specialistsSummed_Table[2,3]))
ITS_foliarPropTest # this shows that air samples from matrix have a higher proportion of foliar-associated reads!
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(ITS_specialistsSummed_Table[1, 1], ITS_specialistsSummed_Table[2, 1]) out of c(ITS_specialistsSummed_Table[1, 3], ITS_specialistsSummed_Table[2, 3])
# X-squared = 25.288, df = 1, p-value = 4.938e-07
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   0.0004220239 0.0009727809
# sample estimates:
#   prop 1    prop 2 
# 0.9968258 0.9961284 

# Test for soil
ITS_specialistsSummed_Table
# x goes soil forest, then soil savanna. N goes forest total then savanna total
ITS_soilPropTest <- prop.test(x = c(ITS_specialistsSummed_Table[1,2], ITS_specialistsSummed_Table[2,2]), 
                              n = c(ITS_specialistsSummed_Table[1,3], ITS_specialistsSummed_Table[2,3]))
ITS_soilPropTest #shows that soil indicators are higher in the patch than in the matrix.
# X-squared = 25.288, df = 1, p-value = 4.938e-07
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.0009727809 -0.0004220239
# sample estimates:
#   prop 1      prop 2 
# 0.003174159 0.003871561 
sum(ITS_specialistsSummed_Table[,3]) #734,381 entered as the total amount (i.e. N)

# T-tests or non-parametric (Wilcoxon signed rank) to compare, ignorning habitat type
head(ITS_specialistsProportion_air_long)
str(ITS_specialistsProportion_air_long)
ITS_specialistsProportion_air_long$proportionType <- as.factor(ITS_specialistsProportion_air_long$proportionType)
str(ITS_specialistsProportion_air_long)
# Test to see if I can do a one-way ANOVA
car::leveneTest(ReadsByHab ~ proportionType, data = ITS_specialistsProportion_air_long) # shows that variances are NOT roughly homogenous, since significant
# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value    Pr(>F)    
# group   1  151.35 < 2.2e-16 ***
#       218    

# Since assumptions aren't met, try (paired) wilcoxon signed-rank

# Pivot wider to get two numbers to line up
ITS_specPropAir_wide <- ITS_specialistsProportion_air_long %>% 
  select(-ReadsByHab) %>%  #remove because next step does not "widen" this
  pivot_wider(names_from = proportionType, values_from= proportionOfReads)
ITS_specPropAir_wide
#  Wilcoxon signed-rank test (paired, one-sided: foliar > soil)
# (wilcox.test ignores exact p when there are ties/large n which is ok!)
wilcoxPropReads <- wilcox.test(
  ITS_specPropAir_wide$phylloSpecProp, ITS_specPropAir_wide$soilSpecProp,
  paired = TRUE, #because samples are not independent
  alternative = "greater", #tests greater foliar than soil
  exact = FALSE,    #avoids warnings with ties/large n
  conf.int = TRUE,  #gives a CI for the shift (one-sided CI)
  conf.level = 0.95
)
wilcoxPropReads

# Wilcoxon signed rank test with continuity correction
# data:  ITS_specPropAir_wide$phylloSpecProp and ITS_specPropAir_wide$soilSpecProp
# V = 6105, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0
# 95 percent confidence interval:
#   0.7759215       Inf
# sample estimates:
#   (pseudo)median 
# 0.7816257 

##############
# WHAT ARE THE MEAN PROPORTIONS OF LEAF AND SOIL INDICATOR TAXA IN EACH AIR SAMPLE?
##############
# This is the ASV table for all samples
ASVs_allITSr_noAirSingsDoubs.ps
ITS_phylloSpecialists #names of foliar surface indicator taxa
length(ITS_phylloSpecialists) #2955
ITS_soilSpecialists #names of soil indicator taxa (755 of these)

# GET READS PER SAMPLE
rownames(ASVs_allITSr_noAirSingsDoubs.ps) #rownames are the ASVs
rownames(ASVs_allITSr_noAirSingsDoubs.ps)[length(rownames(ASVs_allITSr_noAirSingsDoubs.ps))] #this is the last ASV
colnames(ASVs_allITSr_noAirSingsDoubs.ps) #column names are the sample names 
# How many reads are in each sample?
readsPerSample <- colSums(ASVs_allITSr_noAirSingsDoubs.ps)
readsPerSample #good, all are about 8500, rarefying threshold
readsPerSample2 <- as.data.frame(matrix(ncol=2, nrow=length(readsPerSample)))
colnames(readsPerSample2) <- c("sampleName", "readsPerSample")
readsPerSample2[,1] <- names(readsPerSample)
readsPerSample2[,2] <- unname(readsPerSample)
readsPerSample2

head(rownames(ASVs_allITSr_noAirSingsDoubs.ps)) #first one is ASV_1
tail(rownames(ASVs_allITSr_noAirSingsDoubs.ps)) #last one is ASV_16858
ASVs_allITSr_noAirSingsDoubs_long <- ASVs_allITSr_noAirSingsDoubs.ps %>% 
  t() %>% #transpose so that ASV names are columns
  as.data.frame() %>% #re-make a dataframe because t() function un-did this
  rownames_to_column(var= "sampleName") %>% #make the now rows into a new column called sampleName
  pivot_longer(cols = ASV_1:ASV_16858, names_to = "ASV_name", values_to = "readCount")
# View(ASVs_allITSr_noAirSingsDoubs_long)
# Add in reads per sample
ASVs_allITSr_noAirSingsDoubs_longer <- left_join(ASVs_allITSr_noAirSingsDoubs_long, readsPerSample2, by= "sampleName")
# Add in indicator type for each ASV
colnames(ASVs_allITSr_noAirSingsDoubs_longer)
ASVs_allITSr_noAirSingsDoubs_longer <- ASVs_allITSr_noAirSingsDoubs_longer %>%
  mutate(
    indictorType = case_when(
      ASV_name %in% ITS_phylloSpecialists ~ "foliarIndicator",
      ASV_name %in% ITS_soilSpecialists ~ "soilIndicator",
      TRUE ~ NA_character_
    )
  )
# View(ASVs_allITSr_noAirSingsDoubs_longer)
dim(ASVs_allITSr_noAirSingsDoubs_longer) # WOAH! 4,047,408 rows!

# Now get proportions of each ASV (within each sample, automatically with calculation below given the structure of the data frame)
ITS_proportionsEachASV <- ASVs_allITSr_noAirSingsDoubs_longer #make a copy
head(ITS_proportionsEachASV)
ITS_proportionsEachASV$propASV <- ITS_proportionsEachASV$readCount/ITS_proportionsEachASV$readsPerSample
# View(ITS_proportionsEachASV)

# Add in sample type 
soilSampNames <- rownames(sample_data(allITSr_noAirSingsDoubs.ps)[which(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType == "soil")]) 
soilSampNames
airSampNames <- rownames(sample_data(allITSr_noAirSingsDoubs.ps)[which(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType == "air")])
airSampNames
phylloSampNames <- rownames(sample_data(allITSr_noAirSingsDoubs.ps)[which(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType == "phyllosphere")])
phylloSampNames

ITS_proportionsEachASV <- ITS_proportionsEachASV %>%
  mutate(
    sampleType = case_when(
      sampleName %in% airSampNames ~ "airSample",
      sampleName %in% soilSampNames ~ "soilSample",
      sampleName %in% phylloSampNames ~ "leafSample",
      TRUE ~ NA_character_
    )
  )

# Split into foliar surface and soil indicators
ITS_propsITS_soilASVs <- ITS_proportionsEachASV[which(ITS_proportionsEachASV$indictorType == "soilIndicator"),]
unique(ITS_propsITS_soilASVs$indictorType)
ITS_propsFoliarASVs <- ITS_proportionsEachASV[which(ITS_proportionsEachASV$indictorType == "foliarIndicator"),]
unique(ITS_propsFoliarASVs$indictorType)

# Get MEAN ASV abundance for each ASV within each sample, and then sort by this:
# FOLIAR SURFACES (proportions only for foliar surface indicator taxa)
colnames(ITS_propsFoliarASVs)
head(ITS_propsFoliarASVs)
meanASVPerSample_foliar <- ITS_propsFoliarASVs %>% 
  group_by(ASV_name, sampleType) %>% 
  summarize(meanASV_relAbundPerType = mean(propASV)) %>% 
  arrange(sampleType, meanASV_relAbundPerType)
# Merge with taxonomic information
allITSr_noAirSingsDoubsTax <- as.data.frame(phyloseq::tax_table(allITSr_noAirSingsDoubs.ps), stringsAsFactors = F)
colnames(allITSr_noAirSingsDoubsTax)
allITSr_noAirSingsDoubsTax <- allITSr_noAirSingsDoubsTax %>% 
  rownames_to_column(var= "ASV_name")
str(allITSr_noAirSingsDoubsTax) #12492 ASVs
meanASVPerSample_foliar <- left_join(meanASVPerSample_foliar, allITSr_noAirSingsDoubsTax, by = "ASV_name")
length(unique(meanASVPerSample_foliar$ASV_name)) #2955 trimmed out the ASVs not found in the air
# View(meanASVPerSample_foliar)

# Look at only air samples
meanASVPerSample_foliarOnlyAir <- meanASVPerSample_foliar %>% 
  filter(sampleType == "airSample") %>% 
  arrange(desc(meanASV_relAbundPerType))
# View(meanASVPerSample_foliarOnlyAir)

# SOIL (proportions only for soil indicator taxa)
colnames(ITS_propsITS_soilASVs)
meanASVPerSample_soil <- ITS_propsITS_soilASVs %>% 
  group_by(ASV_name, sampleType) %>% 
  summarize(meanASV_relAbundPerType = mean(propASV)) %>% 
  arrange(sampleType, meanASV_relAbundPerType)
# Merge with taxonomic information
meanASVPerSample_soil <- left_join(meanASVPerSample_soil, allITSr_noAirSingsDoubsTax, by = "ASV_name")
# View(meanASVPerSample_soil)
meanASVPerSample_soilOnlyAir <- meanASVPerSample_soil %>% 
  filter(sampleType == "airSample")
# View(meanASVPerSample_soilOnlyAir)

## FINALLY, LOOK AT ALL AT ONCE TO GET TOP AIR ASVs, this will go into supplemental materials
head(ITS_proportionsEachASV)
ITS_proportionsEachASV_justAirMean <- ITS_proportionsEachASV %>% 
  filter(sampleType == "airSample") %>% #consider only air samples
  group_by(ASV_name) %>% #grouping ASVs together
  summarize(meanASV_relAbund_air = mean(propASV)) %>% #get mean propprtion of the ASV (but only across bioaerosol samples)
  arrange(desc(meanASV_relAbund_air)) #arrange so that the one with the highest mean is on top
# View(ITS_proportionsEachASV_justAirMean)

# Get just top 100
ITS_proportionsEachASV_justAir_t100 <- ITS_proportionsEachASV_justAirMean %>% 
  slice_head(n = 100)
# View(ITS_proportionsEachASV_justAir_t100)

# Now merge with taxonomy
ITS_proportionsEachASVtax_Air_100 <- left_join(ITS_proportionsEachASV_justAir_t100, allITSr_noAirSingsDoubsTax, by = "ASV_name")
# View(ITS_proportionsEachASVtax_Air_100) 

# Add in indicator type
head(ITS_proportionsEachASV[,c(2,5)]) #head shows columns 2 and 5, which are ASV_name and indicator type
ITS_proportionsEachASVtax_Air_100 <- left_join(ITS_proportionsEachASV_justAir_t100, ITS_proportionsEachASV[,c(2,5)], by = "ASV_name")
ITS_proportionsEachASVtax_Air_100 <- ITS_proportionsEachASVtax_Air_100 %>% 
  distinct() #remove duplicate rows that were caused by merge above (was because ITS_proportionsEachASV was in long format)
# View(ITS_proportionsEachASVtax_Air_100)

# Get number of air samples that these 100 taxa occur in
colnames(ITS_proportionsEachASVtax_Air_100)
top100Air_ASVs <- unique(ITS_proportionsEachASVtax_Air_100$ASV_name)
# 1. Subset phyloseq object to be only air samples
unique(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType)
ITSairOnly.ps <- subset_samples(allITSr_noAirSingsDoubs.ps, sampleType == "air")
head(otu_table(ITSairOnly.ps))
dim(as.data.frame(as.matrix(otu_table(ITSairOnly.ps)))) #ASVs are rows, as they should be 
# each entry below is TRUE if the value is greater than 0, and FALSE otherwise, effectively counting the occurrences!
occurrencesEachASVair <- rowSums(as.data.frame(as.matrix(otu_table(ITSairOnly.ps))) > 0)
# Edit this to make it mergable with the dataframe above (needs to be a dataframe)
occurrencesEachASVair_2 <- as.data.frame(matrix(ncol=2, nrow= length(occurrencesEachASVair)))
colnames(occurrencesEachASVair_2) <- c("ASV_name", "numAirOccurrences")
occurrencesEachASVair_2$ASV_name <- names(occurrencesEachASVair)
occurrencesEachASVair_2$numAirOccurrences <- unname(occurrencesEachASVair)
# Merge with above dataframe
ITS_proportionsEachASVtax_Air_100_v2 <- left_join(ITS_proportionsEachASVtax_Air_100, occurrencesEachASVair_2, by = "ASV_name")
ITS_proportionsEachASVtax_Air_100_v2$percentOccupancy <- (ITS_proportionsEachASVtax_Air_100_v2$numAirOccurrences/110)*100 #110 bioaerosol samples
# Double check with ASV_11, one of the top bioaerosol ASVs
length(which(otu_table(ITSairOnly.ps)[rownames(otu_table(ITSairOnly.ps)) %in% "ASV_11",]> 0)) #110 every sample!
ITS_proportionsEachASVtax_Air_100_v2$numAirOccurrences[which(ITS_proportionsEachASVtax_Air_100_v2$ASV_name == "ASV_11")] #20 matches above!
# Double check with ASV_169, 59th top bioaerosol ASVs
length(which(otu_table(ITSairOnly.ps)[rownames(otu_table(ITSairOnly.ps)) %in% "ASV_169",]> 0)) #101 times
ITS_proportionsEachASVtax_Air_100_v2$numAirOccurrences[which(ITS_proportionsEachASVtax_Air_100_v2$ASV_name == "ASV_169")] #101 matches above!

# DOUBLE CHECK ITS_proportionsEachASV_justAirMean MADE A BIT ABOVE 
ITS_proportionsEachASV_justAirMean
# DOUBLE CHECK-- IT WORKED!
airASVsTopAirCheck_df1 <- otu_table(ITSairOnly.ps)[rownames(otu_table(ITSairOnly.ps)) %in% top100Air_ASVs,]
# Look at ASV 11. Divide how much each sample has of ASV 11 by its total read count. 
mean(airASVsTopAirCheck_df1[rownames(airASVsTopAirCheck_df1 ) %in% "ASV_11",]/colSums(otu_table(ITSairOnly.ps))) == ITS_proportionsEachASV_justAirMean[which(ITS_proportionsEachASV_justAirMean$ASV_name == "ASV_11"),2]

# GET NUMBER OF SOIL AND FOLIAR SURFACE SAMPLES THAT THESE TOP AIR TAXA OCCUR IN
# SOIL-- i.e., mean prop ASV within a given ASV for each soil sample
ITS_soilOnly.ps <-  subset_samples(allITSr_noAirSingsDoubs.ps, sampleType == "soil")
head(otu_table(ITS_soilOnly.ps))
dim(as.data.frame(as.matrix(otu_table(ITS_soilOnly.ps)))) #ASVs are rows, as they should be 
# each entry below is TRUE if the value is greater than 0, and FALSE otherwise, effectively counting the occurrences!
occurrencesEachASV_soil <- rowSums(as.data.frame(as.matrix(otu_table(ITS_soilOnly.ps))) > 0)
# Edit this to make it mergable with the dataframe above (needs to be a dataframe)
occurrencesEachASV_soil_2 <- as.data.frame(matrix(ncol=2, nrow= length(occurrencesEachASV_soil)))
colnames(occurrencesEachASV_soil_2) <- c("ASV_name", "numSoilOccurrences")
occurrencesEachASV_soil_2$ASV_name <- names(occurrencesEachASV_soil)
occurrencesEachASV_soil_2$numSoilOccurrences <- unname(occurrencesEachASV_soil)
occurrencesEachASV_soil_topAir <- occurrencesEachASV_soil_2 %>% 
  filter(ASV_name %in% top100Air_ASVs) #get only those that are top air ASVs
# Double check with ASV_42, one of the top bioaerosol ASVs
length(which(otu_table(ITS_soilOnly.ps)[rownames(otu_table(ITS_soilOnly.ps)) %in% "ASV_42",]> 0)) #17
occurrencesEachASV_soil_topAir[which(occurrencesEachASV_soil_topAir$ASV_name == "ASV_42"),] #17 matches above!
# Double check with ASV_11, one of the top bioaerosol ASVs
length(which(otu_table(ITS_soilOnly.ps)[rownames(otu_table(ITS_soilOnly.ps)) %in% "ASV_11",]> 0)) #0 times
occurrencesEachASV_soil_topAir[which(occurrencesEachASV_soil_topAir$ASV_name == "ASV_11"),] #0 matches above!
# Add percent occupancy by dividing number of samples in by 155 samples x 100
occurrencesEachASV_soil_topAir$SoilPercentOccupancy <- occurrencesEachASV_soil_topAir$numSoilOccurrences/155*100

# GET NUMBER OF FOLIAR SURFACES AND FOLIAR SURFACE SAMPLES THAT THESE TOP AIR TAXA OCCUR IN
# FOLIAR SURFACES-- i.e., mean prop ASV within a given ASV for each foliar sample
ITS_foliarOnly.ps <-  subset_samples(allITSr_noAirSingsDoubs.ps, sampleType == "phyllosphere")
head(otu_table(ITS_foliarOnly.ps ))
dim(as.data.frame(as.matrix(otu_table(ITS_foliarOnly.ps)))) #ASVs are rows, as they should be. 12492    59
# each entry below is TRUE if the value is greater than 0, and FALSE otherwise, effectively counting the occurrences!
occurrencesEachASV_foliar <- rowSums(as.data.frame(as.matrix(otu_table(ITS_foliarOnly.ps))) > 0)
# Edit this to make it mergable with the dataframe above (needs to be a dataframe)
occurrencesEachASV_foliar_2 <- as.data.frame(matrix(ncol=2, nrow= length(occurrencesEachASV_foliar)))
colnames(occurrencesEachASV_foliar_2) <- c("ASV_name", "numFoliarOccurrences")
occurrencesEachASV_foliar_2$ASV_name <- names(occurrencesEachASV_foliar)
occurrencesEachASV_foliar_2$numFoliarOccurrences <- unname(occurrencesEachASV_foliar)
occurrencesEachASV_foliar_topAir <- occurrencesEachASV_foliar_2 %>% 
  filter(ASV_name %in% top100Air_ASVs) #get only those that are top air ASVs
# Double check with ASV_78, one of the top bioaerosol ASVs
length(which(otu_table(ITS_foliarOnly.ps)[rownames(otu_table(ITS_foliarOnly.ps)) %in% "ASV_78",]> 0)) #1
occurrencesEachASV_foliar_topAir[which(occurrencesEachASV_foliar_topAir$ASV_name == "ASV_78"),] #1 matches above!
# Double check with ASV_11, one of the top bioaerosol ASVs
length(which(otu_table(ITS_foliarOnly.ps)[rownames(otu_table(ITS_foliarOnly.ps)) %in% "ASV_11",]> 0)) #25 times
occurrencesEachASV_foliar_topAir[which(occurrencesEachASV_foliar_topAir$ASV_name == "ASV_11"),] #25 matches above!
# Add percent occupancy by dividing number of samples in by 59 foliar surface samples x 100
occurrencesEachASV_foliar_topAir$FoliarPercentOccupancy <- occurrencesEachASV_foliar_topAir$numFoliarOccurrences/59*100

# PROPORTIONS OF TOP BIOAEROSOL ASVS IN SOIL AND FOLIAR SAMPLES
# SOIL-- i.e., mean prop ASV within a given ASV for each soil sample
ITS_proportionsEachASV_justSoil <- ITS_proportionsEachASV %>%
  filter(sampleType == "soilSample") %>% #only soil
  filter(ASV_name %in% top100Air_ASVs) %>%  #only top air ASVs
  group_by(ASV_name) %>% #group by ASV_name
  summarize(meanASV_relAbund_soil = mean(propASV)) %>% #mean prop ASV within a given ASV for each soil sample
  arrange(desc(meanASV_relAbund_soil))
# View(ITS_proportionsEachASV_justSoil)
# DOUBLE CHECK-- IT WORKED!
soilASVsTopAirCheck_df1 <- otu_table(ITS_soilOnly.ps)[rownames(otu_table(ITS_soilOnly.ps)) %in% top100Air_ASVs,]
# Look at ASV 11. Divide how much each sample has of ASV 11 by its total read count. 
mean(soilASVsTopAirCheck_df1[rownames(soilASVsTopAirCheck_df1) %in% "ASV_11",]/colSums(otu_table(ITS_soilOnly.ps))) == ITS_proportionsEachASV_justSoil[which(ITS_proportionsEachASV_justSoil$ASV_name == "ASV_11"),2]

# FOLIAR
ITS_proportionsEachASV_justLeaves <- ITS_proportionsEachASV %>%
  filter(sampleType == "leafSample") %>%
  filter(ASV_name %in% top100Air_ASVs) %>%
  group_by(ASV_name) %>%
  summarize(meanASV_relAbund_foliar = mean(propASV)) %>%
  arrange(desc(meanASV_relAbund_foliar))
# View(ITS_proportionsEachASV_justLeaves)

# MERGE ALL OF THIS TOGETHER AND ADD IN TAXONOMY
ITS_topAirTable1 <- merge(ITS_proportionsEachASVtax_Air_100_v2, ITS_proportionsEachASV_justSoil, by= "ASV_name")
ITS_topAirTable2 <- merge(ITS_topAirTable1, occurrencesEachASV_soil_topAir,by= "ASV_name")
ITS_topAirTable3 <- merge(ITS_topAirTable2, ITS_proportionsEachASV_justLeaves, by= "ASV_name")
ITS_topAirTable4 <- merge(ITS_topAirTable3, occurrencesEachASV_foliar_topAir, by= "ASV_name")
ITS_topAirTable_final <- left_join(ITS_topAirTable4, allITSr_noAirSingsDoubsTax, by = "ASV_name")
# View(ITS_topAirTable_final)

# Clean up ITS_topAirTable_final
head(ITS_topAirTable_final)
ITS_topAirTable_final <- ITS_topAirTable_final %>% 
  select(-Kingdom) %>% 
  arrange(desc(meanASV_relAbund_air))
# change NA for non-indicators to be "not indicator"
ITS_topAirTable_final$indictorType[which(is.na(ITS_topAirTable_final$indictorType))] <- "not indicator"
# Make all of these proportions percentages
ITS_topAirTable_final$meanASV_relAbund_air <- ITS_topAirTable_final$meanASV_relAbund_air*100
ITS_topAirTable_final$meanASV_relAbund_soil <- ITS_topAirTable_final$meanASV_relAbund_soil*100
ITS_topAirTable_final$meanASV_relAbund_foliar <- ITS_topAirTable_final$meanASV_relAbund_foliar*100
# View(ITS_topAirTable_final)
# Round off all of these to the nearest tenth
ITS_topAirTable_final$meanASV_relAbund_air <- round(ITS_topAirTable_final$meanASV_relAbund_air, 2) 
ITS_topAirTable_final$meanASV_relAbund_soil <- round(ITS_topAirTable_final$meanASV_relAbund_soil, 2) 
ITS_topAirTable_final$meanASV_relAbund_foliar <- round(ITS_topAirTable_final$meanASV_relAbund_foliar, 2) 
ITS_topAirTable_final$percentOccupancy <- round(ITS_topAirTable_final$percentOccupancy, 1)
ITS_topAirTable_final$SoilPercentOccupancy <- round(ITS_topAirTable_final$SoilPercentOccupancy, 1)
ITS_topAirTable_final$FoliarPercentOccupancy <- round(ITS_topAirTable_final$FoliarPercentOccupancy, 1)

# Make sure that indicator makes sense and it does! 
colnames(ITS_topAirTable_final)
# View(ITS_topAirTable_final[,c("ASV_name", "indictorType", "meanASV_relAbund_soil" , "SoilPercentOccupancy",
#                        "meanASV_relAbund_foliar","FoliarPercentOccupancy" )])

# Make a streamlined version for the manuscript supplemental x
colnames(ITS_topAirTable_final)
ITS_topAirTable_forSupp <- ITS_topAirTable_final[,c(1:3,5,8,11:13,15:17)]
ITS_topAirTable_forSupp$indictorType[which(ITS_topAirTable_forSupp$indictorType == "soilIndicator")] <- "soil"
ITS_topAirTable_forSupp$indictorType[which(ITS_topAirTable_forSupp$indictorType == "foliarIndicator")] <- "foliar"
ITS_topAirTable_forSupp$indictorType[which(ITS_topAirTable_forSupp$indictorType == "not indicator")] <- "not"
ITS_topAirTable_forSupp$ASV_name <- gsub(x=ITS_topAirTable_forSupp$ASV_name, pattern= "ASV_", replacement = "")
ITS_topAirTable_forSupp$ASV_name <- as.character(ITS_topAirTable_forSupp$ASV_name)
# View(ITS_topAirTable_forSupp)

# write.csv(ITS_topAirTable_forSupp, file = "Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_topAirTable_forSupp_09-28-2025.csv") #written September 28, 2025

##############
#  A FEW EXTRA CALCULATIONS FOR MANUSCRIPT
##############
### Double-checking this in manuscript:
# A fungal example is the genus Sporobolomyces, which is a dominant colonizer of foliar surfaces of many species of plants 
# (Ruinen 1963, Vacher et al. 2016) and which was detected in 61% of all bioaerosol samples at about 0.2% average relative abundance within each sample 
head(ITS_topAirTable_forSupp)
ITS_topAirTable_forSupp[which(ITS_topAirTable_forSupp$Genus == "Sporobolomyces"),4] # Shows ASVs detected in 60.9, so 61% is correct!
ITS_topAirTable_forSupp[which(ITS_topAirTable_forSupp$Genus == "Sporobolomyces"),2] # Shows ASVs mean abundance 0.18, so 0.2% is correct!

# To fill in these sentences:
# Specifically, the fungal reads in bioaerosol samples consisted of a median of ___% foliar surface- and _____% soil-associated taxa...
head(ASVs_allITSr_noAirSingsDoubs_longer) #made earlier
ITS_allSampsIndicReads <- ASVs_allITSr_noAirSingsDoubs_longer #make a quick copy to edit
# Make non-indicators labeled just to make sure calculations don't get messed up
ITS_allSampsIndicReads$indictorType[which(is.na(ITS_allSampsIndicReads$indictorType))] <- "notIndicator"
which(is.na(ITS_allSampsIndicReads$readCount)) #none are NA, but remove it below in code any
PropsPerAirSampPerIndicType <- ITS_allSampsIndicReads %>% 
  filter(str_detect(sampleName, "air_ITS")) %>%  #get only air samples
  group_by(sampleName, indictorType) %>% 
  summarize(readsInType = sum(readCount), .groups = "drop_last") %>%
  group_by(sampleName) %>%
  mutate(PropSampIndType = readsInType/ sum(readsInType, na.rm = TRUE)) %>%
  ungroup()
# View(PropsPerAirSampPerIndicType)
# DOUBLE CHECK: Test these calculations for a air_ITS_106 (has good amount of reads in all categories)
# 1. What is total number of reads?
air_ITS_106_ASVtab <- otu_table(ITSairOnly.ps)[,colnames(otu_table(ITSairOnly.ps)) %in% "air_ITS_106"] #get ASV table subset
PropsPerAirSampPerIndicType_106 <-  PropsPerAirSampPerIndicType[which(PropsPerAirSampPerIndicType$sampleName == "air_ITS_106"),] #get PropsPerAirSampPerIndicType subset
colSums(air_ITS_106_ASVtab) == sum(PropsPerAirSampPerIndicType_106$readsInType) # TRUE and is #8500 !
# 2. Check proportions for each type (1 and 3 are foliar and soil specialists)
colSums(air_ITS_106_ASVtab[rownames(air_ITS_106_ASVtab) %in% ITS_phylloSpecialists,])/colSums(air_ITS_106_ASVtab) == PropsPerAirSampPerIndicType_106$PropSampIndType[1]
colSums(air_ITS_106_ASVtab[rownames(air_ITS_106_ASVtab) %in% ITS_soilSpecialists,])/colSums(air_ITS_106_ASVtab) == PropsPerAirSampPerIndicType_106$PropSampIndType[3]

# Medians and means of foliar and soil in air samples
PropsPerAirSampPerIndicType %>% 
  group_by(indictorType) %>% 
  reframe(meanPropIndicType = mean(PropSampIndType))
# A tibble: 3 × 2
# indictorType    meanPropIndicType
# <chr>                       <dbl>
#   1 foliarIndicator           0.784  
# 2 notIndicator              0.214  
# 3 soilIndicator             0.00275
PropsPerAirSampPerIndicType %>% 
  group_by(indictorType) %>% 
  reframe(medianPropIndicType = median(PropSampIndType))
# indictorType    medianPropIndicType
# <chr>                         <dbl>
# 1 foliarIndicator             0.785  
# 2 notIndicator                0.213  
# 3 soilIndicator               0.00200
# 78.5% foliar surface- and 0.2% soil-associated taxa added to manuscript

# ##############
# # SAME ANALYSIS, BUT WITH PRESENCE/ABSENCE 0F INDICATOR TAXA (I.E. FOR EACH SAMPLE, OUT OF ALL OF THE 
# # ASVS THAT OCCURRED, WHAT PROPRPTION WERE SOIL OR PHYLLOSPHERE INDICATORS (e.g. if a sample has 
# # a richness of 100 ASVs and 10 were fungi (abundance doesn't matter), then would have value of .1))
# ##############
# 
ITS_specialistsPresAbs <- as.data.frame(matrix(data=NA, nrow=ncol(ASVs_allITSr_noAirSingsDoubs.ps), ncol = 7))
colnames(ITS_specialistsPresAbs) <- c("sampleName", "sampleType", "numberPhylloIndic", "numberSoilIndic", "totalASVsGreaterThanZero", "propNumPhyllo", "propNumSoil")
# Add in sample names 
ITS_specialistsPresAbs[,1] <- rownames(ITSrarefied_metadat) #ITSrarefied_metadat was created from ASVs_allITSr_noAirSingsDoubs.ps earlier in script 
# Add in sample types (i.e., air, soil, or phyllosphere)
ITS_specialistsPresAbs[,2] <- ITSrarefied_metadat$sampleType

# Before starting, confirm once more that these indices are correct (they are!)
unique(rownames(ASVs_allITSr_noAirSingsDoubs.ps)[ITS_phylloSpecialistsBigIndex] == ITS_phylloSpecialists)
unique(rownames(ASVs_allITSr_noAirSingsDoubs.ps)[ITS_soilSpecialistsBigIndex] == ITS_soilSpecialists)

# Fill in the PRESENCE/ABSENCE for each
# This ignores abundance, getting for each sample the count of foliar surface or soil indicator taxa that appear in each sample.
# The proportion is then this number divided by the total number of ASVs (as presence/absence) in each sample
for (k in 1:ncol(ASVs_allITSr_noAirSingsDoubs.ps)){ #k indexes each sample (columns) and there are 324 samples total
  # Next 2 code lines go only on the rows where there are phyllo or soil specialists, respectively, and give each sample a one if it has at least one count of that ASV.
  # In other words, it sums up total presence/absence of phyllosphere and soil specialists in each sample
  ITS_specialistsPresAbs[k,3] <- sum(ASVs_allITSr_noAirSingsDoubs.ps[ITS_phylloSpecialistsBigIndex,k] > 0)  #add up all of the times that there is a count more than 0
  ITS_specialistsPresAbs[k,4] <- sum(ASVs_allITSr_noAirSingsDoubs.ps[ITS_soilSpecialistsBigIndex,k] > 0)  #add up all of the times that there is a count more than 0
  ITS_specialistsPresAbs[k,5] <- sum(ASVs_allITSr_noAirSingsDoubs.ps[,k] > 0) #adds up all of the total number of non-zero ASVs in each sample k
  ITS_specialistsPresAbs[k,6] <- ITS_specialistsPresAbs[k,3]/ITS_specialistsPresAbs[k,5] #divides number of present foliar surface ASVs by number of total present ASVs in sample
  ITS_specialistsPresAbs[k,7] <- ITS_specialistsPresAbs[k,4]/ITS_specialistsPresAbs[k,5] #divides number of present soil ASVs by number of total present ASVs
}

# What is average richness in each? Wow, air is very rich!
ITS_specialistsPresAbs %>% 
  group_by(sampleType) %>% 
  summarize(medianNumASVs = median(totalASVsGreaterThanZero))
# sampleType   medianNumASVs
# <chr>                <dbl>
# 1 air                    724
# 2 phyllosphere           836
# 3 soil                   287

ITS_specialistsPresAbs %>% 
  group_by(sampleType) %>% 
  summarize(sdNumASVs = sd(totalASVsGreaterThanZero)) 
# sampleType   sdNumASVs
# <chr>            <dbl>
# 1 air             273. 
# 2 phyllosphere    185. 
# 3 soil            72.5

# Add in a column that specifies if that air sample is forest or savanna
ITS_specialistsPresAbs_airOnly <- ITS_specialistsPresAbs %>%
  filter(sampleType == "air") %>% #pull out only air samples
  mutate(HabitatAir = ifelse(sampleName %in% ITS_airForestSampNames, "forest", "savanna")) #calls it forest if in ITS_airForestSampNames, otherwise, "savanna"
# View(ITS_specialistsPresAbs_airOnly)

# PLOT IT!
# FINALLY, PLOT (AIR ONLY COMPARISONS)
# SOIL COMPARISON
ITS_airOnlyPresAbs_soil_plot <- ggplot(ITS_specialistsPresAbs_airOnly, aes(x=HabitatAir, y=propNumSoil, fill=HabitatAir)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=c("darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=1, alpha=0.9, height = 0) +
  scale_x_discrete(name= "Underlying vegetation type", labels=c("forest", "savanna")) +
  scale_y_continuous(name= "Number of soil indicator taxa") + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  ggtitle("Soil fungal indicator taxa across air samples")  +
  theme(plot.title = element_text(size=20)) + #increase plot title size
  theme(legend.position = "none")  #remove legend
# quartz()
ITS_airOnlyPresAbs_soil_plot
# 
# PHYLLOSPHERE COMPARISON
ITS_airOnlyPresAbs_phyllo_plot <- ggplot(ITS_specialistsPresAbs_airOnly, aes(x=HabitatAir, y=propNumPhyllo, fill=HabitatAir)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=c("darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=1, alpha=0.9, height = 0) + #height = zero ensures no vertical jitter
  scale_x_discrete(name= "Underlying vegetation type", labels=c("forest", "savanna")) +
  scale_y_continuous(name= "Number of phyllosphere indicator taxa") + # limITS=c(25, 225) add this argument to change limITS
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  ggtitle("Foliar surface fungal indicator taxa across air samples")  +
  theme(plot.title = element_text(size=20)) + #increase plot title size
  theme(legend.position = "none")  #remove legend

#quartz()
ITS_airOnlyPresAbs_phyllo_plot

### PRESENCE/ABSENCE: PLOT WITH AIR FOREST AND AIR SAVANNA AND FOLIAR SURFACES AND PHYLLOSPHERE AND SOIL ####
# Savanna air names (object made above)
ITS_airSavSampNames
# Forest air names (object made above)
ITS_airForestSampNames
# Get soil names 
ITS_soilSampNames
# Get foliar surface names 
ITS_foliarSampNames 
# Add in a column that specifies what kind of sample it is
colnames(ITS_specialistsProportion)
ITS_specialistsPresAbs_allTypes <- ITS_specialistsPresAbs %>% 
  mutate(
    sampleTypeWithAir = case_when(
      sampleName %in% ITS_airSavSampNames ~ "SavAirSample",
      sampleName %in% ITS_airForestSampNames ~ "ForAirSample",
      sampleName %in% ITS_foliarSampNames ~ "foliarSample",
      sampleName %in% ITS_soilSampNames ~ "soilSample",
      TRUE ~ NA_character_
    )
  )
colnames(ITS_specialistsPresAbs_allTypes)
# View(ITS_specialistsPresAbs_allTypes)

# Make dataframe "longer" for faceting purposes
colnames(ITS_specialistsPresAbs_allTypes)
ITS_specialistsPresAbs_allTypes_long <- ITS_specialistsPresAbs_allTypes %>%
  pivot_longer(cols = c(propNumPhyllo, propNumSoil), 
               names_to = "proportionType", 
               values_to = "proportionOfTaxaByPresAbs")

# Re-order x-axis elements
unique(ITS_specialistsPresAbs_allTypes_long$sampleTypeWithAir)
ITS_specialistsPresAbs_allTypes_long$sampleTypeWithAir <- 
  factor(ITS_specialistsPresAbs_allTypes_long$sampleTypeWithAir, 
         levels= c("ForAirSample", "SavAirSample", "foliarSample", "soilSample"))

# Vector for new facet labels
facet_labels_PresAbs <- c("propNumPhyllo" = "Foliar surfaces",
                          "propNumSoil" = "Soil")

# Create faceted plot
colnames(ITS_specialistsPresAbs_allTypes_long)
ITS_allTypesPresAbs_2panels <- ggplot(ITS_specialistsPresAbs_allTypes_long, aes(x=sampleTypeWithAir, y=proportionOfTaxaByPresAbs, fill=sampleTypeWithAir)) + 
  geom_boxplot() + 
  geom_jitter(aes(color = sampleTypeWithAir), size=1, alpha=0.9, height = 0) + #height = 0 insures no horizontal jitter
  facet_wrap(~proportionType, scales = "fixed", labeller = labeller(proportionType = facet_labels_PresAbs)) +
  theme_bw() +
  scale_fill_manual(values=c("cornflowerblue", "cornflowerblue", "chartreuse4", "chocolate4")) +
  scale_color_manual(values = c("black", "black", "black", "black")) + #color for jittered points
  scale_x_discrete(labels=c("bioaerosol\nmatrix", "bioaerosol\npatch", "foliar surface", "soil")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        ) +
  #ggtitle("Bacterial indicator taxa across all samples")  +
  labs(y = "Proportion of leaf or soil indicator taxa",
       x = NULL) +
  theme(plot.title = element_text(size=16)) + #increase plot title size
  theme(legend.position = "none") + #remove legend 
  ylim(0, 1.00) + #make y limits between 0 and 1
  theme(strip.text = element_text(size = 16))
ITS_allTypesPresAbs_2panels

# SAVE PLOT (saved October 25, 2025) 
# saveRDS(ITS_allTypesPresAbs_2panels, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_allTypesPresAbs_2panels_10-25-25.rds")

### PRESENCE/ABSENCE: STATISTCS
# Paired Wilcoxon signed-rank test
head(ITS_specialistsPresAbs_allTypes)
dim(ITS_specialistsPresAbs_allTypes)
# Get dataframe 
ITS_specialistsPresAbs_air <- ITS_specialistsPresAbs_allTypes %>% 
  filter(sampleType == "air")
head(ITS_specialistsPresAbs_air)

#  Wilcoxon signed-rank test (paired, one-sided: foliar > soil)
# (wilcox.test ignores exact p when there are ties/large n which is ok!)
ITS_wilcoxPresAbs <- wilcox.test(
  ITS_specialistsPresAbs_air$propNumPhyllo, ITS_specialistsPresAbs_air$propNumSoil,
  paired = TRUE, #samples are not independent
  alternative = "greater", #foliar > soil is test
  exact = FALSE,    #avoids warnings with ties/large n
  conf.int = TRUE,  #gives a CI for the shift (one-sided CI)
  conf.level = 0.95
)
ITS_wilcoxPresAbs
# data:  ITS_specialistsPresAbs_air$propNumPhyllo and ITS_specialistsPresAbs_air$propNumSoil
# V = 6105, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0
# 95 percent confidence interval:
#   0.4562191       Inf
# sample estimates:
#   (pseudo)median 
# 0.4664759 

# ##############
# MAKE A VENN DIAGRAM
# ##############
# Get phyllosphere only dataset
ITS_phyllo_rarefied_noControls <- subset_samples(allITSr_noAirSingsDoubs.ps, sampleType == "phyllosphere")
# Get soil only dataset
ITS_soil_rarefied_noControls <- subset_samples(allITSr_noAirSingsDoubs.ps, sampleType == "soil")

# Get air ASVs (i.e., any that were detected in bioaerosols)
ITS_airASVs <- ASVs_outta_ps(ITSairOnly.ps)
unique(colSums(ITS_airASVs))
ITS_airASVnames <- names(which(rowSums(ITS_airASVs) > 0)) #get all of the ASVs that have at least one sequence
unique(sort(rowSums(ITS_airASVs[rownames(ITS_airASVs) %in% ITS_airASVnames,]))) #none are zero and no singletons and doubletons
length(ITS_airASVnames) #4673

# Get foliar surface ASVs (i.e., any that were detected in foliar surface samples)
ITS_phylloASVs <- ASVs_outta_ps(ITS_phyllo_rarefied_noControls)
ITS_phylloASVnames <- names(which(rowSums(ITS_phylloASVs) > 0))
unique(sort(rowSums(ITS_phylloASVs[rownames(ITS_phylloASVs) %in% ITS_phylloASVnames,]))) #none are zero
length(ITS_phylloASVnames) #4703

# Get soil ASVs (i.e., any that were detected in soil samples)
ITS_soilASVs <- ASVs_outta_ps(ITS_soil_rarefied_noControls)
ITS_soilASVnames <- names(which(rowSums(ITS_soilASVs) > 0))
unique(sort(rowSums(ITS_soilASVs[rownames(ITS_soilASVs) %in% ITS_soilASVnames,]))) #none are zero
length(ITS_soilASVnames) == sum(4541, 169, 369, 533) #5612 and this equals what's reported in Venn diagram

# Make a list with an element for each of the categories of ASVs
ITS_ASVsForVenn <- list(
  ITS_air <- ITS_airASVnames,
  ITS_soil <- ITS_soilASVnames,
  ITS_foliarSurface <- ITS_phylloASVnames
)
str(ITS_ASVsForVenn)

# Using the vignette as a guide, manually recreate and edit Venn Diagram
vignette("fully-customed", package = "ggVennDiagram")

ITS_venn <- Venn(ITS_ASVsForVenn)
ITS_vennData <- process_data(ITS_venn)
ITS_vennPlot <- ggplot() +
  # 1. region count layer
  geom_polygon(aes(X, Y, fill = id, group = id), 
               data = venn_regionedge(ITS_vennData)) +
  # 2. set edge layer
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(ITS_vennData),
            linewidth = 3.5,
            show.legend = FALSE) +
  # 3. set label layer -- remove labels 
  # geom_text(aes(X, Y, label = name), 
  #       data = venn_setlabel(ITS_vennData)) +
  # 4. region label layer-- changing this from geom_label to geom_text removes label outlines
  geom_text(aes(X, Y, label = count), 
            data = venn_regionlabel(ITS_vennData),
            size=5.5) +
  # 5. Change colors so that all of the filled in areas are white
  scale_color_manual(values = c("cornflowerblue","chocolate4","chartreuse4")) +
  scale_fill_manual(values = c("white", "white", "white","white", "white", "white", "white")) +
  coord_equal() +
  theme_void() + 
  theme(legend.position='none') #remove legend

ITS_vennPlot #excellent, matches original Venn diagram
# saveRDS(ITS_vennPlot, file = "RobjectsSaved/ITS_vennPlot_Jan8_2024.rds") #saved January 14, 2024
# saveRDS(ITS_vennPlot, file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_vennPlot_09-29-2025.rds") #saved Oct 25, 2025

# Get an all-white version without text for editing if need be:
vennPlotNOTEXT <- ggplot() +
  # 1. region count layer
  geom_polygon(aes(X, Y, fill = id, group = id), 
               data = venn_regionedge(ITS_vennData)) +
  # 2. set edge layer
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(ITS_vennData),
            linewidth = 3.5,
            show.legend = FALSE) +
  # 3. set label layer -- remove labels 
  # geom_text(aes(X, Y, label = name), 
  #       data = venn_setlabel(ITS_vennData)) +
  # 4. region label layer-- changing this from geom_label to geom_text removes label outlines
  geom_text(aes(X, Y, label = count), 
            data = venn_regionlabel(ITS_vennData),
            size=5.5, color= "white") +
  # 5. Change colors so that all of the filled in areas are white
  scale_color_manual(values = c("cornflowerblue","chocolate4","chartreuse4")) +
  scale_fill_manual(values = c("white", "white", "white","white", "white", "white", "white")) +
  coord_equal() +
  theme_void() + 
  theme(legend.position='none') #remove legend
vennPlotNOTEXT
# saveRDS(vennPlotNOTEXT, file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/vennPlotNOTEXT_09-29-2025.rds") #saved October 25, 2025

# Double check that the calculations that were manually made in the Venn Diagrams are correct and get ASV names for subsets:
# i. Only soil
# Union gets all of the ASVs that are found in the air and phyllo datasets and sees how this differs
# from those that are found only in the soil
ITS_unique_soil <- setdiff(ITS_soilASVnames, union(ITS_airASVnames, ITS_phylloASVnames))
length(ITS_unique_soil) == 4541 #great, this shows that calculations were correct!
# ii. Only foliar surface
ITS_unique_foliar <- setdiff(ITS_phylloASVnames, union(ITS_airASVnames, ITS_soilASVnames))
length(ITS_unique_foliar) == 2745 #great, this shows that calculations were correct!
# iii. Only in bioaerosols
ITS_unique_air <- setdiff(ITS_airASVnames, union(ITS_soilASVnames, ITS_phylloASVnames))
length(ITS_unique_air) #3079, correct, matches Venn calculations
# iv. Overlap between air and soil
length(which(ITS_airASVnames %in% ITS_soilASVnames==TRUE)) == 169 + 369 #overlap between air and soil
airSoilOverlapNames <- ITS_airASVnames[which(ITS_airASVnames %in% ITS_soilASVnames==TRUE)]
# v. Overlap between air and foliar surfaces
length(which(ITS_airASVnames %in% ITS_phylloASVnames==TRUE)) == 1056 + 369 #overlap between air and foliar surfaces
airFoliarOverlapNames <- ITS_airASVnames[which(ITS_airASVnames %in% ITS_phylloASVnames==TRUE)]
# vi. Shared among all!
allOverlapNames <- intersect(airFoliarOverlapNames, ITS_soilASVnames)
length(allOverlapNames) == 369
# vii. Overlap between air and soil NO FOLIAR
airSoilOverLapOnlyNames <- setdiff(airSoilOverlapNames, allOverlapNames) #takes the 169+369 shared by air and soil and subtracts the 369
length(airSoilOverLapOnlyNames) == 169
# viii. Overlap between air and foliar (NO SOIL)
airFoliarOverLapOnlyNames <- setdiff(airFoliarOverlapNames, allOverlapNames) #takes the 1056 + 369 shared by air and phyllo and subtracts the 369 shared by all
length(airFoliarOverLapOnlyNames) == 1056

#####################
# SOME ADDITIONAL CALCULATIONS FOR MANUSCRIPT AND VENN DIAGRAMS (TO ADD IN POWERPOINT OR INKSCAPE)
#####################
# In manuscript: what percentage of total foliar surface ASVs overlapped with soil for fungi (compared with for bacteria, which
# was calculated in a different script)?
((369+533)/length(ITS_phylloASVnames))*100 #19.17925. Added 19.2% to manuscript

# Set up of air reads
ITS_airASVnames #these are the ASVs present in air 
3079+169+369+1056 == length(ITS_airASVnames) #shows that it matches Venn Diagram
ITS_airASVs #Made above and gets air ASV table
head(ITS_airASVs)
# For manuscript: Fill in text here: Second, there were some taxa present in the bioaerosols that were not found in either the sampled soils or foliar surfaces, but the contribution of these taxa was minimal
# (xxx% of fungal and xxx% of all bacterial ASVs detected in the air samples
length(ITS_unique_air)/length(ITS_airASVnames)*100 #65.88915, so 65.9%

# To add within blue circles: What percentage do all of these ASVs (3079, 169, 369, and 1056)
# make up of the total air community? (Thinking of all air samples merged as one air sample)
# BUILD A DATAFRAME TO GET PERCENTAGE OF EACH SUBSET

head(ITS_airASVs) #this is the ASV table with ASVs as rows and samples as columns
colnames(ITS_airASVs)
# i. Get a longer dataframe with ASV name, sample name, and ASV abundance in that sample
ITS_airASVs_df <- ITS_airASVs %>% 
  rownames_to_column(var="ASVname") %>% 
  pivot_longer(air_ITS_1:air_ITS_99, names_to = "sampleName", values_to = "ASVabundance")
head(ITS_airASVs_df)
# ii. Grouping by each ASV, get the total abundance of that ASV (across only air samples)
unique(ITS_airASVs_df$sampleName)
ITSASVabund_AirSamps <- ITS_airASVs_df %>% 
  group_by(ASVname) %>% 
  summarize(ASVtotalAbund = sum(ASVabundance))
sum(ITSASVabund_AirSamps$ASVtotalAbund) #933866
# iii. Grouping by each ASV, proportion that each ASV makes up of total community is that ASV total 
# abundance divided by the total ASV abundance of all reads in the air samples
ITSASVabund_AirSamps_2 <- ITSASVabund_AirSamps %>% 
  group_by(ASVname) %>% 
  mutate(propASVabund = ASVtotalAbund/sum(ITSASVabund_AirSamps$ASVtotalAbund))
# iv. For the ASVs that are only in air (3,079 ASVs), what percentage is this of total reads?
# This code pulls out the proportional abundances of all of the ASVs unique to the air and then 
# adds these together
sumOf_ITSuniqueAir <- sum(ITSASVabund_AirSamps_2$propASVabund[ITSASVabund_AirSamps_2$ASVname %in% ITS_unique_air])
sumOf_ITSuniqueAir #0.1561027, so 15.6% goes below 3,079, i.e., bioaerosol-only taxa
# v. For the ASVs that are shared between soil and air (169 ASVs), what percentage is this of total reads?
sumOf_ITSairSoil <- sum(ITSASVabund_AirSamps_2$propASVabund[ITSASVabund_AirSamps_2$ASVname %in% airSoilOverLapOnlyNames])
sumOf_ITSairSoil #0.00805683, so 0.8% goes below 169
# vi. For the ASVs that are shared between phyllo and air (1056 ASVs), what percentage is this of total reads?
sumOf_ITSairPhyllo <- sum(ITSASVabund_AirSamps_2$propASVabund[ITSASVabund_AirSamps_2$ASVname %in% airFoliarOverLapOnlyNames])
sumOf_ITSairPhyllo  #0.5250304%, so 52.5% goes below 1056, i.e., bioaerosol and foliar surface taxa
# vii. Of the ASVs that are shared among all taxa, what proportion are these among air samples?
sumOf_ITSall <- sum(ITSASVabund_AirSamps_2$propASVabund[ITSASVabund_AirSamps_2$ASVname %in% allOverlapNames])
sumOf_ITSall #0.3108101, so 31.1% goes below 369
# DO THESE ADD UP TO 1, AS THEY SHOULD?
sum(sumOf_ITSuniqueAir, sumOf_ITSairSoil, sumOf_ITSairPhyllo, sumOf_ITSall) #YES!

