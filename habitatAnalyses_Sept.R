# habitatAnalyses_Sept.R
# TAXONOMY BARPLOTS AND QPCR RESULTS FOR AIR SAMPLES WITHIN EACH SAMPLING DAY

# Script description: The goal of this script is to investigate the question: Does the near-surface atmosphere
# differ in the forest and the savanna? This script follows full16S_EDArarefied_Part2.R (and fungal equivalent),
# because from these scripts airOnly_16Sr_noAirSingsDoubs_phyloseq.rds and airOnly_ITS_noAirSingsDoubs_phyloseq.rds
# are used. 

# ALL below repeated for 16S and ITS datasets
# 1. Exploration of compositional differences:
## a. Plots:
### i. Stacked barplots showing taxonomic differences b/w open patch and forest matrix air
### ii. NMDSes:
##### 1) showing differences among air, soil, and phyllosphere samples
##### 2) between air microbiome in the open patch and in the forested matrix
### iii. Taxonomic barplot by day
### iv. Heatmap showing top taxa (all groups)-- note that some of this code was repeated from EDA_rarefied_Part2.R 
## b. Analyses:
# PERMANOVAs that test BCdists among air samples ~ EU (random effect) + Habitat (with permutations restricted) 

# 2. qPCR 
## a. Plots :
### i. boxplots of qPCR results by sampling day comparing forest and savanna
### ii. with all days combined comparing forest and savanna
## b. Analyses: generalized linear models with negative binomial error distribution

# Many of these plots went in my PowerPoint presented to BROADN May 2, 2024.

#######################################
# I.  SCRIPT SET UP
#######################################
##### 1. LOAD LIBRARIES #####
library(tidyverse)
library(ggplot2); packageVersion("ggplot2")
library(indicspecies)
library(vegan)
library(phyloseq)
library(stringr)
library(Polychrome) #for color palette
library(grid) #for making a legend
library(gridExtra)
library(ggh4x)
library(emmeans)
library(car)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(bbmle)
library(performance) #for check_dispersion function
library(svglite)
library(ggnewscale)

options(scipen = 999) #I don't like scientific notation here!

##### 2. SET WORKING DIRECTORY AND READ IN DATA #####
# Read in phyloseq objects FOR JUST AIR DATA 
# on server: 
## I6S ###
airOnly_16Sr_noAirSingsDoubs_phyloseq <- readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/air16S_noSingDoubs_ps_sept25.rds")
## ITS ###
airOnly_ITS_noAirSingsDoubs_phyloseq <- readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/airOnly_ITS_noAirSingsDoubs_ps.rds")

# Read in earlier versions of all phyloseq objects for soil, phyllosphere, and control information
## I6S ###
all16Sr_noAirSingsDoubs.ps <- readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/all16Sr_noAirSingsDoubs_ps_sept25.rds")
all16Sr_noAirSingsDoubs.ps
load("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_dcAP_rarefied_ps_Sept15") 
I6S_dcAP_rarefied.ps 

## ITS ###
allITSr_noAirSingsDoubs.ps <- readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/sallITSr_noAirSingsDoubs_ps.rds") #made in fullITS_EDArarefied_Part2_Sept.R
allITSr_noAirSingsDoubs.ps
load(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/allITS_rarefied_phyloseq") #Made in fullITS_EDA_part1.R
allITS_rarefied.ps

# Get unique EUs
uniqueEUs <- unique(sample_data(airOnly_16Sr_noAirSingsDoubs_phyloseq)$EU)
uniqueEUs
unique(sample_data(airOnly_16Sr_noAirSingsDoubs_phyloseq)$daysOut)
# double check, once more -- looks correct!
cbind(sample_data(airOnly_16Sr_noAirSingsDoubs_phyloseq)$daysOut, sample_data(airOnly_16Sr_noAirSingsDoubs_phyloseq)$DateSetOut)
colnames(sample_data(airOnly_16Sr_noAirSingsDoubs_phyloseq))

# This lists the daysOut by EU, with the four dates for EU 52, EU 53S, EU54S, 
# and EU 8, respectively. 
daysOutByEU <- c("16-Jun-2022", "18-Jun-2022", "26-Jun-2022", "1-Jul-2022", "22-Jun-2022", "27-Jun-2022", "2-Jul-2022",
                 "6-Jul-2022", "24-Jun-2022", "30-Jun-2022", "5-Jul-2022", "8-Jul-2022", "20-Jun-2022", "23-Jun-2022",
                 "29-Jun-2022", "4-Jul-2022")

###############################################
# II.  EXPLORATION OF COMPOSITIONAL DIFFERENCES 
###############################################
############# A. PLOTS ##################
###### i. Stacked barplots showing taxonomic differences b/w open patch and forest matrix air ######
# These plots have only a few taxa and were made for the ISME presentation
############# BACTERIAL CLASSES #############
# ***ALMOST THE SAME COLOR PALETTE AS ANCOM_traits_I6SandITS.R EXCEPT ANCOM HAS MORE CLASESS AND CLOSTRIDIA NOT INCLUDED AND COLOR RECYCLED
# From https://github.com/Nowosad/rcartocolor (except for "#FE8F42", which I think will work)
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#FE8F42")
scales::show_col(safe_colorblind_palette)


colsForClassesI6S <- c("Bacilli" = "#88CCEE","Alphaproteobacteria" = "#CC6677","Actinobacteria"= "#DDCC77","Bacteroidia" = "#117733","Deinococci" = "#332288",
                             "Clostridia" ="#AA4499", "Gammaproteobacteria" = "#44AA99", "Acidobacteriae" =  "#999933","Verrucomicrobiae" = "#882255","Planctomycetes" =  "#661100")



airOnly_16Sr_class.glom <-  phyloseq::tax_glom(airOnly_16Sr_noAirSingsDoubs_phyloseq, taxrank = "Class") 
length(unique(tax_table(airOnly_16Sr_class.glom)[,3])) #67

# Transform sample counts on just glommed samples 
airOnly_16Sr_noAirSingsDoubs_phyloseq.Class.0 <- transform_sample_counts(airOnly_16Sr_class.glom, function(x) x / sum(x) )
#ASVs in this above are just a representative from each class
length(unique(tax_table(airOnly_16Sr_noAirSingsDoubs_phyloseq.Class.0)[,3])) #still 67!

# Merge samples so that we only have combined abundances for open patch and matrix
airOnly_16Sr_relabun.Class.1 <- merge_samples(airOnly_16Sr_noAirSingsDoubs_phyloseq.Class.0 , group = "HabitatAir") 
sample_data(airOnly_16Sr_relabun.Class.1) #sample categories are correct

# Convert to proportions again b/c total abundance of each site will equal number of species that were merged
airOnly_16Sr_relabun.Class.2 <- transform_sample_counts(airOnly_16Sr_relabun.Class.1, function(x) x / sum(x))
sample_data(airOnly_16Sr_relabun.Class.2)

# Now get only taxa that comprise a certain proportion of the abundance (the rest will be grouped together as one color)
airOnly_16Sr_relabun.Class.df <-psmelt(airOnly_16Sr_relabun.Class.2)
colnames(airOnly_16Sr_relabun.Class.df) #metadata variables 
airOnly_16Sr_relabun.Class.top99 <- airOnly_16Sr_relabun.Class.df 
airOnly_16Sr_relabun.Class.top98 <- airOnly_16Sr_relabun.Class.df 
airOnly_16Sr_relabun.Class.top97 <- airOnly_16Sr_relabun.Class.df 

airOnly_16Sr_relabun.Class.top99$Class[airOnly_16Sr_relabun.Class.top99$Abundance < 0.01] <- "< 1% abundance"
airOnly_16Sr_relabun.Class.top98$Class[airOnly_16Sr_relabun.Class.top98$Abundance < 0.02] <- "< 2% abundance"
airOnly_16Sr_relabun.Class.top97$Class[airOnly_16Sr_relabun.Class.top97$Abundance < 0.03] <- "< 3% abundance"

airOnly_16Sr_top_99p_Class <- unique(airOnly_16Sr_relabun.Class.top99$Class)
length(airOnly_16Sr_top_99p_Class) #16
airOnly_16Sr_top_98p_Class <- unique(airOnly_16Sr_relabun.Class.top98$Class) #try this one first?
length(airOnly_16Sr_top_98p_Class) #12
airOnly_16Sr_top_97p_Class <- unique(airOnly_16Sr_relabun.Class.top97$Class)
length(airOnly_16Sr_top_97p_Class) #6

# Tidy up the dataframe a little bit before plotting
airOnly_16Sr_top_98p_Class
airOnly_16Sr_relabun.Class.top98$Class <- gsub(pattern= "NA", x=airOnly_16Sr_relabun.Class.top98$Class, replacement = "Unclassified")
unique(airOnly_16Sr_relabun.Class.top98$Class) #fixed

# Get good, sequential colors for different classes
setdiff(names(colsForClassesI6S), unique(airOnly_16Sr_relabun.Class.top98$Class)) #Verrucomicrobiae"
setdiff(unique(airOnly_16Sr_relabun.Class.top98$Class), names(colsForClassesI6S)) #"Myxococcia"   "Unclassified" "< 2% abundance"
# These above show that need to add colors for #"Myxococcia"   "Unclassified" "< 2% abundance" adn remove for #Verrucomicrobiae"
colsForClassesI6S_3 <- colsForClassesI6S[-which(names(colsForClassesI6S) %in% "Verrucomicrobiae" == TRUE)] #remove Verrucomicrobiae
# Add #"Myxococcia"   "Unclassified" "< 2% abundance"
toAdd <- c("Myxococcia" =  "#6699CC",
           "Unclassified" = "thistle4",
           "< 2% abundance" = "lightgray")

colsForClassesI6S_3 <- c(colsForClassesI6S_3, toAdd) #make sure these are compatible with ANCOM_traits_I6SandITS.R

# PLOT!
airBactClassAir.plot <- ggplot(data=airOnly_16Sr_relabun.Class.top98, aes(x=Sample, y=Abundance, fill=Class)) + theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold"))
# quartz()
airBactClassAir.plot <- airBactClassAir.plot + geom_bar(aes(), stat="identity", position="fill") + theme_bw() +
  scale_fill_manual(values = colsForClassesI6S_3) + 
  theme(legend.position="bottom",
        legend.text = element_text(size = 14),
        axis.text.x=element_text(size=18),
        axis.text.y =element_text(size=18),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18)
  ) +
  guides(fill=guide_legend(nrow=4)) +
  labs( y = "Relative abundance"
  )
airBactClassAir.plot #exported as airBactClassAirPlot.pdf

# Get these proportions of each class in each habitat type
colnames(airOnly_16Sr_relabun.Class.df)
str(airOnly_16Sr_relabun.Class.df)
airOnly_16Sr_relabun.Class.df$Abundance #relative abundance of each class in both air habitat types
airOnly_16Sr_top_Class <- airOnly_16Sr_relabun.Class.df %>%
  group_by(Sample, Class) %>%
  summarise(Mean = mean(Abundance), na.rm = TRUE) %>%
  arrange(desc(Mean))
# View(airOnly_16Sr_top_Class)

# GET MEAN PERCENTAGES OF CLASS IN EACH SAMPLE FOR MANUSCRIPT
# This is where this is reported in MS: "The bacterial assemblages were dominated by classes Alphaproteobacteria (27%), Bacilli (15%), and Actinobacteria (14%), as shown in Fig. 2b and Table S4."
# Produced earlier in the script, "gloms" at the class level and produces a new phyloseq 
airOnly_16Sr_class.glom 
unique(tax_table(airOnly_16Sr_class.glom)[,3]) #has 67 unique classes (some unclassified)
# Produced earlier in the script and transformed sample counts on just glommed samples 
airOnly_16Sr_noAirSingsDoubs_phyloseq.Class.0 
# View(as.data.frame(as.matrix(otu_table(airOnly_16Sr_noAirSingsDoubs_phyloseq.Class.0))))
# Check, yes, all samples sum to 1 as they should
colSums(as.data.frame(as.matrix(otu_table(airOnly_16Sr_noAirSingsDoubs_phyloseq.Class.0))))
# "Melt" together different parts of the phyloseq 
airSamp_16Sr_relabun_class.df <-psmelt(airOnly_16Sr_noAirSingsDoubs_phyloseq.Class.0)
#View(airSamp_16Sr_relabun_class.df)
# Get the mean relative abundance in each sample!
colnames(airSamp_16Sr_relabun_class.df)
# View(classPerSample)
meanClassPerSample <- airSamp_16Sr_relabun_class.df %>% 
  group_by(Class) %>% 
  summarize(meanClassAbundance = mean(Abundance))
# View(meanClassPerSample)

########## STACKED BAR PLOT FOR SAMPLES X HABITAT (FOR RENDEZVOUS 2025 POSTER)-- looks nicer
airOnly_16Sr_class.glom #this has has all the sample numbers for each class
# View(as.data.frame(otu_table(airOnly_16Sr_class.glom)))
head(as.data.frame(otu_table(airOnly_16Sr_class.glom)))
airOnly_16Sr_class.glom_t <- as.data.frame(t(as.data.frame(otu_table(airOnly_16Sr_class.glom))))
rownames(airOnly_16Sr_class.glom_t)
airOnly_16Sr_class.glom_t <- rownames_to_column(airOnly_16Sr_class.glom_t, var = "sampleName")
head(airOnly_16Sr_class.glom_t)
#View(airOnly_16Sr_class.glom_t)
colnames(airOnly_16Sr_class.glom_t)
airOnly_16Sr_class.glomLong <- pivot_longer(airOnly_16Sr_class.glom_t, cols= colnames(airOnly_16Sr_class.glom_t)[2]:colnames(airOnly_16Sr_class.glom_t)[ncol(airOnly_16Sr_class.glom_t)], names_to = "ASV_name", values_to = "count")
# View(airOnly_16Sr_class.glomLong)
# ADD BACK IN CLASS INFO
airOnly_16Sr_class.glomTAXA <- rownames_to_column(as.data.frame(tax_table(airOnly_16Sr_class.glom)), var= "ASV_name")
head(airOnly_16Sr_class.glomTAXA)
head(airOnly_16Sr_class.glomLong)
# Merge the two dataframes, keeping only Kingdom, Phylum, and Class
air16S_classDFlong <- merge(airOnly_16Sr_class.glomLong, airOnly_16Sr_class.glomTAXA[,1:4], by="ASV_name")
head(air16S_classDFlong) #looks good!
unique(air16S_classDFlong$ASV_name) #73 unique, as before
unique(air16S_classDFlong$Class)
# Add in habitat data
air16S_classMetaDat <- as.data.frame(sample_data(airOnly_16Sr_class.glom))
rownames(air16S_classMetaDat)
air16S_classMetaDat_2 <- rownames_to_column(air16S_classMetaDat, var= "sampleName")
# View(air16S_classMetaDat_2)
colnames(air16S_classMetaDat_2)
air16S_classMetaDat_2 <- as.data.frame(as.matrix(air16S_classMetaDat_2[,c(1,41,42,43)])) #trim down some of the columns
air16S_classDFlongBySamp <- merge(air16S_classDFlong, air16S_classMetaDat_2, by = "sampleName")

# Calculate proportions within each sample
air16S_classDF_BySamp <- air16S_classDFlongBySamp %>%
  group_by(HabitatAir, sampleName) %>%
  mutate(Proportion = count / sum(count)) %>%
  ungroup()

# Combine Habitat and Sample for spacing
air16S_classDF_BySamp$HabitatSample <- interaction(air16S_classDF_BySamp$HabitatAir, air16S_classDF_BySamp$sampleName, sep = "_")
#View(air16S_classDF_BySamp)
colnames(air16S_classDF_BySamp)

# Show only those classes that were in the original barplots.
unique(airOnly_16Sr_relabun.Class.top98$Class) #these are the classes that I want to keep
unique(air16S_classDF_BySamp$Class) #here are all the classes still present
# Change any NAs to unclassified
air16S_classDF_BySamp$Class <- gsub(pattern= "NA", x=air16S_classDF_BySamp$Class, replacement = "Unclassified")
sort(unique(air16S_classDF_BySamp$Class)) #great, now we have no NAs!
# Get those classes that are less than 2% of reads in either all savanna or all forest air samples
classesLess2percent <- unique(air16S_classDF_BySamp$Class)[which(unique(air16S_classDF_BySamp$Class) %in% unique(airOnly_16Sr_relabun.Class.top98$Class)== FALSE)]
# Name these <2% abundance
air16S_classDF_BySamp$Class[which(air16S_classDF_BySamp$Class %in% classesLess2percent== TRUE)] <- "< 2% abundance"
sort(unique(air16S_classDF_BySamp$Class))

# They all add up to one, as expected
doubleCheckSampProps <- air16S_classDF_BySamp %>% 
  group_by(sampleName) %>% 
  summarize(sampTotal = sum(Proportion))
unique(air16S_classDF_BySamp$Class)

airBacterialClassesBySamplePlot <- ggplot(air16S_classDF_BySamp, aes(x = HabitatSample, y = Proportion, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = function(x) gsub(".*_", "", x)) + # Just show Sample
  scale_fill_manual(values = colsForClassesI6S_3) + 
  facet_grid(. ~ HabitatAir, scales = "free_x", space = "free_x") +
  ylab("Relative Abundance") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, colour = "black"),  # x-axis tick size
    axis.text.y = element_text(size = 14, colour = "black"),                         # y-axis tick size
    axis.title.x = element_text(size = 12, colour = "black"),                        # x-axis label
    axis.title.y = element_text(size = 16, colour = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5, size =16)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1)) + #this instead of setting limits prevents cutting off samples slightly > 1 due to floating point
  guides(fill = guide_legend(ncol = 3)) +
  xlab(NULL) +
  labs(fill = "Bacterial Class")

airBacterialClassesBySamplePlot #(saved as Figures/BactClassesBySample_Sept18.pdf on computer)

# Save it -- again Sept 18
# ggsave("~/Desktop/CU_Research/SRS_Aeromicrobiome/Figures/airBacterialClassesBySamplePlot_Sept2025.svg", plot = airBacterialClassesBySamplePlot, width = 7.5, height = 4.5, units = "in")
# Save another version without the legend
airBacterialClassesBySamplePlotNoLeg <- airBacterialClassesBySamplePlot + theme(legend.position = "none")

# SAVE VERSION WITH RIGHT LEGEND
airBacterialClassesBySamplePlotRightLeg <- ggplot(air16S_classDF_BySamp, aes(x = HabitatSample, y = Proportion, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = function(x) gsub(".*_", "", x)) + # Just show Sample
  scale_fill_manual(values = colsForClassesI6S_3) + 
  facet_grid(. ~ HabitatAir, scales = "free_x", space = "free_x") +
  ylab("Relative Abundance") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, colour = "black"),  # x-axis tick size
    axis.text.y = element_text(size = 14, colour = "black"),                         # y-axis tick size
    axis.title.x = element_text(size = 12, colour = "black"),                        # x-axis label
    axis.title.y = element_text(size = 16, colour = "black"),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5, size =16)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1)) + 
  guides(fill = guide_legend(ncol = 1)) +
  xlab(NULL) +
  labs(fill = "Bacterial Class")

airBacterialClassesBySamplePlotRightLeg

############# FUNGI #############
##### ORDER #####
airOnly_ITS_order.glom <-  phyloseq::tax_glom(airOnly_ITS_noAirSingsDoubs_phyloseq, taxrank = "Order") 
length(unique(tax_table(airOnly_ITS_order.glom)[,4])) #78

# Transform sample counts on just glommed samples 
airOnly_ITS_noAirSingsDoubs_phyloseq.Order.0 <- transform_sample_counts(airOnly_ITS_order.glom, function(x) x / sum(x) )
#ASVs in this above are just a representative from each order
length(unique(tax_table(airOnly_ITS_noAirSingsDoubs_phyloseq.Order.0)[,4])) #still 78!

# Merge samples so that we only have combined abundances for open patch and matrix
airOnly_ITS_relabun.Order.1 <- merge_samples(airOnly_ITS_noAirSingsDoubs_phyloseq.Order.0 , group = "HabitatAir") 
sample_data(airOnly_ITS_relabun.Order.1) #sample categories are correct

# Convert to proportions again b/c total abundance of each site will equal number of species that were merged
airOnly_ITS_relabun.Order.2 <- transform_sample_counts(airOnly_ITS_relabun.Order.1, function(x) x / sum(x))
sample_data(airOnly_ITS_relabun.Order.2)

# Now get only taxa that comprise a certain proportion of the abundance (the rest will be grouped together as one color)
airOnly_ITS_relabun.Order.df <-psmelt(airOnly_ITS_relabun.Order.2)
colnames(airOnly_ITS_relabun.Order.df) #metadata variables 
airOnly_ITS_relabun.Order.top99 <- airOnly_ITS_relabun.Order.df 

# 99%
airOnly_ITS_relabun.Order.top99$Order[airOnly_ITS_relabun.Order.top99$Abundance < 0.01] <- "< 1% abundance"

airOnly_ITS_top_99p_Order <- unique(airOnly_ITS_relabun.Order.top99$Order)
length(airOnly_ITS_top_99p_Order) #12 #go with this one!

# Tidy up the dataframe a little bit before plotting
airOnly_ITS_top_99p_Order
airOnly_ITS_relabun.Order.top99$Order <- gsub(pattern= "NA", x=airOnly_ITS_relabun.Order.top99$Order, replacement = "Unclassified")
unique(airOnly_ITS_relabun.Order.top99$Order) #fixed

# Get good, sequential colors for different orders
# NOTE THAT THESE DO NOT MATCH ANCOM_traits_ITSandITS_Sept.R
# Make sure colors match plots in ANCOM_traits_ITSandITS.R (here specifically color_mapping_ITSANCOMplot_2)!!
scales::show_col(safe_colorblind_palette)
# [1] "< 1% abundance"  "Agaricales"      "Auriculariales"  "Cantharellales"  "Capnodiales"     "Corticiales"    
# [7] "Hymenochaetales" "Pleosporales"    "Polyporales"     "Russulales"      "Trechisporales"  "Unclassified"   

color_mapping_ITSOrders <- c(
  "< 1% abundance" = "lightgray", #same as in bacterial plot
  "Agaricales" = "#88CCEE",
  "Auriculariales" = "#CC6677",
  "Cantharellales" = "#DDCC77",
  "Capnodiales" =  "#117733",
  "Corticiales" = "#332288",
  "Hymenochaetales" = "#AA4499",
  "Pleosporales" = "#999933",
  "Polyporales" = "#44AA99",
  "Russulales" = "#882255",
  "Trechisporales" = "#661100",
  "Unclassified" = "thistle4"  #same as in bacterial plot
)

# PLOT!
airFungOrderAir.plot <- ggplot(data=airOnly_ITS_relabun.Order.top99, aes(x=Sample, y=Abundance, fill=Order)) + theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold"))
# quartz()
airFungOrderAir.plot <- airFungOrderAir.plot + geom_bar(aes(), stat="identity", position="fill") + theme_bw() +
  scale_fill_manual(values = color_mapping_ITSOrders) + 
  theme(legend.position="bottom",
        legend.text = element_text(size = 14),
        axis.text.x=element_text(size=18),
        axis.text.y =element_text(size=18),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18)
  ) +
  guides(fill=guide_legend(nrow=4)) +
  labs( y = "Relative abundance"
  )
airFungOrderAir.plot  

# Get these proportions of each order in each habitat type
colnames(airOnly_ITS_relabun.Order.df)
str(airOnly_ITS_relabun.Order.df)
airOnly_ITS_relabun.Order.df$Abundance
airOnly_ITS_top_Order <- airOnly_ITS_relabun.Order.df %>%
  group_by(Sample, Order) %>% #Sample here since "glomming" means that sample is air habitat
  summarise(Mean = mean(Abundance), na.rm = TRUE) %>%
  arrange(desc(Mean))
airOnly_ITS_top_Order

# GET MEAN PERCENTAGES OF ORDER IN EACH SAMPLE FOR MANUSCRIPT
# Here is sentence in MS: "The top fungal orders and their mean relative abundances in bioaerosol samples included Polyporales (49%), Russulales (14%), and Hymenochaetales"
# Produced earlier in the script, "gloms" at the order level and produces a new phyloseq 
airOnly_ITS_order.glom #has 90, 
# Produced earlier in the script and transformed sample counts on just glommed samples 
airOnly_ITS_noAirSingsDoubs_phyloseq.Order.0 
# View(as.data.frame(as.matrix(otu_table(airOnly_ITS_noAirSingsDoubs_phyloseq.Order.0 ))))
# Check, yes, all samples sum to 1 as they should
colSums(as.data.frame(as.matrix(otu_table(airOnly_ITS_noAirSingsDoubs_phyloseq.Order.0 ))))
# "Melt" together different parts of the phyloseq 
airSamp_ITSr_relabun_order.df <-psmelt(airOnly_ITS_noAirSingsDoubs_phyloseq.Order.0)
#View(airSamp_ITSr_relabun_order.df)
# Get the mean relative abundance in each sample!
colnames(airSamp_ITSr_relabun_order.df)
# View(orderPerSample)
meanOrderPerSample <- airSamp_ITSr_relabun_order.df %>% 
  group_by(Order) %>% 
  summarize(meanOrderAbundance = mean(Abundance))
# View(meanOrderPerSample)
meanOrderPerSample

########## STACKED BAR PLOT FOR SAMPLES X HABITAT (FOR RENDEZVOUS 2025 POSTER)
airOnly_ITS_order.glom #this has has all the sample numbers for each order
# View(as.data.frame(otu_table(airOnly_ITS_order.glom)))
head(as.data.frame(otu_table(airOnly_ITS_order.glom)))
# Switch to make samples rows and ASVs columns 
airOnly_ITS_order.glom_t <- as.data.frame(t(as.data.frame(otu_table(airOnly_ITS_order.glom))))
rownames(airOnly_ITS_order.glom_t) #samples
airOnly_ITS_order.glom_t <- rownames_to_column(airOnly_ITS_order.glom_t, var = "sampleName")
head(airOnly_ITS_order.glom_t)
# View(airOnly_ITS_order.glom_t)
colnames(airOnly_ITS_order.glom_t) #column one is ASVs
airOnly_ITS_order.glomLong <- pivot_longer(airOnly_ITS_order.glom_t, cols= colnames(airOnly_ITS_order.glom_t)[2]:colnames(airOnly_ITS_order.glom_t)[ncol(airOnly_ITS_order.glom_t)], names_to = "ASV_name", values_to = "count")
# View(airOnly_ITS_order.glomLong)
# ADD BACK IN ORDER INFO
# Get tax table formatted
airOnly_ITS_order.glomTAXA <- rownames_to_column(as.data.frame(tax_table(airOnly_ITS_order.glom)), var= "ASV_name")
head(airOnly_ITS_order.glomTAXA) #looks good
head(airOnly_ITS_order.glomLong)
# Merge the two dataframes, keeping only Kingdom, Phylum, order, and Order
airITS_orderDFlong <- merge(airOnly_ITS_order.glomLong, airOnly_ITS_order.glomTAXA[,1:5], by="ASV_name")
head(airITS_orderDFlong) #looks good!
# View(airITS_orderDFlong)
sort(unique(airITS_orderDFlong$ASV_name)) #90 unique
sort(unique(airITS_orderDFlong$Order)) #78 unique
# Add in habitat data
airITS_orderMetaDat <- as.data.frame(as.matrix(sample_data(airOnly_ITS_order.glom)))
head(airITS_orderMetaDat)
rownames(airITS_orderMetaDat)
airITS_orderMetaDat_2 <- rownames_to_column(airITS_orderMetaDat, var= "sampleName")
# View(airITS_orderMetaDat_2)
colnames(airITS_orderMetaDat_2)
airITS_orderMetaDat_2 <- as.data.frame(as.matrix(airITS_orderMetaDat_2[,c(1,2,41,42,43)])) #trim down some of the columns
airITS_orderDFlongBySamp <- merge(airITS_orderDFlong, airITS_orderMetaDat_2, by = "sampleName")
head(airITS_orderDFlongBySamp)
unique(airITS_orderDFlongBySamp$HabitatAir)

# Calculate proportions within each sample
colnames(airITS_orderDFlongBySamp)
airITS_orderDF_BySamp <- airITS_orderDFlongBySamp %>%
  group_by(HabitatAir, sampleName) %>%
  mutate(Proportion = count / sum(count)) %>%
  ungroup()

# Combine Habitat and Sample for spacing
airITS_orderDF_BySamp$HabitatSample <- interaction(airITS_orderDF_BySamp$HabitatAir, airITS_orderDF_BySamp$sampleName, sep = "_")
# View(airITS_orderDF_BySamp)
colnames(airITS_orderDF_BySamp)

# Show only those orderes that were in the original barplots.
unique(airOnly_ITS_relabun.Order.top99$Order) #these are the orders that I want to keep
unique(airITS_orderDF_BySamp$Order) #here are all the orders still present
# Change any NAs to unclassified
airITS_orderDF_BySamp$Order <- gsub(pattern= "NA", x=airITS_orderDF_BySamp$Order, replacement = "Unclassified")
sort(unique(airITS_orderDF_BySamp$Order)) #great, now we have no NAs!
# Get those orders that are less than 1% of reads in either all savanna or all forest air samples
ordersLess1percent <- unique(airITS_orderDF_BySamp$Order)[which(unique(airITS_orderDF_BySamp$Order) %in% unique(airOnly_ITS_relabun.Order.top99$Order)== FALSE)]
# Name these <1% abundance
airITS_orderDF_BySamp$Order[which(airITS_orderDF_BySamp$Order %in% ordersLess1percent== TRUE)] <- "< 1% abundance"
sort(unique(airITS_orderDF_BySamp$Order))

airFungalOrdersBySamplePlot <- ggplot(airITS_orderDF_BySamp, aes(x = HabitatSample, y = Proportion, fill = Order)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = function(x) gsub(".*_", "", x)) + # Just show sample
  scale_fill_manual(values = color_mapping_ITSOrders) + 
  facet_grid(. ~ HabitatAir, scales = "free_x", space = "free_x") +
  ylab("Relative Abundance") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, colour = "black"),  # x-axis tick size
    axis.text.y = element_text(size = 14, colour = "black"),                         # y-axis tick size
    axis.title.x = element_text(size = 12, colour = "black"),                        # x-axis label
    axis.title.y = element_text(size = 16, colour = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5, size =16)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1)) + 
  guides(fill = guide_legend(ncol = 3)) +
  xlab(NULL) +
  labs(fill = "Fungal Order")

airFungalOrdersBySamplePlot

###### ii. ORDINATIONS #####
############# BACTERIA #############
# NOTE: Although these are useful, they are NOT on Hellinger-transformed data
#### 1. Ordination air based on habitat type
sample_data(airOnly_16Sr_noAirSingsDoubs_phyloseq)$HabitatAir

# 16S
set.seed(19)
air_16S.ord <- ordinate(airOnly_16Sr_noAirSingsDoubs_phyloseq, "NMDS", "bray") #stress is 0.3218862

I6S_habitatAirOrd <- plot_ordination(airOnly_16Sr_noAirSingsDoubs_phyloseq, air_16S.ord, type="samples", color="HabitatAir") +
  scale_color_manual(values =c("savanna" = "goldenrod", "forest"  = "forestgreen")) +
  geom_point(size=2) +
  theme_bw() +
  theme(panel.grid = element_blank()) # remove ALL gridlines
I6S_habitatAirOrd
# saveRDS(I6S_habitatAirOrd, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_habitatAirOrd_10-17-25.rds")

# ITS
set.seed(190)
air_ITS.ord <- ordinate(airOnly_ITS_noAirSingsDoubs_phyloseq, "NMDS", "bray") #stress is  0.1009676

ITS_habitatAirOrd <- plot_ordination(airOnly_ITS_noAirSingsDoubs_phyloseq, air_ITS.ord, type="samples", color="HabitatAir") +
  scale_color_manual(values =c("savanna" = "goldenrod", "forest"  = "forestgreen")) +
  geom_point(size=2) +
  theme_bw() +
  theme(panel.grid = element_blank()) # remove ALL gridlines
ITS_habitatAirOrd
# saveRDS(ITS_habitatAirOrd, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_habitatAirOrd_10-17-25.rds")

#### 2. "Big" ordination for I6S across air, soil, and phyllosphere
# Grayed out since in earlier EDArarefied_part2_Sept.R script and that script had where these were saved for MS.
# RESULTS WERE THE SAME THOUGH
# unique(as.data.frame((as.matrix(sample_data(all16Sr_noAirSingsDoubs.ps))))$sampleType)
# all16S.ord <- ordinate(all16Sr_noAirSingsDoubs.ps, "NMDS", "bray") #stress is 0.1148462 
# 
# I6S_habitatOrd <- plot_ordination(all16Sr_noAirSingsDoubs.ps, all16S.ord, type="samples", color="sampleType") +
#   scale_color_manual(values =c("soil" = "chocolate", "air" = "cornflowerblue", "phyllosphere" = "forestgreen"))
# 
# sampDat16S <- as.data.frame(as.matrix(sample_data(all16Sr_noAirSingsDoubs.ps)))
# length(unique(sampDat16S$PlantSpecies)) #30
# length(which(sampDat16S$sampleType == "soil")) #157
# length(which(sampDat16S$sampleType == "phyllosphere")) #58
# length(which(sampDat16S$sampleType == "air")) #84

# quartz()
# I6S_habitatOrd + geom_point(size=5) + theme_bw()

#### 3. PCoA just comparing air in open patch and forested matrix
# (NMDS had a stress of 0.32 so it's not good!!!)
sample_data(airOnly_16Sr_noAirSingsDoubs_phyloseq)$HabitatAir
# Get ordination data
onlyAir16S.ord <- ordinate(airOnly_16Sr_noAirSingsDoubs_phyloseq, "PCoA", "bray") 
propVar <- onlyAir16S.ord$values$Relative_eig
cumulativeVar <- cumsum(propVar)
cumulativeVar #first two axes only explain 7% of the variation, so this is also a terrible way
# to show the data

############# FUNGI #############
unique(as.data.frame((as.matrix(sample_data(allITSr_noAirSingsDoubs.ps))))$sampleType)
# Grayed out since in earlier EDArarefied_part2_Sept.R script and that script had where these were saved for MS
# allITS.ord <- ordinate(allITSr_noAirSingsDoubs.ps, "NMDS", "bray") # stress is 0.09624066
# 
# ITS_habitatOrd <- plot_ordination(allITSr_noAirSingsDoubs.ps, allITS.ord , type="samples", color="sampleType") 
# 
# # quartz()
# ITS_habitatOrd + geom_point(size=5) + theme_bw() + 
#   scale_color_manual(values =c("soil" = "chocolate", "air" = "cornflowerblue", "phyllosphere" = "forestgreen"))
# 
# sampDatITS <- as.data.frame(as.matrix(sample_data(allITSr_noAirSingsDoubs.ps)))
# length(unique(sampDatITS$PlantSpecies)) #31
# length(which(sampDatITS$sampleType == "soil")) #155
# length(which(sampDatITS$sampleType == "phyllosphere")) #59
# length(which(sampDatITS$sampleType == "air")) #110

###### iii TAXONOMIC BARPLOTS BY DAY ######
# SEPT. 17, 2025: RETAINED IN SCRIPT IN CASE REVIEWERS WANT IT BUT NOT RE-CHECKED AGAIN SINCE NOT USED IN MANUSCRIPT
# ############# BACTERIA #############
# First split phyloseq object by sampling day and remove ASVs that don't occur in subset
# This lists the daysOut by EU, with the four dates for EU 52, EU 53S, EU54S,
# and EU 8, respectively.
daysOutByEU # (MADE EARLIER ON IN SCRIPT)
I6S_bySamplingDay_psList <- vector("list", length=length(daysOutByEU)) # pre-allocate list to store results, with 1 spot per EU
for (i in 1:length(I6S_bySamplingDay_psList)){
  I6S_bySamplingDay_psList[[i]] <- subset_samples(airOnly_16Sr_noAirSingsDoubs_phyloseq, DateSetOut == daysOutByEU[i])
  I6S_bySamplingDay_psList[[i]] <- prune_taxa(taxa_sums(I6S_bySamplingDay_psList[[i]]) > 0, I6S_bySamplingDay_psList[[i]]) #remove non-occurring ASVs
  names(I6S_bySamplingDay_psList)[[i]] <- unique(sample_data(I6S_bySamplingDay_psList[[i]])$DateSetOut) #this should match the daysOutByEU
}
# Double check!
names(I6S_bySamplingDay_psList) == daysOutByEU

# MAKE PHYLUM LEVEL TAXONOMY PLOTS IN A FOR LOOP
I6S_phylumPlotInfo_list <- vector("list", length=length(I6S_bySamplingDay_psList)) # pre-allocate list to store results, with 1 spot per day
names(I6S_phylumPlotInfo_list) <- names(I6S_bySamplingDay_psList)
I6SuniquePhyla_byDay <- vector("list", length=length(I6S_bySamplingDay_psList)) # pre-allocate list to store unique phyla, with 1 spot per EU
names(I6SuniquePhyla_byDay) <- names(I6S_bySamplingDay_psList)

for (j in 1:length(I6S_phylumPlotInfo_list)){
  # Get all of the standardized relative abundances of each phylum in each EU
  I6S_phylumPlotInfo_list[[j]] <- I6S_bySamplingDay_psList[[j]] %>%
    tax_glom(taxrank = "Phylum") %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    merge_samples(group = "HabitatAir") %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    psmelt()
  # Make the phyla that comprise less than 5% be "< 5% abundance"
  I6S_phylumPlotInfo_list[[j]]$Phylum[I6S_phylumPlotInfo_list[[j]]$Abundance < 0.05] <- "< 5% abundance"
  # Make phyla a factor for better control later of aesthetics in plot
  I6S_phylumPlotInfo_list[[j]]$Phylum <- as.factor(I6S_phylumPlotInfo_list[[j]]$Phylum)
  # Store all of these unique phyla for ease of later plotting (they may differ by EU)
  I6SuniquePhyla_byDay[[j]] <- unique(I6S_phylumPlotInfo_list[[j]]$Phylum)
}

# Get pretty colors
I6SuniquePhyla_byDay #most that a given day has is 9 phyla
I6SuniquePhyla_byDayAll <- sort(unique(unlist(I6SuniquePhyla_byDay)))
I6SuniquePhyla_byDayAll #12 unique phyla, so I need 12 unique colors + <5% abundance
colorsI6SPhyla <- glasbey.colors(15) #do less to avoid black-looking colors
swatch(colorsI6SPhyla)
swatch(colorsI6SPhyla[c(2:4,6:15)]) #looks good
colorsI6SPhyla <- colorsI6SPhyla[c(2:4,6:15)]
swatch(colorsI6SPhyla)

# Master color mapping for all phyla
color_mapping_I6Sphyla <- c(
  "< 5% abundance" = "#766C95",
  "Actinobacteriota" = "#0000FF",
  "Firmicutes" = "#783FC1",
  "Chloroflexi"= "#00FF00",
  "Myxococcota" = "#FF00B6",
  "Cyanobacteria" = "#005300",
  "Proteobacteria" = "#FFD300",
  "Verrucomicrobiota" = "#FF0000",
  "Acidobacteriota" = "#009FFF",
  "Bacteroidota" = "#9A4D42",
  "Deinococcota" = "#00FFBE",
  "Planctomycetota" = "#B1CC71",
  "Bdellovibrionota" = "#1F9698"
)

# Shows that these coloring mapping match!
setdiff(names(color_mapping_I6Sphyla),I6SuniquePhyla_byDayAll)

# Filter the master color mapping to include only the phyla found on each date
I6S_phylumPlotColors <- vector("list", length=length(I6S_bySamplingDay_psList)) # pre-allocate list to store results, with 1 spot per day
names(I6S_phylumPlotColors) <- names(I6S_bySamplingDay_psList)
for (j in 1:length(I6S_phylumPlotColors)){
  I6S_phylumPlotColors[[j]] <- color_mapping_I6Sphyla[names(color_mapping_I6Sphyla) %in% I6SuniquePhyla_byDay[[j]]]
}

# PLOT ALL OF THESE (can't really easily be done in a for loop since number of phyla differs by day).
# Store them in a list
I6S_phylum_PlotList <- vector("list", length=length(I6S_bySamplingDay_psList))
names(I6S_phylum_PlotList) <- names(I6S_phylumPlotInfo_list)

# Can do a little in a for loop:
for (k in 1:length(I6S_phylum_PlotList)) {
  I6S_phylum_PlotList[[k]] <- ggplot(data=I6S_phylumPlotInfo_list[[k]], aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold")) +
    geom_bar(aes(), stat="identity", position="fill") +
    theme_bw() +
    theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5)) +
    theme(legend.text = element_text(colour="black", size = 16))  +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values = I6S_phylumPlotColors[[k]])
}
names(I6S_phylumPlotInfo_list) == daysOutByEU
I6S_phylum_PlotList[[1]]
# Another for loop for plots without legends (plotted with legends above to check everything)
I6S_phylum_PlotListNoLeg <- I6S_phylum_PlotList
for (k in 1:length(I6S_phylum_PlotList)) {
  I6S_phylum_PlotListNoLeg[[k]] <- ggplot(data=I6S_phylumPlotInfo_list[[k]], aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold")) +
    geom_bar(aes(), stat="identity", position="fill") +
    theme_bw() + scale_fill_manual(values = I6S_phylumPlotColors[[k]]) +
    theme(legend.position="none") +
    #remove all x axis labels
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    #remove all y axis labels
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )

}

# quartz()
do.call(grid.arrange, c(I6S_phylum_PlotListNoLeg, ncol = 4, nrow = 4))

### Make a new legend with ALL colors, to be used in figure ##
# Define the number of rows and columns for the grid
n_rows <- 5
n_cols <- 3

# Create a grid newpage
grid.newpage()

# Define the width of the color squares relative to the viewport size
square_size <- unit(0.8, "npc")

# Loop over each color and factor name
pushViewport(viewport(layout = grid.layout(n_rows, n_cols * 2)))  # Multiply cols by 2 for squares and labels
for (i in 1:length(color_mapping_I6Sphyla)) {
  row <- ((i-1) %% n_rows) + 1
  col <- ((i-1) %/% n_rows) * 2 + 1  # Determine column for the square (even cols are labels)

  # Subview for the square
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
  grid.rect(gp = gpar(fill = color_mapping_I6Sphyla[i], col = NA), width = square_size, height = square_size)
  popViewport()

  # Subview for the label
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col + 1))
  grid.text(names(color_mapping_I6Sphyla)[i], x = 0, just = "left")
  popViewport()
}

# How many samples in each day?
sampsByDate <- as.data.frame(matrix(ncol=2, nrow=length(daysOutByEU)))
sampsByDate$V1 <- daysOutByEU
for (j in 1:length(I6S_bySamplingDay_psList)){
  sampsByDate[j,2] <- ncol(otu_table(I6S_bySamplingDay_psList[[j]]))
}
sampsByDate

############# FUNGI #############
# First split phyloseq object by sampling day and remove ASVs that don't occur in subset
# This lists the daysOut by EU, with the four dates for EU 52, EU 53S, EU54S,
# and EU 8, respectively.
daysOutByEU # (MADE EARLIER ON IN SCRIPT)
ITS_bySamplingDay_psList <- vector("list", length=length(daysOutByEU)) # pre-allocate list to store results, with 1 spot per EU
for (i in 1:length(ITS_bySamplingDay_psList)){
  ITS_bySamplingDay_psList[[i]] <- subset_samples(airOnly_ITS_noAirSingsDoubs_phyloseq, DateSetOut == daysOutByEU[i])
  ITS_bySamplingDay_psList[[i]] <- prune_taxa(taxa_sums(ITS_bySamplingDay_psList[[i]]) > 0, ITS_bySamplingDay_psList[[i]]) #remove non-occurring ASVs
  names(ITS_bySamplingDay_psList)[[i]] <- unique(sample_data(ITS_bySamplingDay_psList[[i]])$DateSetOut) #this should match the daysOutByEU
}

# Double check!
names(ITS_bySamplingDay_psList) == daysOutByEU

# # MAKE order LEVEL TAXONOMY PLOTS IN A FOR LOOP
# 9/21/2025: Not double-checked since not used in manuscript (but I have no reason to think it isn't functional)

# ITS_OrderPlotInfo_list <- vector("list", length=length(ITS_bySamplingDay_psList)) # pre-allocate list to store results, with 1 spot per day
# names(ITS_OrderPlotInfo_list) <- names(ITS_bySamplingDay_psList)
# ITSuniqueOrder_byDay <- vector("list", length=length(ITS_bySamplingDay_psList)) # pre-allocate list to store unique phyla, with 1 spot per EU
# names(ITSuniqueOrder_byDay) <- names(ITS_bySamplingDay_psList)
# 
# for (j in 1:length(ITS_OrderPlotInfo_list)){
#   # Get all of the standardized relative abundances of each order in each EU
#   ITS_OrderPlotInfo_list[[j]] <- ITS_bySamplingDay_psList[[j]] %>%
#     tax_glom(taxrank = "Order") %>%
#     transform_sample_counts(function(x) x / sum(x)) %>%
#     merge_samples(group = "HabitatAir") %>%
#     transform_sample_counts(function(x) x / sum(x)) %>%
#     psmelt()
#   # Make the orders that comprise less than 5% be "< 5% abund
#   ITS_OrderPlotInfo_list[[j]]$Order[ITS_OrderPlotInfo_list[[j]]$Abundance < 0.01] <- "< 1% abundance"
#   # Make order a factor for better control later of aesthetics in plot
#   ITS_OrderPlotInfo_list[[j]]$Order <- as.factor(ITS_OrderPlotInfo_list[[j]]$Order)
#   # Store all of these unique phyla for ease of later plotting (they may differ by EU)
#   ITSuniqueOrder_byDay[[j]] <- unique(ITS_OrderPlotInfo_list[[j]]$Order)
# }
# 
# 
# # Get pretty colors
# ITSuniqueOrder_byDay #most that a given day has is 13 orders
# ITSuniqueOrder_byDayAll <- sort(unique(unlist(ITSuniqueOrder_byDay)))
# ITSuniqueOrder_byDayAll #16 unique orders
# colorsITSOrders <- glasbey.colors(18) #do less to avoid white and black
# swatch(colorsITSOrders)
# swatch(colorsITSOrders[c(2:4,6:18)]) #looks good
# colorsITSOrders <- colorsITSOrders[c(2:4,6:18)]
# swatch(colorsITSOrders)
# 
# # Clean up list before plotting
# ITS_OrderPlotInfo_list <- lapply(ITS_OrderPlotInfo_list, function(df) { #creates function without name, which takes one argument df (represents each data frame in the list as lapply iterates over the list)
#   df$Order <- gsub("Agaricomycetes_ord_Incertae_sedis", "Agaricomycetes (inc. sed.)", df$Order)
#   df$Order <- gsub("NA", "Unclassified", df$Order)
#   return(df)
# })
# 
# colorsITSOrders
# # Master color mapping for all orders
# color_mapping_ITS <- c(
#   "< 1% abundance" = "#766C95",
#   "Agaricales" = "#0000FF",
#   "Capnodiales" = "#005300",
#   "Hymenochaetales"= "#FFD300",
#   "Pleosporales" = "#009FFF",
#   "Polyporales" = "#FE8F42",
#   "Russulales" = "#1F9698",
#   "Trechisporales" = "#FF0000",
#   "Cantharellales" = "#00FFBE",
#   "Helotiales" = "#783FC1",
#   "Auriculariales" = "#9A4D42",
#   "Corticiales" = "#FFACFD",
#   "Atheliales" = "#B1CC71",
#   "Sporidiobolales" = "#F1085C",
#   "Agaricomycetes (inc. sed.)" = "#00FF00",
#   "Unclassified" = "#DD00FF"
# )
# 
# # Filter the master color mapping to include only the orders found on each date
# ITS_orderPlotColors <- vector("list", length=length(ITS_bySamplingDay_psList)) # pre-allocate list to store results, with 1 spot per day
# names(ITS_orderPlotColors) <- names(ITS_bySamplingDay_psList)
# 
# for (j in 1:length(ITS_orderPlotColors)){
#   ITS_orderPlotColors[[j]] <- color_mapping_ITS[names(color_mapping_ITS) %in% ITSuniqueOrder_byDay[[j]]]
# }
# 
# # PLOT ALL OF THESE
# # Store them in a list
# ITS_Order_PlotList <- vector("list", length=length(ITS_bySamplingDay_psList))
# names(ITS_Order_PlotList) <- names(ITS_OrderPlotInfo_list)
# daysOutByEU == names(ITS_OrderPlotInfo_list) #confirms that this is order that we wanted!
# 
# # Can do a little in a for loop:
# for (k in 1:length(ITS_Order_PlotList)) {
#   ITS_Order_PlotList[[k]] <- ggplot(data=ITS_OrderPlotInfo_list[[k]], aes(x=Sample, y=Abundance, fill=Order)) + theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) +
#     theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold")) +
#     geom_bar(aes(), stat="identity", position="fill") +
#     theme_bw() +
#     theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5)) +
#     theme(legend.text = element_text(colour="black", size = 16))  +
#     theme(legend.title = element_blank()) +
#     scale_fill_manual(values = ITS_orderPlotColors[[k]])
# }
# 
# ITS_Order_PlotList[[1]]
# ITS_Order_PlotList[[4]]
# # Another for loop for plots without legends (plotted with legends above to check everything)
# names(ITS_Order_PlotList)
# ITS_Order_PlotListNoLeg <- ITS_Order_PlotList
# for (k in 1:length(ITS_Order_PlotList)) {
#   ITS_Order_PlotListNoLeg[[k]] <- ggplot(data=ITS_OrderPlotInfo_list[[k]], aes(x=Sample, y=Abundance, fill=Order)) + theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) +
#     theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold")) +
#     geom_bar(aes(), stat="identity", position="fill") +
#     theme_bw() + scale_fill_manual(values = ITS_orderPlotColors[[k]]) +
#     theme(legend.position="none") +
#     #remove all x axis labels
#     theme(
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.title.x = element_blank()
#     ) +
#     #remove all y axis labels
#     theme(
#       axis.text.y = element_blank(),
#       axis.ticks.y = element_blank(),
#       axis.title.y = element_blank()
#     )
# 
# }
# 
# # quartz()
# do.call(grid.arrange, c(ITS_Order_PlotListNoLeg, ncol = 4, nrow = 4))
# 
# ### Make a new legend with ALL colors, to be used in figure ##
# ITSuniqueOrder_byDayAll
# # Define the number of rows and columns for the grid
# n_rows_ITS <- 5
# n_cols_ITS <- 3
# 
# # Create a grid newpage
# grid.newpage()
# 
# # Define the width of the color squares relative to the viewport size
# square_size <- unit(0.8, "npc")
# 
# # Loop over each color and factor name
# pushViewport(viewport(layout = grid.layout(n_rows_ITS, n_cols_ITS * 2)))  # Multiply cols by 2 for squares and labels
# for (i in 1:length(color_mapping_ITS)) {
#   row <- ((i-1) %% n_rows_ITS) + 1
#   col <- ((i-1) %/% n_rows_ITS) * 2 + 1  # Determine column for the square (even cols are labels)
# 
#   # Subview for the square
#   pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
#   grid.rect(gp = gpar(fill = color_mapping_ITS[i], col = NA), width = square_size, height = square_size)
#   popViewport()
# 
#   # Subview for the label
#   pushViewport(viewport(layout.pos.row = row, layout.pos.col = col + 1))
#   grid.text(names(color_mapping_ITS)[i], x = 0, just = "left")
#   popViewport()
# }
# 
# # How many samples in each day?
# sampsByDate <- as.data.frame(matrix(ncol=2, nrow=length(daysOutByEU)))
# sampsByDate$V1 <- daysOutByEU
# for (j in 1:length(ITS_bySamplingDay_psList)){
#   sampsByDate[j,2] <- ncol(otu_table(ITS_bySamplingDay_psList[[j]]))
# }
# sampsByDate

###### iv. HEAT MAP #####
############# BACTERIA #############
# (CLASS LEVEL OF TAXONOMIC RESOLUTION)
length(which(sample_data(all16Sr_noAirSingsDoubs.ps)$sampleType== "air"))
length(which(sample_data(all16Sr_noAirSingsDoubs.ps)$sampleType== "phyllosphere"))
length(which(sample_data(all16Sr_noAirSingsDoubs.ps)$sampleType== "soil"))
unique(sample_data(all16Sr_noAirSingsDoubs.ps)$PlantSpecies) #29 plant species.
# "EUPCOM" is in ITS data but not 16S

all16Sr_noAirSingsDoubs.ps.class.glom <-  phyloseq::tax_glom(all16Sr_noAirSingsDoubs.ps, taxrank = "Class") 
length(unique(tax_table(all16Sr_noAirSingsDoubs.ps.class.glom)[,3])) #103 classes
mean(unique(colSums(as.data.frame(as.matrix(otu_table(all16Sr_noAirSingsDoubs.ps.class.glom)))))) # 5483.66, ~5500, what rarefied to. Very minimum is 5424, so v. close
# Some lower because I removed singletons and doubletons, but most still 5500

# Transform sample counts on just glommed samples (unlike what we did at first)
all16Sr_noAirSingsDoubs.ps.class.0 <- transform_sample_counts(all16Sr_noAirSingsDoubs.ps.class.glom , function(x) x / sum(x) )
#ASVs in this above are just a representative from each order
unique(colSums(as.data.frame(as.matrix(otu_table(all16Sr_noAirSingsDoubs.ps.class.0))))) #all 1 now!

# Make it all into one giant dataframe
allI6Sr_relabunAllSamples.class.df <-psmelt(all16Sr_noAirSingsDoubs.ps.class.0)
#View(allI6Sr_relabunAllSamples.order.df)
sort(unique(allI6Sr_relabunAllSamples.class.df$Class)) #103 unique classes

# Which are top classes to include?
colnames(allI6Sr_relabunAllSamples.class.df)
I6S_classAbundAcrossAll <- allI6Sr_relabunAllSamples.class.df %>% 
  group_by(Class) %>% 
  summarise(classSum = sum(Abundance)) %>% 
  arrange(desc(classSum))
I6S_classAbundAcrossAll
# View(I6S_classAbundAcrossAll)
# Across air
I6S_classAbundAcrossSampleTypes <- allI6Sr_relabunAllSamples.class.df %>% 
  group_by(Class, sampleType) %>% #figure out the abundance of each class by sample type (air, phyllo, soil)
  summarise(classSumType = sum(Abundance)) %>% 
  arrange(sampleType, desc(classSumType))
I6S_classAbundAcrossSampleTypes
# View(I6S_classAbundAcrossSampleTypes)

# Will include top 35 based on abundance in air
# First, need to remove NAs, so take top 36
top36classes_byAir_I6S <- I6S_classAbundAcrossSampleTypes[1:36,]
# View(top36classes_byAir_I6S)
top36classes_byAir_I6Snames <- top36classes_byAir_I6S$Class
# remove the NA name (2.46426173 is the abundance)
top36classes_byAir_I6Snames <- top36classes_byAir_I6Snames[-which(top36classes_byAir_I6Snames== "NA")]

# Trim allI6Sr_relabunAllSamples.class.df to only have these top 35
I6S_relabunAllSamples.ClassTop35.df <- allI6Sr_relabunAllSamples.class.df[allI6Sr_relabunAllSamples.class.df$Class %in% top36classes_byAir_I6Snames,]
unique(I6S_relabunAllSamples.ClassTop35.df$Class %in% top36classes_byAir_I6Snames) #great, these now match!
# Arrange alphabetically by phylum and then within phylum, class
I6S_relabunAllSamples.ClassTop35.df <- I6S_relabunAllSamples.ClassTop35.df %>% 
  arrange(desc(Phylum), desc(Class)) #descending so that alphabetical from top of plot, not from bottom
# View(I6S_relabunAllSamples.ClassTop35.df)

# Make class into a factor. It's ordered, reverse alphabetically, by phylum and then by class within each phylum. Reverse
# alphabetically so that the alphabetically first taxa appear at the top, i.e., further up on y-axis.
I6S_relabunAllSamples.ClassTop35.df$Class <- factor(I6S_relabunAllSamples.ClassTop35.df$Class, levels = unique(I6S_relabunAllSamples.ClassTop35.df$Class))
uniqBactClass <- as.character(unique(I6S_relabunAllSamples.ClassTop35.df$Class)) #remove factor and get names 
uniqBactClass

# What phyla are each of these classes in?
bactPhylaClasses.df <- as.data.frame(matrix(ncol=2, nrow = length(uniqBactClass)))
colnames(bactPhylaClasses.df) <- c("Phylum", "Class")
for (i in 1:length(uniqBactClass)){
  bactPhylaClasses.df[i,1] <- unique(I6S_relabunAllSamples.ClassTop35.df$Phylum[which(I6S_relabunAllSamples.ClassTop35.df$Class %in% uniqBactClass[i])])
  bactPhylaClasses.df[i,2] <- uniqBactClass[i]
}
bactPhylaClasses.df #this will be added to the plot to show the phylum that each class is in!

# HEAT MAP NOT DONE RE-CHECKING, BUT NEED A BREAK, SO GRAYED OUT UNTIL CAN CHECK AND FIX
# Create a named vector for label colors
sort(unique(I6S_relabunAllSamples.ClassTop35.df$sampleType)) #"air", "phyllosphere", "soil"
I6S_relabunAllSamples.ClassTop35.df$sampleType <- factor(I6S_relabunAllSamples.ClassTop35.df$sampleType, levels = sort(unique(I6S_relabunAllSamples.ClassTop35.df$sampleType)))
colorsSampType <- c("cornflowerblue", "forestgreen", "chocolate")
colorsSampTypeDarker <- c("#004CFF", "darkgreen", "chocolate4") #darker color scheme to make colors pop more

# Make heat map!
# Comparison of different transformations
log(c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) #log
exp(c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) #exponent
sqrt(c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) #square root

#View(I6S_relabunAllSamples.ClassTop35.df)
I6S_relabunAllSamples.ClassTop35.df$AbundanceBigger <- I6S_relabunAllSamples.ClassTop35.df$Abundance + 1e-6

# Get corrected order of the samples, for example, by sampling day
colnames(I6S_relabunAllSamples.ClassTop35.df)
I6S_relabunAllSamples.ClassTop35.df$Sample
# Air
I6S_AirsamplesByDay <- I6S_relabunAllSamples.ClassTop35.df %>%
  filter(sampleType == "air") %>%
  arrange(daysOut)
dim(cbind(I6S_AirsamplesByDay$Sample,I6S_AirsamplesByDay$daysOut))
I6S_AirsamplesByDay$Sample
# Soil
I6S_SoilSamplesByEU <- I6S_relabunAllSamples.ClassTop35.df %>%
  filter(sampleType == "soil") %>%
  arrange(EU)
cbind(I6S_SoilSamplesByEU$Sample,I6S_SoilSamplesByEU$EU)
I6S_SoilSamplesByEU$Sample
# Leaf
I6S_leafSamplesBySpeciesEU <- I6S_relabunAllSamples.ClassTop35.df %>%
  filter(sampleType == "phyllosphere") %>%
  arrange(PlantSpecies, EU)
cbind(I6S_leafSamplesBySpeciesEU$Sample,I6S_leafSamplesBySpeciesEU$EU, I6S_leafSamplesBySpeciesEU$PlantSpecies)
I6S_leafSamplesBySpeciesEU$Sample

I6S_sampleOrder <- c(unique(I6S_AirsamplesByDay$Sample), unique(I6S_leafSamplesBySpeciesEU$Sample), unique(I6S_SoilSamplesByEU$Sample))
I6S_relabunAllSamples.ClassTop35.df$Sample <- factor(I6S_relabunAllSamples.ClassTop35.df$Sample, levels = c(I6S_sampleOrder))
# Edit this to make it smaller to work with
I6S_relabun_class35_small <- I6S_relabunAllSamples.ClassTop35.df[,colnames(I6S_relabunAllSamples.ClassTop35.df) %in%
                                                                   c("Sample", "Phylum", "Class", "Abundance", "sampleType")]
# Save main object so that can edit heat maps quickly if need be.
# saveRDS(I6S_relabun_class35_small, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_relabun_class35_small.df.RData") #saved 09/18/2025 on laptop
# I6S_relabunAllSamples.ClassTop35.df <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_relabun_class35_small.df.RData")
# find duplicates for the exact tile positions used by the plot
I6S_relabun_class35_small %>% count(sampleType, Sample, Class) %>% filter(n > 1)

is.factor(I6S_relabun_class35_small$Sample)
sampLabs <- c(
  air  = "bioaerosols",
  phyllosphere = "foliar surfaces",
  soil        = "soil")

# Square root -- presented in BROADN meeting AND in storyboardingPlots_Summer_2024
I6SClass_heatmap_faceted_sqrt <- ggplot(I6S_relabun_class35_small, aes(Sample, Class, fill = sampleType, alpha= Abundance)) +
  geom_tile(width = 1, height = 1,                    #fill cell (defaults to ~0.9 category)
            colour = "white", linewidth = 0.05,        #more defined grid
            linejoin = "mitre") +
  theme_bw() +
  facet_grid(. ~ sampleType, scales = "free_x", space = "free_x", labeller = labeller(sampleType = sampLabs)) +
  scale_y_discrete(expand = expansion(add = 0.5))  +
  scale_alpha_continuous(trans= "sqrt", range = c(0.01,1.0)) +
  scale_fill_manual(values = colorsSampTypeDarker, name = "sampleType") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),  #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    strip.text.x = element_text(size = 12), #increase font size of facet labels
    legend.position = "none"
  )
# quartz()
I6SClass_heatmap_faceted_sqrt #now ordered by sampling day for air, plant species for phyllosphere,
# and EU for soil
# 
# 
# soilVerrucosIndex <- intersect(which(I6S_relabun_class35_small$sampleType == "soil"), which(I6S_relabun_class35_small$Class == "Verrucomicrobiae"))
# soilKtedonoIndex <- intersect(which(I6S_relabun_class35_small$sampleType == "soil"), which(I6S_relabun_class35_small$Class == "Ktedonobacteria"))
# sort(I6S_relabun_class35_small$Abundance[soilVerrucosIndex])
# sort(I6S_relabun_class35_small$Abundance[soilVerrucosIndex])
# sort(I6S_relabun_class35_small$Abundance[soilKtedonoIndex])
# 
# 
# airVerrucosIndex <- intersect(which(I6S_relabun_class35_small$sampleType == "air"), which(I6S_relabun_class35_small$Class == "Verrucomicrobiae"))
# sort(I6S_relabun_class35_small$Abundance[airVerrucosIndex])
# 
# saveRDS(I6SClass_heatmap_faceted_sqrt, file = "Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6SClass_heatmap_faceted_sqrt.rds")
# 
# ############# FUNGI #############
# # (ORDER LEVEL OF TAXONOMIC RESOLUTION)
# length(which(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType== "air"))
# length(which(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType== "phyllosphere"))
# length(which(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType== "soil"))
# length(which(sample_data(allITSr_noAirSingsDoubs.ps)$sampleType== "soil"))
# unique(sample_data(allITSr_noAirSingsDoubs.ps)$PlantSpecies) #30 plant species
# 
# allITSr_noAirSingsDoubs.ps.order.glom <-  phyloseq::tax_glom(allITSr_noAirSingsDoubs.ps, taxrank = "Order") 
# length(unique(tax_table(allITSr_noAirSingsDoubs.ps.order.glom)[,3])) #133 orders
# unique(colSums(as.data.frame(as.matrix(otu_table(allITSr_noAirSingsDoubs.ps.order.glom))))) # ~8500, what rarefied to. 
# # Some lower because I removed singletons and doubletons
# 
# # Transform sample counts on just glommed samples 
# allITSr_noAirSingsDoubs.ps.order.0 <- transform_sample_counts(allITSr_noAirSingsDoubs.ps.order.glom , function(x) x / sum(x) )
# #ASVs in this above are just a representative from each order
# unique(colSums(as.data.frame(as.matrix(otu_table(allITSr_noAirSingsDoubs.ps.order.0))))) #all 1 now!
# 
# # Make it all into one dataframe
# allITSr_relabunAllSamples.order.df <-psmelt(allITSr_noAirSingsDoubs.ps.order.0)
# #View(allITSr_relabunAllSamples.order.df)
# unique(allITSr_relabunAllSamples.order.df$Order) #133 unique orders
# 
# # Tidy up the dataframe a little bit before plotting
# unique(allITSr_relabunAllSamples.order.df$Order) #fixed
# # Clean up name of one of the orders' names 
# allITSr_relabunAllSamples.order.df$Order <- gsub(pattern= "_ord_Incertae_sedis", x=allITSr_relabunAllSamples.order.df$Order, replacement = ", inc. sed.")
# # Remove NAs
# allITSr_relabunAllSamples.order.df <- allITSr_relabunAllSamples.order.df[-which(allITSr_relabunAllSamples.order.df$Order == "NA"),]
# sort(unique(allITSr_relabunAllSamples.order.df$Order)) #fixed
# 
# # Which are top orders to include?
# colnames(allITSr_relabunAllSamples.order.df)
# ITS_orderAbundAcrossAll <- allITSr_relabunAllSamples.order.df %>% 
#   group_by(Order) %>% 
#   summarise(orderSum = sum(Abundance)) %>% 
#   arrange(desc(orderSum))
# ITS_orderAbundAcrossAll
# #View(ITS_orderAbundAcrossAll)
# # Across air
# ITS_orderAbundAcrossSampleTypes <- allITSr_relabunAllSamples.order.df %>% 
#   group_by(Order, sampleType) %>% 
#   summarise(orderSumType = sum(Abundance)) %>% 
#   arrange(sampleType, desc(orderSumType))
# ITS_orderAbundAcrossSampleTypes
# # View(ITS_orderAbundAcrossSampleTypes)
# 
# # Will include top 35 based on abundance in air
# ITS_top35Ords_byAir <- ITS_orderAbundAcrossSampleTypes[1:35,]
# # View(ITS_top35Ords_byAir)
# ITS_top35Ords_byAirNames <- ITS_top35Ords_byAir$Order
# ITS_top35Ords_byAirNames
# 
# # Trim allITSr_relabunAllSamples.order.df to only have these top 35
# class(allITSr_relabunAllSamples.order.df$Order)
# ITS_relabunAllSamples.orderTop35.df <- allITSr_relabunAllSamples.order.df[allITSr_relabunAllSamples.order.df$Order %in% ITS_top35Ords_byAirNames,]
# unique(ITS_relabunAllSamples.orderTop35.df$Order %in% ITS_top35Ords_byAirNames) #great, these now match!
# # Arrange alphabetically by phylum and then within class, order
# ITS_relabunAllSamples.orderTop35.df <- ITS_relabunAllSamples.orderTop35.df %>% 
#   arrange(desc(Class), desc(Order))
# # View(ITS_relabunAllSamples.orderTop35.df)
# 
# # Make order into a factor. It's ordered, reverse alphabetically, by Order and then by classes within each order. Reverse
# # alphabetically so that the alphabetically first taxa appear at the top, i.e., further up on y-axis.
# ITS_relabunAllSamples.orderTop35.df$Order <- factor(ITS_relabunAllSamples.orderTop35.df$Order, levels = unique(ITS_relabunAllSamples.orderTop35.df$Order))
# uniqFungOrders <- as.character(unique(ITS_relabunAllSamples.orderTop35.df$Order)) #remove factor and get names 
# uniqFungOrders
# 
# # What classes are each of these orders in?
# fungClassOrders.df <- as.data.frame(matrix(ncol=2, nrow = length(uniqFungOrders)))
# colnames(fungClassOrders.df) <- c("Class", "Order")
# for (i in 1:length(uniqFungOrders)){
#   fungClassOrders.df[i,1] <- unique(ITS_relabunAllSamples.orderTop35.df$Class[which(ITS_relabunAllSamples.orderTop35.df$Order %in% uniqFungOrders[i])])
#   fungClassOrders.df[i,2] <- uniqFungOrders[i]
# }
# fungClassOrders.df #this will be added to the plot to show class!
# length(unique(fungClassOrders.df$Order))
# 
# ## Create a named vector for label colors
# ITS_relabunAllSamples.orderTop35.df$sampleType <- factor(ITS_relabunAllSamples.orderTop35.df$sampleType, levels = sort(unique(ITS_relabunAllSamples.orderTop35.df$sampleType)))
# 
# # Get corrected order of the samples, for example, by sampling day
# colnames(ITS_relabunAllSamples.orderTop35.df)
# ITS_relabunAllSamples.orderTop35.df$Sample
# # Air
# ITS_AirsamplesByDay <- ITS_relabunAllSamples.orderTop35.df %>% 
#   filter(sampleType == "air") %>% 
#   arrange(daysOut)
# dim(cbind(ITS_AirsamplesByDay$Sample,ITS_AirsamplesByDay$daysOut))
# ITS_AirsamplesByDay$Sample
# # Soil
# ITS_SoilSamplesByEU <- ITS_relabunAllSamples.orderTop35.df %>% 
#   filter(sampleType == "soil") %>% 
#   arrange(EU)
# cbind(ITS_SoilSamplesByEU$Sample,ITS_SoilSamplesByEU$EU)
# ITS_SoilSamplesByEU$Sample
# # Leaf
# ITS_leafSamplesBySpeciesEU <- ITS_relabunAllSamples.orderTop35.df %>% 
#   filter(sampleType == "phyllosphere") %>% 
#   arrange(PlantSpecies, EU)
# cbind(ITS_leafSamplesBySpeciesEU$Sample,ITS_leafSamplesBySpeciesEU$EU, ITS_leafSamplesBySpeciesEU$PlantSpecies)
# ITS_leafSamplesBySpeciesEU$Sample
# 
# ITS_sampleOrder <- c(unique(ITS_AirsamplesByDay$Sample), unique(ITS_leafSamplesBySpeciesEU$Sample), unique(ITS_SoilSamplesByEU$Sample))
# ITS_relabunAllSamples.orderTop35.df$Sample <- factor(ITS_relabunAllSamples.orderTop35.df$Sample, levels = c(ITS_sampleOrder))
# 
# # Save main object so that can edit heat maps quickly if need be. 
# # saveRDS(ITS_relabunAllSamples.orderTop35.df, file = "ITS_relabunAllSamples.orderTop35.df.RData") #saved 12/2/2024
# 
# # Make heat map!
# # Square Root
# ITSOrder_heatmap_faceted_Sqrt <- ggplot(ITS_relabunAllSamples.orderTop35.df, aes(Sample, Order, fill = sampleType, alpha= Abundance)) +
#   geom_tile() +
#   theme_bw() + 
#   facet_grid(. ~ sampleType, scales = "free_x", space = "free") +
#   scale_alpha_continuous(trans= "sqrt", range = c(0.0,1.0)) +
#   scale_fill_manual(values = colorsSampTypeDarker, name = "sampleType") +
#   #scale_fill_brewer(palette = "Set1"
#   #       , name = "sampleType") +
#   theme(
#     axis.text.x = element_blank(),  
#     axis.ticks.x = element_blank(), 
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     legend.title = element_blank(),
#     panel.grid.major = element_blank(),  #remove major gridlines
#     panel.grid.minor = element_blank(), #remove minor gridlines
#     strip.text.x = element_text(size = 12), #increase font size of facet labels
#     legend.position = "none"
#     )
# 
# # quartz()
# ITSOrder_heatmap_faceted_Sqrt

############# BACTERIA #############
airOnly_16Sr_noAirSingsDoubs_phyloseq #phyloseq object to use that only has air data 
# First, to make subsetting easier, add a column with sampleID (same as the rownames in the metadata)
air_I6S_IDs.ps <- airOnly_16Sr_noAirSingsDoubs_phyloseq #make a copy just in case so it doesn't mess anything up downstream
sample_data(air_I6S_IDs.ps)$sampleID <- rownames(sample_data(air_I6S_IDs.ps))
sample_data(air_I6S_IDs.ps) #looks good!
# Because of needing to restrict permutations, will need to make data within days equal

# CHECK NUMBERS IN EACH
I6Sair_meta <- as.data.frame(as.matrix(sample_data(air_I6S_IDs.ps)))
colnames(I6Sair_meta)
rownames(I6Sair_meta) #sample IDs are rows 
table(I6Sair_meta$HabitatAir, I6Sair_meta$DateSetOut) #shows def not enough to do it by day
#          1-Jul-2022 16-Jun-2022 18-Jun-2022 2-Jul-2022 20-Jun-2022 22-Jun-2022 23-Jun-2022 24-Jun-2022 26-Jun-2022 27-Jun-2022
# forest           1           2           4          4           4           4           3           2           2           2
# savanna          1           3           3          1           2           2           3           2           3           2
# 
#         29-Jun-2022 30-Jun-2022 4-Jul-2022 5-Jul-2022 6-Jul-2022 8-Jul-2022
# forest            3           4          2          4          2          3
# savanna           3           3          4          2          2          2

table(I6Sair_meta$HabitatAir, I6Sair_meta$EU)
#         EU_52 EU_53S EU_54S EU_8
# forest      9     12     13   12
# savanna    10      7      9   12

# RANDOMLY SUBSET SAMPLES FOR 7 SAMPLES FROM EACH HABITAT FROM EACH OF THE 4 EUS
set.seed(1121) #to make reproducible
I6Sair_meta_7hab <- I6Sair_meta %>%
  group_by(HabitatAir, EU) %>%
  sample_n(size = 7) %>%
  ungroup()
# View(I6Sair_meta_7hab) 56 entries = 7 SAMPLES *2 HABTIAT AIRS *4 EUS
# DOUBLE CHECK
table(I6Sair_meta_7hab$HabitatAir, I6Sair_meta_7hab$EU) #all 7

# MAKE NEW PHYLOSEQ OBJECTs WITH JUST THESE SAMPLES 
colnames(I6Sair_meta_7hab)
I6Sair_meta_7hab$sampleID; length(I6Sair_meta_7hab$sampleID) #sample names as expected and there are 56
colnames(sample_data(air_I6S_IDs.ps))
# Make new phyloseq object
airOnly_16S_7hab.ps  <- air_I6S_IDs.ps %>%
  subset_samples(sampleID %in% I6Sair_meta_7hab$sampleID)
# Double check-- these all look good!
table(sample_data(airOnly_16S_7hab.ps)$EU, sample_data(airOnly_16S_7hab.ps)$HabitatAir)

# GET THE NECESSARY PIECES OF DATA FROM PHYLOSEQ OBJECT:
# SET UP ASV TABLE
# First, need to get data out of phyloseq so that it can be processed in vegan
airOnly_16S_7hab_ASVs <- t(as.data.frame(as.matrix(otu_table(airOnly_16S_7hab.ps)))) #get ASV table
dim(airOnly_16S_7hab_ASVs) #56 rows and 4105 columns, which means that samples are rows, taxa are columns as they should be for vegan
length(which(colSums(airOnly_16S_7hab_ASVs)== 0)) #741 were dropped/lost by subsetting these samples
# Remove these non-occurring ASVs
airOnly_16S_7hab_ASVs <- airOnly_16S_7hab_ASVs[,-which(colSums(airOnly_16S_7hab_ASVs)== 0)]
length(which(colSums(airOnly_16S_7hab_ASVs)== 0))
# View(airOnly_16S_7hab_ASVs)
ncol(airOnly_16S_7hab_ASVs) == (4105- 741) #number of expected ASVs removed

# SET UP METADATA
airOnly_16S_7hab_sampNames <- rownames(sample_data(airOnly_16S_7hab.ps))
unique(rownames(airOnly_16S_7hab_ASVs) == airOnly_16S_7hab_sampNames) #since these match, we are good to go!
colnames(sample_data(airOnly_16S_7hab.ps))
# Make a new metadata object
airOnly_16S_7hab_metadata <- as.data.frame(cbind(sample_data(airOnly_16S_7hab.ps)$sampleID, sample_data(airOnly_16S_7hab.ps)$HabitatAir, 
                                                 sample_data(airOnly_16S_7hab.ps)$EU, sample_data(airOnly_16S_7hab.ps)$DateSetOut))
colnames(airOnly_16S_7hab_metadata) <-c("sampleID", "HabitatAir", "EU", "DateSetOut")
str(airOnly_16S_7hab_metadata)
head(airOnly_16S_7hab_metadata)
(airOnly_16S_7hab_metadata$sampleID == airOnly_16S_7hab_sampNames) #also looks correct!
# Make variables factors as needed:
airOnly_16S_7hab_metadata$HabitatAir <- as.factor(airOnly_16S_7hab_metadata$HabitatAir)
airOnly_16S_7hab_metadata$EU <- as.factor(airOnly_16S_7hab_metadata$EU)
dim(airOnly_16S_7hab_metadata) #56, 4

# Next get Bray-Curtis dissimilarities
airOnly_16S_7hab_16SBC <- vegdist(airOnly_16S_7hab_ASVs, method= "bray")

# MAKE A CONTROL OBJECT TO USE TO SET PERMUTATIONS- this permutes only within each EU
CTRL_I6SairHab_7 <- how(within = Within(type = "free"),
                        plots = Plots(strata = airOnly_16S_7hab_metadata$EU, type = "none"),
                        nperm = 9999,
                        observed = TRUE)

set.seed(0719) #set seed so that results are reproducible!
I6S_airHab_7_permanova <- adonis2(airOnly_16S_7hab_16SBC ~ airOnly_16S_7hab_metadata$EU + airOnly_16S_7hab_metadata$HabitatAir, permutations= CTRL_I6SairHab_7, by = "terms")
I6S_airHab_7_permanova

# NOT SIGNIFICANT, NO EFFECT OF EU EITHER
# adonis2(formula = airOnly_16S_7hab_16SBC ~ airOnly_16S_7hab_metadata$EU + airOnly_16S_7hab_metadata$HabitatAir, permutations = CTRL_I6SairHab, by = "terms")
#                                      Df SumOfSqs      R2      F Pr(>F)
# airOnly_16S_7hab_metadata$EU          3   1.7082 0.06602 1.2246 0.5707
# airOnly_16S_7hab_metadata$HabitatAir  1   0.4549 0.01758 0.9784 0.5707
# Residual                             51  23.7125 0.91640              
# Total                                55  25.8756 1.00000  


# # WHAT IF WE DID 9 AND DROPPED EU EU_53S?
# NOT USED IN SUBMISSION, SO GRAYED OUT AND NOT CLOSELY RE-CHECKED
# set.seed(1121) #to make reproducible
# I6Sair_meta_9hab <- I6Sair_meta %>%
#   filter(EU != "EU_53S") %>% 
#   group_by(HabitatAir, EU) %>%
#   sample_n(size = 9) %>%
#   ungroup()
# # View(I6Sair_meta_9hab)
# # DOUBLE CHECK
# table(I6Sair_meta_9hab$HabitatAir, I6Sair_meta_9hab$EU) #all 9
# # View(I6Sair_meta_9hab_no53S)
# 
# # MAKE NEW PHYLOSEQ OBJECTs WITH JUST THESE SAMPLES 
# colnames(I6Sair_meta_9hab)
# I6Sair_meta_9hab$sampleID #sample names as expected
# colnames(sample_data(air_I6S_IDs.ps))
# # Make new phyloseq object
# airOnly_16S_9hab.ps  <- air_I6S_IDs.ps %>%
#   subset_samples(sampleID %in% I6Sair_meta_9hab$sampleID)
# # Double check-- these all look good!
# table(sample_data(airOnly_16S_9hab.ps)$EU, sample_data(airOnly_16S_9hab.ps)$HabitatAir)
# 
# # GET THE NECESSARY PIECES OF DATA FROM PHYLOSEQ OBJECT:
# # SET UP ASV TABLE
# # First, need to get data out of phyloseq so that it can be processed in vegan
# airOnly_16S_9hab_ASVs <- t(as.data.frame(as.matrix(otu_table(airOnly_16S_9hab.ps)))) #get ASV table
# dim(airOnly_16S_9hab_ASVs) #54 rows and 4105 columns, which means that samples are rows, taxa are columns as they should be for vegan
# length(which(colSums(airOnly_16S_9hab_ASVs)== 0)) #785 were dropped/lost by subsetting these samples
# # Remove these ASVs
# airOnly_16S_9hab_ASVs <- airOnly_16S_9hab_ASVs[,-which(colSums(airOnly_16S_9hab_ASVs)== 0)]
# length(which(colSums(airOnly_16S_9hab_ASVs)== 0))
# # View(airOnly_16S_9hab_ASVs)
# ncol(airOnly_16S_9hab_ASVs) == (4105- 785) #number of expected ASVs removed
# 
# # SET UP METADATA
# airOnly_16S_9hab_sampNames <- rownames(sample_data(airOnly_16S_9hab.ps))
# unique(rownames(airOnly_16S_9hab_ASVs) == airOnly_16S_9hab_sampNames) #since these match, we are good to go!
# colnames(sample_data(airOnly_16S_9hab.ps))
# # Make a new metadata object
# airOnly_16S_9hab_metadata <- as.data.frame(cbind(sample_data(airOnly_16S_9hab.ps)$sampleID, sample_data(airOnly_16S_9hab.ps)$HabitatAir, 
#                                                  sample_data(airOnly_16S_9hab.ps)$EU, sample_data(airOnly_16S_9hab.ps)$DateSetOut))
# colnames(airOnly_16S_9hab_metadata) <-c("sampleID", "HabitatAir", "EU", "DateSetOut")
# str(airOnly_16S_9hab_metadata)
# head(airOnly_16S_9hab_metadata)
# # Make variables factors as needed:
# airOnly_16S_9hab_metadata$HabitatAir <- as.factor(airOnly_16S_9hab_metadata$HabitatAir)
# airOnly_16S_9hab_metadata$EU <- as.factor(airOnly_16S_9hab_metadata$EU)
# dim(airOnly_16S_9hab_metadata) #54, 4
# 
# # Next get Bray-Curtis dissimilarities
# airOnly_16S_9hab_16SBC <- vegdist(airOnly_16S_9hab_ASVs, method= "bray")
# 
# # MAKE A CONTROL OBJECT TO USE TO SET PERMUTATIONS- this permutes only within each EU
# CTRL_I6SairHab_9 <- how(within = Within(type = "free"),
#                         plots = Plots(strata = airOnly_16S_9hab_metadata$EU, type = "none"),
#                         nperm = 9999,
#                         observed = TRUE)
# 
# set.seed(0719) #set seed so that results are reproducible!
# I6S_airHab_9_permanova <- adonis2(airOnly_16S_9hab_16SBC ~ airOnly_16S_9hab_metadata$EU + airOnly_16S_9hab_metadata$HabitatAir, permutations= CTRL_I6SairHab_9, by = "terms")
# I6S_airHab_9_permanova
# # adonis2(formula = airOnly_16S_9hab_16SBC ~ airOnly_16S_9hab_metadata$EU + airOnly_16S_9hab_metadata$HabitatAir, permutations = CTRL_I6SairHab_9, by = "terms")
# # Df SumOfSqs      R2      F Pr(>F)
# # airOnly_16S_9hab_metadata$EU          2   1.1614 0.04668 1.2468 0.8555
# # airOnly_16S_9hab_metadata$HabitatAir  1   0.4311 0.01733 0.9256 0.8555
# # Residual                             50  23.2868 0.93599              
# # Total                                53  24.8793 1.00000    


############# FUNGI #############
airOnly_ITS_noAirSingsDoubs_phyloseq #phyloseq object to use.
# First, to make subsetting easier, add a column with sampleID (same as the rownames in the metadata)
air_ITS_IDs.ps <- airOnly_ITS_noAirSingsDoubs_phyloseq #make a copy just in case so it doesn't mess anything up downstream
sample_data(air_ITS_IDs.ps)$sampleID <- rownames(sample_data(air_ITS_IDs.ps))
sample_data(air_ITS_IDs.ps) #looks good!
# Because of needing to restrict permutations, will need to make data within days equal

# CHECK NUMBERS IN EACH
ITSair_meta <- as.data.frame(as.matrix(sample_data(air_ITS_IDs.ps)))
colnames(ITSair_meta)
rownames(ITSair_meta) #sample IDs are rows 
table(ITSair_meta$HabitatAir, ITSair_meta$DateSetOut) #could probably do it within day if omitting all days with less than 3 samples in 
# matrix or patch (omit 16-Jun-2022, 22-Jun-2022, 27-Jun-2022, and 5-Jul-2022)
#         1-Jul-2022 16-Jun-2022 18-Jun-2022 2-Jul-2022 20-Jun-2022 22-Jun-2022 23-Jun-2022 24-Jun-2022 26-Jun-2022 27-Jun-2022 29-Jun-2022 30-Jun-2022 4-Jul-2022 5-Jul-2022
# forest           4           2           4          4           4           4           3           4           3           3           4           4          4          4
# savanna          4           3           3          4           3           2           3           3           4           2           3           4          4          2
# 
# 6-Jul-2022 8-Jul-2022
# forest           4          4
# savanna          4          3

table(ITSair_meta$HabitatAir, ITSair_meta$EU) #shows minimum any day/habitat type is 12
# EU_52 EU_53S EU_54S EU_8
# forest     13     15     16   15
# savanna    14     12     12   13

# RANDOMLY SUBSET SAMPLES FOR 12 SAMPLES FROM EACH HABITAT
set.seed(1121) #to make reproducible
ITSair_meta_12hab <- ITSair_meta %>%
  group_by(HabitatAir, EU) %>%
  sample_n(size = 12) %>%
  ungroup()
# View(ITSair_meta_12hab)
# DOUBLE CHECK
table(ITSair_meta_12hab$HabitatAir, ITSair_meta_12hab$EU) #all 12

# MAKE NEW PHYLOSEQ OBJECTs WITH JUST THESE SAMPLES 
colnames(ITSair_meta_12hab)
ITSair_meta_12hab$sampleID #sample names as expected
colnames(sample_data(air_ITS_IDs.ps))
# Make new phyloseq object
airOnly_ITS_12hab.ps  <- air_ITS_IDs.ps %>%
  subset_samples(sampleID %in% ITSair_meta_12hab$sampleID)
# Double check-- these all look good!
table(sample_data(airOnly_ITS_12hab.ps)$EU, sample_data(airOnly_ITS_12hab.ps)$HabitatAir)
unique(sort(sample_data(airOnly_ITS_12hab.ps)$sampleID) == sort(ITSair_meta_12hab$sampleID))

# GET THE NECESSARY PIECES OF DATA FROM PHYLOSEQ OBJECT:
# SET UP ASV TABLE
# First, need to get data out of phyloseq so that it can be processed in vegan
airOnly_ITS_12hab_ASVs <- t(as.data.frame(as.matrix(otu_table(airOnly_ITS_12hab.ps)))) #get ASV table
dim(airOnly_ITS_12hab_ASVs) #96 rows and 4673 columns, which means that samples are rows, taxa are columns as they should be for vegan
length(which(colSums(airOnly_ITS_12hab_ASVs)== 0)) #90 ASVs were dropped/lost by subsetting these samples
# Remove these ASVs with no reads
airOnly_ITS_12hab_ASVs <- airOnly_ITS_12hab_ASVs[,-which(colSums(airOnly_ITS_12hab_ASVs)== 0)]
length(which(colSums(airOnly_ITS_12hab_ASVs)== 0)) #all good
# View(airOnly_ITS_12hab_ASVs)
ncol(airOnly_ITS_12hab_ASVs) == (4673- 90) #number of expected ASVs removed

# SET UP METADATA
airOnly_ITS_12hab_sampNames <- rownames(sample_data(airOnly_ITS_12hab.ps))
unique(rownames(airOnly_ITS_12hab_ASVs) == airOnly_ITS_12hab_sampNames) #since these match, we are good to go!
colnames(sample_data(airOnly_ITS_12hab.ps))
# Make a new metadata object
airOnly_ITS_12hab_metadata <- as.data.frame(cbind(sample_data(airOnly_ITS_12hab.ps)$sampleID, sample_data(airOnly_ITS_12hab.ps)$HabitatAir, 
                                                  sample_data(airOnly_ITS_12hab.ps)$EU, sample_data(airOnly_ITS_12hab.ps)$DateSetOut))
colnames(airOnly_ITS_12hab_metadata) <-c("sampleID", "HabitatAir", "EU", "DateSetOut")
str(airOnly_ITS_12hab_metadata)
head(airOnly_ITS_12hab_metadata)
# Make variables factors as needed:
airOnly_ITS_12hab_metadata$HabitatAir <- as.factor(airOnly_ITS_12hab_metadata$HabitatAir)
airOnly_ITS_12hab_metadata$EU <- as.factor(airOnly_ITS_12hab_metadata$EU)
dim(airOnly_ITS_12hab_metadata) #96, 4, as expected

# Next get Bray-Curtis dissimilarities
airOnly_ITS_12hab_ITSBC <- vegdist(airOnly_ITS_12hab_ASVs, method= "bray")

# MAKE A CONTROL OBJECT TO USE TO SET PERMUTATIONS- this permutes only within each EU
CTRL_ITSairHab_12 <- how(within = Within(type = "free"),
                         plots = Plots(strata = airOnly_ITS_12hab_metadata$EU, type = "none"),
                         nperm = 9999,
                         observed = TRUE)

set.seed(0719) #set seed so that results are reproducible!
ITS_airHab_12_permanova <- adonis2(airOnly_ITS_12hab_ITSBC ~ airOnly_ITS_12hab_metadata$EU + airOnly_ITS_12hab_metadata$HabitatAir, permutations= CTRL_ITSairHab_12, by = "terms")
ITS_airHab_12_permanova

# adonis2(formula = airOnly_ITS_12hab_ITSBC ~ airOnly_ITS_12hab_metadata$EU + airOnly_ITS_12hab_metadata$HabitatAir, permutations = CTRL_ITSairHab_12, by = "terms")
#                                         Df SumOfSqs      R2      F Pr(>F)
# airOnly_ITS_12hab_metadata$EU          3   3.8737 0.18571 6.9629 0.8077
# airOnly_ITS_12hab_metadata$HabitatAir  1   0.1091 0.00523 0.5883 0.8077
# Residual                              91  16.8753 0.80905              
# Total                                 95  20.8581 1.00000   


###############################################
# III.  QPCR PLOTS AND STATS
###############################################
prok1 <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/claire_qPCR/QuantSummariesUsed/ProkPlate1qPCRForR.csv")
prok2 <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/claire_qPCR/QuantSummariesUsed/ProkPlate2qPCRForR.csv")
fung1 <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/claire_qPCR/QuantSummariesUsed/fungiPlate1_qPCRforR.csv")
fung2 <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/claire_qPCR/QuantSummariesUsed/fungiPlate2_qPCRforR.csv")
head(prok1)
head(prok2)
head(fung1)
head(fung2)

################# 16S #################
### ONLY SAMPLES (no controls) ###
# Get sample data for samples used so that can subset qPCR samples to match
I6SsampDat <- as.data.frame(as.matrix(sample_data(all16Sr_noAirSingsDoubs.ps)))
head(I6SsampDat)
rownames(I6SsampDat) #sample names
unique(I6SsampDat$isControl) #only has samples (no controls were retained after rarefying)
I6SsampNames <- rownames(I6SsampDat)
I6SairNames <- I6SsampNames[grepl(x=I6SsampNames, pattern="air")] #84 as expected

# Filter the I6S qPCR by the air samples that were retained after rarefying
I6SairNamesNumbers <- gsub(x=I6SairNames, pattern = "air_16S_", replacement= "")
# Plate 1
I6SqPCR_Plate1_trimmed <- prok1[prok1$sampleName %in% I6SairNamesNumbers,]
# Plate 2
I6SqPCR_Plate2_trimmed <- prok2[prok2$sampleName %in% I6SairNamesNumbers,]

# Check to make sure that the column names are the same for the two I6S plates
colnames(I6SqPCR_Plate1_trimmed[1:7]) == colnames(I6SqPCR_Plate2_trimmed)
# Bind plate 1 and plate 2 together to make 84 rows representing all remaining bacterial samples
I6SqPCRdata <- rbind(I6SqPCR_Plate1_trimmed[,1:7], I6SqPCR_Plate2_trimmed)
nrow(I6SqPCRdata)

# Make a copy and re-organize so that it can be put in table in the supplementary data
I6SqPCRdata_forTable <- I6SqPCRdata
# Make sample name numeric so that it can be arranged!
I6SqPCRdata_forTable$sampleName <- as.numeric(I6SqPCRdata_forTable$sampleName)
# Round these to the nearest whole number (which was ultimately also done in the statistical analyses!)
I6SqPCRdata_forTable$GenomeEquiv <- round(I6SqPCRdata_forTable$GenomeEquiv)
I6SqPCRdata_forTable <- I6SqPCRdata_forTable %>% 
  arrange(sampleName)  #arrange by sample name
# View(I6SqPCRdata_forTable)
# write.csv(I6SqPCRdata_forTable, file= "I6SqPCRdata_forTable.csv") #written Nov 13, 2024 and exported to own computer 
# to make supplemental figure (but note not ultimate form used since I will bind bacterial and fungal later)
# Plot the non-control samples (only those that made it through filters employed 
# earlier)
#quartz()
plot(sort(I6SqPCRdata$GenomeEquiv))

# ### BRIEFLY CHECK OUT CONTROLS ### to fix this section, save  all16S_rarefied.ps in full16S_EDA_part1_Sept.R and then load here
# all16Sr_Controls.ps <- subset_samples(all16S_rarefied.ps, isControl == "control") #new phyloseq object without controls/blanks
# 
# # Get control sample names-- 
# I6SCONTROLSsampDat <- as.data.frame(as.matrix(sample_data(all16Sr_Controls.ps)))
# head(I6SCONTROLSsampDat)
# rownames(I6SCONTROLSsampDat)
# I6SsampNamesCONTROLS <- rownames(I6SCONTROLSsampDat)
# I6SairCONTROLNames <- I6SsampNamesCONTROLS[grepl(x=I6SsampNamesCONTROLS, pattern="air")]
# I6SairNamesCONTROLSNumbers <- gsub(x=I6SairCONTROLNames, pattern = "air_16S_", replacement= "")
# I6SqPCR_Plate1_trimmedCONTROLS <- prok1[prok1$sampleName %in% I6SairNamesCONTROLSNumbers,]
# I6SqPCR_Plate1_trimmedCONTROLS #only plate 1 has some controls
# I6SqPCR_Plate2_trimmedCONTROLS <- prok2[prok2$sampleName %in% I6SairNamesCONTROLSNumbers,] 
# I6SqPCR_Plate2_trimmedCONTROLS #zero, as assumed
# # Most look pretty low!
# plot(sort(I6SqPCR_Plate1_trimmedCONTROLS$GenomeEquiv)) #controls were very low, but some did remain

# ### PLOT 16S QPCR DATA BY DAY! (to match the taxonomy plots above) ### to use this again, bring back in by day analyses grayed out
# earlier in this script
head(I6SqPCRdata)
# This list, made earlier in this script, has all of the samples by day. Need to get the sample
# names from each and use this to organize the qPCR samples
I6S_bySamplingDay_psList
# pre-allocate a list to store all of the sample names and HabitatAir info by day
I6S_byDaySampNamesHabitat <- vector("list", length=length(I6S_bySamplingDay_psList))
names(I6S_byDaySampNamesHabitat) <- names(I6S_bySamplingDay_psList)
# Get sample names for each day
for (j in 1:length(I6S_byDaySampNamesHabitat)){
  # get the sample names and their habitat air information
  I6S_byDaySampNamesHabitat[[j]] <- cbind(rownames(sample_data(I6S_bySamplingDay_psList[[j]])), sample_data(I6S_bySamplingDay_psList[[j]])$HabitatAir)
  # Replace the full sample names with the shorter ones that match the qPCR labels
  I6S_byDaySampNamesHabitat[[j]] <- gsub(x=I6S_byDaySampNamesHabitat[[j]], pattern = "air_16S_", replacement = "")
  I6S_byDaySampNamesHabitat[[j]] <- as.data.frame(as.matrix(I6S_byDaySampNamesHabitat[[j]]))
  colnames(I6S_byDaySampNamesHabitat[[j]]) <- c("sampleName", "HabitatAir")
}
# Tests to show above for loop is working as expected with one day
test1 <- cbind(rownames(sample_data(I6S_bySamplingDay_psList[[9]])), sample_data(I6S_bySamplingDay_psList[[9]])$HabitatAir)
test2 <- gsub(x=test1, pattern = "air_16S_", replacement = "")
test2 <-as.data.frame(as.matrix(test2))
test2 ==I6S_byDaySampNamesHabitat[[9]]
colnames(test2) <- c("sampleName", "HabitatAir")
test2 == I6S_byDaySampNamesHabitat[[9]]

I6S_byDaySampNamesHabitat[[9]]

# GET qPCR DATA SEPARATED BY DAY
# Pre-allocate a new list to hold each days qPCR info
I6S_qPCR_dayList <- vector("list", length(I6S_byDaySampNamesHabitat))
names(I6S_qPCR_dayList) <- names(I6S_byDaySampNamesHabitat)

# Loop through each element in list1
for (i in 1:length(I6S_byDaySampNamesHabitat)) {
  # Extract the sample names from the current element in the list
  sample_names <- I6S_byDaySampNamesHabitat[[i]]$sampleName

  # Subset df1 to only include rows where sampleName matches those in the current list element
  day_subset <- I6SqPCRdata[I6SqPCRdata$sampleName %in% sample_names, ]

  # Merge the subset of df1 with the current element of list1 based on sampleName
  # Using all.x = TRUE to keep all rows from the list1 data frame in case there are unmatched sampleNames
  I6S_qPCR_dayList[[i]] <- merge(I6S_byDaySampNamesHabitat[[i]], day_subset, by = "sampleName", all.x = TRUE)
}

# Spot check
I6S_qPCR_dayList[[9]] #looks correct!
I6S_byDaySampNamesHabitat[[9]] #looks correct!

# # Re-order these based on the order that I used for the earlier plots
# # Reorder merged_dfs according to the desired order
# I6S_qPCR_dayList <- I6S_qPCR_dayList[names(ITS_Order_PlotList)]
# names(I6S_qPCR_dayList) == names(ITS_Order_PlotList)
# 
# PLOT BOXPLOTS!!!
# Store them in a list
I6S_qPCR_PlotList <- vector("list", length=length(I6S_qPCR_dayList))
names(I6S_qPCR_dayList) <- names(I6S_qPCR_dayList)
length(I6S_qPCR_dayList)

# Can do a little in a for loop:
for (k in 1:length(I6S_qPCR_PlotList)) {
  I6S_qPCR_PlotList[[k]] <- ggplot(data=I6S_qPCR_dayList[[k]], aes(x=HabitatAir, y=GenomeEquiv, fill=HabitatAir)) +
    geom_boxplot() +
    theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) +
    scale_fill_manual(values = c("forest" = "forestgreen", "savanna" = "goldenrod")) +
    geom_point(color="black", size=2, alpha=0.9) +
    theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold")) +
    geom_bar(aes(), stat="identity", position="fill") +
    theme_bw() +
    theme(legend.position="none")
}

# PLOT THEM ALL
I6S_qPCR_PlotList[[11]]
I6S_qPCR_dayList[[11]]
# quartz()
do.call(gridExtra::grid.arrange, c(I6S_qPCR_PlotList, ncol = 4, nrow = 4))

aggregate(GenomeEquiv ~ HabitatAir ,data = I6S_qPCR_dayList[[1]], mean)
aggregate(GenomeEquiv ~ HabitatAir ,data = I6S_qPCR_dayList[[2]], mean)
aggregate(GenomeEquiv ~ HabitatAir ,data = I6S_qPCR_dayList[[5]], mean)
aggregate(GenomeEquiv ~ HabitatAir ,data = I6S_qPCR_dayList[[8]], mean)
aggregate(GenomeEquiv ~ HabitatAir ,data = I6S_qPCR_dayList[[10]], mean)
I6S_qPCR_dayList[[13]]
aggregate(GenomeEquiv ~ HabitatAir ,data = I6S_qPCR_dayList[[13]], mean)

# # TEST #
# ggplot(data=I6S_qPCR_dayList[[1]], aes(x=HabitatAir, y=GenomeEquiv, fill=HabitatAir)) +
#   geom_boxplot() +
#   theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) +
#   scale_fill_manual(values =c("forest" = "goldenrod", "savanna" = "forestgreen")) +
#   geom_point(color="black", size=2, alpha=0.9) +
#   theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold")) +
#   geom_bar(aes(), stat="identity", position="fill") +
#   theme_bw() +
#   theme(legend.position="none")


################# ITS #################
### ONLY SAMPLES (no controls) ###
# Get sample data for samples used so that can subset qPCR samples to match
ITSsampDat <- as.data.frame(as.matrix(sample_data(allITSr_noAirSingsDoubs.ps)))
head(ITSsampDat)
rownames(ITSsampDat)
unique(ITSsampDat$isControl) #only has samples (no controls were retained after rarefying)
ITSsampNames <- rownames(ITSsampDat)
ITSairNames <- ITSsampNames[grepl(x=ITSsampNames, pattern="air")]

# Filter the ITS qPCR by the air samples that were retained after rarefying
ITSairNamesNumbers <- gsub(x=ITSairNames, pattern = "air_ITS_", replacement= "")
# Plate 1
ITSqPCR_Plate1_trimmed <- fung1[fung1$sampleName %in% ITSairNamesNumbers,]
# Plate 2
ITSqPCR_Plate2_trimmed <- fung2[fung2$sampleName %in% ITSairNamesNumbers,]

# Check to make sure that the column names are the same for the two ITS plates
colnames(ITSqPCR_Plate1_trimmed[1:7]) == colnames(ITSqPCR_Plate2_trimmed)
ITSqPCRdata <- rbind(ITSqPCR_Plate1_trimmed[,1:7], ITSqPCR_Plate2_trimmed)
nrow(ITSqPCRdata) #110 like number of samples

# Make a copy and re-organize so that it can be put in table in the supplementary data
ITSqPCRdata_forTable <- ITSqPCRdata
# Make sample name numeric so that it can be arranged!
ITSqPCRdata_forTable$sampleName <- as.numeric(ITSqPCRdata_forTable$sampleName)
# Round these to the nearest whole number (which was ultimtely also done in the statistical analyses!)
ITSqPCRdata_forTable$GenomeEquiv <- round(ITSqPCRdata_forTable$GenomeEquiv)
ITSqPCRdata_forTable <- ITSqPCRdata_forTable %>% 
  arrange(sampleName)  #arrange by sample name
# View(ITSqPCRdata_forTable)
# write.csv(ITSqPCRdata_forTable, file= "ITSqPCRdata_forTable.csv") #written Nov 13, 2024 and exported to own computer 
# to make supplemental figure (but note not final way of doing it since I ended up making one plot for fungal and bacterial metadata)

# Plot the non-control samples (only those that made it through filters employed 
# earlier)
#quartz()
plot(sort(ITSqPCRdata$GenomeEquiv))
plot(sort(log(ITSqPCRdata$GenomeEquiv)))

### PLOT ITS QPCR DATA BY DAY! (to match the taxonomy plots above) ###
head(ITSqPCRdata)
# This list, made earlier in this script, has all of the samples by day. Will use ITS since this has all bioaerosol samples
# that were in 16S and then some. 
# Need to get the sample names from each and use this to organize the qPCR samples 
ITS_bySamplingDay_psList
# pre-allocate a list to store all of the sample names and HabitatAir info by day
ITS_byDaySampNamesHabitat <- vector("list", length=length(ITS_bySamplingDay_psList))
names(ITS_byDaySampNamesHabitat) <- names(ITS_bySamplingDay_psList)
# Get sample names for each day
for (j in 1:length(ITS_byDaySampNamesHabitat)){
  # get the sample names (rownames extracted below) and their habitat air information
  ITS_byDaySampNamesHabitat[[j]] <- cbind(rownames(sample_data(ITS_bySamplingDay_psList[[j]])), sample_data(ITS_bySamplingDay_psList[[j]])$HabitatAir)
  # Replace the full sample names with the shorter ones that match the qPCR labels 
  ITS_byDaySampNamesHabitat[[j]] <- gsub(x=ITS_byDaySampNamesHabitat[[j]], pattern = "air_ITS_", replacement = "")
  ITS_byDaySampNamesHabitat[[j]] <- as.data.frame(as.matrix(ITS_byDaySampNamesHabitat[[j]]))
  colnames(ITS_byDaySampNamesHabitat[[j]]) <- c("sampleName", "HabitatAir")
}
# Tests to show above for loop is working as expected with one day
test1 <- cbind(rownames(sample_data(ITS_bySamplingDay_psList[[9]])), sample_data(ITS_bySamplingDay_psList[[9]])$HabitatAir)
test2 <- gsub(x=test1, pattern = "air_ITS_", replacement = "")
test2 <-as.data.frame(as.matrix(test2))
test2 ==ITS_byDaySampNamesHabitat[[9]]
colnames(test2) <- c("sampleName", "HabitatAir")
test2 == ITS_byDaySampNamesHabitat[[9]]

ITS_byDaySampNamesHabitat[[9]]

# GET qPCR DATA SEPARATED BY DAY 
# Pre-allocate a new list to hold each days qPCR info
ITS_qPCR_dayList <- vector("list", length(ITS_byDaySampNamesHabitat))
names(ITS_qPCR_dayList) <- names(ITS_byDaySampNamesHabitat)

# Loop through each element in ITS_byDaySampNamesHabitat
for (i in 1:length(ITS_byDaySampNamesHabitat)) {
  # Extract the sample names from the current element in the list
  sample_names <- ITS_byDaySampNamesHabitat[[i]]$sampleName
  
  # Subset ITSqPCRdata to only include rows where sampleName matches those in the current list element
  day_subset <- ITSqPCRdata[ITSqPCRdata$sampleName %in% sample_names, ]
  
  # Merge the subset of df1 with the current element of list1 based on sampleName
  # Using all.x = TRUE to keep all rows from the list1 data frame in case there are unmatched sampleNames
  ITS_qPCR_dayList[[i]] <- merge(ITS_byDaySampNamesHabitat[[i]], day_subset, by = "sampleName", all.x = TRUE)
}

# Spot check
ITS_qPCR_dayList[[9]] #looks correct!
ITS_byDaySampNamesHabitat[[9]] #looks correct!

# # Re-order these based on the order that I used for the earlier plots
# # Reorder merged_dfs according to the desired order
# ITS_qPCR_dayList <- ITS_qPCR_dayList[names(ITS_Order_PlotList)]
# names(ITS_qPCR_dayList) == names(ITS_Order_PlotList)

# PLOT BOXPLOTS!!!
# Make new list for plots
ITS_qPCR_PlotList <- vector("list", length=length(ITS_qPCR_dayList))
names(ITS_qPCR_dayList) <- names(ITS_qPCR_dayList)

# Can do a little in a for loop:
for (k in 1:length(ITS_qPCR_PlotList)) {
  ITS_qPCR_PlotList[[k]] <- ggplot(data=ITS_qPCR_dayList[[k]], aes(x=HabitatAir, y=GenomeEquiv, fill=HabitatAir)) + 
    geom_boxplot() +
    theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) + 
    scale_fill_manual(values = c("forest" = "forestgreen", "savanna" = "goldenrod")) + 
    geom_point(color="black", size=2, alpha=0.9) +
    theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold")) +
    geom_bar(aes(), stat="identity", position="fill") +
    theme_bw() +
    theme(legend.position="none") 
}

# PLOT THEM ALL
# quartz()
do.call(grid.arrange, c(ITS_qPCR_PlotList, ncol = 4, nrow = 4))

aggregate(GenomeEquiv ~ HabitatAir ,data = ITS_qPCR_dayList[[1]], mean)
aggregate(GenomeEquiv ~ HabitatAir ,data = ITS_qPCR_dayList[[2]], mean)
aggregate(GenomeEquiv ~ HabitatAir ,data = ITS_qPCR_dayList[[5]], mean)
aggregate(GenomeEquiv ~ HabitatAir ,data = ITS_qPCR_dayList[[8]], mean)
aggregate(GenomeEquiv ~ HabitatAir ,data = ITS_qPCR_dayList[[10]], mean)
ITS_qPCR_dayList[[13]]
aggregate(GenomeEquiv ~ HabitatAir ,data = ITS_qPCR_dayList[[13]], mean)

##### SECOND WAY OF MAKING QPCR PLOTS ######
###### 16S ######
# Box and whisker plots with all samples, forest one column and savanna the other
# y-axis should be in log SQ, so like 10^1, 10^2, 10^3 and so on 

# MAKE OBJECT THAT HAS ALL DATA IN ONE:
allI6S_qPCR <- do.call(rbind, I6S_qPCR_dayList)
aggregate(GenomeEquiv ~ HabitatAir, data = allI6S_qPCR, mean)
mean(allI6S_qPCR$GenomeEquiv)
#View(allI6S_qPCR)
allI6S_qPCR <- allI6S_qPCR %>% 
  as.data.frame() %>% 
  rownames_to_column(var= "samplingDay") 
# Remove numbers and dots that signifed numbers in dataframe 
allI6S_qPCR$samplingDay <- gsub(x= allI6S_qPCR$samplingDay, pattern = "\\.\\d+", replacement = "")
# View(allI6S_qPCR)

# Add in EU data. First get which EU was was sampled when
dayEU_df <- as.data.frame(matrix(ncol =2, nrow= length(ITS_bySamplingDay_psList)))
for (i in 1:length(ITS_bySamplingDay_psList)){
  dayEU_df[i,] <- as.data.frame(cbind(names(ITS_bySamplingDay_psList)[[i]], unique(sample_data(ITS_bySamplingDay_psList[[i]])$EU)))
}
dayEU_df

# Add each of days per EU to an objects below. Cross referenced, just to make
# super super sure, with "SamplingDetailsByDate Excel sheet"
EU_52days <- c("16-Jun-2022", "18-Jun-2022", "26-Jun-2022", "1-Jul-2022")
EU_53Sdays <- c("22-Jun-2022", "27-Jun-2022", "2-Jul-2022", "6-Jul-2022")
EU_54Sdays <- c("24-Jun-2022", "30-Jun-2022", "5-Jul-2022", "8-Jul-2022")
EU_8days <- c("20-Jun-2022", "23-Jun-2022", "29-Jun-2022", "4-Jul-2022")

# EDIT 16S DATA FRAME
allI6S_qPCR <- allI6S_qPCR %>%
  mutate(
    EU = case_when(
      samplingDay %in% EU_52days ~ "EU_52",
      samplingDay %in% EU_8days ~ "EU_8",
      samplingDay %in% EU_53Sdays  ~ "EU_53S",
      samplingDay %in% EU_54Sdays ~ "EU_54S",
      TRUE ~ NA_character_
    )
  )
# View(allI6S_qPCR)

# PLOT IT:
# Make plot:
allI6S_qPCR_byHabitat <- ggplot(data=allI6S_qPCR, aes(x=HabitatAir, y=GenomeEquiv, fill=HabitatAir)) + 
  geom_boxplot() +
  scale_y_log10(name = "Bacterial genome equivalents", breaks = c(1e+03, 1e+04),  labels = c("1e+03", "1e+04")) +  # Apply logarithmic scale
  labs(x= NULL) +
  scale_fill_manual(values = c("forest" = "forestgreen", "savanna" = "goldenrod")) + 
  geom_jitter(color="black", size=2, alpha=0.9, height = 0, width = 0.35) + 
  theme_bw() +
  theme(
    legend.position="none",
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 16),  
    axis.text.x = element_text(colour = "black", size = 14)  #x-axis tick labels
  )

# quartz()
allI6S_qPCR_byHabitat

###### ITS ######
# MAKE OBJECT THAT HAS ALL DATA IN ONE, WITH Y-AXIS ON LOG SCALE:
allITS_qPCR <- do.call(rbind, ITS_qPCR_dayList)
mean(allITS_qPCR$GenomeEquiv)

#View(allITS_qPCR)
allITS_qPCR <- allITS_qPCR %>% 
  as.data.frame() %>% 
  rownames_to_column(var= "samplingDay") 
# Remove numbers and dots that signifed numbers in dataframe 
allITS_qPCR$samplingDay <- gsub(x= allITS_qPCR$samplingDay, pattern = "\\.\\d+", replacement = "")
#View(allITS_qPCR)

# EDIT ITS DATA FRAME
allITS_qPCR <- allITS_qPCR %>%
  mutate(
    EU = case_when(
      samplingDay %in% EU_52days ~ "EU_52",
      samplingDay %in% EU_8days ~ "EU_8",
      samplingDay %in% EU_53Sdays  ~ "EU_53S",
      samplingDay %in% EU_54Sdays ~ "EU_54S",
      TRUE ~ NA_character_
    )
  )
# View(allITS_qPCR)

# PLOT IT:
allITS_qPCR_byHabitat <- ggplot(data=allITS_qPCR, aes(x=HabitatAir, y=GenomeEquiv, fill=HabitatAir)) + 
  geom_boxplot() +
  scale_y_log10(name = "Fungal genome equivalents", labels = scales::scientific) +  # Apply logarithmic scale
  labs(x= NULL) +
  scale_fill_manual(values = c("forest" = "forestgreen", "savanna" = "goldenrod")) + 
  geom_jitter(color="black", size=2, alpha=0.9, height = 0, width = 0.35) + 
  theme_bw() +
  theme(
    legend.position="none",
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 16),  
    axis.text.x = element_text(colour = "black", size = 14)  #x-axis tick labels
  )

# quartz()
allITS_qPCR_byHabitat

# PLOT BACTERIAL AND FUNGAL PLOTS TOGETHER
# quartz()
grid.arrange(allITS_qPCR_byHabitat, allI6S_qPCR_byHabitat) #saved these plots (from Quartz window) as "qPCR_allDaysBoxplots"

# Save and export 
# saveRDS(allITS_qPCR_byHabitat, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/allITS_qPCR_byHabitat_09-21-2025.rds")
# saveRDS(allI6S_qPCR_byHabitat, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/allI6S_qPCR_byHabitat_09-21-2025.rds")

################# QPCR STATISTICS #################
######### 16S #########
### 1. SET UP AND CHECKING
head(allI6S_qPCR)
unique(is.na(allI6S_qPCR$GenomeEquiv)) #no NAs
unique(is.na(allI6S_qPCR$HabitatAir)) #no NAs
# Make everything that needs to be a factor into a factor (especially important for the GLMMs below)
colnames(allI6S_qPCR)
allI6S_qPCR$HabitatAir <- as.factor(allI6S_qPCR$HabitatAir)
allI6S_qPCR$samplingDay <- as.factor(allI6S_qPCR$samplingDay)
allI6S_qPCR$EU <- as.factor(allI6S_qPCR$EU)
allI6S_qPCR$GenomeEquiv <- as.numeric(allI6S_qPCR$GenomeEquiv)

### 2. SIMPLE T-TEST/WILCOXON (NOT REPORTED IN MANUSCRIPT SINCE ENDED UP USING MORE COMPLICATED MODELS BELOW ###
# Check to see if volume is normally distributed with Shapiro-Wilk test
shapiro.test(x = allI6S_qPCR$GenomeEquiv) #W = 0.7224, p-value = 0.00000000002626 these data are NOT normal (reject null that is normal)
# Does log transforming these make them normal (this also matches the likely plots...). YES!
# Q-Q plot - looking much better and shapiro tests confirms that normal once log-transformed!
qqnorm(allI6S_qPCR$GenomeEquiv)
qqline(allI6S_qPCR$GenomeEquiv)
qqnorm(log(allI6S_qPCR$GenomeEquiv))
qqline(log(allI6S_qPCR$GenomeEquiv))
shapiro.test(x = log(allI6S_qPCR$GenomeEquiv)) # W = 0.98876, p-value = 0.6853
# Given that these basically normal, prepare Welch's t-test on log-transformed values
t.test(log(GenomeEquiv) ~ HabitatAir, data = allI6S_qPCR) #t = -0.29295, df = 81.886, p-value = 0.7703
# But for good measure, non-parametric test (Wilcoxon rank sum test/Mann-Whitney test)
qPCR_I6S_Wilcoxon <- wilcox.test(log(GenomeEquiv) ~ HabitatAir, data = allI6S_qPCR)
qPCR_I6S_Wilcoxon # W = 840, p-value = 0.7646 Shows that there is no difference among groups

### 3. MODELS! ###
# GENERALIZED LINEAR MODEL TO TEST FOR THE IMPORTANCE OF HABITAT WITH DAY AND EU AS RANDOM EFFECTS ###
# using guide here: https://entnemdept.ufl.edu/Hahn/generalized-linear-mixed-models.html
# First, round all Genome Equivalents to the nearest whole number, so I can use a Poisson distribution
allI6S_qPCR$GenomeEquiv_rounded <- round(allI6S_qPCR$GenomeEquiv)
# Get natural log-transformed data, since that was better in first models above
allI6S_qPCR$NatLogGenomeEquivRound <- log(allI6S_qPCR$GenomeEquiv_rounded)

# CONSTRUCT SEVERAL POSSIBLE MODELS:
colnames(allI6S_qPCR)
#### 1. Use NON log-transformed data, just as a first pass. First pass is bad!
I6SqPCRgaussianBasic_notLog <- glmmTMB(GenomeEquiv ~ HabitatAir, data=allI6S_qPCR, family="gaussian") 
I6SqPCRgaussianBasic_notLog
summary(I6SqPCRgaussianBasic_notLog)
# AIC      BIC   logLik deviance df.resid 
# 1779.8   1787.0   -886.9   1773.8       81 
I6SqPCRgaussianBasic_notLog_resid <- resid(I6SqPCRgaussianBasic_notLog)
qqnorm(I6SqPCRgaussianBasic_notLog_resid)
qqline(I6SqPCRgaussianBasic_notLog_resid) #ooh this looks BAD

#### 2. Gaussian GLM w/o random effects just to check it out. But will let structure of the collected
# data inform the model
# Will use log-transformed data, since that made data closer to normal in simple linear models above
I6SqPCRgaussianBasic <- glmmTMB(NatLogGenomeEquivRound ~ HabitatAir, data=allI6S_qPCR, family="gaussian") 
I6SqPCRgaussianBasic
# Data: allI6S_qPCR
# AIC       BIC    logLik -2*log(L)  df.resid 
# 272.8113  280.1037 -133.4056  266.8113        81 
# Dispersion estimate for gaussian family (sigma^2):  1.4 
# CHECKS 
# 1. Check for overdispersion
check_overdispersion(I6SqPCRgaussianBasic)
# # Overdispersion test
# 
# dispersion ratio = 1.011
# p-value = 0.904
# 
# No overdispersion detected.
summary(I6SqPCRgaussianBasic)
# 2. Confirm that residuals are normal (they are)
I6SqPCRgaussianBasic_Log_resid <- resid(I6SqPCRgaussianBasic)
qqnorm(I6SqPCRgaussianBasic_Log_resid)
qqline(I6SqPCRgaussianBasic_Log_resid)
shapiro.test(I6SqPCRgaussianBasic_Log_resid) #W = 0.98992, p-value = 0.7644

#### 3. Gaussian GLM with random effects. Will also use log-transformed data
I6SqPCRgaussianRandom_DayEU <- glmmTMB(NatLogGenomeEquivRound ~ HabitatAir + (1|EU/samplingDay), data=allI6S_qPCR, family = "gaussian")
I6SqPCRgaussianRandom_DayEU
summary(I6SqPCRgaussianRandom_DayEU)
# AIC      BIC   logLik deviance df.resid 
# 272.1    284.3   -131.1    262.1       79 
# Conditional model:
#   Estimate Std. Error z value            Pr(>|z|)    
# (Intercept)        8.28276    0.20779   39.86 <0.0000000000000002 ***
#   HabitatAirsavanna  0.03342    0.23762    0.14               0.888   
# Because original data are log-transformed,
exp(8.28277) #3955.135 is the estimate for the intercept
exp(0.03342) #1.033985 is the estimate for the HabitatAirSavanna
# CHECKS 
# 1. Check for overdispersion
check_overdispersion(I6SqPCRgaussianRandom_DayEU)
# # Overdispersion test -- none detected!
# dispersion ratio = 1.014
# p-value = 0.88

# 2. Confirm that residuals are normal
I6SqPCRgaussianRandom_DayEU_resid <- resid(I6SqPCRgaussianRandom_DayEU)
qqnorm(I6SqPCRgaussianRandom_DayEU_resid)
qqline(I6SqPCRgaussianRandom_DayEU_resid)
shapiro.test(I6SqPCRgaussianRandom_DayEU_resid) #data:  I6SqPCRgaussianRandom_DayEU_resid
# W = 0.99408, p-value = 0.9707.. THEY ARE NORMAL!

#### 4. # Basic Poisson GLM Used rounded genome equivalent data (Poisson requires integers)
I6SqPCR_poissonBasic <- glmmTMB(GenomeEquiv_rounded ~ HabitatAir, data=allI6S_qPCR, family="poisson") #poisson glm
I6SqPCR_poissonBasic
summary(I6SqPCR_poissonBasic)

# Data: allI6S_qPCR
# AIC       BIC    logLik  df.resid 
# 699273.0  699277.8 -349634.5        82 
check_overdispersion(I6SqPCR_poissonBasic)
# # Overdispersion test -- OVERDISPERSED!!!

# dispersion ratio =  11039.211
# Pearson's Chi-Squared = 905215.296
#                 p-value =    < 0.001

#### 5. Poisson GLM with both possible random effects. Used rounded genome equivalent data (Poisson requires integers)
I6SqPCR_poissonRandom_DayEU <- glmmTMB(GenomeEquiv_rounded ~ HabitatAir + (1|EU/samplingDay), data=allI6S_qPCR, family="poisson") #poisson glm
I6SqPCR_poissonRandom_DayEU
summary(I6SqPCR_poissonRandom_DayEU)
# Conditional model:
#   Estimate Std. Error z value            Pr(>|z|)    
# (Intercept)        8.802988   0.192724   45.68 <0.0000000000000002 ***
#   HabitatAirsavanna -0.306270   0.002601 -117.73 <0.0000000000000002 ***

check_overdispersion(I6SqPCR_poissonRandom_DayEU)
# # Overdispersion test -- oh no, OVERDISPERSED!
# 
# dispersion ratio =   5675.455
# Pearson's Chi-Squared = 454036.372
#                 p-value =    < 0.001

# 6. Negative binomial model with all predictors
I6SqPCR_NB_Random_DayEU <- glmmTMB(GenomeEquiv_rounded ~ HabitatAir + (1|EU/samplingDay), data=allI6S_qPCR, family="nbinom2")
Anova(I6SqPCR_NB_Random_DayEU)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: GenomeEquiv_rounded
# Chisq Df Pr(>Chisq)
# HabitatAir 1.3937  1     0.2378
dim(allI6S_qPCR)

# Set seed since these are simulations!-- ALL LOOK GOOD!
set.seed(1121)
plot(simulateResiduals(I6SqPCR_NB_Random_DayEU,  n=1000))
set.seed(1121)
hist(simulateResiduals(I6SqPCR_NB_Random_DayEU, n=1000))

######### ITS #########
head(allITS_qPCR)
dim(allITS_qPCR)
unique(is.na(allITS_qPCR$GenomeEquiv)) #no NAs
unique(is.na(allITS_qPCR$HabitatAir)) #no NAs
# Make everything that needs to be a factor into a factor (especially important for the GLMMs below)
colnames(allITS_qPCR)
allITS_qPCR$HabitatAir <- as.factor(allITS_qPCR$HabitatAir)
allITS_qPCR$samplingDay <- as.factor(allITS_qPCR$samplingDay)
allITS_qPCR$EU <- as.factor(allITS_qPCR$EU)
allITS_qPCR$GenomeEquiv <- as.numeric(allITS_qPCR$GenomeEquiv)

### SIMPLE T-TEST/WILCOXON ###
# Check to see if volume is normally distributed with Shapiro-Wilk test
shapiro.test(x = allITS_qPCR$GenomeEquiv) #W = 0.60261, p-value = 0.0000000000000008461 these data are NOT normal (reject null that is normal)
# Does log transforming these make them normal (this also matches the likely plots...). 
# Q-Q plot - looking much better and shapiro tests confirms that normal once log-transformed!
qqnorm(log(allITS_qPCR$GenomeEquiv))
qqline(log(allITS_qPCR$GenomeEquiv))
shapiro.test(x = log(allITS_qPCR$GenomeEquiv)) #W = 0.98351, p-value = 0.1934
# Given that these basically normal, prepare Welch's t-test on log-transformed values
t.test(log(GenomeEquiv) ~ HabitatAir, data = allITS_qPCR) #t = -0.88981, df = 103.13, p-value = 0.3756
# But for good measure, non-parametric test (Wilcoxon rank sum test/Mann-Whitney test)
allITS_qPCR_Wilcoxon <- wilcox.test(log(GenomeEquiv) ~ HabitatAir, data = allITS_qPCR)
allITS_qPCR_Wilcoxon #W = 1340, p-value = 0.3256 Shows that there is no difference among groups

### MODELS! ###
# GENERALIZED LINEAR MODEL TO TEST FOR THE IMPORTANCE OF HABITAT WITH DAY AND EU AS RANDOM EFFECTS ###
# using guide here: https://entnemdept.ufl.edu/Hahn/generalized-linear-mixed-models.html
# First, round all Genome Equivalents to the nearest whole number, so I can use a Poisson distribution
allITS_qPCR$GenomeEquiv_rounded <- round(allITS_qPCR$GenomeEquiv)
# Get natural log-transformed data, since that was better in first models above
allITS_qPCR$NatLogGenomeEquivRound <- log(allITS_qPCR$GenomeEquiv_rounded)

# CONSTRUCT SEVERAL POSSIBLE MODELS:
colnames(allITS_qPCR)
#### 1. Use NON log-transformed data, just as a first pass
ITSqPCRgaussianBasic_notLog <- glmmTMB(GenomeEquiv ~ HabitatAir, data=allITS_qPCR, family="gaussian") 
ITSqPCRgaussianBasic_notLog
summary(ITSqPCRgaussianBasic_notLog)
# AIC      BIC   logLik deviance df.resid 
# 3052.6   3060.7  -1523.3   3046.6      107 
# Check out residuals
ITSqPCRgaussianBasic_notLog_resid <- resid(ITSqPCRgaussianBasic_notLog)
qqnorm(ITSqPCRgaussianBasic_notLog_resid)
qqline(ITSqPCRgaussianBasic_notLog_resid) #ooh this looks not great, as expected

#### 2. Gaussian GLM w/o random effects just to check it out. But will let structure of the collected
# data inform the model
# Will use log-transformed data, since that made data closer to normal in simple linear models above
ITSqPCRgaussianBasic <- glmmTMB(NatLogGenomeEquivRound ~ HabitatAir, data=allITS_qPCR, family="gaussian") 
ITSqPCRgaussianBasic
# Data: allITS_qPCR
# AIC       BIC    logLik -2*log(L)  df.resid 
# 356.9538  365.0553 -175.4769  350.9538       107 
# CHECKS 
# 1. Check for overdispersion
check_overdispersion(ITSqPCRgaussianBasic)
# # Overdispersion test
# 
# dispersion ratio = 1.005
# p-value = 0.912
# 
# No overdispersion detected.
summary(ITSqPCRgaussianBasic)
# 2. Confirm that residuals are normal
ITSqPCRgaussianBasic_Log_resid <- resid(ITSqPCRgaussianBasic)
qqnorm(ITSqPCRgaussianBasic_Log_resid)
qqline(ITSqPCRgaussianBasic_Log_resid)
shapiro.test(ITSqPCRgaussianBasic_Log_resid) #W = 0.98198, p-value = 0.1437

#### 3. Gaussian GLM with random effects. Will also use log-transformed data
ITSqPCRgaussianRandom_DayEU <- glmmTMB(NatLogGenomeEquivRound ~ HabitatAir + (1|EU/samplingDay), data=allITS_qPCR, family = "gaussian")
ITSqPCRgaussianRandom_DayEU
summary(ITSqPCRgaussianRandom_DayEU)

# AIC      BIC   logLik deviance df.resid 
# 326.9    340.4   -158.4    316.9      105 
# 
# Random effects:
#   
# Conditional model:
#   Groups         Name        Variance      Std.Dev.  
# samplingDay:EU (Intercept) 0.61891575681 0.78671199
# EU             (Intercept) 0.00000000177 0.00004207
# Residual                   0.79910723376 0.89392798
# Number of obs: 110, groups:  samplingDay:EU, 16; EU, 4
# 
# Dispersion estimate for gaussian family (sigma^2): 0.799 
# 
# Conditional model:
#   Estimate Std. Error z value            Pr(>|z|)    
# (Intercept)        11.1248     0.2290   48.59 <0.0000000000000002 ***
#   HabitatAirsavanna   0.2226     0.1722    1.29               0.196    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Because original data are log-transformed,
exp(11.1248) #67832.72 is the estimate for the intercept
exp(0.2226) #1.249321 is the estimate for the HabitatAirSavanna
# CHECKS 
# 1. Check for overdispersion
check_overdispersion(ITSqPCRgaussianRandom_DayEU)
# Overdispersion test
# dispersion ratio = 1.049
# p-value = 0.712
# 
# No overdispersion detected.
# 2. Confirm that residuals are normal-- but they are NOT!
ITSqPCRgaussianRandom_DayEU_resid <- resid(ITSqPCRgaussianRandom_DayEU)
qqnorm(ITSqPCRgaussianRandom_DayEU_resid)
qqline(ITSqPCRgaussianRandom_DayEU_resid)
shapiro.test(ITSqPCRgaussianRandom_DayEU_resid) #data:  ITSqPCRgaussianRandom_DayEU_resid

# ... OOOOO these are NOT normal, so need to try the Poisson
# Shapiro-Wilk normality test
# 
# data:  ITSqPCRgaussianRandom_DayEU_resid
# W = 0.96991, p-value = 0.01366


#### 4. # Basic Poisson GLM Used rounded genome equivalent data (Poisson requires integers)
ITSqPCR_poissonBasic <- glmmTMB(GenomeEquiv_rounded ~ HabitatAir, data=allITS_qPCR, family="poisson") #poisson glm
ITSqPCR_poissonBasic
summary(ITSqPCR_poissonBasic)

# Data: allI6S_qPCR
# AIC      BIC   logLik deviance df.resid 
# 18186890 18186896 -9093443 18186886      108 

# Conditional model:
#   Estimate Std. Error z value            Pr(>|z|)    
# (Intercept)       11.7597956  0.0003639   32317 <0.0000000000000002 ***
#   HabitatAirsavanna  0.2843356  0.0004977     571 <0.0000000000000002 ***

#### 5. Poisson GLM with both possible random effects. Used rounded genome equivalent data (Poisson requires integers)
ITSqPCR_poissonRandom_DayEU <- glmmTMB(GenomeEquiv_rounded ~ HabitatAir + (1|EU/samplingDay), data=allITS_qPCR, family="poisson") #poisson glm
ITSqPCR_poissonRandom_DayEU
summary(ITSqPCR_poissonRandom_DayEU)
check_overdispersion(ITSqPCR_poissonRandom_DayEU) #OVERDISPERSED!!!!
# # Overdispersion test
# 
# dispersion ratio =   65057.393
# Pearson's Chi-Squared = 6896083.646
#                 p-value =     < 0.001


# Formula:          GenomeEquiv_rounded ~ HabitatAir + (1 | samplingDay) + (1 | EU)
# Data: allI6S_qPCR
# AIC       BIC    logLik  df.resid 
# 437684.4  437694.1 -218838.2        80 

# 6. Poisson GLM with only EU random effects. Used rounded genome equivalent data (Poisson requires integers)
ITSqPCR_poissonRandom_EU <- glmmTMB(GenomeEquiv_rounded ~ HabitatAir + (1|EU), data=allITS_qPCR, family="poisson") #poisson glm
ITSqPCR_poissonRandom_EU
check_overdispersion(ITSqPCR_poissonRandom_EU) #OVERDISPERSED!!!!
# Overdispersion test
# 
# dispersion ratio =   203586.550
# Pearson's Chi-Squared = 21783760.872
#                 p-value =      < 0.001

# 7. Negative binomial model with all predictors
ITSqPCR_NB_Random_DayEU <- glmmTMB(GenomeEquiv_rounded ~ HabitatAir + (1|EU/samplingDay), data=allITS_qPCR, family="nbinom2")
Anova(ITSqPCR_NB_Random_DayEU)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: GenomeEquiv_rounded
# Chisq Df Pr(>Chisq)
# HabitatAir 1.484  1     0.2232
dim(allITS_qPCR) #110 samples

# Set seed since these are simulations! And all model assumptions look good!
set.seed(1121)
plot(simulateResiduals(ITSqPCR_NB_Random_DayEU,  n=1000))
set.seed(1121)
hist(simulateResiduals(ITSqPCR_NB_Random_DayEU, n=1000))

# Compare AICs -- as expected, NB sooooo much lower
AIC(ITSqPCR_poissonRandom_EU,ITSqPCR_NB_Random_DayEU)
# df          AIC
# ITSqPCR_poissonRandom_EU  3 16793158.172
# ITSqPCR_NB_Random_DayEU   5     2796.129
