# ANCOM_traits_I6SandITS_Sept.R
# Old version was called ANCOM_traits_I6SandITS.R. Old, old version was called DiffAbundTests_I6SandITS.R

# This script follows ANCOMBC_robust_foliarAir_pt0.05 (own computer) and bacteriaSporesAndPigment.R (run on server, but saved in both locations)

# Differential Abundance between savanna and forest air
# PROCESSING ANCOM-BC RESULTS FOR BIOAEROSOLS

# Script description: This script follows I6S_sourceTracking.R (and fungal equilavent) and ancom_Sept.R on own computer
# It also requires "ANCOMspore.RData" from bacterialSporeFormers.R

# Objects saved in this script: UPDATE THIS

#######################################
# I.  SCRIPT SET UP
#######################################
##### 1. LOAD LIBRARIES #####
library(tidyverse)
library(indicspecies)
library(vegan)
library(phyloseq)
library(stringr)
library(Polychrome) #for color palette
library(robustbase) #for colMeans() function
library(grDevices)
library(rcartocolor)
library(gridExtra)

# 1. READ IN PHYLOSEQ DATA (NON-RAREFIED DATA):
# i. BIOAEROSOLS ONLY
# 1. BACTERIA: (From full16S_EDArarefied_part2_Sept.R, has singletons and doubletons and samples with < 5500 reads removed)
I6Sall_d_5.5K_justAir.ps <- readRDS(file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6Sall_d_5.5K_justAir_phyloseq.rds")
# 2. FUNGI:(From fullITS_EDArarefied_part2_Sept.R, has singletons and doubletons and samples with < 8500 reads removed)
airITS_8.5K.ps <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airITS_8.5K_phyloseq.rds")
# ii. ALL data
# 1.BACTERIALoad unrarefied data for all samples (air singletons and doubletons removed and samples lower than 5.5K reads also removed)
I6Sall_5.5Kfiltered.ps <-  readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airI6S_5.5Kfiltered_phyloseq.rds")
# 2. FUNGI: Load unrarefied data for all samples (air singletons and doubletons removed and samples lower than 8.5K reads alos removed)
ITSall_8.5Kfiltered.ps <-  readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airITS_8.5Kfiltered_phyloseq.rds")

# 2. READ IN ANCOM RESULTS (FOR BIOAEROSOLS/FOLIAR SURFACES ONLY, made in ANCOMBC_robust_foliarAir_pt0.05.R)
I6S_ANCOMall_.05pct_df <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_ANCOMall_.05pct_df.rds")
ITS_ANCOMall_.05pct_df <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_ANCOMall_.05pct_df.rds")

# # 3. READ IN SPORULATION AND PIGMENTATION INFORMATION (BOTH made in bacterialSporeFormers_Sept.R on server)
ANCOMspore <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airLeafANCOM_withSpore_Oct13_2025.rds")
ANCOMgenPig <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ANCOM_ijsemGenusPig_10-14-2025.rds")

########################################################################
#   II. PHYLOSEQ OBJECT DOUBLE CHECKING 
########################################################################
I6Sall_5.5Kfiltered.ps #32744 taxa and 299 samples
sort(colSums(as.data.frame(as.matrix(otu_table(I6Sall_5.5Kfiltered.ps))))) #varied, as they should be 
ITSall_8.5Kfiltered.ps #13778 taxa and 324 samples
sort(colSums(as.data.frame(as.matrix(otu_table(ITSall_8.5Kfiltered.ps))))) #varied, as they should be 

########### i. BACTERIA ###########
# subset no soil
unique(sample_data(I6Sall_5.5Kfiltered.ps)$sampleType)
I6SairPhylloOnly.ps <- subset_samples(I6Sall_5.5Kfiltered.ps, sampleType != "soil")
unique(sample_data(I6SairPhylloOnly.ps)$sampleType)
# Remove ASVs that do not have taxa present
I6SairPhylloOnly.psZeros <- which(rowSums(otu_table(I6SairPhylloOnly.ps))==0)
unique(rowSums(otu_table(I6SairPhylloOnly.ps)[I6SairPhylloOnly.psZeros,])) #all zeros
I6SairPhylloOnly.psZerosNames <- names(I6SairPhylloOnly.psZeros)
I6SairPhylloOnly.psASVsToKeep <- setdiff(rownames(otu_table(I6SairPhylloOnly.ps)), I6SairPhylloOnly.psZerosNames)
length(I6SairPhylloOnly.psASVsToKeep) #10,540
I6SairPhylloOnly.ps <- prune_taxa(taxa= I6SairPhylloOnly.psASVsToKeep, x=I6SairPhylloOnly.ps) #9603 taxa and 142 samples
I6SairPhylloOnly.ps #10,540 taxa and 142 samples
which(rowSums(otu_table(I6SairPhylloOnly.ps)) == 0) #none, as expected

########### ii. FUNGI ###########
ITSall_8.5Kfiltered.ps
unique(sample_data(ITSall_8.5Kfiltered.ps)$sampleType)
ITSairPhylloOnly.ps <- subset_samples(ITSall_8.5Kfiltered.ps, sampleType != "soil")
unique(sample_data(ITSairPhylloOnly.ps)$sampleType)

# remove ASVs that do not have taxa present
ITSairPhylloOnly.psZeros <- which(rowSums(otu_table(ITSairPhylloOnly.ps))==0)
unique(rowSums(otu_table(ITSairPhylloOnly.ps)[ITSairPhylloOnly.psZeros,])) #all zeros
ITSairPhylloOnly.psZerosNames <- names(ITSairPhylloOnly.psZeros)
ITSairPhylloOnly.psASVsToKeep <- setdiff(rownames(otu_table(ITSairPhylloOnly.ps)), ITSairPhylloOnly.psZerosNames)
length(ITSairPhylloOnly.psASVsToKeep) #9,012
ITSairPhylloOnly.ps <- prune_taxa(taxa= ITSairPhylloOnly.psASVsToKeep, x=ITSairPhylloOnly.ps) #7,951 taxa and 142 samples
ITSairPhylloOnly.ps #9,012 taxa and 169 samples
which(rowSums(otu_table(ITSairPhylloOnly.ps)) == 0) #none, as expected

########### i. BACTERIA ###########
I6SairPhylloOnly.ps_ASVs <- t(as.data.frame(as.matrix(otu_table(I6SairPhylloOnly.ps))))
dim(I6SairPhylloOnly.ps_ASVs) #ASVs are columms and samples are rows
# Get sample type "cluster"
I6SairPhylloOnly.ps_meta <- as.data.frame(as.matrix(sample_data(I6SairPhylloOnly.ps)))

########### ii. FUNGI ###########
ITSairPhylloOnly.ps_ASVs <- t(as.data.frame(as.matrix(otu_table(ITSairPhylloOnly.ps))))
dim(ITSairPhylloOnly.ps_ASVs) #ASVs are columms (9,012) and samples are rows
# Get metadata
ITSairPhylloOnly.ps_meta <- as.data.frame(as.matrix(sample_data(ITSairPhylloOnly.ps)))

########################################################################
#   III. BACTERIA: EXPLORE ANCOM RESULTS (BIOAEROSOLS AND FOLIAR SURFACE) 
########################################################################
########### i. BACTERIA ###########
length(which(I6S_ANCOMall_.05pct_df$ANCOMcat == "bioaerosol")) #94
length(which(I6S_ANCOMall_.05pct_df$ANCOMcat == "foliar surface")) #126

##### FACETED DOT PLOT ####
#### TRY FIRST WITH FAMILY ####
# As calculated below, shows only 55/184 bioaerosol ASVs and 91/126 foliar surface ASVs, 
# so perhaps not the most representative.
# Get number of ASVs in each family
# Try at family level, as Noah suggested
numberFamsANCOM_I6S <- I6S_ANCOMall_.05pct_df %>% 
  group_by(ANCOMcat, Family) %>% 
  summarise(ASVsInFamByHab = n_distinct(ASV_name), .groups = "drop") %>%  #.groups = "drop" no more grouping
  complete(ANCOMcat, Family, fill = list(ASVsInFamByHab = 0)) #add in zero if not found in that ANCOM cat
# View(numberFamsANCOM_I6S)
table(I6S_ANCOMall_.05pct_df$Family, I6S_ANCOMall_.05pct_df$ANCOMcat) #matches calculations above!
print(arrange(numberFamsANCOM_I6S, -ASVsInFamByHab), n= nrow(numberFamsANCOM_I6S))
length(I6S_ANCOMall_.05pct_df$ASV_name) #220
# Get names of families that have at least 3 representatives in one family or the other
fams3_I6S <- numberFamsANCOM_I6S$Family[which(numberFamsANCOM_I6S$ASVsInFamByHab >=3)]
length(fams3_I6S)

# Looks good!
# View(numberFamsANCOM_I6S[numberFamsANCOM_I6S$Family %in% fams3_I6S,])
I6S_3fam_df <- numberFamsANCOM_I6S[numberFamsANCOM_I6S$Family %in% fams3_I6S,]
I6S_3fam_df <- I6S_3fam_df[-which(I6S_3fam_df$Family == "NA"),] #drop NA families
dim(I6S_3fam_df)

# Add phylum info to this plot by merging with subset of I6S_ANCOMall_.05pct_df that has phylum info! DO distinct because otherwise subset is v long and get many-to-many relationship
famPhylaI6S <- distinct(I6S_ANCOMall_.05pct_df[,colnames(I6S_ANCOMall_.05pct_df) %in% c("Family", "Phylum")])
I6S_3fam_df <- left_join(I6S_3fam_df, famPhylaI6S, by = "Family")
# View(I6S_3fam_df) # looks good!
sum(I6S_3fam_df$ASVsInFamByHab) #shows 146/(184+126 = 310 diff abund ASVs)
sum(I6S_3fam_df$ASVsInFamByHab[which(I6S_3fam_df$ANCOMcat == "foliar surface")]) #91/126 foliar surface ASVs
sum(I6S_3fam_df$ASVsInFamByHab[which(I6S_3fam_df$ANCOMcat == "bioaerosol")]) #55/184 biaoerosol ASVs

# Re-order families (making it a factor) based on alphabetical in phylum than in family
I6S_3fam_df <- I6S_3fam_df %>%
  arrange(desc(Phylum), desc(Family)) %>%  #alphabetic by Phylum, then Family, descending so that it goes up y-axis
  mutate(Family = factor(Family, levels = unique(Family))) 
# View(I6S_3fam_df) #Looks good!
I6S_3fam_df$Family
unique(I6S_3fam_df$Family)
levels(I6S_3fam_df$Family)

# Now remove those with 0 so that there are no dots for zero
I6S_3fam_df_no0s <- I6S_3fam_df[-which(I6S_3fam_df$ASVsInFamByHab == 0),]
# View(I6S_3fam_df_no0s)
levels(I6S_3fam_df_no0s$Family) #great kept levels!
# Make this numeric to squeeze left and right dots closer together
I6S_3fam_df_no0s$xpos <- as.numeric(factor(I6S_3fam_df_no0s$ANCOMcat))

# Make the plot!
colnames(I6S_3fam_df_no0s)
I6S_3famANCOM_bubPlot <- ggplot(I6S_3fam_df_no0s,
                               aes(x = xpos, y = Family, size = ASVsInFamByHab, fill = ANCOMcat)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_size(range = c(1, 6)) +   
  facet_grid(Family ~ ., scales = "free_y", space = "free_y", switch = "y", as.table = FALSE) + #as.table = FALSE orders by factor
  scale_x_continuous(
    breaks = 1:2,  # positions of the categories
    labels = levels(factor(I6S_3fam_df$ANCOMcat)),
    expand = expansion(mult = c(0.0, 0.0)),
    limits = c(0.2, 3)  # pull the two categories closer together
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "white", fill = NA, linewidth = 0.5),
    strip.placement = "outside",                 # put strips on the outer side (left)
    strip.text.y.right = element_text(angle = 0, hjust = 1, size =9, color= "black"),  # control orientation
    strip.background = element_blank(),
    panel.spacing.y = unit(0.01, "lines"),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size =5, color= "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    # Legend layout (stack the two guides)
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.box.just = "left",
    axis.text = element_text(color = "black"),
    legend.title = element_text(size = 8, color= "black"),
    legend.text  = element_text(size = 8, color= "black"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.1, "lines"),
    # Center titles above legend
    legend.title.position = "top",
    legend.title.align = 0.5,
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Number of ASVs in Family"
  ) +
  # Change bubble colors
  scale_fill_manual(
    values = c("cornflowerblue", "forestgreen") #color bioaerosol bubbles blue and foliar surface green
  ) +
  guides(fill = "none") #+ #remove fill legend for ANCOMcat
  # # SIZE legend (bottom)
  # scale_size(
  #   range = c(1, 10),
  #   guide = guide_legend(
  #     order = 1,
  #     title.position = "top",
  #     title.hjust = 0.5,
  #     override.aes = list(
  #       shape  = 21,
  #       stroke = 0.3,
  #       colour = "black"
  #     )
  #   )
  # )
#quartz(width = 6, height = 7)
I6S_3famANCOM_bubPlot

I6S_3famANCOM_bubPlot #
# NOTE: LEGEND WILL NOT LEFT JUSTIFY NOR CAN I REMOVE Y-AXIS TICKS FROM THE RIGHT SIDE, SO I WILL
# MANUALLY DO THESE THINGS IN POWERPOINT


#### CLASS ####
# (Benefits versus Class: shows more of the diff abund. ASVs and matches heat maps)
# 1. Get number of ASVs in each Class
numberClassesANCOM_I6S <- I6S_ANCOMall_.05pct_df %>% 
  group_by(ANCOMcat, Class) %>% 
  summarise(ASVsInClassByHab = n_distinct(ASV_name), .groups = "drop") %>%  #.groups = "drop" no more grouping
  complete(ANCOMcat, Class, fill = list(ASVsInClassByHab = 0)) #add in zero if not found in that ANCOM cat
# View(numberClassesANCOM_I6S)
# Check calculations
table(I6S_ANCOMall_.05pct_df$Class, I6S_ANCOMall_.05pct_df$ANCOMcat) #these match the calculations above!
unique(I6S_ANCOMall_.05pct_df$Class) #only 17, which maybe could fit except some are very ASV-rich so bubbles likely too big
sum(numberClassesANCOM_I6S$ASVsInClassByHab) #220
sum(numberClassesANCOM_I6S$ASVsInClassByHab) == length(unique(I6S_ANCOMall_.05pct_df$ASV_name)) #220

allnumberClassesANCOM_I6S <- numberClassesANCOM_I6S[-which(numberClassesANCOM_I6S$Class == "NA"),] #drop NA classes
dim(allnumberClassesANCOM_I6S)
sum(allnumberClassesANCOM_I6S$ASVsInClassByHab) #shows 215/220
sum(allnumberClassesANCOM_I6S$ASVsInClassByHab[which(allnumberClassesANCOM_I6S$ANCOMcat == "foliar surface")]) #123/126 foliar surface ASVs
sum(allnumberClassesANCOM_I6S$ASVsInClassByHab[which(allnumberClassesANCOM_I6S$ANCOMcat == "bioaerosol")]) #92/94 biaoerosol ASVs

# Add phylum info to this plot by merging with subset of I6S_ANCOMall_.05pct_df that has phylum info! DO distinct because otherwise subset is v long and get many-to-many relationship
classPhylaI6S <- distinct(I6S_ANCOMall_.05pct_df[,colnames(I6S_ANCOMall_.05pct_df) %in% c("Class", "Phylum")])
allnumberClassesANCOM_I6S <- left_join(allnumberClassesANCOM_I6S, classPhylaI6S, by = "Class")
# View(allnumberClassesANCOM_I6S) # looks good!

# Re-order classes (making it a factor) based on alphabetical in phylum than in Class
allnumberClassesANCOM_I6S <- allnumberClassesANCOM_I6S %>%
  arrange(desc(Phylum), desc(Class)) %>%  #alphabetic by Phylum, then Class, descending so that it goes up y-axis
  mutate(Class = factor(Class, levels = unique(Class))) 
# View(allnumberClassesANCOM_I6S) #Looks good!
allnumberClassesANCOM_I6S$Class
unique(allnumberClassesANCOM_I6S$Class)
levels(allnumberClassesANCOM_I6S$Class)

# Now remove those with 0 so that there are no dots for zero (had to have these zeros first to avoid arranging by ANCOMcat)
allnumberClassesANCOM_I6S_no0s <- allnumberClassesANCOM_I6S[-which(allnumberClassesANCOM_I6S$ASVsInClassByHab == 0),]
# View(allnumberClassesANCOM_I6S_no0s)
levels(allnumberClassesANCOM_I6S_no0s$Class) #great kept levels!
# Make this numeric to squeeze left and right dots closer together
allnumberClassesANCOM_I6S_no0s$xpos <- as.numeric(factor(allnumberClassesANCOM_I6S_no0s$ANCOMcat))

# Make the plot!
colnames(allnumberClassesANCOM_I6S_no0s)
I6S_allClassANCOM_bubPlot <- ggplot(allnumberClassesANCOM_I6S_no0s,
                                  aes(x = xpos, y = Class, size = ASVsInClassByHab, fill = ANCOMcat)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_size(range = c(1, 6)) +   
  facet_grid(Class ~ ., scales = "free_y", space = "free_y", switch = "y", as.table = FALSE) + #as.table = FALSE orders by factor
  scale_x_continuous(
    breaks = 1:2,  # positions of the categories
    labels = levels(factor(allnumberClassesANCOM_I6S$ANCOMcat)),
    expand = expansion(mult = c(0.0, 0.0)),
    limits = c(0.2, 3)  # pull the two categories closer together
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "white", fill = NA, linewidth = 0.5),
    strip.placement = "outside",                 # put strips on the outer side (left)
    strip.text.y.right = element_text(angle = 0, hjust = 1, size =9, color= "black"),  # control orientation
    strip.background = element_blank(),
    panel.spacing.y = unit(0.01, "lines"),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size =5, color= "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    # Legend layout (stack the two guides)
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.box.just = "left",
    axis.text = element_text(color = "black"),
    legend.title = element_text(size = 8, color= "black"),
    legend.text  = element_text(size = 8, color= "black"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.1, "lines"),
    # Center titles above legend
    legend.title.position = "top",
    legend.title.align = 0.5,
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Number of ASVs in Class"
  ) +
  # Change bubble colors
  scale_fill_manual(
    values = c("cornflowerblue", "forestgreen") #color bioaerosol bubbles blue and foliar surface green
  ) + 
  guides(fill = "none") #remove fill legend for ANCOMcat
#quartz(width = 6, height = 7)
I6S_allClassANCOM_bubPlot
# NOTE: LEGEND WILL NOT LEFT JUSTIFY NOR CAN I REMOVE 2nd set of FACET LABELS ALONG Y-AXIS, SO I WILL
# MANUALLY DO THESE THINGS IN POWERPOINT when I add phyla
# saveRDS(I6S_allClassANCOM_bubPlot, file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_allClassANCOM_bubPlot_10-16-2025.rds")

# WHAT ARE TOP CLASSES FOR BIOAEROSOL AND FOLIAR SURFACE ASVS? (ADDED TO MANUSCRIPT)
head(I6S_ANCOMall_.05pct_df)
countsDistinctClassI6SANCOM <- I6S_ANCOMall_.05pct_df %>%
  group_by(ANCOMcat, Class) %>% 
  summarise(nASVs = n_distinct(ASV_name), .groups = "drop") %>% #number of distinct ASVs per group then drop group
  group_by(ANCOMcat) %>%
  mutate(totalInANCOMcat = sum(nASVs),
         prop = nASVs / totalInANCOMcat,
         pct  = 100 * prop) %>% #pct is percentage
  ungroup()
# View(countsDistinctClassI6SANCOM)
# Do these add up to one? YES!
sum(countsDistinctClassI6SANCOM$prop[which(countsDistinctClassI6SANCOM$ANCOMcat == "bioaerosol")])
sum(countsDistinctClassI6SANCOM$prop[which(countsDistinctClassI6SANCOM$ANCOMcat == "foliar surface")])

# A check of the code above
AlphaProFS_check <- length(intersect(which(I6S_ANCOMall_.05pct_df$Class == "Alphaproteobacteria"), 
                                 which(I6S_ANCOMall_.05pct_df$ANCOMcat == "foliar surface")))/
  length(which(I6S_ANCOMall_.05pct_df$ANCOMcat == "foliar surface")) 
# True so this is working as expected
AlphaProFS_check == countsDistinctClassI6SANCOM$prop[intersect(which(countsDistinctClassI6SANCOM$ANCOMcat == "foliar surface"), which(countsDistinctClassI6SANCOM$Class == "Alphaproteobacteria"))]

########################################################################
#   III. FUNGI: EXPLORE ANCOM RESULTS (BIOAEROSOLS AND FOLIR SURFACE) 
########################################################################
length(which(ITS_ANCOMall_.05pct_df$ANCOMcat == "bioaerosol")) #152
length(which(ITS_ANCOMall_.05pct_df$ANCOMcat == "foliar surface")) #243

########### FUNGI ###########
##### FACETED DOT PLOT ####
#### TRY FIRST WITH FAMILY ####
# As calculated below, shows only 109/152 bioaerosol ASVs and 117/243 foliar surface ASVs, 
# so perhaps not the most representative.
# Get number of ASVs in each family
# Try at family level, as Noah suggested
numberFamsANCOM_ITS <- ITS_ANCOMall_.05pct_df %>% 
  group_by(ANCOMcat, Family) %>% 
  summarise(ASVsInFamByHab = n_distinct(ASV_name), .groups = "drop") %>%  #.groups = "drop" no more grouping
  complete(ANCOMcat, Family, fill = list(ASVsInFamByHab = 0)) #add in zero if not found in that ANCOM cat
# View(numberFamsANCOM_ITS)
table(ITS_ANCOMall_.05pct_df$Family, ITS_ANCOMall_.05pct_df$ANCOMcat) #matches calculations above!
length(ITS_ANCOMall_.05pct_df$ASV_name) #395 = 152 BA + 243 FS
# How many families should I try?
length(numberFamsANCOM_ITS$Family[which(numberFamsANCOM_ITS$ASVsInFamByHab >=3)]) #39 is tooo many!
length(numberFamsANCOM_ITS$Family[which(numberFamsANCOM_ITS$ASVsInFamByHab >=4)]) #32 is tooo many!
length(numberFamsANCOM_ITS$Family[which(numberFamsANCOM_ITS$ASVsInFamByHab >=5)]) #25 is probably possible, especially once NA removed
# Get names of families that have at least 3 representatives in one family or the other
fams5_ITS <- numberFamsANCOM_ITS$Family[which(numberFamsANCOM_ITS$ASVsInFamByHab >=5)]
fams5_ITS #two NAs in there
length(unique(fams5_ITS)) #24 unique families

# Looks good!
# View(numberFamsANCOM_ITS[numberFamsANCOM_ITS$Family %in% fams5_ITS,])
ITS_5fam_df <- numberFamsANCOM_ITS[numberFamsANCOM_ITS$Family %in% fams5_ITS,]
dim(ITS_5fam_df)
ITS_5fam_df <- ITS_5fam_df[-which(ITS_5fam_df$Family == "NA"),] #drop NA families
dim(ITS_5fam_df) #46, so 23 unique fams

# Add phylum info to this plot by merging with subset of ITS_ANCOMall_.05pct_df that has phylum info! DO distinct because otherwise subset is v long and get many-to-many relationship
famPhylaITS <- distinct(ITS_ANCOMall_.05pct_df[,colnames(ITS_ANCOMall_.05pct_df) %in% c("Family", "Phylum")])
ITS_5fam_df <- left_join(ITS_5fam_df, famPhylaITS, by = "Family")
# View(ITS_5fam_df) # looks good!
sum(ITS_5fam_df$ASVsInFamByHab) #shows 226/(152+243 = 395 diff abund ASVs)
sum(ITS_5fam_df$ASVsInFamByHab[which(ITS_5fam_df$ANCOMcat == "foliar surface")]) #117/243 foliar surface ASVs
sum(ITS_5fam_df$ASVsInFamByHab[which(ITS_5fam_df$ANCOMcat == "bioaerosol")]) #109/152 biaoerosol ASVs

# Re-order families (making it a factor) based on alphabetical in phylum than in family
ITS_5fam_df <- ITS_5fam_df %>%
  arrange(desc(Phylum), desc(Family)) %>%  #alphabetic by Phylum, then Family, descending so that it goes up y-axis
  mutate(Family = factor(Family, levels = unique(Family))) 
# View(ITS_5fam_df) #Looks good!
ITS_5fam_df$Family
unique(ITS_5fam_df$Family)
levels(ITS_5fam_df$Family)
# Change _fam_Incertae_sedis and make it 
levels(ITS_5fam_df$Family) <- gsub(x=levels(ITS_5fam_df$Family), pattern="_fam_Incertae_sedis", replacement = ", inc. sed.")

# Now remove those with 0 so that there are no dots for zero
ITS_5fam_df_no0s <- ITS_5fam_df[-which(ITS_5fam_df$ASVsInFamByHab == 0),]
# View(ITS_5fam_df_no0s)
levels(ITS_5fam_df_no0s$Family) #great kept levels!
# Make this numeric to squeeze left and right dots closer together
ITS_5fam_df_no0s$xpos <- as.numeric(factor(ITS_5fam_df_no0s$ANCOMcat))

# Make the plot!
colnames(ITS_5fam_df_no0s)
ITS_5famANCOM_bubPlot <- ggplot(ITS_5fam_df_no0s,
                                aes(x = xpos, y = Family, size = ASVsInFamByHab, fill = ANCOMcat)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_size(range = c(1, 6)) +   
  facet_grid(Family ~ ., scales = "free_y", space = "free_y", switch = "y", as.table = FALSE) + #as.table = FALSE orders by factor
  scale_x_continuous(
    breaks = 1:2,  # positions of the categories
    labels = levels(factor(ITS_5fam_df$ANCOMcat)),
    expand = expansion(mult = c(0.0, 0.0)),
    limits = c(0.2, 3)  # pull the two categories closer together
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "white", fill = NA, linewidth = 0.5),
    strip.placement = "outside",                 # put strips on the outer side (left)
    strip.text.y.right = element_text(angle = 0, hjust = 1, size =9, color= "black"),  # control orientation
    strip.background = element_blank(),
    panel.spacing.y = unit(0.01, "lines"),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size =5, color= "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    # Legend layout (stack the two guides)
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.box.just = "left",
    axis.text = element_text(color = "black"),
    legend.title = element_text(size = 8, color= "black"),
    legend.text  = element_text(size = 8, color= "black"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.1, "lines"),
    # Center titles above legend
    legend.title.position = "top",
    legend.title.align = 0.5,
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Number of ASVs in Family"
  ) +
  # Change bubble colors
  scale_fill_manual(
    values = c("cornflowerblue", "forestgreen") #color bioaerosol bubbles blue and foliar surface green
  ) +
  guides(fill = "none") #remove fill legend for ANCOMcat

#quartz(width = 6, height = 7)
ITS_5famANCOM_bubPlot
#
# NOTE: LEGEND WILL NOT LEFT JUSTIFY NOR CAN I REMOVE Y-AXIS TICKS FROM THE RIGHT SIDE, SO I WILL
# MANUALLY DO THESE THINGS IN POWERPOINT

#### ORDER ####
# As calculated below, shows 144/152 bioaerosol ASVs and 201/243 foliar surface ASVs
# Get number of ASVs in each order
numberOrdersANCOM_ITS <- ITS_ANCOMall_.05pct_df %>% 
  group_by(ANCOMcat, Order) %>% 
  summarise(ASVsInOrderByHab = n_distinct(ASV_name), .groups = "drop") %>%  #.groups = "drop" no more grouping
  complete(ANCOMcat, Order, fill = list(ASVsInOrderByHab = 0)) #add in zero if not found in that ANCOM cat
# View(numberOrdersANCOM_ITS)
table(ITS_ANCOMall_.05pct_df$Order, ITS_ANCOMall_.05pct_df$ANCOMcat) #matches calculations above!
# How many orders should I try?
length(numberOrdersANCOM_ITS$Order[which(numberOrdersANCOM_ITS$ASVsInOrderByHab >=3)]) #23 would work!!
length(numberOrdersANCOM_ITS$Order[which(numberOrdersANCOM_ITS$ASVsInOrderByHab >=2)]) #27 would probably work too!!
# Get names of orders that have at least 2 representatives in one order or the other
orders2_ITS <- numberOrdersANCOM_ITS$Order[which(numberOrdersANCOM_ITS$ASVsInOrderByHab >=2)]
orders2_ITS #two NAs in there
length(unique(orders2_ITS)) #26 unique orders

# Looks good!
# View(numberOrdersANCOM_ITS[numberOrdersANCOM_ITS$Order %in% orders2_ITS,])
ITS_2order_df <- numberOrdersANCOM_ITS[numberOrdersANCOM_ITS$Order %in% orders2_ITS,]
dim(ITS_2order_df)
ITS_2order_df <- ITS_2order_df[-which(ITS_2order_df$Order == "NA"),] #drop NA orders
dim(ITS_2order_df) #50, so 25 unique orders

# Add phylum info to this plot by merging with subset of ITS_ANCOMall_.05pct_df that has phylum info! DO distinct because otherwise subset is v long and get many-to-many relationship
orderPhylaITS <- distinct(ITS_ANCOMall_.05pct_df[,colnames(ITS_ANCOMall_.05pct_df) %in% c("Order", "Phylum")])
ITS_2order_df <- left_join(ITS_2order_df, orderPhylaITS, by = "Order")
# View(ITS_2order_df) # looks good!
sum(ITS_2order_df$ASVsInOrderByHab) #shows 345/395 diff abund ASVs
sum(ITS_2order_df$ASVsInOrderByHab[which(ITS_2order_df$ANCOMcat == "foliar surface")]) #201/243 foliar surface ASVs
sum(ITS_2order_df$ASVsInOrderByHab[which(ITS_2order_df$ANCOMcat == "bioaerosol")]) #144/152 biaoerosol ASVs

# Re-order orders (making it a factor) based on alphabetical in phylum than in order
ITS_2order_df <- ITS_2order_df %>%
  arrange(desc(Phylum), desc(Order)) %>%  #alphabetic by Phylum, then Order, descending so that it goes up y-axis
  mutate(Order = factor(Order, levels = unique(Order))) 
# View(ITS_2order_df) #Looks good!
ITS_2order_df$Order
unique(ITS_2order_df$Order)
levels(ITS_2order_df$Order)
# Remove _ord_Incertae_sedis. Note that Cystobasidiomycetes is inc. sed in caption
levels(ITS_2order_df$Order) <- gsub(x=levels(ITS_2order_df$Order), pattern="_ord_Incertae_sedis", replacement = "")
levels(ITS_2order_df$Order)

# Now remove those with 0 so that there are no dots for zero
ITS_2order_df_no0s <- ITS_2order_df[-which(ITS_2order_df$ASVsInOrderByHab == 0),]
# View(ITS_2order_df_no0s)
levels(ITS_2order_df_no0s$Order) #great kept levels!
# Make this numeric to squeeze left and right dots closer together
ITS_2order_df_no0s$xpos <- as.numeric(factor(ITS_2order_df_no0s$ANCOMcat))

# Make the plot!
colnames(ITS_2order_df_no0s)
ITS_2ordANCOM_bubPlot <- ggplot(ITS_2order_df_no0s,
                                aes(x = xpos, y = Order, size = ASVsInOrderByHab, fill = ANCOMcat)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_size(range = c(1, 6)) +   
  facet_grid(Order ~ ., scales = "free_y", space = "free_y", switch = "y", as.table = FALSE) + #as.table = FALSE orders by factor
  scale_x_continuous(
    breaks = 1:2,  # positions of the categories
    labels = levels(factor(ITS_2order_df$ANCOMcat)),
    expand = expansion(mult = c(0.0, 0.0)),
    limits = c(0.2, 3)  # pull the two categories closer together
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "white", fill = NA, linewidth = 0.5),
    strip.placement = "outside",                 # put strips on the outer side (left)
    strip.text.y.right = element_text(angle = 0, hjust = 1, size =9, color= "black"),  # control orientation
    strip.background = element_blank(),
    panel.spacing.y = unit(0.01, "lines"),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size =5, color= "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    # Legend layout (stack the two guides)
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.box.just = "left",
    axis.text = element_text(color = "black"),
    legend.title = element_text(size = 8, color= "black"),
    legend.text  = element_text(size = 8, color= "black"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.1, "lines"),
    # Center titles above legend
    legend.title.position = "top",
    legend.title.align = 0.5,
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Number of ASVs in Order"
  ) +
  # Change bubble colors
  scale_fill_manual(
    values = c("cornflowerblue", "forestgreen") #color bioaerosol bubbles blue and foliar surface green
  ) +
  guides(fill = "none") #remove fill legend for ANCOMcat

# quartz(width = 6, height = 7)
ITS_2ordANCOM_bubPlot

# saveRDS(ITS_2ordANCOM_bubPlot, file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_2ordANCOM_bubPlot_10-16-2025.rds")
# NOTE: LEGEND WILL NOT LEFT JUSTIFY NOR CAN I REMOVE Y-AXIS TICKS FROM THE RIGHT SIDE, SO I WILL
# MANUALLY DO THESE THINGS IN POWERPOINT

# WHAT ARE TOP ORDERS FOR BIOAEROSOL AND FOLIAR SURFACE ASVS? (ADDED TO MANUSCRIPT)
head(ITS_ANCOMall_.05pct_df)
countsDistinctOrderITSANCOM <- ITS_ANCOMall_.05pct_df %>%
  group_by(ANCOMcat, Order) %>% 
  summarise(nASVs = n_distinct(ASV_name), .groups = "drop") %>% #number of distinct ASVs per group then drop group
  group_by(ANCOMcat) %>%
  mutate(totalInANCOMcat = sum(nASVs),
         prop = nASVs / totalInANCOMcat,
         pct  = 100 * prop) %>% #pct is percentage
  ungroup()
# View(countsDistinctOrderITSANCOM)
# Do these add up to one? YES!
sum(countsDistinctOrderITSANCOM$prop[which(countsDistinctOrderITSANCOM$ANCOMcat == "bioaerosol")])
sum(countsDistinctOrderITSANCOM$prop[which(countsDistinctOrderITSANCOM$ANCOMcat == "foliar surface")])

# A check of the code above
pleoFS_check <- length(intersect(which(ITS_ANCOMall_.05pct_df$Order == "Pleosporales"), 
                 which(ITS_ANCOMall_.05pct_df$ANCOMcat == "foliar surface")))/
  length(which(ITS_ANCOMall_.05pct_df$ANCOMcat == "foliar surface")) 
# True so this is working as expected
pleoFS_check == countsDistinctOrderITSANCOM$prop[intersect(which(countsDistinctOrderITSANCOM$ANCOMcat == "foliar surface"), which(countsDistinctOrderITSANCOM$Order == "Pleosporales"))]

########################################################################
#   V. BACTERIAL SPORULATION
########################################################################
# Read in data and clean up for use
head(ANCOMspore)
str(ANCOMspore)
which(is.na(ANCOMspore$spore4cats)==FALSE)
# Keep only those ASVs with a match in Madin, i.e., those that DON'T have NA as their sporulation category
ANCOMspore <- ANCOMspore[which(is.na(ANCOMspore$spore4cats)==FALSE),]
# View(ANCOMspore)
nrow(ANCOMspore) #111 had data
111/220*100 #55.5% matches percentage reported that I calculated in bacteriaSporesAndPigment

# Split these into separate air and phyllosphere dataframes
ANCOMspore_air <- ANCOMspore[which(ANCOMspore$ANCOMcat == "bioaerosol"),]
unique(ANCOMspore_air$ANCOMcat)
nrow(ANCOMspore_air) #61
ANCOMspore_foliar <- ANCOMspore[which(ANCOMspore$ANCOMcat == "foliar surface"),]
unique(ANCOMspore_foliar$ANCOMcat)
nrow(ANCOMspore_foliar) #50 with information

# INVESTIGATE PATTERNS AND MAKE PIE CHARTS -- tentatively supports hypothesis!
table(ANCOMspore$ANCOMcat, ANCOMspore$spore4cats)
#                 no possible yes
# bioaerosol     41       13   7
# foliar surface 43        5   2

yesOrPossibleSporeIndex <- c(which(ANCOMspore$spore4cats == "yes"), which(ANCOMspore$spore4cats == "possible"))
sporeYesOrPossible <- ANCOMspore[yesOrPossibleSporeIndex,]
# View(sporeYesOrPossible)
unique(sporeYesOrPossible$Order) #11
unique(sporeYesOrPossible$Class) #  "Bacilli"              "Thermoanaerobacteria" "Clostridia"           "Actinobacteria"       "Myxococcia"          
# [6] "Gammaproteobacteria"  "Alphaproteobacteria" 
# Some Ramlibacter (Alphaproteobacteria, produce cysts!!!)
# Confirmed too that some Pseudomonas and Sphingomonadaceae in Madin database are classified as spore producing
head(sporeYesOrPossible)
sporeYesOrPossible_order <- sporeYesOrPossible %>% 
  group_by(ANCOMcat, Order) %>% 
  summarise(n_ASVs = n_distinct(ASV_name), .groups = "drop")

sporeAirPhyllo_barPlot <- ggplot(sporeYesOrPossible_order) +
  geom_bar(aes(x = ANCOMcat, y = n_ASVs, fill = Order), 
           position = "stack", stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 14) 
  ) +
  labs(
    fill = "Order",  #change legend title
    x = NULL,
    y = "Number of ASVs"  
  )


##### STACKED BARPLOTS:
# Make dataframe
sporeLeafAir_df <- as.data.frame(matrix(ncol=3, nrow=6))
colnames(sporeLeafAir_df) <- c("ANCOMcat", "numberASVs", "spore4cats")
sporeLeafAir_df[,1] <- c(rep("bioaerosol",3), rep("foliar surface", 3)) 
sporeLeafAir_df
sporeLeafAir_df[,3] <- c(rep(c("no", "possible", "yes"), 2))
sporeLeafAir_df[,2] <- c(unname(table(ANCOMspore_air$spore4cats)), unname(table(ANCOMspore_foliar$spore4cats)))
sporeLeafAir_df #looks good,  but double check by inspecting that these match the tables below:
table(ANCOMspore_air$spore4cats)
table(ANCOMspore_foliar$spore4cats)

sporeAirPhyllo_barPlot <- ggplot(sporeLeafAir_df) +
  geom_bar(aes(x = ANCOMcat, y = numberASVs, fill = spore4cats), 
           position = "dodge", stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 14) 
  ) +
  scale_fill_manual(values = c("#6C244C", "#ffe200", "#02AD24")) +
  labs(
    fill = "Sporulation Ability",  #change legend title
    x = NULL,
    y = "Number of ASVs"  
  )

sporeAirPhyllo_barPlot 

# With the legend on the bottom:
sporeAirPhyllo_barPlot_LegBott <- ggplot(sporeLeafAir_df) +
  geom_bar(aes(x = ANCOMcat, y = numberASVs, fill = spore4cats), 
           position = "dodge", stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.position = "bottom",
    legend.justification = "center" 
  ) +
  scale_fill_manual(values = c("#6C244C", "#ffe200", "#02AD24")) +
  scale_x_discrete(labels = c("air", "foliar surfaces")) +
  labs(
    fill = "Sporulation Ability",  #change legend title
    x = NULL,
    y = "Number of ASVs "  
  )

# October 15, 2025
# saveRDS(sporeAirPhyllo_barPlot_LegBott, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/sporeAirFoliar_barPlot_10-15-2025.rds")

# INVESTIGATION OF TOP TAXA IN BIOAEROSOLS AND SPORULATION ABILITY
# Are some of the top taxa in Table S5 that are found in bioaerosols but not on foliar surfaces, spore-formers?
# View(ANCOMspore_air)
# Note that ASV 126, Deinococcus, had no sporulation data
ANCOMspore_air[which(ANCOMspore_air$ASV_name == "ASV_10"),] #ASV_10, Staphylococcus, not spore former as expected
ANCOMspore_air[which(ANCOMspore_air$ASV_name == "ASV_48"),] #48, Cystobacter, IS A SPORE FORMER
ANCOMspore_air[which(ANCOMspore_air$ASV_name == "ASV_129"),] #129, Ligilactobacillus, is not a spore former
ANCOMspore_air[which(ANCOMspore_air$ASV_name == "ASV_142"),] #not found
ANCOMspore_air[which(ANCOMspore_air$ASV_name == "ASV_239"),] #Anoxybacillus is possible. Low in soil, not on foliar surfaces
ANCOMspore_air[which(ANCOMspore_air$ASV_name == "ASV_198"),] #not found
# "ASV_1728" "ASV_1855" "ASV_2180"
ANCOMspore_air$ASV_name[which(ANCOMspore_air$Genus == "[Ruminococcus] torques group")] #all of these are "possible" spore formerand are bioaerosol
ANCOMspore_air$ASV_name[which(ANCOMspore_air$Genus == "Pseudomonas")] #all of these are "possible" spore formers

# ADD SPORULATION DATA TO TABLE S5;
I6S_topAirTable_forSupp <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_topAirTable_03-24-26.rds")

colnames(ANCOMspore_air)
# Make a copy and add back "ASV_" for merging
I6S_topAirTable_forSupp2 <- I6S_topAirTable_forSupp
I6S_topAirTable_forSupp2$ASV_name <- paste0("ASV_", I6S_topAirTable_forSupp2$ASV_name)
top16S_plusANCOMspores <- left_join(I6S_topAirTable_forSupp2, ANCOMspore_air[,c(1, 148,149, 156, 157, 158)], by = "ASV_name")
# View(top16S_plusANCOMspores)

top16S_plusANCOMspores[which(top16S_plusANCOMspores$Genus == "Pseudomonas"),]

# Look at all that are only in bioaerosols
View(top16S_plusANCOMspores[intersect(which(top16S_plusANCOMspores$SoilPercentOccupancy == 0.0), which(top16S_plusANCOMspores$FoliarPercentOccupancy == 0.0)),])

#### PERFORM chi-squared tests to SEE IF AIR AND PHYLLOSPHERE VARY IN SPORULATION 
# (TESTED NO AGAINST PROPORTION OF YES/POTENTIALLY)

# i. SET UP
ANCOMspore
contTableANCOMgroups <- table(ANCOMspore$ANCOMcat, ANCOMspore$spore4cats)
contTableANCOMgroups #looks like before, excellent

# ii. # FISHER EXACT TEST:
# According to Jim Frost: "Fisher’s exact test doesn’t use the chi-square statistic and sampling distribution. 
# Instead, it calculates the number of all possible contingency tables with the same row and column totals 
# (i.e., marginal distributions) as the observed table. Then it calculates the probability for the p-value
# by finding the proportion of possible tables that are more extreme than the observed table. 

# Requires a 2 x 2 contingency table, so will combine possible and yes
contTableSporesCombined <- contTableANCOMgroups #make a copy 
contTableSporesCombined[,2] <- contTableANCOMgroups[,2] + contTableANCOMgroups[,3]
contTableSporesCombined <- contTableSporesCombined[,-3]
colnames(contTableSporesCombined)[2] <- "possibleAndYes"
class(contTableSporesCombined)
contTableSporesCombined

mosaicplot(contTableSporesCombined,
           main = "Mosaic plot",
           color = TRUE
)

# Can flip this to be more natural (odds ratio in R defaults to reporting odds
# ratio from first cell-perspective)
contTableSporesCombinedReordered <- cbind(contTableSporesCombined[,2], contTableSporesCombined[,1])
colnames(contTableSporesCombinedReordered) <- c("possibleAndYes", "no")
contTableSporesCombinedReordered
contTableANCOMgroups #double check that this adds up like the original

set.seed(18)
sporesFisherTest <- fisher.test(contTableSporesCombinedReordered)
sporesFisherTest

# Fisher's Exact Test for Count Data
# 
# data:  contTableSporesCombinedReordered
# p-value = 0.02663
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.064311 9.224660
# sample estimates:
# odds ratio 
#   2.967765 

########################################################################
#   VI. BACTERIAL PIGMENTATION!
########################################################################
ANCOMgenPig

# Make stacked barplot
# Make dataframe
pigLeafAir_df <- as.data.frame(matrix(ncol=3, nrow=6))
colnames(pigLeafAir_df) <- c("ANCOMcat", "numberASVs", "PigConsensus")
table(ANCOMgenPig$ANCOMcat, ANCOMgenPig$PigConsensus)
#                 no yes yesAndNo
# bioaerosol      3  16       33
# foliar surface  2  21       18

# How many available?
length(unique(ANCOMgenPig$ASV_name))/length(unique(I6S_ANCOMall_.05pct_df$ASV_name))*100 #42.27273, so 42.3% in MS

pigLeafAir_df[,1] <- c(rep("bioaerosol",3), rep("foliar surface", 3))
pigLeafAir_df
pigLeafAir_df[,3] <- c(rep(c("no", "possible", "yes"), 2))
pigLeafAir_df
table(ANCOMgenPig$ANCOMcat, ANCOMgenPig$PigConsensus)
pigLeafAir_df[,2] <- c(3, 33, 16, 2, 18, 21)
pigLeafAir_df #looks good,  but double check by inspecting that these match the tables below:
table(ANCOMgenPig$ANCOMcat, ANCOMgenPig$PigConsensus)

# Make plot
colnames(ANCOMgenPig)
ANCOMgenPigPlot <- ggplot(pigLeafAir_df) +
  geom_bar(aes(x = ANCOMcat, y = numberASVs, fill = PigConsensus),
           position = "dodge", stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.position = "bottom",
    legend.justification = "center" 
  ) +
  scale_fill_manual(values = c("#6C244C", "#ffe200", "#02AD24")) +
  labs(
    fill = "Pigmentation Ability",  #change legend title
    x = NULL,
    y = "Number of ASVs"
  )

ANCOMgenPigPlot
# October 15, 2025
# saveRDS(ANCOMgenPigPlot, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ANCOMgenPigPlot_10-15-2025.rds")

#### PERFORM chi-squared tests to SEE IF AIR AND PHYLLOSPHERE VARY IN PIGMENTATION PRODUCTION ABILITY
# (TESTED NO AGAINST PROPORTION OF YES/POTENTIALLY)

# i. SET UP
pigContigencyTable <- table(ANCOMgenPig$ANCOMcat, ANCOMgenPig$PigConsensus)

# ii. # FISHER EXACT TEST:
# According to Jim Frost: "Fisher’s exact test doesn’t use the chi-square statistic and sampling distribution.
# Instead, it calculates the number of all possible contingency tables with the same row and column totals
# (i.e., marginal distributions) as the observed table. Then it calculates the probability for the p-value
# by finding the proportion of possible tables that are more extreme than the observed table.

# Requires a 2 x 2 contingency table, so will combine possible/yesAndNo and yes
pigContigencyTableCombined <- pigContigencyTable #make a copy
pigContigencyTableCombined[,2] <- pigContigencyTable[,2] + pigContigencyTable[,3]
pigContigencyTableCombined <- pigContigencyTableCombined[,-3] #remove third column
colnames(pigContigencyTableCombined)[2] <- "possibleAndYes" #rename second column
class(pigContigencyTableCombined)
pigContigencyTableCombined

# Check our mosaic plot
mosaicplot(pigContigencyTableCombined,
           main = "Mosaic plot",
           color = TRUE
)

# Can flip this to be more natural (odds ratio in R defaults to reporting odds
# ratio from first cell-perspective)
pigContigencyTableCombinedReordered <- cbind(pigContigencyTableCombined[,2], pigContigencyTableCombined[,1])
colnames(pigContigencyTableCombinedReordered) <- c("possibleAndYes", "no")
pigContigencyTableCombinedReordered
pigContigencyTable  #double check that this adds up like the original and it does

set.seed(21)
pigmentFisherTest <- fisher.test(pigContigencyTableCombinedReordered)
pigmentFisherTest

# Fisher's Exact Test for Count Data
# 
# data:  pigContigencyTableCombinedReordered
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.06702525 7.70672979
# sample estimates:
# odds ratio 
#  0.8391869 

############# SPORE SIZE NOW IN SCRIPT sporeSizeSept5.R ########

########################################################################
#   VI. FUNGAL FRUITING BODIES
########################################################################
########
#         FORMATTING FOR FUNGuild 
########
# Make a dataframe to feed into FUNGuild 
# format is tsv, where first row should be OTU_IOTU_ID(tab)sample1(tab)sample2(tab)sample3(tab)taxonomy(return)
# taxonomic levels should be separated by semicolons
# Remove ASVs with no reads (i.e. were only found in soil)
ITSairPhylloOnly.ps <- prune_taxa(taxa_sums(ITSairPhylloOnly.ps) > 0, ITSairPhylloOnly.ps)
min(taxa_sums(ITSairPhylloOnly.ps)) #no non-present taxa!
airFoliar_ITS_ASVs <- as.data.frame(as.matrix(otu_table(ITSairPhylloOnly.ps))) # ASV table out of phyloseq made above; ASVs are rows
dim(airFoliar_ITS_ASVs) #9,012  169
airFoliar_ITS_tax <- as.data.frame(as.matrix(tax_table(ITSairPhylloOnly.ps)))# tax table out of phyloseq made above; ASVs are also rows!
dim(airFoliar_ITS_tax) #9012    7
colnames(airFoliar_ITS_tax)
unique(rownames(airFoliar_ITS_tax) == rownames(airFoliar_ITS_ASVs)) #TRUE
airFoliarFUNGuildTable1 <- merge(airFoliar_ITS_ASVs, airFoliar_ITS_tax, by= "row.names", all=TRUE)
colnames(airFoliarFUNGuildTable1)[1] <- "OTU_ID"
# Collapse taxonomy into one column 
airFoliarFUNGuildTable1 <- tidyr::unite(airFoliarFUNGuildTable1, sep=";",col= "taxonomy", Kingdom:Species)
# View(airFoliarFUNGuildTable1)
head(airFoliarFUNGuildTable1) #has species info as we want as one big taxonomy with ; separation!
colnames(airFoliarFUNGuildTable1)
dim(airFoliarFUNGuildTable1)

#### HERE: SHOULD CHECK THAT SAME OUTPUT AS ON SERVER, i.e.
# Written to file October 16, 2025
# write.table(airFoliarFUNGuildTable1, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airFoliarFUNGuildTable1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# Code above should make it so that no quotation marks and rownames appear, but I confirmed this by opening in Excel anyway
# Note: double checked here: https://github.com/UMNFuN/FUNGuild/blob/master/README.md that the version I have, Guilds_v1.1.py, is the most current
# SEE THE FILE Desktop/CU_Research/SoilEdgeEffectsResearch/EdgeEffectsOnSoil/FUNGuildStepByStep/FUNGuildStepByStep for a more detailed guide to how I ran FUNGuild I ran FUNGuild
# Ran the following code in Terminal on personal MacBook:
# python3 ~/Desktop/CU_Research/SoilEdgeEffectsResearch/FUNGuild/Guilds_v1.1.py -otu ~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airFoliarFUNGuildTable1.txt -m -db fungi

# This created the output file in the same folder where the input ASV table was found:
# ~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airFoliarFUNGuildTable1.guilds_matched.txt !

# Load raw results file (note that this has ALL fungal ASVs in air and foliar surfaces, not just ANCOM):
rawFungResults <- read.delim(file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airFoliarFUNGuildTable1.guilds_matched.txt", header = T)
# View(rawFungResults) 
dim(rawFungResults)
#  6422 out of 9012 had results in fungGuild OTUs!!
# How many fungal taxa could we match?
6422/9012 *100 #71.26054, or 71.3% of ASVs had a match, but this is overall, not ANCOM yet
# Compare with old results:  3485/4673*100 #74.58% were matched
colnames(rawFungResults)[1] <- "ASV_name"
dim(rawFungResults)
unique(rawFungResults$Confidence.Ranking) #"Probable"        "Highly Probable" "Possible" 

# Get only those ASVs with probable or highly probable results (i.e. NO Possible)
highProbIndex <- which(rawFungResults$Confidence.Ranking == "Highly Probable")
probIndex <- which(rawFungResults$Confidence.Ranking == "Probable") 
fungAssignedIndex <- c(highProbIndex, probIndex)
fungResults <- rawFungResults[fungAssignedIndex,]
dim(fungResults) # 5710  180
head(fungResults)
colnames(fungResults)
5710/9012*100 #63.35996% or 63.4% of all fungal ASVs had high probably or highly probable info!!
# OLD RESULTS: 3074/4673*100 #65.78% were matched

# INVESTIGATE WITH ANCOM RESULTS
head(ITS_ANCOMall_.05pct_df)
# Get subsetted columns because do not need whole ASV table twice!!!
subsetFungCols <- colnames(fungResults) %in% c("ASV_name" , "taxonomy", "Taxon", "Taxon.Level", "Trophic.Mode", "Guild", "Growth.Morphology","Trait","Confidence.Ranking")
length(unique(ITS_ANCOMall_.05pct_df$ASV_name)) #395
# Join with ANCOM results results with an inner join so that only ANCOM pulled out
airLeafANCOM_FUNGUILD <- inner_join(ITS_ANCOMall_.05pct_df, x=fungResults[,subsetFungCols], by="ASV_name")
length(unique(airLeafANCOM_FUNGUILD$ASV_name)) #300
# This shows that 512/692 differentially abundant were matched!
length(unique(airLeafANCOM_FUNGUILD$ASV_name))/length(unique(ITS_ANCOMall_.05pct_df$ASV_name))*100 #75.94937, or 75.9%, but see below because this is before taxa w/o growth morphology removed

# Make a smaller dataset without specific sample abundances 
colnames(airLeafANCOM_FUNGUILD)
airLeafANCOM_FUNGUILD_Smaller <- airLeafANCOM_FUNGUILD %>% 
  select(ASV_name, ANCOMcat, Trophic.Mode, Growth.Morphology, 
         Guild, Order, Family, Genus, Species, foliarSurfaceOcc, bioaerosolOcc, 
         pctBioaerosol, pctFoliar)
# View(airLeafANCOM_FUNGUILD_Smaller)

length(which(airLeafANCOM_FUNGUILD_Smaller$Growth.Morphology == "NULL")) #14 were unassigned
dim(airLeafANCOM_FUNGUILD_Smaller[-which(airLeafANCOM_FUNGUILD_Smaller$Growth.Morphology == "NULL"),]) #286  10
dim(airLeafANCOM_FUNGUILD_Smaller[-which(airLeafANCOM_FUNGUILD_Smaller$Growth.Morphology == "NULL"),])[1]/length(unique(ITS_ANCOMall_.05pct_df$ASV_name))*100 #72.40506 or 72.4% had growth morphology data!!!

# Remove those with no growth morph info
airLeafANCOM_FG_NoNullMorph <- airLeafANCOM_FUNGUILD_Smaller[-which(airLeafANCOM_FUNGUILD_Smaller$Growth.Morphology == "NULL"),]
dim(airLeafANCOM_FG_NoNullMorph)

# What proportion of the original air or phyllo species could be matched for growth morphology?
airMorphNumber <- length(which(airLeafANCOM_FG_NoNullMorph$ANCOMcat == "bioaerosol")) #number of ASVs enriched in air with morph data = 130
airANCOMNumber <- length(which(ITS_ANCOMall_.05pct_df$ANCOMcat == "bioaerosol")) #152 original bioaerosol taxa
airMorphNumber/airANCOMNumber*100 #85.52632% or 85.5% of air enriched taxa has morphological data

foliarMorphNumber <- length(which(airLeafANCOM_FG_NoNullMorph$ANCOMcat == "foliar surface")) #number of ASVs enriched in foliar data with morph data = 156
foliarANCOMNumber <- length(which(ITS_ANCOMall_.05pct_df$ANCOMcat == "foliar surface")) #number of ASVs enriched in foliar = 243
foliarMorphNumber/foliarANCOMNumber*100 #64.19753% or 64.2% of foliar-enriched taxa could be matched :)

MorphAirLeavesTable <- table(airLeafANCOM_FG_NoNullMorph$Growth.Morphology, airLeafANCOM_FG_NoNullMorph$ANCOMcat)
MorphAirLeaves_df <- as.data.frame(MorphAirLeavesTable)
colnames(MorphAirLeaves_df) <- c("morphology", "ANCOMcat", "freq")
# View(MorphAirLeaves_df)

# Add Unclassified to show number that could be matched for bioaerosols and foliar surface
# 1. Add a new level, "unclassified"
levels(MorphAirLeaves_df$morphology) <- c(levels(MorphAirLeaves_df$morphology), "Unclassified")
airANCOMNumber - airMorphNumber #22 could not be matched (22/243+152 could not be matched)
airMorphNumber/airANCOMNumber*100 #85.52632% could be classified (matches above)
# View(MorphAirLeaves_df)
nrow(MorphAirLeaves_df) 
# Make last row unclassified bioaerosol and add in number unclassified
MorphAirLeaves_df[nrow(MorphAirLeaves_df) +1,] <- c("Unclassified", "bioaerosol", airANCOMNumber - airMorphNumber)
foliarANCOMNumber - foliarMorphNumber #87 (out of 243 foliar surface ANCOM cat) could NOT be matched
foliarMorphNumber/foliarANCOMNumber*100 #64.19753 could be matched, matches above!
MorphAirLeaves_df[nrow(MorphAirLeaves_df) +1,] <- c("Unclassified", "foliar surface", foliarANCOMNumber - foliarMorphNumber)
MorphAirLeaves_df$freq <- as.numeric(MorphAirLeaves_df$freq)
class(MorphAirLeaves_df$freq)
# View(MorphAirLeaves_df)

# WHAT PERCENTAGE OF BIOAEROSOL AND FOLIAR SURFACE ASSOCIATED TAXA belonged to each type of
# "growth morphology"?
# AIR
MorphAir_only <- MorphAirLeaves_df %>% 
  filter(morphology != "Unclassified") %>% #filter out unclassified
  filter(ANCOMcat == "bioaerosol") %>% 
  mutate(percentMorph = (freq/sum(freq))*100) %>% 
  arrange(desc(percentMorph))
MorphAir_only
# View(MorphAir_only)
sum(MorphAir_only$freq[-nrow(MorphAir_only)]) # 130 matches which is correct (matches number of ANCOM taxa that 
# I was able to assign morphology information), Did -nrow(MorphAir_only) to subtract out Unclassified

# FOLIAR
MorphLEAF_only <- MorphAirLeaves_df %>% 
  filter(morphology != "Unclassified") %>% #filter out unclassified
  filter(ANCOMcat == "foliar surface") %>% 
  mutate(percentMorph = (freq/sum(freq))*100) %>% 
  arrange(desc(percentMorph))
MorphLEAF_only
# View(MorphLEAF_only)
sum(MorphLEAF_only$freq[-nrow(MorphLEAF_only)]) #156 matches which is correct (matches number of ANCOM taxa that
# I was able to assign morphology information), Did -nrow(MorphAir_only) to subtract out Unclassified

# MAKE A PLOT THAT SHOWS ANCOM CATEGORY AND MORPHOLOGY
# 1. Get which morphologies to include
length(unique(MorphAirLeaves_df$morphology)) #16 unique, but this is too many
numMorphCat <- MorphAirLeaves_df %>% 
  group_by(morphology) %>% 
  summarize(n= sum(freq)) %>% 
  arrange(n)
print(numMorphCat, n= nrow(numMorphCat))
# View(numMorphCat)
# Will use only those morphologies with at least 2 representative ASVs
manyMorphs <- numMorphCat$morphology[which(numMorphCat$n >= 2)]
manyMorphs <- as.character(manyMorphs) #remove factor levels
# View(numMorphCat) #checked that these match!

# 2. Get a smaller dataframe just with these morphologies
MorphAirLeaves_2_df <- MorphAirLeaves_df #make a copy to edit or stupid factors will mess up
MorphAirLeaves_2_df$morphology <- as.character(MorphAirLeaves_2_df$morphology) #remove factor levels
is.factor(MorphAirLeaves_2_df$morphology) #good this is false
# Remove rarer morphologies
MorphAirLeaves_2_df <- MorphAirLeaves_2_df[which(MorphAirLeaves_2_df$morphology %in% manyMorphs),]
# Remove unclassified
MorphAirLeaves_2_df <- MorphAirLeaves_2_df[-which(MorphAirLeaves_2_df$morphology == "Unclassified"),]
# View(MorphAirLeaves_2_df)
setdiff(unique(MorphAirLeaves_2_df$morphology), manyMorphs) #no differences, yay!
# View(MorphAirLeaves_2_df)
length(unique(MorphAirLeaves_2_df$morphology)) #12

# 3. Set colors
carto_pal(12, "Safe")
morphsToPlot <- unique(MorphAirLeaves_2_df$morphology)
class(morphsToPlot); is.vector(morphsToPlot)
length(morphsToPlot)
colsForANCOMmorphs <- setNames(
  c("#FE8F42","#CC6677","#117733","#DDCC77","#332288",
    "#44AA99","#999933","#882255","#88CCEE","#661100",
    "#6699CC", "#AA4499"),
  morphsToPlot[1:12]
)
# saveRDS(colsForANCOMmorphs, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/colsForANCOMmorphs.rds")
# 4. Make the plot!
airLeafANCOM_FUNGUILD_2_barPlot <- ggplot(MorphAirLeaves_2_df) +
  geom_bar(aes(x = ANCOMcat, y = freq, fill = morphology), 
           position = "stack", stat = "identity") +
  theme_bw() +
  #scale_x_discrete(labels = c("air"= "air", "phyllo" = "foliar surfaces")) + #option to change colors
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 11) 
  ) +
  scale_fill_manual(values = colsForANCOMmorphs) +
  labs(
    fill = "Growth Morphology",  #change legend title
    x = NULL,
    y = "Number of ASVs"  
  )

airLeafANCOM_FUNGUILD_2_barPlot 
# saveRDS(airLeafANCOM_FUNGUILD_2_barPlot, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airLeafANCOM_FUNGUILD_2_barPlot_Oct16_2025.rds")
# CRUST FUNGI (Corticoid): https://www.crustfungi.com/html/sidebar/introduction.html#:~:text=Morphological%20Definition,or%20hydnoid%20hymenophore;%20and%20holobasidia.

# EXPLORING SPECIFIC ASVS
# Capnodiales, the leaf specialists!
# View(airLeafANCOM_FG_NoNullMorph[airLeafANCOM_FG_NoNullMorph$Order=="Capnodiales",])
# Take a look at all matches for leaf surfaces
# View(airLeafANCOM_FG_NoNullMorph[airLeafANCOM_FG_NoNullMorph$ANCOMcat=="foliar surface",])

# View(airLeafANCOM_FG_NoNullMorph[airLeafANCOM_FG_NoNullMorph$ANCOMcat=="bioaerosol",])

# TRAMETES EXAMPLE (in manuscript). #Shows 5 ASVs, all with polyporoid growth morphology
airLeafANCOM_FG_NoNullMorph[airLeafANCOM_FG_NoNullMorph$Genus=="Trametes",]
# Get % of total fungal reads these ANCOM-identified Trametes comprised
sum(airLeafANCOM_FG_NoNullMorph$pctBioaerosol[airLeafANCOM_FG_NoNullMorph$Genus=="Trametes"]) #3.377989, or 3.4% of all biaoerosol reads
# Get occupancy in bioaerosol samples. This is the average that each ASV is in in bioaerosol samples
mean(airLeafANCOM_FG_NoNullMorph$bioaerosolOcc[airLeafANCOM_FG_NoNullMorph$Genus=="Trametes"] )/110 #0.9036364, 90% added to paper

