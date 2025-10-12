# ANCOM_traits_I6SandITS_Sept.R
# Old version was called ANCOM_traits_I6SandITS.R. Old, old version was called DiffAbundTests_I6SandITS.R

# This script follows ANCOMBC_robust_foliarAir_pt0.05 (own computer) and bacterialSporeFormers_Sept.R (run on server, but saved in both locations)

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

# Not read in right now since making bubble plots first
# # 3. READ IN SPORULATION AND PIGMENTATION INFORMATION (BOTH made in bacterialSporeFormers_Sept.R on server)
# ANCOMspore <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airLeafANCOM_withSpore_Oct2_2025.RData")
# ANCOMgenPig <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ANCOM_ijsemGenusPig_10-03-2025.rds")

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

####################################
##### ANOTHER WAY OF SHOWING THE DATA: FACETED DOT PLOT ####
# Get number of ASVs in each family
# Try at family level, as Noah suggested
numberFamsANCOM_I6S <- I6S_ANCOMall_.05pct_df %>% 
  group_by(ANCOMcat, Family) %>% 
  summarise(ASVsInFamByHab = n_distinct(ASV_name), .groups = "drop") %>%  #.groups = "drop" no more grouping
  complete(ANCOMcat, Family, fill = list(ASVsInFamByHab = 0)) #add in zero if not found in that ANCOM cat
# View(numberFamsANCOM_I6S)
print(arrange(numberFamiliesANCOM_I6S, -ASVsInFamByHab), n= nrow(numberFamsANCOM_I6S))
# Get names of families that have at least 3 representatives in one family or the other
fams3_I6S <- numberFamsANCOM_I6S$Family[which(numberFamsANCOM_I6S$ASVsInFamByHab >=3)]

# Looks good!
# View(numberFamsANCOM_I6S[numberFamsANCOM_I6S$Family %in% fams3_I6S,])
I6S_3fam_df <- numberFamsANCOM_I6S[numberFamsANCOM_I6S$Family %in% fams3_I6S,]
I6S_3fam_df <- I6S_3fam_df[-which(I6S_3fam_df$Family == "NA"),] #drop NA families

# Addphylum info to this plot by merging with subset of I6S_ANCOMall_.05pct_df that has phylum info! DO distinct because otherwise subset is v long and get many-to-many relationship
famPhylaI6S <- distinct(I6S_ANCOMall_.05pct_df[,colnames(I6S_ANCOMall_.05pct_df) %in% c("Family", "Phylum")])
I6S_3fam_df <- left_join(I6S_3fam_df, famPhylaI6S, by = "Family")
# View(I6S_3fam_df) # looks good!

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
    labels = levels(factor(I6S_2fam_df$ANCOMcat)),
    expand = expansion(mult = c(0.0, 0.0)),
    limits = c(0.2, 3)  # pull the two categories closer together
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "white", fill = NA, linewidth = 0.5),
    strip.placement = "outside",                 # put strips on the outer side (left)
    strip.text.y.right = element_text(angle = 0, hjust = 1, size =9, color= "black"),  # control orientation
    strip.background = element_blank(),
    panel.spacing.y = unit(0.1, "lines"),
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
quartz(width = 6, height = 7)
I6S_3famANCOM_bubPlot

I6S_3famANCOM_bubPlot #
# NOTE: LEGEND WILL NOT LEFT JUSTIFY NOR CAN I REMOVE Y-AXIS TICKS FROM THE RIGHT SIDE, SO I WILL
# MANUALLY DO THESE THINGS IN POWERPOINT


########################################################################
#   IV. FUNGI: EXPLORE ANCOM RESULTS (BIOAEROSOLS AND FOLIR SURFACE) 
########################################################################
ITS_leafAir_ANCOM_ASVlvl
colnames(ITS_leafAir_ANCOM_ASVlvl) 
# Make dataframe smaller
colsToKeepANCOM_ITS <- c("taxon", "lfc_sampleTypephyllosphere", "W_sampleTypephyllosphere", "q_sampleTypephyllosphere", "diff_sampleTypephyllosphere", "Phylum", "Class",
                         "Order","Family","Genus","Species")
ITS_leafAir_ANCOM_small <- ITS_leafAir_ANCOM_ASVlvl[,colnames(ITS_leafAir_ANCOM_ASVlvl) %in% colsToKeepANCOM_ITS]
colnames(ITS_leafAir_ANCOM_small)[colnames(ITS_leafAir_ANCOM_small) %in% "taxon"] <- "ASV_name" #change to ASV name for merging
head(ITS_leafAir_ANCOM_small)
# Add in category for bioaerosol or foliar surface
ITS_leafAir_ANCOM_small$ANCOMcat[which(ITS_leafAir_ANCOM_small$lfc_sampleTypephyllosphere > 0)] <- "foliar surface" #if lfc greater than 0, then foliar surface
ITS_leafAir_ANCOM_small$ANCOMcat[which(ITS_leafAir_ANCOM_small$lfc_sampleTypephyllosphere < 0)] <- "bioaerosol" #if lfc less than 0, then bioaerosol
# Double check
sort(ITS_leafAir_ANCOM_small$lfc_sampleTypephyllosphere[which(ITS_leafAir_ANCOM_small$ANCOMcat == "bioaerosol")]) #great, all are negative!
sort(ITS_leafAir_ANCOM_small$lfc_sampleTypephyllosphere[which(ITS_leafAir_ANCOM_small$ANCOMcat == "foliar surface")]) #great, all are positive!

length(which(ITS_leafAir_ANCOM_small$ANCOMcat == "bioaerosol")) #422
bioaerosolITSANCOMASVsName <- ITS_leafAir_ANCOM_small$ASV_name[which(ITS_leafAir_ANCOM_small$ANCOMcat == "bioaerosol")]
length(which(ITS_leafAir_ANCOM_small$ANCOMcat == "foliar surface")) #270
foliarITSANCOMASVsName <- ITS_leafAir_ANCOM_small$ASV_name[which(ITS_leafAir_ANCOM_small$ANCOMcat == "foliar surface")]
# Saved October 5, 2025

# WHAT ARE TOP ORDERS FOR BIOAEROSOL AND FOLIAR SURFACE ASVS? (ADDED TO MANUSCRIPT)
head(ITS_leafAir_ANCOM_small)
countsDistinctOrderITSANCOM <- ITS_leafAir_ANCOM_small %>%
  group_by(ANCOMcat, Order) %>%
  summarise(n = n_distinct(ASV_name), .groups = "drop") %>%
  group_by(ANCOMcat) %>%
  mutate(totalInANCOMcat = sum(n),
         prop = n / totalInANCOMcat,
         pct  = 100 * prop) %>% #pct is percentage
  ungroup()
# View(countsDistinctOrderITSANCOM)
# Do these add up to one? YES!
sum(countsDistinctOrderITSANCOM$prop[which(countsDistinctOrderITSANCOM$ANCOMcat == "bioaerosol")])
sum(countsDistinctOrderITSANCOM$prop[which(countsDistinctOrderITSANCOM$ANCOMcat == "foliar surface")])

# saveRDS(ITS_leafAir_ANCOM_small, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_leafAir_ANCOM_small.rds")

##### INVESTIGATE HOW ANCOM MODEL AND ACTUAL DATA COMPARE ######
# 1. Create ASV table with only ANCOM ASVs in air and leaf samples. Can use this to compare how ANCOM model and actual data compare
ITS_leafAir_ANCOM_ASVs <- ITS_leafAir_ANCOM_small$ASV_name
length(ITS_leafAir_ANCOM_ASVs) #421 different ASVs
# Pull out only ASV table for significantly different ASVs (ASVs are columns!)
ITS_airLeaf_sigASVtab <- ITSairPhylloOnly.ps_ASVs[,colnames(ITSairPhylloOnly.ps_ASVs) %in% ITS_leafAir_ANCOM_ASVs]
# View(ITS_airLeaf_sigASVtab)
# But which samples are air versus leaves?
#View(ITS_airLeaf_sigASVtab)
ITS_airLeaf_sigASVtab <- as.data.frame(ITS_airLeaf_sigASVtab) %>% 
  rownames_to_column(var= "sampleName")
dim(ITS_airLeaf_sigASVtab) #rows are samples, columns are ASVs
# View(ITS_airLeaf_sigASVtab)

# Add in air or leaf, depending on what sample type it is
ITS_airLeaf_sigASVtab <- ITS_airLeaf_sigASVtab %>%
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "air",
      str_detect(sampleName, "phyllo") ~ "phyllo",
      TRUE ~ NA_character_
    )
  )

head(ITS_airLeaf_sigASVtab)
cbind(ITS_airLeaf_sigASVtab$sampleName, ITS_airLeaf_sigASVtab$sampleType) #quick check to make sure code above worked... and it did!
# View(ITS_airLeaf_sigASVtab)
# One last check before calculating...
all.equal(ITS_leafAir_ANCOM_ASVs, colnames(ITS_airLeaf_sigASVtab)[2:(ncol(ITS_airLeaf_sigASVtab)-1)]) #okay these are still the significant ones

# 2. Add in number of reads per sample
# 1. Get dataframe with these from big, non-ANCOM set
readNumberITSairFoliar <- rowSums(ITSairPhylloOnly.ps_ASVs) #get read numbers from whole dataset
readNumberITSairFoliar_df <- as.data.frame(matrix(nrow=length(readNumberITSairFoliar), ncol=2))
colnames(readNumberITSairFoliar_df) <- c("sampleName", "TotalReads")
readNumberITSairFoliar_df[,1] <- names(readNumberITSairFoliar)
readNumberITSairFoliar_df[,2] <- unname(readNumberITSairFoliar)
# 2. Add to working dataframe
head(ITS_airLeaf_sigASVtab)
ITS_airLeaf_sigASVs <- merge(ITS_airLeaf_sigASVtab, readNumberITSairFoliar_df, by = "sampleName")
# View(ITS_airLeaf_sigASVs)

# 3. MAKE A "LONGER" DATAFRAME AND ADD IN TAXONOMY
head(ITS_airLeaf_sigASVs) 
# First, pivot longer so that ASVs are column 1 and Sample is second column 
colnames(ITS_airLeaf_sigASVs)
ITS_airLeaf_sigASVs_long <- ITS_airLeaf_sigASVs %>% 
  pivot_longer(cols = ASV_1:ASV_10516, names_to= "ASV_name", values_to = "ASV_abundance")
head(ITS_airLeaf_sigASVs_long)
# View(ITS_airLeaf_sigASVs_long)
dim(ITS_airLeaf_sigASVs_long) #116948      5
unique(ITS_airLeaf_sigASVs_long$sampleType) #great, only air and phyllo!
ITSairPhylloOnly_tax <- as.data.frame(as.matrix(tax_table(ITSairPhylloOnly.ps)))
ITSairPhylloOnly_tax <- ITSairPhylloOnly_tax %>% #make ASV_name column
  rownames_to_column(var = "ASV_name")
colnames(ITS_airLeaf_sigASVs_long)
# Result is "long" object with the abundance of each ANCOM ASV in each air and phyllo sample, plus taxonomy!
ITS_airLeaf_ANCOM_long <- left_join(ITS_airLeaf_sigASVs_long, ITSairPhylloOnly_tax, by = "ASV_name")
head(ITS_airLeaf_ANCOM_long)
# View(ITS_airLeaf_ANCOM_long)
# Add in ANCOM cat of each 
colnames(ITS_leafAir_ANCOM_small)
ITS_airLeaf_ANCOM_long <- merge(ITS_airLeaf_ANCOM_long, ITS_leafAir_ANCOM_small[,c(1,12)], by ="ASV_name")

# 4. Get mean prop of each ASV in bioaerosol and foliar
ITS_airLeaf_ANCOM_long <- ITS_airLeaf_ANCOM_long %>% 
  mutate(propASVsamp = ASV_abundance/TotalReads) #get each ASVs proportion in each sample. Double checked that this worked!
# View(ITS_airLeaf_ANCOM_long)
ITS_airLeaf_ANCOM_long2 <- ITS_airLeaf_ANCOM_long %>% 
  group_by(ASV_name, ANCOMcat, sampleType) %>% 
  summarize(meanASVPctSampType = mean(propASVsamp)*100)
# View(ITS_airLeaf_ANCOM_long2)

# Make this wider for direct comparison
colnames(ITS_airLeaf_ANCOM_long2)
ITS_airLeaf_ANCOM_wide <- ITS_airLeaf_ANCOM_long2 %>% 
  pivot_wider(names_from= sampleType, values_from = meanASVPctSampType)
head(ITS_airLeaf_ANCOM_wide)
# View(ITS_airLeaf_ANCOM_wide)

# When is foliar greater than air and vice versa?
ITS_airLeaf_ANCOM_wide2 <- ITS_airLeaf_ANCOM_wide %>%
  mutate(whichGreater = case_when(
    air  > phyllo ~ "air",
    phyllo > air  ~ "foliar",
    air == phyllo ~ "ties",
    TRUE          ~ NA_character_     # any NA comparisons land here
  ))
# View(ITS_airLeaf_ANCOM_wide2)
ITSbioFoliarBiggerIndex <- intersect(which(ITS_airLeaf_ANCOM_wide2$ANCOMcat == "bioaerosol"), which(ITS_airLeaf_ANCOM_wide2$whichGreater == "foliar"))
length(ITSbioFoliarBiggerIndex)/422 # 3.7% of the time, foliar more abundant in a biaoerosol indicator
ITSfoliarAirBiggerIndex <- intersect(which(ITS_airLeaf_ANCOM_wide2$ANCOMcat == "foliar surface"), which(ITS_airLeaf_ANCOM_wide2$whichGreater == "air"))
length(ITSfoliarAirBiggerIndex)/270 #3.3% of the time, bioaerosol more abundant in bioaerosol indicator

# Add this to plot
ITS_airLeaf_ANCOM_wide2$discrepencyANCOM <- NA
ITS_airLeaf_ANCOM_wide2$discrepencyANCOM[ITSbioFoliarBiggerIndex] <- "foliar more abundant in bioaerosol indicator!"
ITS_airLeaf_ANCOM_wide2$discrepencyANCOM[ITSfoliarAirBiggerIndex] <- "bioaerosol more abundant in foliar indicator!"

######## ANCOM CHECKING PLOTS ########
head(ITS_airLeaf_ANCOM_long)
# View(ITS_airLeaf_ANCOM_long)
# GET THE PROPORTION OF READS IN EACH SAMPLE THAT COME FROM BIOAEROSOL AND FOLIAR ANCOM TAXA
ITS_airLeaf_ANCOM_longBySamp <- ITS_airLeaf_ANCOM_long %>% 
  group_by(sampleName, ANCOMcat) %>% 
  summarise(propSampANCOMcat = sum(propASVsamp)) #propASVsamp is proportion
# View(ITS_airLeaf_ANCOM_longBySamp)
colnames(ITS_airLeaf_ANCOM_longBySamp)
# Add back in sample type
# Add in air or leaf, depending on what sample type it is
ITS_airLeaf_ANCOM_longBySamp <- ITS_airLeaf_ANCOM_longBySamp %>%
  mutate(
    sampleType = case_when(
      str_detect(sampleName, "air") ~ "bioaerosol",
      str_detect(sampleName, "phyllo") ~ "foliar surface",
      TRUE ~ NA_character_
    )
  )
ITS_airLeaf_ANCOM_longBySamp
# Double check
# air_ITS_105 and phyllo_ITS_43 proportion of reads from two types of ANCOM taxa
ASVs_air105andPhyllo43_ITS <- ITSairPhylloOnly.ps_ASVs[which(rownames(ITSairPhylloOnly.ps_ASVs) %in% c("air_ITS_105", "phyllo_ITS_43")),]
totalReadsair105andPhyllo43_ITS <- rowSums(ASVs_air105andPhyllo43_ITS)
# Bioaerosol ASVs
bioaerosolANCOMASVs_air105andPhyllo43_ITS <- ASVs_air105andPhyllo43_ITS[,colnames(ASVs_air105andPhyllo43_ITS) %in% bioaerosolITSANCOMASVsName]
# Foliar surface ASVs
foliarANCOMASVs_air105andPhyllo43_ITS <- ASVs_air105andPhyllo43_ITS[,colnames(ASVs_air105andPhyllo43_ITS) %in% foliarITSANCOMASVsName]
propBioaerosol_air105andPhyllo43_ITS <- rowSums(bioaerosolANCOMASVs_air105andPhyllo43_ITS)/totalReadsair105andPhyllo43_ITS
propBioaerosol_air105andPhyllo43_ITS
propFoliar_air105andPhyllo43_ITS <- rowSums(foliarANCOMASVs_air105andPhyllo43_ITS)/totalReadsair105andPhyllo43_ITS
propFoliar_air105andPhyllo43_ITS
# Get indices
ITSair105bioIndex <- intersect(which(ITS_airLeaf_ANCOM_longBySamp$sampleName == "air_ITS_105"),which(ITS_airLeaf_ANCOM_longBySamp$ANCOMcat == "bioaerosol"))
ITSair105FoliarIndex <- intersect(which(ITS_airLeaf_ANCOM_longBySamp$sampleName == "air_ITS_105"),which(ITS_airLeaf_ANCOM_longBySamp$ANCOMcat == "foliar surface"))
ITSfoliar43bioIndex <- intersect(which(ITS_airLeaf_ANCOM_longBySamp$sampleName == "phyllo_ITS_43"),which(ITS_airLeaf_ANCOM_longBySamp$ANCOMcat == "bioaerosol"))
ITSfoliar43FoliarIndex <- intersect(which(ITS_airLeaf_ANCOM_longBySamp$sampleName == "phyllo_ITS_43"),which(ITS_airLeaf_ANCOM_longBySamp$ANCOMcat == "foliar surface"))

# These are correct, so calculation was correct
all.equal(ITS_airLeaf_ANCOM_longBySamp$propSampANCOMcat[ITSair105bioIndex],unname(propBioaerosol_air105andPhyllo43_ITS[2]))
all.equal(ITS_airLeaf_ANCOM_longBySamp$propSampANCOMcat[ITSair105FoliarIndex],unname(propFoliar_air105andPhyllo43_ITS[2]))
all.equal(ITS_airLeaf_ANCOM_longBySamp$propSampANCOMcat[ITSfoliar43bioIndex],unname(propBioaerosol_air105andPhyllo43_ITS[1]))
all.equal(ITS_airLeaf_ANCOM_longBySamp$propSampANCOMcat[ITSfoliar43FoliarIndex],unname(propFoliar_air105andPhyllo43_ITS[1]))

# 1. For each sample, what proportion of reads are bioaerosol versus foliar surface taxa?
colnames(ITS_airLeaf_ANCOM_longBySamp)

propANCOMByReadsPlot_ITS <- ggplot(ITS_airLeaf_ANCOM_longBySamp, aes(x = sampleType, y = propSampANCOMcat, fill=sampleType)) +
  geom_boxplot(outlier.shape = NA) +           
  scale_fill_manual(values=c("cornflowerblue", "forestgreen")) +
  geom_jitter(width = 0.15, size = 1.8, height =0) + #each sample point
  facet_wrap(~ ANCOMcat, nrow = 1)   +          # ANCOM categoey by sample
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(name= "Proportion bioaerosol or \nfoliar surface ANCOM taxa")
propANCOMByReadsPlot_ITS
# saveRDS(propANCOMByReadsPlot_ITS, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/propANCOMByReadsPlot_ITS_10-6-25.rds")

###################### PLOTS OF ANCOM DATA (NO STRUCTURAL ZEROS) ###################
ITS_leafAir_ANCOM_small <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/ITS_leafAir_ANCOM_small.rds")

##### LOG FOLD CHANGE PLOT, USING LFC ESTIMATIONS FROM ANCOMBC2 ####
ITS_leafAir_ANCOM_small
length(unique(which(ITS_leafAir_ANCOM_small$ANCOMcat == "foliar surface"))) #270
length(unique(which(ITS_leafAir_ANCOM_small$ANCOMcat == "bioaerosol"))) #422
# Set up colors
length(unique(ITS_leafAir_ANCOM_small$Order))
table(ITS_leafAir_ANCOM_small$ANCOM, ITS_leafAir_ANCOM_small$Order)
# Can use these 20 with at least 8 differentially abundant
table(ITS_leafAir_ANCOM_small$Order)[which(unname(table(ITS_leafAir_ANCOM_small$Order)) >= 12)]

# Keep only orders with at least 5 ANCOM-identified ASVs (not including NAs, which will report)
ITS_leafAir_ANCOM_small$Class[which(ITS_leafAir_ANCOM_small$Order == "NA")]

ANCOMITSOrdersToKeep <- c("Agaricales", "Auriculariales", "Cantharellales", "Capnodiales", "Chaetothyriales", "Dothideales", "Hymenochaetales", "Pleosporales", "Polyporales",
                          "Russulales", "Trechisporales", "Tremellales")
length(ANCOMITSOrdersToKeep) #12 orders

# ALMOST MATCHES COLORS IN HABITATANALYSES_SEPT.R EXCEPT Clostridia AND SOME EXTRAS
# From https://github.com/Nowosad/rcartocolor (except for "#FE8F42", which I think will work)
display_carto_all(colorblind_friendly = TRUE)
display_carto_pal(12, "Safe")
carto_pal(12, "Safe")


setdiff(carto_pal(12, "Safe"), unname(colsForOrderesITS))

colsForOrdersITS_ANCOM <- c("Agaricales" = "#88CCEE","Auriculariales" = "#CC6677","Cantharellales"= "#DDCC77","Capnodiales" = "#117733","Chaetothyriales" = "#332288",
                            "Dothideales" = "#44AA99", "Hymenochaetales" =  "#999933","Pleosporales" = "#882255","Polyporales" =  "#661100", "Russulales" = "#FE8F42",
                            "Trechisporales" = "#6699CC",  "Tremellales" = "#AA4499")

# Get a smaller dataframe only for the orders I'll plot
ITS_leafAir_ANCOM_forPlot <- ITS_leafAir_ANCOM_small[which(ITS_leafAir_ANCOM_small$Order %in% ANCOMITSOrdersToKeep),]
unique(ITS_leafAir_ANCOM_forPlot$Order)
ITS_leafAir_ANCOM_forPlot$Order <- as.factor(ITS_leafAir_ANCOM_forPlot$Order) #make class a factor
# Make a new column that SWITCHES these so that positive is biaoerosol!!
ITS_leafAir_ANCOM_forPlot$lfc_sampleTypeBIOAEROSOL <- (ITS_leafAir_ANCOM_forPlot$lfc_sampleTypephyllosphere)*-1

# THE PLOT

colnames(ITS_leafAir_ANCOM_forPlot)
head(ITS_leafAir_ANCOM_forPlot)
ITS_ANCOM_ASVsLog2FoldPlot <- ggplot(ITS_leafAir_ANCOM_forPlot, aes(x = reorder(ASV_name, -lfc_sampleTypeBIOAEROSOL), y = lfc_sampleTypeBIOAEROSOL, color = Order)) +
  geom_point(size = 3.5, shape = 21, color = "black", stroke = 0.3, aes(fill = Order)) +  # Add points
  scale_fill_manual(values= colsForOrdersITS_ANCOM) +
  theme_bw() +
  labs(title = "",
       x = NULL,
       y = expression(Log[e] ~ "Fold Change (bioaerosol \nrelative to foliar surface)"),
       fill = "Bacterial Order") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Add a horizontal line at y=0
  theme(axis.title.x = element_blank(),  #remove x-axis title
        axis.text.x = element_blank(), #remove tick labels
        axis.ticks.x = element_blank(), #remove x ticks
        panel.grid.major.x = element_blank(),  # remove major vertical gridlines
        panel.grid.minor.x = element_blank(), # remove minor vertical gridlines)
        legend.text = element_text(size = 9),  # Increase legend text size
        legend.title = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  expand_limits(x = c(-4, length(ITS_leafAir_ANCOM_forPlot$ASV_name) + 4.5)) ##  Add padding to the x-axis limits to make it so points aren't cut off!

# For a horizontal legend:
ITS_ANCOM_ASVsLog2FoldPlot <- ITS_ANCOM_ASVsLog2FoldPlot +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  theme(
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 18),
    axis.title.y = element_text(
      margin = margin(t = 6, r = 8, b = 6, l = 8)  # t,r,b,l (before rotation)
    ),
    legend.position      = "bottom",
    legend.title.position= "top",
    # make keys smaller
    legend.key.height    = unit(0.30, "lines"),
    legend.key.width     = unit(0.60, "lines"),
    legend.key           = element_rect(fill = NA, colour = NA),
    # reduce gaps between rows/columns of keys
    legend.spacing.y     = unit(0.10, "lines"),
    legend.spacing.x     = unit(0.15, "lines"),
    legend.box.spacing   = unit(0.10, "lines"),
    legend.margin        = margin(0, 0, 0, 0),
    # tighten text
    legend.text          = element_text(size = 8, lineheight = 0.9),
    panel.grid           = element_blank()
  )
quartz(width = 7, height = 5)
ITS_ANCOM_ASVsLog2FoldPlot

# Saved October 6, 2025
# saveRDS(ITS_ANCOM_ASVsLog2FoldPlot, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_ANCOM_ASVsLog2FoldPlot_10-6-25.rds")
# quartz()
ITS_ANCOM_ASVsLog2FoldPlot

##### MAKE PIE CHART OF TAXA ####
# (NOTE THAT ASVS UNCLASSFIED AT CLASS LEVEL NOT SHOWN)
########## AIR BACTERIA ##########
# Make data
head(ITS_leafAir_ANCOM_forPlot)
ITS_ANCOM_justAir <- ITS_leafAir_ANCOM_forPlot[which(ITS_leafAir_ANCOM_forPlot$ANCOMcat == "bioaerosol"),]
unique(ITS_ANCOM_justAir$ANCOMcat)
# pull out only the orders with at least one ANCOM air ASV
ITSairPieDat <- table(ITS_ANCOM_justAir$Order)

# Pie Chart with colors for orders consistent with larger plot
ITSairPieDat_copy <- ITSairPieDat
# ensure it's a named numeric vector
ITSairPieDat_copy <- as.numeric(ITSairPieDat_copy); names(ITSairPieDat_copy) <- names(ITSairPieDat)
# drop zeros/negatives so colors stay aligned
keep <- ITSairPieDat_copy > 0 & is.finite(ITSairPieDat_copy)
ITSairPieDat_copy <- ITSairPieDat_copy[keep]
# align colors to slice names explicitly
cols <- colsForOrdersITS_ANCOM[names(ITSairPieDat_copy)]
# Plot it
pie(ITSairPieDat_copy,
    labels = names(ITSairPieDat_copy),
    col = cols)

# Pie chart without labels (but above chart ensures that they are as they should be)
pie(unname(ITSairPieDat_copy), labels = NA, col = cols)
sum(unname(ITSairPieDat_copy)) #358
length(which(ITS_leafAir_ANCOM_small$lfc_sampleTypephyllosphere < 0)) #422
# 358 out of 422 ASVs shown
# Saved October 6, 2025
# saveRDS(ITSairPieDat, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITSairPieDat_10-6-25.rds")

########## LEAF BACTERIA ##########
# Make data
head(ITS_leafAir_ANCOM_forPlot)
ITS_ANCOM_justFoliar <- ITS_leafAir_ANCOM_forPlot[which(ITS_leafAir_ANCOM_forPlot$ANCOMcat == "foliar surface"),]
unique(ITS_ANCOM_justFoliar$ANCOMcat) #foliar surface, yay!
# pull out only the orders with at least one ANCOM phyllosphere ASV

ITSfoliarPieDat <- table(ITS_ANCOM_justFoliar$Order)
sort(ITSfoliarPieDat) #remove Ktedonobacteria since it only has one ASV


# Pie Chart with colors for orders consistent with larger plot
ITSfoliarPieDat_copy <- ITSfoliarPieDat
# ensure it's a named numeric vector
ITSfoliarPieDat_copy <- as.numeric(ITSfoliarPieDat_copy); names(ITSfoliarPieDat_copy) <- names(ITSfoliarPieDat)
# drop zeros/negatives so colors stay aligned
keep <- ITSfoliarPieDat_copy > 0 & is.finite(ITSfoliarPieDat_copy)
ITSfoliarPieDat_copy <- ITSfoliarPieDat_copy[keep]
# align colors to slice names explicitly
cols <- colsForOrdersITS_ANCOM[names(ITSfoliarPieDat_copy)]
# Plot it
pie(ITSfoliarPieDat_copy,
    labels = names(ITSfoliarPieDat_copy),
    col = cols)

# Pie chart without labels (but above chart ensures that they are as they should be)
pie(unname(ITSfoliarPieDat_copy), labels = NA, col = cols)
sum(unname(ITSfoliarPieDat)) #168 ASVs in the air ANCOM category (goes next to plot)
length(which(ITS_leafAir_ANCOM_small$lfc_sampleTypephyllosphere > 0)) #170
# 168 out of 170 plotted

# ## WHAT ARE TOP CLASSES FOR AIR AND FOLIAR SURFACE SPECIALISTS (REPORTED IN PAPER)
# # View(ITS_leafAir_ANCOM_forPlot)
# ITS_MEAN_airLeafANCOM_OrderesByType <- ITS_leafAir_ANCOM_forPlot %>%
#   group_by(ANCOMcat) %>%
#   count(Order) %>%
#   arrange(desc(n))
# ITS_MEAN_airLeafANCOM_OrderesByType
# # View(ITS_MEAN_airLeafANCOM_OrderesByType)
# 
# # Get percentages of ASVs
# ITS_MEAN_airOrderes <- ITS_MEAN_airLeafANCOM_OrderesByType %>%
#   filter(ANCOMcat == "air") %>%
#   mutate(classPercentage = n/sum(n)*100)
# ITS_MEAN_airOrderes
# # Check
# 22/37 # looks correct!
# 
# ITS_MEAN_leafOrderes <- ITS_MEAN_airLeafANCOM_OrderesByType %>%
#   filter(ANCOMcat == "phyllo") %>%
#   mutate(classPercentage = n/sum(n)*100)
# ITS_MEAN_leafOrderes


####################################
##### ANOTHER WAY OF SHOWING THE DATA: FACETED DOT PLOT ####
# Get number of taxa in each phylum  
colnames(ITS_leafAir_ANCOM_small)
numberOrdersANCOM_ITS <- ITS_leafAir_ANCOM_small %>% 
  group_by(ANCOMcat, Order) %>% 
  summarise(ASVsInOrderByHab = n_distinct(ASV_name))
# View(numberOrdersANCOM_ITS)
print(arrange(numberOrdersANCOM_ITS, -ASVsInOrderByHab), n= nrow(numberOrdersANCOM_ITS))

# Get the top 8 Orderes for air and for leaf surfaces (do top 9 b/c 9 gets ten orders hence this)
topFungalOrders_10 <- numberOrdersANCOM_ITS %>%
  #filter(Order != NA) %>% #drop NA
  group_by(ANCOMcat) %>%
  top_n(n = 7, wt = ASVsInOrderByHab) %>% # wt is the "weighting variable for selecting the top n rows
  arrange(ANCOMcat, desc(ASVsInOrderByHab))
unique(topFungalOrders_10$Order) #10
topFungalOrders_10 <- topFungalOrders_10[-which(topFungalOrders_10$Order == "NA"),]

## For air: 
topFungalOrders_10$Order[which(topFungalOrders_10$ANCOMcat == "bioaerosol")] #8 for bioaerosols
topFungalOrders_10$Order[which(topFungalOrders_10$ANCOMcat == "foliar surface")] #8 for foliar

topFungalOrders_10
# get these top 10 Orderes into a vector 
topFungalOrders_10_vec <- as.character(unique(topFungalOrders_10$Order))

# Put counts back in original data frame
ITS_airLeafANCOM_OrderCount <- ITS_leafAir_ANCOM_small %>%
  left_join(numberOrdersANCOM_ITS, by = c("ANCOMcat", "Order"))
# View(ITS_airLeafANCOM_OrderCount)

# MAKE DOT/BUBBLE PLOTS -- TOP 10 Orderes IN AIR AND PHYLLO 
# Get colors
colsForOrdersITS_ANCOM 
setdiff(unique(topFungalOrders_10$Order),names(colsForOrdersITS_ANCOM)) 

colsForOrdersITS_ANCOM_bubble <- c(colsForOrdersITS_ANCOM, "Xylariales" =  "#888888")


# Re-Order this based on highest to lowest # ASVs in air then leaves
# Make Order a factor so it can be colored (Ordered from highest # ASVs to lowest in air, then leaves)
topFungalOrders_10$Order <- factor(topFungalOrders_10$Order, levels = topFungalOrders_10_vec)
unique(topFungalOrders_10$Order) #no NA
colnames(topFungalOrders_10)
ITS_ANCOM_bubblePlot <- ggplot(topFungalOrders_10,
                               aes(x = ANCOMcat, y = Order, size = ASVsInOrderByHab, fill = Order)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  facet_grid(Order ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0),
    axis.text.y = element_text(size = 10),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    # Legend layout (stack the two guides)
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.box.just = "left",
    
    # Make room so bubbles aren't clipped
    legend.key.height = unit(8, "mm"),
    legend.key.width  = unit(8, "mm"),
    legend.spacing.y  = unit(0.1, "lines"),
    legend.box.spacing= unit(0.1, "lines"),
    
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 10),
    
    # Center titles above each legend
    legend.title.position = "top",
    legend.title.align = 0.5,
    
    panel.grid = element_blank()
  ) +
  labs(
    x = "Bioaerosol or Foliar Surface Habitat",
    y = "Order",
    size = "Number of ASVs in Order"
  ) +
  # 1) FILL legend (first/top): bigger bubbles + thin outline IN THE LEGEND ONLY
  scale_fill_manual(
    values = colsForOrdersITS_ANCOM_bubble,
    guide = guide_legend(
      order = 1, byrow = TRUE, nrow = 4,
      override.aes = list(
        shape  = 21,
        size   = 3,      # <-- bigger legend bubbles
        stroke = 0.25,   # <-- thinner outline so colors show
        alpha  = 1
      )
    )
  ) +
  # 2) SIZE legend (second/bottom)
  scale_size(
    range = c(1, 10),
    guide = guide_legend(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(
        shape  = 21,
        stroke = 0.3,
        colour = "black"
      )
    )
  )
quartz()
ITS_ANCOM_bubblePlot

ITS_ANCOM_bubblePlot + theme(legend.position = "none") # exported as RobjectsSaved/ITS_ANCOM_bubblePlotNoLeg.pdf on August 7.

ITS_ANCOM_bubblePlot # exported as RobjectsSaved/ITS_ANCOM_bubblePlot.pdf on August 7.
# NOTE: LEGEND WILL NOT LEFT JUSTIFY NOR CAN I REMOVE Y-AXIS TICKS FROM THE RIGHT SIDE, SO I WILL
# MANUALLY DO THESE THINGS IN POWERPOINT




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
nrow(ANCOMspore) #173 had data
# Add in ANCOM category
ANCOMspore$ANCOMcat <- NA
ANCOMspore$ANCOMcat[which(ANCOMspore$lfc_sampleTypephyllosphere > 0)] <- "foliar surface"
ANCOMspore$ANCOMcat[which(ANCOMspore$lfc_sampleTypephyllosphere < 0)] <- "bioaerosol"

# Split these into separate air and phyllosphere dataframes
ANCOMspore_air <- ANCOMspore[which(ANCOMspore$ANCOMcat == "bioaerosol"),]
unique(ANCOMspore_air$ANCOMcat)
nrow(ANCOMspore_air) #80
ANCOMspore_foliar <- ANCOMspore[which(ANCOMspore$ANCOMcat == "foliar surface"),]
unique(ANCOMspore_foliar$ANCOMcat)
nrow(ANCOMspore_foliar) #93 with information

# INVESTIGATE PATTERNS AND MAKE PIE CHARTS -- tentatively supports hypothesis!
table(ANCOMspore$ANCOMcat, ANCOMspore$spore4cats)
#                no possible yes
# bioaerosol     47       22  11
# foliar surface 76        9   8

yesOrPossibleSporeIndex <- c(which(ANCOMspore$spore4cats == "yes"), which(ANCOMspore$spore4cats == "possible"))
sporeYesOrPossible <- ANCOMspore[yesOrPossibleSporeIndex,]
# View(sporeYesOrPossible)
unique(sporeYesOrPossible$Order) #14
unique(sporeYesOrPossible$Class) # "Actinobacteria", "Bacilli", "Gammaproteobacteria", "Alphaproteobacteria
# Some Ramlibacter (Alphaproteobacteria, produce cysts!!!)
# Confirmed too that some Pseudomonas in Madin database produce
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


# What percentage was there sporulation information for?
nrow(ANCOMspore)/nrow(I6S_leafAir_ANCOM_ASVlvl)*100 #41.09264, so 41.1%!!
length(unique(ANCOMspore$ASV_name))/nrow(I6S_leafAir_ANCOM_ASVlvl)*100 #matches above
# Previously: 46/(37+49)*100 #53.48837%

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
    y = "Number of differentially-\nabundant ASVs "  
  )

# October 3, 2025
# saveRDS(sporeAirPhyllo_barPlot_LegBott, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/sporeAirFoliar_barPlot_10-03-2025.rds")

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

sporesFisherTest <- fisher.test(contTableSporesCombinedReordered)
sporesFisherTest

# Fisher's Exact Test for Count Data -- MUCH MORE SIGNIFICANT
# 
# data:  contTableSporesCombinedReordered
# p-value = 0.001286
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.498123 6.681042
# sample estimates:
# odds ratio 
#   3.117245 


########################################################################
#   VI. BACTERIAL PIGMENTATION!
########################################################################
ANCOMgenPig
# Add in ANCOM category
ANCOMgenPig$ANCOMcat <- NA
ANCOMgenPig$ANCOMcat[which(ANCOMgenPig$lfc_sampleTypephyllosphere > 0)] <- "foliar surface"
ANCOMgenPig$ANCOMcat[which(ANCOMgenPig$lfc_sampleTypephyllosphere < 0)] <- "bioaerosol"

# Make stacked barplot
# Make dataframe
pigLeafAir_df <- as.data.frame(matrix(ncol=3, nrow=6))
colnames(pigLeafAir_df) <- c("ANCOMcat", "numberASVs", "PigConsensus")
table(ANCOMgenPig$ANCOMcat, ANCOMgenPig$PigConsensus)
#                 no yes yesAndNo
# bioaerosol      3  47       38
# foliar surface  5  36       37

# How many available?
length(unique(ANCOMgenPig$ASV_name))/length(unique(I6S_leafAir_ANCOM_small$ASV_name))*100 #39.42993, so 39% in MS

pigLeafAir_df[,1] <- c(rep("bioaerosol",3), rep("foliar surface", 3))
pigLeafAir_df
pigLeafAir_df[,3] <- c(rep(c("no", "possible", "yes"), 2))
pigLeafAir_df
table(ANCOMgenPig$ANCOMcat, ANCOMgenPig$PigConsensus)
pigLeafAir_df[,2] <- c(3, 38, 47, 5, 37, 36)
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
# October 6, 2025
# saveRDS(ANCOMgenPigPlot, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ANCOMgenPigPlot_10-06-2025.rds")

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

pigmentFisherTest <- fisher.test(pigContigencyTableCombinedReordered)
pigmentFisherTest

# Fisher's Exact Test for Count Data
# 
# data:  pigContigencyTableCombinedReordered
# p-value = 0.4766
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.3619491 12.8720344
# sample estimates:
# odds ratio 
#   1.932947 



############# SPORE SIZE NOW IN SCRIPT sporeSizeSept5.R ########


########################################################################
#   VI. FUNGAL FRUITING BODIES
########################################################################
##### RE DO THIS SECTION ##### 
# ITS_MEAN_airLeafANCOM
# # What are the top fungal orders by percentage of the ANCOM results (to be 
# # included in the paper with fruiting bodies)
# which(ITS_MEAN_airLeafANCOM$ANCOMcat == "air")
# unique(ITS_MEAN_airLeafANCOM$Order[which(ITS_MEAN_airLeafANCOM$ANCOMcat == "air")])
# 
# colnames(ITS_MEAN_airLeafANCOM)
# ITS_MEAN_airLeafANCOM_orderByType <- ITS_MEAN_airLeafANCOM %>% 
#   group_by(ANCOMcat) %>% 
#   count(Order) %>% 
#   arrange(desc(n)) %>% 
#   print(n= 33)
# ITS_MEAN_airLeafANCOM_orderByType %>% 
#   print(n= 33)
# 
# # This is the sum of all of the distinct ASVs (which matches other calculations)
# sum(ITS_MEAN_airLeafANCOM_orderByType$n[ITS_MEAN_airLeafANCOM_orderByType$ANCOMcat == "air"])
# 
# # IN AIR
# ITS_MEAN_airLeafANCOM_orderAir <- ITS_MEAN_airLeafANCOM_orderByType %>% 
#   filter(ANCOMcat == "air") %>% 
#   mutate(percentOrder = n/sum(ITS_MEAN_airLeafANCOM_orderByType$n[ITS_MEAN_airLeafANCOM_orderByType$ANCOMcat == "air"])*100)
# ITS_MEAN_airLeafANCOM_orderAir
# 
# # Check a few to confirm above is working... and yes!
# # Polyporales in air
# 85/150*100 #56.7%
# # Hymenochaetales in air
# 27/150*100 #18
# 
# # ON FOLIAR SURFACES
# ITS_MEAN_airLeafANCOM_orderFoliar <- ITS_MEAN_airLeafANCOM_orderByType %>% 
#   filter(ANCOMcat == "phyllo") %>% 
#   mutate(percentOrder = n/sum(ITS_MEAN_airLeafANCOM_orderByType$n[ITS_MEAN_airLeafANCOM_orderByType$ANCOMcat == "phyllo"])*100)
# ITS_MEAN_airLeafANCOM_orderFoliar


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
head(airFoliarFUNGuildTable1) #has species info as we want!
colnames(airFoliarFUNGuildTable1)
dim(airFoliarFUNGuildTable1)

#### HERE: SHOULD CHECK THAT SAME OUTPUT AS ON SERVER, i.e.
# Written to file October 4, 2025
# write.table(airFoliarFUNGuildTable1, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airFoliarFUNGuildTable1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# Code above should make it so that no quotation marks and rownames appear, but I confirmed this by opening in Excel anyway
# Note: double checked here: https://github.com/UMNFuN/FUNGuild/blob/master/README.md that the version I have, Guilds_v1.1.py, is the most current
# SEE THE FILE Desktop/CU_Research/SoilEdgeEffectsResearch/EdgeEffectsOnSoil/FUNGuildStepByStep/FUNGuildStepByStep for a more detailed guide to how I ran FUNGuild I ran FUNGuild
# Ran the following code in Terminal on personal MacBook:
# python3 ~/Desktop/CU_Research/SoilEdgeEffectsResearch/FUNGuild/Guilds_v1.1.py -otu ~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airFoliarFUNGuildTable1.txt -m -db fungi

# This created the output file in the same folder where the input ASV table was found:
# ~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airFoliarFUNGuildTable1.guilds_matched.txt !

# load results file:
rawFungResults <- read.delim(file= "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airFoliarFUNGuildTable1.guilds_matched.txt", header = T)
# View(rawFungResults) 
dim(rawFungResults)
#  6422 out of 9012 had results in fungGuild OTUs!!
# How many fungal taxa could we match?
6422/9012 *100 #71.26054, or 71.3% of ASVs had a match, but that's a lot more ASVs and foliar surfaces and air!!!
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
5710/9012*100 #63.35996% or 63.4% of all ANCOM-identified fungal ASVs had high probably or highly probable info!!
# OLD RESULTS: 3074/4673*100 #65.78% were matched

# INVESTIGATE with ANCOM results
ITS_leafAir_ANCOM_small
length(unique(ITS_leafAir_ANCOM_small$ASV_name)) #692
# Join with ANCOM results results with an inner join so that only ANCOM pulled out
ITS_leafAir_ANCOM_small
airLeafANCOM_FUNGUILD <- inner_join(ITS_leafAir_ANCOM_small, x=fungResults, by="ASV_name")
length(unique(airLeafANCOM_FUNGUILD$ASV_name)) #512
# This shows that 512/692 differentially abundant were matched!
512/629*100 #81.39905, or 81.4%
# OLD RESULTS:184/387*100 #47.5%

# Make a smaller dataset without specific sample abundances 
colnames(airLeafANCOM_FUNGUILD)
airLeafANCOM_FUNGUILD_Smaller <- airLeafANCOM_FUNGUILD %>% 
  select(ANCOMcat, lfc_sampleTypephyllosphere, diff_sampleTypephyllosphere, Trophic.Mode, Growth.Morphology, 
         Guild, Order, Family, Genus, Species)
# View(airLeafANCOM_FUNGUILD_Smaller)

length(which(airLeafANCOM_FUNGUILD_Smaller$Growth.Morphology == "NULL")) #28 were unassigned
dim(airLeafANCOM_FUNGUILD_Smaller[-which(airLeafANCOM_FUNGUILD_Smaller$Growth.Morphology == "NULL"),]) #484  10
484/629*100 #76.94754 or 76.9% had growth morphology data!!!
# OLD RESULTS: 121/387*100 #Only 31.27% had growth morphology data

# Remove those with no growth morph info
airLeafANCOM_FG_NoNullMorph <- airLeafANCOM_FUNGUILD_Smaller[-which(airLeafANCOM_FUNGUILD_Smaller$Growth.Morphology == "NULL"),]
dim(airLeafANCOM_FG_NoNullMorph)

# What proportion of the original air or phyllo species could be matched for growth morphology?
airMorphNumber <- length(which(airLeafANCOM_FG_NoNullMorph$ANCOMcat == "bioaerosol")) #number of ASVs enriched in air with morph data = 324
airANCOMNumber <- length(which(ITS_leafAir_ANCOM_small$ANCOMcat == "bioaerosol")) #422 original bioaerosol taxa
airMorphNumber/airANCOMNumber*100 #76.77725% or 76.8% of air enriched taxa has morphological data

foliarMorphNumber <- length(which(airLeafANCOM_FG_NoNullMorph$ANCOMcat == "foliar surface")) #number of ASVs enriched in foliar data with morph data = 23
foliarANCOMNumber <- length(which(ITS_leafAir_ANCOM_small$ANCOMcat == "foliar surface")) #number of ASVs enriched in foliar = 237
foliarMorphNumber/foliarANCOMNumber*100 #59.25926% or 59.26% of foliar-enriched taxa could be matched :)

MorphAirLeavesTable <- table(airLeafANCOM_FG_NoNullMorph$Growth.Morphology, airLeafANCOM_FG_NoNullMorph$ANCOMcat)
MorphAirLeaves_df <- as.data.frame(MorphAirLeavesTable)
colnames(MorphAirLeaves_df) <- c("morphology", "ANCOMcat", "freq")
# View(MorphAirLeaves_df)

# Add Unclassified to show number that could be matched for bioaerosols and foliar surface
# 1. Add a new level, "unclassified"
levels(MorphAirLeaves_df$morphology) <- c(levels(MorphAirLeaves_df$morphology), "Unclassified")
airANCOMNumber - airMorphNumber #98 could not be matched (98/422 could not be matched)
airMorphNumber/airANCOMNumber*100 #76.77725% could be classified (matches above)
# View(MorphAirLeaves_df)
nrow(MorphAirLeaves_df) 
# Make last row unclassified bioaerosol and add in number unclassified
MorphAirLeaves_df[nrow(MorphAirLeaves_df) +1,] <- c("Unclassified", "bioaerosol", airANCOMNumber - airMorphNumber)
foliarANCOMNumber - foliarMorphNumber #110 (out of 270 foliar surface ANCOM cat) could NOT be matched
foliarMorphNumber/foliarANCOMNumber*100 #59.25926 could be matched, matches above!
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
sum(MorphAir_only$freq[-nrow(MorphAir_only)]) # 324 matches which is correct (matches number of ANCOM taxa that 
# I was able to assign morphology information), Did -nrow(MorphAir_only) to subtract out Unclassified

# FOLIAR
MorphLEAF_only <- MorphAirLeaves_df %>% 
  filter(morphology != "Unclassified") %>% #filter out unclassified
  filter(ANCOMcat == "foliar surface") %>% 
  mutate(percentMorph = (freq/sum(freq))*100) %>% 
  arrange(desc(percentMorph))
MorphLEAF_only
# View(MorphLEAF_only)
sum(MorphLEAF_only$freq[-nrow(MorphLEAF_only)]) # 160 matches which is correct (matches number of ANCOM taxa that
# I was able to assign morphology information), Did -nrow(MorphAir_only) to subtract out Unclassified

# MAKE A PLOT THAT SHOWS ANCOM CATEGORY AND MORPHOLOGY
# 1. Get which morphologies to include
length(unique(MorphAirLeaves_df$morphology)) #21 unique, but this is too many
numMorphCat <- MorphAirLeaves_df %>% 
  group_by(morphology) %>% 
  summarize(n= sum(freq)) %>% 
  arrange(n)
print(numMorphCat, n= nrow(numMorphCat))
# View(numMorphCat)
# Will use only those morphologies with at least 5 representative ASVs
manyMorphs <- numMorphCat$morphology[which(numMorphCat$n >= 5)]
manyMorphs <- as.character(manyMorphs) #remove factor levels
# View(numMorphCat) #checked that these match!

# 2. Get a smaller dataframe just with these morphologies
MorphAirLeaves_5_df <- MorphAirLeaves_df #make a copy to edit or stupid factors will mess up
MorphAirLeaves_5_df$morphology <- as.character(MorphAirLeaves_5_df$morphology) #remove factor levels
is.factor(MorphAirLeaves_5_df$morphology) #good this is false
# Remove rarer morphologies
MorphAirLeaves_5_df <- MorphAirLeaves_5_df[which(MorphAirLeaves_5_df$morphology %in% manyMorphs),]
# Remove unclassified
MorphAirLeaves_5_df <- MorphAirLeaves_5_df[-which(MorphAirLeaves_5_df$morphology == "Unclassified"),]
# View(MorphAirLeaves_5_df)
rownames(MorphAirLeaves_5_df) <- NULL #remove pesky, unintentional rownames
setdiff(unique(MorphAirLeaves_5_df$morphology), manyMorphs) #no differences, yay!
# View(MorphAirLeaves_5_df)
length(unique(MorphAirLeaves_5_df$morphology)) 11

# 3. Set colors
carto_pal(12, "Safe")
morphsToPlot <- unique(MorphAirLeaves_5_df$morphology)
class(morphsToPlot); is.vector(morphsToPlot)
length(morphsToPlot)
colsForANCOMmorphs <- setNames(
  c("#88CCEE","#CC6677","#117733","#DDCC77","#332288",
    "#44AA99","#999933","#882255","#FE8F42","#661100",
    "#6699CC"),
  morphsToPlot[1:11]
)

# 4. Make the plot!
airLeafANCOM_FUNGUILD_5_barPlot <- ggplot(MorphAirLeaves_5_df) +
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

airLeafANCOM_FUNGUILD_5_barPlot 
# saveRDS(airLeafANCOM_FUNGUILD_5_barPlot, "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airLeafANCOM_FUNGUILD_5_barPlot_Oct6_2025.rds")
# CRUST FUNGI (Corticoid): https://www.crustfungi.com/html/sidebar/introduction.html#:~:text=Morphological%20Definition,or%20hydnoid%20hymenophore;%20and%20holobasidia.

# tiff(filename="airLeafANCOM_FUNGUILD_barPlotApril15.tiff",height=5600,width=5200,units="px",res=800,compression="lzw")
# airLeafANCOM_FUNGUILD_barPlot
# dev.off()

# What did some these specific ASVs look like?
# Capnodiales, the leaf specialists!
# View(airLeafANCOM_FG_NoNullMorph[airLeafANCOM_FG_NoNullMorph$Order=="Capnodiales",])
# Take a look at all matches for leaf surfaces
#View(airLeafANCOM_FG_NoNullMorph[airLeafANCOM_FG_NoNullMorph$ANCOMcat=="phyllo",])

#View(airLeafANCOM_FG_NoNullMorph[airLeafANCOM_FG_NoNullMorph$ANCOMcat=="air",])

