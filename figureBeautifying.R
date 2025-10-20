# Just Figures
# September 3, 2025


# This script takes all main figures for the manuscript and re-formats them together
# THINGS TO DO:
# 1. REMOVE ALLL GRIDLINES ON THESE PLOTS!!!!!!!!!!!!!!
# 2. Genome equivalents should be 16S COPY NUMBER?

library(gridExtra)
library(ggplot2); packageVersion("ggplot2") #‘3.5.1’
library(dplyr)
library(scales)
library(patchwork)
library(grid)

###############################
# SCRIPT SET UP 
###############################
##### Objects used throughout:
spacer_h <- unit(0.4, "in")  #gap for grid arrange



##### 1. IMPORT FIGURES MADE IN VARIOUS SCRIPTS #####
#### MAIN FIGURES ####
# FIGURE 2: "Figure 2.  Heat maps showing the relative abundance of the top 35 fungal orders (a)
# and bacterial phyla (b) in air samples."


# FIGURE 4: "Figure 3.  Associations among taxa in the air and the potential source environments
# of foliar surfaces and soil"
# From: I6S_sourceTracking_Sept.R and ITS_sourceTracking_Sept.R
I6S_allTypesAir_2panels <- readRDS(file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_allTypesAir_2panels_I6S_09-22-25")
I6S_allTypesAir_2panels2 <- I6S_allTypesAir_2panels +
  theme(panel.grid = element_blank()) + # remove ALL gridlines
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.text = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_discrete(labels=c(
    "bioaerosol\nmatrix", "bioaerosol\npatch",
    "foliar\nsurfaces", "soil")
  )

I6S_allTypesAir_2panels2 
ITS_allTypesAir_2panels <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_allTypesAir_2panels_ITS_09-28-25")
ITS_allTypesAir_2panels2 <- ITS_allTypesAir_2panels +
  theme(panel.grid = element_blank()) + # remove ALL gridlines
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.text = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_discrete(labels=c(
    "bioaerosol\nmatrix", "bioaerosol\npatch",
    "foliar\nsurfaces", "soil")
  )
ITS_allTypesAir_2panels2
# Venn Diagram
I6S_vennPlot <- readRDS(file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_vennPlot_09-26-2025.rds")
# 50.3% goes below 2869;  4.7% goes below 229; 25.9% goes below 528; 19.0% goes below 479
ITS_vennPlot <- readRDS(file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_vennPlot_09-29-2025.rds") #saved September 29, 2025
# 5.6% goes below 3,079; 0.8% goes below 169; 52.5% goes below 1056; 31.1% goes below 369

# A no-text Venn to fill in manually:
vennPlotNOTEXT <- readRDS(file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/vennPlotNOTEXT_09-29-2025.rds")

grid.arrange(
  ITS_allTypesAir_2panels2, 
  rectGrob(gp = gpar(col = NA, fill = NA), height = spacer_h),
  I6S_allTypesAir_2panels2,
  ncol = 1,
  heights = unit.c(unit(1, "null"), spacer_h, unit(1, "null"))
)

# FIGURE 4: "Figure 5. Comparison of traits of fungi between the foliar surface
# source environment and the air."
# a. 
ITS_2ordANCOM_bubPlot <- readRDS(file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_2ordANCOM_bubPlot_10-16-2025.rds")
#quartz()
ITS_2ordANCOM_bubPlot + theme(
  strip.text = element_text(size = 11),
  strip.text.y.left = element_text(angle = 0, hjust = 1, size =.5, color= "black"),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 10),
  axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
  legend.justification = "left" 
) +
  scale_x_discrete(labels=c(
  "bioaerosol", "foliar surfaces")
)

# b. Description: "volume distributions of fungal uninucleate sexual spores by ASV in the air and foliar surface samples"
# From sporeSizeAug26_2025.R
sexSpore_ViolinPlot <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/sexSpore_ViolinPlot_Oct16.rds")
sexSpore_ViolinPlot +
  theme(panel.grid = element_blank()) + # remove ALL gridlines
  scale_x_discrete(labels=c(
    "bioaerosol", "foliar surfaces")
  )

# c. # Morphology, made in ANCOM_traits_I6SandITS_.05pct.R
colsForANCOMmorphs <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/colsForANCOMmorphs.rds")
airLeafANCOM_FUNGUILD_5_barPlot <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/airLeafANCOM_FUNGUILD_2_barPlot_Oct16_2025.rds")
quartz()
airLeafANCOM_FUNGUILD_5_barPlot + 
  theme(panel.grid = element_blank(), # remove ALL gridlines
        axis.text.y = element_text(size = 10, color= "black"),
        axis.title.y = element_text(size = 10), color= "black",
        axis.text.x = element_text(size = 10, color = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size =10, color = "black"),
        legend.justification = "left" ,
        legend.title.position = "top")  +
  guides(fill = guide_legend(nrow = 6)) +
  scale_x_discrete(labels=c(
    "bioaerosol", "foliar surfaces")) +
  scale_fill_manual(
    values = colsForANCOMmorphs,
    labels = c(
      "Agaricoid" = "Agaricoid",
      "Corticoid" = "Corticoid",
      "Corticoid-Polyporoid" = "Corticoid-Polyporoid",
      "Dimorphic" = "Dimorphic",
      "Dimorphic-Facultative Yeast" = "Dimorphic-Fac. Yeast",
      "Dimorphic-Facultative Yeast-Microfungus-Tremelloid" =
        "Dimorphic-Fac. Yeast\nMicrofungus-Tremelloid",
      "Dimorphic-Microfungus" = "Dimorphic-Microfungus",
      "Dimorphic-Tremelloid" = "Dimorphic-Tremelloid",
      "Dimorphic-Yeast" = "Dimorphic-Yeast",
      "Microfungus" = "Microfungus",
      "Polyporoid" = "Polyporoid",
      "Yeast" = "Yeast"
    )
  )
# FIGURE 5: Comparison of traits of bacteria between the foliar surface

I6S_allClassANCOM_bubPlot <- readRDS(file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_allClassANCOM_bubPlot_10-16-2025.rds")
#quartz()
I6S_allClassANCOM_bubPlot +
  theme(
  strip.text = element_text(size = 11),
  strip.text.y.left = element_text(angle = 0, hjust = 1, size =.5, color= "black"),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 10),
  axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
  legend.justification = "left" 
) +
  scale_x_discrete(labels=c(
    "bioaerosol", "foliar surfaces")
  )

# d. Sporulation, made in ANCOM_traits_I6SandITS_.05pct.R
spore_barPlot <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/sporeAirFoliar_barPlot_10-15-2025.rds")
# quartz()
spore_barPlot + 
  theme(panel.grid = element_blank(), # remove ALL gridlines
        axis.text.y = element_text(size = 10, color= "black"),
        axis.title.y = element_text(size = 10), color= "black",
        axis.text.x = element_text(size = 10, color = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size =10, color = "black"),
        legend.justification = "center" ,
        legend.title.position = "top")  +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_discrete(labels=c(
    "bioaerosol", "foliar surfaces")) 

# e. Pigmentation, made in ANCOM_traits_I6SandITS_.05pct.R
ANCOMgenPigPlot <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ANCOMgenPigPlot_10-15-2025.rds")
quartz()
ANCOMgenPigPlot + 
  theme(panel.grid = element_blank(), # remove ALL gridlines
        axis.text.y = element_text(size = 10, color= "black"),
        axis.title.y = element_text(size = 10), color= "black",
        axis.text.x = element_text(size = 10, color = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size =10, color = "black"),
        legend.justification = "center" ,
        legend.title.position = "top")  +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_discrete(labels=c(
    "bioaerosol", "foliar surfaces"))




#### SUPPLEMENTARY FIGURES ####
# Fig S2. Non-metric multidimensional scaling (NMDS) ordinations based on Bray-Curtis dissimilarities among air, foliar surface, and soil samples for fungal (a) and bacterial (b) datasets
# I6S (from full16S_EDArarefied_Part2_Sept.R) Stress is still 0.1148366 , so 0.11 as currently on figure
I6S_bySampTypeOrd <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_bySampTypeOrd_09-17-2025.rds")
I6S_bySampTypeOrd +
  theme(panel.grid = element_blank()) # remove ALL gridlines

# ITS (from fullITS_EDArarefied_Part2_Sept.R). Stress is  0.09624302, so matches current figure
ITS_bySampTypeOrd <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_bySampTypeOrd_09-17-2025.rds")
ITS_bySampTypeOrd +
  theme(panel.grid = element_blank()) # remove ALL gridlines

# Presence/absence figure ????
# From I6S_sourceTracking_Sept.R
I6S_allTypesPresAbs_2panels <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_allTypesPresAbs_2panels_09-25-25.rds")
# Median number of ASVs for each group
# sampleType   medianNumASVs

#   1 air                   106.
# 2 phyllosphere          592.
# 3 soil                 1974 

# SD of number of ASVs for each group
# sampleType   sdNumASVs
#   1 air               82.8
# 2 phyllosphere     385. 
# 3 soil             353. 

# From I6S_sourceTracking_Sept.R
ITS_allTypesPresAbs_2panels <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_allTypesPresAbs_2panels_09-29-25")
# Median number of ASVs per group
# sampleType   medianNumASVs
# <chr>                <dbl>
# 1 air                    724
# 2 phyllosphere           836
# 3 soil                   287

# SD of number of ASVs per group
# sampleType   sdNumASVs
# <chr>            <dbl>
# 1 air             273. 
# 2 phyllosphere    185. 
# 3 soil            72.5

# Fig S4. Genome equivalents and NMDS
# FIGURE 3: "Figure 3. Box and whisker plots indicating the number of fungal (a) and bacterial 
# (b) genome equivalents detected per meter cubed of sampled air in each air sample, separated by air samples collected in the forested matrix and open patch"
# From: habitatAnalyses_Sept.R
# CHANGE Y AXIS TO GENE COPY NUMBER!!
allITS_qPCR_byHabitatPlot <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/allITS_qPCR_byHabitat_09-21-2025.rds")
allITS_qPCR_byHabitatPlot2 <- allITS_qPCR_byHabitatPlot +
  theme(panel.grid = element_blank()) + # remove ALL gridlines
  scale_y_log10(
    name = expression(atop("Fungal genome", "equivalents" ~ symbol("\xb7") ~ m^{-3})),
    breaks = c(1e3, 1e4, 1e5, 1e6),
    labels = parse(text = c("10^3", "10^4", "10^5", "10^6")),
    minor_breaks = NULL,
    limits = c(1e3, 1e6) 
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10)
  ) +
  scale_x_discrete(labels = c( #change x-axis labels
    forest = "forested matrix",
    savanna = "open patch"
  ))
allITS_qPCR_byHabitatPlot2

allI6S_qPCR_byHabitatPlot <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/allI6S_qPCR_byHabitat_09-21-2025.rds")
allI6S_qPCR_byHabitatPlot2 <- allI6S_qPCR_byHabitatPlot  +
  theme(panel.grid = element_blank()) + # remove ALL gridlines
  scale_y_log10(
    name = expression(atop("Bacterial genome", "equivalents" ~ symbol("\xb7") ~ m^{-3})),
    breaks = c(1e3, 1e4, 1e5, 1e6),
    labels = parse(text = c("10^3", "10^4", "10^5", "10^6")),
    minor_breaks = NULL,
    limits = c(1e3, 1e6) 
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10)
  ) +
  scale_x_discrete(labels = c( #change x-axis labels
    forest = "forested matrix",
    savanna = "open patch"
  ))
allI6S_qPCR_byHabitatPlot2

# Ordinations
I6S_habitatAirOrd <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_habitatAirOrd_10-17-25.rds")
I6S_habitatAirOrd2 <- I6S_habitatAirOrd +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    legend.position = "none"
  )
I6S_habitatAirOrd2
  
ITS_habitatAirOrd <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_habitatAirOrd_10-17-25.rds")
ITS_habitatAirOrd2 <- ITS_habitatAirOrd +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    legend.position = "none"
  )
ITS_habitatAirOrd2

grid.arrange(allITS_qPCR_byHabitatPlot2, allI6S_qPCR_byHabitatPlot2, ITS_habitatAirOrd2, I6S_habitatAirOrd2, nrow=2)
###################################
# MAKING MAIN FIGURES PRETTY! 
###################################
# REMOVE ALLL GRIDLINES ON THESE PLOTS!!!!!!!!!!!!!!
# REMOVE ALLL GRIDLINES ON THESE PLOTS!!!!!!!!!!!!!!
# REMOVE ALLL GRIDLINES ON THESE PLOTS!!!!!!!!!!!!!!
# REMOVE ALLL GRIDLINES ON THESE PLOTS!!!!!!!!!!!!!!
# BELOW FORMAT EXAMPLE FROM EDGE EFFECTS PAPER!!!!!
  # From EnvDataPlotsStats.R (on own computer), as a starting place to standardizing all figures
  # Horizontal color bar legend
  #quartz() # Horizontal color bar legend - note that on March 10, 2025, I edited this to add a black outline around the legend color bar
  #quartz() 
  canCoverEU_10_grid + guides(fill = guide_colourbar(
    title.position = "top",  # Moves title above color bar
    title.hjust = 0.5,       # Centers title
    ticks = FALSE, 
    label.position = "bottom",  # Moves labels below the bar
    nbin = 100, 
    title = "Canopy cover (%)",
    direction = "horizontal", 
    barwidth = 6, 
    barheight = 1.2,
    frame.colour = "black"  # Adds a black outline to the color bar
  )) + 
  theme_bw() +
  theme(legend.text = element_text(size = 8),  # Adjust legend text size
        legend.title = element_text(size = 8, face = "bold") # Adjust legend title size
  )  # Moves the title above the color bar

# MARCH 10, 2025 -- saving this to include color bar legends:
tiff(filename="manscript/Final_MS_Version2/Figures/Fig2_transectTurnoverCanCover_March8_2025.tiff",height=5600,width=5200,units="px",res=800,compression="lzw")
canCoverEU_10_grid + guides(fill = guide_colourbar(
  title.position = "top",  # Moves title above color bar
  title.hjust = 0.5,       # Centers title
  ticks = FALSE, 
  label.position = "bottom",  # Moves labels below the bar
  nbin = 100, 
  title = "Canopy cover (%)",
  direction = "horizontal", 
  barwidth = 6, 
  barheight = 1.2,
  frame.colour = "black"  # Adds a black outline to the color bar
)) + 
  theme_bw() +
  theme(legend.text = element_text(size = 8),  # Adjust legend text size
        legend.title = element_text(size = 8, face = "bold") # Adjust legend title size
  )  # Moves the title above the color bar
dev.off()


