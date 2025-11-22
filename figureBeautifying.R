# Just Figures
# September 3, 2025


# This script takes all main figures for the manuscript and re-formats them together

library(gridExtra)
library(ggplot2); packageVersion("ggplot2") #‘4.0.0’
library(dplyr)
library(scales)
library(patchwork)
library(grid)
library(cowplot)

###############################
# SCRIPT SET UP 
###############################
options(scipen = 999) #turn off awful scientific notation
##### Objects used throughout:
spacer_h <- unit(0.4, "in")  #gap for grid arrange

##### 1. IMPORT FIGURES MADE IN VARIOUS SCRIPTS #####
#### MAIN FIGURES ####
# FIGURE 2: "Figure 2.  Heat maps showing the relative abundance of the top 35 fungal orders (a)
# and bacterial phyla (b) in air samples."
# Made in habitatAnalyses_Sept.R
ITSOrder_heatmap_faceted_sqrt <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITSOrder_heatmap_faceted_sqrt-10-26-2025.rds")
I6SClass_heatmap_faceted_sqrt <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6SClass_heatmap_faceted_sqrt_10-26-2025.rds")
ITSOrder_heatmap_faceted_sqrt <- ITSOrder_heatmap_faceted_sqrt +
  theme(
    axis.text.y = element_text(size= 7, color="black"))

I6SClass_heatmap_faceted_sqrt <- I6SClass_heatmap_faceted_sqrt  +
  theme(
    axis.text.y = element_text(size= 7, color="black"))

# Have to re-arrange manually in the quartz window
quartz(width = 5.5, height = 8.6, family = "Helvetica")
(ITSOrder_heatmap_faceted_sqrt  / I6SClass_heatmap_faceted_sqrt) + plot_layout(heights = c(1, 1))


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
##############################
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
##############################

# b. Description: "volume distributions of fungal uninucleate sexual spores by ASV in the air and foliar surface samples"
# From sporeSizeAug26_2025.R
sexSpore_ViolinPlot <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/sexSpore_ViolinPlot_Oct16.rds")
sexSpore_ViolinPlot +
  theme(panel.grid = element_blank()) + # remove ALL gridlines
  scale_x_discrete(labels=c(
    "bioaerosol", "foliar surfaces")
  )

##############################
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
##############################
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
##############################
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
##############################
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
# Figure S2. Non-metric multidimensional scaling (NMDS) ordinations based on Bray-Curtis dissimilarities among bioaerosol, foliar surface, and soil
# Made in fullITS_EDArarefied_Part2_Sept.R
ITS_bySampTypeOrd <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_bySampTypeOrd_10-25-2025.rds")

ITS_bySampTypeOrd <- ITS_bySampTypeOrd +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
  )
# Made in full1TS_EDArarefied_Part2_Sept.R
I6S_bySampTypeOrd <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_bySampTypeOrd_10-25-2025.rds")
I6S_bySampTypeOrd <- I6S_bySampTypeOrd +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    legend.position = "none"
  )

# Convert plots to grobs
grobITS_bySampTypeOrd <- ggplotGrob(ITS_bySampTypeOrd)
grobI6S_bySampTypeOrd <- ggplotGrob(I6S_bySampTypeOrd)

# Ensure matching panel widths(so that, without legend, bacterial plot is not wider than fungal plot)
grobI6S_bySampTypeOrd$widths <- grobITS_bySampTypeOrd$widths

# quartz() #exported as S2_sampleTypeOrdinations.pdf 
grid.arrange(grobITS_bySampTypeOrd, grobI6S_bySampTypeOrd, ncol=2)

##############################
# Figure S3. Proportions of taxa in bioaerosols and the potential source environments of foliar surfaces and soil
# From I6S_sourceTracking_Sept.R
I6S_allTypesPresAbs_2panels <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_allTypesPresAbs_2panels_10-25-25.rds")
I6S_allTypesPresAbs_2panels <- 
  I6S_allTypesPresAbs_2panels + theme(panel.grid = element_blank()) + # remove ALL gridlines
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.text = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
    legend.position = "none" #remove legend
  ) +
  scale_x_discrete(labels=c(
    "bioaerosol\nmatrix", "bioaerosol\npatch",
    "foliar\nsurfaces", "soil")
  ) +
  ylab("Proportion of foliar surface\nor soil indicator taxa")
I6S_allTypesPresAbs_2panels
# From ITS_sourceTracking_Sept.R
ITS_allTypesPresAbs_2panels <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITS_allTypesPresAbs_2panels_10-25-25.rds")
ITS_allTypesPresAbs_2panels <- ITS_allTypesPresAbs_2panels + theme(panel.grid = element_blank()) + # remove ALL gridlines
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.text = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
    legend.position = "none" #remove legend
  ) +
  scale_x_discrete(labels=c(
    "bioaerosol\nmatrix", "bioaerosol\npatch",
    "foliar\nsurfaces", "soil")
  ) +
  ylab("Proportion of foliar surface\nor soil indicator taxa")
ITS_allTypesPresAbs_2panels
quartz()
grid.arrange(
  ITS_allTypesPresAbs_2panels, 
  rectGrob(gp = gpar(col = NA, fill = NA), height = spacer_h),
  I6S_allTypesPresAbs_2panels,
  ncol = 1,
  heights = unit.c(unit(1, "null"), spacer_h, unit(1, "null"))
)

##############################
# Fig S4. Gene copies and NMDS
# From: (dataframe) habitatAnalyses_Sept.R 
allI6S_qPCR <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allI6S_qPCR.rds")
max(allI6S_qPCR$I6Scopies)
min(allI6S_qPCR$I6Scopies)

allI6S_qPCR_byHabitat <- ggplot(
  data = allI6S_qPCR,
  aes(x = HabitatAir, y = I6Scopies, fill = HabitatAir)
) +
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  geom_boxplot(outlier.shape = NA) +
  # jittered points in black
  geom_jitter(
    color  = "black",
    size   = 2,
    alpha  = 0.9,
    width  = 0.25,
    height = 0
  ) +
  # log10 y axis 
  scale_y_continuous(
    trans  = "log10",
    name   = "16S copies (bacteria)",
    breaks = c(1e3, 1e4, 1e5, 1e6),
    limits = c(1e3, 1e6)
  ) +
  scale_fill_manual(
    values = c(
      forest  = "forestgreen",
      savanna = "goldenrod"
    )
  ) +
  scale_x_discrete(name = NULL,
                   labels = c( #change x-axis labels
                     forest = "forested matrix",
                     savanna = "open patch"
                   )) +
 theme(
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      legend.position = "none"
    ) 

allI6S_qPCR_byHabitat

## ITS ####
allITS_qPCR <- readRDS(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/allITS_qPCR.rds")
max(allITS_qPCR$ITScopies)
min(allITS_qPCR$ITScopies)
allITS_qPCR_byHabitat <- ggplot(
  data = allITS_qPCR,
  aes(x = HabitatAir, y = ITScopies, fill = HabitatAir)
) +
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  geom_boxplot(outlier.shape = NA) +
  # jittered points in black
  geom_jitter(
    color  = "black",
    size   = 2,
    alpha  = 0.9,
    width  = 0.25,
    height = 0
  ) +
  # log10 y axis 
  scale_y_continuous(
    trans  = "log10",
    name   = "ITS copies (fungi)",
    breaks = c(1e4, 1e5, 1e6, 1e7, 1e8),
    limits = c(1e4, 1e8)
  ) +
  scale_fill_manual(
    values = c(
      forest  = "forestgreen",
      savanna = "goldenrod"
    )
  ) +
  scale_x_discrete(name = NULL,
                   labels = c( #change x-axis labels
                     forest = "forested matrix",
                     savanna = "open patch"
                   )) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  ) 

allITS_qPCR_byHabitat


# Ordinations
I6S_habitatAirOrd <- readRDS("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6S_habitatAirOrd_10-17-25.rds")
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

quartz()
grid.arrange(allITS_qPCR_byHabitat, allI6S_qPCR_byHabitat, ITS_habitatAirOrd2, I6S_habitatAirOrd2, nrow=2)



