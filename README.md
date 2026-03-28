# BioaerosolSourcesTraits

Scripts in order of processing for bacteria and fungal datasets:

1. bioinformatics_16S.R -- processing of raw 16S rRNA sequences. Mostly follows tutorial found at https://github.com/fiererlab/dada2_fiererlab and on B. Callahan (author of dada2) website.
2. bioinformatics_final_ITS.R -- processing of raw ITS region sequences. Nearly identical pipeline to 16S above, but with some steps modified to reflect ITS primers and sequences

3. full16S_EDA_part1_Sept.R -- bacteria/archaea (16S) data cleaning, archaea and contaminant checking/removal, rarefying
4. fullITS_EDA_part1_Sept.R -- fungi (ITS) data cleaning and contaminant checking/removal, rarefying

5. full16S_EDArarefied_part2_Sept.R -- bacteria (16S) exploratory data analysis and sample type analysis and plots
6. fullITS_EDArarefied_Part2_Sept.R -- fungi (ITS) exploratory data analysis and sample type analysis and plots

7. habitatAnalyses_Sept.R -- bacterial and fungal analyses exploring differences in open patch and forested matrix, including compositional and biomass potential differences

8. I6S_sourceTracking_Sept.R -- analyses to determine origin of bacteria in bioaerosols
9. ITS_sourceTracking_Sept.R -- analyses to determine origin of bacteria in fungi

10. ANCOMBC_robust_foliarAir_pt0.05.R -- performs ANCOMBC-2 analyses to determine differentially abundant taxa between bioaerosols and foliar surfaces for both bacterial and fungal datasets.
11. bacteriaSporesAndPigment.R -- infers pigmentation and sporulation status for ANCOM-identified ASVs using Madin et al. 2020 and ijsem databases
12. ANCOM_traits_I6SandITS_.05pct.R -- stats for how foliar surface and biaoerosol enriched taxa differ in traits. Both bacterial traits, just morphology (inference from FUNGuild and stats) for fungi.
13. sporeSize.05pct.R -- infers spore size for fungi (using Aguilar- Trigueros et al. 2023 dataset) and associated stats.

14. weather_qPCRtable -- creates formatted meteorological and qPCR plots for supplement.
15. figureBeautifying.R --beautifies paper figures together
16. NCBIuploads.R -- formatting metadata for Bioproject and SRA

