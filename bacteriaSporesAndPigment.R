# bacteriaSporesAndPigment.R
# October 13, 2025
# Earlier version was bacterialSporeFormers.R (July 2, 2024)

# Description: This script:
# 1) attempts to match bacterial and archaeal ASVs from SRS air and phyllosphere samples 
# to those in the Madin et al. 2020 dataset, in order to investigate the hypothesis that the air microbiome is
# enriched in bacterial spore-formers. Code was adapted from code written by Josep Ramoneda and Michael Hoffert.
# 2) Explores pigmentation in bacteria (by matching ANCOM taxa that differed between foliar surfaces and air to
# ijsem database by Barbaran et al. 2017 by genus)

# Important objects created in this script and their saved locations for quick access.
# "RobjectsSaved/airLeafANCOM_withSpore_Oct2_2025.RData" -- has ANCOM results, taxonomy, and sporulation information from Madin!

###### I. SET-UP OF SCRIPT ######
library(car)
library(ggplot2)
library(tidyr)
library(dplyr)

# i. SET WORKING DIRECTORY:
setwd("/data/winfreyc/SRS_aeromicrobiome_2022")
list.files() #list all of the files in this working directory

# ii. BRING IN IMPORTANT FILES:
# 1. # Dataframe made in "MATCH MY 16S SEQUENCES TO THOSE IN THE GTDB DATABASE" below (reading in now though since below is Python code run in Terminal)
# Note that the fasta file that was matched was made in bioinformatics_16S.R  April 11, 2023
raw97AllMatchedWithGTDB.gtdb <- read.csv("madinMatching/16S_GTDB_aligned_ssu_97.tsv", sep = '\t', header=FALSE)
# 2. GTDB metadata file... it's big!
GTDBmeta <- read.csv("/data/winfreyc/GTDB_databases/bac120_metadata_r207.tsv", sep= '\t', header=TRUE)
# 3. Read in formatted bacterial ANCOM results (object itself made October 1, 2025 in ancom_Sept.R on own computer)
I6S_airLeafANCOM <- readRDS(file="importsFromLaptop/I6S_ANCOMall_.05pct_df.rds") 
str(I6S_airLeafANCOM) # just has differentially abundant taxa with new, 0.05% threshold
length(unique(I6S_airLeafANCOM$ASV_name)) #220 ANCOM taxa

# 4. read in Madin data 
madinDat <- read.csv("/data/winfreyc/madinDataset/Metadata_Madin2020_new.csv", header=TRUE)
# Do a few checks in this original database of taxa that are supposed to be able to sporulate, but that seem a little suss to me..
unique(madinDat$sporulation[which(madinDat$order== "Pseudomonadales")]) #there are some with yes!
sporePseudoIndex <- intersect(which(madinDat$order== "Pseudomonadales"), which(madinDat$sporulation== "yes")) #there are some with yes!
# View(madinDat[sporePseudoIndex, ])
unique(madinDat$sporulation[which(madinDat$family== "Sphingomonadaceae")]) #there are some with yes!, so that is why a Sphingomonas is appearing as a possible spore-former


###### II. MATCH MY 16S SEQUENCES TO THOSE IN THE GTDB DATABASE, USING THE CODE BELOW ###### 
# (Based on instructions and code in vsearch_align_16S_genomes.txt from Josep. I originally did steps 
# i and ii below for my soil edge effects paper (so code for that repeated in madinMatchingCode_Dec18_2023 
# script)
# 
# i. INSTALL VSEARCH (location is /data/winfreyc/miniconda3/bin/vsearch)
# # Version is v2.15.2_linux_x86_64, 251.8GB RAM, 32 cores
# conda install -c bioconda vsearch
# 
# ii. USE THE "MAKEUDB_USEARCH" COMMAND TO CREATE A UDB DATABASE FILE FROM SEQUENCES IN THE GTDB METADATA FILE
# # (can grab it from file on microbe, as is done here). NOTE: since this code was run, "bac120_ssu_reps_r207.fna"
# is now ALSO in the folder winfreyc/GTDB_databases on microbe!
# vsearch --makeudb_usearch /data/ramonedaj/bac120_ssu_reps_r207.fna -output /data/winfreyc/madinDataset/bac120_ssu_reps_r207.udb
# 
# iii. ALIGN THE READS FROM MY DATASET TO THE DATABASE OF SSUS (AND PUT THIS INTO THE FOLDER CALLED MADINMATCHINGDEC2023 IN /DATA/WINFREYC/SRS_MAY_2021/ FOLDER):
# # FIRST THROW ON A SCREEN!!
# # Note: the --id 0.97 flag set the identity threshold at 97%, meaning only sequences that are at least 97% identical over the aligned region will be considered as matches.
# vsearch --usearch_global /data/winfreyc/SRS_aeromicrobiome_2022/all16S_repset.fasta --db /data/winfreyc/madinDataset/bac120_ssu_reps_r207.udb --strand both --notrunclabels --iddef 0 --id 0.97 --maxrejects 100 --maxaccepts 100 --blast6out /data/winfreyc/SRS_aeromicrobiome_2022/madinMatching/16S_GTDB_aligned_ssu_97.tsv --threads 16
#####################

###### III. BRING IN THE MATCHED ASVs MADE ABOVE AND RE-FORMAT FOR LATER MATCHING ###### 
# i. EXPLORE DATAFRAME MADE ABOVE
head(raw97AllMatchedWithGTDB.gtdb)
nrow(raw97AllMatchedWithGTDB.gtdb) #202814 rows (compare with 36,183 rows at 99% match) 
length(unique(raw97AllMatchedWithGTDB.gtdb$V2)) #there are 14,896 unique accession numbers in this, meaning ofc that 
# some of my ASVs matched to the same GTDB taxa. 

# ii. FIX COLUMN NAMES. COLUMN ONE IS THE ASV INFORMATION BASED ON MY REPSET FILE. V2 IS THE ACCESSION. V3 IS THE % SIMILARITY THAT THE GTDB
# ASV/taxon has with my ASVs, V4 is the length of this overlapping region. V5 is the number of mismatches. After that, the column names 
# correspond to those here: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
colnames(raw97AllMatchedWithGTDB.gtdb) <- c("ASVandTaxonomy", "accession", "percSimilarity", "overlapRegionLength", "numbMismatch", "gapopen",
                                            "qstart",  "qend", "sstart", "send", "evalue", "bitscore")
head(raw97AllMatchedWithGTDB.gtdb)

# iii. Get the metadata from the GTDB from GTDB meta, but only for the taxa matched in raw97AllMatchedWithGTDB.gtdb
# "ncbi_taxid" should be the column name that corresponds with the Madin data one
raw97AllMatchedWithGTDBwithMeta <- merge(raw97AllMatchedWithGTDB.gtdb, GTDBmeta, by= "accession", all.x=TRUE)
# View(raw97AllMatchedWithGTDBwithMeta) #202,814, same as matches above

# iv. TRIM DOWN FILE ABOVE TO ONLY BE ANCOM- SELECTED ASVs
I6S_airLeafANCOM
head(I6S_airLeafANCOM) # is the ASV table with info on occupancy, ANCOM category, resType (whether structural zero or main 0), 
# percent of each ASV in each sample, taxonomy
str(I6S_airLeafANCOM) #can use taxon or Species to match below

# Make ASV that appears at the beginning of raw97AllMatchedWithGTDBwithMeta$ASVandTaxonomy into its own column
colnames(raw97AllMatchedWithGTDBwithMeta)
raw97AllMatchedWithGTDBwithMeta[,2] #this has all ASV names and their big long taxonomy
raw97AllMatchedWithGTDBwithMeta2 <- tidyr::separate(raw97AllMatchedWithGTDBwithMeta, ASVandTaxonomy, into = c("ASV_name", "SILVAtaxonomy"), sep = " ")
head(raw97AllMatchedWithGTDBwithMeta2)
colnames(raw97AllMatchedWithGTDBwithMeta2)
# Pull out only those GTDB ASVs that are in the ANCOM list of names
matchedASVsToKeep <- which(raw97AllMatchedWithGTDBwithMeta2$ASV_name %in% I6S_airLeafANCOM$ASV_name == TRUE) #get index of ASV names to keep from big matched dataset
length(matchedASVsToKeep) #matched 7,709 ASVs (so that means multiple matches in many cases)
ANCOM_ASVs_WithGTDB <- raw97AllMatchedWithGTDBwithMeta2[matchedASVsToKeep,]
colnames(ANCOM_ASVs_WithGTDB)

# EXPLORE MATCHES
dim(ANCOM_ASVs_WithGTDB) # 7,709   122
length(which(ANCOM_ASVs_WithGTDB$ASV_name %in% I6S_airLeafANCOM$ASV_name == TRUE)) #7709, as expected
length(which(ANCOM_ASVs_WithGTDB$ASV_name %in% I6S_airLeafANCOM$ASV_name == TRUE)) - length(I6S_airLeafANCOM$ASV_name) # 7489 "extra matches" woah
length(I6S_airLeafANCOM$ASV_name) - length(unique(ANCOM_ASVs_WithGTDB$ASV_name)) #59/220 were unable to be matched. So 161/220 could be matched

# Which ASVs could not be matched?
notMatchedASVs <- I6S_airLeafANCOM$ASV_name[which(I6S_airLeafANCOM$ASV_name %in% unique(ANCOM_ASVs_WithGTDB$ASV_name)==FALSE)]
notMatchedASVs #59, as expected
# Check out these unmatched ASVs
View(I6S_airLeafANCOM[I6S_airLeafANCOM$ASV_name %in% notMatchedASVs,])

# Which ASVs could be matched?
length(unique(ANCOM_ASVs_WithGTDB$ASV_name)) #270 unique ASVs were able to be matched to the database (this used to be only 69!!!)
# The fact that there are more GTDB ASVs than names in I6S_airLeafANCOM$ASV_name (and more importantly, than in the 
# line above) implies that some ASVs matched to multiple GTDB entries.
matchedASVs <- I6S_airLeafANCOM$ASV_name[which(I6S_airLeafANCOM$ASV_name %in% unique(ANCOM_ASVs_WithGTDB$ASV_name)==TRUE)]
# Check out these matched ASVs
# View(I6S_airLeafANCOM[I6S_airLeafANCOM$ASV_name %in% matchedASVs,])

# v. FILTER MORE SO THAT I ONLY INCLUDE GENOMES THAT ARE AT LEAST 90% COMPLETE AND ONLY UP TO 7 MISMATCHES WITH OUR 16S RRNA GENE READS (~95% SEQUENCE SIMILARITY)
colnames(ANCOM_ASVs_WithGTDB)
ANCOM_ASVs_WithGTDB.90 <- subset(ANCOM_ASVs_WithGTDB, checkm_completeness >= 90) #90% complete or more!
nrow(ANCOM_ASVs_WithGTDB) - nrow(ANCOM_ASVs_WithGTDB.90) #removed 74 sequences
max(ANCOM_ASVs_WithGTDB.90$numbMismatch) #there are currently up to 7 mismatches (which is about equal to 95%), so I don't need to do anything else
# Calculation for 7 mismatches == 95% similarity. Thinking through getting 97% similarity: 
# 515 and 806 primers; 515 is 19 bp long and 806 is 20 bp long 
# So! 515 primer would go to 533 and 806 would go to 787 
# 254 * 0.97 = 246.38 or 8 mismatches 
# 254 - 7 = 247. 247/254 = 0.972
min(ANCOM_ASVs_WithGTDB.90$numbMismatch) #0 mismatches
unique(ANCOM_ASVs_WithGTDB.90$ASV_name) #161
length(unique(ANCOM_ASVs_WithGTDB.90$ASV_name)) #161 ASVs (so removed one I think!)
# What percentage of ASVs could I match? 
length(unique(ANCOM_ASVs_WithGTDB.90$ASV_name))/length(I6S_airLeafANCOM$ASV_name)*100 #73.18182, so able to match 73.2%
# Were there any Archaea?
which(grepl(ANCOM_ASVs_WithGTDB.90$SILVAtaxonomy, pattern= "Archaea")== TRUE) #no Archaea!
unique(ANCOM_ASVs_WithGTDB.90$SILVAtaxonomy) #only 152 unique here, even though 161 ASVs

# vi. HOW MANY MATCHES DOES EACH ASV HAVE?
ASVfreqTable <- table(ANCOM_ASVs_WithGTDB.90$ASV_name)
min(ASVfreqTable) #1 (duh)
max(ASVfreqTable) #200 One ASV had 200 matches :(
mean(ASVfreqTable) #47.42236
median(ASVfreqTable) #median is 19
length(which(ASVfreqTable == 1)) #17 ASVs had only 1 match
length(which(ASVfreqTable > 1)) #144 ASVs had more than 1 match
length(which(ASVfreqTable == 2)) #16 ASVs had 2 matches

# viii. FINALLY, TO MAKE GTDB EASIER TO WORK WITH, KEEP ONLY COLUMNS OF INTEREST AND SPLIT UP TAXONOMY COLUMN
colnames(ANCOM_ASVs_WithGTDB.90)
ANCOM_ASVs_WithGTDB.90_df <- ANCOM_ASVs_WithGTDB.90[,c(1:6,15,25:26, 29, 59, 86, 90, 91)]
# View(ANCOM_ASVs_WithGTDB.90_df)
ANCOM_ASVs_WithGTDB.90_df <- tidyr::separate(ANCOM_ASVs_WithGTDB.90_df, SILVAtaxonomy, into = c("Kingdom", "Phylum",
                                                                                                "Class", "Order", "Family", "Genus", "ASVnameAgain"), sep = ";")
# Warning message:
#   Expected 7 pieces. Missing pieces filled with `NA` in 97 rows [9, 196, 199, 724, 884, 885, 952, 956, 1027, 1028, 1029, 1030,
#                                                                  1053, 1059, 1615, 1990, 1991, 1992, 1993, 2696, ...]. 
colnames(ANCOM_ASVs_WithGTDB.90_df)
unique(ANCOM_ASVs_WithGTDB.90_df$ASVname == ANCOM_ASVs_WithGTDB.90_df$ASVnameAgain) #since this is always true, drop ASVnameAgain
ANCOM_ASVs_WithGTDB.90_df$ASVnameAgain <- NULL
# View(ANCOM_ASVs_WithGTDB.90_df)

###### IV. WHAT ARE GENOME SIZES OF THE TAX ABOVE? ######
dim(ANCOM_ASVs_WithGTDB.90_df) #7635   19
head(ANCOM_ASVs_WithGTDB.90_df)

# i. MERGE BELOW TO GET INFORMATION ON ASVs (ANCOM, etc.) FOR EACH ASV AND THEN CLEAN UP
head(I6S_airLeafANCOM)
head(ANCOM_ASVs_WithGTDB.90_df)

ANCOM_GTDBMadin_df <- merge(ANCOM_ASVs_WithGTDB.90_df, I6S_airLeafANCOM, by = "ASV_name", all.x = TRUE)
# View(ANCOM_GTDBMadin_df)
colnames(ANCOM_GTDBMadin_df)
# KEEP UP BY DROPPING REDUNDANT COLUMNS
unique(ANCOM_GTDBMadin_df$Kingdom.x == ANCOM_GTDBMadin_df$Kingdom.y) #true so drop
ANCOM_GTDBMadin_df$Kingdom.y <- NULL
unique(ANCOM_GTDBMadin_df$Phylum.x == ANCOM_GTDBMadin_df$Phylum.y) #true so drop
ANCOM_GTDBMadin_df$Phylum.y <- NULL
unique(ANCOM_GTDBMadin_df$Class.x == ANCOM_GTDBMadin_df$Class.y) #true so drop
ANCOM_GTDBMadin_df$Class.y <- NULL
unique(ANCOM_GTDBMadin_df$Order.x == ANCOM_GTDBMadin_df$Order.y) #true so drop
ANCOM_GTDBMadin_df$Order.y <- NULL
unique(ANCOM_GTDBMadin_df$Family.x == ANCOM_GTDBMadin_df$Family.y) #FALSE, which differ?
# View(ANCOM_GTDBMadin_df[which(ANCOM_GTDBMadin_df$Family.x %in% ANCOM_GTDBMadin_df$Family.y ==FALSE),])
colnames(ANCOM_GTDBMadin_df) #use to pull out Family.x and Family.y
ANCOM_GTDBMadin_df[which(ANCOM_GTDBMadin_df$Family.x %in% ANCOM_GTDBMadin_df$Family.y ==FALSE),colnames(ANCOM_GTDBMadin_df) %in% c("Family.x", "Family.y")]
# When I ran SILVA, it made this Acidobacteriaceae (Subgroup 1), but was Acidobacteriaceae in NCBI
# So these are basically the same, but will keep both for now
unique(ANCOM_GTDBMadin_df$Genus.x == ANCOM_GTDBMadin_df$Genus.y) #true and false! Which differ?
colnames(ANCOM_GTDBMadin_df) 
ANCOM_GTDBMadin_df[which(ANCOM_GTDBMadin_df$Genus.x %in% ANCOM_GTDBMadin_df$Genus.y ==FALSE),colnames(ANCOM_GTDBMadin_df) %in% c("Genus.x", "Genus.y")]
# This shows that they are basically the same, but that, in some cases, I6S_airLeafANCOM has more information
# likely because taxonomy is more updated. Keep both columns for now

#Re-name some columns! 
colnames(ANCOM_GTDBMadin_df)
colnames(ANCOM_GTDBMadin_df)[c(3:6)] <- c("Kingdom", "Phylum",
                                             "Class", "Order")
# Drop unneeded columns (a lot of GTDB information and the ASV table)
colnames(ANCOM_GTDBMadin_df)
ANCOM_GTDBMadin_df <- ANCOM_GTDBMadin_df[,c(1:8,11,12,14,18, 162:171)] #keep 
dim(ANCOM_GTDBMadin_df) #10848    22
# View(ANCOM_GTDBMadin_df)

# Save to the server 
# saveRDS(ANCOM_GTDBMadin_df, file="RobjectsSaved/ANCOM_GTDBMadin_df_Oct13_2025.RData") #saved October 13, 2025
# # ii. SINCE THERE ARE SOME ASVS THAT ARE MATCHED TO MULTIPLE GTDB ENTRIES, I'LL TAKE THE MEAN GENOME SIZE OF THESE.
#### GRAYED OUT SINCE WE DON'T USE THIS ANALYSIS IN PAPER ####
# meanGenomeSize <- ANCOM_GTDBMadin_df %>%
#   dplyr::group_by(ASV_name) %>%
#   dplyr::summarise(meanGenSize = mean(genome_size))
# # View(meanGenomeSize)
# dim(meanGenomeSize) 
# # Merge these back to keep the original information
# ANCOM_ASVs_GTDB_meanGenSize_df <- left_join(ANCOM_GTDBMadin_df, meanGenomeSize, by="ASV_name")
# head(ANCOM_ASVs_GTDB_meanGenSize_df)
# # Finally, make it so that I only have unique ASVs (so that some ASVs are not double counted):
# ancomGTDBgenMatches_final <- dplyr::distinct(ANCOM_ASVs_GTDB_meanGenSize_df, ASV_name, .keep_all = TRUE)
# # View(ancomGTDBgenMatches_final)
# nrow(ancomGTDBgenMatches_final) == nrow(meanGenomeSize)
# sort(colnames(ancomGTDBgenMatches_final))

# Save to the server 
# saveRDS(ancomGTDBgenMatches_final, file="RobjectsSaved/ancomGTDBgenMatches_final.RData") #saved July 3, 2024 

###### IV. FORMAT MADIN DATA, MATCH WITH THE GTDB NCBI TAXON IDS FROM ABOVE, AND THEN GET SPORE INFO FROM MATCHING ######
head(madinDat)
dim(madinDat) #154,739     36
# View(madinDat)

# i. GET SUBSET OF MADIN DATAFRAME THAT HAS SPORE-FORMING INFORMATION
# filter Madin data to only have those ASVs/taxa with sporulation information
sporeInfoIndex <- which(is.na(madinDat$sporulation) == FALSE) #get index for the rows that do NOT have NA in their sporulation column
length(sporeInfoIndex) #17,153
madinSporeInfo <- madinDat[which(is.na(madinDat$sporulation) == FALSE),] #pull out only non-NA rows for sporulation
dim(madinSporeInfo) # 17,153 taxa have sporulation information!
# View(madinSporeInfo)
head(madinSporeInfo)
colnames(madinSporeInfo)
# For ease below, keep only some of the columns
madinSporeInfoStep1 <- madinSporeInfo[,c(2:7,10,11, 31:36)]
head(madinSporeInfoStep1)
#View(madinSporeInfoStep1)
colnames(madinSporeInfoStep1)[1] <- "ncbi_taxid" #rename this column so that I can do merge below
# saveRDS(madinSporeInfoStep1, file= "RobjectsSaved/madinSporeInfoStep1.RData") #saved here July 3, 2024

# ii. FURTHER SUBSET MADIN SO THAT IT ONLY HAS THE ASVS THAT ARE IN MY "ANCOM" DATASET USING GTDB THING ABOVE
madinSporeInfoStep1$ncbi_taxid[which(madinSporeInfoStep1$ncbi_taxid %in% ANCOM_GTDBMadin_df$ncbi_taxid == TRUE)]
madinANCOMMatch <- merge(ANCOM_GTDBMadin_df, madinSporeInfoStep1, by= "ncbi_taxid")
head(madinANCOMMatch)
# View(madinANCOMMatch)
dim(madinANCOMMatch) #2522   35
length(unique(madinANCOMMatch$ASV_name)) #111 unique ASVs had either yes or no data... SOOOO 111/220 SO 50.5%!
unique(madinANCOMMatch$sporulation[which(is.na(madinANCOMMatch$sporulation) == FALSE)]) #double check, yes there are all 'yes' or 'no'
# View(madinANCOMMatch)
# What percentage of our ASVs had sporulation data in Madin?
length(unique(madinANCOMMatch$ASV_name))/length(I6S_airLeafANCOM$ASV_name) *100 #50.45455 So, about 50.5%!, or 111/220 ASVs (as I said above!)
#View(madinANCOMMatch)
# Is it ever the case that an ASV (with multiple matches) has yes and no?
colnames(madinANCOMMatch)
head(madinANCOMMatch)
ASVsporeMatches <- madinANCOMMatch %>%
  group_by(ASV_name, sporulation) %>%
  summarise(count = n(), .groups = 'drop')
print(ASVsporeMatches, n= nrow(ASVsporeMatches)) #yes, many do

# Add in a column that counts the number of unique spore states for each ASV
ASVsporeMatches_extra <- ASVsporeMatches %>% 
  group_by(ASV_name) %>% 
  mutate(uniqueSporeStates = n()) %>% 
  ungroup()
print(ASVsporeMatches_extra, n= nrow(ASVsporeMatches_extra))
# View(ASVsporeMatches_extra)
nrow(ASVsporeMatches_extra[which(ASVsporeMatches_extra$uniqueSporeStates == 2),])/2 #divide by two since if two rows if 2 states!
# This shows that 18 of these ASVs had matches to both 'yes' and 'no' taxa, so these will be counted as "possible"

# COUNT ASVs WITH YES AND NOS AS POSSIBLE AND FORMAT TO BE MERGED BACK IN 
head(ASVsporeMatches_extra)
possibleIndex <- which(ASVsporeMatches_extra$uniqueSporeStates > 1) #index of those with yeses and nos
ASVsPossible <- ASVsporeMatches_extra$ASV_name[which(ASVsporeMatches_extra$uniqueSporeStates > 1)] #names of ASVs with multiple spore types
ASVsporeMatches_extra$spore4cats <- ASVsporeMatches_extra$sporulation #make this new column match "sporulation"
# NOW write over the spore4cats designation for only the ASVs that had yeses and nos matches and make them "potential"
ASVsporeMatches_extra$spore4cats[possibleIndex] <- "possible"
# For matching purposes, below, get only unique rows (so no duplicate ASVs), considering only columns 1 (ASV name) and 5 (consensus sporulation, i.e., spore4cats)
colnames(ASVsporeMatches_extra)
ASVsSporUnique_df <- ASVsporeMatches_extra[,c(1,5)] %>% distinct(ASV_name, .keep_all = TRUE)
print(ASVsSporUnique_df, n=nrow(ASVsSporUnique_df)) #great, shows the 111 matches!
head(ASVsSporUnique_df)
table(ASVsSporUnique_df$spore4cats) 
# no possible      yes 
# 84       18        9 

# PUT THE SPORULATION INFORMATION INTO A NEW DATAFRAME
colnames(madinANCOMMatch)
airLeafANCOM_withSpore <- merge(I6S_airLeafANCOM, ASVsSporUnique_df, by = "ASV_name", all.x = TRUE)
# View(airLeafANCOM_withSpore)
table(airLeafANCOM_withSpore$ANCOMcat, airLeafANCOM_withSpore$spore4cats)
# no possible yes
# bioaerosol     41       13   7
# foliar surface 43        5   2

#### SAVE THIS DATAFRAME ####
# saveRDS(airLeafANCOM_withSpore, file = "RobjectsSaved/airLeafANCOM_withSpore_Oct13_2025.rds") #October 13, 2025

###### V. INVESTIGATE PIGMENTATION DATABASE ######
# October 28, 2024: Madin 2020 does not have pigmnentation data. However, one of the source databases for Madin,
# Barberán et al., (2017), does. Downloaded database from (https://doi.org/10.6084/m9 .figshare.4272392. Read into microbe.

# Note that the following commented out lines were just investigating the "fierer" source taxa in the 
# Madin et al. database, which are a subset of the Barbaran dataset, ijsem, ultimately used below.
# # How many of the madinDat matches were from the Fierer database? 
# madinFiererOnly <- madinDat[which(madinDat$data_source== "fierer"),] #Fierer is the Barbaran database with pigmentation
# nrow(madinFiererOnly) #2648 taxa here
# madinDat2 <- madinDat #make a copy to manipulate
# colnames(madinDat2)[2] <- "ncbi_taxid"
# # only 2648 entries ....
# # rename column 2 to be "ncbi_taxid" for merging below
# colnames(madinFiererOnly)[2] <- "ncbi_taxid"
# # SUBSET THESE TO JUST BE "ANCOM" taxa...
# madinFiererOnly$ncbi_taxid[which(madinFiererOnly$ncbi_taxid %in% ANCOM_GTDBMadin_df$ncbi_taxid == TRUE)] #just take a look at the matches
# madinFiererANCOMMatch <- merge(ANCOM_GTDBMadin_df, madinFiererOnly, by= "ncbi_taxid")
# head(madinFiererANCOMMatch)
# dim(madinFiererANCOMMatch) #475 of these were all from Madin!
# length(unique(madinFiererANCOMMatch$ASV_name)) #79 unique ANCOM ASVs
# # View(madinFiererANCOMMatch)
# 
# # Madin dat that match ANCOM (but all Madin)
# madinANCOMMatch <- merge(ANCOM_GTDBMadin_df, madinDat2, by= "ncbi_taxid")
# madinANCOMMatch

# Open the Barbaran/fierer database stuff here:
# (Next 40 or so lines from the README.txt downloaded with the pigmentation database):
# SOME CURATION STEPS IN R:

#read table
ijsem <-read.delim("IJSEM_pheno_db_v1.0.txt", sep="\t", header=T, check.names=F, fill=T,
                   na.strings=c("NA", "", "Not indicated", " Not indicated","not indicated", "Not Indicated", "n/a", "N/A", "Na", "Not given", "not given","Not given for yeasts", "not indicated, available in the online version", "Not indicated for yeasts", "Not Stated", "Not described for yeasts", "Not determined", "Not determined for yeasts"))

#simplify column names
colnames(ijsem)<-c("Habitat", "Year", "DOI", "rRNA16S", "GC", "Oxygen",
                   "Length", "Width", "Motility", "Spore", "MetabAssays", "Genus", "Species", "Strain", "pH_optimum", "pH_range", "Temp_optimum", "Temp_range", "Salt_optimum", "Salt_range", "Pigment", "Shape", "Aggregation", "FirstPage", "CultureCollection", "CarbonSubstrate", "Genome", "Gram", "Subhabitat", "Biolog")

#clean Habitat column
levels(ijsem$Habitat)[levels(ijsem$Habitat)=="freshwater (river, lake, pond)"]<-"freshwater"
levels(ijsem$Habitat)[levels(ijsem$Habitat)=="freshwater sediment (river, lake, pond"]<-"freshwater sediment"

#clean Oxygen column
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="aerobic"]<-"obligate aerobe"
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="anerobic"]<-"obligate anerobe"
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="microerophile"]<-"microaerophile"

#clean pH_optimum column
ijsem$pH_optimum<-as.character(ijsem$pH_optimum)
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
ijsem$pH_optimum < -sapply(ijsem$pH_optimum, simplify=T, function(x){mean(as.numeric(unlist(strsplit(x, split="-", fixed=T))))})
#remove pH values <0 and >10
ijsem$pH_optimum[ijsem$pH_optimum<0 | ijsem$pH_optimum>10]<-NA

#clean Temp_optimum column
ijsem$Temp_optimum<-as.character(ijsem$Temp_optimum)
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
ijsem$Temp_optimum<-sapply(ijsem$Temp_optimum, simplify=T, function(x){mean(as.numeric(unlist(strsplit(x, split="-", fixed=T))))})

#clean Salt_optimum column
ijsem$Salt_optimum<-as.character(ijsem$Salt_optimum)
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
ijsem$Salt_optimum<-sapply(ijsem$Salt_optimum, simplify=T, function(x){mean(as.numeric(unlist(strsplit(x, split="-", fixed=T))))})
#there are some formatting issues that should be solved

# View(ijsem)
colnames(ijsem)
colnames(madinANCOMMatch)

# SHOWS THAT THERE ARE NO OVERLAPS BETWEEN THESE DATA TYPES
# See if any of these match
unique(ijsem$rRNA16S %in% madinANCOMMatch$ncbi_taxid)
unique(ijsem$rRNA16S %in% madinANCOMMatch$accession)
unique(ijsem$rRNA16S %in% madinANCOMMatch$species_tax_id)
unique(ijsem$rRNA16S %in% madinANCOMMatch$rRNA16S_genes)
unique(ijsem$rRNA16S %in% madinANCOMMatch$rRNA16S_genes)

unique(ijsem$Genome%in% madinANCOMMatch$ncbi_taxid)
unique(ijsem$Genome%in% madinANCOMMatch$accession)
unique(ijsem$Genome%in% madinANCOMMatch$species_tax_id)
unique(ijsem$Genome%in% madinANCOMMatch$rRNA16S_genes)
unique(ijsem$Genome%in% madinANCOMMatch$rRNA16S_genes)

####### WILL TRY TO GET MY ANCOM TAXA AND MATCH BY SPECIES AND IF NOT SPECIES, GENUS #####
dim(ANCOM_GTDBMadin_df)
length(which(is.na(ANCOM_GTDBMadin_df$Genus.y)==TRUE)) #0
# None do not have genus information (Genus.y because this is the original SILVA taxonomy not the GTDB one)

# 1. Which ANCOM ASVs have matches at the genus level in the Barbaran data?
colnames(ijsem)
colnames(ANCOM_GTDBMadin_df)
length(which(ANCOM_GTDBMadin_df$Genus.y %in% ijsem$Genus == TRUE)) #4721!
ijsem$Pigment[which(is.na(ijsem$Pigment)==TRUE)] #isolate the NAs (this step just shows that one below will work)
# Get dataframe with only those ASVs with pigmentation information
ijsemWithPig <- ijsem[which(is.na(ijsem$Pigment)==FALSE),] 
# View(ijsemWithPig)
# Do all of these have Genus information?
unique(is.na(ijsemWithPig$Genus)) #all of these have genus information, so good to move to the next step!
unique(ijsemWithPig$Pigment) #contains yeses and nos!
str(ijsemWithPig)

# Summarize pigmentation by genus within the ijsemWithPig dataframe
str(ijsemWithPig)
ijsemByGenusPig <- ijsemWithPig %>%
  group_by(Genus) %>%
  summarize(PigConsensus = case_when(
    all(Pigment == "yes") ~ "yes",
    all(Pigment == "no") ~ "no",
    TRUE ~ "yesAndNo"
  )) %>%
  ungroup()
# View(ijsemByGenusPig)
unique(ijsemByGenusPig$PigConsensus) #"yes", "no", "yesAndNo", as expected

# Now, merge with inner_join to keep only those with ANCOM info. Crucially, use I6S_airLeafANCOM NOT
# Madin pre-matched dataframe, since there may be some ASVs that match my ANCOM ones in ijsem that
# were not included in Madin
# (unmatched rows in either input are not included in the result)
ANCOM_ijsemGenusPig <- inner_join(ijsemByGenusPig, I6S_airLeafANCOM, by ="Genus")
# View(ANCOM_ijsemGenusPig)
length(unique(ANCOM_ijsemGenusPig$ASV_name)) #93 unique ASVs! Better than the 40 from last time (and better than
# using jsut those fierer matches in madin)
dim(I6S_airLeafANCOM) #220
length(unique(ANCOM_ijsemGenusPig$ASV_name))/dim(I6S_airLeafANCOM)[1]*100 #42.27273, 42.3% (added to manuscript)
length(intersect(ANCOM_ijsemGenusPig$ASV_name, I6S_airLeafANCOM$ASV_name))/220
# Said another way, I matched 93/220 differentially abundant taxa
# compare with: OLD 40/86*100 #46.51163! Lower percentage but more taxa!!!

# OCTOBER 3, 2025
# saveRDS(ANCOM_ijsemGenusPig, "RobjectsSaved/ANCOM_ijsemGenusPig_10-14-2025.rds")