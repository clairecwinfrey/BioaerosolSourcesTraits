# ITS BIOINFORMATICS PIPELINE FOR SRS AEROMICROBIOME PROJECT (air, phyllosphere, and soil samples run together)
# started March 13, 2023
# Claire Winfrey

# First, this script demultiplexes raw reads from the MiSeq (separately for air and phyllosphere
# since they were on different MiSeq runs). Next, for relevant 2021 soil samples, air samples, and 
# phyllosphere samples, I use dada2 to trim, denoise, dereplicate, and merge reads, and the assign 
# taxonomy (using UNITE database, version released October 5, 2021).
# Finally, this script formats the data for use in R.

# *NOTE* Some of the important diagnostic plots (learning error rates) were not able to be plotted 
# using our server, so I scp'd them to my personal computer and did all the checks before continuing.

# This script is based on a combination of 1) Fierer Lab DADA2 tutorial created by Angela 
# Oliverio and Hannah Holland-Moritz (Mar 2, 2020 version, accessed from: 
# https://github.com/fiererlab/dada2_fiererlab) AND 
# 2) Benjamin Callahan's DADA2 ITS pipeline tutorial (last accessed Jan. 19, 2022 at 
# https://benjjneb.github.io/dada2/ITS_workflow.html)

####################
#    SET UP 
####################

# Load DADA2 and required packages (versions are as of March 13, 2023)
library(dada2); packageVersion("dada2") # the dada2 pipeline, version ‘1.22.0’
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this ‘1.52.0’
library(dplyr); packageVersion("dplyr") # for manipulating data #dplyr is working ‘1.0.10’
library(tidyr); packageVersion("tidyr") # for creating the final graph at the end of the pipeline ‘1.2.0’
library(Hmisc); packageVersion("Hmisc") # for creating the final graph at the end of the pipeline ‘4.6.0’
library(ggplot2); packageVersion("ggplot2") # for creating the final graph at the end of the pipeline ‘3.3.5’
library(plotly); packageVersion("plotly") # enables creation of interactive graphs, especially helpful for quality plots ‘4.10.0’

# Once the packages are installed, you can check to make sure the auxillary
# software is working and set up some of the variables that you will need 
# along the way.

# Set up pathway to idemp (demultiplexing tool) and test
idemp <- "/usr/bin/idemp"
system2(idemp) # Check that idemp is in your path and you can run shell commands from R

# Set up pathway to cutadapt (primer trimming tool) and test
cutadapt <- "/usr/local/Python27/cutadapt" #for microbe server
system2(cutadapt, args = "--version") # Check by running shell command from R
# 03/13/23: version is 1.8.1

### SET UP FILE PATHS FOR DATA AND BARCODES ###
# 1. AIR ITS SAMPLES
# Set path to shared data folder and contents
airITSdata.fp <- "/data/shared/MiSeq/02.08.2023_Mixed_ITS_Claire_Havrilla_McNorvell"
# List all files in shared folder to check path
list.files(airITSdata.fp) #Data looks correct
# Set file paths for barcodes file, map file, and fastqs
airITSbarcode.fp <- file.path(airITSdata.fp, "airITS_mappingFileNoNs") #used sequences from #Barcodes_demultiplexing column in 
# Barcodes_DADA2_demultiplex_FiererLab_Feb2019.xlsx Excel sheet (Plates 3 & 8) and then reverse-complimented them
airITS_I1.fp <- file.path(airITSdata.fp, "Undetermined_S0_L001_I1_001.fastq.gz") 
airITS_R1.fp <- file.path(airITSdata.fp, "Undetermined_S0_L001_R1_001.fastq.gz") 
airITS_R2.fp <- file.path(airITSdata.fp, "Undetermined_S0_L001_R2_001.fastq.gz") 

# 2. PHYLLOSPHERE ITS SAMPLES
# Set path to shared data folder and contents
phylloITSdata.fp <- "/data/shared/MiSeq/9.28.2022_Mixed_ITS_Claire_Wilkins_Havrilla"
# List all files in shared folder to check path
list.files(phylloITSdata.fp) 
# Set file paths for barcodes file, map file, and fastqs
phylloITSbarcode.fp <- file.path(phylloITSdata.fp, "phylloITS_mappingFileNoNs") #used sequences from #Barcodes_demultiplexing column in 
# Barcodes_DADA2_demultiplex_FiererLab_Feb2019.xlsx Excel sheet and then reverse-complimented them
phylloITS_I1.fp <- file.path(phylloITSdata.fp, "Undetermined_S0_L001_I1_001.fastq.gz") 
phylloITS_R1.fp <- file.path(phylloITSdata.fp, "Undetermined_S0_L001_R1_001.fastq.gz") 
phylloITS_R2.fp <- file.path(phylloITSdata.fp, "Undetermined_S0_L001_R2_001.fastq.gz") 

# SET UP FILE PATHS FOR WHERE RESULTS WILL GO
project.fp <- "/data/winfreyc/SRS_aeromicrobiome_2022" #this is where all of the results will go! 
list.files(project.fp) #when first doing this script, this only had copy of two mapping files

# Set up file paths/names of sub directories to stay organized
#all of these are now in SRS_aeromicrobiome_2022 folder
## For these, the processing has to be done separately for air and phyllo samples, since they
## were on separate runs.
airITS_preprocess.fp <- file.path(project.fp, "airITS_01_preprocess") 
airITS_demultiplex.fp <- file.path(airITS_preprocess.fp, "airITS_demultiplexed")
airITS_filtN.fp <- file.path(airITS_preprocess.fp, "airITS_filtN") 
airITS_trimmed.fp <- file.path(airITS_preprocess.fp, "airITS_trimmed")
phylloITS_preprocess.fp <- file.path(project.fp, "phylloITS_01_preprocess") 
phylloITS_demultiplex.fp <- file.path(phylloITS_preprocess.fp, "phylloITS_demultiplexed")
phylloITS_filtN.fp <- file.path(phylloITS_preprocess.fp, "phylloITS_filtN") 
phylloITS_trimmed.fp <- file.path(phylloITS_preprocess.fp, "phylloITS_trimmed")

########################################################
#    # PRE-PROCESSING DATA FOR DADA2 - 
# DEMULTIPLEX, REMOVE SEQUENCES WITH Ns, CUTADAPT
########################################################
# This is done separately for air and phyllosphere samples, since they were run on different ITS runs.
###### AIR ITS SAMPLES ######
### 1. DEMULTIPLEX ###
# #### Call the demultiplexing script
# Demultiplexing splits your reads out into separate files based on the barcodes associated with each sample. 
airITS_Demultiplex <- paste("-b", airITSbarcode.fp, "-I1", airITS_I1.fp, "-R1", airITS_R1.fp, "-R2", airITS_R2.fp, "-o", airITS_demultiplex.fp) 
system2(idemp, args = airITS_Demultiplex)

# Look at output of demultiplexing
list.files(airITS_demultiplex.fp) 
sort(list.files(airITS_demultiplex.fp))
list.files(airITS_demultiplex.fp, pattern = "R1") #166 R1 files
list.files(airITS_demultiplex.fp, pattern = "R2") #166 R2 files
list.files(airITS_demultiplex.fp, pattern = "unsigned") #2 "unsigned file

### 2. CLEAN UP THE OUTPUT FROM IDEMP ###
# Change names of unassignable reads so they are not included in downstream processing
airITS_unassigned_1 <- paste0("mv", " ", airITS_demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_unsigned.fastq.gz",
                              " ", airITS_demultiplex.fp, "/Unassigned_reads1.fastq.gz") #this moves undetermined read 1 sequences into a new .gz file 
airITS_unassigned_2 <- paste0("mv", " ", airITS_demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_unsigned.fastq.gz", 
                              " ", airITS_demultiplex.fp, "/Unassigned_reads2.fastq.gz") #likewise, this moves undetermined R2 reads to their own folder
system(airITS_unassigned_1)
system(airITS_unassigned_2)

# Rename files - use gsub to get names in order!
airITS_R1_names <- gsub(paste0(airITS_demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_"), "", 
                        list.files(airITS_demultiplex.fp, pattern="R1", full.names = TRUE))
file.rename(list.files(airITS_demultiplex.fp, pattern="R1", full.names = TRUE), 
            paste0(airITS_demultiplex.fp, "/R1_", airITS_R1_names))

airITS_R2_names <- gsub(paste0(airITS_demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_"), "", 
                        list.files(airITS_demultiplex.fp, pattern="R2", full.names = TRUE))
file.rename(list.files(airITS_demultiplex.fp, pattern="R2", full.names = TRUE),
            paste0(airITS_demultiplex.fp, "/R2_", airITS_R2_names))

# Get full paths for all files and save them for downstream analyses
# Forward and reverse fastq filenames have format: (fnFs and frRs stand for "full name forward" and
# "full name reverse" )
airITS_fnFs <- sort(list.files(airITS_demultiplex.fp, pattern="R1_", full.names = TRUE))
airITS_fnRs <- sort(list.files(airITS_demultiplex.fp, pattern="R2_", full.names = TRUE))

### 3. PRE-FILTER TO REMOVE SEQUENCE READS WITH NS ###
# Ambiguous bases will make it hard for cutadapt to find short primer sequences in the reads.
# To solve this problem, we will remove sequences with ambiguous bases (Ns)

# Name the N-filtered files to put them in filtN/ subdirectory
airITS_fnFs.filtN <- file.path(airITS_preprocess.fp, "filtN", basename(airITS_fnFs))
airITS_fnRs.filtN <- file.path(airITS_preprocess.fp, "filtN", basename(airITS_fnRs))

# Filter Ns from reads and put them into the filtN directory (file path is "SRS_aeromicrobiome_2022/airITS_01_preprocess/filtN")
filterAndTrim(airITS_fnFs, airITS_fnFs.filtN, airITS_fnRs, airITS_fnRs.filtN, maxN = 0, multithread = TRUE) 

### 4. CUTADAPT ###
# The setting up of the primer sequences and the custom functions are the same for 
# the air and phyllosphere processing!
## i. Prepare the primers sequences and custom functions for analyzing the results from cutadapt ##
# Set up the primer sequences to pass along to cutadapt
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  ## ITS forward primer sequence (ITS1F)
REV <- "GCTGCGTTCTTCATCGATGC"  ## ITS reverse primer sequence (ITS2)

# Write a function that creates a list of all orientations of the primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# Save the primer orientations to pass to cutadapt
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Write a function that counts how many time primers appear in a sequence
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Save the reverse complements of the primers to variables
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

#  Create the cutadapt flags #
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 50") 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 50") 

## ii. Prepare the primers sequences and custom functions for analyzing the results from cutadapt ##
# Before running cutadapt, we will look at primer detection for the first sample, as a check. There may be some primers here, we will remove them below using cutadapt.
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = airITS_fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = airITS_fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = airITS_fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = airITS_fnRs.filtN[[1]]))
# #                 Forward Complement Reverse RevComp
# FWD.ForwardReads       0          0       0       0
# FWD.ReverseReads       0          0       0    3944
# REV.ForwardReads       0          0       0   15902
# REV.ReverseReads       0          0       0       0

## iii. Remove primers with cutadapt and assess the output ##

# Create directory to hold the output from cutadapt
if (!dir.exists(airITS_trimmed.fp)) dir.create(airITS_trimmed.fp)
airITS_fnFs.cut <- file.path(airITS_trimmed.fp, basename(airITS_fnFs))
airITS_fnRs.cut <- file.path(airITS_trimmed.fp, basename(airITS_fnRs))
# An example file path from the reverses here is:
# "/data/winfreyc/SRS_aeromicrobiome_2022/airITS_01_preprocess/airITS_trimmed/R2_air_ITS_PCR_NTC_2.fastq.gz

# Run Cutadapt
for (i in seq_along(airITS_fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", airITS_fnFs.cut[i], "-p", airITS_fnRs.cut[i], # output files
                             airITS_fnFs.filtN[i], airITS_fnRs.filtN[i])) # input files
}

# As a sanity check, we will check for primers in the first cutadapt-ed sample:
airITS_CutAdaptCheck <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = airITS_fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = airITS_fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = airITS_fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = airITS_fnRs.cut[[1]]))

airITS_CutAdaptCheck
####################################################################################

###### PHYLLOSPHERE ITS SAMPLES ######
### 1. DEMULTIPLEX PHYLLOSPHERE ITS samples ###
phylloITS_Demultiplex <- paste("-b", phylloITSbarcode.fp, "-I1", phylloITS_I1.fp, "-R1", phylloITS_R1.fp, "-R2", phylloITS_R2.fp, "-o", phylloITS_demultiplex.fp) 
system2(idemp, args = phylloITS_Demultiplex)

# Look at output of demultiplexing
list.files(phylloITS_demultiplex.fp) 
sort(list.files(phylloITS_demultiplex.fp))
list.files(phylloITS_demultiplex.fp) #140 files
list.files(phylloITS_demultiplex.fp, pattern = "unsigned") #2 "unsigned" files

### 2. CLEAN UP THE OUTPUT FROM IDEMP ###
# Change names of unassignable reads so they are not included in downstream processing
phylloITS_unassigned_1 <- paste0("mv", " ", phylloITS_demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_unsigned.fastq.gz",
                              " ", phylloITS_demultiplex.fp, "/Unassigned_reads1.fastq.gz") #this moves undetermined read 1 sequences into a new .gz file 
phylloITS_unassigned_2 <- paste0("mv", " ", phylloITS_demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_unsigned.fastq.gz", 
                              " ", phylloITS_demultiplex.fp, "/Unassigned_reads2.fastq.gz") #likewise, this moves undetermined R2 reads to their own folder
system(phylloITS_unassigned_1)
system(phylloITS_unassigned_2)

# Rename files - use gsub to get names in order!
phylloITS_R1_names <- gsub(paste0(phylloITS_demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_"), "", 
                        list.files(phylloITS_demultiplex.fp, pattern="R1", full.names = TRUE))
file.rename(list.files(phylloITS_demultiplex.fp, pattern="R1", full.names = TRUE), 
            paste0(phylloITS_demultiplex.fp, "/R1_", phylloITS_R1_names))

phylloITS_R2_names <- gsub(paste0(phylloITS_demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_"), "", 
                        list.files(phylloITS_demultiplex.fp, pattern="R2", full.names = TRUE))
file.rename(list.files(phylloITS_demultiplex.fp, pattern="R2", full.names = TRUE),
            paste0(phylloITS_demultiplex.fp, "/R2_", phylloITS_R2_names))

# Get full paths for all files and save them for downstream analyses
# Forward and reverse fastq filenames have format: (fnFs and frRs stand for "full name forward" and
# "full name reverse" )
phylloITS_fnFs <- sort(list.files(phylloITS_demultiplex.fp, pattern="R1_", full.names = TRUE))
phylloITS_fnRs <- sort(list.files(phylloITS_demultiplex.fp, pattern="R2_", full.names = TRUE))

#### 3. PRE-FILTER TO REMOVE SEQUENCE READS WITH NS ####
# Name the N-filtered files to put them in filtN/ subdirectory (so file path is "SRS_aeromicrobiome_2022/phylloITS_01_preprocess/filtN")
phylloITS_fnFs.filtN <- file.path(phylloITS_preprocess.fp, "filtN", basename(phylloITS_fnFs))
phylloITS_fnRs.filtN <- file.path(phylloITS_preprocess.fp, "filtN", basename(phylloITS_fnRs))

# Filter Ns from reads and put them into the filtN directory
filterAndTrim(phylloITS_fnFs, phylloITS_fnFs.filtN, phylloITS_fnRs, phylloITS_fnRs.filtN, maxN = 0, multithread = TRUE) 

#### 4. CUTADAPT ####
# (the set up is the same as in air samples, see above)
## ii. Prepare the primers sequences and custom functions for analyzing the results from cutadapt ####
# Before running cutadapt, we will look at primer detection for the first sample, as a check. There may be some primers here, we will remove them below using cutadapt.
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = phylloITS_fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = phylloITS_fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = phylloITS_fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = phylloITS_fnRs.filtN[[1]]))

# Forward Complement Reverse RevComp
# FWD.ForwardReads       0          0       0       0
# FWD.ReverseReads       0          0       0   15856
# REV.ForwardReads       0          0       0   33023
# REV.ReverseReads       0          0       0       0

## iii. Remove primers with cutadapt and assess the output ##

# Create directory to hold the output from cutadapt
if (!dir.exists(phylloITS_trimmed.fp)) dir.create(phylloITS_trimmed.fp)
phylloITS_fnFs.cut <- file.path(phylloITS_trimmed.fp, basename(phylloITS_fnFs))
phylloITS_fnRs.cut <- file.path(phylloITS_trimmed.fp, basename(phylloITS_fnRs))
# An example file path from the reverses here is:
# "/data/winfreyc/SRS_aeromicrobiome_2022/phylloITS_01_preprocess/phylloITS_trimmed/R2_phyllo_ITS_PCR_NTC_2.fastq.gz

# Run Cutadapt
for (i in seq_along(phylloITS_fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", phylloITS_fnFs.cut[i], "-p", phylloITS_fnRs.cut[i], # output files
                             phylloITS_fnFs.filtN[i], phylloITS_fnRs.filtN[i])) # input files
}

# As a sanity check, we will check for primers in the first cutadapt-ed sample:
# got all zeros!
phylloITS_CutAdaptCheck <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = phylloITS_fnFs.cut[[1]]), 
                                 FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = phylloITS_fnRs.cut[[1]]), 
                                 REV.ForwardReads = sapply(REV.orients, primerHits, fn = phylloITS_fnFs.cut[[1]]), 
                                 REV.ReverseReads = sapply(REV.orients, primerHits, fn = phylloITS_fnRs.cut[[1]]))

phylloITS_CutAdaptCheck #all zeros, showing this worked

######################################################
#### DADA2 PIPELINE ####
######################################################
#### PRE-PIPELINE. GET ORGANIZED ####
## Create different sub directories for each different kind of sample
# (already processed soil data will be brought in later, when using 
# subset of each sample type, separately, to infer error rates.)

airITS_filter.fp <- file.path(project.fp, "airITS__02_filter") 
airITS_table.fp <- file.path(project.fp, "airITS__03_tabletax") 
phylloITS_filter.fp <- file.path(project.fp, "phylloITS_02_filter") 
phylloITS_table.fp <- file.path(project.fp, "phylloITS_03_tabletax")

# Put filtered reads into separate sub-directories for big data workflow
###### AIR ######
dir.create(airITS_filter.fp)
airITS_subF.fp <- file.path(airITS_filter.fp, "preprocessed_F") 
airITS_subR.fp <- file.path(airITS_filter.fp, "preprocessed_R") 
dir.create(airITS_subF.fp)
dir.create(airITS_subR.fp)
### PHYLLO ###
dir.create(phylloITS_filter.fp)
phylloITS_subF.fp <- file.path(phylloITS_filter.fp, "preprocessed_F") 
phylloITS_subR.fp <- file.path(phylloITS_filter.fp, "preprocessed_R") 
dir.create(phylloITS_subF.fp)
dir.create(phylloITS_subR.fp)

# Move R1 and R2 from trimmed to separate forward/reverse sub-directories
###### AIR ######
airITS_fnFs.Q <- file.path(airITS_subF.fp,  basename(airITS_fnFs)) 
airITS_fnRs.Q <- file.path(airITS_subR.fp,  basename(airITS_fnRs))
file.rename(from = airITS_fnFs.cut, to = airITS_fnFs.Q)
file.rename(from = airITS_fnRs.cut, to = airITS_fnRs.Q)
### PHYLLO ###
phylloITS_fnFs.Q <- file.path(phylloITS_subF.fp,  basename(phylloITS_fnFs)) 
phylloITS_fnRs.Q <- file.path(phylloITS_subR.fp,  basename(phylloITS_fnRs))
file.rename(from = phylloITS_fnFs.cut, to = phylloITS_fnFs.Q)
file.rename(from = phylloITS_fnRs.cut, to = phylloITS_fnRs.Q)

# File parsing; create file names and make sure that forward and reverse files match
###### AIR ######
airITS_filtpathF <- file.path(airITS_subF.fp, "filtered") # files go into preprocessed_F/filtered/
airITS_filtpathR <- file.path(airITS_subR.fp, "filtered") # ...
airITS_fastqFs <- sort(list.files(airITS_subF.fp, pattern="fastq.gz"))
airITS_fastqRs <- sort(list.files(airITS_subR.fp, pattern="fastq.gz"))
if(length(airITS_fastqFs) != length(airITS_fastqRs)) stop("Forward and reverse files do not match.")
### PHYLLO ###
phylloITS_filtpathF <- file.path(phylloITS_subF.fp, "filtered") # files go into preprocessed_F/filtered/
phylloITS_filtpathR <- file.path(phylloITS_subR.fp, "filtered") # ...
phylloITS_fastqFs <- sort(list.files(phylloITS_subF.fp, pattern="fastq.gz"))
phylloITS_fastqRs <- sort(list.files(phylloITS_subR.fp, pattern="fastq.gz"))
if(length(phylloITS_fastqFs) != length(phylloITS_fastqRs)) stop("Forward and reverse files do not match.")

# ### 1. FILTER AND TRIM FOR QUALITY
###### AIR ######
# #### Inspect read quality profiles
# It's important to get a feel for the quality of the data that we are using. To do this, we will plot the quality of some of the samples.
# If the number of samples is 20 or less, plot them all, otherwise, just plot 20 randomly selected samples
if( length(airITS_fastqFs) <= 20) {
  plotQualityProfile(paste0(airITS_subF.fp, "/", airITS_fastqFs))
  plotQualityProfile(paste0(airITS_subR.fp, "/", airITS_fastqRs))
} else {
  rand_samples <- sample(size = 20, 1:length(airITS_fastqFs)) # grab 20 random samples to plot
  airITS_fwd_qual_plots <- plotQualityProfile(paste0(airITS_subF.fp, "/", airITS_fastqFs[rand_samples]))
  airITS_rev_qual_plots <- plotQualityProfile(paste0(airITS_subR.fp, "/", airITS_fastqRs[rand_samples]))
}

airITS_fwd_qual_plots 
airITS_rev_qual_plots

# Or, to make these quality plots interactive, just call the plots through plotly
ggplotly(airITS_fwd_qual_plots)
ggplotly(airITS_rev_qual_plots)

# write plots to disk
saveRDS(airITS_fwd_qual_plots, paste0(airITS_filter.fp, "/airITS_fwd_qual_plots.rds"))
saveRDS(airITS_rev_qual_plots, paste0(airITS_filter.fp, "/airITS_rev_qual_plots.rds"))

ggsave(plot = airITS_fwd_qual_plots, filename = paste0(airITS_filter.fp, "/airITS_fwd_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = airITS_rev_qual_plots, filename = paste0(airITS_filter.fp, "/airITS_rev_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")

# #### Filter the data

airITS_filt_out <- filterAndTrim(fwd=file.path(airITS_subF.fp, airITS_fastqFs), filt=file.path(airITS_filtpathF, airITS_fastqFs), #removed trunclen parameter for the ITS data
                                 rev=file.path(airITS_subR.fp, airITS_fastqRs), filt.rev=file.path(airITS_filtpathR, airITS_fastqRs),
                                 maxEE=c(2,2), truncQ=2, maxN=0, rm.phix=TRUE, minLen = 50, 
                                 compress=TRUE, verbose=TRUE, multithread=TRUE) #I added minLen = 50 to remove spurrious, short reads and I removed the trunclen parameter (as recommended in ITS DADA2 tutorial). Nothing else was changed from original tutorial

dim(airITS_filt_out) 

# summary of samples in airITS_filt_out by percentage
airITS_filt_out %>% 
  data.frame() %>% 
  mutate(Samples = rownames(.),
         percent_kept = 100*(reads.out/reads.in)) %>%
  select(Samples, everything())

# Plot the quality of the filtered fastq files.
# figure out which samples, if any, have been filtered out (out of the random ones)
airITS_remaining_samplesF <-  airITS_fastqFs[rand_samples][
  which(airITS_fastqFs[rand_samples] %in% list.files(airITS_filtpathF))] # keep only samples that haven't been filtered out
airITS_remaining_samplesR <-  airITS_fastqRs[rand_samples][
  which(airITS_fastqRs[rand_samples] %in% list.files(airITS_filtpathR))] # keep only samples that haven't been filtered out
airITS_fwd_qual_plots_filt <- plotQualityProfile(paste0(airITS_filtpathF, "/", airITS_remaining_samplesF))
airITS_rev_qual_plots_filt <- plotQualityProfile(paste0(airITS_filtpathR, "/", airITS_remaining_samplesR))

airITS_fwd_qual_plots_filt
airITS_rev_qual_plots_filt

# write plots to disk
saveRDS(airITS_fwd_qual_plots_filt, paste0(airITS_filter.fp, "/airITS_fwd_qual_plots_filt.rds"))
saveRDS(airITS_rev_qual_plots_filt, paste0(airITS_filter.fp, "/airITS_rev_qual_plots_filt.rds"))

ggsave(plot = airITS_fwd_qual_plots_filt, filename = paste0(airITS_filter.fp, "/airITS_fwd_qual_plots_filt.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = airITS_rev_qual_plots_filt, filename = paste0(airITS_filter.fp, "/airITS_rev_qual_plots_filt.png"), 
       width = 10, height = 10, dpi = "retina")

###### PHYLLO ######
# #### Inspect read quality profiles
# It's important to get a feel for the quality of the data that we are using. To do this, we will plot the quality of some of the samples.
# If the number of samples is 20 or less, plot them all, otherwise, just plot 20 randomly selected samples
if( length(phylloITS_fastqFs) <= 20) {
  plotQualityProfile(paste0(phylloITS_subF.fp, "/", phylloITS_fastqFs))
  plotQualityProfile(paste0(phylloITS_subR.fp, "/", phylloITS_fastqRs))
} else {
  rand_samples <- sample(size = 20, 1:length(phylloITS_fastqFs)) # grab 20 random samples to plot
  phylloITS_fwd_qual_plots <- plotQualityProfile(paste0(phylloITS_subF.fp, "/", phylloITS_fastqFs[rand_samples]))
  phylloITS_rev_qual_plots <- plotQualityProfile(paste0(phylloITS_subR.fp, "/", phylloITS_fastqRs[rand_samples]))
}

phylloITS_fwd_qual_plots 
phylloITS_rev_qual_plots

# Or, to make these quality plots interactive, just call the plots through plotly
ggplotly(phylloITS_fwd_qual_plots)
ggplotly(phylloITS_rev_qual_plots)

# write plots to disk
saveRDS(phylloITS_fwd_qual_plots, paste0(phylloITS_filter.fp, "/phylloITS_fwd_qual_plots.rds"))
saveRDS(phylloITS_rev_qual_plots, paste0(phylloITS_filter.fp, "/phylloITS_rev_qual_plots.rds"))

ggsave(plot = phylloITS_fwd_qual_plots, filename = paste0(phylloITS_filter.fp, "/phylloITS_fwd_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = phylloITS_rev_qual_plots, filename = paste0(phylloITS_filter.fp, "/phylloITS_rev_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")

# #### Filter the data

phylloITS_filt_out <- filterAndTrim(fwd=file.path(phylloITS_subF.fp, phylloITS_fastqFs), filt=file.path(phylloITS_filtpathF, phylloITS_fastqFs), #removed trunclen parameter for the ITS data
                                    rev=file.path(phylloITS_subR.fp, phylloITS_fastqRs), filt.rev=file.path(phylloITS_filtpathR, phylloITS_fastqRs),
                                    maxEE=c(2,2), truncQ=2, maxN=0, rm.phix=TRUE, minLen = 50, 
                                    compress=TRUE, verbose=TRUE, multithread=TRUE) #I added minLen = 50 to remove spurrious, short reads and I removed the trunclen parameter (as recommended in ITS DADA2 tutorial). Nothing else was changed from original tutorial

# summary of samples in phylloITS_filt_out by percentage
phylloITS_filt_out %>% 
  data.frame() %>% 
  mutate(Samples = rownames(.),
         percent_kept = 100*(reads.out/reads.in)) %>%
  select(Samples, everything())

# Plot the quality of the filtered fastq files.
# figure out which samples, if any, have been filtered out (out of the random ones)
phylloITS_remaining_samplesF <-  phylloITS_fastqFs[rand_samples][
  which(phylloITS_fastqFs[rand_samples] %in% list.files(phylloITS_filtpathF))] # keep only samples that haven't been filtered out
phylloITS_remaining_samplesR <-  phylloITS_fastqRs[rand_samples][
  which(phylloITS_fastqRs[rand_samples] %in% list.files(phylloITS_filtpathR))] # keep only samples that haven't been filtered out
phylloITS_fwd_qual_plots_filt <- plotQualityProfile(paste0(phylloITS_filtpathF, "/", phylloITS_remaining_samplesF))
phylloITS_rev_qual_plots_filt <- plotQualityProfile(paste0(phylloITS_filtpathR, "/", phylloITS_remaining_samplesR))

phylloITS_fwd_qual_plots_filt
phylloITS_rev_qual_plots_filt

# write plots to disk
saveRDS(phylloITS_fwd_qual_plots_filt, paste0(phylloITS_filter.fp, "/phylloITS_fwd_qual_plots_filt.rds"))
saveRDS(phylloITS_rev_qual_plots_filt, paste0(phylloITS_filter.fp, "/phylloITS_rev_qual_plots_filt.rds"))

ggsave(plot = phylloITS_fwd_qual_plots_filt, filename = paste0(phylloITS_filter.fp, "/phylloITS_fwd_qual_plots_filt.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = phylloITS_rev_qual_plots_filt, filename = paste0(phylloITS_filter.fp, "/phylloITS_rev_qual_plots_filt.png"), 
       width = 10, height = 10, dpi = "retina")

# ### 2. INFER sequence variants
# In this part of the pipeline dada2 will learn to distinguish error from biological 
# differences using a subset of our data as a training set. After it understands the 
# error rates, we will reduce the size of the dataset by combining all identical 
# sequence reads into "unique sequences". Then, using the dereplicated data and 
# error rates, dada2 will infer the sequence variants (OTUs) in our data. Finally, 
# we will merge the coresponding forward and reverse reads to create a list of the 
# fully denoised sequences and create a sequence table from the result.

# In this part of the pipeline, I infer sequence variants *separately* for air, phyllo,
# and soil, since they were run on different MiSeq runs. (although these soils have
# been processed before, I do so here again since I am only working with a subset of the 
# soils from last time)

##### SOILS ######
# first, need to bring in soils that were processed in 2022 (collected 2021)
soilProcessedFfiltered.fp <- "~/ITS_SRS_May_2021/02_filter/preprocessed_F/filtered"
soilProcessedRfiltered.fp <- "~/ITS_SRS_May_2021/02_filter/preprocessed_R/filtered"
list.files(soilProcessedFfiltered.fp)
list.files(soilProcessedRfiltered.fp)
# Get names of ExtBlank samples to add to processing
extBlanks <- list.files(soilProcessedRfiltered.fp)[grep(list.files(soilProcessedFfiltered.fp), pattern="ExtBlank_")]
extBlanksNames <- gsub(extBlanks, pattern = ".fastq.gz", replacement = "")
extBlanksNames <- gsub(extBlanksNames, pattern = "R2_", replacement = "")
extBlanksNames

# Get names of water blanks samples to add to processing
waterBlanks <- list.files(soilProcessedRfiltered.fp)[grep(list.files(soilProcessedFfiltered.fp), pattern="ExtControlWater_")]
waterBlanksNames <- gsub(waterBlanks, pattern = ".fastq.gz", replacement = "")
waterBlanksNames <- gsub(waterBlanksNames, pattern = "R2_", replacement = "")
waterBlanksNames

# Get names of PCR blanks samples to add to processing
PCRBlanks <- list.files(soilProcessedRfiltered.fp)[grep(list.files(soilProcessedFfiltered.fp), pattern="ITS_NTC_")]
PCRBlanksNames <- gsub(PCRBlanks, pattern = ".fastq.gz", replacement = "")
PCRBlanksNames <- gsub(PCRBlanksNames, pattern = "R2_", replacement = "")
PCRBlanksNames

# These sample names were verified in the script called soilSamplesOrganizing.R (on my personal computer).
# I then added the controls whose names I found above
soilSamples <- c("1", "10", "103", "104", "105", "106", "108", "110", "111", "112", "113", "114", "115", "116", "117", "118", "119", "12", "120", "121", "123", "124", "126", "127", "128", "129",
           "13", "130", "131", "132", "133", "134", "135", "136", "137", "140", "141", "142", "143", "144", "145", "147", "148", "15", "150", "152", "154", "155", "157", "158", "159", "16", 
           "160", "161", "162", "166", "168", "169", "17", "170", "171", "172", "174", "176", "177", "178", "179", "18", "180", "182", "183", "185", "186", "187", "188", "189", "19", "190",
           "191", "192", "194", "196", "197", "198", "199", "20", "200", "202", "203", "204", "209", "21", "211", "212", "213", "217", "218", "219", "22", "220", "221", "224", "225", "227",
           "228", "229", "23", "233", "239", "240", "28", "29", "30", "31", "33", "35", "36", "38", "4", "40", "41", "44", "45", "46", "47", "48", "49", "51", "52", "54", 
           "55", "56", "57", "58", "6", "60", "61", "63", "65", "66", "67", "7", "74", "75", "76", "80", "81", "85", "86", "88", "89", "9", "91", "92", "93", "95", 
           "96", "97", "98", "99", extBlanksNames, waterBlanksNames, PCRBlanksNames)

length(soilSamples)
setdiff(seq_along(1:243), soilSamples) #as a final check, these are the soils that SHOULD NOT be included. 

# Get just the subset of the soil samples found above to process:
# First, make a new folders
soilsfilteredPathF.fp <- file.path("~/SRS_aeromicrobiome_2022/soilsfilteredPathF")
if (!dir.exists(soilsfilteredPathF.fp)) dir.create(soilsfilteredPathF.fp)
soilsfilteredPathR.fp <- file.path("~/SRS_aeromicrobiome_2022/soilsfilteredPathR")
if (!dir.exists(soilsfilteredPathR.fp)) dir.create(soilsfilteredPathR.fp)

# Copy all of the R1 and R2 files for the relevant soil samples into the ~/ITS_SRS_May_2021/02_filter/preprocessed_F/filtered folder
for (i in 1:length(soilSamples)) {
  system(paste0('cp ', soilProcessedFfiltered.fp, '/R1_', soilSamples[i], ".fastq.gz ", soilsfilteredPathF.fp))
  system(paste0('cp ', soilProcessedRfiltered.fp, '/R2_', soilSamples[i], ".fastq.gz ", soilsfilteredPathR.fp))
}

# Check to make sure that all of the files are in the correct location:
list.files(soilsfilteredPathF.fp) #looks good!
list.files(soilsfilteredPathR.fp) 

# #### Housekeeping step - set up and verify the file names for the output:
soilITS_filtFs <- list.files(soilsfilteredPathF.fp, pattern="fastq.gz", full.names = TRUE)
length(soilITS_filtFs) #180
soilITS_filtRs <- list.files(soilsfilteredPathR.fp, pattern="fastq.gz", full.names = TRUE)
length(soilITS_filtRs) #180

# Sample names in order
soilITS_sample.names <- substring(basename(soilITS_filtFs), regexpr("_", basename(soilITS_filtFs)) + 1) # doesn't drop fastq.gz
soilITS_sample.names <- gsub(".fastq.gz", "", soilITS_sample.names)
soilITS_sample.namesR <- substring(basename(soilITS_filtRs), regexpr("_", basename(soilITS_filtRs)) + 1) # doesn't drop fastq.gz
soilITS_sample.namesR <- gsub(".fastq.gz", "", soilITS_sample.namesR)

# Double check
if(!identical(soilITS_sample.names, soilITS_sample.namesR)) stop("Forward and reverse files do not match.")
names(soilITS_filtFs) <- soilITS_sample.names
names(soilITS_filtRs) <- soilITS_sample.names

# #### Learn the error rates
set.seed(100) # set seed to ensure that randomized steps are replicatable 
# Learn forward error rates (Notes: randomize default is FALSE)
soilITS_errF <- learnErrors(soilITS_filtFs, nbases = 1e8, multithread = TRUE, randomize = TRUE) 
# Learn reverse error rates
soilITS_errR <- learnErrors(soilITS_filtRs, nbases = 1e8, multithread = TRUE, randomize = TRUE)

# #### Plot Error Rates
# We want to make sure that the machine learning algorithm is learning the error rates properly. In the plots below,
# the red line represents what we should expect the learned error rates to look like for each of the 16 possible base
# transitions (A->A, A->C, A->G, etc.) and the black line and grey dots represent what the observed error rates are. 
# If the black line and the red lines are very far off from each other, it may be a good idea to increase the 
# ```nbases``` parameter. This alows the machine learning algorthim to train on a larger portion of your data and may
# help improve the fit.
soilITS_errF_plot <- plotErrors(soilITS_errF, nominalQ = TRUE) 
soilITS_errR_plot <- plotErrors(soilITS_errR, nominalQ = TRUE)

soilITS_errF_plot
soilITS_errR_plot
#
# write to disk
saveRDS(soilITS_errF_plot, paste0(soilsfilteredPathF.fp, "/soilITS_errF_plot.rds"))
saveRDS(soilITS_errR_plot, paste0(soilsfilteredPathF.fp, "/soilITS_errR_plot.rds"))
# saved to microbe and evaluated on own computer (since microbe does not like ggplot :< )

# #### Dereplication, sequence inference, and merging of paired-end reads
# In this part of the pipeline, dada2 will make decisions about assigning sequences to ASVs (called "sequence inference"). There is a major parameter option in the core function dada() that changes how samples are handled during sequence inference. The parameter ```pool = ``` can be set to: ```pool = FALSE``` (default), ```pool = TRUE```, or ```pool = psuedo```. For details on parameter choice, please see below, and further information on this blogpost [http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/](http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/), and explanation on the dada2 tutorial [https://benjjneb.github.io/dada2/pool.html](https://benjjneb.github.io/dada2/pool.html).

# # SAMPLES POOLED 
# For complex communities when you want to preserve rare taxa
# same steps, not in loop

# Dereplicate forward reads
soilITS_derepF.p <- derepFastq(soilITS_filtFs)
names(soilITS_derepF.p) <- soilITS_sample.names
# Infer sequences for forward reads
soilITS_dadaF.p <- dada(soilITS_derepF.p, err = soilITS_errF, multithread = TRUE, pool = TRUE)
names(soilITS_dadaF.p) <- soilITS_sample.names

# Dereplicate reverse reads
soilITS_derepR.p <- derepFastq(soilITS_filtRs)
names(soilITS_derepR.p) <- soilITS_sample.names
# Infer sequences for reverse reads
soilITS_dadaR.p <- dada(soilITS_derepR.p, err = soilITS_errR, multithread = TRUE, pool = TRUE)
names(soilITS_dadaR.p) <- soilITS_sample.names

# Merge reads together
soilITS_mergers <- mergePairs(soilITS_dadaF.p, soilITS_derepF.p, soilITS_dadaR.p, soilITS_derepR.p)
head(soilITS_mergers)

# #### Construct sequence table
soilITS_seqtab <- makeSequenceTable(soilITS_mergers)

# Save table as an r data object file
dir.create("~/SRS_aeromicrobiome_2022/soilsfilteredPathF/soilITS_table")
saveRDS(soilITS_seqtab, paste0("~/SRS_aeromicrobiome_2022/soilsfilteredPathF/soilITS_table", "/soilITS_seqtab.rds"))

###################################################################################
##### AIR ######

# #### Housekeeping step - set up and verify the file names for the output:
airITS_filtFs <- list.files(airITS_filtpathF, pattern="fastq.gz", full.names = TRUE)
length(airITS_filtFs) #164
airITS_filtRs <- list.files(airITS_filtpathR, pattern="fastq.gz", full.names = TRUE)
length(airITS_filtRs) #164
# which ones are we missing? (should have 164)
airITS_mappingFileNoNs <- read.delim("SRS_aeromicrobiome_2022/airITS_mappingFileNoNs", header= FALSE)
sampleNamesMapFile_airITS <- airITS_mappingFileNoNs$V2
airITS_filtFsSampleNames <- airITS_filtFs
airITS_filtFsSampleNames <-gsub("/data/winfreyc/SRS_aeromicrobiome_2022/airITS__02_filter/preprocessed_F/filtered/R1_","", airITS_filtFsSampleNames)
airITS_filtFsSampleNames <-gsub(".fastq.gz","", airITS_filtFsSampleNames)
setdiff(sampleNamesMapFile_airITS, airITS_filtFsSampleNames) #this shows that air_ITS_21 (wash blank #2), air_ITS_PCR_NTC_1,
# and air_ITS_PCR_NTC_1_H12 were not found in the processed data. 

# Sample names in order
airITS_sample.names <- substring(basename(airITS_filtFs), regexpr("_", basename(airITS_filtFs)) + 1) # doesn't drop fastq.gz
airITS_sample.names <- gsub(".fastq.gz", "", airITS_sample.names)
airITS_sample.namesR <- substring(basename(airITS_filtRs), regexpr("_", basename(airITS_filtRs)) + 1) # doesn't drop fastq.gz
airITS_sample.namesR <- gsub(".fastq.gz", "", airITS_sample.namesR)

# Double check
if(!identical(airITS_sample.names, airITS_sample.namesR)) stop("Forward and reverse files do not match.")
names(airITS_filtFs) <- airITS_sample.names
names(airITS_filtRs) <- airITS_sample.names

# #### Learn the error rates
set.seed(100) # set seed to ensure that randomized steps are replicatable 

# Learn forward error rates (Notes: randomize default is FALSE)
airITS_errF <- learnErrors(airITS_filtFs, nbases = 1e8, multithread = TRUE, randomize = TRUE) 

# Learn reverse error rates
airITS_errR <- learnErrors(airITS_filtRs, nbases = 1e8, multithread = TRUE, randomize = TRUE)

# #### Plot Error Rates
# We want to make sure that the machine learning algorithm is learning the error rates properly. In the plots below, the red line represents what we should expect the learned error rates to look like for each of the 16 possible base transitions (A->A, A->C, A->G, etc.) and the black line and grey dots represent what the observed error rates are. If the black line and the red lines are very far off from each other, it may be a good idea to increase the ```nbases``` parameter. This alows the machine learning algorthim to train on a larger portion of your data and may help imporve the fit.

airITS_errF_plot <- plotErrors(airITS_errF, nominalQ = TRUE) 
airITS_errF_plot
# Error: Warning messages:
# 1: Transformation introduced infinite values in continuous y-axis 
# 2: Transformation introduced infinite values in continuous y-axis 

# Message from Ben Callahan (dada2 creator) on an online forum
# This isn't an error, just a message from the plotting function to let you know that there were some zero values 
# in the data plotted (which turn into infinities on the log-scale). That is completely expected, it results 
# from the fact that not every combination of error type (e.g. A->C) and quality score (e.g. 33) is observed 
# in your data which is normal.

airITS_errR_plot <- plotErrors(airITS_errR, nominalQ = TRUE)
airITS_errR_plot
#
# write to disk
saveRDS(airITS_errF_plot, paste0(airITS_filtpathF, "/airITS_errF_plot.rds"))
saveRDS(airITS_errR_plot, paste0(airITS_filtpathR, "/airITS_errR_plot.rds"))

# #### Dereplication, sequence inference, and merging of paired-end reads
# In this part of the pipeline, dada2 will make decisions about assigning sequences to ASVs (called "sequence inference"). There is a major parameter option in the core function dada() that changes how samples are handled during sequence inference. The parameter ```pool = ``` can be set to: ```pool = FALSE``` (default), ```pool = TRUE```, or ```pool = psuedo```. For details on parameter choice, please see below, and further information on this blogpost [http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/](http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/), and explanation on the dada2 tutorial [https://benjjneb.github.io/dada2/pool.html](https://benjjneb.github.io/dada2/pool.html).

# # SAMPLES POOLED 
# For complex communities when you want to preserve rare taxa

# Dereplicate forward reads
airITS_derepF.p <- derepFastq(airITS_filtFs)
names(airITS_derepF.p) <- airITS_sample.names
# Infer sequences for forward reads
airITS_dadaF.p <- dada(airITS_derepF.p, err = airITS_errF, multithread = TRUE, pool = TRUE)
names(airITS_dadaF.p) <- airITS_sample.names

# Dereplicate reverse reads
airITS_derepR.p <- derepFastq(airITS_filtRs)
names(airITS_derepR.p) <- airITS_sample.names
# Infer sequences for reverse reads
airITS_dadaR.p <- dada(airITS_derepR.p, err = airITS_errR, multithread = TRUE, pool = TRUE)
names(airITS_dadaR.p) <- airITS_sample.names

# Merge reads together
airITS_mergers <- mergePairs(airITS_dadaF.p, airITS_derepF.p, airITS_dadaR.p, airITS_derepR.p)
head(airITS_mergers)

# #### Construct sequence table
airITS_seqtab <- makeSequenceTable(airITS_mergers)

# Save table as an r data object file
dir.create(airITS_table.fp)
saveRDS(airITS_seqtab, paste0(airITS_table.fp, "/airITS_seqtab.rds"))

###################################################################################
##### PHYLLO ######

# #### Housekeeping step - set up and verify the file names for the output:
phylloITS_filtFs <- list.files(phylloITS_filtpathF, pattern="fastq.gz", full.names = TRUE)
length(phylloITS_filtFs) #68
phylloITS_filtRs <- list.files(phylloITS_filtpathR, pattern="fastq.gz", full.names = TRUE)
length(phylloITS_filtRs) #68
# which ones are we missing? (should have 164)
phylloITS_mappingFileNoNs <- read.delim("SRS_aeromicrobiome_2022/phylloITS_mappingFileNoNs", header= FALSE)
sampleNamesMapFile_phylloITS <- phylloITS_mappingFileNoNs$V2
phylloITS_filtFsSampleNames <- phylloITS_filtFs
phylloITS_filtFsSampleNames <-gsub("/data/winfreyc/SRS_aeromicrobiome_2022/phylloITS_02_filter/preprocessed_F/filtered/R1_","", phylloITS_filtFsSampleNames)
phylloITS_filtFsSampleNames <-gsub(".fastq.gz","", phylloITS_filtFsSampleNames)
setdiff(sampleNamesMapFile_phylloITS, phylloITS_filtFsSampleNames) #none are missing

# Sample names in order
phylloITS_sample.names <- substring(basename(phylloITS_filtFs), regexpr("_", basename(phylloITS_filtFs)) + 1) # doesn't drop fastq.gz
phylloITS_sample.names <- gsub(".fastq.gz", "", phylloITS_sample.names)
phylloITS_sample.namesR <- substring(basename(phylloITS_filtRs), regexpr("_", basename(phylloITS_filtRs)) + 1) # doesn't drop fastq.gz
phylloITS_sample.namesR <- gsub(".fastq.gz", "", phylloITS_sample.namesR)

# Double check
if(!identical(phylloITS_sample.names, phylloITS_sample.namesR)) stop("Forward and reverse files do not match.")
names(phylloITS_filtFs) <- phylloITS_sample.names
names(phylloITS_filtRs) <- phylloITS_sample.names

# #### Learn the error rates
set.seed(100) # set seed to ensure that randomized steps are replicatable 

# Learn forward error rates (Notes: randomize default is FALSE)
phylloITS_errF <- learnErrors(phylloITS_filtFs, nbases = 1e8, multithread = TRUE, randomize = TRUE) 

# Learn reverse error rates
phylloITS_errR <- learnErrors(phylloITS_filtRs, nbases = 1e8, multithread = TRUE, randomize = TRUE)

# #### Plot Error Rates
# We want to make sure that the machine learning algorithm is learning the error rates properly. In the plots below, the red line represents what we should expect the learned error rates to look like for each of the 16 possible base transitions (A->A, A->C, A->G, etc.) and the black line and grey dots represent what the observed error rates are. If the black line and the red lines are very far off from each other, it may be a good idea to increase the ```nbases``` parameter. This alows the machine learning algorthim to train on a larger portion of your data and may help imporve the fit.

phylloITS_errF_plot <- plotErrors(phylloITS_errF, nominalQ = TRUE) 
phylloITS_errF_plot
# Error: Warning messages:
# 1: Transformation introduced infinite values in continuous y-axis 
# 2: Transformation introduced infinite values in continuous y-axis 

# Message from Ben Callahan (dada2 creator) on an online forum
# This isn't an error, just a message from the plotting function to let you know that there were some zero values 
# in the data plotted (which turn into infinities on the log-scale). That is completely expected, it results 
# from the fact that not every combination of error type (e.g. A->C) and quality score (e.g. 33) is observed 
# in your data which is normal.

phylloITS_errR_plot <- plotErrors(phylloITS_errR, nominalQ = TRUE)
phylloITS_errR_plot
#
# write to disk
saveRDS(phylloITS_errF_plot, paste0(phylloITS_filtpathF, "/phylloITS_errF_plot.rds"))
saveRDS(phylloITS_errR_plot, paste0(phylloITS_filtpathR, "/phylloITS_errR_plot.rds"))

# #### Dereplication, sequence inference, and merging of paired-end reads
# In this part of the pipeline, dada2 will make decisions about assigning sequences to ASVs (called "sequence inference"). There is a major parameter option in the core function dada() that changes how samples are handled during sequence inference. The parameter ```pool = ``` can be set to: ```pool = FALSE``` (default), ```pool = TRUE```, or ```pool = psuedo```. For details on parameter choice, please see below, and further information on this blogpost [http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/](http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/), and explanation on the dada2 tutorial [https://benjjneb.github.io/dada2/pool.html](https://benjjneb.github.io/dada2/pool.html).

# # SAMPLES POOLED 
# For complex communities when you want to preserve rare taxa

# Dereplicate forward reads
phylloITS_derepF.p <- derepFastq(phylloITS_filtFs)
names(phylloITS_derepF.p) <- phylloITS_sample.names
# Infer sequences for forward reads
phylloITS_dadaF.p <- dada(phylloITS_derepF.p, err = phylloITS_errF, multithread = TRUE, pool = TRUE)
names(phylloITS_dadaF.p) <- phylloITS_sample.names

# Dereplicate reverse reads
phylloITS_derepR.p <- derepFastq(phylloITS_filtRs)
names(phylloITS_derepR.p) <- phylloITS_sample.names
# Infer sequences for reverse reads
phylloITS_dadaR.p <- dada(phylloITS_derepR.p, err = phylloITS_errR, multithread = TRUE, pool = TRUE)
names(phylloITS_dadaR.p) <- phylloITS_sample.names

# Merge reads together
phylloITS_mergers <- mergePairs(phylloITS_dadaF.p, phylloITS_derepF.p, phylloITS_dadaR.p, phylloITS_derepR.p)
head(phylloITS_mergers)

# #### Construct sequence table
phylloITS_seqtab <- makeSequenceTable(phylloITS_mergers)

# Save table as an r data object file
dir.create(phylloITS_table.fp)
saveRDS(phylloITS_seqtab, paste0(phylloITS_table.fp, "/phylloITS_seqtab.rds"))

###############################################################################
# ### 3. REMOVE Chimeras and ASSIGN Taxonomy
# Although dada2 has searched for indel errors and subsitutions, there may still be chimeric
# sequences in our dataset (sequences that are derived from forward and reverse sequences from 
# two different organisms becoming fused together during PCR and/or sequencing). To identify 
# chimeras, we will search for rare sequence variants that can be reconstructed by combining
# left-hand and right-hand segments from two more abundant "parent" sequences. After removing chimeras, we will use a taxonomy database to train a classifier-algorithm
# to assign names to our sequence variants.
# 

# Read in all of the RDS for the different sample types made earlier:
airSeqTab <- readRDS(paste0(airITS_table.fp, "/airITS_seqtab.rds"))
phylloSeqTab <- readRDS(paste0(phylloITS_table.fp, "/phylloITS_seqtab.rds"))
soilSeqtab <- readRDS("~/SRS_aeromicrobiome_2022/soilsfilteredPathF/soilITS_table/soilITS_seqtab.rds")

# Combine all these seqtabs
ITS_st.all <- mergeSequenceTables(airSeqTab, phylloSeqTab, soilSeqtab)

# Remove chimeras
ITSall_seqtab.nochim <- removeBimeraDenovo(ITS_st.all, method="consensus", multithread=TRUE)

# Print percentage of our sequences that were not chimeric.
100*sum(ITSall_seqtab.nochim)/sum(ITS_st.all) #99.17567

# Assign taxonomy
ITSall_tax <- assignTaxonomy(ITSall_seqtab.nochim, "/db_files/dada2/sh_general_release_dynamic_10.05.2021.fasta", tryRC = TRUE,
                             multithread=TRUE)
summary(ITSall_tax)  # in the reference database it looks like some of the fungal sequences have all the way to species information

# Write results to disk
#saveRDS(ITSall_seqtab.nochim, paste0(project.fp, "/ITSall_seqtab_final.rds")) #saved April 11, 2023
#saveRDS(ITSall_tax, paste0(project.fp, "/ITSall_tax_final.rds")) #saved April 11, 2023

# ### 4. Optional - FORMAT OUTPUT to obtain ASV IDs and repset, and input for mctoolsr
# For convenience sake, we will now rename our ASVs with numbers, output our 
# results as a traditional taxa table, and create a matrix with the representative
# sequences for each ASV. 

ITSall_tax <- readRDS("/data/winfreyc/SRS_aeromicrobiome_2022/ITSall_tax_final.rds")
ITSall_seqtab.nochim <- readRDS("/data/winfreyc/SRS_aeromicrobiome_2022/ITSall_seqtab_final.rds")
# Flip table
allITS_seqtab.t <- as.data.frame(t(ITSall_seqtab.nochim))

# Pull out ASV repset
allITS_rep_set_ASVs <- as.data.frame(rownames(allITS_seqtab.t))
allITS_rep_set_ASVs <- mutate(allITS_rep_set_ASVs, ASV_ID = 1:n())
allITS_rep_set_ASVs$ASV_ID <- sub("^", "ASV_", allITS_rep_set_ASVs$ASV_ID)
allITS_rep_set_ASVs$ASV <- allITS_rep_set_ASVs$`rownames(allITS_seqtab.t)` 
allITS_rep_set_ASVs$`rownames(allITS_seqtab.t)` <- NULL

# Add ASV numbers to table
rownames(allITS_seqtab.t) <- allITS_rep_set_ASVs$ASV_ID

# Add ASV numbers to taxonomy
allITS_taxonomy <- as.data.frame(ITSall_tax) 
allITS_taxonomy$ASV <- as.factor(rownames(allITS_taxonomy))
allITS_taxonomy <- merge(allITS_rep_set_ASVs, allITS_taxonomy, by = "ASV")
rownames(allITS_taxonomy) <- allITS_taxonomy$ASV_ID
allITS_taxonomy_for_mctoolsr <- unite(allITS_taxonomy, "allITS_taxonomy", 
                                       c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ASV_ID"),
                                       sep = ";")

# Write repset to fasta file
# create a function that writes fasta sequences
writeRepSetFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# Arrange the allITS_taxonomy dataframe for the writeRepSetFasta function
allITS_taxonomy_for_fasta <- allITS_taxonomy %>%
  unite("TaxString", c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ASV_ID"), 
        sep = ";", remove = FALSE) %>%
  unite("name", c("ASV_ID", "TaxString"), 
        sep = " ", remove = TRUE) %>%
  select(ASV, name) %>%
  rename(seq = ASV)

# write fasta file- written April 11, 2023
writeRepSetFasta(allITS_taxonomy_for_fasta, "/data/winfreyc/SRS_aeromicrobiome_2022/allITS_repset.fasta")

# Merge allITS_taxonomy and table
allITS_seqtab_wTax <- merge(allITS_seqtab.t, allITS_taxonomy_for_mctoolsr, by = 0)
allITS_seqtab_wTax$ASV <- NULL 

# Set name of table in mctoolsr format and save
names(allITS_seqtab_wTax)[1] = "#ASV_ID"
write("#Exported for mctoolsr", "~/SRS_aeromicrobiome_2022/allITS_seqtab_wTax_mctoolsr.txt")
suppressWarnings(write.table(allITS_seqtab_wTax, "~/SRS_aeromicrobiome_2022/allITS_seqtab_wTax_mctoolsr.txt", sep = "\t", row.names = FALSE, append = TRUE))

# Also export files as .txt
write.table(allITS_seqtab.t, file = "~/SRS_aeromicrobiome_2022/allITS_seqtab_final.txt",
            sep = "\t", row.names = TRUE, col.names = NA)
write.table(ITSall_tax, file = "~/SRS_aeromicrobiome_2022/allITS_tax_final.txt", 
            sep = "\t", row.names = TRUE, col.names = NA)

# ### Summary of output files:
# 1. allITS_seqtab_final.txt - A tab-delimited sequence-by-sample (i.e. OTU) table 
# 2. allITS_tax_final.txt - a tab-demilimited file showing the relationship between ASVs, ASV IDs, and their taxonomy 
# 3. allITS_seqtab_wTax_mctoolsr.txt - a tab-delimited file with ASVs as rows, samples as columns and the final column showing the taxonomy of the ASV ID 
# 4. allITS_repset.fasta - a fasta file with the representative sequence of each ASV. Fasta headers are the ASV ID and taxonomy string.  
#

# ### 5. Summary of reads throughout pipeline
# Here we track the reads throughout the pipeline to see if any step is resulting in a greater-than-expected loss of reads. If a step is showing a greater than expected loss of reads, it is a good idea to go back to that step and troubleshoot why reads are dropping out. The dada2 tutorial has more details about what can be changed at each step. 
# 
getN <- function(x) sum(getUniques(x)) # function to grab sequence counts from output objects
# tracking reads by counts
allITS_filt_out_track <- allITS_filt_out %>%
  data.frame() %>%
  mutate(Sample = gsub("(R1\\_)(.{1,})(\\.fastq\\.gz)","\\2",rownames(.))) %>%
  rename(input = reads.in, filtered = reads.out)
rownames(allITS_filt_out_track) <- allITS_filt_out_track$Sample

# ###### APRIL 11, 2023 Here on does not work because ddf is not a defined object (and is needed for next line of code)
# #### Will go back to DADA2 pipeline and try to figure out the issue
# ddF_track <- data.frame(denoisedF = sapply(ddF[allITS_sample.names], getN)) %>%
#   mutate(Sample = row.names(.))
# ddR_track <- data.frame(denoisedR = sapply(ddR[allITS_sample.names], getN)) %>%
#   mutate(Sample = row.names(.))
# merge_track <- data.frame(merged = sapply(allITS_mergers, getN)) %>%
#   mutate(Sample = row.names(.))
# chim_track <- data.frame(nonchim = rowSums(allITS_seqtab.nochim)) %>%
#   mutate(Sample = row.names(.))
# 
# 
# track <- left_join(allITS_filt_out_track, ddF_track, by = "Sample") %>%
#   left_join(ddR_track, by = "Sample") %>%
#   left_join(merge_track, by = "Sample") %>%
#   left_join(chim_track, by = "Sample") %>%
#   replace(., is.na(.), 0) %>%
#   select(Sample, everything())
# row.names(track) <- track$Sample
# head(track)
# 
# # tracking reads by percentage
# track_pct <- track %>% 
#   data.frame() %>%
#   mutate(Sample = rownames(.),
#          filtered_pct = ifelse(filtered == 0, 0, 100 * (filtered/input)),
#          denoisedF_pct = ifelse(denoisedF == 0, 0, 100 * (denoisedF/filtered)),
#          denoisedR_pct = ifelse(denoisedR == 0, 0, 100 * (denoisedR/filtered)),
#          merged_pct = ifelse(merged == 0, 0, 100 * merged/((denoisedF + denoisedR)/2)),
#          nonchim_pct = ifelse(nonchim == 0, 0, 100 * (nonchim/merged)),
#          total_pct = ifelse(nonchim == 0, 0, 100 * nonchim/input)) %>%
#   select(Sample, ends_with("_pct"))
# 
# # summary stats of tracked reads averaged across samples
# track_pct_avg <- track_pct %>% summarize_at(vars(ends_with("_pct")), 
#                                             list(avg = mean))
# head(track_pct_avg)
# 
# track_pct_med <- track_pct %>% summarize_at(vars(ends_with("_pct")), 
#                                             list(avg = stats::median))
# head(track_pct_avg)
# head(track_pct_med)
# 
# # Plotting each sample's reads through the pipeline
# track_plot <- track %>% 
#   data.frame() %>%
#   mutate(Sample = rownames(.)) %>%
#   gather(key = "Step", value = "Reads", -Sample) %>%
#   mutate(Step = factor(Step, 
#                        levels = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"))) %>%
#   ggplot(aes(x = Step, y = Reads)) +
#   geom_line(aes(group = Sample), alpha = 0.2) +
#   geom_point(alpha = 0.5, position = position_jitter(width = 0)) + 
#   stat_summary(fun.y = median, geom = "line", group = 1, color = "steelblue", size = 1, alpha = 0.5) +
#   stat_summary(fun.y = median, geom = "point", group = 1, color = "steelblue", size = 2, alpha = 0.5) +
#   stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), 
#                geom = "ribbon", group = 1, fill = "steelblue", alpha = 0.2) +
#   geom_label(data = t(track_pct_avg[1:5]) %>% data.frame() %>% 
#                rename(Percent = 1) %>%
#                mutate(Step = c("filtered", "denoisedF", "denoisedR", "merged", "nonchim"),
#                       Percent = paste(round(Percent, 2), "%")),
#              aes(label = Percent), y = 1.1 * max(track[,2])) +
#   geom_label(data = track_pct_avg[6] %>% data.frame() %>%
#                rename(total = 1),
#              aes(label = paste("Total\nRemaining:\n", round(track_pct_avg[1,6], 2), "%")), 
#              y = mean(track[,6]), x = 6.5) +
#   expand_limits(y = 1.1 * max(track[,2]), x = 7) +
#   theme_classic()
# 
# track_plot
# #
# # Write results to disk
# saveRDS(track, paste0(project.fp, "/tracking_reads.rds"))
# saveRDS(track_pct, paste0(project.fp, "/tracking_reads_percentage.rds"))
# saveRDS(track_plot, paste0(project.fp, "/tracking_reads_summary_plot.rds"))
# 
# # NEXT STEPS:
# # This output of this script is processed in ITSExploratoryDataAnalysis.R before further downstream analysis
