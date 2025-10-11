# 16S BIOINFORMATICS PIPELINE FOR SRS AEROMICROBIOME PROJECT (air, phyllosphere, and soil samples run together)
# started May 1, 2023
# Claire Winfrey

# First, this script demultiplexes raw reads from the MiSeq (separately for air and phyllosphere
# since they were on different MiSeq runs). Next, working with the air, phyllosphere, and soil (only 
# matched samples to air collection) data TOGETHER, I use dada2 to trim, denoise, dereplicate,
# and merge reads, and the assign taxonomy. 
# Finally, this script formats the data for use in R.

# *NOTE* Some of the important diagnostic plots (learning error rates) were not able to be plotted 
# using our server, so I scp'd them to my personal computer and did all the checks before continuing.

# This script is based on a combination of 1) Fierer Lab DADA2 tutorial created by Angela 
# Oliverio and Hannah Holland-Moritz (Mar 2, 2020 version, accessed from: 
# https://github.com/fiererlab/dada2_fiererlab) AND 
# 2) Benjamin Callahan's DADA2 16S pipeline tutorial (last accessed Jan. 19, 2022 at 
# https://benjjneb.github.io/dada2/16S_workflow.html)

####################
#    SET UP 
####################

library(dada2); packageVersion("dada2") # the dada2 pipeline, version ‘1.22.0’, May 19- 1.26.0
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this ‘1.52.0’-- May 19- 1.56.1
library(dplyr); packageVersion("dplyr") # for manipulating data #dplyr is working ‘1.0.10’, May 19- 1.1.1
library(tidyr); packageVersion("tidyr") # for creating the final graph at the end of the pipeline ‘1.2.0’, May 19- 1.3.0
library(Hmisc); packageVersion("Hmisc") # for creating the final graph at the end of the pipeline ‘4.6.0’, May 19- 5.0.1
library(ggplot2); packageVersion("ggplot2") # for creating the final graph at the end of the pipeline ‘3.3.5’, May 19 -3.4.1
library(plotly); packageVersion("plotly") # enables creation of interactive graphs, especially helpful for quality plots ‘4.10.0’, May 19 4.10.1
# Once the packages are installed, you can check to make sure the auxillary
# software is working and set up some of the variables that you will need 
# along the way.

# Set up pathway to idemp (demultiplexing tool) and test
idemp <- "/usr/bin/idemp"
system2(idemp) # Check that idemp is in your path and you can run shell commands from R

# Set up pathway to cutadapt (primer trimming tool) and test
cutadapt <- "/usr/local/Python27/cutadapt" #for microbe server
system2(cutadapt, args = "--version") # Check by running shell command from R
# 05/01/23: version is 1.8.1

### SET UP FILE PATHS FOR DATA AND BARCODES ###
# 1. AIR 16S SAMPLES
# Set path to shared data folder and contents
air16Sdata.fp <- "/data/shared/MiSeq/03.22.2023_Mixed16S_Claire_Nick"
# List all files in shared folder to check path
list.files(air16Sdata.fp) #Data looks correct
# Set file paths for barcodes file, map file, and fastqs
air16Sbarcode.fp <- file.path(air16Sdata.fp, "air16S_MappingFile")
air16S_I1.fp <- file.path(air16Sdata.fp, "Undetermined_S0_L001_I1_001.fastq.gz") 
air16S_R1.fp <- file.path(air16Sdata.fp, "Undetermined_S0_L001_R1_001.fastq.gz") 
air16S_R2.fp <- file.path(air16Sdata.fp, "Undetermined_S0_L001_R2_001.fastq.gz") 

# 2. PHYLLOSPHERE 16S SAMPLES
# Set path to shared data folder and contents
phyllo16Sdata.fp <- "/data/shared/MiSeq/11.11.2022_MixRun16S_try2"
# List all files in shared folder to check path
list.files(phyllo16Sdata.fp) 
# Set file paths for barcodes file, map file, and fastqs
phyllo16Sbarcode.fp <- file.path(phyllo16Sdata.fp, "phyllo16S_mappingFile") #used sequences from #Barcodes_demultiplexing column in 

phyllo16S_I1.fp <- file.path(phyllo16Sdata.fp, "Undetermined_S0_L001_I1_001.fastq.gz") 
phyllo16S_R1.fp <- file.path(phyllo16Sdata.fp, "Undetermined_S0_L001_R1_001.fastq.gz") 
phyllo16S_R2.fp <- file.path(phyllo16Sdata.fp, "Undetermined_S0_L001_R2_001.fastq.gz") 

# SET UP FILE PATHS FOR WHERE RESULTS WILL GO
project.fp <- "/data/winfreyc/SRS_aeromicrobiome_2022" #this is where all of the results will go! 
list.files(project.fp) #when first doing this script, this only had copy of two mapping files

# Set up file paths/names of sub directories to stay organized
#all of these are now in SRS_aeromicrobiome_2022 folder
## For these, the processing has to be done separately for air and phyllo samples, since they
## were on separate runs.
air16S_preprocess.fp <- file.path(project.fp, "air16S_01_preprocess") 
air16S_demultiplex.fp <- file.path(air16S_preprocess.fp, "air16S_demultiplexed")
air16S_filtN.fp <- file.path(air16S_preprocess.fp, "air16S_filtN") 
air16S_trimmed.fp <- file.path(air16S_preprocess.fp, "air16S_trimmed")
phyllo16S_preprocess.fp <- file.path(project.fp, "phyllo16S_01_preprocess") 
phyllo16S_demultiplex.fp <- file.path(phyllo16S_preprocess.fp, "phyllo16S_demultiplexed")
phyllo16S_filtN.fp <- file.path(phyllo16S_preprocess.fp, "phyllo16S_filtN") 
phyllo16S_trimmed.fp <- file.path(phyllo16S_preprocess.fp, "phyllo16S_trimmed")

########################################################
#    # PRE-PROCESSING DATA FOR DADA2 - 
# DEMULTIPLEX, REMOVE SEQUENCES WITH Ns, CUTADAPT
########################################################
# This is done separately for air and phyllosphere samples, since they were run on different 16S runs.
###### AIR 16S SAMPLES ######
### 1. DEMULTIPLEX ###
# #### Call the demultiplexing script
# Demultiplexing spl16S your reads out into separate files based on the barcodes associated with each sample. 
air16S_Demultiplex <- paste("-b", air16Sbarcode.fp, "-I1", air16S_I1.fp, "-R1", air16S_R1.fp, "-R2", air16S_R2.fp, "-o", air16S_demultiplex.fp) 
system2(idemp, args = air16S_Demultiplex)

# Look at output of demultiplexing
list.files(air16S_demultiplex.fp) 
sort(list.files(air16S_demultiplex.fp))
list.files(air16S_demultiplex.fp, pattern = "R1") 
list.files(air16S_demultiplex.fp, pattern = "R2") 
list.files(air16S_demultiplex.fp, pattern = "unsigned") #2 "unsigned file

### 2. CLEAN UP THE OUTPUT FROM IDEMP ###
# Change names of unassignable reads so they are not included in downstream processing
air16S_unassigned_1 <- paste0("mv", " ", air16S_demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_unsigned.fastq.gz",
                              " ", air16S_demultiplex.fp, "/Unassigned_reads1.fastq.gz") #this moves undetermined read 1 sequences into a new .gz file 
air16S_unassigned_2 <- paste0("mv", " ", air16S_demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_unsigned.fastq.gz", 
                              " ", air16S_demultiplex.fp, "/Unassigned_reads2.fastq.gz") #likewise, this moves undetermined R2 reads to their own folder
system(air16S_unassigned_1)
system(air16S_unassigned_2)

# Rename files - use gsub to get names in order!
air16S_R1_names <- gsub(paste0(air16S_demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_"), "", 
                        list.files(air16S_demultiplex.fp, pattern="R1", full.names = TRUE))
file.rename(list.files(air16S_demultiplex.fp, pattern="R1", full.names = TRUE), 
            paste0(air16S_demultiplex.fp, "/R1_", air16S_R1_names))

air16S_R2_names <- gsub(paste0(air16S_demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_"), "", 
                        list.files(air16S_demultiplex.fp, pattern="R2", full.names = TRUE))
file.rename(list.files(air16S_demultiplex.fp, pattern="R2", full.names = TRUE),
            paste0(air16S_demultiplex.fp, "/R2_", air16S_R2_names))

# Get full paths for all files and save them for downstream analyses
# Forward and reverse fastq filenames have format: (fnFs and frRs stand for "full name forward" and
# "full name reverse" )
air16S_fnFs <- sort(list.files(air16S_demultiplex.fp, pattern="R1_", full.names = TRUE))
air16S_fnRs <- sort(list.files(air16S_demultiplex.fp, pattern="R2_", full.names = TRUE))

### 3. PRE-FILTER TO REMOVE SEQUENCE READS WITH NS ###
# Ambiguous bases will make it hard for cutadapt to find short primer sequences in the reads.
# To solve this problem, we will remove sequences with ambiguous bases (Ns)

# Name the N-filtered files to put them in filtN/ subdirectory
air16S_fnFs.filtN <- file.path(air16S_preprocess.fp, "filtN", basename(air16S_fnFs))
air16S_fnRs.filtN <- file.path(air16S_preprocess.fp, "filtN", basename(air16S_fnRs))

# Filter Ns from reads and put them into the filtN directory (file path is "SRS_aeromicrobiome_2022/air16S_01_preprocess/filtN")
filterAndTrim(air16S_fnFs, air16S_fnFs.filtN, air16S_fnRs, air16S_fnRs.filtN, maxN = 0, multithread = TRUE) 

### 4. CUTADAPT ###
# The setting up of the primer sequences and the custom functions are the same for 
# the air and phyllosphere processing!
# #### Prepare the primers sequences and custom functions for analyzing the results from cutadapt
# Assign the primers you used to "FWD" and "REV" below. Note primers should be not be reverse complemented ahead of time. Our tutorial data uses 515f and 806br those are the primers below. Change if you sequenced with other primers

# Set up the primer sequences to pass along to cutadapt
FWD <- "GTGYCAGCMGCCGCGGTAA"  # this is 515f
REV <- "GGACTACNVGGGTWTCTAAT"  # this is 806Br

# Save the reverse complements of the primers to variables
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

##  Create the cutadapt flags ##
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 50") 

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 50") 

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
FWD.orients

# Write a function that counts how many time primers appear in a sequence
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Before running cutadapt, we will look at primer detection for the first sample, as a check. There may be some primers here, we will remove them below using cutadapt.
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = air16S_fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = air16S_fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = air16S_fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = air16S_fnRs.filtN[[1]]))
# Forward Complement Reverse RevComp
# FWD.ForwardReads       0          0       0       0
# FWD.ReverseReads       0          0       0       2
# REV.ForwardReads       0          0       0       8
# REV.ReverseReads       2          0       0       0

## iii. Remove primers with cutadapt and assess the output ##

# Create directory to hold the output from cutadapt
if (!dir.exists(air16S_trimmed.fp)) dir.create(air16S_trimmed.fp)
air16S_fnFs.cut <- file.path(air16S_trimmed.fp, basename(air16S_fnFs))
air16S_fnRs.cut <- file.path(air16S_trimmed.fp, basename(air16S_fnRs))
# An example file path from the reverses here is:
# "/data/winfreyc/SRS_aeromicrobiome_2022/air16S_01_preprocess/air16S_trimmed/R2_air_16S_PCR_NTC_2.fastq.gz

# Run Cutadapt
for (i in seq_along(air16S_fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", air16S_fnFs.cut[i], "-p", air16S_fnRs.cut[i], # output files
                             air16S_fnFs.filtN[i], air16S_fnRs.filtN[i])) # input files
}

# As a sanity check, we will check for primers in the first cutadapt-ed sample:
air16S_CutAdaptCheck <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = air16S_fnFs.cut[[1]]), 
                              FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = air16S_fnRs.cut[[1]]), 
                              REV.ForwardReads = sapply(REV.orients, primerHits, fn = air16S_fnFs.cut[[1]]), 
                              REV.ReverseReads = sapply(REV.orients, primerHits, fn = air16S_fnRs.cut[[1]]))

air16S_CutAdaptCheck

####################################################################################

###### PHYLLOSPHERE 16S SAMPLES ######
### 1. DEMULTIPLEX PHYLLOSPHERE 16S samples ###
phyllo16S_Demultiplex <- paste("-b", phyllo16Sbarcode.fp, "-I1", phyllo16S_I1.fp, "-R1", phyllo16S_R1.fp, "-R2", phyllo16S_R2.fp, "-o", phyllo16S_demultiplex.fp) 
system2(idemp, args = phyllo16S_Demultiplex)

# Look at output of demultiplexing
list.files(phyllo16S_demultiplex.fp) 
sort(list.files(phyllo16S_demultiplex.fp))
list.files(phyllo16S_demultiplex.fp) #136 files
list.files(phyllo16S_demultiplex.fp, pattern = "unsigned") #2 "unsigned" files

### 2. CLEAN UP THE OUTPUT FROM IDEMP ###
# Change names of unassignable reads so they are not included in downstream processing
phyllo16S_unassigned_1 <- paste0("mv", " ", phyllo16S_demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_unsigned.fastq.gz",
                                 " ", phyllo16S_demultiplex.fp, "/Unassigned_reads1.fastq.gz") #this moves undetermined read 1 sequences into a new .gz file 
phyllo16S_unassigned_2 <- paste0("mv", " ", phyllo16S_demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_unsigned.fastq.gz", 
                                 " ", phyllo16S_demultiplex.fp, "/Unassigned_reads2.fastq.gz") #likewise, this moves undetermined R2 reads to their own folder
system(phyllo16S_unassigned_1)
system(phyllo16S_unassigned_2)

# Rename files - use gsub to get names in order!
phyllo16S_R1_names <- gsub(paste0(phyllo16S_demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_"), "", 
                           list.files(phyllo16S_demultiplex.fp, pattern="R1", full.names = TRUE))
file.rename(list.files(phyllo16S_demultiplex.fp, pattern="R1", full.names = TRUE), 
            paste0(phyllo16S_demultiplex.fp, "/R1_", phyllo16S_R1_names))

phyllo16S_R2_names <- gsub(paste0(phyllo16S_demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_"), "", 
                           list.files(phyllo16S_demultiplex.fp, pattern="R2", full.names = TRUE))
file.rename(list.files(phyllo16S_demultiplex.fp, pattern="R2", full.names = TRUE),
            paste0(phyllo16S_demultiplex.fp, "/R2_", phyllo16S_R2_names))

# Get full paths for all files and save them for downstream analyses
# Forward and reverse fastq filenames have format: (fnFs and frRs stand for "full name forward" and
# "full name reverse" )
phyllo16S_fnFs <- sort(list.files(phyllo16S_demultiplex.fp, pattern="R1_", full.names = TRUE))
phyllo16S_fnRs <- sort(list.files(phyllo16S_demultiplex.fp, pattern="R2_", full.names = TRUE))

#### 3. PRE-FILTER TO REMOVE SEQUENCE READS WITH NS ####
# Name the N-filtered files to put them in filtN/ subdirectory (so file path is "SRS_aeromicrobiome_2022/phyllo16S_01_preprocess/filtN")
phyllo16S_fnFs.filtN <- file.path(phyllo16S_preprocess.fp, "filtN", basename(phyllo16S_fnFs))
phyllo16S_fnRs.filtN <- file.path(phyllo16S_preprocess.fp, "filtN", basename(phyllo16S_fnRs))

# Filter Ns from reads and put them into the filtN directory (filtN directory is also made during this step)
filterAndTrim(phyllo16S_fnFs, phyllo16S_fnFs.filtN, phyllo16S_fnRs, phyllo16S_fnRs.filtN, maxN = 0, multithread = TRUE) 

#### 4. CUTADAPT ####
# (the set up is the same as in air samples, see above)
## ii. Prepare the primers sequences and custom functions for analyzing the results from cutadapt ####
# Before running cutadapt, we will look at primer detection for the first sample, as a check. There may be some primers here, we will remove them below using cutadapt.
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = phyllo16S_fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = phyllo16S_fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = phyllo16S_fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = phyllo16S_fnRs.filtN[[1]]))

# Forward Complement Reverse RevComp
# FWD.ForwardReads       0          0       0       0
# FWD.ReverseReads       0          0       0      20
# REV.ForwardReads       0          0       0      21
# REV.ReverseReads       0          0       0       0

## iii. Remove primers with cutadapt and assess the output ##

# Create directory to hold the output from cutadapt
if (!dir.exists(phyllo16S_trimmed.fp)) dir.create(phyllo16S_trimmed.fp)
phyllo16S_fnFs.cut <- file.path(phyllo16S_trimmed.fp, basename(phyllo16S_fnFs))
phyllo16S_fnRs.cut <- file.path(phyllo16S_trimmed.fp, basename(phyllo16S_fnRs))
# An example file path from the reverses here is:
# "/data/winfreyc/SRS_aeromicrobiome_2022/phyllo16S_01_preprocess/phyllo16S_trimmed/R2_phyllo_16S_PCR_NTC_2.fastq.gz

# Run Cutadapt
for (i in seq_along(phyllo16S_fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", phyllo16S_fnFs.cut[i], "-p", phyllo16S_fnRs.cut[i], # output files
                             phyllo16S_fnFs.filtN[i], phyllo16S_fnRs.filtN[i])) # input files
}

# As a sanity check, we will check for primers in the first cutadapt-ed sample:
# got all zeros!
phyllo16S_CutAdaptCheck <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = phyllo16S_fnFs.cut[[1]]), 
                                 FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = phyllo16S_fnRs.cut[[1]]), 
                                 REV.ForwardReads = sapply(REV.orients, primerHits, fn = phyllo16S_fnFs.cut[[1]]), 
                                 REV.ReverseReads = sapply(REV.orients, primerHits, fn = phyllo16S_fnRs.cut[[1]]))

phyllo16S_CutAdaptCheck #all zeros, showing this worked

######################################################
#### DADA2 PIPELINE ####
######################################################
#### PRE-PIPELINE. GET ORGANIZED ####
## Create different sub directories for each different kind of sample
# (already processed soil data will be brought in later, when using 
# subset of each sample type, separately, to infer error rates.)

air16S_filter.fp <- file.path(project.fp, "air16S_02_filter") 
air16S_table.fp <- file.path(project.fp, "air16S_03_tabletax") 
phyllo16S_filter.fp <- file.path(project.fp, "phyllo16S_02_filter") 
phyllo16S_table.fp <- file.path(project.fp, "phyllo16S_03_tabletax")

# Put filtered reads into separate sub-directories for big data workflow
###### AIR ######
dir.create(air16S_filter.fp)
air16S_subF.fp <- file.path(air16S_filter.fp, "preprocessed_F") 
air16S_subR.fp <- file.path(air16S_filter.fp, "preprocessed_R") 
dir.create(air16S_subF.fp)
dir.create(air16S_subR.fp)
### PHYLLO ###
dir.create(phyllo16S_filter.fp)
phyllo16S_subF.fp <- file.path(phyllo16S_filter.fp, "preprocessed_F") 
phyllo16S_subR.fp <- file.path(phyllo16S_filter.fp, "preprocessed_R") 
dir.create(phyllo16S_subF.fp)
dir.create(phyllo16S_subR.fp)

# Move R1 and R2 from trimmed to separate forward/reverse sub-directories
###### AIR ######
air16S_fnFs.Q <- file.path(air16S_subF.fp,  basename(air16S_fnFs)) 
air16S_fnRs.Q <- file.path(air16S_subR.fp,  basename(air16S_fnRs))
file.rename(from = air16S_fnFs.cut, to = air16S_fnFs.Q)
file.rename(from = air16S_fnRs.cut, to = air16S_fnRs.Q)
### PHYLLO ###
phyllo16S_fnFs.Q <- file.path(phyllo16S_subF.fp,  basename(phyllo16S_fnFs)) 
phyllo16S_fnRs.Q <- file.path(phyllo16S_subR.fp,  basename(phyllo16S_fnRs))
file.rename(from = phyllo16S_fnFs.cut, to = phyllo16S_fnFs.Q)
file.rename(from = phyllo16S_fnRs.cut, to = phyllo16S_fnRs.Q)

# File parsing; create file names and make sure that forward and reverse files match
###### AIR ######
air16S_filtpathF <- file.path(air16S_subF.fp, "filtered") # files go into preprocessed_F/filtered/
air16S_filtpathR <- file.path(air16S_subR.fp, "filtered") # ...
air16S_fastqFs <- sort(list.files(air16S_subF.fp, pattern="fastq.gz"))
air16S_fastqRs <- sort(list.files(air16S_subR.fp, pattern="fastq.gz"))
if(length(air16S_fastqFs) != length(air16S_fastqRs)) stop("Forward and reverse files do not match.")
### PHYLLO ###
phyllo16S_filtpathF <- file.path(phyllo16S_subF.fp, "filtered") # files go into preprocessed_F/filtered/
phyllo16S_filtpathR <- file.path(phyllo16S_subR.fp, "filtered") # ...
phyllo16S_fastqFs <- sort(list.files(phyllo16S_subF.fp, pattern="fastq.gz"))
phyllo16S_fastqRs <- sort(list.files(phyllo16S_subR.fp, pattern="fastq.gz"))
if(length(phyllo16S_fastqFs) != length(phyllo16S_fastqRs)) stop("Forward and reverse files do not match.")

# ### 1. FILTER AND TRIM FOR QUALITY
###### AIR ######
# #### Inspect read quality profiles
# It's important to get a feel for the quality of the data that we are using. To do this, we will plot the quality of some of the samples.
# If the number of samples is 20 or less, plot them all, otherwise, just plot 20 randomly selected samples
if( length(air16S_fastqFs) <= 20) {
  plotQualityProfile(paste0(air16S_subF.fp, "/", air16S_fastqFs))
  plotQualityProfile(paste0(air16S_subR.fp, "/", air16S_fastqRs))
} else {
  rand_samples <- sample(size = 20, 1:length(air16S_fastqFs)) # grab 20 random samples to plot
  air16S_fwd_qual_plots <- plotQualityProfile(paste0(air16S_subF.fp, "/", air16S_fastqFs[rand_samples]))
  air16S_rev_qual_plots <- plotQualityProfile(paste0(air16S_subR.fp, "/", air16S_fastqRs[rand_samples]))
}

air16S_fwd_qual_plots 
air16S_rev_qual_plots

# Or, to make these quality plots interactive, just call the plots through plotly
ggplotly(air16S_fwd_qual_plots)
ggplotly(air16S_rev_qual_plots)

# write plots to disk- saved May 24, 2023 and then scp'd to my own computer (~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/bioinformaticsOutput)
# saveRDS(air16S_fwd_qual_plots, paste0(air16S_filter.fp, "/air16S_fwd_qual_plots.rds"))
# saveRDS(air16S_rev_qual_plots, paste0(air16S_filter.fp, "/air16S_rev_qual_plots.rds"))
# 
# ggsave(plot = air16S_fwd_qual_plots, filename = paste0(air16S_filter.fp, "/air16S_fwd_qual_plots.png"), 
#        width = 10, height = 10, dpi = "retina")
# ggsave(plot = air16S_rev_qual_plots, filename = paste0(air16S_filter.fp, "/air16S_rev_qual_plots.png"), 
#        width = 10, height = 10, dpi = "retina")

# #### Filter the data
air16S_filt_out <- filterAndTrim(fwd=file.path(air16S_subF.fp, air16S_fastqFs), filt=file.path(air16S_filtpathF, air16S_fastqFs),
              rev=file.path(air16S_subR.fp, air16S_fastqRs), filt.rev=file.path(air16S_filtpathR, air16S_fastqRs),
              truncLen=c(150,140), maxEE=1, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

dim(air16S_filt_out) 

# summary of samples in air16S_filt_out by percentage
air16S_filt_out %>% 
  data.frame() %>% 
  mutate(Samples = rownames(.),
         percent_kept = 100*(reads.out/reads.in)) %>%
  select(Samples, everything())

# Plot the quality of the filtered fastq files.
# figure out which samples, if any, have been filtered out (out of the random ones)
air16S_remaining_samplesF <-  air16S_fastqFs[rand_samples][
  which(air16S_fastqFs[rand_samples] %in% list.files(air16S_filtpathF))] # keep only samples that haven't been filtered out
air16S_remaining_samplesR <-  air16S_fastqRs[rand_samples][
  which(air16S_fastqRs[rand_samples] %in% list.files(air16S_filtpathR))] # keep only samples that haven't been filtered out
air16S_fwd_qual_plots_filt <- plotQualityProfile(paste0(air16S_filtpathF, "/", air16S_remaining_samplesF))
air16S_rev_qual_plots_filt <- plotQualityProfile(paste0(air16S_filtpathR, "/", air16S_remaining_samplesR))

air16S_fwd_qual_plots_filt
air16S_rev_qual_plots_filt

# write plots to disk- 
# saved May 24, 2023 and then scp'd to my own computer (~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/bioinformaticsOutput)
# saveRDS(air16S_fwd_qual_plots_filt, paste0(air16S_filter.fp, "/air16S_fwd_qual_plots_filt.rds"))
# saveRDS(air16S_rev_qual_plots_filt, paste0(air16S_filter.fp, "/air16S_rev_qual_plots_filt.rds"))
# 
# ggsave(plot = air16S_fwd_qual_plots_filt, filename = paste0(air16S_filter.fp, "/air16S_fwd_qual_plots_filt.png"), 
#        width = 10, height = 10, dpi = "retina")
# ggsave(plot = air16S_rev_qual_plots_filt, filename = paste0(air16S_filter.fp, "/air16S_rev_qual_plots_filt.png"), 
#        width = 10, height = 10, dpi = "retina")

###### PHYLLO ######
# #### Inspect read quality profiles
# It's important to get a feel for the quality of the data that we are using. To do this, we will plot the quality of some of the samples.
# If the number of samples is 20 or less, plot them all, otherwise, just plot 20 randomly selected samples
if( length(phyllo16S_fastqFs) <= 20) {
  plotQualityProfile(paste0(phyllo16S_subF.fp, "/", phyllo16S_fastqFs))
  plotQualityProfile(paste0(phyllo16S_subR.fp, "/", phyllo16S_fastqRs))
} else {
  rand_samples <- sample(size = 20, 1:length(phyllo16S_fastqFs)) # grab 20 random samples to plot
  phyllo16S_fwd_qual_plots <- plotQualityProfile(paste0(phyllo16S_subF.fp, "/", phyllo16S_fastqFs[rand_samples]))
  phyllo16S_rev_qual_plots <- plotQualityProfile(paste0(phyllo16S_subR.fp, "/", phyllo16S_fastqRs[rand_samples]))
}

# Interactive plots (plotly library)
ggplotly(phyllo16S_fwd_qual_plots)
ggplotly(phyllo16S_rev_qual_plots)

# write plots to disk
saveRDS(phyllo16S_fwd_qual_plots, paste0(phyllo16S_filter.fp, "/phyllo16S_fwd_qual_plots.rds"))
saveRDS(phyllo16S_rev_qual_plots, paste0(phyllo16S_filter.fp, "/phyllo16S_rev_qual_plots.rds"))

ggsave(plot = phyllo16S_fwd_qual_plots, filename = paste0(phyllo16S_filter.fp, "/phyllo16S_fwd_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = phyllo16S_rev_qual_plots, filename = paste0(phyllo16S_filter.fp, "/phyllo16S_rev_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")

# #### Filter the data

phyllo16S_filt_out <- filterAndTrim(fwd=file.path(phyllo16S_subF.fp, phyllo16S_fastqFs), filt=file.path(phyllo16S_filtpathF, phyllo16S_fastqFs),
              rev=file.path(phyllo16S_subR.fp, phyllo16S_fastqRs), filt.rev=file.path(phyllo16S_filtpathR, phyllo16S_fastqRs),
              truncLen=c(150,140), maxEE=1, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

# summary of samples in phyllo16S_filt_out by percentage
phyllo16S_filt_out %>% 
  data.frame() %>% 
  mutate(Samples = rownames(.),
         percent_kept = 100*(reads.out/reads.in)) %>%
  select(Samples, everything())

# Plot the quality of the filtered fastq files.
# figure out which samples, if any, have been filtered out (out of the random ones)
phyllo16S_remaining_samplesF <-  phyllo16S_fastqFs[rand_samples][
  which(phyllo16S_fastqFs[rand_samples] %in% list.files(phyllo16S_filtpathF))] # keep only samples that haven't been filtered out
phyllo16S_remaining_samplesR <-  phyllo16S_fastqRs[rand_samples][
  which(phyllo16S_fastqRs[rand_samples] %in% list.files(phyllo16S_filtpathR))] # keep only samples that haven't been filtered out
phyllo16S_fwd_qual_plots_filt <- plotQualityProfile(paste0(phyllo16S_filtpathF, "/", phyllo16S_remaining_samplesF))
phyllo16S_rev_qual_plots_filt <- plotQualityProfile(paste0(phyllo16S_filtpathR, "/", phyllo16S_remaining_samplesR))

phyllo16S_fwd_qual_plots_filt
phyllo16S_rev_qual_plots_filt

# write plots to disk (then scp and view on own computer (SRS_Aeromicrobiome/BioinformaticsAndMetadata/bioinformaticsOutput folder) 
# since ggplot not working on server) -- saved May 24, 2023

#saveRDS(phyllo16S_fwd_qual_plots_filt, paste0(phyllo16S_filter.fp, "/phyllo16S_fwd_qual_plots_filt.rds"))
#saveRDS(phyllo16S_rev_qual_plots_filt, paste0(phyllo16S_filter.fp, "/phyllo16S_rev_qual_plots_filt.rds"))

# ggsave(plot = phyllo16S_fwd_qual_plots_filt, filename = paste0(phyllo16S_filter.fp, "/phyllo16S_fwd_qual_plots_filt.png"), 
#        width = 10, height = 10, dpi = "retina")
# ggsave(plot = phyllo16S_rev_qual_plots_filt, filename = paste0(phyllo16S_filter.fp, "/phyllo16S_rev_qual_plots_filt.png"), 
#        width = 10, height = 10, dpi = "retina")

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
soilProcessedFfiltered.fp <- "~/SRS_May_2021/02_filter/preprocessed_F/filtered"
soilProcessedRfiltered.fp <- "~/SRS_May_2021/02_filter/preprocessed_R/filtered"
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
PCRBlanks <- list.files(soilProcessedRfiltered.fp)[grep(list.files(soilProcessedFfiltered.fp), pattern="NTC")]
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
soils_16S_filteredPathF.fp <- file.path("~/SRS_aeromicrobiome_2022/soils_16S_filteredPathF")
if (!dir.exists(soils_16S_filteredPathF.fp)) dir.create(soils_16S_filteredPathF.fp)
soils_16S_filteredPathR.fp <- file.path("~/SRS_aeromicrobiome_2022/soils_16S_filteredPathR")
if (!dir.exists(soils_16S_filteredPathR.fp)) dir.create(soils_16S_filteredPathR.fp)

# Copy all of the R1 and R2 files for the relevant soil samples into the ~/SRS_May_2021/02_filter/preprocessed_F/filtered folder
for (i in 1:length(soilSamples)) {
  system(paste0('cp ', soilProcessedFfiltered.fp, '/R1_', soilSamples[i], ".fastq.gz ", soils_16S_filteredPathF.fp))
  system(paste0('cp ', soilProcessedRfiltered.fp, '/R2_', soilSamples[i], ".fastq.gz ", soils_16S_filteredPathR.fp))
}

# Check to make sure that all of the files are in the correct location:
list.files(soils_16S_filteredPathF.fp) #looks good!
list.files(soils_16S_filteredPathR.fp) 

# #### Housekeeping step - set up and verify the file names for the output:
soil16S_filtFs <- list.files(soils_16S_filteredPathF.fp, pattern="fastq.gz", full.names = TRUE)
length(soil16S_filtFs) #186
soil16S_filtRs <- list.files(soils_16S_filteredPathR.fp, pattern="fastq.gz", full.names = TRUE)
length(soil16S_filtRs) #186

# Sample names in order
soil16S_sample.names <- substring(basename(soil16S_filtFs), regexpr("_", basename(soil16S_filtFs)) + 1) # doesn't drop fastq.gz
soil16S_sample.names <- gsub(".fastq.gz", "", soil16S_sample.names)
soil16S_sample.namesR <- substring(basename(soil16S_filtRs), regexpr("_", basename(soil16S_filtRs)) + 1) # doesn't drop fastq.gz
soil16S_sample.namesR <- gsub(".fastq.gz", "", soil16S_sample.namesR)

# Double check
if(!identical(soil16S_sample.names, soil16S_sample.namesR)) stop("Forward and reverse files do not match.")
names(soil16S_filtFs) <- soil16S_sample.names
names(soil16S_filtRs) <- soil16S_sample.names

# #### Learn the error rates
set.seed(100) # set seed to ensure that randomized steps are replicatable 
# Learn forward error rates (Notes: randomize default is FALSE)
soil16S_errF <- learnErrors(soil16S_filtFs, nbases = 1e8, multithread = TRUE, randomize = TRUE) 
# Learn reverse error rates
soil16S_errR <- learnErrors(soil16S_filtRs, nbases = 1e8, multithread = TRUE, randomize = TRUE)

# #### Plot Error Rates
# We want to make sure that the machine learning algorithm is learning the error rates properly. In the plots below,
# the red line represents what we should expect the learned error rates to look like for each of the 16 possible base
# transitions (A->A, A->C, A->G, etc.) and the black line and grey dots represent what the observed error rates are. 
# If the black line and the red lines are very far off from each other, it may be a good idea to increase the 
# ```nbases``` parameter. This alows the machine learning algorthim to train on a larger portion of your data and may
# help improve the fit.
soil16S_errF_plot <- plotErrors(soil16S_errF, nominalQ = TRUE) 
soil16S_errR_plot <- plotErrors(soil16S_errR, nominalQ = TRUE)

soil16S_errF_plot
soil16S_errR_plot
#
# write to disk- saved May 24/25, 2023
saveRDS(soil16S_errF_plot, paste0(soils_16S_filteredPathF.fp, "/soil16S_errF_plot.rds"))
saveRDS(soil16S_errR_plot, paste0(soils_16S_filteredPathF.fp, "/soil16S_errR_plot.rds"))
# saved to microbe and evaluated on own computer (since microbe does not like ggplot :< )
# # How to copy these to my computer
# scp winfreyc@microbe.colorado.edu:/data/winfreyc/SRS_aeromicrobiome_2022/soils_16S_filteredPathF/soil16S_errF_plot.rds ~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/bioinformaticsOutput
# scp winfreyc@microbe.colorado.edu:/data/winfreyc/SRS_aeromicrobiome_2022/soils_16S_filteredPathR/soil16S_errR_plot.rds ~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/bioinformaticsOutput

# #### Dereplication, sequence inference, and merging of paired-end reads
# In this part of the pipeline, dada2 will make decisions about assigning sequences to ASVs (called "sequence inference"). There is a major parameter option in the core function dada() that changes how samples are handled during sequence inference. The parameter ```pool = ``` can be set to: ```pool = FALSE``` (default), ```pool = TRUE```, or ```pool = psuedo```. For details on parameter choice, please see below, and further information on this blogpost [http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/](http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/), and explanation on the dada2 tutorial [https://benjjneb.github.io/dada2/pool.html](https://benjjneb.github.io/dada2/pool.html).

# # SAMPLES POOLED 
# For complex communities when you want to preserve rare taxa
# same steps, not in loop

# Dereplicate forward reads
soil16S_derepF.p <- derepFastq(soil16S_filtFs)
names(soil16S_derepF.p) <- soil16S_sample.names
# Infer sequences for forward reads
soil16S_dadaF.p <- dada(soil16S_derepF.p, err = soil16S_errF, multithread = TRUE, pool = TRUE)
names(soil16S_dadaF.p) <- soil16S_sample.names

# Dereplicate reverse reads
soil16S_derepR.p <- derepFastq(soil16S_filtRs)
names(soil16S_derepR.p) <- soil16S_sample.names
# Infer sequences for reverse reads
soil16S_dadaR.p <- dada(soil16S_derepR.p, err = soil16S_errR, multithread = TRUE, pool = TRUE)
names(soil16S_dadaR.p) <- soil16S_sample.names

# Merge reads together
soil16S_mergers <- mergePairs(soil16S_dadaF.p, soil16S_derepF.p, soil16S_dadaR.p, soil16S_derepR.p)
head(soil16S_mergers)

# #### Construct sequence table
soil16S_seqtab <- makeSequenceTable(soil16S_mergers)

# Save table as an r data object file
dir.create("~/SRS_aeromicrobiome_2022/soils_16S_filteredPathF/soil16S_table")
saveRDS(soil16S_seqtab, paste0("~/SRS_aeromicrobiome_2022/soils_16S_filteredPathF/soil16S_table", "/soil16S_seqtab.rds"))

###################################################################################
##### AIR ######

# #### Housekeeping step - set up and verify the file names for the output:
air16S_filtFs <- list.files(air16S_filtpathF, pattern="fastq.gz", full.names = TRUE)
length(air16S_filtFs) #165
air16S_filtRs <- list.files(air16S_filtpathR, pattern="fastq.gz", full.names = TRUE)
length(air16S_filtRs) #165
# which ones are we missing? (should have 164)
air16S_mappingFile <- read.delim("~/SRS_aeromicrobiome_2022/air16S_MappingFile", header= FALSE)
sampleNamesMapFile_air16S <- air16S_mappingFile$V2
air16S_filtFsSampleNames <- air16S_filtFs

air16S_filtFsSampleNames <-gsub("/data/winfreyc/SRS_aeromicrobiome_2022/air16S_02_filter/preprocessed_F/filtered/R1_","", air16S_filtFsSampleNames)
air16S_filtFsSampleNames <-gsub(".fastq.gz","", air16S_filtFsSampleNames)
setdiff(sampleNamesMapFile_air16S, air16S_filtFsSampleNames) #air_16S_PCR_NTC_1" is found in mapping file but not files

# Sample names in order
air16S_sample.names <- substring(basename(air16S_filtFs), regexpr("_", basename(air16S_filtFs)) + 1) # doesn't drop fastq.gz
air16S_sample.names <- gsub(".fastq.gz", "", air16S_sample.names)
air16S_sample.namesR <- substring(basename(air16S_filtRs), regexpr("_", basename(air16S_filtRs)) + 1) # doesn't drop fastq.gz
air16S_sample.namesR <- gsub(".fastq.gz", "", air16S_sample.namesR)

# Double check
if(!identical(air16S_sample.names, air16S_sample.namesR)) stop("Forward and reverse files do not match.")
names(air16S_filtFs) <- air16S_sample.names
names(air16S_filtRs) <- air16S_sample.names

# #### Learn the error rates
set.seed(100) # set seed to ensure that randomized steps are replicatable 

# Learn forward error rates (Notes: randomize default is FALSE)
air16S_errF <- learnErrors(air16S_filtFs, nbases = 1e8, multithread = TRUE, randomize = TRUE) 

# Learn reverse error rates
air16S_errR <- learnErrors(air16S_filtRs, nbases = 1e8, multithread = TRUE, randomize = TRUE)

# #### Plot Error Rates
# We want to make sure that the machine learning algorithm is learning the error rates properly. In the plots below, the red line represents what we should expect the learned error rates to look like for each of the 16 possible base transitions (A->A, A->C, A->G, etc.) and the black line and grey dots represent what the observed error rates are. If the black line and the red lines are very far off from each other, it may be a good idea to increase the ```nbases``` parameter. This alows the machine learning algorthim to train on a larger portion of your data and may help imporve the fit.

air16S_errF_plot <- plotErrors(air16S_errF, nominalQ = TRUE) 
air16S_errF_plot
# Error: Warning messages:
# 1: Transformation introduced infinite values in continuous y-axis 
# 2: Transformation introduced infinite values in continuous y-axis 

# Message from Ben Callahan (dada2 creator) on an online forum
# This isn't an error, just a message from the plotting function to let you know that there were some zero values 
# in the data plotted (which turn into infinities on the log-scale). That is completely expected, it results 
# from the fact that not every combination of error type (e.g. A->C) and quality score (e.g. 33) is observed 
# in your data which is normal.

air16S_errR_plot <- plotErrors(air16S_errR, nominalQ = TRUE)
air16S_errR_plot
#
# write to disk - saved May 24/24, 2023 (attempted overnight)
saveRDS(air16S_errF_plot, paste0(air16S_filtpathF, "/air16S_errF_plot.rds"))
saveRDS(air16S_errR_plot, paste0(air16S_filtpathR, "/air16S_errR_plot.rds"))

# #### Dereplication, sequence inference, and merging of paired-end reads
# In this part of the pipeline, dada2 will make decisions about assigning sequences to ASVs (called "sequence inference"). There is a major parameter option in the core function dada() that changes how samples are handled during sequence inference. The parameter ```pool = ``` can be set to: ```pool = FALSE``` (default), ```pool = TRUE```, or ```pool = psuedo```. For details on parameter choice, please see below, and further information on this blogpost [http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/](http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/), and explanation on the dada2 tutorial [https://benjjneb.github.io/dada2/pool.html](https://benjjneb.github.io/dada2/pool.html).

# # SAMPLES POOLED 
# For complex communities when you want to preserve rare taxa

# Dereplicate forward reads
air16S_derepF.p <- derepFastq(air16S_filtFs)
names(air16S_derepF.p) <- air16S_sample.names
# Infer sequences for forward reads
air16S_dadaF.p <- dada(air16S_derepF.p, err = air16S_errF, multithread = TRUE, pool = TRUE)
names(air16S_dadaF.p) <- air16S_sample.names

# Dereplicate reverse reads
air16S_derepR.p <- derepFastq(air16S_filtRs)
names(air16S_derepR.p) <- air16S_sample.names
# Infer sequences for reverse reads
air16S_dadaR.p <- dada(air16S_derepR.p, err = air16S_errR, multithread = TRUE, pool = TRUE)
names(air16S_dadaR.p) <- air16S_sample.names

# Merge reads together
air16S_mergers <- mergePairs(air16S_dadaF.p, air16S_derepF.p, air16S_dadaR.p, air16S_derepR.p)
head(air16S_mergers)

# #### Construct sequence table
air16S_seqtab <- makeSequenceTable(air16S_mergers)

# Save table as an r data object file
dir.create(air16S_table.fp)
saveRDS(air16S_seqtab, paste0(air16S_table.fp, "/air16S_seqtab.rds")) #attempted to save May 24/25, 2023

###################################################################################
##### PHYLLO ######

# #### Housekeeping step - set up and verify the file names for the output:
phyllo16S_filtFs <- list.files(phyllo16S_filtpathF, pattern="fastq.gz", full.names = TRUE)
length(phyllo16S_filtFs) #68
phyllo16S_filtRs <- list.files(phyllo16S_filtpathR, pattern="fastq.gz", full.names = TRUE)
length(phyllo16S_filtRs) #68
# which ones are we missing? (should have 164)
phyllo16S_mappingFile <- read.delim("/data/shared/MiSeq/11.11.2022_MixRun16S_try2/phyllo16S_mappingFile", header= FALSE)
sampleNamesMapFile_phyllo16S <- phyllo16S_mappingFile$V2
phyllo16S_filtFsSampleNames <- phyllo16S_filtFs
phyllo16S_filtFsSampleNames <-gsub("/data/winfreyc/SRS_aeromicrobiome_2022/phyllo16S_02_filter/preprocessed_F/filtered/R1_","", phyllo16S_filtFsSampleNames)
phyllo16S_filtFsSampleNames <-gsub(".fastq.gz","", phyllo16S_filtFsSampleNames)
setdiff(sampleNamesMapFile_phyllo16S, phyllo16S_filtFsSampleNames) # "phyllo_16S_PCR_NTC_1" "phyllo_16S_PCR_NTC_2" are not in filtFs

# Sample names in order
phyllo16S_sample.names <- substring(basename(phyllo16S_filtFs), regexpr("_", basename(phyllo16S_filtFs)) + 1) # doesn't drop fastq.gz
phyllo16S_sample.names <- gsub(".fastq.gz", "", phyllo16S_sample.names)
phyllo16S_sample.namesR <- substring(basename(phyllo16S_filtRs), regexpr("_", basename(phyllo16S_filtRs)) + 1) # doesn't drop fastq.gz
phyllo16S_sample.namesR <- gsub(".fastq.gz", "", phyllo16S_sample.namesR)

# Double check
if(!identical(phyllo16S_sample.names, phyllo16S_sample.namesR)) stop("Forward and reverse files do not match.")
names(phyllo16S_filtFs) <- phyllo16S_sample.names
names(phyllo16S_filtRs) <- phyllo16S_sample.names

# #### Learn the error rates
set.seed(100) # set seed to ensure that randomized steps are replicatable 

# Learn forward error rates (Notes: randomize default is FALSE)
phyllo16S_errF <- learnErrors(phyllo16S_filtFs, nbases = 1e8, multithread = TRUE, randomize = TRUE) 

# Learn reverse error rates
phyllo16S_errR <- learnErrors(phyllo16S_filtRs, nbases = 1e8, multithread = TRUE, randomize = TRUE)

# #### Plot Error Rates
# We want to make sure that the machine learning algorithm is learning the error rates properly. In the plots below, the red line represents what we should expect the learned error rates to look like for each of the 16 possible base transitions (A->A, A->C, A->G, etc.) and the black line and grey dots represent what the observed error rates are. If the black line and the red lines are very far off from each other, it may be a good idea to increase the ```nbases``` parameter. This alows the machine learning algorthim to train on a larger portion of your data and may help imporve the fit.

phyllo16S_errF_plot <- plotErrors(phyllo16S_errF, nominalQ = TRUE) 
phyllo16S_errF_plot
# Error: Warning messages:
# 1: Transformation introduced infinite values in continuous y-axis 
# 2: Transformation introduced infinite values in continuous y-axis 

# Message from Ben Callahan (dada2 creator) on an online forum
# This isn't an error, just a message from the plotting function to let you know that there were some zero values 
# in the data plotted (which turn into infinities on the log-scale). That is completely expected, it results 
# from the fact that not every combination of error type (e.g. A->C) and quality score (e.g. 33) is observed 
# in your data which is normal.

phyllo16S_errR_plot <- plotErrors(phyllo16S_errR, nominalQ = TRUE)
phyllo16S_errR_plot
#
# write to disk
saveRDS(phyllo16S_errF_plot, paste0(phyllo16S_filtpathF, "/phyllo16S_errF_plot.rds"))
saveRDS(phyllo16S_errR_plot, paste0(phyllo16S_filtpathR, "/phyllo16S_errR_plot.rds"))

# #### Dereplication, sequence inference, and merging of paired-end reads
# In this part of the pipeline, dada2 will make decisions about assigning sequences to ASVs (called "sequence inference"). There is a major parameter option in the core function dada() that changes how samples are handled during sequence inference. The parameter ```pool = ``` can be set to: ```pool = FALSE``` (default), ```pool = TRUE```, or ```pool = psuedo```. For details on parameter choice, please see below, and further information on this blogpost [http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/](http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/), and explanation on the dada2 tutorial [https://benjjneb.github.io/dada2/pool.html](https://benjjneb.github.io/dada2/pool.html).

# # SAMPLES POOLED 
# For complex communities when you want to preserve rare taxa

# Dereplicate forward reads
phyllo16S_derepF.p <- derepFastq(phyllo16S_filtFs)
names(phyllo16S_derepF.p) <- phyllo16S_sample.names
# Infer sequences for forward reads
phyllo16S_dadaF.p <- dada(phyllo16S_derepF.p, err = phyllo16S_errF, multithread = TRUE, pool = TRUE)
names(phyllo16S_dadaF.p) <- phyllo16S_sample.names

# Dereplicate reverse reads
phyllo16S_derepR.p <- derepFastq(phyllo16S_filtRs)
names(phyllo16S_derepR.p) <- phyllo16S_sample.names
# Infer sequences for reverse reads
phyllo16S_dadaR.p <- dada(phyllo16S_derepR.p, err = phyllo16S_errR, multithread = TRUE, pool = TRUE)
names(phyllo16S_dadaR.p) <- phyllo16S_sample.names

# Merge reads together
phyllo16S_mergers <- mergePairs(phyllo16S_dadaF.p, phyllo16S_derepF.p, phyllo16S_dadaR.p, phyllo16S_derepR.p)
head(phyllo16S_mergers)

# #### Construct sequence table
phyllo16S_seqtab <- makeSequenceTable(phyllo16S_mergers)

# Save table as an r data object file
dir.create(phyllo16S_table.fp)
saveRDS(phyllo16S_seqtab, paste0(phyllo16S_table.fp, "/phyllo16S_seqtab.rds"))


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
airSeqTab <- readRDS(paste0(air16S_table.fp, "/air16S_seqtab.rds"))
phylloSeqTab <- readRDS(paste0(phyllo16S_table.fp, "/phyllo16S_seqtab.rds"))
soilSeqtab <- readRDS("~/SRS_aeromicrobiome_2022/soilsfilteredPathF/soil16S_table/soil16S_seqtab.rds")

# Combine all these seqtabs
seqTab16S.all <- mergeSequenceTables(air16S_seqtab, phyllo16S_seqtab, soil16S_seqtab)

# Remove chimeras
all_16S_seqtab.nochim <- removeBimeraDenovo(seqTab16S.all, method="consensus", multithread=TRUE)

# Print percentage of our sequences that were not chimeric.
100*sum(all_16S_seqtab.nochim)/sum(seqTab16S.all) 

# Assign taxonomy
all_16S_tax <- assignTaxonomy(all_16S_seqtab.nochim, "/db_files/dada2/silva_nr99_v138.1_train_set.fa", tryRC = TRUE,
                             multithread=TRUE)
summary(all_16S_tax)  # 

# Write results to disk #saved May 25, 2023
saveRDS(all_16S_seqtab.nochim, paste0(project.fp, "/all_16S_seqtab_final.rds")) 
saveRDS(all_16S_tax, paste0(project.fp, "/all_16S_tax_final.rds")) 


# ### 4. Optional - FORMAT OUTPUT to obtain ASV IDs and repset, and input for mctoolsr
# For convenience sake, we will now rename our ASVs with numbers, output our 
# results as a traditional taxa table, and create a matrix with the representative
# sequences for each ASV. 

# Flip table
all16S_seqtab.t <- as.data.frame(t(all_16S_seqtab.nochim))

# Pull out ASV repset
all16S_rep_set_ASVs <- as.data.frame(rownames(all16S_seqtab.t))
all16S_rep_set_ASVs <- mutate(all16S_rep_set_ASVs, ASV_ID = 1:n())
all16S_rep_set_ASVs$ASV_ID <- sub("^", "ASV_", all16S_rep_set_ASVs$ASV_ID)
all16S_rep_set_ASVs$ASV <- all16S_rep_set_ASVs$`rownames(all16S_seqtab.t)` 
all16S_rep_set_ASVs$`rownames(all16S_seqtab.t)` <- NULL

# Add ASV numbers to table
rownames(all16S_seqtab.t) <- all16S_rep_set_ASVs$ASV_ID

# Add ASV numbers to taxonomy
all16S_taxonomy <- as.data.frame(all_16S_tax) 
all16S_taxonomy$ASV <- as.factor(rownames(all16S_taxonomy))
all16S_taxonomy <- merge(all16S_rep_set_ASVs, all16S_taxonomy, by = "ASV")
rownames(all16S_taxonomy) <- all16S_taxonomy$ASV_ID
all16S_taxonomy_for_mctoolsr <- unite(all16S_taxonomy, "all16S_taxonomy", 
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

# Arrange the all16S_taxonomy dataframe for the writeRepSetFasta function
all16S_taxonomy_for_fasta <- all16S_taxonomy %>%
  unite("TaxString", c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ASV_ID"), 
        sep = ";", remove = FALSE) %>%
  unite("name", c("ASV_ID", "TaxString"), 
        sep = " ", remove = TRUE) %>%
  select(ASV, name) %>%
  rename(seq = ASV)

# write fasta file- written April 11, 2023
writeRepSetFasta(all16S_taxonomy_for_fasta, "/data/winfreyc/SRS_aeromicrobiome_2022/all16S_repset.fasta")

# Merge all16S_taxonomy and table
all16S_seqtab_wTax <- merge(all16S_seqtab.t, all16S_taxonomy_for_mctoolsr, by = 0)
all16S_seqtab_wTax$ASV <- NULL 

# Set name of table in mctoolsr format and save
names(all16S_seqtab_wTax)[1] = "#ASV_ID"
write("#Exported for mctoolsr", "~/SRS_aeromicrobiome_2022/all16S_seqtab_wTax_mctoolsr.txt")
suppressWarnings(write.table(all16S_seqtab_wTax, "~/SRS_aeromicrobiome_2022/all16S_seqtab_wTax_mctoolsr.txt", sep = "\t", row.names = FALSE, append = TRUE))

# Also export files as .txt
write.table(all16S_seqtab.t, file = "~/SRS_aeromicrobiome_2022/all16S_seqtab_final.txt",
            sep = "\t", row.names = TRUE, col.names = NA)
write.table(all_16S_tax, file = "~/SRS_aeromicrobiome_2022/all16S_tax_final.txt", 
            sep = "\t", row.names = TRUE, col.names = NA)

# ### Summary of output files:
# 1. all16S_seqtab_final.txt - A tab-delimited sequence-by-sample (i.e. OTU) table 
# 2. all16S_tax_final.txt - a tab-demilimited file showing the relationship between ASVs, ASV IDs, and their taxonomy 
# 3. all16S_seqtab_wTax_mctoolsr.txt - a tab-delimited file with ASVs as rows, samples as columns and the final column showing the taxonomy of the ASV ID 
# 4. all16S_repset.fasta - a fasta file with the representative sequence of each ASV. Fasta headers are the ASV ID and taxonomy string.  
#

# 
# # NEXT STEPS:
# 