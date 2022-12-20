# Set working directory for this test
setwd("~/workspace/Incliva_wkspace/TTV_estela")

### Load packages
library(dada2) # Info on sequence processing
library(DECIPHER) # Sequence alignment
library(phangorn) # Phylogenetic tree generation
library(ggplot2)
library(phyloseq)

# Specify path with data files
path = "./raw_reads" # Files must be unzipped (We have fastq R1 and R2)
list.files(path)

# # Sort the files to ensure forward/reverse are in the same order
fnFs = sort(list.files(path, pattern = '_R1_merged.fastq'))
fnRs = sort(list.files(path, pattern = '_R2_merged.fastq'))

# Extract sample names, assuming filenames have format: SAMPLENAME_Rx_001.fastq
sample.names = sapply(strsplit(fnFs, "_"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs = file.path(path, fnFs)
fnRs = file.path(path, fnRs)

# Plot quality scores and determine trimming cutoffs (here we are looking at the
# first two files)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# We want to filter so we need to set a path to store them
filt_path = file.path(path, "filtered") # Place the filtered files in filtered/ directory
dir.create(filt_path)

# Rename filtered files
filtFs = file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs = file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Quality filtering and trimming
# truncLen corresponds to the number of the base in which the quality drops below
# 30, so we will be trimming bases from 240 onwards in the Fwd and 160 onwards in Rvs

# Keep in mind how long is our amplicon and how much overlap we need, since
# DADA2 needs 20nt by default
# For V4 normally we have almost complete overlap and can trim reads aggressively
# maxN = 0 means that we dont want any ambiguous bases allowed
# maxEE, it is the maximum number of estimated errors allowed for an individual read,
# reads with more EE than this number will be discarded
# truncQ = truncates reads at the first instance of a Q score less than or equal to 
# the value specified (so a double check on quality)

out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(130,130),
                    maxN=0, maxEE = c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=TRUE)
head(out)

#### Since some files did not have reads that passed quality control, we have 
# to call for files again
filtFs = file.path(filt_path, sort(list.files(filt_path, pattern = '_F_filt.fastq.gz')))
filtRs = file.path(filt_path, sort(list.files(filt_path, pattern = '_R_filt.fastq.gz')))

# Estimate the error model for DADA2 algorithm
# Remember that every batch of sequencing will have a different error rate
# This alternates error rate estimation and sample composition inference until
# they converge at a consistent solution

errF = learnErrors(filtFs, multithread=TRUE)
errR = learnErrors(filtRs, multithread = TRUE)

# Plot error rates for all possible base transitions as a function of quality score
# Black line are observed error rates
# Red line are expected error rate under the nominal definition of the Q-value
# Remember that frequency of errors decrease as quality score increases
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
