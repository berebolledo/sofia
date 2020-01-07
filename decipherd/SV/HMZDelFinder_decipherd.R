################################################################################

# Running HMZDelFinder
# https://github.com/berebolledo/HMZDelFinder

################################################################################

# Load packages
library(RCurl)
library(data.table)
library(gdata)
library(parallel)
library(Hmisc)
library(matrixStats)
library(DNAcopy)
library(GenomicRanges)
library(Rsubread)

################################################################################

# Directories

proj <- "/Users/boris/GoogleDrive/UDD/research/bioinformatics/SABIO/projects/"
wd <- paste0(proj, "01_DECIPHERD/baylor/LABWORK/")
setwd(wd)
workDir <- getwd()
mainDir <- paste0(workDir ,"/HMZDelFinder/")
if (!file.exists(mainDir)){dir.create(mainDir)}

dataDir <- paste0(mainDir, "data/" , sep="")
if (!file.exists(dataDir)){dir.create(dataDir)}

outputDir <- paste0(mainDir, "out/" , sep="")
if (!file.exists(outputDir)){dir.create(outputDir)}

plotsDir <- paste0(mainDir, "plots/" , sep="")
if (!file.exists(plotsDir)){dir.create(plotsDir)}
################################################################################

# Files
# BED file
# Tab delimited file without header and four columns:
# Chromosome Start Stop Gene symbol
# RPKM files
# Tab delimited file with a header and two columns:
#    count // the number of reads overlapping with each capture target
#    RPKM // the RPKM value for each capture target

bedFile <- paste0(dataDir, "refGene.hg19.sorted.nochr.4col.bed")
countsDir <- "baylor_counts_tpm"
rpkmFiles <- list.files(paste0(dataDir, countsDir ), ".txt")
rpkmFids <- gsub("_counts.txt", "", rpkmFiles)
rpkmPaths <- paste0(paste0(dataDir, countsDir, "/"), rpkmFiles)

outCalls <- paste0(outputDir, countsDir, "_hmzCalls.csv")
################################################################################

# Load HMZDelFinder source code
gitHub <- "https://raw.githubusercontent.com/"
Lupski <- "BCM-Lupskilab/HMZDelFinder/master/src/HMZDelFinder.R"
sourceCode <- paste0(gitHub, Lupski)
eval(expr = parse(text = getURL(sourceCode)))

# THRESHOLDS
# See description of HMZDelFinder function for details

is_cmg <- FALSE 		     # only for CMG project - otherwhise use FALSE
lowRPKMthreshold <- 0.65 # RPKM threshold  
maxFrequency <- 0.005	   # max frequncy of HMZ deletion; default =0.005
minAOHsize <- 1000		   # min AOH size
minAOHsig <- 0.45		     # min AOH signal threshold
mc.cores <- 4 				   # number of cores
vR_id<- "VR"				     # ID from VCF FORMAT variant reads
tR_id<- "DP"				     # ID from VCF FORMAT indicating the number total reads 
filter <- "PASS"		     # for other variant callers be  '.'

# running HMZDelFinder
results <- runHMZDelFinder(
  NULL,		    # vcfPaths - paths to VCF files for AOH analysis
  NULL,		    # vcfFids - sample identifiers corresponding to VCF files
  rpkmPaths, 	# paths to RPKM files 
  rpkmFids,	  # samples identifiers corresponding to RPKM files
  mc.cores,	  # number of CPU cores
  aohRDataOut,# temp file to store AOH data
  bedFile,	  # bed file with target 
  lowRPKMthreshold, #  RPKM threshold 
  minAOHsize, # min AOH size
  minAOHsig,	# min AOH signal threshold
  is_cmg,		  # flag used for CMG specific annotations
  vR_id, 		  # ID for 'the number of variants reads' in VCF FORMAT
  tR_id,		  # ID for 'the number of total reads' in VCF FORMAT
  filter)		  # only variants with this value for AOH analysis 


# saving results in csv files
write.csv(results$filteredCalls, outCalls, row.names=F)

# plotting deletions

plots <- function(i){
  plotDeletion (results$filteredCalls, i, results$bedOrdered, 
                results$rpkmDtOrdered,  lowRPKMthreshold, 
                plotsDir, mainText="")}

lapply(1:nrow(results$filteredCalls), plots)
