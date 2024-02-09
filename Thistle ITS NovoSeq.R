setwd("Analysis/amplicon/ITS")

rm(list = ls()) # Clears data from global environment


library(ggplot2)
library(knitr)
library(gridExtra)
library(dada2)
library(phyloseq)
library(ggpubr)
library(lme4)
library(dplyr)
library(plyr)
library(vegan)
library(microbiomeutilities)
library(microbiome)
library(magrittr)
set.seed(1234)


# Set path to fastq files

path <- "/Thistle_ITS/fastq/00.RawData"

# Sort ensures forward/reverse reads are in same order
# Novogene sequencing results include raw and processed fastq files. The following uses the raw only
fnFs <- sort(list.files(path, pattern="raw_1.fastq", recursive = TRUE)) # Use recursive if fastq files are in subdirectories
fnRs <- sort(list.files(path, pattern="raw_2.fastq", recursive = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
sampleNames <- sapply(strsplit(sampleNames, "/"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

# Specify the full path to the fnFs and fnRs
out <- filterAndTrim(fwd=file.path(fnFs), filt=file.path(filtFs),
                     rev=file.path(fnRs), filt.rev=file.path(filtRs), maxN = 0, multithread = TRUE)
head(out)


# Primer removal with cutadapt
FWD <- "GCATCGATGAAGAACGCAGC"
REV <- "TCCTCCGCTTATTGATATGC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

# NOTE
# cutadapt cannot handle spaces in path name. If sequences are on a drive with spaces in name
# Move filtered sequences to a different directory (e.g., desktop) to run cut adapt.
cutadapt <- "mambaforge/envs/cutadaptenv/bin/cutadapt"

filt_path <- "sequencing/filtered"
path <- "sequencing"

fnFs.filtN <- sort(list.files(filt_path, pattern="F_filt.fastq", recursive = TRUE)) # Use recursive if fastq files are in subdirectories
fnRs.filtN <- sort(list.files(filt_path, pattern="R_filt.fastq", recursive = TRUE))

fnFs.filtN <- file.path(filt_path, fnFs.filtN)
fnRs.filtN <- file.path(filt_path, fnRs.filtN)

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, paste0(sampleNames, "_F_cut.fastq.gz"))
fnRs.cut <- file.path(path.cut, paste0(sampleNames, "_R_cut.fastq.gz"))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
       primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
       fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Move cutadapt folder back to the initial working directory
# Continue with dada2 workflow

path.cut <- "/Thistle_2023/Thistle_ITS/cutadapt"

cutFs <- sort(list.files(path.cut, pattern = "_F_cut.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R_cut.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(cutFs, "_F"), `[`, 1)
sampleNames <- sapply(strsplit(sampleNames, "/"), `[`, 8)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path.cut, cutFs)
fnRs <- file.path(path.cut, cutRs)

# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2]) 

miseq_path <- "/Thistle_2023/Thistle_ITS"
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))


out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(4, 10), truncQ = 2, verbose = TRUE,
                     minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE, matchIDs = TRUE) 
# Note: matchIDs=TRUE is needed to match forward and reverse reads with NovaSeq reads
head(out)

# Set names to match sample names
names(filtFs) <- sampleNames
names(filtRs) <- sampleNames

# Run NovaSeq Error Function.R before the following
errF <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e8,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)

errR <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e8,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, verbose = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, verbose = TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)

seqtabAll <- makeSequenceTable(mergers) 
table(nchar(getSequences(seqtabAll)))
saveRDS(seqtabAll, "Sequence Table ITS.RData")

#Remove chimeras
seqtabNoC <- removeBimeraDenovo(seqtabAll, multithread=TRUE, verbose=TRUE)
saveRDS(seqtabNoC, "Sequence Table ITS No Chimaras.RData") #Save sequence table
seqtabNoC <- readRDS("Sequence Table ITS No Chimaras.RData")

# Track number of reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtabNoC))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sampleNames
head(track)

#Assign taxonomy
fastaRef <- "/Unite/sh_general_release_dynamic_10.05.2021.fasta"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
saveRDS(taxTab, "taxTab ITS.RData")
taxTab <- readRDS("taxTab ITS.RData")

# Make phyloseq object
sd <- read.csv("sample-data-ITS.csv", header = TRUE)
rownames(sd) <- sd[,1]
OTU <- otu_table(seqtabNoC, taxa_are_rows = FALSE)
TAX <- tax_table(taxTab)
sampledata <- sample_data(data.frame(sd))

ps <- phyloseq(OTU, TAX, sampledata)
ps

# Add ASV short names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Add tree 
library("ape")
random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))

ps <- merge_phyloseq(ps, random_tree)
ps 
saveRDS(ps, "ps_ITS_Nova.RDS")

