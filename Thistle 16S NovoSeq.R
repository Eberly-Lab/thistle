setwd("Analysis/amplicon/16S")

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

path <- "/Thistle_16S/fastq/00.RawData"

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

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2]) 

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

out <- filterAndTrim(fwd=file.path(fnFs), filt=file.path(filtFs),
                     rev=file.path(fnRs), filt.rev=file.path(filtRs), minLen = 50,
                     maxEE=c(2,4), truncQ=2, maxN=0, trimLeft=c(20, 20), trimRight = 5, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE, matchIDs=TRUE)
# Note: matchIDs=TRUE is needed to match forward and reverse reads with NovaSeq reads
head(out)

derepFs <- derepFastq(fnFs, verbose=TRUE, n = 1e+08) # Set n = 1e+03 if < 32 GB RAM
derepRs <- derepFastq(filtRs, verbose=TRUE, n = 1e+08) 
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

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

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, verbose = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, verbose = TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

seqtabAll <- makeSequenceTable(mergers) 

table(nchar(getSequences(seqtabAll)))
saveRDS(seqtabAll, "Sequence Table.RData")

#Remove chimeras
seqtabNoC <- removeBimeraDenovo(seqtabAll, multithread=TRUE, verbose=TRUE)
saveRDS(seqtabNoC, "Sequence Table No Chimaras.RData") #Save sequence table
seqtabNoC <- readRDS("Sequence Table No Chimaras.RData")

# Track number of reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtabNoC))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sampleNames
head(track)

# Plotting each sample's reads through the pipeline
track_plot <- track %>% 
  data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  gather(key = "Step", value = "Reads", -Sample) %>%
  mutate(Step = factor(Step, 
                       levels = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"))) %>%
  ggplot(aes(x = Step, y = Reads)) +
  geom_line(aes(group = Sample), alpha = 0.2) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0)) + 
  stat_summary(fun.y = median, geom = "line", group = 1, color = "steelblue", size = 1, alpha = 0.5) +
  stat_summary(fun.y = median, geom = "point", group = 1, color = "steelblue", size = 2, alpha = 0.5) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), 
               geom = "ribbon", group = 1, fill = "steelblue", alpha = 0.2) +
  geom_label(data = t(track_pct_avg[1:5]) %>% data.frame() %>% 
               rename(Percent = 1) %>%
               mutate(Step = c("filtered", "denoisedF", "denoisedR", "merged", "nonchim"),
                      Percent = paste(round(Percent, 2), "%")),
             aes(label = Percent), y = 1.1 * max(track[,2])) +
  geom_label(data = track_pct_avg[6] %>% data.frame() %>%
               rename(total = 1),
             aes(label = paste("Total\nRemaining:\n", round(track_pct_avg[1,6], 2), "%")), 
             y = mean(track[,6]), x = 6.5) +
  expand_limits(y = 1.1 * max(track[,2]), x = 7) +
  theme_classic()

track_plot

#Assign taxonomy

# Taxonomy assignment with Silva database with species included
fastaRef2 <- "/Users/Jed/Silva_138.1/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
taxTabSp <- assignTaxonomy(seqtabNoC, refFasta = fastaRef2, multithread=TRUE, verbose = TRUE)
saveRDS(taxTabSp, "taxTabSp.RDS")
taxTabSp <- readRDS("taxTabSp.RDS")

# Make phyloseq object
sd <- read.csv("sample-data-16S.csv", header = TRUE)
rownames(sd) <- sd[,1]
OTU <- otu_table(seqtabNoC, taxa_are_rows = FALSE)
TAX <- tax_table(taxTabSp)
sampledata <- sample_data(data.frame(sd))

ps <- phyloseq(OTU, TAX, sampledata)
ps

# Add refseqs for predictive functional profiling with PICRUSt
ps_ref <- add_refseq(ps, tag = "ASV")
ps_ref
saveRDS(ps_ref, "ps_ref_seq.RDS")

# Add ASV short names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(pss, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps <- pss

saveRDS(ps, "ps_16S.RDS")

# Add tree 
library("ape")
random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))

ps <- merge_phyloseq(ps, random_tree)
ps 
saveRDS(ps, "ps_16S_Nova.RDS")


