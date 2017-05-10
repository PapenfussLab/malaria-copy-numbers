## This reads the reference genome of P. falciparum in BioStrings format,
## reads bam files (first time) or binned read counts if available, 
## for the parent and resistant clones, then uses functions from 
## QDNAseq to calculate read depths and copynumbers.
## A 'mappability file' is required for correcting the counts
## Copied from malvacR to maldrugR Feb 2017 and modified for new problem and improved QDNASeq
## Jocelyn Sietsma Penington  
## 
## Current version of QDNAseq is 1.10.0

library("BSgenome")
library("QDNAseq")
library("Biobase")

source("file_paths.R")
alignDir <- file.path(baseDir, "alignment")

## For first use of new version of reference genome, need to forge the genome data package:
# library(devtools)
# forgeBSgenomeDataPkg(file.path(refDir, "BioStringsFormatSeedFile.txt"), 
#                      seqs_srcdir = refDir, destdir = refDir)
# install_local(path = file.path(refDir, 
#                         "BSgenome.Pfalciparum3D7.PlasmoDBv29"))

library("BSgenome.Pfalciparum3D7.PlasmoDBv29")
pfg <- BSgenome.Pfalciparum3D7.PlasmoDBv29
seqinfo(pfg)
str(pfg[[1]])
head(pfg[[1]])
## BSgenome for Pfalciparum3D7 successfully loaded

## Statement below was executed for 1k, then edited and re-run for 5k.
bin_in_kbases <- "10"
## It creates an annotated dataframe of bins, listing their start, end, GC percentage, etc.
if (file.exists(file.path(refDir, paste0("PfBins", bin_in_kbases, "kb.rds")))) {
  pfBins <- readRDS(file.path(refDir, paste0("PfBins", bin_in_kbases, "kb.rds")))
} else {
    pfBins <- createBins(pfg, as.numeric(bin_in_kbases))
    bwfile <- file.path(refDir, "Pfalciparum3D7_Genome50mer.bigWig")
    pfBins$mappability <- calculateMappability(pfBins, 
                                        bigWigFile=bwfile, 
                                        bigWigAverageOverBed='/usr/local/bioinf/bin/bigWigAverageOverBed',
                                        chrPrefix='')
  ## Use a blacklist BEDfile to mask telomeres and centromeres. 
    pfBins$blacklist <- calculateBlacklist(pfBins, 
                                             bedFiles= file.path(
                                               "/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects",
                                               "reference_genomes/plasmodium/", 
                                               '12.0/CentromereTelomereRegions.bed') )
    pfBins <- AnnotatedDataFrame(pfBins,
                                 varMetadata=data.frame(labelDescription=c(
                                   'Chromosome name',
                                   'Base pair start position',
                                   'Base pair end position',
                                   'Percentage of non-N nucleotides (of full bin size)',
                                   'Percentage of C and G nucleotides (of non-N nucleotides)',
                                   'Average mappability of 50mers',
                                   'Percentage overlap of bin with blacklist region'),
                                   row.names=colnames(pfBins)) )
  saveRDS(pfBins, file.path(refDir,   paste0("PfBins", bin_in_kbases, "kb.rds")))
}

## If bin counts calculated previously, read saved data.
## Otherwise, using the bins above, calculate read depths from files
if (file.exists(
  file.path(cnDir, paste0("Counts", bin_in_kbases, "k_mergedS1S2_Set1.rds"))
  )) {
  counts_set1 <- readRDS(file.path(cnDir, paste0("Counts", bin_in_kbases, "k_mergedS1S2_Set1.rds")))
    } else {
    pfBins <- readRDS(file.path(refDir, paste0("PfBins", bin_in_kbases, "kb.rds")))
    counts_set1 <- binReadCounts(pfBins, 
                            bamfiles=c(file.path(alignDir,'3D7-merge-B3_S1-C5_S2.bam'), 
                                       file.path(alignDir,'754-1-D2_S3_recalQ.bam'), 
                                       file.path(alignDir,'754-1-E8_S4_recalQ.bam'), 
                                       file.path(alignDir,'754-1-G7_S5_recalQ.bam')
                                       ),
                            bamnames=c('ParentsS1S2', 'S3', 'S4', 'S5'))
  
  saveRDS(counts_set1, file.path(cnDir, paste0("Counts", bin_in_kbases, "k_mergedS1S2_Set1.rds")))  
    }

if (file.exists(
  file.path(cnDir, paste0("Counts", bin_in_kbases, "k_mergedS1S2_Set2.rds")))
  ) {
  counts_set2 <- readRDS(
    file.path(cnDir, paste0("Counts", bin_in_kbases, "k_mergedS1S2_Set2.rds"))
    ) 
  } else {
  pfBins <- readRDS(file.path(refDir, paste0("PfBins", bin_in_kbases, "kb.rds")))
  counts_set2 <- binReadCounts(pfBins, 
                                 bamfiles=c(file.path(alignDir,'3D7-merge-B3_S1-C5_S2.bam'), 
                                            file.path(alignDir,'754-2-B3_S6_recalQ.bam'), 
                                            file.path(alignDir,'754-2-C11_S7_recalQ.bam'), 
                                            file.path(alignDir,'754-2-F5_S8_recalQ.bam'), 
                                            file.path(alignDir,'754-2-G7_S9_recalQ.bam')
                                 ),
                                 bamnames=c('ParentsS1S2', 'S6', 'S7', 'S8', 'S9'))
  
  saveRDS(counts_set2, 
          file.path(cnDir, paste0("Counts", bin_in_kbases, "k_mergedS1S2_Set2.rds")))  
  }

## Plot read counts for set 1 (Default is log2 with log2(0) set to -1022):
plot(counts_set1)

## Filtering with min-mappability=50% , blacklist and exclude non-nuclear
## I didn't calculate residuals - see version in malvacR for code if wanted
counts1Filtered <- applyFilters(counts_set1, mappability=50, blacklist=TRUE, 
                    chromosomes=c('Pf3D7_API_v3', 'Pf_M76611'))
counts2Filtered <- applyFilters(counts_set2, mappability=50, blacklist=TRUE, 
                                  chromosomes=c('Pf3D7_API_v3', 'Pf_M76611'))

## 3-d contour plot of read depth vs GC and mappability
isobarPlot(counts_set1)

## Estimate the correction for GC content and mappability. 
## This calculates a loess fit to the counts (as a function of GC and mappability)
## and adds it to the data structure
counts1Corrected <- estimateCorrection(counts1Filtered)
noisePlot(counts1Corrected) 
counts2Corrected <- estimateCorrection(counts2Filtered)

## Use the values added by estimateCorrection to correct the data:
## Copy numbers are counts divided by 'fit'  
copyNums_set1 <- correctBins(counts1Corrected) 
copyNums_set2 <- correctBins(counts2Corrected) 
## Scale copy numbers by median (default)
copyNums1Normed <- normalizeBins(copyNums_set1)
copyNums2Normed <- normalizeBins(copyNums_set2)
plot(copyNums2Normed)

# copyNums1k1Smoothed <- smoothOutlierBins(copyNums1k1Normed, logTransform = FALSE)
# copyNums1k2Smoothed <- smoothOutlierBins(copyNums1k2Normed, logTransform = FALSE)
# 
# exportBins(copyNums1k1Smoothed, 
#            file=file.path(cnDir, "Set1_1kFiltCorrectNormSmoothBins.txt") )
## By default, exported values are log2 of copynumber values, with log2(0) set to -1022

Set1ScaledByParents <- compareToReference(copyNums1Normed, c(FALSE, 1, 1, 1))
plot(Set1ScaledByParents)
Set2ScaledByParents <- compareToReference(copyNums2Normed, c(FALSE, 1, 1, 1, 1))
plot(Set2ScaledByParents)
#################################################################################

## manipulating results in ggplot2

convertQDNAtoDF <- function(qobject) {
  if ("counts" %in% assayDataElementNames(qobject) ) {
    scoreDF <- as.data.frame( assayDataElement(qobject, "counts"))
  } else if ("copynumber" %in% assayDataElementNames(qobject) ) {
    scoreDF <- as.data.frame( assayDataElement(qobject, "copynumber"))
  }
  scoreDF$chrom <- sapply(strsplit(rownames(scoreDF), ':'), function(x) (unlist(x)[1]) )
  scoreDF$range <- sapply(strsplit(rownames(scoreDF), ':'), function(x) (unlist(x)[2]) )
  scoreDF$start <- sapply(strsplit(scoreDF$range, '-'), function(x) (as.numeric(unlist(x)[1]) ) )
  scoreDF$end <- sapply(strsplit(scoreDF$range, '-'), function(x) (as.numeric(unlist(x)[2]) ) )
  
  return(scoreDF)
}

library(ggplot2)
theme_set(theme_bw())

plotScaledCounts <- function(binCounts, orgCol, bin_size) {
  ## Plot nuclear chromosomes in 2 rows
  library(gridExtra)  # to arrange the 14 chromosomes in 2 rows
  library(scales)     # Need the scales package for log2 transform
  row1 <- paste0('Pf3D7_0', 1:9, '_v3')
  nonNuc <- c('Pf3D7_API_v3', 'Pf_M76611')
  binCounts <- na.omit(subset(binCounts, !(binCounts$chrom %in% nonNuc)))
  binCounts$score <- binCounts[, orgCol]
  maxCN  <- max(binCounts$score)
#  minCN <- min(binCounts[which(binCounts$score > 0), 'score'])
  minCN <- 2^-5
  ## Changing to fixed minimum. QDNASeq uses 2^-3, I will use 2^-5
  ## Previous y-min was 0.9*minCN
  binCounts$score <- pmax(binCounts$score, minCN)
  p1 <- ggplot(binCounts[ which(binCounts$chrom %in% row1), ] , 
               aes(end/1E5, score)) + 
    geom_point(size=0.5) + 
    facet_grid(~chrom, space="free_x", scales="free_x"
    ) + 
    scale_x_continuous(breaks=seq(0,max(binCounts$end/1E5),5)) +
    labs(title = paste(orgCol, "scores for", bin_size, "kb bins"), x='') +
    # scale_y_continuous(trans=log2_trans(), limits=c(minCN, 1.1*maxCN), 
    #                 breaks = trans_breaks("log2", function(x) 2^x),
    #                 labels = trans_format("log2", math_format(2^.x))) ## removing log-transform
    scale_y_continuous(limits=c(minCN, 1.1*maxCN))
  
  p2 <- ggplot(na.omit(binCounts[ which( ! binCounts$chrom %in% row1), ] ), 
               aes(end/1E5, score)) + 
    geom_point(size=0.5) + 
    facet_grid(~chrom, space="free_x", scales="free_x"
    ) + 
    scale_x_continuous(breaks=seq(0,max(binCounts$end/1E5),5)) +
    labs(x="bin-end location in chromosome (bp x 100,000)") +
    # scale_y_continuous(trans=log2_trans(), limits=c(minCN, 1.1*maxCN), 
    #                breaks = trans_breaks("log2", function(x) 2^x),
    #                labels = trans_format("log2", math_format(2^.x)))  ## removing log-transform
    scale_y_continuous(limits=c(minCN, 1.1*maxCN))
  
  arrangeGrob(p1,p2, ncol=1)
  
}

Set1scaled_df <- convertQDNAtoDF(Set1ScaledByParents)
saveRDS(Set1scaled_df, 
        file.path(cnDir, paste0("CN_filt_correct_norm_compare_df_", bin_in_kbases, "k_754.1.rds")))  

p <- plotScaledCounts(Set1scaled_df, 'S3 vs. ParentsS1S2', bin_in_kbases)
plot(p)
ggsave(filename = 
         file.path(cnDir, 
                   paste0("S3vsS1S2_", bin_in_kbases, "k_CN_no_smoothing.pdf"))
       , units = "mm", width = 297, height = 210, p)
p <- plotScaledCounts(Set1scaled_df, 'S4 vs. ParentsS1S2', bin_in_kbases)
ggsave(filename = 
         file.path(cnDir, 
                   paste0("S4vsS1S2_", bin_in_kbases, "k_CN_no_smoothing.pdf"))
       , units = "mm", width = 297, height = 210, p)
p <- plotScaledCounts(Set1scaled_df, 'S5 vs. ParentsS1S2', bin_in_kbases)
ggsave(filename = 
         file.path(cnDir, 
                   paste0("S5vsS1S2_", bin_in_kbases, "k_CN_no_smoothing.pdf"))
       , units = "mm", width = 297, height = 210, p)

Set2scaled_df <- convertQDNAtoDF(Set2ScaledByParents)
saveRDS(Set2scaled_df, 
        file.path(cnDir, paste0("CN_filt_correct_norm_compare_df_", bin_in_kbases, "k_754.2.rds")))  
p <- plotScaledCounts(Set2scaled_df, 'S6 vs. ParentsS1S2', bin_in_kbases)
ggsave(filename = 
         file.path(cnDir, 
                   paste0("S6vsS1S2_", bin_in_kbases, "k_CN_no_smoothing.pdf"))
       , units = "mm", width = 297, height = 210
       , p)
p <- plotScaledCounts(Set2scaled_df, 'S7 vs. ParentsS1S2', bin_in_kbases)
ggsave(filename = 
         file.path(cnDir, 
                   paste0("S7vsS1S2_", bin_in_kbases, "k_CN_no_smoothing.pdf"))
       , units = "mm", width = 297, height = 210
       , p)
p <- plotScaledCounts(Set2scaled_df, 'S8 vs. ParentsS1S2', bin_in_kbases)
ggsave(filename = 
         file.path(cnDir, 
                   paste0("S8vsS1S2_", bin_in_kbases, "k_CN_no_smoothing.pdf"))
       , units = "mm", width = 297, height = 210
       , p)
p <- plotScaledCounts(Set2scaled_df, 'S9 vs. ParentsS1S2', bin_in_kbases)
ggsave(filename = 
         file.path(cnDir, 
                   paste0("S9vsS1S2_", bin_in_kbases, "k_CN_no_smoothing.pdf"))
       , units = "mm", width = 297, height = 210
       , p)

###
## Zoom in on Chrom 8
## facet_wrap the samples in a set
require(reshape2)
set1_chrom8 <- melt(Set1scaled_df[Set1scaled_df$chrom=="Pf3D7_08_v3", ],
                   id.vars = c("chrom", "range", "start", "end"), 
                   variable.name = "sample", value.name = "copynum")
set1_chrom8$copynum[is.na(set1_chrom8$copynum)] <- 1  # set undefined copynum (probably 0/0) to 1
plotSet1_chrom8 <- ggplot(data = set1_chrom8, 
                          aes(x=end/1E5, y = copynum)) + 
  geom_point(size=0.5) + facet_wrap(~sample) + 
  labs(title = paste("Chrom 8 in", bin_in_kbases, "kb bins"), x='bin end position in 100kb') + 
  scale_y_continuous(breaks=seq(0,max(set1_chrom8$copynum)+1,1))
plotSet1_chrom8

set1_chrom8_half <- set1_chrom8[set1_chrom8$end < 600500, ]
ggplot(data = set1_chrom8_half, 
       aes(x=end/1E5, y = copynum)) + 
  geom_point(size=0.5) + facet_wrap(~sample) + 
  labs(title = paste("Chrom 8:1-600kb in", bin_in_kbases, "kb bins"), x='bin end position in 100kb') + 
  scale_y_continuous(breaks=seq(0,max(set1_chrom8$copynum)+1,1))

ggsave(filename = 
         file.path(cnDir, 
                   paste0("set1_chrom8_", bin_in_kbases, "k_CN_1to600.pdf"))
       , units = "mm", width = 297, height = 210)

set2_chrom8 <- melt(Set2scaled_df[Set2scaled_df$chrom=="Pf3D7_08_v3", ],
                    id.vars = c("chrom", "range", "start", "end"), 
                    variable.name = "sample", value.name = "copynum")
set2_chrom8$copynum[is.na(set2_chrom8$copynum)] <- 1  # set undefined copynum (probably 0/0) to 1
plotSet2_chrom8 <- ggplot(data = set2_chrom8, 
                          aes(x=end/1E5, y = copynum)) + 
  geom_point(size=0.5) + facet_wrap(~sample) + 
  labs(title = paste("Chrom 8 in", bin_in_kbases, "kb bins"), x='bin end position in 100kb') + 
  scale_y_continuous(breaks=seq(0,max(set2_chrom8$copynum)+1,1))
plotSet2_chrom8
ggsave(filename = 
         file.path(cnDir, 
                   paste0("set1_chrom8_", bin_in_kbases, "k_CN_no_smoothing.pdf"))
       , units = "mm", width = 297, height = 210)

set2_chrom8_half <- set2_chrom8[set2_chrom8$end < 600500, ]
ggplot(data = set2_chrom8_half, 
       aes(x=end/1E5, y = copynum)) + 
  geom_point(size=0.5) + facet_wrap(~sample) + 
  labs(title = paste("Chrom 8:1-600kb in", bin_in_kbases, "kb bins"), x='bin end position in 100kb') + 
  scale_y_continuous(breaks=seq(0,max(set2_chrom8$copynum)+1,1))

ggsave(filename = 
         file.path(cnDir, 
                   paste0("set2_chrom8_", bin_in_kbases, "k_CN_1to600.pdf"))
       , units = "mm", width = 297, height = 210)
#######
## Plot the parent as well
set2norm_df <- convertQDNAtoDF(copyNums2Normed)
p <- plotScaledCounts(set2norm_df, 'ParentsS1S2', as.numeric(bin_in_kbases))
plot(p)
ggsave(filename = 
         file.path(cnDir, 
                   paste0("S1S2_", bin_in_kbases, "k_CN_normed.pdf"))
       , units = "mm", width = 297, height = 210, p)

###############################################################################
## Segmentation from example in QDNAseq intro section 1.3. Not complete or useful

copyNumbersSegmented <- segmentBins(copyNums5k1Smoothed, transformFun="sqrt")

copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
plot(copyNumbersSegmented)
