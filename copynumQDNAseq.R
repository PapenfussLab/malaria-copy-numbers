## This reads the reference genome of P. falciparum in BioStrings format,
## reads bam files (first time) or binned read counts if available, 
## for the parent and resistant clones, then uses functions from 
## QDNAseq to calculate read depths and copynumbers.
## A 'mappability file' is required for correcting the counts

## Jocelyn Sietsma Penington  
## 
## Current version of QDNAseq is 1.10.0

## Modifying for github to make slightly less specific, 
## and include plot function not used with latest. 7 Nov 2017

library("BSgenome")
library("QDNAseq")
library("Biobase")

source("file_paths.R")
# defines currDir, refDir, cnDir 

alignDir <- file.path(currDir, "alignment")

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

## Statement below has been executed for 10k, 5k and 1k.
bin_in_kbases <- "5"
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
                                             bedFiles= file.path(refDir, "..",
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
######################################################
### Binned counts for samples
## If bin counts calculated previously, read saved data.
## Otherwise, using the bins above, calculate read depths from files
## Example file names for samples, by clonal sets:-
parentbam <- '3D7-merge-B2_S1-F4_S4.bam'
parentname <- 'mergedS1toS4'

strain <- 'E'
sampleL <- paste0('S', as.character(c(15:18)))

strainbamL <- paste0(strain, '_', sampleL, '_nodup.bam')
bamnameL <- c(parentname, sampleL)

countfilen <- paste0("Counts", bin_in_kbases, "k_", parentname, "_Strain", strain, ".rds")
if (file.exists(
  file.path(cnDir, countfilen))
  ) {
  countsinbins <- readRDS(file.path(cnDir, countfilen)) 
  } else {
  pfBins <- readRDS(file.path(refDir, paste0("PfBins", bin_in_kbases, "kb.rds")))
  countsinbins <- binReadCounts(pfBins, 
                                 bamfiles=file.path(alignDir,c(parentbam, strainbamL)), 
                                 bamnames=bamnameL)
  
  saveRDS(countsinbins, 
          file.path(cnDir, countfilen))  
  }

## Plot read counts (Default is log2 with log2(0) set to -1022):
plot(countsinbins)

## Filtering with min-mappability=50% , blacklist and exclude non-nuclear
## I didn't calculate residuals - see version in malvacR for code if wanted
countsFiltered <- applyFilters(countsinbins, mappability=50, blacklist=TRUE, 
                    chromosomes=c('Pf3D7_API_v3', 'Pf_M76611'))

## Estimate the correction for GC content and mappability. 
## This calculates a loess fit to the counts (as a function of GC and mappability)
## and adds it to the data structure
countsCorrected <- estimateCorrection(countsFiltered)
# noisePlot(countsCorrected) 

## Use the values added by estimateCorrection to correct the data:
## Copy numbers are counts divided by 'fit'  
copyNums <- correctBins(countsCorrected) 
# copyNums_set2 <- correctBins(counts2Corrected) 
## Scale copy numbers by median (default)
copyNumsNormed <- normalizeBins(copyNums)
# copyNums2Normed <- normalizeBins(copyNums_set2)
# plot(copyNumsANormed)

StrainScaledByParents <- compareToReference(copyNumsNormed, 
                                            c(FALSE, rep(1, times=length(strainbamL)))
                                           )
plot(StrainScaledByParents)

## manipulating results in ggplot2 requires converting to a data-frame ##
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

scaled_df <- convertQDNAtoDF(StrainScaledByParents)
saveRDS(scaled_df, 
        file.path(cnDir, paste0("CN_filt_correct_norm_compare_df_", bin_in_kbases, 
                                strain, ".rds"))) 

#################################################################################
### Making plots with ggplot ###

pdf.options(useDingbats=FALSE)  
## This should remove need for 'useDingbats=FALSE' in individual ggsave commands.
# Alternative solution is to edit fonts used by Illustrator
library(ggplot2)
theme_set(theme_bw())
library(gridExtra)  # to arrange the 14 chromosomes in 2 rows
library(scales)     # Need the scales package for log2 transform
require(reshape2)   # for melt

plotDir <- cnDir  ## option to change for final version

## Functions to draw various plot styles ###
nonNuc <- c('Pf3D7_API_v3', 'Pf_M76611')

plotScaledCounts <- function(binCounts, orgCol, bin_size) {
  ## Plot nuclear chromosomes in 2 rows
  row1 <- paste0('Pf3D7_0', 1:9, '_v3')
  binCounts <- na.omit(subset(binCounts, !(binCounts$chrom %in% nonNuc)))
  binCounts$score <- binCounts[, orgCol]
  maxCN  <- max(binCounts$score)
  minCN <- 2^-5
  ## Changing to fixed minimum. QDNASeq uses 2^-3, I will use 2^-5
  ## Previous y-min was 0.9 * minCN where minCN <- min(binCounts[which(binCounts$score > 0), 'score'])
  binCounts$score <- pmax(binCounts$score, minCN)
  p1 <- ggplot(binCounts[ which(binCounts$chrom %in% row1), ] , 
               aes(end/1E5, score)) + 
    geom_point(size=0.5) + 
    facet_grid(~chrom, space="free_x", scales="free_x"
    ) + 
    scale_x_continuous(breaks=seq(0,max(binCounts$end/1E5),5)) +
    labs(title = paste(orgCol, "scores for", bin_size, "kb bins"), x='') +
    scale_y_continuous(limits=c(minCN, 1.1*maxCN))
  
  p2 <- ggplot(na.omit(binCounts[ which( ! binCounts$chrom %in% row1), ] ), 
               aes(end/1E5, score)) + 
    geom_point(size=0.5) + 
    facet_grid(~chrom, space="free_x", scales="free_x"
    ) + 
    scale_x_continuous(breaks=seq(0,max(binCounts$end/1E5),5)) +
    labs(x="bin-end location in chromosome (bp x 100,000)") +
    scale_y_continuous(limits=c(minCN, 1.1*maxCN))
  
  arrangeGrob(p1,p2, ncol=1)
  
}

plotLogRow <- function(binCounts, orgCol, bin_size, maxCN, minCN) {
  ## Plot nuclear chromosomes in single row, log-y scale
  ## ymax and ymin can be passed in or calculated from data
  binCounts <- na.omit(subset(binCounts, !(binCounts$chrom %in% nonNuc)))
  binCounts$score <- binCounts[, orgCol]
  binCounts$shortChrom <- sapply(binCounts$chrom, function(x){strsplit(x, split='_')[[1]][2]})
  if (missing(maxCN)) maxCN <-  max(binCounts$score)
  if (missing(minCN)) minCN <-  min(binCounts[which(binCounts$score > 0), 'score'])
  binCounts$score <- pmax(binCounts$score, minCN)
  ggplot(binCounts , aes(end/1E5, score)) + 
    geom_point(size=0.5) + 
    facet_grid(~shortChrom, space="free_x", scales="free_x"
    ) + 
    scale_x_continuous(breaks=seq(0,max(binCounts$end/1E5),5)) +
    labs(title = paste(orgCol, "scores for", bin_size, "kb bins"), x='') +
    scale_y_continuous(trans=log2_trans(), limits=c(minCN, 1.1*maxCN),
                    breaks = trans_breaks("log2", function(x) 2^x),
                    labels = trans_format("log2", math_format(2^.x))) ## removing log-transform
}
  
plotWholeCN <- function(bin_df, bin_size){
  # plot all copynumbers in a single panel, coloured by Sample
  nonNuc <- c('API', 'Pf_M76611')
  bin_df$chrom <- gsub("Pf3D7_|_v3", "",bin_df$chrom) # shorten chromosome names
  bin_df <- na.omit(subset(bin_df, !(bin_df$chrom %in% nonNuc)))
  bin_df <- reshape2::melt(bin_df, id.vars = c("chrom", "range", "start", "end"),
                           variable.name = "sample", value.name = "copynum")
  bin_df$sample <-  sub(relativeText, "", bin_df$sample)
  bin_df$pos <- paste(bin_df$chrom, bin_df$start, sep=":")
  chromlist <- unique(bin_df$chrom)
  # counting on all chromosomes having a value for copy number with start==100001, for x-axis labelling
  ymax <- ceiling(max(bin_df$copynum)) 
  ggplot(bin_df, aes(x=pos, y=copynum, colour=sample)) + 
    geom_point(position=position_dodge(width = 0.5), size=0.5) + 
    scale_colour_brewer(palette="Set1") +
    scale_x_discrete(breaks=paste(chromlist, "100001", sep=":"), 
                     labels=chromlist) + 
    scale_y_continuous(limits=c(0, ymax),
                       breaks=seq(0, ymax, 2)) + 
    labs(title= paste(bin_size,"kb copy numbers"), x= "Chromosome", y="Relative copy numbers")
}

plotZoomedROI <- function(bin_df, chrom_2ch, startkb, endkb) {
  bin_df$chrom <- gsub("Pf3D7_|_v3", "",bin_df$chrom) # shorten chromosome names
  chromOI <- melt(bin_df[bin_df$chrom==chrom_2ch, ],
                 id.vars = c("chrom", "range", "start", "end"), 
                 variable.name = "sample", value.name = "copynum")
  chromOI$sample <-  sub(relativeText, "", chromOI$sample)
  roi <- na.omit(chromOI[which(chromOI$start > startkb*1000 & chromOI$end < endkb*1000), ])
  ymax <- ceiling(max(roi$copynum)) 
  ggplot(roi, aes(x = end/1000, y = copynum)) + geom_point(aes(color = sample)) + 
    facet_grid(sample ~ ., scales = "free_y") + theme(legend.position = "none") + 
    scale_colour_brewer(palette="Dark2") + 
    scale_y_continuous(limits=c(0, ymax),
                       breaks=seq(0, ymax, 2)) + 
    labs( x= "Position in chromosome ",chrom_2ch," (kb)", y="Relative copy numbers")
}

#### end of plotting functions ####

### Load saved data frame from part 1 if needed ###
## If not in environment already, will need to define bin_in_kbases, strain and sampleL

scaled_df <- readRDS(file.path(
  plotDir, paste0("CN_filt_correct_norm_compare_df_", bin_in_kbases, 
                strain, ".rds"))) 
relativeText <- paste('vs.', bamnameL[1])  # or "vs. ParentsS1S2": column used for scaling

## Make and save a separate plot for each sample 
for (samplen in sampleL) {
  p <- plotScaledCounts(scaled_df, paste(samplen, relativeText), bin_in_kbases )
  plot(p)
  ggsave(filename = 
           file.path(plotDir, 
                     paste0(samplen, 'vs', bamnameL[1], "_", bin_in_kbases, "k_CN_no_smoothing.pdf"))
         , units = "mm", width = 297, height = 210, p)
}

## Make a plot for all the samples in a single column
## Option to fix values for ymax and ymin to give better comparisons
plotset <- arrangeGrob(grobs = lapply(sampleL, function(samplen)
  {plotLogRow(scaled_df, paste(samplen, relativeText), bin_in_kbases 
              , maxCN = 2^3, minCN = 2^-2
              )}),
  ncol = 1)
plot(plotset)
ggsave(filename = 
         file.path(plotDir, 
                   paste0(strain, "_", bin_in_kbases, "k_CN.pdf"))
       , units = "mm", width = 297, height = 210, plotset)

###
## Zoom in on Chrom 4
## facet_wrap the samples in a set
require(reshape2)
A_chrom4 <- melt(Ascaled_df[Ascaled_df$chrom=="Pf3D7_04_v3", ],
                   id.vars = c("chrom", "range", "start", "end"), 
                   variable.name = "sample", value.name = "copynum")
A_chrom4$copynum[is.na(A_chrom4$copynum)] <- 1  # set undefined copynum (probably 0/0) to 1
plotA_chrom4 <- ggplot(data = A_chrom4, 
                          aes(x=end/1E5, y = copynum)) + 
  geom_point(size=0.5) + facet_wrap(~sample) + 
  labs(title = paste("Chrom 4 in", bin_in_kbases, "kb bins"), x='bin end position in 100kb') + 
  scale_y_continuous(breaks=seq(0,max(A_chrom4$copynum)+1,1))
plotA_chrom4

## Plot all the samples for a set in a single panel, coloured by sample
plotWholeCN(scaled_df, bin_in_kbases)
ggsave(filename = file.path(plotDir, paste0(strain, "_", bin_in_kbases, "k_CN_colours.pdf"))
       , units = "mm", width = 160, height = 60
)

## Plot all samples in set zoomed in on Region of Interest
# scaled_df[which(scaled_df$start==100001),] # checking that there is an x-axis-tick for each chromosome
# panelplot <- plotZoomedROI(scaled_df, chrom_2ch ="04", startkb = 400, endkb = 500)
# panelplot + labs(title=paste(bin_in_kbases, 
#                              "kb copy numbers zoomed in on Chrom 4 region of interest"))

