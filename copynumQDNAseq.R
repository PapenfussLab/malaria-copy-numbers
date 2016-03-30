## This reads the reference genome of P. falciparum in BioStrings format,
## reads three bam files, for the WT and two clones, then uses functions from 
## QDNAseq to calculate read depths and copynumbers.
## A 'mappability file' is required for correcting the counts
## Jocelyn Sietsma Penington  August 2015

## Small alterations for re-running with WTsolo in March 2016. QDNAseq has updated to 1.6.0
## Manipulation of chromosome names is due to a bug in QDNAseq which has 
## been fixed in QDNAseq version 1.6.0, (although not completely) but I'm not removing the code yet
library("BSgenome")
library("QDNAseq")
library("Biobase")

## BioStringsFormatSeedFile.txt has for original chromosome names
## First forge the genome data package:
#forgeBSgenomeDataPkg("BioStringsFormatSeedFile.txt")

## then at the unix shell command line build, check and install.
## Resulting reference genome library is in 
## /home/users/allstaff/penington.j/R/x86_64-unknown-linux-gnu-library/3.2

library("BSgenome.Pfalciparum3D7.PlasmoDB.3D7v3")
pfg <- BSgenome.Pfalciparum3D7.PlasmoDB.3D7v3
seqinfo(pfg)
str(pfg[[1]])
head(pfg[[1]])
## BSgenome and DNAstring item successfully loaded

genomesDir <- '/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/malaria/vaccine/genomes/'
outputDir <- '/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/malaria/vaccine/output/bincounts_copynums/'

## QDNAseq version 1.6.0 has a bug on calculateBlacklist(), so I cut-and-pasted code from 
## https://github.com/ccagc/QDNAseq/blob/master/R/createBinAnnotations.R 
## in March 2016 to edit and replace . Simplified!
## Code has been fixed on Github, so this code for calculateBlacklist won't be 
## needed from next version
calculateBlacklist <- function(bins, bedFiles, ...) {
  # JSP removed message code
  beds <- list()
  for (bed in bedFiles)
    beds[[bed]] <- read.table(bed, sep="\t", as.is=TRUE)
  combined <- beds[[1L]]
  if (length(beds) >= 2L)
    for (i in 2:length(beds))
      combined <- rbind(combined, beds[[i]])
  combined <- combined[, 1:3]
  colnames(combined) <- c("chromosome", "start", "end")
  combined <- combined[combined$chromosome %in% unique(bins$chromosome), ]
  combined <- combined[!is.na(combined$chromosome), ]
  combined$start <- combined$start + 1
  ## define correct sorting order of chromosomes as the order in which they
  ## are in the bins
  chromosomes <- unique(bins$chromosome)
  chromosomeOrder <- factor(combined$chromosome, levels=chromosomes,
                            ordered=TRUE)
  combined <- combined[order(chromosomeOrder, combined$start), ]
  # JSP removed code which combined overlapping regions
  overlap.counter <- function(x, joined) {
    chr <- x["chromosome"]
    start <- as.integer(x["start"])
    end <- as.integer(x["end"])
    overlaps <- joined[joined$chromosome == chr &
                         joined$start      <= end &
                         joined$end        >= start, ]
    bases <- 0
    for (i in rownames(overlaps))
      bases <- bases + min(end, overlaps[i, "end"]) -
      max(start, overlaps[i, "start"]) + 1
    bases / (end - start + 1) * 100
  }
  blacklist <- apply(bins, MARGIN=1L, FUN=overlap.counter, combined)
  return(blacklist)
}

## If run previously, read saved data
if (file.exists(paste0(outputDir, "Counts10k_WTsolo_C5_F6.rds"))) {
  counts10k <- readRDS(paste0(outputDir, "Counts10k_WTsolo_C5_F6.rds"))
  
  } else {
    setwd(genomesDir)
    pfBins10k <- createBins(pfg, 10, ignoreMitochondria=FALSE)
  ## default of ignoreMitochondria=TRUE ignores all chromosomes 
  
  ## Steps to make bigWig file for calculating mappability are in 'mappabilityInBigWig.sh'
  ## Length of kmer for calculations was 50 bp
    bwfile <- "Pfalciparum3D7_Genome50mer.bigWig"
    pfBins10k$mappability <- calculateMappability(pfBins10k, 
                                        bigWigFile=bwfile, 
                                        bigWigAverageOverBed='/usr/local/bioinf/bin/bigWigAverageOverBed',
                                        chrPrefix='')
  ## Mark nuclear chromosomes which are not all N as use=TRUE, others as use=FALSE
  ## (For this reference genome, all nuclear chromosome names start 'Pf3D7_', but mito and apico don't)
    pfBins10k$use <-  (substr(pfBins10k$chromosome,1,6)=="Pf3D7_") & (pfBins10k$bases > 0)
  ## Use a blacklist BEDfile to mask telomeres and centromeres. March 2016 made local copy of function.
    pfBins10k$blacklist <- calculateBlacklist(pfBins10k, 
                                             bedFiles='CentromereTelomereRegions.bed')
    pfBins10k <- AnnotatedDataFrame(pfBins10k,
                                 varMetadata=data.frame(labelDescription=c(
                                   'Chromosome name',
                                   'Base pair start position',
                                   'Base pair end position',
                                   'Percentage of non-N nucleotides (of full bin size)',
                                   'Percentage of C and G nucleotides (of non-N nucleotides)',
                                   'Average mappability of 50mers',
                                   'Whether the bin should be used in subsequent analysis steps',
                                   'Percentage overlap of bin with blacklist region'),
                                   row.names=colnames(pfBins10k)) )
  
  setwd("~/Pfvaccine_ln/output")
  
  ## Calculate median residuals of the LOESS fit from the control (WT) dataset
  ## I doubt the value of this step - these values are used by applyFilters() to 
  ## filter out bins where the count for WT is unusual. 
  ## It also means I read the file twice.
  wtCounts <- binReadCounts(pfBins10k, bamfiles='3D7WTsoloPhg_s.bam')
  pfBins10k$residual <- iterateResiduals(wtCounts)
  
  counts10k <- binReadCounts(pfBins10k, 
                            bamfiles=c('3D7WTsoloPhg_s.bam','Clone-C5_s.bam', 'Clone-F6_s.bam'),
                            bamnames=c('3D7wt', 'cloneC5', 'cloneF6'))
  
  ## Save a local copy
  saveRDS(counts10k, paste0(outputDir, "Counts10k_WTsolo_C5_F6.rds"))  
}


## Convert chromosome names so that 'chromosome' method does not return NA. 
## Need to change rownames of data as well.
## This was needed due to bug in QDNAseq that required chromosome names to be integer or 'X' or 'Y'
## The bug was fixed in QDNASeq 1.6.0, but I like the shorter names
newnames <- sapply(seqnames(pfg), function(x) {
  if (x=="M76611") {'Mito'} else
    if (x=="PFC10_API_IRAB") {'Apico'} else
      if (strsplit(x, '_')[[1]][1]=="Pf3D7") {strsplit(x, '_')[[1]][2]} 
   } )
counts10kchrNew <- counts10k
counts10kchrNew@featureData@data$chromosome <- newnames[counts10kchrNew@featureData@data$chromosome]
rownames(assayDataElement(counts10kchrNew, 'counts')) <- sub('Pf3D7_', '', 
                                   rownames(assayDataElement(counts10kchrNew, 'counts')))
rownames(assayDataElement(counts10kchrNew, 'counts')) <- sub('_v3', '', 
                                   rownames(assayDataElement(counts10kchrNew, 'counts')))
rownames(assayDataElement(counts10kchrNew, 'counts')) <- sub('M76611', 'Mito', 
                                   rownames(assayDataElement(countskchrNew, 'counts')))
rownames(assayDataElement(counts10kchrNew, 'counts')) <- sub('PFC10_API_IRAB', 'Apico', 
                                  rownames(assayDataElement(counts10kchrNew, 'counts')))
featureNames(featureData(counts10kchrNew)) <- rownames(assayDataElement(counts10kchrNew, 'counts'))

## Plot read counts for WT and clones (as log2 with log2(0) set to -1022):
plot(counts10kchrNew)

## Filter out bins based on residual > 4.0 standard deviations:
## Filtering with mappability=min-mappability-cutoff is an option.
## Blacklist is available internally March 2016
counts10kFiltered <- applyFilters(counts10kchrNew, residual=TRUE,
                    mappability=50, blacklist=TRUE, chromosomes=c('Apico', 'Mito'))

## 3-d contour plot of read depth vs GC and mappability
isobarPlot(counts10kchrNew)

## Estimate the correction for GC content and mappability. 
## This calculates a loess fit to the counts (as a function of GC and mappability)
## and adds it to the data structure
counts10kCorrected <- estimateCorrection(counts10kFiltered)
#counts1kCorrected <- estimateCorrection(counts1kchrNew)

noisePlot(counts10kCorrected) 

## Use the values added by estimateCorrection to correct the data:
## Copy numbers are counts divided by 'fit'  
copyNums10k <- correctBins(counts10kCorrected) 
## Scale copy numbers by median (default)
copyNums10kNormed <- normalizeBins(copyNums10k)

copyNums10kSmoothed <- smoothOutlierBins(copyNums10kNormed, logTransform = FALSE)
plot(copyNums10kSmoothed)

exportBins(copyNums10kSmoothed, 
           file=paste0(outputDir, "clones10kFiltCorrectNormSmoothBins.txt") )
## By default, exported values are log2 of copynumber values, with log2(0) set to -1022

clonesScaledByWT <- compareToReference(copyNums10kSmoothed, c(FALSE, 1, 1))
plot(clonesScaledByWT)
exportBins(clonesScaledByWT, file=paste0(outputDir, 'clones10kScaledByWT.txt') )
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

plotScaledCounts <- function(binCounts, orgCol) {
  ## Plot nuclear chromosomes in 2 rows
  library(gridExtra)  # to arrange the 14 chromosomes in 2 rows
  library(scales)     # Need the scales package for log2 transform
  row1 <- paste0('0', 1:9)
  nonNuc <- c('Mito', 'Apico')
  binCounts <- na.omit(subset(binCounts, !(binCounts$chrom %in% nonNuc)))
  
  binCounts$score <- binCounts[, orgCol]
  maxCN  <- max(binCounts$score)
  minCN <- min(binCounts[which(binCounts$score > 0), 'score'])
  p1 <- ggplot(binCounts[ which(binCounts$chrom %in% row1), ] , 
               aes(end/1E5, score)) + 
    geom_point() + 
    facet_grid(~chrom, space="free_x", scales="free_x"
    ) + 
    scale_x_continuous(breaks=seq(0,max(binCounts$end/1E5),5)) +
    labs(title = paste(orgCol, "bin score"), x='') +
    scale_y_continuous(trans=log2_trans(), limits=c(0.9*minCN, 1.1*maxCN), 
                    breaks = trans_breaks("log2", function(x) 2^x),
                    labels = trans_format("log2", math_format(2^.x)))
  
  p2 <- ggplot(na.omit(binCounts[ which( ! binCounts$chrom %in% row1), ] ), 
               aes(end/1E5, score)) + 
    geom_point() + 
    facet_grid(~chrom, space="free_x", scales="free_x"
    ) + 
    scale_x_continuous(breaks=seq(0,max(binCounts$end/1E5),5)) +
    labs(x="bin-end location in chromosome (bp x 100,000)") +
    scale_y_continuous(trans=log2_trans(), limits=c(0.9*minCN, 1.1*maxCN), 
                   breaks = trans_breaks("log2", function(x) 2^x),
                   labels = trans_format("log2", math_format(2^.x)))
  
  grid.arrange(p1,p2, ncol=1)
  
}

copynumdf <- convertQDNAtoDF(copyNums10kSmoothed)
plotScaledCounts(copynumdf, 'cloneC5')

plotScaledCounts(convertQDNAtoDF(clonesScaledByWT), 'cloneC5 vs. 3D7wt')
