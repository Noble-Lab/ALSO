######################################################################################################################
######################################################################################################################
#
# Expression Analysis from HTseq-count output
# Author: GB
# 20160510
#
######################################################################################################################
######################################################################################################################

sessionInfo()

require(RColorBrewer)
require(gplots)

require(ggplot2)
# Fonts: http://www.cookbook-r.com/Graphs/Fonts/
require(gridExtra)
require(grid)

# http://seananderson.ca/2013/10/19/reshape.html
require(reshape2)

require(TeachingDemos)
require(GGally)
require(colorspace)

# # install.packages("xlsx")
# require(xlsx)
require(openxlsx)

######################################################################################################################
######################################################################################################################
# Important functions for later

# version of vioplot with colour
sourceDir = "<scriptDirectory>"
source(file.path(sourceDir, "colVioplot.R"))

# The normalize.quantiles function in the Bolstad 'preprocessCore' R package 
# does not appear to retain the order of 'probes' of equal value.
# So wrote a custom version that does.
source(file.path(sourceDir, "quantNormGB.R"))

source(file.path(sourceDir, "corChart.R"))

# For loggin RPKMs in a list
addThenLog = function(x, add=1, base=2) {
  return(log(x+add, base))
}

# read key
readkey <- function() {
  cat ("Press [enter] to continue")
  line <- readline()
}


######################################################################################################################
######################################################################################################################
# Preliminaries

projectDir = "<projectDirectory>"
subDir = "20170620_RNAseq_InvDxz4abAndWTab"
workDir = "Patski.combo.AllData"

outputDir = file.path(projectDir, subDir, "expressionAnalysis")
if ( ! file.exists( outputDir ) ) {
  if ( ! dir.create( outputDir, showWarnings=FALSE, recursive=TRUE )  ) outputDir = getwd()
}
setwd(outputDir)
getwd()

#save.image() # To .RData
#load(file = ".RData")


######################################################################################################################
# Data

refdataDir="<refdataDirectory>"

dataDir="<dataDirectory>"
inputDir=file.pathj(dataDir, "HTSeq")

# Read in raw counts table and remove the last 5 lines
countFile = file.path(inputDir, paste(workDir, "ReadCounts.tsv.gz", sep="."))
rpkmFile = file.path(inputDir, paste(workDir, "RPKMs.tsv.gz", sep="."))
tpmFile = file.path(inputDir, paste(workDir, "TPMs.tsv.gz", sep="."))


###########################################################
###########################################################
# Load expression files

##########################################
#
# If uploading individual .out files from HTseq-count, use the code below:
#
# Create a merged dataframe from HTseq-count output files
# file_list        <- list.files()
# file_list
# countData       <- do.call("cbind", lapply(file_list, FUN=function(files){read.table(files,header=F, row.names=1, sep="\t")}))
# countData       <- head(countData, -5) 
# names(countData) <- list.files()
#
##########################################

# Read in raw counts table and remove the last 5 lines
countData <- read.table(gzfile(countFile), sep="\t", header=TRUE)
countData <- head(countData, -5) 
rownames(countData) = countData[,1]
countData = countData[,-1]
countData = data.matrix(countData)
head(countData)
dim(countData)

# # Remove genes with 0's across the board for DE/clustering
# countData.filtered = countData
# filterIdx = which(rowSums(countData.filtered) == 0)
# countData.filtered <- countData.filtered[-filterIdx,] 
# cat(sprintf("Filtered out %s genes with 0's across the board!\n", length(filterIdx)))


# Read in raw counts table and remove the last 5 lines
RPKMdata = read.table(file=rpkmFile, sep="\t", header=TRUE)
rownames(RPKMdata) = RPKMdata[,1]
RPKMdata = RPKMdata[,-1]
head(RPKMdata)
dim(RPKMdata)


# Read in raw counts table and remove the last 5 lines
TPMdata = read.table(file=tpmFile, sep="\t", header=TRUE)
rownames(TPMdata) = TPMdata[,1]
TPMdata = TPMdata[,-1]
head(TPMdata)
dim(TPMdata)


# Sample colors
cellLines = colnames(RPKMdata)
sampleCols = c("dodgerblue", "blue1", "blue1", "blue4", "blue4", 
               "orangered", "violetred", "red1", "red1", "red4", "red4",
               "orange", "green",
               "purple", "indianred", "greenyellow",
               "lightblue1", "lightblue1", "lightblue2", "lightblue4",
               "darkseagreen3", "darkseagreen3",
               "deepskyblue1", "deepskyblue1",
               "darkslategray", "darkslategray" 
               )
names(sampleCols) = cellLines
# barplot(rep(1, length(sampleCols)), col=sampleCols)
#
# relatedCols = grep("gray", colors(), value=TRUE)
# barplot(rep(1, length(relatedCols)), col=relatedCols, name=relatedCols, las=2, cex.names=.75) 


# http://stackoverflow.com/questions/30219738/is-there-a-way-to-programmatically-darken-the-color-given-rgb-values/30412528
# Lighter versions of colors
sampleColsLight <- readhex(file = textConnection(paste(rgb(t(col2rgb(sampleCols)), maxColorValue=255), collapse = "\n")),
                           class = "RGB")

#additive decrease of lightness
#transform to hue/lightness/saturation colorspace
sampleColsLightS25 = as(sampleColsLight, "HLS")
sampleColsLightS25@coords[, "L"] <- pmax(0, sampleColsLightS25@coords[, "L"] - 0.25)
#going via rgb seems to work better  
sampleColsLightS25 <- hex(as(sampleColsLightS25, "RGB"))

#multiplicative decrease of lightness
#transform to hue/lightness/saturation colorspace
sampleColsLightD2 <- as(sampleColsLight, "HLS")
sampleColsLightD4 <- as(sampleColsLight, "HLS")
sampleColsLightD8 <- as(sampleColsLight, "HLS")

sampleColsLightD2@coords[, "L"] <- sampleColsLightD2@coords[, "L"] * 0.50
sampleColsLightD4@coords[, "L"] <- sampleColsLightD4@coords[, "L"] * 0.25
sampleColsLightD8@coords[, "L"] <- sampleColsLightD8@coords[, "L"] * 0.125

#going via rgb seems to work better  
sampleColsLightD2 <- hex(as(sampleColsLightD2, "RGB"))
sampleColsLightD4 <- hex(as(sampleColsLightD4, "RGB"))
sampleColsLightD8 <- hex(as(sampleColsLightD8, "RGB"))

names(sampleColsLightD2) = names(sampleColsLightD4) =  names(sampleColsLightD8) = names(sampleCols)
names(sampleColsLightS25) = names(sampleCols)

# barplot(rep(10, 4), col=rbind(sampleCols, sampleColsLightD2, sampleColsLightD4, sampleColsLightD8) [,1])



###########################################################
###########################################################
# Chromosome-specific genes
# Must reference GTF using for genes counts to get chromosomes

###
# *** NOTE: THESE ARE REFSEQ GENES FROM refFlat:
# ***       ftp://ftp.illumina.com/README.txt
# ***       https://genome.ucsc.edu/FAQ/FAQdownloads.html#download33

# First, import the GTF-file that you have also used as input for htseq-count

# GTF file.
GTFfile = file.path(refdataDir, "iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf")

require(GenomicFeatures)
# txdb <- makeTranscriptDbFromGFF(GTFfile, format="gtf")
# # Error: 'makeTranscriptDbFromGFF' is defunct.
# # Use 'makeTxDbFromGFF' instead.
# # See help("Defunct")
# *** THIS IS VERY SLOW TO LOAD!!! ***
txdb <- makeTxDbFromGFF(GTFfile, format="gtf")

# attributes(txdb)
# str(txdb)
# seqlevelsStyle(txdb)
# seqinfo(txdb)
# seqlevels(txdb)
# genes = genes(txdb)
# transcripts = transcripts(txdb)
# cds = cds(txdb)

# # TOOOOO SLOW!
# getGeneChrom <- function(x){as.character(seqnames(x)[1])}
# # # Test
# # transcripts.list.per.gene[[1]]
# # getGeneChrom(transcripts.list.per.gene[[1]])
# genesChroms = lapply(X=transcripts.list.per.gene,FUN=getGeneChrom)

# CAN COLLAPSE TO CHROMS
mm10genes = data.frame(transcriptsBy(txdb, by="gene"))
mm10geneNames = unique(mm10genes[,c("group_name")])
length(mm10geneNames)
mm10geneChrs= unique(mm10genes[,c("group_name", "seqnames")])
nrow(mm10geneChrs)
whichDup = duplicated(mm10geneChrs$group_name)
sum((whichDup))
# 79
sum(!(whichDup))
# 24421
mm10geneChrs = mm10geneChrs[!whichDup,]
# First has no gene ID
mm10geneChrs = mm10geneChrs[-1,]  
nrow(mm10geneChrs)
# 24420
colnames(mm10geneChrs) = c("geneId", "chr")



###
# Ensure the gene ids match up: unfiltered
matchRes = match(rownames(countData), as.character(mm10geneChrs$geneId))
# all(1:nrow(countData) ==  matchRes)

# # Get differnces
# matchDiffs = which(1:nrow(countData) != matchRes)
# rownames(countData)[matchDiffs]
# as.character(mm10geneChrs$geneId)[matchRes[matchDiffs]]
# all(rownames(countData)[matchDiffs] == as.character(mm10geneChrs$geneId)[matchRes[matchDiffs]])

# Re-sort
mm10geneChrs = mm10geneChrs[matchRes, ]

# Check
matchResAfter = match(rownames(countData), as.character(mm10geneChrs$geneId))
all(1:nrow(countData) ==  matchResAfter)
# TRUE
matchResAfter = match(rownames(RPKMdata), as.character(mm10geneChrs$geneId))
all(1:nrow(RPKMdata) ==  matchResAfter)
all(1:nrow(TPMdata) ==  matchResAfter)
# TRUE





###
# Drop ChrY and ChrOther genes

chrYidx.orig = which(as.character(mm10geneChrs$chr) == "chrY")
chrXidx.orig = which(as.character(mm10geneChrs$chr) == "chrX")
chrAidx.orig = which(as.character(mm10geneChrs$chr) %in% paste("chr", 1:22, sep=""))
chrOther.orig = which(!as.character(mm10geneChrs$chr) %in% c("chrX", "chrY", paste("chr", 1:22, sep="")) )
length(chrYidx.orig)
length(chrXidx.orig)
length(chrAidx.orig)
length(chrOther.orig)
# Check
sum(length(chrYidx.orig), length(chrXidx.orig), length(chrAidx.orig), length(chrOther.orig))
countData[chrOther.orig,]
mm10geneChrs[chrOther.orig,]

# chrIndices = list(chrAidx, chrXidx, chrYidx)
# names(chrIndices) = c("autosomes", "chrX", "chrY")

# Perform dropping
countData = countData[-c(chrYidx.orig, chrOther.orig),]
RPKMdata = RPKMdata[-c(chrYidx.orig, chrOther.orig),]
TPMdata = TPMdata[-c(chrYidx.orig, chrOther.orig),]
mm10geneChrs = mm10geneChrs[-c(chrYidx.orig, chrOther.orig),]

# Check
chrYidx = which(as.character(mm10geneChrs$chr) == "chrY")
chrXidx = which(as.character(mm10geneChrs$chr) == "chrX")
chrAidx = which(as.character(mm10geneChrs$chr) %in% paste("chr", 1:22, sep=""))
chrOther = which(!as.character(mm10geneChrs$chr) %in% c("chrX", "chrY", paste("chr", 1:22, sep="")) )
length(chrYidx)
length(chrXidx)
length(chrAidx)
length(chrOther)
# Check
sum(length(chrYidx), length(chrXidx), length(chrAidx), length(chrOther))
# 24360




###
# New indices

# chrIndices = list(chrAidx, chrXidx)
# names(chrIndices) = c("autosomes", "chrX")
chrAidx = which(as.character(mm10geneChrs$chr) %in% paste("chr", 1:22, sep=""))
chrSafeAidx = which(as.character(mm10geneChrs$chr) %in% paste("chr", c(1:2,5:22), sep="")) # Remove known aneuploid chromosomes
chrXidx = which(as.character(mm10geneChrs$chr) == "chrX")
chr3idx = which(as.character(mm10geneChrs$chr) == "chr3")
chr4idx = which(as.character(mm10geneChrs$chr) == "chr4")
chr5idx = which(as.character(mm10geneChrs$chr) == "chr5")
chr12idx = which(as.character(mm10geneChrs$chr) == "chr12")
chr13idx = which(as.character(mm10geneChrs$chr) == "chr13")
chr14idx = which(as.character(mm10geneChrs$chr) == "chr14")
chrIndices = list(chrAidx, chrSafeAidx, chrXidx, 
                  chr3idx, chr4idx,
                  chr5idx, chr12idx,
                  chr13idx, chr14idx)
names(chrIndices) = c("autosomes", "diplosomes","chrX",
                      "chr3", "chr4",
                      "chr5", "chr12",
                      "chr13", "chr14")


###
# Boolean indices for all chromosomes
chrBoolindices = list()
theChroms = paste("chr", c(seq(1, 19), "X"), sep="")
for (theCurr in theChroms) {
  chrBoolindices[[theCurr]] = as.character(mm10geneChrs$chr) == theCurr
}
chrBoolindices[["autosomes"]] = as.character(mm10geneChrs$chr) %in% paste("chr", 1:22, sep="")

# ###
# # Ensure the gene ids match up: filtered
# matchRes = match(rownames(countData.filtered), as.character(mm10geneChrs$geneId))
# # all(1:nrow(countData.filtered) ==  matchRes)
#
# # # Get differnces
# # matchDiffs = which(1:nrow(countData.filtered) != matchRes)
# # rownames(countData.filtered)[matchDiffs]
# # as.character(mm10geneChrs$geneId)[matchRes[matchDiffs]]
# # all(rownames(countData.filtered)[matchDiffs] == as.character(mm10geneChrs$geneId)[matchRes[matchDiffs]])
#
# # Re-sort
# mm10geneChrs = mm10geneChrs[matchRes, ]
#
# # Check
# matchResAfter = match(rownames(countData.filtered), as.character(mm10geneChrs$geneId))
# all(1:nrow(countData.filtered) ==  matchResAfter)
# 
# chrYidx.filtered = which(as.character(mm10geneChrs$chr) == "chrY")
# chrXidx.filtered = which(as.character(mm10geneChrs$chr) == "chrX")
# chrAidx.filtered = which(as.character(mm10geneChrs$chr) %in% paste("chr", 1:22, sep=""))
# chrOther.filtered = which(!as.character(mm10geneChrs$chr) %in% c("chrX", "chrY", paste("chr", 1:22, sep="")) )
# length(chrYidx.filtered)
# length(chrXidx.filtered)
# length(chrAidx.filtered)
# length(chrOther.filtered)
# sum(length(chrYidx.filtered), length(chrXidx.filtered), length(chrAidx.filtered), length(chrOther.filtered))
# countData.filtered[chrOther.filtered,]
# mm10geneChrs[chrOther.filtered,]
# 
# chrIndices.filtered = list(chrAidx.filtered, chrXidx.filtered, chrYidx.filtered)
# names(chrIndices.filtered) = c("autosomes", "chrX", "chrY")

























######################################################################################################################
######################################################################################################################
# Allele-specific coverage data

# SNP level data not necessary here
# SNPlevelASEfile = file.path(dataDir, "snpPileups/Patski.combo.SNPpileup.annotated.tsv.gz")
# SNPlevelASEdata = read.table(gzfile(SNPlevelASEfile), header=TRUE, sep="\t")
# head(SNPlevelASEdata[,1:10])
# dim(SNPlevelASEdata)

# gene level data
geneLevelASEfile = file.path(dataDir, "snpPileups/Patski.combo.SNPcoverage.sumOverGenes.tsv.gz")
geneLevelASEdata = read.table(gzfile(geneLevelASEfile), header=TRUE, sep="\t")
head(geneLevelASEdata[,1:11])
head(geneLevelASEdata[,c(1:6, grep("snpCountOverGenes", colnames(geneLevelASEdata)))])
head(geneLevelASEdata[,c(1:6, grep("sumOverGenes", colnames(geneLevelASEdata)))])
head(geneLevelASEdata[,c(1:6, seq(7, 214, 8))])
tail(geneLevelASEdata[,c(1:6, seq(7, 214, 8))])
dim(geneLevelASEdata)

# Get largest isoforms for each gene
# Fast: No 'rbinding' to tables
geneLengths = geneLevelASEdata$end-geneLevelASEdata$start+1
geneLevelASEdata.uniquedIdx = c()
tally = 0
for (geneId in unique(geneLevelASEdata[,"geneID"])) {
  tally = tally + 1
  if (tally%%1000 == 0 ) {
    cat(tally, " genes processed\n", sep="")
  }
  geneIdx = which(as.character(geneLevelASEdata[,"geneID"]) == geneId)
  if (length(geneIdx)==1) { 
    geneLevelASEdata.uniquedIdx = c(geneLevelASEdata.uniquedIdx, geneIdx)
  } else {
    geneLevelASEdata.uniquedIdx = c(geneLevelASEdata.uniquedIdx, 
                                    geneIdx[ which.max(geneLengths[geneIdx]) ])
  }
}
length(geneLevelASEdata.uniquedIdx)

geneLevelASEdata.uniqued = geneLevelASEdata[geneLevelASEdata.uniquedIdx, ]
geneLevelASEdata[1:10,1:5]
geneLevelASEdata.uniqued[1:10,1:5]
head(geneLevelASEdata.uniqued[,1:11])


# Sort genes to match count table
matchRes = match(rownames(countData), as.character(geneLevelASEdata.uniqued$geneID))
# all(1:nrow(countData) ==  matchRes)

# # Get differnces
# matchDiffs = which(1:nrow(countData) != matchRes)
# rownames(countData)[matchDiffs]
# as.character(geneLevelASEdata.uniqued$geneID)[matchRes[matchDiffs]]
# all(rownames(countData)[matchDiffs] == as.character(geneLevelASEdata.uniqued$geneID)[matchRes[matchDiffs]])

# Re-sort 
geneLevelASEdata.uniqued = geneLevelASEdata.uniqued[matchRes, ]

# Check
matchResAfter = match(rownames(countData), as.character(geneLevelASEdata.uniqued$geneID))
all(1:nrow(countData) ==  matchResAfter)
# TRUE

countData[1:10,]
geneLevelASEdata.uniqued[1:10,1:5]
dim(geneLevelASEdata.uniqued)
tail(geneLevelASEdata.uniqued[,4])

# Check all chroms present:
unique(as.character(geneLevelASEdata$chr))
unique(as.character(geneLevelASEdata.uniqued$chr))

currChr = "chr4"
currChrIdx = which(as.character(geneLevelASEdata$chr) == currChr)
head(geneLevelASEdata[currChrIdx,])
currChrIdx = which(as.character(geneLevelASEdata.uniqued$chr) == currChr)
head(geneLevelASEdata.uniqued[currChrIdx,])
currChrIdx = which(as.character(mm10geneChrs$chr) == currChr)
head(geneLevelASEdata.uniqued[currChrIdx,])


###########################################################
###########################################################
# New gene lengths
geneLengths = geneLevelASEdata.uniqued$end-geneLevelASEdata.uniqued$start+1
names(geneLengths) = as.character(geneLevelASEdata.uniqued$geneID)





###########################################################
###########################################################
# SNP counts per gene (coverage)
# Total SNP coverage per cell line
# allelic SNP proportion
# allelic prop correction factor
# allelic RPKMs

geneLevelASEdata.uniqued.SNPcounts = geneLevelASEdata.uniqued[,10] # Same for all genes: geneLevelASEdata.uniqued[1:10,c(10,14,18)]
names(geneLevelASEdata.uniqued.SNPcounts) = geneLevelASEdata.uniqued[,4]

geneLevelASEdata.uniqued.autosomeAltRefRatio = c()

geneLevelASEdata.uniqued.altCov = c()
geneLevelASEdata.uniqued.refCov = c()
geneLevelASEdata.uniqued.totCov = c()
geneLevelASEdata.uniqued.altProp = c()
geneLevelASEdata.uniqued.refProp = c()

geneLevelASEdata.uniqued.adjAltCov = c()
geneLevelASEdata.uniqued.adjRefCov = c()
geneLevelASEdata.uniqued.adjTotCov = c()
geneLevelASEdata.uniqued.adjAltProp = c()
geneLevelASEdata.uniqued.adjRefProp = c()

altRefRatiosTable = c()

# 170417
geneLevelAScovAndRatios = list()

for (cellLine in cellLines) {
  # cellLine = "PatskiWT"
  altCol = paste(cellLine, ".altCov.sumOverGenes", sep="")
  refCol= paste(cellLine, ".refCov.sumOverGenes", sep="")

  altCov = geneLevelASEdata.uniqued[, altCol]   
  refCov = geneLevelASEdata.uniqued[, refCol]
  totCov = altCov + refCov
  altProp = altCov / totCov
  refProp = 1-altProp
  
  geneLevelASEdata.uniqued.altCov = cbind( geneLevelASEdata.uniqued.altCov, altCov )
  geneLevelASEdata.uniqued.refCov = cbind( geneLevelASEdata.uniqued.refCov, refCov )
  geneLevelASEdata.uniqued.totCov = cbind( geneLevelASEdata.uniqued.totCov, totCov )
  geneLevelASEdata.uniqued.altProp = cbind( geneLevelASEdata.uniqued.altProp, altProp )
  geneLevelASEdata.uniqued.refProp = cbind( geneLevelASEdata.uniqued.refProp, 1-altProp)
  
  ######
  # Adjust Reference counts
  
  # Expect this to be 1:1 for autosomes 
  # But not always the case if there are aneuploidies
  obsAutosomeAltRefRatioGenome = sum(altCov, na.rm=TRUE) / sum(refCov, na.rm=TRUE)
  obsAutosomeAltRefRatio = sum(altCov[ chrIndices[["autosomes"]] ], na.rm=TRUE) / sum(refCov[ chrIndices[["autosomes"]] ], na.rm=TRUE)
  obsAutosomeAltRefRatioSafe = sum(altCov[ chrIndices[["diplosomes"]] ], na.rm=TRUE) / sum(refCov[ chrIndices[["diplosomes"]] ], na.rm=TRUE)
  
  altRefRatiosTable = rbind(altRefRatiosTable, c(obsAutosomeAltRefRatioGenome, obsAutosomeAltRefRatio, obsAutosomeAltRefRatioSafe))
  
  # Correction of refcov
  if (obsAutosomeAltRefRatio < 1) {
    adjRefCov = refCov * obsAutosomeAltRefRatio
    adjAltCov = altCov
  } else {
    adjAltCov = altCov * (1/obsAutosomeAltRefRatio)
    adjRefCov = refCov 
  }
  adjTotCov = adjAltCov + adjRefCov
  adjAltProp = adjAltCov / adjTotCov
  adjRefProp = 1-adjAltProp
  
  geneLevelASEdata.uniqued.adjAltCov = cbind( geneLevelASEdata.uniqued.adjAltCov, adjAltCov )
  geneLevelASEdata.uniqued.adjRefCov = cbind( geneLevelASEdata.uniqued.adjRefCov, adjRefCov )
  geneLevelASEdata.uniqued.adjTotCov = cbind( geneLevelASEdata.uniqued.adjTotCov, adjTotCov )
  geneLevelASEdata.uniqued.adjAltProp = cbind( geneLevelASEdata.uniqued.adjAltProp, adjAltProp )
  geneLevelASEdata.uniqued.adjRefProp = cbind( geneLevelASEdata.uniqued.adjRefProp, 1-adjAltProp)
  geneLevelASEdata.uniqued.autosomeAltRefRatio = c(geneLevelASEdata.uniqued.autosomeAltRefRatio, obsAutosomeAltRefRatio)
  
  # 170417
  geneLevelAScovAndRatios[[cellLine]] = cbind(geneLevelASEdata.uniqued[,1:3],
                                              geneLevelASEdata.uniqued.SNPcounts, 
                                              altCov, refCov, totCov, altProp, refProp,
                                              adjAltCov, adjRefCov, adjTotCov, adjAltProp, adjRefProp)
  colnames(geneLevelAScovAndRatios[[cellLine]]) = c("chr", "start", "end", 
                                                   "SNPcounts", 
                                                   "altCov", "refCov", "totCov", "altProp", "refProp",
                                                   "adjAltCov", "adjRefCov", "adjTotCov", "adjAltProp", "adjRefProp")
  rownames(geneLevelAScovAndRatios[[cellLine]]) = rownames(geneLevelAScovAndRatios[[cellLine]])
  
}

names(geneLevelASEdata.uniqued.autosomeAltRefRatio) = cellLines

colnames(geneLevelASEdata.uniqued.altCov) = colnames(geneLevelASEdata.uniqued.refCov) = 
  colnames(geneLevelASEdata.uniqued.totCov) = colnames(geneLevelASEdata.uniqued.altProp) = colnames(geneLevelASEdata.uniqued.refProp) = 
  colnames(geneLevelASEdata.uniqued.adjAltCov) = colnames(geneLevelASEdata.uniqued.adjRefCov) = 
  colnames(geneLevelASEdata.uniqued.adjTotCov) = colnames(geneLevelASEdata.uniqued.adjAltProp) = colnames(geneLevelASEdata.uniqued.adjRefProp) = cellLines

rownames(geneLevelASEdata.uniqued.altCov) = rownames(geneLevelASEdata.uniqued.refCov) = 
  rownames(geneLevelASEdata.uniqued.totCov) = rownames(geneLevelASEdata.uniqued.altProp) = rownames(geneLevelASEdata.uniqued.refProp) = 
  rownames(geneLevelASEdata.uniqued.adjAltCov) = rownames(geneLevelASEdata.uniqued.adjRefCov) = 
  rownames(geneLevelASEdata.uniqued.adjTotCov) = rownames(geneLevelASEdata.uniqued.adjAltProp) = rownames(geneLevelASEdata.uniqued.adjRefProp) = 
  as.character(geneLevelASEdata.uniqued$geneID)

rownames(altRefRatiosTable) = cellLines
colnames(altRefRatiosTable) = c("genome", "autosomes", "diplosomes")

altRefRatiosTableOut = apply(altRefRatiosTable, 2, function(x) round(as.numeric(as.character(x)), 3))
rownames(altRefRatiosTableOut) = cellLines


altRefRatioDir = file.path(outputDir, "altRefRatios")
if ( ! file.exists( altRefRatioDir ) ) {
  if ( ! dir.create( altRefRatioDir, showWarnings=FALSE, recursive=TRUE )  ) altRefRatioDir = outputDir
}

altRefRatiosTableOutFilename = file.path(altRefRatioDir, paste("altRefRatios.tsv", sep = ""))

write.table(altRefRatiosTableOut, altRefRatiosTableOutFilename, 
            append = FALSE, sep = "\t", 
            col.names = NA, row.names = TRUE, quote = FALSE)

# try(system(paste(binDir, "/rdb2html -noformatline  \"", sub("~/", "$HOME/", altRefRatiosTableOutFilename), "\" > \"", 
#                  sub(".tsv", ".html", sub("~/", "$HOME/", altRefRatiosTableOutFilename)), "\"", sep = "")))

# Check
altRefRatiosTableOut
geneLevelASEdata.uniqued.autosomeAltRefRatio


######  
# Determine allelic RPKMs based on adjusted allelic proportion

expressionTables = list()
expressionTables[["RPKM"]] = RPKMdata
expressionTables[["altRPKM"]] = expressionTables[["RPKM"]] * geneLevelASEdata.uniqued.adjAltProp 
expressionTables[["refRPKM"]] = expressionTables[["RPKM"]] * geneLevelASEdata.uniqued.adjRefProp

# Check
colSums(expressionTables[["RPKM"]], na.rm=TRUE)

colSums(expressionTables[["altRPKM"]], na.rm=TRUE)

colSums(expressionTables[["refRPKM"]], na.rm=TRUE)

all(expressionTables[["altRPKM"]] + expressionTables[["refRPKM"]]  == expressionTables[["RPKM"]])
# FALSE
# *** NOT ALL GENES WILL BE ASSOCIATED WITH SNPS ***

colSums(expressionTables[["altRPKM"]], na.rm=TRUE) / colSums(expressionTables[["refRPKM"]], na.rm=TRUE)
# NOTE: altRPKM / refRPKM != genome-wide alt:ref SNP coverage, 
#                       but instead close to 1 B/C adjusted

colSums(geneLevelASEdata.uniqued.altCov, na.rm=TRUE) / colSums(geneLevelASEdata.uniqued.refCov, na.rm=TRUE)


######  
# Determine allelic TPMs based on adjusted allelic proportion

expressionTables[["TPM"]] = TPMdata 
expressionTables[["altTPM"]] = expressionTables[["TPM"]] * geneLevelASEdata.uniqued.adjAltProp
expressionTables[["refTPM"]] = expressionTables[["TPM"]] * geneLevelASEdata.uniqued.adjRefProp

# Check
colSums(expressionTables[["TPM"]], na.rm=TRUE)

colSums(expressionTables[["altTPM"]], na.rm=TRUE)

colSums(expressionTables[["refTPM"]], na.rm=TRUE)

colSums(expressionTables[["altTPM"]] + expressionTables[["refTPM"]], na.rm=TRUE)
# NOTE: alt + ref != TPM (niether for RPKM)
#     This is b/c some genes have no SNPs (or SNP coverage), so have NAs for proportion.

all(expressionTables[["altTPM"]] + expressionTables[["refTPM"]]  == expressionTables[["TPM"]])
# FALSE
# *** NOT ALL GENES WILL BE ASSOCIATED WITH SNPS ***

colSums(expressionTables[["altTPM"]], na.rm=TRUE) / colSums(expressionTables[["refTPM"]], na.rm=TRUE)
# NOTE: altTPM / refTPM != genome-wide alt:ref SNP coverage, 
#                       but instead close to 1 B/C adjusted

colSums(geneLevelASEdata.uniqued.adjAltCov, na.rm=TRUE) / colSums(geneLevelASEdata.uniqued.adjRefCov, na.rm=TRUE)


###
# ChromosomeByChrom normalized values (to get around effects of aneuploidy)
expressionTables.chromNorm = list()
expressionTables.chromNorm[["PCRPKM"]] = PCRPKMtable
expressionTables.chromNorm[["altPCRPKM"]] = expressionTables.chromNorm[["PCRPKM"]] * geneLevelASEdata.uniqued.adjAltProp
expressionTables.chromNorm[["refPCRPKM"]] = expressionTables.chromNorm[["PCRPKM"]] * geneLevelASEdata.uniqued.adjRefProp

expressionTables.chromNorm[["PCTPM"]] = PCTPMtable
expressionTables.chromNorm[["altPCTPM"]] = expressionTables.chromNorm[["PCTPM"]] * geneLevelASEdata.uniqued.adjAltProp
expressionTables.chromNorm[["refPCTPM"]] = expressionTables.chromNorm[["PCTPM"]] * geneLevelASEdata.uniqued.adjRefProp



#save.image() # To .RData
#load(file = ".RData")






###########################################################
###########################################################
# Save expression tables with expression for all samples

tableDir = file.path(outputDir, "expressionTables")
if ( ! file.exists( tableDir ) ) {
  if ( ! dir.create( tableDir, showWarnings=FALSE, recursive=TRUE )  ) tableDir = outputDir
}

datasetSelections = colnames(countData)

datasetRenames = colnames(countData)

names(datasetRenames) = datasetSelections

for (expDataType in names(expressionTables)) {
  
  bigOldTable = cbind( as.character(geneLevelASEdata.uniqued[,4]),
                       as.character(geneLevelASEdata.uniqued[,1]),
                       format(geneLevelASEdata.uniqued[,2], scientific=FALSE, trim=TRUE, justify="none"),
                       format(geneLevelASEdata.uniqued[,3], scientific=FALSE, trim=TRUE, justify="none"),      
                       expressionTables[[expDataType]][,datasetSelections]
                      )
  colnames(bigOldTable) = c("geneId", "chr", "start", "stop", datasetRenames)
    
  bigOldFilename = file.path(tableDir, paste(expDataType, ".tsv", sep = ""))

  write.table(bigOldTable, bigOldFilename, 
              append = FALSE, sep = "\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  try(system(paste("gzip -f \"", sub("~/", "$HOME/", bigOldFilename), "\"", sep = "")))
}


##########
# 20171205
# Save with strand info
tableDir = file.path(outputDir, "expressionTables.withStrand")
if ( ! file.exists( tableDir ) ) {
  if ( ! dir.create( tableDir, showWarnings=FALSE, recursive=TRUE )  ) tableDir = outputDir
}

datasetSelections = colnames(countData)

datasetRenames = colnames(countData)

names(datasetRenames) = datasetSelections

for (expDataType in names(expressionTables)) {
  
  bigOldTable = cbind( as.character(geneLevelASEdata.uniqued[,4]),
                       as.character(geneLevelASEdata.uniqued[,1]),
                       format(geneLevelASEdata.uniqued[,2], scientific=FALSE, trim=TRUE, justify="none"),
                       format(geneLevelASEdata.uniqued[,3], scientific=FALSE, trim=TRUE, justify="none"),
                       as.character(geneLevelASEdata.uniqued[,6]),
                       expressionTables[[expDataType]][,datasetSelections]
  )
  colnames(bigOldTable) = c("geneId", "chr", "start", "stop", "strand", datasetRenames)
  
  bigOldFilename = file.path(tableDir, paste(expDataType, ".tsv", sep = ""))
  
  write.table(bigOldTable, bigOldFilename, 
              append = FALSE, sep = "\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  try(system(paste("gzip -f \"", sub("~/", "$HOME/", bigOldFilename), "\"", sep = "")))
}


####################
# chrX-specific data
# See '20170620_RNAseq_InvDxz4abAndWTab_expressionAnalysis_HTSeqAndExpressionAnalysis_workflow.sh'
#
# cd ~/proj/2017_RNAseq_analysis/results/gbonora/20170620_RNAseq_InvDxz4abAndWTab/expressionAnalysis/expressionTables
# cd ~/proj/2017_RNAseq_analysis/results/gbonora/20170620_RNAseq_InvDxz4abAndWTab/expressionAnalysis/expressionTables.withStrand
# 
# # while read FILE; do echo "${FILE/.tsv.gz/.chrX.tsv.gz}"; done < <( ls -1 *.tsv.gz )
# while read FILE; do 
# # FILE="counts.tsv.gz"
# zcat "${FILE}" | head -n 1 > "${FILE/.tsv.gz/.chrX.tsv}"; 
# zcat "${FILE}" | grep "chrX" >> "${FILE/.tsv.gz/.chrX.tsv}"; 
# gzip -f "${FILE/.tsv.gz/.chrX.tsv}"; 
# done < <( ls -1 *.tsv.gz | grep -v "chrX" )









###########################################################
###########################################################
# Save expression bedgraphs

bedDir = file.path(outputDir, "bedgraphs")
if ( ! file.exists( bedDir ) ) {
  if ( ! dir.create( bedDir, showWarnings=FALSE, recursive=TRUE )  ) bedDir = outputDir
}

datasetSelections = colnames(countData)

datasetRenames = colnames(countData)

names(datasetRenames) = datasetSelections

for (expDataType in names(expressionTables)) {

  bedDir = file.path(outputDir, "bedgraphs", expDataType)
  if ( ! file.exists( bedDir ) ) {
    if ( ! dir.create( bedDir, showWarnings=FALSE, recursive=TRUE )  ) bedDir = outputDir
  }
  
  for( subSetName in datasetSelections) {
      
      cat(expDataType, " ", subSetName, "\n", sep="")
    
    #       expValues = expressionTables[[expDataType]][,subSetName]
    #       expValues[is.na(expValues)] = "NaN"
    
      bedTable = cbind( as.character(geneLevelASEdata.uniqued[,1]),
                             format(geneLevelASEdata.uniqued[,2], scientific=FALSE, trim=TRUE, justify="none"),
                             format(geneLevelASEdata.uniqued[,3], scientific=FALSE, trim=TRUE, justify="none"),      
                             expressionTables[[expDataType]][,subSetName],
                             as.character(geneLevelASEdata.uniqued[,4]) )
                             
      bedFilename = file.path(bedDir, paste(expDataType, "_", datasetRenames[subSetName], ".bedgraph", sep = ""))
      bed2Filename = file.path(bedDir, paste(expDataType, "_", datasetRenames[subSetName], ".noGeneId.bedgraph", sep = ""))
      
      cat("track type=bedGraph ",
          "name=\"", datasetRenames[subSetName], "\"", 
          "description=\"", datasetRenames[subSetName], " ", expDataType, "\"", 
          "visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0:1000 autoScale=off\n", 
          file=bedFilename, sep="")
      cat("track type=bedGraph ",
          "name=\"", datasetRenames[subSetName], "\"", 
          "description=\"", datasetRenames[subSetName], " ", expDataType, "\"", 
          "visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0:1000 autoScale=off\n", 
          file=bed2Filename, sep="")
      
      write.table(bedTable, bedFilename, 
                  append = TRUE, sep = "\t", 
                  col.names = FALSE, row.names = FALSE, quote = FALSE)
      try(system(paste("gzip -f \"", sub("~/", "$HOME/", bedFilename), "\"", sep = "")))
      
      write.table(bedTable[,c(1:4)], bed2Filename, 
                  append = TRUE, sep = "\t", 
                  col.names = FALSE, row.names = FALSE, quote = FALSE)
      try(system(paste("gzip -f \"", sub("~/", "$HOME/", bed2Filename), "\"", sep = "")))
      
 
  }
}











###########################################################
###########################################################
# Find expressed genes

######################
# One line over 1 RPKM
# *** USE QNORMED RPKMS? NOPE ***

minRPKM = 1
expressedGenesBoolIdx = apply(RPKMdata, 1, max, na.rm=TRUE) >= minRPKM
sum(expressedGenesBoolIdx)/length(expressedGenesBoolIdx) * 100 
# # 56.4% >= 1 RPKM
# 56.62151

sum(apply(RPKMdata, 1, max, na.rm=TRUE) > 0)/nrow(RPKMdata)
# 86% > 0 RPKM
# 86.42%

head(expressedGenesBoolIdx)
tail(RPKMdata)

tail(expressedGenesBoolIdx)
tail(RPKMdata)



######################
# One line over 1 TPM

minTPM = 1
TPMexpressedGenesBoolIdx = apply(TPMdata, 1, max, na.rm=TRUE) >= minTPM
sum(TPMexpressedGenesBoolIdx)/length(TPMexpressedGenesBoolIdx) * 100 
#  63% >= 1 TPM

sum(apply(TPMdata, 1, max, na.rm=TRUE) > 0)/nrow(TPMdata)
# 86% > 0 TPM as with RPKM










###########################################################
###########################################################
# SNP count and coverage

SNPcovDir = file.path(outputDir, "SNPcoverageData")
if ( ! file.exists( SNPcovDir ) ) {
  if ( ! dir.create( SNPcovDir, showWarnings=FALSE, recursive=TRUE )  ) SNPcovDir = outputDir
}

###
# SNP count distributions

PDFfilename = file.path(SNPcovDir, 
                        paste(workDir, "_hist_SNPcountsPerGene.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
# par(mfrow=c(2,2), oma=c(0,0,2.5,0))
layout(matrix(c(1,2,3,4),2,2,TRUE))
par(oma=c(0,0,2.5,0))

for (expressed in c(FALSE, TRUE)) {
  #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
  for (whichChrs in c("autosomes", "chrX") ) {
    
    if (expressed) {
      if (whichChrs == "genome-wide") {
        SNPcountVector = geneLevelASEdata.uniqued.SNPcounts[expressedGenesBoolIdx]
      } else {
        SNPcountVector = geneLevelASEdata.uniqued.SNPcounts[intersect(which(expressedGenesBoolIdx),
                                                                     chrIndices[[whichChrs]])]
      }
    } else {
      if (whichChrs == "genome-wide") {
        SNPcountVector = geneLevelASEdata.uniqued.SNPcounts
      } else {
        SNPcountVector = geneLevelASEdata.uniqued.SNPcounts[chrIndices[[whichChrs]]]
      }
    }
    
    ret = hist(SNPcountVector,
               freq=TRUE, xlim =c(0, 180), 
               breaks=c(-0.5, 0.5, 20.5, 40.5, 60.5, 80.5, 100.5, 
                        120.5, 140.5, 160.5, 
                        max(SNPcountVector, na.rm=TRUE)+0.5),
               main=paste(ifelse(expressed, "expressed", "all"), " genes\n", 
                          whichChrs, " (n=", length(SNPcountVector), ")", sep=""),
               xlab="SNPs per gene")
    
    pieValues = c(ret$counts[1], sum(ret$counts[-1]) )
    cols <- c("grey90", "grey50")
    percentlabels<- round(100*pieValues/sum(pieValues), 1)
    pielabels<- paste(percentlabels, "%", sep="")
    #     subPlotRef = subplot(pie(pieValues, main="Genes with SNPs", col=cols, labels=pielabels, cex=0.8), 
    #                          "top", size=c(3,3), inset=c(0, 0.2))
    #     op <- par(no.readonly=TRUE)
    #     par(subPlotRef)
    #     legend("topright", c("0 SNPs", "1+ SNPs"), 
    #            cex=0.8, fill=cols, inset=c(-0.1,0), xpd=TRUE)
    #     par(op)
    subplot( pie(pieValues, main="Genes with SNPs", col=cols, labels=pielabels, cex=0.8), 
             "top", size=c(2,2), inset=c(0, 0.2))
    legend("topright", c("0 SNPs", "1+ SNPs"), cex=0.8, fill=cols, inset=c(0.15, 0.15))
    

  }
}
title(main=paste("Distribution of SNP counts within genes", sep=""), outer=TRUE)
dev.off()




###
# SNP coverage distributions
for (expressed in c(FALSE, TRUE)) {
  for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
    # cellLines = c("PatskiDel1", "PatskiDel2", "PatskiDel5", "PatskiWT")
    for ( cellLine in cellLines) {
                   
      if (expressed) {
        if (whichChrs == "genome-wide") {
          SNPcovVector = geneLevelASEdata.uniqued.totCov[expressedGenesBoolIdx, cellLine]
        } else {
          SNPcovVector = geneLevelASEdata.uniqued.totCov[intersect(which(expressedGenesBoolIdx),
                                                                    chrIndices[[whichChrs]]), cellLine]
        }
      } else {
        if (whichChrs == "genome-wide") {
          SNPcovVector = geneLevelASEdata.uniqued.totCov[, cellLine]
        } else {
          SNPcovVector = geneLevelASEdata.uniqued.totCov[chrIndices[[whichChrs]], cellLine]
        }
      }
      
      # plot distribution
      PDFfilename = file.path(SNPcovDir, 
                              paste(workDir, "_hist_SNPtotalReadCoveragePerGene_", 
                                    whichChrs, 
                                    ifelse(expressed, "_expressedGenes_", "_allGenes_"), 
                                    cellLine,
                                    ".pdf", sep=""))
      pdf(PDFfilename, height=6, width=6)
      ret = hist(SNPcovVector, col=sampleCols[cellLine],
                 freq=TRUE, xlim =c(0, 110), 
                 breaks=c(-0.5, 9.5, 19.5, 29.5, 39.5, 49.5, 
                          59.5, 69.5, 79.5, 89.5, 99.5,  
                          max(SNPcovVector, na.rm=TRUE)+0.5),
                 main=paste(cellLine, ": distribution of total SNP read coverage for ", 
                            ifelse(expressed, "expressed", "all"), " genes\n", 
                            whichChrs, " (n=", length(SNPcovVector), ")", sep=""),
                 xlab="total SNP read coverage")
      
      pieValues = c(ret$counts[1], sum(ret$counts[-1]) )
      cols <- c("grey90", "grey50")
      percentlabels<- round(100*pieValues/sum(pieValues), 1)
      pielabels<- paste(percentlabels, "%", sep="")
      subPlotRef = subplot(pie(pieValues, main="Total SNP read coverage per gene", col=cols, labels=pielabels, cex=0.8), 
                           "top", size=c(3,3), inset=c(0, 0.2))
      op <- par(no.readonly=TRUE)
      par(subPlotRef)
      legend("topright", c("<10X", ">=10X"), 
             cex=0.8, fill=cols, inset=c(-0.1,0))
      par(op)
      
      dev.off()
    }
  }     
}


###########################################################
###########################################################
# Allelic Proportions

ASEdir = file.path(outputDir, "allelicProportions")
if ( ! file.exists( ASEdir ) ) {
  if ( ! dir.create( ASEdir, showWarnings=FALSE, recursive=TRUE )  ) ASEdir = outputDir
}

###
# Allelic proportion distributions
for (expressed in c(FALSE, TRUE)) {
  for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
    # cellLines = c("PatskiDel1", "PatskiDel2", "PatskiDel5", "PatskiWT")
    for ( cellLine in cellLines) {
      
      if (expressed) {
        if (whichChrs == "genome-wide") {
          altProp = geneLevelASEdata.uniqued.altProp[expressedGenesBoolIdx, cellLine]
        } else {
          altProp = geneLevelASEdata.uniqued.altProp[intersect(which(expressedGenesBoolIdx),
                                                                     chrIndices[[whichChrs]]), cellLine]
        }
      } else {
        if (whichChrs == "genome-wide") {
          altProp = geneLevelASEdata.uniqued.altProp[, cellLine]
        } else {
          altProp = geneLevelASEdata.uniqued.altProp[chrIndices[[whichChrs]], cellLine]
        }
      }
      
      # plot distribution
      PDFfilename = file.path(ASEdir, 
                              paste(workDir, "_hist_AllelicProportionPerGene_", 
                                    whichChrs, 
                                    ifelse(expressed, "_expressedGenes_", "_allGenes_"), 
                                    cellLine,
                                    ".pdf", sep=""))
      pdf(PDFfilename, height=6, width=6)
      ret = hist(altProp, col=sampleCols[cellLine],
                 freq=TRUE, xlim =c(0, 1.05), 
                 breaks=seq(-0.05, 1.05, 0.1),
                 main=paste(cellLine, ": distribution of allelic proportion for ", 
                            ifelse(expressed, "expressed", "all"), " genes\n", 
                            whichChrs, " (n=", length(altProp), ")", sep=""),
                 xlab="Allelic proportion (Xa/[Xa+Xi])")

      dev.off()  
    }
  }     
}


###
# Adjusted allelic proportion distributions
for (expressed in c(FALSE, TRUE)) {
  for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
    # cellLines = c("PatskiDel1", "PatskiDel2", "PatskiDel5", "PatskiWT")
    for ( cellLine in cellLines) {
      
      if (expressed) {
        if (whichChrs == "genome-wide") {
          altProp = geneLevelASEdata.uniqued.adjAltProp[expressedGenesBoolIdx, cellLine]
        } else {
          altProp = geneLevelASEdata.uniqued.adjAltProp[intersect(which(expressedGenesBoolIdx),
                                                                         chrIndices[[whichChrs]]), cellLine]
        }
      } else {
        if (whichChrs == "genome-wide") {
          altProp = geneLevelASEdata.uniqued.adjAltProp[, cellLine]
        } else {
          altProp = geneLevelASEdata.uniqued.adjAltProp[chrIndices[[whichChrs]], cellLine]
        }
      }
      
      # plot distribution
      PDFfilename = file.path(ASEdir, 
                              paste(workDir, "_hist_AdjustedAllelicProportionPerGene_", 
                                    whichChrs, 
                                    ifelse(expressed, "_expressedGenes_", "_allGenes_"), 
                                    cellLine,
                                    ".pdf", sep=""))
      pdf(PDFfilename, height=6, width=6)
      ret = hist(altProp, col=sampleCols[cellLine],
                 freq=TRUE, xlim =c(0, 1.05), 
                 breaks=seq(-0.05, 1.05, 0.1),
                 main=paste(cellLine, ": distribution of adjusted allelic proportion for ", 
                            ifelse(expressed, "expressed", "all"), " genes\n", 
                            whichChrs, " (n=", length(altProp), ")", sep=""),
                 xlab="Allelic proportion (Xa/[Xa+Xi])")
      legend("topleft", paste("autosomes\nalt/ref ratio: ", round(geneLevelASEdata.uniqued.autosomeAltRefRatio[cellLine], 2), sep=""), cex=0.8)
      
      dev.off()  
    }
  }     
}


######################################################################################################################
######################################################################################################################
# Allelic Ratios

ASEdir = file.path(outputDir, "allelicRatios")
if ( ! file.exists( ASEdir ) ) {
  if ( ! dir.create( ASEdir, showWarnings=FALSE, recursive=TRUE )  ) ASEdir = outputDir
}

minCov = 5

minCovIndices = list()
for (cellLine in cellLines) {
  minCovIndices[[cellLine]] = which(geneLevelAScovAndRatios[[cellLine]]$totCov >= minCov) 
}


###
# Allelic proportion distributions - Adjusted vs unadjusted - overlaid density plots

# vector of spretus minima
autosomalMinima = list()
autosomalProps = list()
chrXSubProps = list()

for (cellLine in cellLines) {
  # cellLine=cellLines[1]
  print(cellLine)
  
  allelicPropVector=list()
  ggplotList = list()
 
  autosomalMinima[[cellLine]] = list()
  autosomalProps[[cellLine]] = list()
  chrXSubProps[[cellLine]] = list()
  
  # Density   
  #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
  for (whichChrs in c("autosomes", "diplosomes", "chrX") ) {
    # whichChrs="autosomes"
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[["unadjusted"]] = as.data.frame(geneLevelAScovAndRatios[[cellLine]]$altProp[ minCovIndices[[cellLine]] ])
      allelicPropVector[["adjusted"]] = as.data.frame(geneLevelAScovAndRatios[[cellLine]]$adjAltProp[ minCovIndices[[cellLine]] ])
    } else {
      allelicPropVector[["unadjusted"]] = as.data.frame(geneLevelAScovAndRatios[[cellLine]]$altProp[ intersect(chrIndices[[whichChrs]], minCovIndices[[cellLine]]) ])
      allelicPropVector[["adjusted"]] = as.data.frame(geneLevelAScovAndRatios[[cellLine]]$adjAltProp[ intersect(chrIndices[[whichChrs]], minCovIndices[[cellLine]]) ])
    }
    names(allelicPropVector[["unadjusted"]]) = "ratio"
    names(allelicPropVector[["adjusted"]]) = "ratio"
    
    list.names <- names(allelicPropVector)
    lns <- sapply(allelicPropVector, nrow)
    dat <- as.data.frame(do.call("rbind", allelicPropVector))
    colnames(dat)[1] = "ratio"
    dat$group <- rep(list.names, lns)
    
    theCols2use = c(sampleCols[cellLine], sampleColsLightD8[cellLine])
    names(theCols2use) = list.names
    
    ggplotList[[whichChrs]] <- ggplot(dat, aes(ratio, fill=group)) + 
      geom_density(alpha=0.6, size=0.25) +
      ggtitle(paste(whichChrs, "\n(n=", lns[1], ")", sep="")) +
      scale_fill_manual(values=theCols2use) +
      labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
      theme(legend.position="top") +
      theme(legend.title=element_blank())

    # Get B6 and Sp-specifc autosomal minima & proportions below and above, resp.
    if ( whichChrs == "autosomes" || whichChrs == "diplosomes"  ) {

      # autosomal minima
      ggplotData = print(ggplotList[[whichChrs]])
      # names(ggplotData)
      # str(ggplotData[["data"]])
      # attributes(ggplotData[["data"]][[1]])
      yVals = ggplotData[["data"]][[1]]$y
      xVals = ggplotData[["data"]][[1]]$x
      idx1 = which.min(yVals[xVals<0.5])
      yMin1 = yVals[xVals<0.5][idx1]
      xMin1 = xVals[xVals<0.5][idx1]
      idx2 = which.min(yVals[xVals>0.5])
      yMin2 = yVals[xVals>0.5][idx2]
      xMin2 = xVals[xVals>0.5][idx2]
      
      autosomalMinima[[cellLine]][[whichChrs]] = list()
      autosomalMinima[[cellLine]][[whichChrs]][["xMin1"]] = xMin1
      autosomalMinima[[cellLine]][[whichChrs]][["yMin1"]] = yMin1
      autosomalMinima[[cellLine]][[whichChrs]][["xMin2"]] = xMin2
      autosomalMinima[[cellLine]][[whichChrs]][["yMin1"]] = yMin2
      
      # B6 and Sp-specific proportions 
      B6autoProp = sum(ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$x<xMin1]) / sum(ggplotData[["data"]][[1]]$y)
      SpAutoProp = sum(ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$x>xMin2]) / sum(ggplotData[["data"]][[1]]$y)
      
      autosomalProps[[cellLine]][[whichChrs]] = list()
      autosomalProps[[cellLine]][[whichChrs]][["B6autoProp"]] = B6autoProp
      autosomalProps[[cellLine]][[whichChrs]][["SpAutoProp"]] = SpAutoProp
      
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
      
      ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
        geom_vline(xintercept=xMin1, linetype="dashed") +
        geom_vline(xintercept=xMin2, linetype="dashed") +
        annotate(geom="text", label=paste("(", round(xMin1, 2), "; ", round(yMin1, 2), "); Pr(Y|X<", round(xMin1, 2), ") = ", round(B6autoProp, 2), sep=""), x=xMin1, y=yMid,  angle = 90, vjust=-1) +
        annotate(geom="text", label=paste("(", round(xMin2, 2), "; ", round(yMin2, 2), "); Pr(Y|X>", round(xMin2, 2), ") = ", round(SpAutoProp, 2), sep=""), x=xMin2, y=yMid,  angle = 90, vjust=+1.5)
    }    

    if ( whichChrs == "chrX" ) {
      xMin1 = autosomalMinima[[cellLine]][["diplosomes"]][["xMin1"]]
      yMin1 = autosomalMinima[[cellLine]][["diplosomes"]][["yMin1"]] 
      xMin2 = autosomalMinima[[cellLine]][["diplosomes"]][["xMin2"]] 
      yMin2 = autosomalMinima[[cellLine]][["diplosomes"]][["yMin1"]] 
      
      ggplotData = print(ggplotList[[whichChrs]])
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65

      # B6 and Sp-sub proportions 
      B6subProp = sum(ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$x<xMin1]) / sum(ggplotData[["data"]][[1]]$y)
      SpSubProp = sum(ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$x<xMin2]) / sum(ggplotData[["data"]][[1]]$y)
      
      chrXSubProps[[cellLine]][["B6subProp"]] = B6subProp
      chrXSubProps[[cellLine]][["SpSubProp"]] = SpSubProp
      
      ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
        geom_vline(xintercept=xMin1, linetype="dashed") +
        geom_vline(xintercept=xMin2, linetype="dashed") +
        annotate(geom="text", label=paste("(", round(xMin1, 2), "; ", round(yMin1, 2), "); Pr(Y|X<", round(xMin1, 2), ") = ", round(B6subProp, 2), sep=""), x=xMin1, y=yMid,  angle = 90, vjust=-1) +
        annotate(geom="text", label=paste("(", round(xMin2, 2), "; ", round(yMin2, 2), "); Pr(Y|X<", round(xMin2, 2), ") = ", round(SpSubProp, 2), sep=""), x=xMin2, y=yMid,  angle = 90, vjust=-1)

      ggplotList[["chrXzoom"]] = ggplot(dat, aes(ratio, fill=group)) + 
        geom_density(alpha=0.6, size=0.25) +
        ggtitle(paste("chrX zoomed\n(n=", lns[1], ")", sep="")) +
        scale_fill_manual(values=theCols2use) +
        labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
        theme(legend.position="top") +
        theme(legend.title=element_blank()) + 
        coord_cartesian(ylim = c(0, 5))

      ggplotData = print(ggplotList[["chrXzoom"]])
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
      
      ggplotList[["chrXzoom"]] =  ggplotList[["chrXzoom"]] + 
        geom_vline(xintercept=xMin1, linetype="dashed") +
        geom_vline(xintercept=xMin2, linetype="dashed") +
        annotate(geom="text", label=paste("(", round(xMin1, 2), "; ", round(yMin1, 2), "); Pr(Y|X<", round(xMin1, 2), ") = ", round(B6subProp, 2), sep=""), x=xMin1, y=yMid,  angle = 90, vjust=-1) +
        annotate(geom="text", label=paste("(", round(xMin2, 2), "; ", round(yMin2, 2), "); Pr(Y|X<", round(xMin2, 2), ") = ", round(SpSubProp, 2), sep=""), x=xMin2, y=yMid,  angle = 90, vjust=-1)
    }    
  }

  # Histograms
  #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
  for (whichChrs in c("autosomes", "diplosomes", "chrX") ) {
    # whichChrs = "chrX"
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[["unadjusted"]] = as.data.frame(geneLevelAScovAndRatios[[cellLine]]$altProp[ minCovIndices[[cellLine]] ])
      allelicPropVector[["adjusted"]] = as.data.frame(geneLevelAScovAndRatios[[cellLine]]$adjAltProp[ minCovIndices[[cellLine]] ])
    } else {
      allelicPropVector[["unadjusted"]] = as.data.frame(geneLevelAScovAndRatios[[cellLine]]$altProp[ intersect(chrIndices[[whichChrs]], minCovIndices[[cellLine]]) ])
      allelicPropVector[["adjusted"]] = as.data.frame(geneLevelAScovAndRatios[[cellLine]]$adjAltProp[ intersect(chrIndices[[whichChrs]], minCovIndices[[cellLine]]) ])
    }
    names(allelicPropVector[["unadjusted"]]) = "ratio"
    names(allelicPropVector[["adjusted"]]) = "ratio"
    
    list.names <- names(allelicPropVector)
    lns <- sapply(allelicPropVector, nrow)
    dat <- as.data.frame(do.call("rbind", allelicPropVector))
    colnames(dat)[1] = "ratio"
    dat$group <- rep(list.names, lns)
    
    theCols2use = c(sampleCols[cellLine], sampleColsLightD8[cellLine])
    names(theCols2use) = list.names
    
    ggplotList[[paste(whichChrs, "_hist", sep="")]] <- ggplot(dat, aes(ratio, fill=group)) + 
      geom_histogram(breaks=seq(0, 1, by=0.05), alpha=0.6, position='dodge',
                     aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                             ..count..[..group..==2]/sum(..count..[..group..==2]))*100)) +
      ggtitle(paste("chrX zoomed\n(n=", lns[1], ")", sep="")) +
      scale_fill_manual(values=theCols2use) +
      labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])"), y="%") + 
      theme(legend.position="top") +
      theme(legend.title=element_blank())

    # Get B6 and Sp-specifc autosomal minima
    if ( whichChrs == "autosomes" || whichChrs == "diplosomes"  ) {
      xMin1 = autosomalMinima[[cellLine]][[whichChrs]][["xMin1"]]
      yMin1 = autosomalMinima[[cellLine]][[whichChrs]][["yMin1"]] 
      xMin2 = autosomalMinima[[cellLine]][[whichChrs]][["xMin2"]] 
      yMin2 = autosomalMinima[[cellLine]][[whichChrs]][["yMin1"]] 
      
      B6autoProp = autosomalProps[[cellLine]][[whichChrs]][["B6autoProp"]]
      SpAutoProp = autosomalProps[[cellLine]][[whichChrs]][["SpAutoProp"]]
      
      ggplotData = print(ggplotList[[paste(whichChrs, "_hist", sep="")]])
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
      
      ggplotList[[paste(whichChrs, "_hist", sep="")]] = ggplotList[[paste(whichChrs, "_hist", sep="")]] +
        geom_vline(xintercept=xMin1, linetype="dashed") +
        geom_vline(xintercept=xMin2, linetype="dashed") +
        annotate(geom="text", label=paste("(", round(xMin1, 2), "; ", round(yMin1, 2), "); Pr(Y|X<", round(xMin1, 2), ") = ", round(B6autoProp, 2), sep=""), x=xMin1, y=yMid,  angle = 90, vjust=-1) +
        annotate(geom="text", label=paste("(", round(xMin2, 2), "; ", round(yMin2, 2), "); Pr(Y|X>", round(xMin2, 2), ") = ", round(SpAutoProp, 2), sep=""), x=xMin2, y=yMid,  angle = 90, vjust=+1.5)
    }
    
    if ( whichChrs == "chrX" ) {
      xMin1 = autosomalMinima[[cellLine]][["diplosomes"]][["xMin1"]]
      yMin1 = autosomalMinima[[cellLine]][["diplosomes"]][["yMin1"]] 
      xMin2 = autosomalMinima[[cellLine]][["diplosomes"]][["xMin2"]] 
      yMin2 = autosomalMinima[[cellLine]][["diplosomes"]][["yMin1"]] 
  
      B6subProp = chrXSubProps[[cellLine]][["B6subProp"]]
      SpSubProp = chrXSubProps[[cellLine]][["SpSubProp"]]          

      ggplotData = print(ggplotList[[paste(whichChrs, "_hist", sep="")]])
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
      
      ggplotList[[paste(whichChrs, "_hist", sep="")]] = ggplotList[[paste(whichChrs, "_hist", sep="")]] +
        geom_vline(xintercept=xMin1, linetype="dashed") +
        geom_vline(xintercept=xMin2, linetype="dashed") +
        annotate(geom="text", label=paste("(", round(xMin1, 2), "; ", round(yMin1, 2), "); Pr(Y|X<", round(xMin1, 2), ") = ", round(B6subProp, 2), sep=""), x=xMin1, y=yMid,  angle = 90, vjust=-1) +
        annotate(geom="text", label=paste("(", round(xMin2, 2), "; ", round(yMin2, 2), "); Pr(Y|X<", round(xMin2, 2), ") = ", round(SpSubProp, 2), sep=""), x=xMin2, y=yMid,  angle = 90, vjust=-1)

      ggplotList[["chrXzoom_hist"]] <- ggplot(dat, aes(ratio, fill=group)) + 
        geom_histogram(breaks=seq(0, 1, by=0.05), alpha=0.6, position='dodge',
                       aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                               ..count..[..group..==2]/sum(..count..[..group..==2]))*100)) +
        ggtitle(paste("chrX zoomed\n(n=", lns[1], ")", sep="")) +
        scale_fill_manual(values=theCols2use) +
        labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])"), y="%") + 
        theme(legend.position="top") +
        theme(legend.title=element_blank()) + 
        coord_cartesian(ylim = c(0, 5))

      ggplotData = print(ggplotList[["chrXzoom_hist"]])
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
      
      ggplotList[["chrXzoom_hist"]] =  ggplotList[["chrXzoom_hist"]] + 
        geom_vline(xintercept=xMin1, linetype="dashed") +
        geom_vline(xintercept=xMin2, linetype="dashed") +
        annotate(geom="text", label=paste("(", round(xMin1, 2), "; ", round(yMin1, 2), "); Pr(Y|X<", round(xMin1, 2), ") = ", round(B6subProp, 2), sep=""), x=xMin1, y=yMid,  angle = 90, vjust=-1) +
        annotate(geom="text", label=paste("(", round(xMin2, 2), "; ", round(yMin2, 2), "); Pr(Y|X<", round(xMin2, 2), ") = ", round(SpSubProp, 2), sep=""), x=xMin2, y=yMid,  angle = 90, vjust=-1)
    }    
  }
  
  PDFfilename = file.path(ASEdir, 
                          paste(workDir, ".", cellLine, ".densityOverlaid.unadjustedVsAdjsutedAllelicRatios.minCov=", minCov, ".pdf", sep=""))
  pdf(PDFfilename, height=12, width=20)
  #par(mfrow=c(1,2), oma=c(0,0,2.5,0))
  
  grid.arrange(ggplotList[[1]], ggplotList[[2]], ggplotList[[3]], ggplotList[[4]], 
               ggplotList[[5]], ggplotList[[6]], ggplotList[[7]], ggplotList[[8]], 
               nrow=2, ncol=4, 
               top=textGrob(paste(cellLine, " distributions of ", list.names[1], ", and ", list.names[2], " allelic ratios (minCov=", minCov, ")", sep=""),
                            gp=gpar(fontsize=16)))
  dev.off()
  
}









###########################################################
###########################################################
# Dendrograms and heat maps

ClusterDir = file.path(outputDir, "sampleClustering")
if ( ! file.exists( ClusterDir ) ) {
  if ( ! dir.create( ClusterDir, showWarnings=FALSE, recursive=TRUE )  ) ClusterDir = outputDir
}

subSetSelections = subSetRenames = list()

for (expDataType in names(expressionTables)) {
  
  for( subSetName in c("all", names(subSetSelections))) {
    if (subSetName != "all") {
      thisNormExpdata = expressionTables[[expDataType]][,subSetSelections[[subSetName]]]
      colnames(thisNormExpdata) = subSetRenames[[subSetName]]
      thisCountData = countData[,subSetSelections[[subSetName]]]
      colnames(thisCountData) = subSetRenames[[subSetName]]
      theseSampleCols = sampleCols[subSetSelections[[subSetName]]]
      ClusterDir = file.path(outputDir, paste("sampleClustering/", subSetName, sep="") )
      if ( ! file.exists( ClusterDir ) ) {
        if ( ! dir.create( ClusterDir, showWarnings=FALSE, recursive=TRUE )  ) ClusterDir = outputDir
      }
    } else {
      thisNormExpdata = expressionTables[[expDataType]]
      thisCountData = countData
      theseSampleCols = sampleCols
      ClusterDir = file.path(outputDir, "sampleClustering")
      if ( ! file.exists( ClusterDir ) ) {
        if ( ! dir.create( ClusterDir, showWarnings=FALSE, recursive=TRUE )  ) ClusterDir = outputDir
      }    
    }
    logNormExpdata <- addThenLog(thisNormExpdata)

    #for (whichChrs in c("genome-wide", "chrX", "autosomes") ) {
    for (whichChrs in c( "chrX", "autosomes") ) {
      
      if (whichChrs == "genome-wide") {
        expTable = logNormExpdata
        countTable = thisCountData
      } else {
        expTable = logNormExpdata[chrIndices[[whichChrs]],]
        countTable = thisCountData[chrIndices[[whichChrs]],]
      }
      
      #######################
      # Correlation across samples
      cor_matrix  <- cor(expTable, use="pairwise.complete")
      
      #     # Dendrogram V1
      #     hclust_results_asDistCor  <- hclust(as.dist(1-cor_matrix))
      #     PDFfilename = file.path(ClusterDir, paste(workDir, "_hclust_asdistCorr_log2", expDataType, "_", whichChrs, ".pdf", sep=""))
      #     pdf(PDFfilename, height=6, width=6)
      #     plot(hclust_results_asDistCor,hang=-1,
      #          main=paste("Hierarchical clust by correlation of expression (log2 ", expDataType, ")\n", whichChrs, " genes", sep=""))
      #     dev.off()
      #     
      #     # Dendrogram V2
      #     hclust_results_distCor  <- hclust(dist(1-cor_matrix))
      #     PDFfilename = file.path(ClusterDir, paste(workDir, "_hclust_distCorr_log2", expDataType, "_", whichChrs, ".pdf", sep=""))
      #     pdf(PDFfilename, height=10, width=10)
      #     plot(hclust_results_distCor,hang=-1, 
      #          main=paste("Hierarchical clust by correlation of expression (log2 ", expDataType, ")\n", whichChrs, " genes", sep=""))
      #     dev.off()
      #     
      #     # Cut-Tree Analysis
      #     cutree_results  <- cutree(hclust_results_distCor, k=3)
      #     # Cluster groups colours
      #     cols <- brewer.pal(max(cutree_results), "Accent")
      
      # total count colors
      hmcol       <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
      newcols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
      W2Rpalette  <- colorRampPalette(c("white", "red"))
      W2Rcols     <- W2Rpalette(ncol(expTable))
      totalCounts <- colSums(countTable)
      sizeColors  <- W2Rcols[rank(totalCounts)]
      sizeText    <- as.character(totalCounts)
      sizeText[totalCounts>=100000]  <- paste(as.character(round(totalCounts[totalCounts>=100000]/1000)), "K", sep="")
      sizeText[totalCounts>=1000000] <- paste(as.character(round(totalCounts[totalCounts>=1000000]/1000000, 2)), "M", sep="")
      
      PDFfilename = file.path(ClusterDir, paste(workDir, "_corrMatrix_log2", expDataType, "_", whichChrs, ".pdf", sep=""))
      pdf(PDFfilename, height=10, width=10)
      heatmap.2(cor_matrix, distfun = function(x) dist(1-x), trace="none", col=(hmcol), margin=c(20, 20),
                ColSideColors=theseSampleCols, # cols[cutree_results], 
                RowSideColors=sizeColors, 
                labRow=paste(colnames(expTable), sizeText), labCol=paste(colnames(expTable), sizeText),
                main=paste("Correlation heat map of expression (log2 ", expDataType, ")\n", 
                           whichChrs, " (Samples: ", ncol(expTable)," Genes: ",nrow(expTable),")", sep=""))
      dev.off()
    }
  }
}



###########################################################
###########################################################
# CDFs

###########################################################
# ChrByChr CDFs: WT, Del1, and Del2

CDFdir = file.path(outputDir, "CDFs")
if ( ! file.exists( CDFdir ) ) {
  if ( ! dir.create( CDFdir, showWarnings=FALSE, recursive=TRUE )  ) CDFdir = outputDir
}

subSetSelections = subSetRenames = list()

subSetSelections[["WT-vs-Del1"]] = c(
  "PatskiWT",
  "PatskiWT_0uMaza_rep1",
  "PatskiWT_0uMaza_rep2",  
  "PatskiWT_4uMaza_rep2",  
  "PatskiWT_4uMaza_rep1",
  "PatskiDel1",
  "PatskiDel1b",
  "PatskiDel1_0uMaza_rep1",
  "PatskiDel1_0uMaza_rep2",
  "PatskiDel1_4uMaza_rep2",
  "PatskiDel1_4uMaza_rep1") 

subSetRenames[["WT-vs-Del1"]] = c("WT",
                                 "WT_0uMaza_rep1", 
                                 "WT_0uMaza_rep2",
                                 "WT_4uMaza_rep1",
                                 "WT_4uMaza_rep2",
                                 "Del1",
                                 "Del1b",
                                 "Del1_0uMaza_rep1", 
                                 "Del1_0uMaza_rep2",
                                 "Del1_4uMaza_rep1",
                                 "Del1_4uMaza_rep2")

#for( subSetName in c("all", names(subSetSelections))) {
for( subSetName in names(subSetSelections)) {
#for( subSetName in "newData_NoDel5") {
  # subSetName = "newData_NoDel5"
  # subSetName = "newData"
  thisExpressionTables = list()
  if (subSetName != "all") {
    for (expDataType in names(expressionTables) )  {
      thisExpressionTables[[expDataType]] = expressionTables[[expDataType]][,subSetSelections[[subSetName]]]
      colnames(thisExpressionTables[[expDataType]]) = subSetRenames[[subSetName]]
    }
    theseSampleCols = sampleCols[subSetSelections[[subSetName]]]
    CDFdir = file.path(outputDir, paste("CDFs/", subSetName, sep="") )
    if ( ! file.exists( CDFdir ) ) {
      if ( ! dir.create( CDFdir, showWarnings=FALSE, recursive=TRUE )  ) CDFdir = outputDir
    }
  } else {
    thisExpressionTables[[expDataType]] = expressionTables[[expDataType]]
    theseSampleCols = sampleCols
    CDFdir = file.path(outputDir, "CDFs")
    if ( ! file.exists( CDFdir ) ) {
      if ( ! dir.create( CDFdir, showWarnings=FALSE, recursive=TRUE )  ) CDFdir = outputDir
    }    
  }
  
  # RPKMS
  for (expNormType in c("RPKM", "TPM")) {
    print(expNormType)
    
    # Quantile normalize? 
    for (Qnorm in c(FALSE, TRUE)) {
      # ONLY makes (dubious) sense for RPKMs 
      if (Qnorm && expNormType == "TPM") next
      print(Qnorm)
      
      for (whichChrs in names(chrBoolindices) ) {
        
        print(whichChrs)
        
        # CDFs between samples
        PDFfilename = file.path(CDFdir, paste(workDir, "_CDFs_log2_", whichChrs, 
                                              ".", expNormType, ifelse(Qnorm, ".Qnorm",""), ".pdf", sep=""))
        pdf(PDFfilename, height=6, width=12)
        par(mfrow=c(1,3), oma=c(0,0,2.5,0))
        
        for (expDataType in grep(expNormType, names(thisExpressionTables), value=TRUE)) {
          
          entireExpTable = thisExpressionTables[[expDataType]]
          
          # Quant norm?
          if (Qnorm) {
            # *** Qnorm NOT WORKING B/C altRPKM and refRPKM have NAs
            NonNArows = apply(entireExpTable, 1, function(x){ !any(is.na(x))})
            preprocExpTable = entireExpTable[NonNArows, ]
            preprocExpTable = quantNormGB(as.matrix(preprocExpTable))
          } else {
            NonNArows = rep(TRUE, nrow(entireExpTable)) # Needed for subsetting chrBoolindices[[whichChrs]]
            preprocExpTable = entireExpTable
          }
          
          expTable = preprocExpTable[chrBoolindices[[whichChrs]][NonNArows], ]
          
          # Log data      
          expTable <- addThenLog(expTable)
          
          if ( grepl("RPKM", expDataType) ) {
            maxExp = log2(64) # 6
            xlims = c(0, maxExp)
            if ( whichChrs == "chrX" && expDataType == "refRPKM") {
              ylims=c(0.6, 1)        
            } else {
              ylims=c(0.15, 1)        
            }
          } else { # TPM
            maxExp = log2(10^3)
            xlims = c(0, maxExp)
            if ( whichChrs == "chrX" && expDataType == "refTPM") {
              ylims=c(0.6, 1)        
            } else {
              ylims=c(0.15, 1)        
            }
          }
          
          colName = colnames(expTable)[1]
          plot(ecdf(expTable[,1]), verticals=TRUE, do.points=FALSE,
               main=expDataType, 
               ylab=paste("Percentage of genes"), xlab=paste("Expression level (log2[", expDataType, " + 1])", sep=""),
               lwd=3, col=theseSampleCols[1], xaxt="n", yaxt="n",
               xlim=xlims, ylim=ylims)
          axis(1, at=pretty(c(0, maxExp)), labels=round(2^pretty(c(0, maxExp))))  
          axis(2, at=seq(0,1,0.1), labels=seq(0,100,10), las=2)  
          for (colidx in 2:ncol(expTable)) {
            colName = colnames(expTable)[colidx]
            plot(ecdf(expTable[,colidx]), verticals=TRUE, do.points=FALSE,
                 lwd=3, col=theseSampleCols[colidx], 
                 xaxt="n", yaxt="n", add=TRUE)
          }
          legend( "bottomright", legend=colnames(expTable), fill=theseSampleCols)
        } # expDataType
        title( main=paste( "CDFs of ", ifelse(Qnorm, "Quantile-normalized ","raw " ), "expression", "\n",
                           whichChrs, " (Genes: ", nrow(expTable),")", sep="" ), outer=TRUE ) 
        dev.off()   
      } # whichChrs
    } # Qnorm
  } # RPKM or TPM  
  
} # subselection






###########################################################
###########################################################
# Chromosome-wide expression:
# Lots of evidence for aneuploidy seen by CDFs
# This confirms it by looking at total/average chromosomal expression.

aneuploidyDir = file.path(outputDir, "aneuploidyCheck")
if ( ! file.exists( aneuploidyDir ) ) {
  if ( ! dir.create( aneuploidyDir, showWarnings=FALSE, recursive=TRUE )  ) aneuploidyDir = outputDir
}

theChroms = paste("chr", c(seq(1, 19), "X"), sep="")

sampleExpByChrTable = list() 
sampleExpByChrByWTtable = list() 
sampleExpByChrByWTlog2table = list() 
sampleExpPropByChrTable = list()
sampleExpPropByChrByWTtable = list()
sampleExpPropByChrByWTlog2table = list()
for (ExpNormType in names(expressionTables)) {
  for (cellLine in cellLines) {
    sampleExpByChr=c()
    for (currChr in theChroms) {
      currChrIdx = which(as.character(mm10geneChrs$chr) == currChr)
      chrExp = sum(expressionTables[[ExpNormType]][currChrIdx, cellLine], na.rm=TRUE)
      sampleExpByChr = c( sampleExpByChr, chrExp )
    }
    sampleExpByChrTable[[ExpNormType]] = cbind( sampleExpByChrTable[[ExpNormType]], sampleExpByChr )
  }  
  colnames(sampleExpByChrTable[[ExpNormType]]) = cellLines
  rownames(sampleExpByChrTable[[ExpNormType]]) = theChroms
  
  sampleExpByChrByWTtable[[ExpNormType]] = round(sampleExpByChrTable[[ExpNormType]] / sampleExpByChrTable[[ExpNormType]][,1], 2)
  sampleExpByChrByWTlog2table[[ExpNormType]] = round(log2(sampleExpByChrTable[[ExpNormType]] / sampleExpByChrTable[[ExpNormType]][,1]), 2)
  sampleExpPropByChrTable[[ExpNormType]] = round(sampleExpByChrTable[[ExpNormType]] / colSums(sampleExpByChrTable[[ExpNormType]]) * 100, 2)
  sampleExpPropByChrByWTtable[[ExpNormType]] = round(sampleExpPropByChrTable[[ExpNormType]] / sampleExpPropByChrTable[[ExpNormType]][,1], 2)
  sampleExpPropByChrByWTlog2table[[ExpNormType]] = round(log2(sampleExpPropByChrTable[[ExpNormType]] / sampleExpPropByChrTable[[ExpNormType]][,1]), 2)
  
  fileName = file.path(aneuploidyDir, paste(workDir, "_exprPerChrPerSample_", ExpNormType, ".tsv", sep=""))
  write.table(sampleExpByChrTable[[ExpNormType]], fileName, row.names=TRUE, col.names=NA, sep="\t")  
  
  fileName = file.path(aneuploidyDir, paste(workDir, "_exprPerChrPerSampleRelativeToWT_", ExpNormType, ".tsv", sep=""))
  write.table(sampleExpByChrByWTtable[[ExpNormType]], fileName, row.names=TRUE, col.names=NA, sep="\t")  
  
  fileName = file.path(aneuploidyDir, paste(workDir, "_exprPerChrPerSampleRelativeToWTlog2_", ExpNormType, ".tsv", sep=""))
  write.table(sampleExpByChrByWTlog2table[[ExpNormType]], fileName, row.names=TRUE, col.names=NA, sep="\t")  
  
  fileName = file.path(aneuploidyDir, paste(workDir, "_propOfExprPerChrPerSample_", ExpNormType, ".tsv", sep=""))
  write.table(sampleExpPropByChrTable[[ExpNormType]], fileName, row.names=TRUE, col.names=NA, sep="\t")  
  
  fileName = file.path(aneuploidyDir, paste(workDir, "_propOfExprPerChrPerSampleRelativeToWT_", ExpNormType, ".tsv", sep=""))
  write.table(sampleExpPropByChrByWTtable[[ExpNormType]], fileName, row.names=TRUE, col.names=NA, sep="\t")  
  
  fileName = file.path(aneuploidyDir, paste(workDir, "_propOfExprPerChrPerSampleRelativeToWTlog2_", ExpNormType, ".tsv", sep=""))
  write.table(sampleExpPropByChrByWTlog2table[[ExpNormType]], fileName, row.names=TRUE, col.names=NA, sep="\t")  
}


###
# RPKM haploid/diploid ratios per chromosome

###
# alt/2n and ref/2n: RPKM
sampleExpByChrRefVsDiploidTable = round(sampleExpByChrTable[["refRPKM"]] / sampleExpByChrTable[["RPKM"]], 2)
sampleExpByChrAltVsDiploidTable = round(sampleExpByChrTable[["altRPKM"]] / sampleExpByChrTable[["RPKM"]], 2)

fileName = file.path(aneuploidyDir, paste(workDir, "_exprPerChrRefVsDiploidPerSample_RPKM.tsv", sep=""))
write.table(sampleExpByChrRefVsDiploidTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  

fileName = file.path(aneuploidyDir, paste(workDir, "_exprPerChrAltVsDiploidPerSample_RPKM.tsv", sep=""))
write.table(sampleExpByChrAltVsDiploidTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  



###
# alt/(alt+ref) and ref/(alt+ref): RPKM
# This masks the effect of reduced SNP coverage in some chromosomes.
sampleExpByChrRefVsDiploidTable = round(sampleExpByChrTable[["refRPKM"]] / (sampleExpByChrTable[["altRPKM"]] + sampleExpByChrTable[["refRPKM"]]), 2)
sampleExpByChrAltVsDiploidTable = round(sampleExpByChrTable[["altRPKM"]] / (sampleExpByChrTable[["altRPKM"]] + sampleExpByChrTable[["refRPKM"]]), 2)

fileName = file.path(aneuploidyDir, paste(workDir, "_exprPerChrRefVsAltPlusRefPerSample_RPKM.tsv", sep=""))
write.table(sampleExpByChrRefVsDiploidTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  

fileName = file.path(aneuploidyDir, paste(workDir, "_exprPerChrAltVsAltPlusRefPerSample_RPKM.tsv", sep=""))
write.table(sampleExpByChrAltVsDiploidTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  



###
# TPM haploid/diploid ratios per chromosome should give same resuls as RPKMs b/c they represent the same proportions of total

###
# alt/2n and ref/2n: TPM
sampleExpByChrRefVsDiploidTable = round(sampleExpByChrTable[["refTPM"]] / sampleExpByChrTable[["TPM"]], 2)
sampleExpByChrAltVsDiploidTable = round(sampleExpByChrTable[["altTPM"]] / sampleExpByChrTable[["TPM"]], 2)

fileName = file.path(aneuploidyDir, paste(workDir, "_exprPerChrRefVsDiploidPerSample_TPM.tsv", sep=""))
write.table(sampleExpByChrRefVsDiploidTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  

fileName = file.path(aneuploidyDir, paste(workDir, "_exprPerChrAltVsDiploidPerSample_TPM.tsv", sep=""))
write.table(sampleExpByChrAltVsDiploidTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  


###
# alt/(alt+ref) and ref/(alt+ref): TPM
# This masks the effect of reduced SNP coverage in some chromosomes.
sampleExpByChrRefVsDiploidTable = round(sampleExpByChrTable[["refTPM"]] / (sampleExpByChrTable[["altTPM"]] + sampleExpByChrTable[["refTPM"]]), 2)
sampleExpByChrAltVsDiploidTable = round(sampleExpByChrTable[["altTPM"]] / (sampleExpByChrTable[["altTPM"]] + sampleExpByChrTable[["refTPM"]]), 2)

fileName = file.path(aneuploidyDir, paste(workDir, "_exprPerChrRefVsAltPlusRefPerSample_TPM.tsv", sep=""))
write.table(sampleExpByChrRefVsDiploidTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  

fileName = file.path(aneuploidyDir, paste(workDir, "_exprPerChrAltVsAltPlusRefPerSample_TPM.tsv", sep=""))
write.table(sampleExpByChrAltVsDiploidTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  









###########################################################
###########################################################
# Expression for specific genes

geneIDlist = c(
"1110012L19Rik",
"A830080D01Rik",
"Apoo",
"BC023829",
"Bcorl1",
"Brwd3",
"Ccdc120",
"Cited1",
"Dkc1",
"Fam120c",
"Fgd1",
"Foxo4",
"Gpc3",
"Kctd12b",
"Mbtps2",
"Mcts1",
"Morc4",
"Morf4l2",
"Mtmr1",
"Ndufb11",
"Nhs",
"Nkap",
"Ophn1",
"Prps1",
"Rab9",
"Rnf113a1",
"Slc10a3",
"Tsc22d3",
"Upf3b",
"Zdhhc9",
"Zmat1",
"Apex2",
"Atrx",
"B630019K06Rik",
"Cd99l2",
"Cenpi",
"Diap2",
"Elf4",
"Elk1",
"Gprasp1",
"Mmgt1",
"Mospd1",
"Phka1",
"Ppp1r3f",
"Prps2",
"Psmd10",
"Rap2c",
"Rps6ka3",
"Shroom2",
"Slc9a6",
"Suv39h1",
"Tbc1d25",
"Tmem164",
"Trmt2b",
"Tspyl2",
"Usp11",
"Utp14a",
"Xiap",
"Zfp185",
"Acsl4",
"Ammecr1",
"Bhlhb9",
"Cldn2",
"Gpkow",
"Hccs",
"Med12",
"Ogt",
"Pdzd4",
"Ripply1",
"Wdr45",
"Hcfc1",
"Hprt",
"Rbm10",
"Tfe3",
"Tspan6",
"Prkx",
"Sms",
"Cdk16",
"Las1l",
"Prdx4",
"Taf1",
"C1galt1c1",
"Magee1",
"Rbbp7",
"Pgrmc1",
"Tspan7",
"Zfp275",
"Fhl1",
"Maoa"
)
  
for (geneID in geneIDlist) {
  # geneID = "Apoo"
  print(geneLevelASEdata.uniqued[geneLevelASEdata.uniqued$geneID==geneID,c(c(1:5), 
                                                                           grep("sumOverGenes", colnames(geneLevelASEdata.uniqued)) )])
}
# chr    start      end geneID Null
# 28491 chrX 94367109 94417092   Apoo    

# geneLevelASEdata
tmpTable = geneLevelASEdata[, grepl("sumOverGenes", colnames(geneLevelASEdata))]
shortColnames = gsub(".sumOverGenes", "", colnames(tmpTable))
colnames(tmpTable) = shortColnames
#rownames(tmpTable) = geneLevelASEdata$geneID
head(tmpTable)

# geneLevelASEdata.uniqued
tmpTable.uniqued = geneLevelASEdata.uniqued[, grepl("sumOverGenes", colnames(geneLevelASEdata.uniqued))]
shortColnames = gsub(".sumOverGenes", "", colnames(tmpTable.uniqued))
colnames(tmpTable.uniqued) = shortColnames
rownames(tmpTable.uniqued) = geneLevelASEdata.uniqued$geneID
head(tmpTable.uniqued)

for (geneID in geneIDlist) {
#   print("ASE")
#   print(tmpTable[geneLevelASEdata$geneID==geneID, ])
  print("ASE uniqued")
  print(tmpTable.uniqued[geneID, ])

#   print("altProp")
#   print(geneLevelASEdata.uniqued.altProp[geneID, ])
#   print("adjAltProp")

  print(geneLevelASEdata.uniqued.adjAltProp[geneID, ])
  for (expDataType in grep("RPKM", names(expressionTables), value=TRUE)) {
    print(expDataType)
    thisGeneExpdata = expressionTables[[expDataType]][geneID,]
    print(thisGeneExpdata)
  }

  readkey()
}






###########################################################
###########################################################
# Differentially expressed genes BASED ON GENOME-WIDE NORMALIZATION.
# But because aneuploidy is so bad, one really needs to use chromosome normalized values (see below)


diffExpGenes = file.path(outputDir, "differentialExpression")
if ( ! file.exists( diffExpGenes ) ) {
  if ( ! dir.create( diffExpGenes, showWarnings=FALSE, recursive=TRUE )  ) diffExpGenes = outputDir
}

#FCthresh = log2(1.5) # 50% more i.e. 0.5849625 fold change
FCthresh = log2(2) # 100% more i.e. 1 fold change

diffexpTables = list()
for (expDataType in names(expressionTables)) {
  thisNormExpdata = expressionTables[[expDataType]]
  #logNormExpdata <- addThenLog(thisNormExpdata)  
  
  if (grepl("RPKM", expDataType)) {
    minNormExpValue = 1
  } else { # TPM
    minNormExpValue = 1
  }
  
  diffexpResults = cellCompDescs = c()
  for ( cellLineA in cellLines) {
    for ( cellLineB in cellLines) {
      cellCompDesc = paste(cellLineA, "-", cellLineB, sep="")
      cellCompDescRev = paste(cellLineB, "-", cellLineA, sep="")
      if (cellLineA == cellLineB || 
          ( cellCompDescRev %in% cellCompDescs ) ) {
        next
      }
      diffExpBoolidx = ( abs( log2( (thisNormExpdata[,cellLineA] + 0.001) / 
                                      (thisNormExpdata[,cellLineB] + 0.001) ) ) >= FCthresh ) &
        (thisNormExpdata[,cellLineA] >= minNormExpValue | thisNormExpdata[,cellLineB] >= minNormExpValue)
      diffexpResults = cbind(diffexpResults, diffExpBoolidx)
      cellCompDescs = c(cellCompDescs, cellCompDesc)
      cat(expDataType, "\t", cellCompDesc, "\t", 
          sum(diffExpBoolidx, na.rm=TRUE), "\t", sum(diffExpBoolidx, na.rm=TRUE)/length(diffExpBoolidx), 
          "\n", sep="")
    }
  }
  colnames(diffexpResults) = cellCompDescs
  rownames(diffexpResults) = rownames(thisNormExpdata)
  diffexpTables[[expDataType]] = diffexpResults
}


###
# Stats
for (expDataType in names(expressionTables)) {
  for (cellCompDesc in colnames(diffexpTables[[expDataType]])) {
    #for (theCurr in names(chrIndices)) {
    for (theCurr in c("autosomes", "chrX")) {
      cat(expDataType, "\t", cellCompDesc, "\t", theCurr, "\t", 
          sum(diffexpTables[[expDataType]][,cellCompDesc][chrIndices[[theCurr]]], na.rm=TRUE), "\t", 
          round( sum(diffexpTables[[expDataType]][,cellCompDesc][chrIndices[[theCurr]]], na.rm=TRUE) / 
                   length(chrIndices[[theCurr]]) * 100, 2), 
          "\n", sep="")
    }
  }
}






###
# Xi vs. autosomes

# cellCompDesc="PatskiWT-PatskiDel1"
#cellCompDesc="PatskiDel1_0uMaza_rep1-PatskiDel1_4uMaza_rep1" 

for(expDataType in c("TPM", "altTPM", "refTPM")) {
  
  theDiffGenes = diffexpTables[[expDataType]][,cellCompDesc]
  
  celltypeA = strsplit(cellCompDesc, "-")[[1]][1]
  celltypeB = strsplit(cellCompDesc, "-")[[1]][2]
  theGeneXpressionA = expressionTables[[expDataType]][,celltypeA]
  theGeneXpressionB = expressionTables[[expDataType]][,celltypeB]
  
  for (whichChrs in c("autosomes", "chrX") ) {
    cat(expDataType, "\t", whichChrs, "\t", sep="")
    cat(sum(theDiffGenes & chrBoolindices[[whichChrs]], na.rm=TRUE), "\t", sep="")
    #cat(sum(TPMexpressedGenesBoolIdx & chrBoolindices[[whichChrs]], na.rm=TRUE), "\t", sep="")
    cat(sum((theGeneXpressionA >= 1 | theGeneXpressionB >=1) & chrBoolindices[[whichChrs]], na.rm=TRUE), "\t", sep="")
    cat(sum(theDiffGenes & chrBoolindices[[whichChrs]], na.rm=TRUE) /
          #sum(TPMexpressedGenesBoolIdx & chrBoolindices[[whichChrs]], na.rm=TRUE) 
          sum((theGeneXpressionA >= 1 |  theGeneXpressionB >=1 ) & chrBoolindices[[whichChrs]], na.rm=TRUE)
        * 100, "\n", sep="")
  }
}


###
# Differential gene lists

saveXLSX=TRUE
#saveXLSX=FALSE
saveBED=FALSE

for ( cellLineA in cellLines) {
  for ( cellLineB in cellLines) {
    cellCompDesc = paste(cellLineA, "-", cellLineB, sep="")
    if (cellLineA == cellLineB) {
      next
    }
    print(paste("****** ", cellCompDesc, " ******", sep=""))
    
    #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
    for (whichChrs in c("chrX") ) {
      # whichChrs = "autosomes" 
      # expDataType = "altRPKM"; thisNormExpdata = expressionTables[[expDataType]]; cellLineA="PatskiWT"; cellLineB="PatskiDel1"; cellCompDesc = paste(cellLineA, "-", cellLineB, sep="")
      
      diffTableList = list()  
      for (expDataType in names(expressionTables)) {
        print(paste("       ", expDataType, "       ", sep=""))
        
        thisNormExpdata = expressionTables[[expDataType]]
        #logNormExpdata <- addThenLog(thisNormExpdata)  
        
        if ( ! cellCompDesc %in% colnames(diffexpTables[[expDataType]]) ) {
          print("Next!")
          next
        } else {
          # print("Check!")
        }
        
        cellLineAdata = thisNormExpdata[, cellLineA]
        cellLineBdata = thisNormExpdata[, cellLineB]
        theDiffGenes = diffexpTables[[expDataType]][,cellCompDesc]
        diffDir = rep(cellLineA, length(cellLineAdata))
        diffDir[ cellLineAdata < cellLineBdata ] = cellLineB
        diffTable = cbind(rownames(thisNormExpdata),
                          round(cellLineAdata, 2), round(cellLineBdata, 2), 
                          round(log2( (cellLineAdata + 0.001) / (cellLineBdata + 0.001) ), 2),
                          diffDir )
        if (sum(which(theDiffGenes & chrBoolindices[[whichChrs]])) < 2) {
          next
        } 
        diffTable = diffTable[which(theDiffGenes & chrBoolindices[[whichChrs]]),]
        names(diffTable) = c("gene", cellLineA, cellLineB, "foldChange", "directionOfChange")
        #rownames(diffTable) = names(cellLineAdata)[which(theDiffGenes & chrBoolindices[[whichChrs]])]
        
        diffTableList[[expDataType]] = as.data.frame(diffTable)
      }
      
      if (length(names(diffTableList)) == 0 ) {
        next
      }
      
      if (saveXLSX) {
        # To excel spreadsheet
        #         names(diffTableList) =c("diploid.RPKM", "spretus-specific.RPKM", "black6-specific.RPKM",
        #                                 "diploid.TPM", "spretus-specific.TPM", "black6-specific.TPM")
        ouputfilename = file.path(diffExpGenes, 
                                  paste(workDir,
                                        "_DEgenes",
                                        "_FCthresh", FCthresh,
                                        "_avgMinExp", minNormExpValue,
                                        "_", cellCompDesc, "_",
                                        whichChrs,
                                        ".xlsx", sep=""))
        write.xlsx(diffTableList, file=ouputfilename, colNames=TRUE, rowNames=FALSE)
      }  
      
      
      if (saveBED) {
        # Save differential beds
        
        bedDir = file.path(outputDir, "diffbedgraphs")
        if ( ! file.exists( bedDir ) ) {
          if ( ! dir.create( bedDir, showWarnings=FALSE, recursive=TRUE )  ) bedDir = outputDir
        }
        
        for (expDataType in names(diffTableList)) {
          # diffType = "diploid.RPKM"
          
          cat(cellCompDesc, " ", whichChrs, " ", expDataType, "\n", sep="")
          
          bedDir = file.path(outputDir, "diffbedgraphs", expDataType, whichChrs)
          if ( ! file.exists( bedDir ) ) {
            if ( ! dir.create( bedDir, showWarnings=FALSE, recursive=TRUE )  ) bedDir = outputDir
          }
          
          geneIdx = match(diffTableList[[expDataType]]$gene, 
                          as.character(geneLevelASEdata.uniqued[,4]))
          
          bedTable = cbind( as.character(geneLevelASEdata.uniqued[geneIdx,1]),
                            format(geneLevelASEdata.uniqued[geneIdx,2], scientific=FALSE, trim=TRUE, justify="none"),
                            format(geneLevelASEdata.uniqued[geneIdx,3], scientific=FALSE, trim=TRUE, justify="none"),
                            format(diffTableList[[expDataType]]$foldChange, scientific=FALSE, trim=TRUE, justify="none"),
                            as.character(geneLevelASEdata.uniqued[geneIdx,4]) )
          
          bedFilename = file.path(bedDir, paste("DEgenes",
                                                "_", cellCompDesc,
                                                "_", expDataType,
                                                "_FCthresh", FCthresh,
                                                "_avgMinExp", minNormExpValue,
                                                "_", whichChrs,
                                                ".bedgraph", sep=""))
          bed2Filename = file.path(bedDir, paste("DEgenes",
                                                "_", cellCompDesc,
                                                "_", expDataType,
                                                "_FCthresh", FCthresh,
                                                "_avgMinExp", minNormExpValue,
                                                "_", whichChrs,
                                                ".noGeneId.bedgraph", sep=""))
          
          cat("track type=bed ",
              "name=\"", cellCompDesc, " DE genes\"", 
              "description=\"", cellCompDesc, " DE genes ", expDataType, " FCthresh", FCthresh, " avgMinExp", minNormExpValue, " ", whichChrs, "\"", 
              "visibility=full color=200,0,0 altColor=0,0,200 priority=20 viewLimits=-2:2 autoScale=off\n", 
              file=bedFilename, sep="")
          cat("track type=bed ",
              "name=\"", cellCompDesc, " DE genes\"", 
              "description=\"", cellCompDesc, " DE genes ", expDataType, " FCthresh", FCthresh, " avgMinExp", minNormExpValue, " ", whichChrs, "\"", 
              "visibility=full color=200,0,0 altColor=0,0,200 priority=20 viewLimits=-2:2 autoScale=off\n", 
              file=bed2Filename, sep="")
          
          write.table(bedTable, bedFilename, 
                      append = TRUE, sep = "\t", 
                      col.names = FALSE, row.names = FALSE, quote = FALSE)
          try(system(paste("gzip -f \"", sub("~/", "$HOME/", bedFilename), "\"", sep = "")))

          write.table(bedTable[,c(1:4)], bed2Filename, 
                      append = TRUE, sep = "\t", 
                      col.names = FALSE, row.names = FALSE, quote = FALSE)
          try(system(paste("gzip -f \"", sub("~/", "$HOME/", bed2Filename), "\"", sep = "")))
          
        }
      }
    }
  }     
}











###########################################################
###########################################################
# Scatter plots BASED ON GENOME-WIDE NORMALIZATION.
# But because aneuploidy is so bad, one really needs to use chromosome normalized values (see below)

scatterDir = file.path(outputDir, "scatterPlots")
if ( ! file.exists( scatterDir ) ) {
  if ( ! dir.create( scatterDir, showWarnings=FALSE, recursive=TRUE )  ) scatterDir = outputDir
}
  
###
# Individual autosome vs chrX scatter plots: adjusted RPKMs 

for (expDataType in names(expressionTables)) {
  thisNormExpdata = expressionTables[[expDataType]]
  logNormExpdata <- addThenLog(thisNormExpdata)  
  
  for ( cellLineA in cellLines) {
    for ( cellLineB in cellLines) {
      cellCompDesc = paste(cellLineA, "-", cellLineB, sep="")
      if (cellLineA == cellLineB || 
          ( ! cellCompDesc %in% colnames(diffexpTables[[expDataType]]) ) ) {
        next
      }
      
      # SNP read coverage distribution
      PDFfilename = file.path(scatterDir, 
                              paste(workDir, "_scatterplot_",
                                    expDataType, "_",
                                    cellCompDesc,
                                    ".autosomesVschrX.pdf", sep=""))
      pdf(PDFfilename, height=6, width=12)
      par(mfrow=c(1,2), oma=c(0,0,2.5,0))
      
      #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
      for (whichChrs in c("autosomes", "chrX") ) {
        # whichChrs = "autosomes"
        
        # pointCols = rep(rgb(0, 0, 0, alpha=0.5), length(cellLineAdata))
        # pointCols[ diffexpTables[[expDataType]][,cellCompDesc] ] = rgb(1, 0, 0, alpha=0.5)
        pointCols = rep(rgb(0, 0, 0, alpha=0), length(cellLineAdata))
        pointCols[ diffexpTables[[expDataType]][,cellCompDesc] ] = rgb(1, 0, 0, alpha=1)
        if (whichChrs == "genome-wide") {
          cellLineAdata = logNormExpdata[, cellLineA]
          cellLineBdata = logNormExpdata[, cellLineB]
          nDiffGenes = sum( diffexpTables[[expDataType]][,cellCompDesc], na.rm=TRUE )
        } else {
          cellLineAdata = logNormExpdata[chrIndices[[whichChrs]], cellLineA]
          cellLineBdata = logNormExpdata[chrIndices[[whichChrs]], cellLineB]
          pointCols = pointCols[ chrIndices[[whichChrs]] ]
          nDiffGenes = sum( diffexpTables[[expDataType]][,cellCompDesc][chrIndices[[whichChrs]]], na.rm=TRUE )
        }
        
        plot(cellLineAdata, cellLineBdata,
             #pch=20, col=pointCols,
             pch=22, bg=pointCols,
             main=paste(whichChrs, " (", nDiffGenes, " / ", length(cellLineAdata), " differential)", sep=""),
             xlab=paste(cellLineA, " log2(", expDataType, "+1)", sep=""),
             ylab=paste(cellLineB, " log2(", expDataType, "+1)", sep="")
        )
        legend("topleft", legend="differential", fill=rgb(1, 0, 0, alpha=1), bty="n")
        
      }
      
      title(  main=paste(expDataType, " - ", cellCompDesc,  sep=""), outer=TRUE )
      dev.off()
    }     
  }
}

  



###
# Pair plots: adjusted RPKMs 

for (expDataType in names(expressionTables)) {
  thisNormExpdata = expressionTables[[expDataType]]
  logNormExpdata <- addThenLog(thisNormExpdata)  
  
  for (expressed in c(FALSE, TRUE)) {
    # expressed=FALSE
    
    #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
    for (whichChrs in c("autosomes", "chrX") ) {
      # whichChrs = "autosomes"
      
      # cellLines = c("PatskiDel1", "PatskiDel2", "PatskiDel5", "PatskiWT")
      expTableList = list()
      for ( cellLine in cellLines) {
        
        if (expressed) {
          if (whichChrs == "genome-wide") {
            expTableList[[cellLine]] = as.data.frame(logNormExpdata[expressedGenesBoolIdx, cellLine])
          } else {
            expTableList[[cellLine]] = as.data.frame(logNormExpdata[intersect(which(expressedGenesBoolIdx),
                                                                           chrIndices[[whichChrs]]), cellLine])
          }
        } else {
          if (whichChrs == "genome-wide") {
            expTableList[[cellLine]] = as.data.frame(logNormExpdata[, cellLine])
          } else {
            expTableList[[cellLine]] = as.data.frame(logNormExpdata[chrIndices[[whichChrs]], cellLine])
          }
        }
        # colnames(expTableList[[cellLine]]) = "RPKM"  
        # expTableList[[cellLine]]$cellLine <- cellLine
        colnames(expTableList[[cellLine]]) = cellLine 
      }
      
      #and combine into your new data frame 
      NormExpdataframe = expTableList[[1]]
      for ( cellLineIdx in 2:length(cellLines)) {
        NormExpdataframe = cbind(NormExpdataframe, expTableList[[cellLineIdx]])
      }
      print(dim(NormExpdataframe))
        
      # SNP read coverage distribution
      PDFfilename = file.path(scatterDir, 
                              paste(workDir, "_scatterplotMatrices_",
                                    expDataType, "_",
                                    whichChrs, "_",
                                    ifelse(expressed, "expressedGenes", "allGenes"), 
                                    "_ALLlines.pdf", sep=""))
      pdf(PDFfilename, height=20, width=20)

      chart.Correlation.GB(NormExpdataframe, 
            main=paste(expDataType, " for ", ifelse(expressed, "expressed", "all"), " genes - ", whichChrs,  sep=""))
      
      #       ggpairs(NormExpdataframe)
      #       pm <- ggpairs(NormExpdataframe)
      #       pm <- ggtitle(paste(expDataType, " for ", ifelse(expressed, "expressed", "all"), " genes\n", whichChrs,  sep=""))
      #       pm
      
      #       ggplot(NormExpdataframe, aes(x=id, y=value)) + geom_point() + facet_grid(.~variable) +
      #         ggtitle(paste(expDataType, " for ", ifelse(expressed, "expressed", "all"), " genes\n", whichChrs,  sep=""))
      
      dev.off()
    }
  }     
}















###########################################################
###########################################################
# Escape gene expression gene list

escapeGeneDir = file.path(outputDir, "escapeGenes")
if ( ! file.exists( escapeGeneDir ) ) {
  if ( ! dir.create( escapeGeneDir, showWarnings=FALSE, recursive=TRUE )  ) escapeGeneDir = outputDir
}

SNPlevelASEfile = file.path(dataDir, "snpPileups/Patski.combo.SNPpileup.annotated.tsv.gz")
SNPlevelASEdata = read.table(gzfile(SNPlevelASEfile), header=TRUE, sep="\t")
head(SNPlevelASEdata[,1:15])
dim(SNPlevelASEdata)

# journals.plos.org/plosgenetics/article/asset?unique&id=info:doi/10.1371/journal.pgen.1005079.s010
escapeGenes=c(
  "Shroom4",
  "Otud5",
  "Hdac6",
  "Wdr13",
  "Ebp",
  "Dynlt3",
  "Bcor",
  "1810030O07Rik",
  "Ddx3x",
  "Kdm6a",
  "Il13ra1",
  "Zbtb33",
  "Lamp2",
  "Zfp280c",
  "Firre",
  "Hs6st2",
  "Htatsf1",
  "Fmr1",
  "Ids",
  "Slc6a8",
  "Idh3g",
  "Irak1",
  "Mecp2",
  "Flna",
  "Fam50a",
  "Ubl4a",
  "Fam3a",
  "G6pdx",
  "Ikbkg",
  "Fundc2",
  "Mtcp1",
  "Vbp1",
  "Pls3",
  "Tab3",
  "Gyk",
  "Eif2s3x",
  "Maged1",
  "Msn",
  "Yipf6",
  "Dlg3",
  "Snx12",
  "Zmym3",
  "Xist",
  "Rlim",
  "Abcb7",
  "5530601H04Rik",
  "Pbdc1",
  "Magt1",
  "Atp7a",
  "Sh3bgrl",
  "Pcdh19",
  "Cstf2",
  "Gla",
  "Bhlhb9",
  "Mid2",
  "Col4a5",
  "Amot",
  "Tmem29",
  "Gnl3l",
  "Iqsec2",
  "Kdm5c",
  "Rbbp7",
  "Ctps2",
  "Car5b",
  "Ofd1",
  "Mid1"
)

whichChrs = "chrX"


#############
# Autosomal allelic ratio distributions





#############
# Parameter estimation per cell line/sample
# chromosome level: sum(Xi refCov)/sum(Xi totalCov) for all genes
# gene level: mean(Xi refCov/Xi totalCov) for all genes <------------- *** Probably best, most relavent method ***
# SNP level: mean(Xi refCov/Xi totalCov) for all SNPs

expXiProp.chromLevel = list()

expXiProp.adj.chromLevel = list()

expXiProp.adj.geneLevel = list()
sdXiXa.adj.geneLevel = list()

expXiProp.adj.SNPlevel = list()
sdXiXa.adj.SNPlevel = list()

expXiCov.adj.geneLevel = list()
expXiTPM.adj.geneLevel = list()

quantsXiCov.adj.geneLevel = list()
quantsXiTPM.adj.geneLevel = list()

for(cellLine in cellLines) {
  print(cellLine)
  
  ######
  # chromosome level parameters
  
  ###
  # unadjusted expected Xi proportion: *** DON'T USE THIS. JUST FOR INTEREST ***
  chrXNi0 = sum(geneLevelASEdata.uniqued.refCov[chrIndices[[whichChrs]], cellLine], na.rm=TRUE)
  chrXNi1 = sum(geneLevelASEdata.uniqued.altCov[chrIndices[[whichChrs]], cellLine], na.rm=TRUE)
  chrXNi = sum(geneLevelASEdata.uniqued.totCov[chrIndices[[whichChrs]], cellLine], na.rm=TRUE)
  chrXNi ==  (chrXNi0 + chrXNi1)
  # TRUE
  
  # Expected Xi:Xa ratio based on chrX Sp:Bl6 from total (chromosome-wide) SNP coverage
  # Unadjusted by autosomal Sp:Bl6 SNP ratio
  expXiProp.chromLevel[[cellLine]] = chrXNi0 / chrXNi
  
  ###
  # adjusted expected Xi proportion
  chrXNi0.adj = sum(geneLevelASEdata.uniqued.adjRefCov[chrIndices[[whichChrs]], cellLine], na.rm=TRUE)
  chrXNi0.adj == chrXNi0
  # TRUE
  chrXNi1.adj = sum(geneLevelASEdata.uniqued.adjRefCov[chrIndices[[whichChrs]], cellLine], na.rm=TRUE)
  chrXNi.adj = sum(geneLevelASEdata.uniqued.adjTotCov[chrIndices[[whichChrs]], cellLine], na.rm=TRUE)
  chrXNi.adj ==  (chrXNi0.adj + chrXNi1.adj)
  # TRUE
  
  # Expected Xi:Xa ratio based on chrX Sp:Bl6 coverage ratio from total (chromosome-wide) SNP coverage
  # Adjusted by autosomal Sp:Bl6 SNP ratio 
  expXiProp.adj.chromLevel[[cellLine]] = chrXNi0.adj / chrXNi.adj
  
  
  ######
  # gene level parameters *** USE THIS ***
  
  # Expected Xi:Xa ratio based on Sp:Bl6 SNP coverage per chrX gene
  # Adjusted by autosomal Sp:Bl6 SNP coverage ratio 
  expXiProp.adj.geneLevel[[cellLine]] = mean(geneLevelASEdata.uniqued.adjRefCov[chrIndices[[whichChrs]], cellLine] / 
                                             geneLevelASEdata.uniqued.adjTotCov[chrIndices[[whichChrs]], cellLine], na.rm=TRUE)
  # Can get variance this way too.
  sdXiXa.adj.geneLevel[[cellLine]] = sd(geneLevelASEdata.uniqued.adjRefCov[chrIndices[[whichChrs]], cellLine] / 
                                          geneLevelASEdata.uniqued.adjTotCov[chrIndices[[whichChrs]], cellLine], na.rm=TRUE)
  
  # Expected & Q5 Xi coverage
  expXiCov.adj.geneLevel[[cellLine]] = mean(geneLevelASEdata.uniqued.adjRefCov[chrIndices[[whichChrs]], cellLine], na.rm=TRUE)
  quantsXiCov.adj.geneLevel[[cellLine]] = quantile(geneLevelASEdata.uniqued.adjRefCov[chrIndices[[whichChrs]], cellLine], 
                                                   probs=seq(0.5, 1, 0.05),
                                                   na.rm=TRUE)
  
  # Expected & Q5 Xi TPM
  expXiTPM.adj.geneLevel[[cellLine]] = mean(expressionTables[["refTPM"]][chrIndices[[whichChrs]], cellLine], na.rm=TRUE)
  quantsXiTPM.adj.geneLevel[[cellLine]] = quantile(expressionTables[["refTPM"]][chrIndices[[whichChrs]], cellLine],
                                                   probs=seq(0.5, 1, 0.05),
                                                   na.rm=TRUE)
  
  
  ######
  # SNP level parameters
  
  cellLineString = paste(cellLine, "_refCov", sep="")
  chrXSNPlevelRefCov = as.numeric(as.character(SNPlevelASEdata[SNPlevelASEdata$chr == whichChrs, cellLineString]))
  cellLineString = paste(cellLine, "_altCov", sep="")
  chrXSNPlevelAltCov = as.numeric(as.character(SNPlevelASEdata[SNPlevelASEdata$chr == whichChrs, cellLineString]))
  
  # Expected Xi:Xa ratio based on Sp:Bl6 SNP coverage ratio per SNP
  # Adjusted by autosomal Sp:Bl6 SNP ratio 
  expXiProp.adj.SNPlevel[[cellLine]] = mean(chrXSNPlevelRefCov / (chrXSNPlevelAltCov + chrXSNPlevelRefCov), na.rm=TRUE)
  
  # Can get variance this way too.
  sdXiXa.adj.SNPlevel[[cellLine]] = sd(chrXSNPlevelRefCov / (chrXSNPlevelAltCov + chrXSNPlevelRefCov), na.rm=TRUE)
}

#expXiProp.chromLevel
expXiProp.adj.chromLevel

expXiProp.adj.geneLevel
#sdXiXa.adj.geneLevel
#expXiCov.adj.geneLevel
#quantsXiCov.adj.geneLevel
#expXiTPM.adj.geneLevel
#quantsXiTPM.adj.geneLevel

expXiProp.adj.SNPlevel
#sdXiXa.adj.SNPlevel





#############
# Xi:Xa Sig test based on expected values at different 'parameters levels': chrom, gene, and SNP
# *** Gene-level expected value probably best, most relavent method ***

# alpha
CIalpha = 0.01


####
# For each parameter level for binomial test?
escapePvalueTable = list()
escapeLCILtable = c() # 1% Lower Condidence Interval Limit
escapeUCILtable = c() # 1% Upper Condidence Interval Limit
for (parameterLevel in c("chrom", "gene", "SNP")) {
  escapePvalueTable[[parameterLevel]] = c()
  for(cellLine in cellLines) {
    cat( parameterLevel, ";", cellLine, "\n", sep="")

    expXiProp2Use = ifelse(parameterLevel == "SNP", expXiProp.adj.SNPlevel[[cellLine]], 
                           ifelse(parameterLevel == "gene", expXiProp.adj.geneLevel[[cellLine]], 
                                  expXiProp.adj.chromLevel[[cellLine]]))

    escapePvalueVector = c()
    escapeLCILvector = c() # 1% Lower Condidence Interval Limit
    escapeUCILvector = c() # 1% Upper Condidence Interval Limit
    for(theChrXidx in chrIndices[[whichChrs]]) {
      # theChrXidx = 24298
      geneNi0 = round(geneLevelASEdata.uniqued.adjRefCov[theChrXidx, cellLine])
      geneNi = round(geneLevelASEdata.uniqued.adjTotCov[theChrXidx, cellLine])
      if (!is.numeric(geneNi0) || is.na(geneNi0) || geneNi0==0) {
        escapePvalueVector = c(escapePvalueVector, 1)
        escapeLCILvector = c(escapeLCILvector, 0)
        escapeUCILvector = c(escapeUCILvector, 0)
      } else {
        BTret = binom.test(x=geneNi0, 
                           n=geneNi, 
                           p=expXiProp2Use, 
                           conf.level=1-CIalpha,
                           #alternative="greater")
                           # 170425
                           alternative="two.sided")
        escapePvalueVector = c(escapePvalueVector, BTret$p.value)
        escapeLCILvector = c(escapeLCILvector, BTret$conf.int[1])
        escapeUCILvector = c(escapeUCILvector, BTret$conf.int[2])
        }
    }
    escapePvalueTable[[parameterLevel]] = cbind(escapePvalueTable[[parameterLevel]] , escapePvalueVector)
    escapeLCILtable[[parameterLevel]] = cbind(escapeLCILtable[[parameterLevel]] , escapeLCILvector)
    escapeUCILtable[[parameterLevel]] = cbind(escapeUCILtable[[parameterLevel]] , escapeUCILvector)
  }
  colnames(escapePvalueTable[[parameterLevel]]) = colnames(escapeLCILtable[[parameterLevel]]) = colnames(escapeUCILtable[[parameterLevel]]) = cellLines
  rownames(escapePvalueTable[[parameterLevel]]) =  rownames(escapeLCILtable[[parameterLevel]]) = rownames(escapeUCILtable[[parameterLevel]]) = rownames(geneLevelASEdata.uniqued.adjRefCov)[chrIndices[[whichChrs]]]
  #head(escapePvalueTable[[parameterLevel]])
}

all(names(escapePvalueTable[[parameterLevel]][, cellLine]) == names(geneLevelASEdata.uniqued.adjRefProp[chrIndices[[whichChrs]], cellLine]))
# TRUE









#############
## Escape genes

# Gene-level expected value probably best
parameterLevel = "gene"
#parameterLevel = "chrom"
#parameterLevel = "SNP"

# cellLinesToTest = "PatskiWT"
# cellLinesToTest = cellLines
cellLinesToTest = c( "PatskiWT", 
                     grep(paste("PatskiWT_0uMaza.*", sep=""), cellLines, value=TRUE),
                     "PatskiWTa", "PatskiWTb" )
cellLinesToTest = c( cellLinesToTest,
                     "PatskiDel1", 
                     grep(paste("PatskiDel1_0uMaza.*", sep=""), cellLines, value=TRUE),
                     grep(paste("PatskiInvDxz4.*", sep=""), cellLines, value=TRUE))
cellLinesToTest

####
# The proportion of genes any Xi expression that meet confidence limit
# Some genes not assigned RPKM by HTSeq: exclude these
LCIpropTable = c()
for(cellLine in cellLines) {
  HTseqExpression = expressionTables[["RPKM"]][chrIndices[[whichChrs]],cellLine] > 0 # Some genes not assigned RPKM by HTSeq: exclude these
  for (parameterLevel in c("chrom", "gene", "SNP")) {
    LCIescape = sum( HTseqExpression & escapeLCILtable[[parameterLevel]][,cellLine] > 0, na.rm=TRUE ) 
    XiExpression = sum( HTseqExpression & expressionTables[["refRPKM"]][chrIndices[[whichChrs]], cellLine] > 0, na.rm=TRUE)
    OverlapPct = sum( HTseqExpression & 
                        escapeLCILtable[[parameterLevel]][,cellLine] > 0 & 
                        expressionTables[["refTPM"]][chrIndices[[whichChrs]], cellLine] > 0, na.rm=TRUE ) /
      sum(HTseqExpression & expressionTables[["refTPM"]][chrIndices[[whichChrs]], cellLine] > 0, na.rm=TRUE)  * 100
    LCIpropTable = rbind(LCIpropTable, c( cellLine, parameterLevel, LCIescape, XiExpression, OverlapPct ))
  }
}
colnames(LCIpropTable) = c("cellLine", "parameterLevel", "#genesLCI>0", "#genesXi>0", "overlap%")
LCIpropTable
LCIpropTable[LCIpropTable[,"parameterLevel"] == "gene", c(1,3,4,5)]

fileName = file.path(escapeGeneDir, paste(workDir, "_overlapLCI0andXi0genes.tsv", sep=""))
write.table(LCIpropTable[LCIpropTable[,"parameterLevel"] == "gene", c(1,3,4,5)], fileName, row.names=FALSE, col.names=TRUE, sep="\t")  



####
# Just using the gene-level expected value,
# All genes with LCIescapeIdx > 0 
# This meets the main Berletch-Ma criterion:
# For each RNA-seq experiment, we called a gene "escape" if 
# (1) the 99% lower confidence limit (alpha = 0.01) of the escape probability was greater than zero, indicating significant contribution from the Xi,
for(cellLine in cellLinesToTest) {
  expXiProp2Use = ifelse(parameterLevel == "SNP", expXiProp.adj.SNPlevel[[cellLine]], 
                         ifelse(parameterLevel == "gene", expXiProp.adj.geneLevel[[cellLine]], 
                                expXiProp.adj.chromLevel[[cellLine]]))
  HTseqExpression = expressionTables[["RPKM"]][chrIndices[[whichChrs]],cellLine] > 0 # Some genes not assigned RPKM by HTSeq: exclude these
  LCIescapeIdx = HTseqExpression & escapeLCILtable[[parameterLevel]][,cellLine] > 0
  LCIescapeGeneTable = c()
  for(theChrXidx in chrIndices[[whichChrs]][LCIescapeIdx]) {
    geneNi0 = round(geneLevelASEdata.uniqued.adjRefCov[theChrXidx, cellLine])
    geneNi = round(geneLevelASEdata.uniqued.adjTotCov[theChrXidx, cellLine])
    if (!is.numeric(geneNi0) || is.na(geneNi0) || geneNi0==0) {
      escapePvalue = 1
      escapeLowerCIlimit = 0
      escapeUpperCIlimit = 0
    } else {
      BTret = binom.test(x=geneNi0, 
                         n=geneNi, 
                         p=expXiProp2Use, 
                         conf.level=1-CIalpha,
                         #alternative="greater")
                         # 170425
                         alternative="two.sided")
      escapePvalue = BTret$p.value
      escapeLowerCIlimit = BTret$conf.int[1]
      escapeUpperCIlimit = BTret$conf.int[2]
    }
    geneName = rownames(geneLevelASEdata.uniqued.adjRefCov)[theChrXidx]
    geneLength = geneLengths[theChrXidx]
    SNPcounts = geneLevelASEdata.uniqued.SNPcounts[theChrXidx]
    RPKM = expressionTables[["RPKM"]][theChrXidx, cellLine]
    TPM = expressionTables[["TPM"]][theChrXidx, cellLine]
    refRPKM = expressionTables[["refRPKM"]][theChrXidx, cellLine]
    TPM = expressionTables[["TPM"]][theChrXidx, cellLine]
    refTPM = expressionTables[["refTPM"]][theChrXidx, cellLine]
    LCIescapeGeneTable = rbind(LCIescapeGeneTable,
                               c( geneLength, SNPcounts, 
                                  RPKM, TPM, 
                                  refRPKM, refRPKM*10, 
                                  refTPM, refTPM*10,
                                  geneNi0, geneNi, signif( geneNi0/geneNi, 3), signif( expXiProp2Use, 3), 
                                  signif( escapePvalue, 3), signif( escapeLowerCIlimit, 3), signif( escapeUpperCIlimit, 3)))
  }
  rownames(LCIescapeGeneTable) = rownames(geneLevelASEdata.uniqued.adjRefCov)[ chrIndices[[whichChrs]][LCIescapeIdx] ]
  colnames(LCIescapeGeneTable) = c( "geneLength", "SNPcounts", "RPKM", "TPM", "XiRPKM", "Xi-SRPM", "XiTPM", "Xi-STPM", 
                                    "XiSNPcov", "Xi+XaSNPcov", "observedXiProp", "expectedXiProp",
                                    "P-value", "LowerCIlimit", "UpperCIlimit")
  
  #head(LCIescapeGeneTable[order(LCIescapeGeneTable[,"escapePvalue"]),])
  #tail(LCIescapeGeneTable[order(LCIescapeGeneTable[,"escapePvalue"]),])

  fileName = file.path(escapeGeneDir, paste(workDir, "_genesMeetingBinomialLowerCIlimit0_", 
                                            cellLine, ".updated170425.tsv", sep=""))
  write.table(LCIescapeGeneTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  
}  

# LCIescapeGeneTable


####
# compare published escape genes to my results
BMescapeIdx = match(escapeGenes, rownames(geneLevelASEdata.uniqued.adjRefCov))
for(cellLine in cellLinesToTest) {
  expXiProp2Use = ifelse(parameterLevel == "SNP", expXiProp.adj.SNPlevel[[cellLine]], 
                         ifelse(parameterLevel == "gene", expXiProp.adj.geneLevel[[cellLine]], 
                                expXiProp.adj.chromLevel[[cellLine]]))
  LCIescapeGeneTable = c()
  for(theEscapeeidx in BMescapeIdx) {
    geneNi0 = round(geneLevelASEdata.uniqued.adjRefCov[theEscapeeidx, cellLine])
    geneNi = round(geneLevelASEdata.uniqued.adjTotCov[theEscapeeidx, cellLine])
    if (!is.numeric(geneNi0) || is.na(geneNi0) || geneNi0==0) {
      escapePvalue = 1
      escapeLowerCIlimit = 0
      escapeUpperCIlimit = 0
    } else {
      BTret = binom.test(x=geneNi0, 
                         n=geneNi, 
                         p=expXiProp2Use, 
                         conf.level=1-CIalpha,
                         #alternative="greater")
                         # 170425
                         alternative="two.sided")
      escapePvalue = BTret$p.value
      escapeLowerCIlimit = BTret$conf.int[1]
      escapeUpperCIlimit = BTret$conf.int[2]
    }
    geneName = rownames(geneLevelASEdata.uniqued.adjRefCov)[theEscapeeidx]
    geneLength = geneLengths[theEscapeeidx]
    SNPcounts = geneLevelASEdata.uniqued.SNPcounts[theEscapeeidx]
    RPKM = expressionTables[["RPKM"]][theEscapeeidx, cellLine]
    TPM = expressionTables[["TPM"]][theChrXidx, cellLine]
    refRPKM = expressionTables[["refRPKM"]][theEscapeeidx, cellLine]
    TPM = expressionTables[["TPM"]][theEscapeeidx, cellLine]
    refTPM = expressionTables[["refTPM"]][theEscapeeidx, cellLine]
    LCIescapeGeneTable = rbind(LCIescapeGeneTable,
                               c( geneLength, SNPcounts, RPKM, TPM, refRPKM, refRPKM*10, refTPM, refTPM*10,
                                  geneNi0, geneNi, signif( geneNi0/geneNi, 3), signif( expXiProp2Use, 3), 
                                  signif( escapePvalue, 3), signif( escapeLowerCIlimit, 3), signif( escapeUpperCIlimit, 3)))
  }
  rownames(LCIescapeGeneTable) = rownames(geneLevelASEdata.uniqued.adjRefCov)[ BMescapeIdx ]
  colnames(LCIescapeGeneTable) = c( "geneLength", "SNPcounts", "RPKM", "TPM", "XiRPKM", "Xi-SRPM", "XiTPM", "Xi-STPM", 
                                    "XiSNPcov", "Xi+XaSNPcov", "observedXiProp", "expectedXiProp",
                                    "P-value", "LowerCIlimit", "UpperCIlimit")
  
  #head(LCIescapeGeneTable[order(LCIescapeGeneTable[,"escapePvalue"]),])
  #tail(LCIescapeGeneTable[order(LCIescapeGeneTable[,"escapePvalue"]),])
  
  fileName = file.path(escapeGeneDir, paste(workDir, "_BerletchMaPatskiEscapeGenes_Vs_", 
                                            cellLine, ".updated170425.tsv", sep=""))
  write.table(LCIescapeGeneTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  
}  

# LCIescapeGeneTable


#############
## Escape gene lists

# incDel=FALSE
incDel=TRUE

# incBM=TRUE
incBM=FALSE

incInvDxz4=TRUE

cellLinesToTest = c( "PatskiWT", 
                     grep(paste("PatskiWT_0uMaza.*", sep=""), cellLines, value=TRUE),
                     "PatskiWTa", "PatskiWTb")
if (incDel) {
  cellLinesToTest = c( cellLinesToTest,
                       "PatskiDel1", 
                       grep(paste("PatskiDel1_0uMaza.*", sep=""), cellLines, value=TRUE) )
}
if (incInvDxz4) {
  cellLinesToTest = c( cellLinesToTest,
                       "PatskiInvDxz4a", "PatskiInvDxz4b" )
}

parameterLevel = "gene"

genesToExclude=c('Gm6568')

escapeGeneLists = list()

# For each RNA-seq experiment, we called a gene "escape" if 
# (1) the 99% lower confidence limit (alpha = 0.01) of the escape probability was greater than 0.01 (rather than 0) indicating some contribution from the Xi,
# (2) the diploid gene expression measured by TPM >=1, indicating that the gene was expressed, 
# (3) the Xi-TPM was >= 0.1, representing sufficient reads from the Xi, and
# (4) SNP coverage >= 5
  
SNPcovThresh = 5
minLCI = 0.01 # 1% 
TPMthresh = 1
XiTPMthresh = 0.1
  
sigLabel = paste(
    "minCov", SNPcovThresh, 
    "_minTPKM", TPMthresh,
    "_minXiTPMthresh", XiTPMthresh,
    "_minLCI", minLCI,
    sep="")
  
for(cellLine in cellLinesToTest) {
    theLabel = paste(whichChrs, "_", cellLine, "_", sigLabel, sep="")
    escapeGeneLists[[theLabel]] =
      names(which(
        geneLevelASEdata.uniqued.adjTotCov[chrIndices[[whichChrs]], cellLine] >= SNPcovThresh &
          expressionTables[["TPM"]][chrIndices[[whichChrs]], cellLine] >= TPMthresh &
          # 171202
          # expressionTables[["refTPM"]][chrIndices[[whichChrs]], cellLine]*10 >= XiTPMthresh &
          expressionTables[["refTPM"]][chrIndices[[whichChrs]], cellLine] >= XiTPMthresh &
          escapeLCILtable[[parameterLevel]][,cellLine] > minLCI ))
}

if (incBM) {
  escapeGeneLists[["BerletchMa.EscapeGenes"]] = escapeGenes
}

# Filter troublesome genes
for (dataSet in names(escapeGeneLists)) {
  escapeGeneLists[[dataSet]] = setdiff(escapeGeneLists[[dataSet]], genesToExclude)
  cat(dataSet, "\t", length(escapeGeneLists[[dataSet]]), "\t", 
      length(intersect(escapeGenes, escapeGeneLists[[dataSet]])) / length(escapeGenes), "\t",
      length(intersect(escapeGenes, escapeGeneLists[[dataSet]])) / length(escapeGeneLists[[dataSet]]), "\n", sep="")
}
#length(escapeGenes)

selectSubDir = ""

escapeGeneSubDir = file.path(escapeGeneDir, paste( ifelse(incDel, "_withDel1", ""),
                                                   ifelse(incInvDxz4, "_withInvDxz4", ""),
                                                   ifelse(!incBM, "_sansBM", ""), sep=""))
if ( ! file.exists( escapeGeneSubDir ) ) {
  if ( ! dir.create( escapeGeneSubDir, showWarnings=FALSE, recursive=TRUE )  ) escapeGeneSubDir = escapeGeneDir
}

###############
# All data sets
selectList = names(escapeGeneLists)
altNameList = sub("chrX_Patski", "", sub(paste("_", sigLabel, sep=""), "", selectList))

# selectList = grep("ENCODE", selectList, invert=TRUE, value=TRUE)
# altNameList = grep("ENCODE", altNameList, invert=TRUE, value=TRUE)

###############
# Overlap table

selectEscapeGeneLists = list()
for (selectID in selectList) {
  selectEscapeGeneLists[[selectID]] = escapeGeneLists[[selectID]]
}
names(selectEscapeGeneLists) = altNameList

geneOverlapTable = c()
for (listA in names(selectEscapeGeneLists)) {
  for (listB in names(selectEscapeGeneLists)) {
    geneOverlapTable = c(geneOverlapTable, length(intersect(selectEscapeGeneLists[[listA]], selectEscapeGeneLists[[listB]])) / length(selectEscapeGeneLists[[listA]]))
  }
}
geneOverlapTable = matrix(geneOverlapTable, nrow=length(selectEscapeGeneLists), ncol=length(selectEscapeGeneLists))
colnames(geneOverlapTable) = names(selectEscapeGeneLists)
rownames(geneOverlapTable) = names(selectEscapeGeneLists)

fileName = file.path(escapeGeneSubDir, paste(workDir, "_geneLists_overlapTable_", 
                                        sigLabel, ".tsv", sep=""))
write.table(geneOverlapTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  


#############
# Escape Genes Matrix
escapeGenesMatrix = c()
allEscapeGenes = sort(unique(unlist(selectEscapeGeneLists))) # Sort by gene names first
for (geneListName in names(selectEscapeGeneLists)) {
  escapeGenesMatrix = cbind( escapeGenesMatrix, as.numeric(!is.na(match(allEscapeGenes, selectEscapeGeneLists[[geneListName]]))) )
}
#colnames(escapeGenesMatrix) = c("WT", "Del1", "Del2", "Del5", "BerletchMa")
colnames(escapeGenesMatrix) = altNameList
rownames(escapeGenesMatrix) = allEscapeGenes

overlapCount = rowSums(escapeGenesMatrix)
escapeGenesMatrix = cbind( escapeGenesMatrix, overlapCount )
# sortIdx = order(overlapCount, 
#                 escapeGenesMatrix[,1], 
#                 escapeGenesMatrix[,2], 
#                 escapeGenesMatrix[,3], 
#                 escapeGenesMatrix[,4], 
#                 escapeGenesMatrix[,5], 
#                 decreasing=TRUE)
orderText = "order(overlapCount"
for (colID in 1:ncol(escapeGenesMatrix)) {
  orderText = paste( orderText, ", escapeGenesMatrix[,", colID, "]", sep="")
}
orderText = paste( orderText, ", decreasing=TRUE)", sep="")
sortIdx = eval(parse(text=orderText))
escapeGenesMatrix = escapeGenesMatrix[sortIdx,]

fileName = file.path(escapeGeneSubDir, paste(workDir, "_geneLists_overlapMatrix_", 
                                        sigLabel, ".tsv", sep=""))
write.table(escapeGenesMatrix, fileName, row.names=TRUE, col.names=NA, sep="\t")  


# Other data sources but not sorted:
# geneLevelASEdata
# SNPlevelASEdata

geneData = c()
geneDataCols = c()
geneNames = rownames(escapeGenesMatrix)
geneDataNameIdx = match(geneNames, rownames(geneLevelASEdata.uniqued.adjRefCov[chrIndices[[whichChrs]], ]))
for (cellLine in cellLines) {
  geneData = cbind( geneData, round(geneLevelASEdata.uniqued.adjRefProp[chrIndices[[whichChrs]], cellLine][geneDataNameIdx], 3) )
  geneData = cbind( geneData, round(geneLevelASEdata.uniqued.adjRefCov[chrIndices[[whichChrs]], cellLine][geneDataNameIdx]) )
  geneData = cbind( geneData, geneLevelASEdata.uniqued.adjAltCov[chrIndices[[whichChrs]], cellLine][geneDataNameIdx] )
  geneData = cbind( geneData, round(expressionTables[["RPKM"]][chrIndices[[whichChrs]], cellLine][geneDataNameIdx], 2) )
  geneData = cbind( geneData, round(expressionTables[["refRPKM"]][chrIndices[[whichChrs]], cellLine][geneDataNameIdx], 2) )
  geneData = cbind( geneData, round(expressionTables[["altRPKM"]][chrIndices[[whichChrs]], cellLine][geneDataNameIdx], 2) )
  geneData = cbind( geneData, round(expressionTables[["TPM"]][chrIndices[[whichChrs]], cellLine][geneDataNameIdx], 2) )
  geneData = cbind( geneData, round(expressionTables[["refTPM"]][chrIndices[[whichChrs]], cellLine][geneDataNameIdx], 2) )
  geneData = cbind( geneData, round(expressionTables[["altTPM"]][chrIndices[[whichChrs]], cellLine][geneDataNameIdx], 2) )
  # geneData = cbind( geneData, signif(escapePvalueTable[[parameterLevel]][,cellLine][geneDataNameIdx], 3) )
  geneData = cbind( geneData, signif(escapeLCILtable[[parameterLevel]][,cellLine][geneDataNameIdx], 3) )
  geneData = cbind( geneData, signif(escapeUCILtable[[parameterLevel]][,cellLine][geneDataNameIdx], 3) )
  geneDataCols = c( geneDataCols,
                    paste( sub("Patski", "", cellLine), c("Xi/(Xi+Xa) SNPcov", "Xi SNPcov", "Xa SNPcov", 
                                                          "RPKM", "XiRPKM", "XaRPKM", 
                                                          "TPM", "XiTPM", "XaTPM",
                                                          #"P-value", 
                                                          "CIlowerLimit", "CIupperLimit"),
                                                          sep=" ") )
}
colnames(geneData) = geneDataCols

escapeGenesMatrixAnnotated = cbind(escapeGenesMatrix, geneData)

fileName = file.path(escapeGeneSubDir, paste(workDir, "_geneLists_overlapMatrix_annotated_", 
                                        sigLabel, ".tsv", sep=""))
write.table(escapeGenesMatrixAnnotated, fileName, row.names=TRUE, col.names=NA, sep="\t")  


#############
# WT escapees

WTescapeGenesMatrix = escapeGenesMatrix[,sub("Patski", "", cellLinesToTest)]
colnames(WTescapeGenesMatrix) = sub("Patski", "", cellLinesToTest)

WTescapeGenesMatrix = WTescapeGenesMatrix[rowSums(WTescapeGenesMatrix)>=2,]

# Drop Gm6568
WTescapeGenesMatrix[!rownames(WTescapeGenesMatrix) == "Gm6568",]
nrow(WTescapeGenesMatrix)
# 29

WTescapeeIdx = match(rownames(WTescapeGenesMatrix), geneLevelASEdata.uniqued$geneID)
WTescapees = geneLevelASEdata.uniqued[WTescapeeIdx, c(4,1,2,3,6)]
WTescapees = WTescapees[order(WTescapees$start),]

fileName = file.path(escapeGeneSubDir, paste(workDir, "_WTescapeGenes_", 
                                             sigLabel, ".tsv", sep=""))
write.table(WTescapees, fileName, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)  


#############
# Another way to list genes

# # for every element in the list (`myl`), write that element to a file
# # and append if necessary. also, if list element is a character, write
# # to as many columns as there are characters.
# cat(paste(names(selectEscapeGeneLists), collapse="\t"), "\n", sep="", file=fileName)
# lapply(X = selectEscapeGeneLists, FUN = function(x) { write(x, append = T, file=fileName, ncolumns = length(x), sep="\t")})

maxLen = 0
for (listA in names(selectEscapeGeneLists)) {
  maxLen = ifelse(length(selectEscapeGeneLists[[listA]])>maxLen, length(selectEscapeGeneLists[[listA]]), maxLen)
}

geneListTable = matrix("", nrow=maxLen, ncol=length(selectEscapeGeneLists))
colnames(geneListTable) = names(selectEscapeGeneLists)
for (listA in names(selectEscapeGeneLists)) {
  geneListTable[1:length(selectEscapeGeneLists[[listA]]),listA] = sort(selectEscapeGeneLists[[listA]])
}

fileName = file.path(escapeGeneSubDir, paste(workDir, "_geneLists_", 
                                        sigLabel,".tsv", sep=""))
write.table(geneListTable, fileName, row.names=TRUE, col.names=NA, sep="\t")  


#############
# UpSetR
require(ggplot2); require(plyr); require(gridExtra); require(grid);
require(UpSetR)

# For some reason must reload the data for upsetR to work
# otherwise I get error:
#     Error in colSums(data[sets]) : 
#     'x' must be an array of at least two dimensions
fileName = file.path(escapeGeneSubDir, paste(workDir, "_geneLists_overlapMatrix_", 
                                        sigLabel, ".tsv", sep=""))
UpSetRdata = read.table(file=fileName, header=TRUE, sep="\t")
rownames(UpSetRdata) = UpSetRdata[,1]
UpSetRdata = UpSetRdata[,c(-1,-ncol(UpSetRdata))]


fileName = file.path(escapeGeneSubDir, paste(workDir, "_geneLists_UpSetRed_", 
                                        sigLabel, ".pdf", sep=""))
pdf(fileName, width=6, height=6)
upset(UpSetRdata, sets = colnames(UpSetRdata), sets.bar.color = "#56B4E9", order.by = "freq")
dev.off()

# Only if in two samples
UpSetRdata.filtered = UpSetRdata[rowSums(UpSetRdata, na.rm=TRUE)>1,]
fileName = file.path(escapeGeneSubDir, paste(workDir, "_geneLists_UpSetRed_", 
                                             sigLabel, ".filtered.pdf", sep=""))
pdf(fileName, width=6, height=6)
upset(UpSetRdata.filtered, sets = colnames(UpSetRdata.filtered), sets.bar.color = "#56B4E9", order.by = "freq")
dev.off()
