###############################################################################
# STAT 540: PhotoRec Data
# Read data and design
# Function to select a small subset 
# Date: 05-02-2014
# Author: Shannon Erdelyi
###############################################################################

# data
dir <- "/Users/shannonerdelyi/Dropbox/UBC/W2014/STAT 540/Data"
des <- read.table(paste(dir, "/GSE4051_design.tsv", sep=""), header=T)
dat <- read.table(paste(dir, "/GSE4051_data.tsv", sep=""), header=T)

# put devStage in chrono order
des$devStage <- factor(des$devStage, levels(des$devStage)[c(2, 4, 5, 3, 1)])

###############################################################################
# Functions
###############################################################################

# Function to return a data frame with only the specified genes
miniDF <- function(geneNames){
  miniDat <- dat[geneNames, ]
  miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),
                        gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
                                      levels = rownames(miniDat)))
  miniDat <- suppressWarnings(data.frame(des, miniDat))
  return(miniDat)
}

# Function to make a stripplot using the photoRec data
makeSP <- function(data, ...){
  stripplot(gExp ~ devStage | gene, data,
            group = gType, jitter.data = TRUE,
            scales=list(tck=c(1, 0), rot=c(45, 0)),
            type = c('p', 'a'), 
            par.settings=simpleTheme(col=c("coral2","cadetblue3"), pch=16),
            grid = TRUE, ...)
}