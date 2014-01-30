###############################################################################
# STAT 540: Seminar 4
# Take home problem: DE analysis times 
# Date: 29-01-2014
# Author: Shannon Erdelyi
###############################################################################

# librarys
library(lattice)
library(plyr)

# data
dir <- "/Users/shannonerdelyi/Dropbox/UBC/W2014/STAT 540/Data"
des <- read.table(paste(dir, "/GSE4051_design.tsv", sep=""), header=T)
dat <- read.table(paste(dir, "/GSE4051_data.tsv", sep=""), header=T)

# constants
set.seed(8)
n <- c(5, 10, 50, 100, 500)

###############################################################################
# Function to perform DE analysis for n genes
# input: n
# output: system time
###############################################################################

DErunTime <- function(n){
  # prepare data
  keepGenes <- sample(1:nrow(dat), n)
  miniDat <- dat[keepGenes, ]
  miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),
                        gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
                                      levels = rownames(miniDat)))
  miniDat <- suppressWarnings(data.frame(des, miniDat))
  
  # DE analysis time
  DEtime <- system.time(pvals <-  suppressWarnings(
    ddply(miniDat, ~ gene, function(miniDf) {
      tt <- t.test(gExp ~ gType, miniDf)$p.value
      wt <- wilcox.test(gExp ~ gType, miniDf)$p.value
      kt <- with(miniDf, ks.test(gExp[gType == "NrlKO"], 
                                 gExp[gType == "wt"]))$p.value
      c(t = tt, wilcoxon = wt, ks = kt)
    })))
  
 return(DEtime)
}

###############################################################################
# Simulation
###############################################################################

# run 10 simulations for each value of n
trials <- rep(n, each=10)
sim <- lapply(trials, DErunTime)
times <- data.frame(n=trials, do.call("rbind", sim))
head(times)

# plot the results
stripplot(elapsed ~ factor(n), times, 
          type = c("p", "a"), fun = mean,
          xlab = "Number of Genes",
          ylab = "Elapsed Time (s)")

