###############################################################################
# STAT 540: Seminar 4
# Take home problem: two sample tests with 100 genes 
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
n <- 100
keepGenes <- sample(1:nrow(dat), n)

###############################################################################
# Data Preparation
###############################################################################

miniDat <- dat[keepGenes, ]
miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),
                      gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
                                    levels = rownames(miniDat)))
miniDat <- suppressWarnings(data.frame(des, miniDat))
str(miniDat)

###############################################################################
# Two sample tests
###############################################################################

# get the p-values 
pvals <- suppressWarnings(ddply(miniDat, ~ gene, function(miniDf) {
  tt <- t.test(gExp ~ gType, miniDf)$p.value
  wt <- wilcox.test(gExp ~ gType, miniDf)$p.value
  kt <- with(miniDf, 
             ks.test(gExp[gType == "NrlKO"], gExp[gType == "wt"]))$p.value
  c(t = tt, wilcoxon = wt, ks = kt)
}))

head(pvals)

# compare pvalues across tests
plot(pvals[,-1])
plot(log(pvals[,-1]))

# apply 0.05 significance threshold
sig <- pvals[,-1] <= 0.05

# test agreement by gene
# 0 --> no tests were significant
# 3 --> all three tests were significant 
table(apply(sig, 1, sum))

# how many genes are significant
apply(sig, 2, sum)