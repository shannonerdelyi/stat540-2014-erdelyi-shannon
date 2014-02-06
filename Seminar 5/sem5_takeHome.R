###############################################################################
# STAT 540: Seminar 5
# Low Volume Linear Modeling 
# Take home problem:  Use data aggregation techniques to repeat the 
#                     analysis for more than one gene
# Date: 05-02-2014
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
n <- 10
(keepGenes <- rownames(dat)[sample(1:nrow(dat), n)])

###############################################################################
# Functions
###############################################################################

# Function to return a data frame with only the specified genes
makeDF <- function(geneNames){
  miniDat <- dat[geneNames, ]
  miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),
                        gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
                                      levels = rownames(miniDat)))
  miniDat <- suppressWarnings(data.frame(des, miniDat))
  miniDat$devStage <- factor(miniDat$devStage, 
                             levels(miniDat$devStage)[c(2, 4, 5, 3, 1)])
  return(miniDat)
}

# Function to make a stripplot using the photoRec data
makeSP <- function(data, ...){
  stripplot(gExp ~ devStage | gene, data,
            group = gType, jitter.data = TRUE,
            scales=list(tck=c(1, 0), rot=c(45, 0)),
            auto.key = TRUE, 
            type = c('p', 'a'), 
            par.settings=simpleTheme(col=c("coral2","cadetblue3"), pch=16),
            grid = TRUE, ...)
}

# Function to perform a t-test comparing P2 to 4_week devStages for one gene
doTwoSampleTests <- function(miniDf) {
  tt <- t.test(gExp ~ devStage, 
               subset(miniDf, devStage %in% c("P2", "4_weeks")),
               var.equal=T)
  results <- with(tt, c(estimate, statistic, parameter, pval=p.value))
  return(results)
}

# Function to fit a linear model with devStage (for wt's only)
fitOneWayANOVA <- function(miniDf) {
  lm(gExp ~ devStage, miniDf, subset=gType=="wt")
}

# Function to fit a linear model with devStage and gType
# Test for importance of interaction term
interactionPval <- function(miniDf) {
  fitBig <- lm(gExp ~ gType*devStage, miniDf)
  fitSmall <- lm(gExp ~ gType + devStage, miniDf)
  c(pval=anova(fitSmall, fitBig)[2,6])
}

###############################################################################
# Analysis
###############################################################################

# make a data frame
sDat <- makeDF(keepGenes)
str(sDat)
head(sDat)

# make a stripplot
makeSP(sDat)

# two sample t-tests (P2 vs. 4_weeks)
ddply(sDat, ~ gene, doTwoSampleTests)

# linear model with categorical covariate 
ow <- dlply(sDat, ~ gene, fitOneWayANOVA)
(owCoefs <- ldply(ow, coef))

# inference for a contrast (P2 vs P10)
cont <- matrix(c(0, 1, 0, -1, 0), nrow=1)
(inf <- ldply(ow, function(x){
  est <- cont %*% coef(x)
  se <- cont %*% vcov(x) %*% t(cont)
  t <- est/se
  p <- 2*pt(abs(t), df = df.residual(x), lower.tail = FALSE)
  c(est=est, se=se, pval=p)
  }))

# importance of interaction term
ddply(sDat, ~ gene, interactionPval)
