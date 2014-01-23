###############################################################################
# STAT 540: Seminar 3
# Graphics: Lattice
# Date: 22-01-2014
# Author: Shannon Erdelyi
###############################################################################

library(lattice)

# data
mini <- read.table("GSE4051_MINI.txt", header=T, row.names=1)
design <- read.table("GSE4051_design.tsv", header=T)
fullDat <- read.table("GSE4051_data.tsv", header=T)

###############################################################################
# Scatterplots
###############################################################################

dat2 <- with(mini,
             data.frame(sidChar = rownames(mini), 
                        sidNum = sample, 
                        devStage, gType, crabHammer,
                        probeset = factor(rep(c("eggBomb", "poisonFang"), 
                                              each = nrow(mini))),
                        geneExp = c(eggBomb, poisonFang)))

xyplot(geneExp ~ crabHammer | probeset, dat2,
       grid = TRUE,
       groups = devStage, 
       auto.key = list(space = "right"),
       scales = list(alternating = F, tck=c(1, 0)), 
       pch = 16)

###############################################################################
# Density plots
###############################################################################

dat3 <- with(mini,
             data.frame(sidChar = rownames(mini), 
                        sidNum = sample, 
                        devStage, gType,
                        probeset = factor(rep(c("crabHammer", 
                                                "eggBomb","poisonFang"), 
                                              each = nrow(mini))),
                        geneExp = c(crabHammer, eggBomb, poisonFang)))

densityplot( ~ geneExp | gType, dat3,
             groups = devStage,
             auto.key = list(space = "right"),
             scales = list(alternating = F, tck=c(1, 0)),
             n = 500)

