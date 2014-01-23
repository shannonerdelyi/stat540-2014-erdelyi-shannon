###############################################################################
# STAT 540: Seminar 3
# Take-home problem: look for gene expression differences between genotypes
# Date: 22-01-2014
# Author: Shannon Erdelyi
###############################################################################

# librarys
library(lattice)
library(RColorBrewer)

# data
des <- read.table("GSE4051_design.tsv", header=T)
dat <- read.table("GSE4051_data.tsv", header=T)

# constants
n <- 20  ## gene samples
myCols <- colorRampPalette(brewer.pal(n = 9, "BuPu"))  ## heatmap scale
ttResults <- data.frame(gene = NA, pval = NA)  ## empty df for loop

###############################################################################
# Data Preparation
###############################################################################

# select random subset of probesets/genes
set.seed(8)
samp <- sample(1:nrow(dat), n)
small <- dat[samp, ]

# create long dataframe with experimental design
long <- reshape(small, 
                v.names = "geneExp", varying = list(1:ncol(small)),
                timevar = "sample", times = names(small),
                idvar = "gene", ids = row.names(small),
                direction="long")
long$gene <- as.factor(long$gene)
long$sample <- as.factor(long$sample)
myDat <- merge(long, des, by.x="sample", by.y="sidChar")
str(myDat)

# create matrix format for heatmap plotting
smallM <- as.matrix(t(small))
rownames(smallM) <- with(des, paste(gType, sidChar, sep="_"))

###############################################################################
# Exploratory Analysis
###############################################################################

# order heatmap by genotypes and look for differences
heatmap(smallM[order(rownames(smallM)), order(apply(smallM, 2, mean))], 
        Rowv = NA, Colv = NA, margins=c(6, 5),
        col = myCols(500))
# 1429021_at and 1457045_at are jumping out

# compare gene expression distributions for genotypes
densityplot( ~ geneExp | gene, myDat,
           groups = gType, auto.key = T)
# the means look different for these genes
#   1457045_at
#   1429021_at 
#   1437750_at

# lets do some t tests to find out which genes differ
# this loop is ugly ... there has to be a better way
pairs <- with(myDat, tapply(geneExp, list(gene, gType), function(x) x))
for(i in 1:nrow(pairs)){
  test <- t.test(pairs[i, "NrlKO"][[1]], pairs[i, "wt"][[1]])
  ttResults <- rbind(ttResults, c(gene = rownames(pairs)[i], 
                                  pval = test$p.value))
}
(ttResults <- ttResults[-1, ])
(diffGenes <- ttResults$gene[as.numeric(ttResults$pval) <= 0.05])

# we can look at the data for only those genes with pval <= 0.05
stripplot(geneExp ~ gType | gene, myDat, subset = gene %in% diffGenes,
          group = gType, jitter.data = T,
          scales = list(y = list(relation = "free")))
