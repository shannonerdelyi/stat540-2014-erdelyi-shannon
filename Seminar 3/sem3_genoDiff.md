# Differences in gene expression by genotype
### STAT 540: Seminar 3
### Shannon Erdelyi


```r
# librarys
library(lattice)
library(RColorBrewer)

# data
des <- read.table("GSE4051_design.tsv", header = T)
dat <- read.table("GSE4051_data.tsv", header = T)

# constants
n <- 20  ## gene samples
myCols <- colorRampPalette(brewer.pal(n = 9, "BuPu"))  ## heatmap scale
ttResults <- data.frame(gene = NA, pval = NA)  ## empty df for loop

###############################################################################
############################################################################### Data
############################################################################### Preparation

# select random subset of probesets/genes
set.seed(8)
samp <- sample(1:nrow(dat), n)
small <- dat[samp, ]

# create long dataframe with experimental design
long <- reshape(small, v.names = "geneExp", varying = list(1:ncol(small)), timevar = "sample", 
    times = names(small), idvar = "gene", ids = row.names(small), direction = "long")
long$gene <- as.factor(long$gene)
long$sample <- as.factor(long$sample)
myDat <- merge(long, des, by.x = "sample", by.y = "sidChar")
str(myDat)
```

```
## 'data.frame':	780 obs. of  6 variables:
##  $ sample  : Factor w/ 39 levels "Sample_1","Sample_10",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ geneExp : num  10.54 10.9 11.38 8.26 5.75 ...
##  $ gene    : Factor w/ 20 levels "1415708_at","1419053_at",..: 11 4 18 15 8 16 7 20 17 14 ...
##  $ sidNum  : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ devStage: Factor w/ 5 levels "4_weeks","E16",..: 5 5 5 5 5 5 5 5 5 5 ...
##  $ gType   : Factor w/ 2 levels "NrlKO","wt": 1 1 1 1 1 1 1 1 1 1 ...
```

```r

# create matrix format for heatmap plotting
smallM <- as.matrix(t(small))
rownames(smallM) <- with(des, paste(gType, sidChar, sep = "_"))

###############################################################################
############################################################################### Exploratory
############################################################################### Analysis

# order heatmap by genotypes and look for differences
heatmap(smallM[order(rownames(smallM)), order(apply(smallM, 2, mean))], Rowv = NA, 
    Colv = NA, margins = c(6, 5), col = myCols(500))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-11.png) 

```r
# 1429021_at and 1457045_at are jumping out

# compare gene expression distributions for genotypes
densityplot(~geneExp | gene, myDat, groups = gType, auto.key = T)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-12.png) 

```r
# the means look different for these genes 1457045_at 1429021_at
# 1437750_at

# lets do some t tests to find out which genes differ this loop is ugly
# ... there has to be a better way
pairs <- with(myDat, tapply(geneExp, list(gene, gType), function(x) x))
for (i in 1:nrow(pairs)) {
    test <- t.test(pairs[i, "NrlKO"][[1]], pairs[i, "wt"][[1]])
    ttResults <- rbind(ttResults, c(gene = rownames(pairs)[i], pval = test$p.value))
}
(ttResults <- ttResults[-1, ])
```

```
##            gene                 pval
## 2    1415708_at    0.177755258547863
## 3    1419053_at    0.815330917561272
## 4  1421392_a_at    0.337458631383954
## 5  1424407_s_at    0.011441778417152
## 6    1426778_at   0.0469922530002784
## 7    1427243_at     0.46735158482114
## 8  1427929_a_at    0.781328560426702
## 9    1429021_at   0.0015549737403702
## 10   1434310_at   0.0758980812135523
## 11   1435194_at    0.912370715352243
## 12   1435552_at    0.842958777709302
## 13   1437750_at 0.000381169546848767
## 14   1438784_at    0.696057828277474
## 15   1443950_at    0.321378325751261
## 16   1444288_at    0.935752405224795
## 17   1448113_at    0.714489960072615
## 18   1450193_at    0.350054842160445
## 19   1451485_at    0.191994323268592
## 20   1457045_at 5.56888082228149e-05
## 21   1457272_at    0.519888324436242
```

```r
(diffGenes <- ttResults$gene[as.numeric(ttResults$pval) <= 0.05])
```

```
## [1] "1424407_s_at" "1426778_at"   "1429021_at"   "1437750_at"  
## [5] "1457045_at"
```

```r

# we can look at the data for only those genes with pval <= 0.05
stripplot(geneExp ~ gType | gene, myDat, subset = gene %in% diffGenes, group = gType, 
    jitter.data = T, scales = list(y = list(relation = "free")))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-13.png) 


