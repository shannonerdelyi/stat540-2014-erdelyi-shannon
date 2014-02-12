#####################################################################
# STAT 540: Seminar 6
# High Volume Linear Models
# Author: Shannon Erdelyi
# Date: 12-02-2014
#####################################################################

# libraries
library(limma)
library(lattice)
library(hexbin)

# source code with
#   read 'dat' and 'des'
#   'miniDF' function to subset data based on gene names
#   'makeSP' function to make a stripplot
dir <- "/Users/shannonerdelyi/Dropbox/UBC/W2014/STAT 540"
source(paste(dir, "/Seminars/photoRec_Functions.R", sep=""))


#####################################################################
# Using limma on the wild types
#####################################################################

# select wild types only
wtDes <- subset(des, gType == "wt")
wtDat <- subset(dat, select = des$gType == "wt")

# linear models with limma
wtDesMat <- model.matrix(~devStage, wtDes)
wtFit <- lmFit(wtDat, wtDesMat)
wtEbFit <- eBayes(wtFit)


#####################################################################
# Developmental Stage Hits
#####################################################################

# test the effect of developmental stage
dsHits <- topTable(wtEbFit, 
                   coef = grep("devStage", colnames(coef(wtEbFit))), 
                   number=nrow(wtDat))

# how many have p-val < 0.00001
table(dsHits$adj.P < 0.00001)

# 63rd hit
dsHits[63, c("F", "adj.P.Val", "devStageP6")]


#####################################################################
# P2 and P10 hits
#####################################################################

# test the effect of developmental stage
p2hit <- topTable(wtEbFit, coef = "devStageP2", 
                  number=nrow(wtDat), sort.by="none")
p10hit <- topTable(wtEbFit, coef = "devStageP10", 
                   number=nrow(wtDat), sort.by="none")

# plot test statistics against one another
xyplot(p10hit$t ~ p2hit$t,
       xlab="P2 t-statistic", ylab="P10 t-statistic",
       scales=list(tck=c(1,0), 
                   limits=c(range(c(p2hit$t, p10hit$t)))),
       panel=function(...){
         panel.hexbinplot(...)
         panel.abline(a=0, b=1, col=2)
       })

# density plot of p-vals
# we expect to observe a greater difference 10 weeks after
# baseline compared to 2 weeks after baseline (and we do!)
densityplot(~p10hit$adj.P + p2hit$adj.P, 
            xlab="adjusted p-values",
            plot.points=F, 
            auto.key=list(x=0.95, y=0.95, corner=c(1,1)),
            scales=list(tck=c(1,0)))

# p-values < 0.001
addmargins(table(p2hit$adj.P<0.001, p10hit$adj.P<0.001,
                 dnn=c("P2", "P10")))

# how many have p-val < 0.00001
table(dsHits$adj.P < 0.00001)

# scatterplot of P10 pvals
p10by <- topTable(wtEbFit, coef = "devStageP10", 
                   number=nrow(wtDat), sort.by="none",
                   adjust.method=c("BY"))
pvals <- data.frame(raw=p10hit$P.Val, 
                    BH=p10hit$adj.P, 
                    BY=p10by$adj.P) 
# BH pvals are larger than raw pvals (as expected)
plot(pvals)


#####################################################################
# Contrasts
#####################################################################

cm1 <- makeContrasts(p10vp6=devStageP10-devStageP6, 
                    fourvp10=devStage4_weeks-devStageP10, 
                    levels=wtDesMat)
wtC1 <- contrasts.fit(wtFit, cm1)
wtEbC1 <- eBayes(wtC1)
c1Hits <- topTable(wtEbC1)

# plot top 4 hits
makeSP(miniDF(rownames(c1Hits)[1:4]), 
       subset=miniDF(rownames(c1Hits)[1:4])$gType=="wt")

# adjust p-vals
cutoff <- 0.001
wtResC1 <- decideTests(wtEbC1, p.value=cutoff, method="global")
summary(wtResC1)


#####################################################################
# TAKE HOME PROBLEM
# See if you can find one or more probes that 
# have some expression changes up to P6 
# and then hold steady all the way to 4_weeks
#####################################################################

# get contrast matrix
#   e16 = p2
#   p2 = p6
(cm2 <- makeContrasts(p2ve16=devStageP6,
                      p6vp2=devStageP6-devStageP2,
                      levels=wtDesMat))

# estimate contrasts
wtC2 <- contrasts.fit(wtFit, cm2)
wtEbC2 <- eBayes(wtC2)
c2Hits <- topTable(wtEbC2)

# plot top 4 hits
# there are clear changes in the beginning stages (as expected)
makeSP(miniDF(rownames(c2Hits)),
       subset=miniDF(rownames(c2Hits))$gType=="wt")

# direction of differences
wtResC2 <- decideTests(wtEbC2, p.value=cutoff, method="global")
summary(wtResC2)

# genes with expression changes up to P6
# (i.e. both contrasts in cm2 non-zero)
change <- names(which(wtResC2[, "p6vp2"]!=0 & wtResC2[, "p2ve16"]!=0))

# genes with steady expression after P6
# (i.e. contrasts in cm1 close to zero)
wtResC1 <- decideTests(wtEbC1, p.value=0.05, method="global")
steady <- names(which(wtResC1[, "p10vp6"]==0 & wtResC1[, "fourvp10"]==0))

# genes with expression changes up to P6 and steady afterward
(both <- intersect(change, steady))

# plot the results
makeSP(miniDF(both), subset=miniDF(both)$gType=="wt")





