#####################################################################
# STAT 540: Seminar 6
# High Volume Linear Models
# Author: Shannon Erdelyi
# Date: 12-02-2014
#####################################################################

# source code with
#   read 'dat' and 'des'
#   function 'miniDF(genes)' to subset the data
dir <- "/Users/shannonerdelyi/Dropbox/UBC/W2014/STAT 540"
source(paste0(dir, "/Seminars/photoRec_Functions.R"))

#####################################################################
# Compare results from limma to lm on probe 1440645_at 
#####################################################################

# get wild type data only for probe 1440645_at
wt <- subset(miniDF(rownames(dat)), gType=="wt" & gene=="1440645_at")
head(wt)

# linear model
mod <- lm(gExp ~ devStage, wt)
summary(mod)
# the parameter estimates are the same 
# the F statistic is larger (589 > 425.4)


