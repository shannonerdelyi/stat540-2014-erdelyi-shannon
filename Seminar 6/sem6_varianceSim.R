#####################################################################
# STAT 540: Seminar 6
# Variance Simulation
# Author: Shannon Erdelyi
# Date: 12-02-2014
#####################################################################

library(lattice)

# constants
set.seed(8)
groups <- 2
genes <- 5
n <- 3

# generate data with different group means and gene variances 
realGroupMeans <- rnorm(groups)
realGeneVar <- runif(genes)

data <- rnorm(groups*genes*n, 
              mean=rep(realGroupMeans, each=genes*n),
              sd=rep(sqrt(realGeneVar), each=n, times=groups))

# create a data frame
df <- data.frame(groups=rep(letters[1:groups], each=genes*n),
                 n=rep(paste("n", 1:n, sep=""), times=groups*genes),
                 genes=rep(paste("g", 1:genes, sep=""), each=n, times=groups),
                 data)
head(df)

# check observed group means
# the means are spot on :)
obsMeans <- tapply(df$data, df$groups, mean)
xyplot(obsMeans ~ realGroupMeans, 
       panel=function(...){
         panel.xyplot(...)
         panel.abline(a=0, b=1, lty=2)
       })

# check observed gene variance 
# the variances are all over the place :(
obsVar <- tapply(df$data, df$genes, var)
xyplot(obsVar ~ realGeneVar, 
       panel=function(...){
         panel.xyplot(...)
         panel.abline(a=0, b=1, lty=2)
       })

densityplot(obsVar)