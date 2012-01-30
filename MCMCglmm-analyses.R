# Genetic analyses for Amberley
# Luis A. Apiolaza
# School of Forestry
# University of Canterbury
#

library(MCMCglmm)

#### Setting Working directory
# This directory has to be changed to reflect
# the location of the files in your computer
setwd('~/Documents/Research/2012/amberley/')

# Read data file
amber <- read.csv('data20091016.csv', header = TRUE)
amber$Block <- factor(amber$Block)
amber$tree.id <- factor(10001:(10000+dim(amber)[1]))
head(amber)

# Read pedigree file
ped <- read.csv('pedigree.csv', header = TRUE)
head(ped)

# Rename tree.id to animal so we can use MCMCglmm
names(amber)[14] <- 'animal'
names(ped)[1] <- 'animal'

## Example for univariate analyses
#
# ngvel (normal green velocity)
#

prior = list(R = list(V = 0.007, n = 0),
             G = list(G1 = list(V = 0.002, n = 0), G2 = list(V = 0.001, n = 0)))

# I use high thinning to avoid autocorrelation with animal model
ngvel.u <-  MCMCglmm(ngvel ~ 1, 
                     random = ~ animal + Block, 
                     family = 'gaussian',
                     data = amber,
                     pedigree = ped,
                     prior = prior,
                     verbose = FALSE,
                     pr = TRUE,
                     burnin = 10000,
                     nitt = 200000,
                     thin = 200)

plot(mcmc.list(nvel.bayes$VCV))

# Heritability for normal velocity
h2.nvel <- nvel.bayes$VCV[, 'animal']/(nvel.bayes$VCV[, 'animal'] + nvel.bayes$VCV[, 'Block'] + nvel.bayes$VCV[, 'units'])
posterior.mode(h2.nvel)
# 0.196

HPDinterval(h2.nvel)
# 0.047 0.457

autocorr(nvel.bayes$VCV)
plot(h2.nvel)


