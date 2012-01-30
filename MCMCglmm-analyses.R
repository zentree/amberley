# Genetic analyses for Amberley
# Luis A. Apiolaza
# School of Forestry
# University of Canterbury
#

# Changing default conversion of strings to factors
# (source of extreme pain when debugging)
options(stringsAsFactors = FALSE)

# Loading MCMCglmm for Bayesian analyses
library(MCMCglmm)

# Rename tree.id to animal so we can use MCMCglmm
names(amber)[14] <- 'animal'
names(pedigree)[1] <- 'animal'

## Example for univariate analyses
#
# ngvel (normal green velocity)
#

# Priors for residuals, animal and block effects
prior = list(R = list(V = 0.007, n = 0),
             G = list(G1 = list(V = 0.002, n = 0), G2 = list(V = 0.001, n = 0)))

# I use high thinning to avoid autocorrelation with animal model
ngvel.u <-  MCMCglmm(ngvel ~ 1, 
                     random = ~ animal + Block, 
                     family = 'gaussian',
                     data = amber,
                     pedigree = pedigree,
                     prior = prior,
                     verbose = FALSE,
                     pr = TRUE,
                     burnin = 10000,
                     nitt = 200000,
                     thin = 200)

plot(mcmc.list(ngvel.u$VCV))

# Heritability for normal velocity
h2.ngvel <- ngvel.u$VCV[, 'animal']/(ngvel.u$VCV[, 'animal'] + ngvel.u$VCV[, 'Block'] + ngvel.u$VCV[, 'units'])
posterior.mode(h2.ngvel)
# 0.196

HPDinterval(h2.ngvel)
# 0.047 0.457

plot(h2.ngvel)

autocorr(ngvel.u$VCV)



