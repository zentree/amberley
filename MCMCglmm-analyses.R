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
prior = list(R = list(V = 0.007, nu = 0.002),
             G = list(G1 = list(V = 0.002, nu = 0.002), G2 = list(V = 0.001, nu = 0.002)))

ngvel.u <-  MCMCglmm(ngvel ~ 1, 
                     random = ~ animal + Block, 
                     family = 'gaussian',
                     data = amber,
                     pedigree = pedigree,
                     prior = prior,
                     verbose = FALSE,
                     pr = TRUE,
                     burnin = 10000,
                     nitt = 210000,
                     thin = 200)

plot(mcmc.list(ngvel.u$VCV))

# Heritability for normal velocity
h2.ngvel <- ngvel.u$VCV[, 'animal']/(ngvel.u$VCV[, 'animal'] + ngvel.u$VCV[, 'Block'] + ngvel.u$VCV[, 'units'])
posterior.mode(h2.ngvel)
# 0.1526903

HPDinterval(h2.ngvel)
# 0.03549398 0.4069589

plot(h2.ngvel)

autocorr(ngvel.u$VCV)


#
# nbden (normal basic density)
#

nbden.u <-  MCMCglmm(nbden ~ 1, 
                     random = ~ animal + Block, 
                     family = 'gaussian',
                     data = amber,
                     pedigree = pedigree,
                     prior = prior,
                     verbose = FALSE,
                     pr = TRUE,
                     burnin = 10000,
                     nitt = 210000,
                     thin = 200)

plot(mcmc.list(nbden.u$VCV))

# Heritability for normal velocity
h2.nbden <- nbden.u$VCV[, 'animal']/(nbden.u$VCV[, 'animal'] + nbden.u$VCV[, 'Block'] + nbden.u$VCV[, 'units'])
posterior.mode(h2.nbden)
# 0.1570195

HPDinterval(h2.nbden)
# 0.06695 0.40988

plot(h2.nbden)

autocorr(nbden.u$VCV)


#
# Bivariate example (for correlation)
#

bivar.prior <- list(R = list(V = diag(c(0.01,0.02)), nu = 2),
                    G = list(G1 = list(V = diag(c(0.01,0.02)), nu = 2), G2 = list(V = diag(c(0.02,0.04)), nu = 2)))

ngvel.nbden <- MCMCglmm(cbind(ngvel, nbden) ~ trait - 1, 
                        random = ~ us(trait):Block + us(trait):animal,
                        rcov = ~ us(trait):units, family = c('gaussian', 'gaussian'),
                        data = amber,
                        prior = bivar.prior,
                        pedigree = pedigree,
                        verbose = FALSE,
                        pr = TRUE,
                        burnin = 10000,
                        nitt = 210000,
                        thin = 200)


# Multivariate h2 ngvel (uses information from density as well)
h2m.ngvel <- ngvel.nbden$VCV[,'ngvel:ngvel.animal']/(ngvel.nbden$VCV[,'ngvel:ngvel.animal'] + 
             ngvel.nbden$VCV[,'ngvel:ngvel.Block'] + ngvel.nbden$VCV[,'ngvel:ngvel.units'])

plot(h2m.ngvel)
posterior.mode(h2m.ngvel)
# 0.2551462 
HPDinterval(h2m.ngvel, prob = 0.95)
# var1 0.1500659 0.503896

# Multivariate h2 nbden (uses information from ngvel as well)
h2m.nbden <- ngvel.nbden$VCV[,'nbden:nbden.animal']/(ngvel.nbden$VCV[,'nbden:nbden.animal'] + 
  ngvel.nbden$VCV[,'nbden:nbden.Block'] + ngvel.nbden$VCV[,'nbden:nbden.units'])

plot(h2m.nbden)
posterior.mode(h2m.nbden)
# 0.1690776 
HPDinterval(h2m.nbden, prob = 0.95)
# var1 0.03816083 0.4012651

# Genetic correlation
rg.ngvel.nbden <- ngvel.nbden$VCV[,'ngvel:nbden.animal']/(ngvel.nbden$VCV[,'ngvel:ngvel.animal'] * 
                  ngvel.nbden$VCV[,'nbden:nbden.animal'])


plot(rg.ngvel.nbden)
posterior.mode(rg.ngvel.nbden)
# 0.3971

HPDinterval(rg.ngvel.nbden, prob = 0.95)
#             lower     upper
# var1 -0.004812858 0.9689158
