# How many SNPs
# Adapted from script written by Ian



sharedSNPs <- function( totalSNPs = 13701, sigSitesFirst = 324, sigSitesOther = 245) {
  sites1 <- sample(totalSNPs, size = sigSitesFirst, replace = F)
  sites2 <- sample(totalSNPs, size = sigSitesOther, replace = F)
  intersection_length <- length(intersect(sites1, sites2))
  return(intersection_length)
}

# total SNPs is all dros genes we are checking? and sig sites


# This is the formula for doing this, straight from Ians code

HowManySharedSites <- replicate(10^4, sharedSNPs())

max(HowManySharedSites) # the maximum number of shared sites among all the simulations

quantile(HowManySharedSites, 
         probs = c(0.5, 0.9, 0.95, 0.99))

mean(HowManySharedSites > 0) # what proportion of simulations have any sites in common.


## Low versus High VS Low versus Control

LVH_vs_LVC <- replicate(10^4, sharedSNPs(sigSitesFirst = 324, sigSitesOther = 245))
max(LVH_vs_LVC)
min(LVH_vs_LVC) 
quantile(LVH_vs_LVC, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))

mean(LVH_vs_LVC > 0)

## Low versus High VS Control vs High

LVH_vs_CVH <- replicate(10^4, sharedSNPs(sigSitesFirst = 324, sigSitesOther = 184))
max(LVH_vs_CVH)
min(LVH_vs_CVH) 
quantile(LVH_vs_CVH, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))

mean(LVH_vs_CVH > 0)


## Low versus Control VS Control vs High

LVC_vs_CVH <- replicate(10^4, sharedSNPs(sigSitesFirst = 245, sigSitesOther = 184))
max(LVC_vs_CVH)
min(LVC_vs_CVH) 
quantile(LVC_vs_CVH, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))

mean(LVC_vs_CVH > 0)