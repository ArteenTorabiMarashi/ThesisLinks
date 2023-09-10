# How many SNPs
# Adapted from script written by Ian


#https://inparanoidb.sbc.su.se/
# ^ Found this site to get an estimate of shared genes

## Says bombus is 6306
## Says humans are 5099
## Couldn't find apis mellifera so just going to use the number for bombus and hope thats alright
#, give that we are doing +/- 10%

# Bombus / Apis
## 6306 / 13701 (number of mapped genes we have) ~ 46%
# So we'll do 36% (4932), 46% (6306) and 56% (7672)

#Humans
## 5099 / 13701 ~ 37%
# we'll do 27% (3699), 37% (5099), and 47% (6439)





# These defaults are for humans and drosophila
sharedSNPs <- function( totalSNPs = 5099, sigSitesMine = 327, sigSitesOtherPaper = 56) {
  sites1 <- sample(totalSNPs, size = sigSitesMine, replace = F)
  sites2 <- sample(totalSNPs, size = sigSitesOtherPaper, replace = F)
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





# Bralten (Human)
# We found 0 overlapping

braltenSites_plus10 <- replicate(10^4, sharedSNPs(totalSNPs = 6439))
max(braltenSites_plus10) # 10
min(braltenSites_plus10) # 0
quantile(braltenSites_plus10, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100% 
# 0    0    0    1    3    5    6    7   10
mean(braltenSites_plus10 > 0) # 0.9401

braltenSites <- replicate(10^4, sharedSNPs(totalSNPs = 5099))
max(braltenSites) # 15
min(braltenSites) # 0
quantile(braltenSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100% 
# 0    0    1    1    3    6    7    8   12  
mean(braltenSites > 0) # 0.9749

braltenSites_minus10 <- replicate(10^4, sharedSNPs(totalSNPs = 3699))
max(braltenSites_minus10) # 14
min(braltenSites_minus10) # 0
quantile(braltenSites_minus10, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100% 
# 0    1    2    2    5    8    9   10   14
mean(braltenSites_minus10 > 0) # 0.9934


# Wang et al
# We found 2 overlapping

wangSites_plus10 <- replicate(10^4, sharedSNPs(totalSNPs = 7672, sigSitesOtherPaper = 115))
max(wangSites_plus10) # 16
min(wangSites_plus10) # 0
quantile(wangSites_plus10, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100% 
# 0    1    2    2    5    8    9   10   16
mean(wangSites_plus10 > 0) # 0.994
mean(wangSites_plus10 <= 2) # 0.1256 ~ proportion of simulations with 2 or less sites in common

wangSites <- replicate(10^4, sharedSNPs(totalSNPs = 6306, sigSitesOtherPaper = 115))
max(wangSites) # 17
min(wangSites) # 0
quantile(wangSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100% 
# 0    1    2    3    6    9   10   12   17
mean(wangSites > 0) # 0.9981
mean(wangSites <= 2) # 0.0554 ~ proportion of simulations with 2 or less sites in common

wangSites_minus10 <- replicate(10^4, sharedSNPs(totalSNPs = 4932, sigSitesOtherPaper = 115))
max(wangSites_minus10) # 19
min(wangSites_minus10) # 0
quantile(wangSites_minus10, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100% 
# 0    2    4    4    8   11   12   14   19 
mean(wangSites_minus10 > 0) # 0.9997
mean(wangSites_minus10 <= 2) # 0.0143 ~ proportion of simulations with 2 or less sites in common




## Sphigler et al
# We found 14 overlapping

sphigSites_plus10 <- replicate(10^4, sharedSNPs(totalSNPs = 7672, sigSitesOtherPaper = 1057))
max(sphigSites_plus10) # 72
min(sphigSites_plus10) # 23
quantile(sphigSites_plus10, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100%  
# 23   32   35   37   45   53   55   60   72  
mean(sphigSites_plus10 > 0) # 100%
mean(sphigSites_plus10 <= 14) # 0 ~ proportion of simulations with 14 or less sites in common

shpigSites <- replicate(10^4, sharedSNPs(totalSNPs = 6306, sigSitesOtherPaper = 1057))
max(shpigSites) # 81
min(shpigSites) # 28
quantile(shpigSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100% 
# 28   40   44   46   55   63   66   71   84 
mean(shpigSites > 0) # 100%
mean(shpigSites <= 14) # 0 ~ proportion of simulations with 14 or less sites in common

shpigSites_minus10 <- replicate(10^4, sharedSNPs(totalSNPs = 4932, sigSitesOtherPaper = 1057))
max(shpigSites_minus10) # 96
min(shpigSites_minus10) # 42
quantile(shpigSites_minus10, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100% 
# 42   53   58   61   70   79   82   87   96 
mean(shpigSites_minus10 > 0) # 100%
mean(shpigSites_minus10 <= 14) # 0 ~ proportion of simulations with 14 or less sites in common




## woodard et al
# We found 4 overlapping

woodardSites_plus10 <- replicate(10^4, sharedSNPs(totalSNPs = 7672, sigSitesOtherPaper = 212))
max(woodardSites_plus10) # 21
min(woodardSites_plus10) # 0
quantile(woodardSites_plus10, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100% 
# 0    3    5    5    9   13   14   16   21 
mean(woodardSites_plus10 > 0) # 0.9998
mean(woodardSites_plus10 <= 4) # 0.0474 ~ proportion of sites with 4 or less sites in common


woodardSites <- replicate(10^4, sharedSNPs(totalSNPs = 6306, sigSitesOtherPaper = 212))
max(woodardSites) # 25
min(woodardSites) # 0
quantile(woodardSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100% 
#  0    4    6    7   11   15   16   19   24  
mean(woodardSites > 0) # 0.9999
mean(woodardSites <= 4) # 0.0147 ~ proportion of sites with 4 or less sites in common

woodardSites_minus10 <- replicate(10^4, sharedSNPs(totalSNPs = 4932, sigSitesOtherPaper = 212))
max(woodardSites_minus10) # 30
min(woodardSites_minus10) # 3
quantile(woodardSites_minus10, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
# 0%   1%   5%  10%  50%  90%  95%  99% 100% 
# 3    6    8   10   14   19   20   23   30 
mean(woodardSites_minus10 > 0) # 100%
mean(woodardSites_minus10 <= 4) # 0.0016 ~ proportion of sites with 4 or less sites in common



## Between total DGE and DTU

dgeVdtu <- replicate(10^4, sharedSNPs(totalSNPs = 13701, sigSitesMine = 327, sigSitesOtherPaper = 619))
max(dgeVdtu)
min(dgeVdtu)
quantile(dgeVdtu, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
mean(dgeVdtu > 0)

############# Between the four methods of DGE


## Gauss Salmon - 174
## Gauss STAR - 183
## nbinom Salmon - 142
## nbinom STAR - 154



gSalmonVgSTAR <- replicate(10^4, sharedSNPs(totalSNPs = 11525, sigSitesMine = 174, sigSitesOtherPaper = 183))
max(gSalmonVgSTAR)
min(gSalmonVgSTAR)
quantile(gSalmonVgSTAR, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
mean(gSalmonVgSTAR > 0)

gSalmonVnSalmon <- replicate(10^4, sharedSNPs(totalSNPs = 11525, sigSitesMine = 174, sigSitesOtherPaper = 142))
max(gSalmonVnSalmon)
min(gSalmonVnSalmon)
quantile(gSalmonVnSalmon, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
mean(gSalmonVnSalmon > 0)

gSalmonVnSTAR <- replicate(10^4, sharedSNPs(totalSNPs = 11525, sigSitesMine = 174, sigSitesOtherPaper = 154))
max(gSalmonVnSTAR)
min(gSalmonVnSTAR)
quantile(gSalmonVnSTAR, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
mean(gSalmonVnSTAR > 0)

gSTARVnSalmon <- replicate(10^4, sharedSNPs(totalSNPs = 11525, sigSitesMine = 183, sigSitesOtherPaper = 142))
max(gSTARVnSalmon)
min(gSTARVnSalmon)
quantile(gSTARVnSalmon, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
mean(gSTARVnSalmon > 0)


gSTARVnSTAR <- replicate(10^4, sharedSNPs(totalSNPs = 11525, sigSitesMine = 183, sigSitesOtherPaper = 154))
max(gSTARVnSTAR)
min(gSTARVnSTAR)
quantile(gSTARVnSTAR, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
mean(gSTARVnSTAR > 0)


nSalmonVnSTAR <- replicate(10^4, sharedSNPs(totalSNPs = 11525, sigSitesMine = 142, sigSitesOtherPaper = 154))
max(nSalmonVnSTAR)
min(nSalmonVnSTAR)
quantile(nSalmonVnSTAR, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
mean(nSalmonVnSTAR > 0)

############# Between the our 3 DGE contrasts


## Low V High - 174
## Low V Control - 271
## Control V High - 194



lvhVlvc <- replicate(10^4, sharedSNPs(totalSNPs = 11525, sigSitesMine = 174, sigSitesOtherPaper = 271))
max(lvhVlvc)
min(lvhVlvc)
quantile(lvhVlvc, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
mean(lvhVlvc > 0)

lvhVcvh <- replicate(10^4, sharedSNPs(totalSNPs = 11525, sigSitesMine = 174, sigSitesOtherPaper = 194))
max(lvhVcvh)
min(lvhVcvh)
quantile(lvhVcvh, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
mean(lvhVcvh > 0)


lvcVcvh <- replicate(10^4, sharedSNPs(totalSNPs = 11525, sigSitesMine = 271, sigSitesOtherPaper = 194))
max(lvcVcvh)
min(lvcVcvh)
quantile(lvcVcvh, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))
mean(lvcVcvh > 0)

