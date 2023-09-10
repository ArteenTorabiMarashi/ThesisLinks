library(ggplot2)
library(poolSeq)
library(ACER)


## Adapted from script by Tyler Audet

reps <- c(1:8)
gen <- rep(26,8)

sync <- read_sync(file = "/Users/arteen/Desktop/pop_gen_out/actual_sync/CVU_subsetted_X.sync", 
                  gen=gen,
                  repl=reps, 
                  polarization = "minor", 
                  keepOnlyBiallelic = TRUE)

# pops <- c("Ancestor1","Ancestor2",
#           "Ancestor3","Ancestor4",
#           "Control1","Control2",
#           "Control3","Control4",
#           "Down1", "Down2",
#           "Down3","Down4",
#           "Up1","Up2",
#           "Up3","Up4")

# pops <- c("Down1","Down2",
#           "Down3","Down4",
#           "Up1", "Up2",
#           "Up3","Up4")


# pops <- c("Control1","Control2",
#           "Control3","Control4",
#           "Down1", "Down2",
#           "Down3","Down4")

# pops <- c("Control1","Control2",
#           "Control3","Control4",
#           "Up1", "Up2",
#           "Up3","Up4")

# pops <- c("Ancestor1","Ancestor2",
#           "Ancestor3","Ancestor4",
#           "Up1", "Up2",
#           "Up3","Up4")

# pops <- c("Ancestor1","Ancestor2",
#           "Ancestor3","Ancestor4",
#           "Down1", "Down2",
#           "Down3","Down4")


af.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 8)
colnames(af.mat) <- pops


for (i in 1:ncol(af.mat)){
  tempdat <- af(sync, repl = i, gen = 26)
  af.mat[,i] <- as.matrix(tempdat)
}

af.mat <- na.omit(af.mat)
head(af.mat)
dim(af.mat) # 1.9 million by 4


# Now a coverage one

cov.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 8)
colnames(cov.mat) <- pops

for (i in 1:ncol(cov.mat)){
  tempdat <- poolSeq::coverage(sync, repl = i, gen = 26)
  cov.mat[,i] <- as.matrix(tempdat)
}

crap <- data.frame(cov.mat, sync@alleles[ , 1:2])
crap[crap == 0] <- NA
crap2 <- na.omit(crap)
location <- crap2[ , 9:10]

cov.mat[cov.mat == 0] <- NA
cov.mat <- na.omit(cov.mat)
dim(cov.mat) # still matches up, thats good
head(cov.mat)

ne <- estimateNe(p0 = af.mat[ , "Control1"],
                 
                 pt = af.mat[ , "Up1"], 
                 
                 cov0 = cov.mat[ , "Control1"],
                 
                 covt = cov.mat[ , "Up1"], 
                 
                 t = 26, 
                 
                 method = "P.planII", 
                 
                 poolSize = c(96, 96))  

#Creating the vars for the CMH test 
rep <- c(1,1,1,1,2,2,2,2)
Ne <- as.integer(rep(24,2))
tp <- c(rep(0,4), rep(1, 4))
ps <- c(rep(96,8))


pval <- adapted.cmh.test(freq = af.mat,
                         coverage = cov.mat, 
                         Ne = Ne,
                         gen = tp,
                         repl = rep,
                         poolSize = ps,
                         order = 1)


padj <- p.adjust(pval, method = "BY")

data2 <- cbind(location, pval, padj)





data2 <- na.omit(data2)


data <- data2[data2$padj < quantile(data2$padj, 0.01), ]


data$chromEnd <- data[, 2]


data <- data[, c(1, 2, 5)]


headers <- c("chrom","chromStart","chromEnd")

colnames(data) <- headers
head(data)

X_top1 <- data

############## Now for the autosomes

sync <- read_sync(file = "/Users/arteen/Desktop/pop_gen_out/actual_sync/CVU_subsetted_autosomes.sync", 
                  gen=gen,
                  repl=reps, 
                  polarization = "minor", 
                  keepOnlyBiallelic = TRUE)


af.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 8)
colnames(af.mat) <- pops


for (i in 1:ncol(af.mat)){
  tempdat <- af(sync, repl = i, gen = 26)
  af.mat[,i] <- as.matrix(tempdat)
}

af.mat <- na.omit(af.mat)
head(af.mat)
dim(af.mat) # 1.9 million by 4


# Now a coverage one

cov.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 8)
colnames(cov.mat) <- pops

for (i in 1:ncol(cov.mat)){
  tempdat <- poolSeq::coverage(sync, repl = i, gen = 26)
  cov.mat[,i] <- as.matrix(tempdat)
}

crap <- data.frame(cov.mat, sync@alleles[ , 1:2])
crap[crap == 0] <- NA
crap2 <- na.omit(crap)
location <- crap2[ , 9:10]

cov.mat[cov.mat == 0] <- NA
cov.mat <- na.omit(cov.mat)
dim(cov.mat) # still matches up, thats good
head(cov.mat)

ne <- estimateNe(p0 = af.mat[ , "Control1"],
                 
                 pt = af.mat[ , "Up1"], 
                 
                 cov0 = cov.mat[ , "Control1"],
                 
                 covt = cov.mat[ , "Up1"], 
                 
                 t = 26, 
                 
                 method = "P.planII", 
                 
                 poolSize = c(96, 96))  

#Creating the vars for the CMH test 
rep <- c(1,1,1,1,2,2,2,2)
Ne <- as.integer(rep(24,2))
tp <- c(rep(0,4), rep(1, 4))
ps <- c(rep(96,8))


pval <- adapted.cmh.test(freq = af.mat,
                         coverage = cov.mat, 
                         Ne = Ne,
                         gen = tp,
                         repl = rep,
                         poolSize = ps,
                         order = 1)


padj <- p.adjust(pval, method = "BY")

data2 <- cbind(location, pval, padj)
data2 <- na.omit(data2)



data <- data2[data2$padj < quantile(data2$padj, 0.01), ]




data$chromEnd <- data[, 2]


data <- data[, c(1, 2, 5)]


headers <- c("chrom","chromStart","chromEnd")

colnames(data) <- headers
head(data)

autosome_top1 <- data




###################### Now writing


results <- rbind(autosome_top1, X_top1)


write.table(x = results, file = "/Users/arteen/Desktop/pop_gen_out/revisions/controlVhigh/cmh_top0.1.bed",
            sep = '\t',col.names = F, row.names = FALSE, quote = FALSE)
