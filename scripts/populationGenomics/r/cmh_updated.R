library(ggplot2)
library(poolSeq)
library(ACER)


#reps <- c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4,4))

#reps <- c(1:8)



# "~20 generations" per andrews methods but actually check what generation
# Given that controls are first, the first four are 0, rest are "~20 generations"
#gen <- c (rep(0, 4), rep(26, 12) )

#gen <- c(0,0,0,0,26,26,26,26)

reps <- c(1:8)
gen <- rep(26,8)

sync <- read_sync(file = "/Users/arteen/Desktop/pop_gen_out/actual_sync/AVH_subsetted_unmerged.sync", 
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

pops <- c("Ancestor1","Ancestor2",
          "Ancestor3","Ancestor4",
          "Up1", "Up2",
          "Up3","Up4")




af.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 8)
colnames(af.mat) <- pops

# for (i in 1:ncol(af.mat)){
#   if (i == 1 | i == 2 | i == 3 | i == 4 ) {
#     gen_num = 0
#   }else {gen_num = 26}
#   tempdat <- af(sync, repl = i, gen = gen_num)
#   af.mat[ , i] <- as.matrix(tempdat)
# }

for (i in 1:ncol(af.mat)){
  tempdat <- af(sync, repl = i, gen = 26)
  af.mat[,i] <- as.matrix(tempdat)
}

af.mat <- na.omit(af.mat)
head(af.mat)
dim(af.mat) # 1.9 million by 4

## lets just do ancestor vs Control first

af.mat2 <- af.mat[ , 2:3]
dim(af.mat2)
head(af.mat2)

# Now a coverage one?

cov.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 8)
colnames(cov.mat) <- pops

# for (i in 1:ncol(cov.mat)){
#   if (i == 1) {
#     gen_num = 0
#   }else {gen_num = 26}
#   tempdat <- poolSeq::coverage(sync, repl = i, gen = gen_num)
#   cov.mat[ , i] <- as.matrix(tempdat)
# }
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


cov.mat2 <- cov.mat[ , 2:3]
dim(cov.mat2)
head(cov.mat)

ne <- estimateNe(p0 = af.mat[ , "Ancestor1"],
                 
                 pt = af.mat[ , "Up1"], 
                 
                 cov0 = cov.mat[ , "Ancestor1"],
                 
                 covt = cov.mat[ , "Up1"], 
                 
                 t = 26, 
                 
                 method = "P.planII", 
                 
                 poolSize = c(96, 96))  

#Creating the vars for the CMH test 
rep <- c(1,1,1,1,2,2,2,2)
Ne <- as.integer(rep(53,2))
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

data2$neg.log10 <- -log10(data2$padj)



data2 <- na.omit(data2)
AVL_CMH <- data2

data2 <- UVD_CMH

data <- data2[data2$padj < quantile(data2$padj, 0.01), ]

top_0.01_data <- data2[data2$padj < quantile(data2$padj, 0.01), ]


data$chromEnd <- data[, 2]


data <- data[, c(1, 2, 6)]


headers <- c("chrom","chromStart","chromEnd")

colnames(data) <- headers
head(data)
write.table(x = data, file = "/Users/arteen/Desktop/pop_gen_out/bedtools/ancestorVhigh/cmh_top0.1.bed", sep = '\t',col.names = F, row.names = FALSE, quote = FALSE)


## if you wanna give it the +/- 5 bp
data$temp <- data$chromStart - 5
data$temp2 <- data$chromEnd + 5
data <- data[c(1,4,5)]


colnames(data) <- headers
head(data)

write.table(x = data, file = "/Users/arteen/Desktop/pop_gen_out/bedtools/new_cmh_top0.1.bed", sep = '\t',col.names = F, row.names = FALSE, quote = FALSE)
