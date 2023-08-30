library(tximport)
library(readr)
library("RColorBrewer")
library(stringr)
library(edgeR)
library(limma)
library(lme4)
library(glmmTMB)
library(emmeans)
library(parallel)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(ggrepel)
library(data.table)
library(car)
library(DESeq2)
library(tidyr)
library(klaR)


#To the directory of counts
setwd("/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data")
#setwd("/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/input")
quant_files <- file.path("quants", list.files("quants"), "quant.sf")
file.exists(quant_files)

samples <- quant_files
samples <- substring(samples, 8, 26) #grab just sample name from the file name

names(quant_files) <- samples


txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(x = txdb, keys = k, "GENEID", "TXNAME")

male_quant_filesIndex <- (grep(pattern = "*_M", quant_files))

male_quant_files <- quant_files[c(male_quant_filesIndex)]
txi <- tximport(male_quant_files, type = "salmon", tx2gene = tx2gene)

cts <- txi$counts

selection <- substring(samples, 4, 4)
selection <- gsub("C", "control", selection)
selection <- gsub("D", "down", selection)
selection <- gsub("U", "up", selection)
selection <- as.factor(selection)

lineage <- substring(samples, 7, 7)
lineage <- as.factor(lineage)

sex <- substring(samples, 9, 9)
sex <- gsub("F", "female", sex)
sex <- gsub("M", "male", sex)
sex <- as.factor(sex)

expMatched <- substring(samples, 11, 11)
expMatched <- gsub("S", "socArena", expMatched)
expMatched <- gsub("V", "vial", expMatched)
expMatched <- as.factor(expMatched)

lane <- substring(samples, 19, 19)
lane <- factor(lane)

rna.design <- data.frame(sample = samples,
                         file = quant_files,
                         selection = selection,
                         sex = sex,
                         expMatched = expMatched,
                         lane = lane,
                         lineage = lineage)

rna.design$selection <- relevel(as.factor(rna.design$selection), ref = "down")

rna.design$lineage2 <- paste(rna.design$selection, rna.design$lineage, sep = "")

males_rna.design <- subset(rna.design, rna.design$sex == "male")
males_rna.design <- males_rna.design[,-4]
## For limma Voom
design <- model.matrix(~ lane +
                         lineage +
                         expMatched +
                         selection +
                         sex:expMatched +
                         selection:expMatched +
                         selection:sex,
                       data = rna.design)


design_simple_PCA <- model.matrix(~ lane +
                                    expMatched +
                                    selection,
                                  data = males_rna.design)

dds <- DESeqDataSetFromTximport(txi, males_rna.design,
                                design = design_simple_PCA)

nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

vsd <- vst(dds, blind = FALSE, nsub = 5000) # using 5000 genes to do the vst, but all genes are back

#crap <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "parametric") # all genes

# Extracted and modified from the plot_pca function from RNAseqQC,  version 0.1.4 
extract_df_for_rda <- function (vsd, na_frac = 0.3, n_feats = 500, scale_feats = FALSE ) {
  mat <- tryCatch(as.matrix(assay(vsd)), error = function(e) {
    stop("Assay cannot be converted to a matrix.")
  })
  col_data <- as_tibble(colData(vsd))
  row_data <- as_tibble(rowData(vsd), rownames = "rowname")
  mat[is.nan(mat)] <- NA
  nrow_mat_original <- nrow(mat)
  mat <- mat[rowSums(is.na(mat))/ncol(mat) <= na_frac, ]
  if (nrow(mat) == 0) {
    stop("There are no features passing the 'na_frac' filtering criterion.")
  }
  else if (nrow(mat)/nrow_mat_original < 1/100) {
    warning("Less than 1% of the features pass the 'na'frac filtering criterion.")
  }
  mat_centered <- mat %>% t() %>% scale(scale = scale_feats) %>% 
    t()
  mat_centered[is.na(mat_centered)] <- 0
  mat_centered <- mat_centered[matrixStats::rowVars(mat_centered, 
                                                    na.rm = T) %>% order(decreasing = T) %>% head(n_feats), 
  ]
  mat_centered  
}


vst_extracted_top500 <- extract_df_for_rda(vsd, n_feats = 500)

vst_extracted_top1000 <- extract_df_for_rda(vsd, n_feats = 1000)

vst_extracted_top5000 <- extract_df_for_rda(vsd, n_feats = 5000)


t_vst_extracted_tmp <- t(vst_extracted_top500)



 
## this is from before, but it doesnt NOT work..
row_forSamples <- rownames(t_vst_extracted_tmp)


sample_rows <- (grep(pattern = "*M", row_forSamples))

sample_names <- row_forSamples[sample_rows]

t_vst_extracted <- t_vst_extracted_tmp[sample_rows,]

####Control Lineage 1#######################################################################################

for_test_set <- (grep(pattern = "*C_L1", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)



test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)


klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

#       predicted
# actual  C  D  U
#      C 18  0  0
#      D  0 24  0
#      U  0  0 24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")

# So this can identify males and females, but not really the treatment (in this lineage)
table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

#       predicted
# actual C D U
#      C 0 0 6



####Control Lineage 2#######################################################################################

for_test_set <- (grep(pattern = "*C_L2", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)


test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)

klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

#       predicted
# actual  C  D  U
#      C 18  0  0
#      D  0 24  0
#      U  0  0 24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

#       predicted
# actual C D U
#      C 0 2 4



####Control Lineage 3#######################################################################################

for_test_set <- (grep(pattern = "*C_L3", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)


test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)

klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

#       predicted
# actual  C  D  U
#      C 18  0  0
#      D  0 24  0
#      U  0  0 24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

#       predicted
# actual C D U
#      C 0 4 2


####Control Lineage 4#######################################################################################

for_test_set <- (grep(pattern = "*C_L4", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)

test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)

klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

#       predicted
# actual  C  D  U
#      C 18  0  0
#      D  0 24  0
#      U  0  0 24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

#       predicted
# actual C D U
#      C 2 0 4


####Up Lineage 1#######################################################################################

for_test_set <- (grep(pattern = "*U_L1", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)


test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)

klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

# Ran through a lot of iterations

#       predicted
# actual  C  D  U
#      C 24  0  0
#      D  0 24  0
#      U  0  0 18

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )


#       predicted
# actual C D U
#      U 4 0 2

####Up Lineage 2#######################################################################################

for_test_set <- (grep(pattern = "*U_L2", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)


test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)

klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

#       predicted
# actual  C  D  U
#      C 24  0  0
#      D  0 24  0
#      U  0  0 18

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

#       predicted
# actual C D U
#      U 1 5 0


####Up Lineage 3#######################################################################################

for_test_set <- (grep(pattern = "*U_L3", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)


test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)

klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

#       predicted
# actual  C  D  U
#      C 24  0  0
#      D  0 24  0
#      U  0  0 18

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

#       predicted
# actual C D U
#      U 0 6 0


####Up Lineage 4#######################################################################################

for_test_set <- (grep(pattern = "*U_L4", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)


test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)

klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

# predicted
# actual  C  D  U
# C 24  0  0
# D  0 24  0
# U  0  0 18

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

#       predicted
# actual C D U
#      U 0 6 0


####Down Lineage 1#######################################################################################

for_test_set <- (grep(pattern = "*D_L1", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)


test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)

klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

#       predicted
# actual  C  D  U
#      C 24  0  0
#      D  0 18  0
#      U  0  0 24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

#       predicted
# actual C D U
#      D 0 1 5


####Down Lineage 2#######################################################################################

for_test_set <- (grep(pattern = "*D_L2", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)


test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)

klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

#       predicted
# actual  C  D  U
#      C 24  0  0
#      D  0 18  0
#      U  0  0 24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )



#       predicted
# actual C D U
#      D 0 3 3


####Down Lineage 3#######################################################################################

for_test_set <- (grep(pattern = "*D_L3", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)


test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)

klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

#       predicted
# actual  C  D  U
#      C 24  0  0
#      D  0 18  0
#      U  0  0 24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

#       predicted
# actual C D U
#      D 0 2 4

# Not bad

####Down Lineage 4#######################################################################################

for_test_set <- (grep(pattern = "*D_L4", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

training_classifiers <- substring(row.names(training_set), 4, 4)


test_set <- t_vst_extracted[c(for_test_set), ]

test_classifiers <- substring(row.names(test_set), 4, 4)

klar_rda_out <- klaR::rda(x = training_set,
                          grouping = training_classifiers,
                          crossval = TRUE,
                          fold = 10,
                          output = TRUE)


## training set
klar_predict_out_training_class <- predict(klar_rda_out, 
                                           training_set,
                                           type = "class")

table(actual = training_classifiers,
      predicted = klar_predict_out_training_class$class)

#       predicted
# actual  C  D  U
#      C 24  0  0
#      D  0 18  0
#      U  0  0 24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

#       predicted
# actual C D U
#      D 0 0 6


# same as male + female RDA run






