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

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)

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

## For limma Voom
design <- model.matrix(~ lane +
                         lineage +
                         sex +
                         expMatched +
                         selection +
                         sex:expMatched +
                         selection:expMatched +
                         selection:sex,
                       data = rna.design)


design_simple_PCA <- model.matrix(~ lane +
                                    sex +
                                    expMatched +
                                    selection +
                                    selection:sex,
                                  data = rna.design)

dds <- DESeqDataSetFromTximport(txi, rna.design,
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


t_vst_extracted <- t(vst_extracted_top500)



sample_names <- rownames(t_vst_extracted)

####Control Lineage 1#######################################################################################

for_test_set <- (grep(pattern = "*C_L1", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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
# actual C_F C_M D_F D_M U_F U_M
  # C_F  18   0   0   0   0   0
  # C_M   0  18   0   0   0   0
  # D_F   0   0  24   0   0   0
  # D_M   0   0   0  24   0   0
  # U_F   0   0   0   0  23   0
  # U_M   0   0   0   0   0  24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                  test_set,
                                  type = "class")

# So this can identify males and females, but not really the treatment (in this lineage)
table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # C_F   0   0   0   0   5   0
  # C_M   0   0   0   0   0   6



####Control Lineage 2#######################################################################################

for_test_set <- (grep(pattern = "*C_L2", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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
# actual C_F C_M D_F D_M U_F U_M
  # C_F  17   0   0   0   0   0
  # C_M   0  18   0   0   0   0
  # D_F   0   0  24   0   0   0
  # D_M   0   0   0  24   0   0
  # U_F   0   0   0   0  23   0
  # U_M   0   0   0   0   0  24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # C_F   1   0   4   0   1   0
  # C_M   0   0   0   0   0   6

# Predicted 1 of the controls!

####Control Lineage 3#######################################################################################

for_test_set <- (grep(pattern = "*C_L3", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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
# actual C_F C_M D_F D_M U_F U_M
  # C_F  17   0   0   0   0   0
  # C_M   0  18   0   0   0   0
  # D_F   0   0  24   0   0   0
  # D_M   0   0   0  24   0   0
  # U_F   0   0   0   0  23   0
  # U_M   0   0   0   0   0  24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # C_F   0   0   6   0   0   0
  # C_M   0   0   0   5   0   1


####Control Lineage 4#######################################################################################

for_test_set <- (grep(pattern = "*C_L4", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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
# actual C_F C_M D_F D_M U_F U_M
  # C_F  17   0   0   0   0   0
  # C_M   0  18   0   0   0   0
  # D_F   0   0  24   0   0   0
  # D_M   0   0   0  24   0   0
  # U_F   0   0   0   0  23   0
  # U_M   0   0   0   0   0  24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # C_F   0   0   0   0   6   0
  # C_M   0   0   0   0   0   6


####Up Lineage 1#######################################################################################

for_test_set <- (grep(pattern = "*U_L1", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # C_F  23   0   0   0   0   0
  # C_M   0  22   0   0   0   2
  # D_F   0   0  24   0   0   0
  # D_M   0   0   0  24   0   0
  # U_F   0   0   0   0  18   0
  # U_M   0   0   0   0   0  18

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )


      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # U_F   5   0   0   0   0   0
  # U_M   0   6   0   0   0   0

####Up Lineage 2#######################################################################################

for_test_set <- (grep(pattern = "*U_L2", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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
# actual C_F C_M D_F D_M U_F U_M
  # C_F  23   0   0   0   0   0
  # C_M   0  24   0   0   0   0
  # D_F   0   0  24   0   0   0
  # D_M   0   0   0  24   0   0
  # U_F   0   0   0   0  17   0
  # U_M   0   0   0   0   0  18

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # U_F   0   0   5   0   1   0
  # U_M   0   1   0   5   0   0


####Up Lineage 3#######################################################################################

for_test_set <- (grep(pattern = "*U_L3", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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
# actual C_F C_M D_F D_M U_F U_M
  # C_F  23   0   0   0   0   0
  # C_M   0  24   0   0   0   0
  # D_F   0   0  24   0   0   0
  # D_M   0   0   0  24   0   0
  # U_F   0   0   0   0  17   0
  # U_M   0   0   0   0   0  18

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # U_F   0   0   5   0   1   0
  # U_M   0   0   0   6   0   0


####Up Lineage 4#######################################################################################

for_test_set <- (grep(pattern = "*U_L4", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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
# actual C_F C_M D_F D_M U_F U_M
  # C_F  23   0   0   0   0   0
  # C_M   0  24   0   0   0   0
  # D_F   0   0  24   0   0   0
  # D_M   0   0   0  24   0   0
  # U_F   0   0   0   0  17   0
  # U_M   0   0   0   0   0  18

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # U_F   0   0   6   0   0   0
  # U_M   0   0   0   6   0   0


####Down Lineage 1#######################################################################################

for_test_set <- (grep(pattern = "*D_L1", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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
# actual C_F C_M D_F D_M U_F U_M
  # C_F  23   0   0   0   0   0
  # C_M   0  24   0   0   0   0
  # D_F   0   0  18   0   0   0
  # D_M   0   0   0  18   0   0
  # U_F   0   0   0   0  23   0
  # U_M   0   0   0   0   0  24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # D_F   1   0   0   0   5   0
  # D_M   0   0   0   1   0   5


####Down Lineage 2#######################################################################################

for_test_set <- (grep(pattern = "*D_L2", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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
# actual C_F C_M D_F D_M U_F U_M
  # C_F  23   0   0   0   0   0
  # C_M   0  24   0   0   0   0
  # D_F   0   0  18   0   0   0
  # D_M   0   0   0  18   0   0
  # U_F   0   0   0   0  23   0
  # U_M   0   0   0   0   0  24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

# Actually predicts !!!

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # D_F   0   0   6   0   0   0
  # D_M   0   0   0   5   0   1


####Down Lineage 3#######################################################################################

for_test_set <- (grep(pattern = "*D_L3", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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
# actual C_F C_M D_F D_M U_F U_M
  # C_F  23   0   0   0   0   0
  # C_M   0  24   0   0   0   0
  # D_F   0   0  18   0   0   0
  # D_M   0   0   0  18   0   0
  # U_F   0   0   0   0  23   0
  # U_M   0   0   0   0   0  24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # D_F   1   0   3   0   2   0
  # D_M   0   0   0   5   0   1

# Not bad

####Down Lineage 4#######################################################################################

for_test_set <- (grep(pattern = "*D_L4", sample_names))

training_set <- t_vst_extracted[-c(for_test_set), ]

#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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
# actual C_F C_M D_F D_M U_F U_M
  # C_F  23   0   0   0   0   0
  # C_M   0  24   0   0   0   0
  # D_F   0   0  18   0   0   0
  # D_M   0   0   0  18   0   0
  # U_F   0   0   0   0  23   0
  # U_M   0   0   0   0   0  24

## now for test data

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")


table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )

      # predicted
# actual C_F C_M D_F D_M U_F U_M
  # D_F   0   0   0   0   6   0
  # D_M   0   0   0   0   0   6

























####### Below was me playing around with the functions#####








for_test_set <- (grep(pattern = "*D_L2", sample_names))



training_set <- t_vst_extracted[(for_test_set), ]



#training_classifiers <- substring(row.names(training_set), 4, 4)
training_classifiers_toMod <- substring(row.names(training_set), 4, 9)
training_classifiers1 <- substring(training_classifiers_toMod, 1,1)
training_classifiers2 <- substring(training_classifiers_toMod, 6,6)
training_classifiers <- paste(training_classifiers1, training_classifiers2, sep = "_")


test_set <- t_vst_extracted[-c(for_test_set), ]

#test_classifiers <- substring(row.names(test_set), 4, 4)
test_classifiers_toMod <- substring(row.names(test_set), 4, 9)
test_classifiers1 <- substring(test_classifiers_toMod, 1,1)
test_classifiers2 <- substring(test_classifiers_toMod, 6,6)
test_classifiers <- paste(test_classifiers1, test_classifiers2, sep = "_")

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

klar_predict_out_test_class <- predict(klar_rda_out,
                                       test_set,
                                       type = "class")

table(actual = test_classifiers,
      predicted = klar_predict_out_test_class$class )



















####Trying stuff################################################################################################
## Tutorial that Ian sent

# For control lineage 1, this performs the same as klaR

library(caret)
library(glmnet)
library(mlbench)



cv_5_grid = trainControl(method = "cv", number = 5)

fit_rda_grid = train(x = training_set,
                     y = training_classifiers,
                     method = "rda",
                     trControl = cv_5_grid)






cv_5_rand = trainControl(method = "cv", number = 5, search = "random")
fit_rda_rand = train(x = training_set,
                     y = training_classifiers,
                     method = "rda",
                     trControl = cv_5_rand,
                     tuneLength = 9)


fit_elnet_grid = train(x = training_set,
                       y = training_classifiers,
                       method = "glmnet",
                       trControl = cv_5_grid,
                       tuneLength = 10)

# fit_elnet_int_grid = train(Class ~ . ^ 2, data = Sonar, method = "glmnet",
#                            trControl = cv_5_grid, tuneLength = 10)


fit_rda_grid
fit_rda_rand
fit_elnet_grid

get_best_result = function(caret_fit) {
  best_result = caret_fit$results[as.numeric(rownames(caret_fit$bestTune)), ]
  rownames(best_result) = NULL
  best_result
}
knitr::kable(rbind(
  get_best_result(fit_rda_grid),
  get_best_result(fit_rda_rand)))



## training set
caret_predict_out_training_class <- predict(fit_rda_rand, 
                                           training_set,
                                           type = "raw")

table(actual = training_classifiers,
      predicted = caret_predict_out_training_class)

## now for test data

caret_predict_out_test_class <- predict(fit_rda_rand,
                                       test_set,
                                       type = "raw")

# So this can identify males and females, but not really the treatment (in this lineage)
table(actual = test_classifiers,
      predicted = caret_predict_out_test_class )




## sparsediscrim

rda_out<- rda_high_dim(x = training_set,
                           y = training_classifiers)

predict_out_training_class <- predict(rda_out,
                                      training_set,
                                      type = "class")

predict_out_training_prob <- predict(rda_out,
                                     training_set,
                                     type = "prob")
table(actual = training_classifiers,
      predicted = predict_out_training_class)

## now for test data

predict_out_test_class <- predict(rda_out,
                                  test_set,
                                  type = "class")

table(actual = test_classifiers,
      predicted = predict_out_test_class )

predict_out_test_class <- predict(rda_out,
                                  test_set,
                                  type = "prob")
