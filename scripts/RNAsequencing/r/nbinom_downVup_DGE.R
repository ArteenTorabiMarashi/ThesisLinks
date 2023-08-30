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


#To the directory of counts
#setwd("/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data")
setwd("/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/input")
quant_files <- file.path("quants", list.files("quants"), "quant.sf")
file.exists(quant_files)

samples <- quant_files
samples <- substring(samples, 8, 26) #grab just sample name from the file name

names(quant_files) <- samples


txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(x = txdb, keys = k, "GENEID", "TXNAME")

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)

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

## For Limma Voom
design <- model.matrix(~ lane + 
                         lineage +
                         sex + 
                         expMatched + 
                         selection + 
                         sex:expMatched + 
                         selection:expMatched + 
                         selection:sex,
                       data = rna.design)

raw.counts <- txi$counts
d0 <- DGEList(raw.counts)

#Remove low counts
keep2 <- filterByExpr(d0, design)
d0 <- d0[keep2, ]


d0 <- calcNormFactors(d0)

samples <- colnames(raw.counts)

#without random effect
y <- voom(d0,
          design,
          plot = FALSE)

fit <- lmFit(y, design)
fit <- eBayes(fit)

selectionDown.limma <- topTable(fit, coef = 8, sort.by = "none", n = Inf)

targets <- row.names(selectionDown.limma)



voomCPM <- y$E

nt <- min(parallel::detectCores(), 8)

cpm.genes.check <- list()


## Empty lists and dataframes for emmeans
selectionPairEmmean_df_nbinom <- data.frame()
experiencePairEmmean_df_nbinom <- data.frame()
sexPairEmmean_df_nbinom <- data.frame()
sexContrastEmmean_df_nbinom <- data.frame()
selectionContrastEmmean_df_nbinom <- data.frame()
experienceContrastEmmean_df_nbinom <- data.frame()
sexSelectionPairEmmean_df_nbinom <- data.frame()
expSelectionPairEmmean_df_nbinom <- data.frame()
selectionSexPairEmmean_df_nbinom <- data.frame()

selectionPair_df_nbinom <- data.frame()
experiencePair_df_nbinom <- data.frame()
sexPair_df_nbinom <- data.frame()
sexContrast_df_nbinom <- data.frame()
selectionContrast_df_nbinom <- data.frame()
experienceContrast_df_nbinom <- data.frame()
sexSelectionPair_df_nbinom <- data.frame()
expSelectionPair_df_nbinom <- data.frame()
selectionSexPair_df_nbinom <- data.frame()

sex_anova_df_nbinom <- data.frame(gene = character(),
                           pval = double(),
                           fullModel = character())

exp_anova_df_nbinom <- data.frame(gene = character(),
                           pval = double(),
                           fullModel = character())

sel_anova_df_nbinom <- data.frame(gene = character(),
                           pval = double(),
                           fullModel = character())

sexExp_anova_df_nbinom <- data.frame(gene = character(),
                              pval = double(),
                              fullModel = character())

expSel_anova_df_nbinom <- data.frame(gene = character(),
                              pval = double(),
                              fullModel = character())

sexSel_anova_df_nbinom <- data.frame(gene = character(),
                              pval = double(),
                              fullModel = character())

lost_genes_nbinom <- list()

## NOTE
# experiencePair_df_nbinom will double save the contrasts (because emmeans and contrasts are 12 & 6 so will not match)
# so it will save the contrasts twice, remember to remove the duplicates

load.model <- formula(~ lane + 
                         lineage +
                         sex + 
                         expMatched + 
                         selection + 
                         sex:expMatched + 
                         selection:expMatched + 
                         selection:sex)

data.in <- DESeqDataSetFromTximport(txi,rna.design,
                                    design = load.model)



data.in <- DESeq2::estimateSizeFactors(data.in)

normFactors <- DESeq2::normalizationFactors(data.in)





count <- 1
Sys.time()

for (g in targets) {
  
  count <- count + 1
  if (count%%100 == 0){
    print(count)
    print(Sys.time())
  }
  
  the_counts <- raw.counts[g,]
  the_counts <- round(the_counts, 0)
  the_counts <- unlist(the_counts)
  the_normFactors <- normFactors[g,]
  
  gaussian.model <- glmmTMB(the_counts ~ lane +
                              sex + 
                              expMatched + 
                              selection + 
                              sex:expMatched + 
                              selection:expMatched + 
                              selection:sex + 
                              diag(0 + sex + expMatched | selection:lineage),
                            offset = the_normFactors,
                            family = nbinom2(),
                            data = rna.design,
                            control = glmmTMBControl(parallel = nt))
 
  
  if (anyNA(as.vector(summary(gaussian.model)$coef$cond[,2:4])) == T |
      anyNA(logLik(gaussian.model)) == T) {
    ## Using reduced model
    gaussian.model.reduced <- glmmTMB(the_counts ~ lane +
                                        sex + 
                                        expMatched + 
                                        selection + 
                                        sex:expMatched + 
                                        selection:expMatched + 
                                        selection:sex + 
                                        diag(0 + sex | selection:lineage),
                                      offset = the_normFactors,
                                      family = nbinom2(),
                                      data = rna.design,
                                      control = glmmTMBControl(parallel = nt))
    
    
    if (anyNA(as.vector(summary(gaussian.model.reduced)$coef$cond[,2:4])) == T|
        anyNA(logLik(gaussian.model)) == T) {
      ## Using reduced model
      gaussian.model.reduced <- glmmTMB(the_counts ~ lane +
                                          sex + 
                                          expMatched + 
                                          selection + 
                                          sex:expMatched + 
                                          selection:expMatched + 
                                          selection:sex + 
                                          diag(1  + sex| selection:lineage),
                                        offset = the_normFactors,
                                        family = nbinom2(),
                                        data = rna.design,
                                        control = glmmTMBControl(parallel = nt))
      
      if (anyNA(as.vector(summary(gaussian.model.reduced)$coef$cond[,2:4])) == T|
          anyNA(logLik(gaussian.model)) == T) {
         
        lost_genes_nbinom <- append(lost_genes_nbinom, g)
         next;
        }
    }
    
    cpm.genes.check[[g]] <- list(summary(gaussian.model.reduced)$coef$cond, "Used Reduced Model")
    
    anova_out <- car::Anova(gaussian.model.reduced)
    
    sex_anova_res <- anova_out[2,3]
    sex_anova_insert <- data.frame(gene = geneID,
                                   pval = sex_anova_res,
                                   fullModel = "No")
    sex_anova_df_nbinom <- rbind(sex_anova_df_nbinom, sex_anova_insert)
    
    
    exp_anova_res <- anova_out[3,3]
    exp_anova_insert <- data.frame(gene = geneID,
                                   pval = exp_anova_res,
                                   fullModel = "No")
    exp_anova_df_nbinom <- rbind(exp_anova_df_nbinom, exp_anova_insert)
    
    
    sel_anova_res <- anova_out[4,3]
    sel_anova_insert <- data.frame(gene = geneID,
                                   pval = sel_anova_res,
                                   fullModel = "No")
    sel_anova_df_nbinom <- rbind(sel_anova_df_nbinom, sel_anova_insert)
    
    sexExp_anova_res <- anova_out[5,3]
    sexExp_anova_insert <- data.frame(gene = geneID,
                                      pval = sexExp_anova_res,
                                      fullModel = "No")
    sexExp_anova_df_nbinom <- rbind(sexExp_anova_df_nbinom, sexExp_anova_insert)
    
    expSel_anova_res <- anova_out[6,3]
    expSel_anova_insert <- data.frame(gene = geneID,
                                      pval = expSel_anova_res,
                                      fullModel = "No")
    expSel_anova_df_nbinom <- rbind(expSel_anova_df_nbinom, expSel_anova_insert)
    
    sexSel_anova_res <- anova_out[7,3]
    sexSel_anova_insert <- data.frame(gene = geneID,
                                      pval = sexSel_anova_res,
                                      fullModel = "No")
    sexSel_anova_df_nbinom <- rbind(sexSel_anova_df_nbinom, sexSel_anova_insert)
    
    
    
    geneID <- g
    
    # If the reduced model still fails, the anova wont work causing the entire script to fail, so just skipping it for now
    # ## Anova and extraction
    # 
    # anovaOut <- car::Anova(gaussian.model.reduced)
    # anovaOut_df_nbinom.tmp <- as.data.frame(anovaOut)
    # anovaOut_df_nbinom.tmp <- setDT(anovaOut_df_nbinom.tmp, keep.rownames = "Coefficient")
    # 
    # anovaOut_df_nbinom <- cbind(geneID, anovaOut_df_nbinom.tmp)
    # anova_df_nbinom <- rbind(anova_df_nbinom, anovaOut_df_nbinom)
    # 
    #emmeans contrasts
    selectionPair <- suppressMessages(emmeans(gaussian.model.reduced, pairwise ~ selection | sex + expMatched))
    experiencePair <- suppressMessages(emmeans(gaussian.model.reduced, pairwise ~ expMatched | sex + selection))
    sexPairEmm <- suppressMessages(emmeans(gaussian.model.reduced, pairwise ~ sex | expMatched + selection))
    sexEmm <- suppressMessages(emmeans(gaussian.model.reduced, pairwise ~ sex))
    selectionEmm <- suppressMessages(emmeans(gaussian.model.reduced, pairwise ~ selection))
    experienceEmm <- suppressMessages(emmeans(gaussian.model.reduced, pairwise ~ expMatched))
    sexSelectionPair <- suppressMessages(emmeans(gaussian.model.reduced, pairwise ~ sex | selection))
    expSelectionPair <- suppressMessages(emmeans(gaussian.model.reduced, pairwise ~ expMatched | selection))
    selectionSexEmm <- suppressMessages(emmeans(gaussian.model.reduced, pairwise ~ selection | sex))
    
    
    
    ## one for emmean one for contrast
    
    selectPairEmmean.insert <- selectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectPairEmmean <- cbind(geneID, selectPairEmmean.insert)
    selectionPairEmmean_df_nbinom <- rbind(selectionPairEmmean_df_nbinom, selectPairEmmean)
    
    selectPair.insert <- selectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectPair <- cbind(geneID, selectPair.insert)
    selectionPair_df_nbinom <- rbind(selectionPair_df_nbinom, selectPair)
    
    expPairEmmean.insert <- experiencePair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expPairEmmean <- cbind(geneID, expPairEmmean.insert)
    experiencePairEmmean_df_nbinom <- rbind(experiencePairEmmean_df_nbinom, expPairEmmean)
    
    expPair.insert <- experiencePair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expPair <- cbind(geneID, expPair.insert)
    experiencePair_df_nbinom <- rbind(experiencePair_df_nbinom, expPair)
    
    sexPairEmmean.insert <- sexPairEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexPairEmmean <- cbind(geneID, sexPairEmmean.insert)
    sexPairEmmean_df_nbinom <- rbind(sexPairEmmean_df_nbinom, sexPairEmmean)
    
    sexPair.insert <- sexPairEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexPair <- cbind(geneID, sexPair.insert)
    sexPair_df_nbinom <- rbind(sexPair_df_nbinom, sexPair)
    
    sexContrastEmmean.insert <- sexEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexContrastEmmean <- cbind(geneID, sexContrastEmmean.insert)
    sexContrastEmmean_df_nbinom <- rbind(sexContrastEmmean_df_nbinom, sexContrastEmmean)
    
    sexContrast.insert <- sexEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexContrast <- cbind(geneID, sexContrast.insert)
    sexContrast_df_nbinom <- rbind(sexContrast_df_nbinom, sexContrast)
    
    selectionContrastEmmean.insert <- selectionEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionContrastEmmean <- cbind(geneID, selectionContrastEmmean.insert)
    selectionContrastEmmean_df_nbinom <- rbind(selectionContrastEmmean_df_nbinom, selectionContrastEmmean)
    
    selectionContrast.insert <- selectionEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionContrast <- cbind(geneID, selectionContrast.insert)
    selectionContrast_df_nbinom <- rbind(selectionContrast_df_nbinom, selectionContrast)
    
    experienceContrastEmmean.insert <- experienceEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    experienceContrastEmmean <- cbind(geneID, experienceContrastEmmean.insert)
    experienceContrastEmmean_df_nbinom <- rbind(experienceContrastEmmean_df_nbinom, experienceContrastEmmean)
    
    experienceContrast.insert <- experienceEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    experienceContrast <- cbind(geneID, experienceContrast.insert)
    experienceContrast_df_nbinom <- rbind(experienceContrast_df_nbinom, experienceContrast)
    
    sexSelectionEmmean.insert <- sexSelectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexSelectionEmmean <- cbind(geneID, sexSelectionEmmean.insert)
    sexSelectionPairEmmean_df_nbinom <- rbind(sexSelectionPairEmmean_df_nbinom, sexSelectionEmmean)
    
    sexSelection.insert <- sexSelectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexSelection <- cbind(geneID, sexSelection.insert)
    sexSelectionPair_df_nbinom <- rbind(sexSelectionPair_df_nbinom, sexSelection)
    
    expSelectionEmmean.insert <- expSelectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expSelectionEmmean <- cbind(geneID, expSelectionEmmean.insert)
    expSelectionPairEmmean_df_nbinom <- rbind(expSelectionPairEmmean_df_nbinom, expSelectionEmmean)
    
    expSelection.insert <- expSelectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expSelection <- cbind(geneID, expSelection.insert)
    expSelectionPair_df_nbinom <- rbind(expSelectionPair_df_nbinom, expSelection)
    
    selectionSexPairEmmean.insert <- selectionSexEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionSexPairEmmean <- cbind(geneID, selectionSexPairEmmean.insert)
    selectionSexPairEmmean_df_nbinom <- rbind(selectionSexPairEmmean_df_nbinom, selectionSexPairEmmean)
    
    selectionSexPair.insert <- selectionSexEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionSexPair <- cbind(geneID, selectionSexPair.insert)
    selectionSexPair_df_nbinom <- rbind(selectionSexPair_df_nbinom, selectionSexPair)
  }
  else {
    # If using full model
    cpm.genes.check[[g]] <- list(summary(gaussian.model)$coef$cond, "Used Full Model")

    geneID <- g
    
    anova_out <- car::Anova(gaussian.model)
    
    sex_anova_res <- anova_out[2,3]
    sex_anova_insert <- data.frame(gene = geneID,
                                   pval = sex_anova_res,
                                   fullModel = "Yes")
    sex_anova_df_nbinom <- rbind(sex_anova_df_nbinom, sex_anova_insert)
    
    
    exp_anova_res <- anova_out[3,3]
    exp_anova_insert <- data.frame(gene = geneID,
                                   pval = exp_anova_res,
                                   fullModel = "Yes")
    exp_anova_df_nbinom <- rbind(exp_anova_df_nbinom, exp_anova_insert)
    
    
    sel_anova_res <- anova_out[4,3]
    sel_anova_insert <- data.frame(gene = geneID,
                                   pval = sel_anova_res,
                                   fullModel = "Yes")
    sel_anova_df_nbinom <- rbind(sel_anova_df_nbinom, sel_anova_insert)
    
    sexExp_anova_res <- anova_out[5,3]
    sexExp_anova_insert <- data.frame(gene = geneID,
                                      pval = sexExp_anova_res,
                                      fullModel = "Yes")
    sexExp_anova_df_nbinom <- rbind(sexExp_anova_df_nbinom, sexExp_anova_insert)
    
    expSel_anova_res <- anova_out[6,3]
    expSel_anova_insert <- data.frame(gene = geneID,
                                      pval = expSel_anova_res,
                                      fullModel = "Yes")
    expSel_anova_df_nbinom <- rbind(expSel_anova_df_nbinom, expSel_anova_insert)
    
    sexSel_anova_res <- anova_out[7,3]
    sexSel_anova_insert <- data.frame(gene = geneID,
                                      pval = sexSel_anova_res,
                                      fullModel = "Yes")
    sexSel_anova_df_nbinom <- rbind(sexSel_anova_df_nbinom, sexSel_anova_insert)
    
    
    # ## Anova and extraction
    # anovaOut <- car::Anova(gaussian.model)
    # anovaOut_df_nbinom.tmp <- as.data.frame(anovaOut)
    # anovaOut_df_nbinom.tmp <- setDT(anovaOut_df_nbinom.tmp, keep.rownames = "Coefficient")
    # 
    # anovaOut_df_nbinom <- cbind(geneID, anovaOut_df_nbinom.tmp)
    # anova_df_nbinom <- rbind(anova_df_nbinom, anovaOut_df_nbinom)
    # 
    #emmeans contrasts
    selectionPair <- suppressMessages(emmeans(gaussian.model, pairwise ~ selection | sex + expMatched))
    experiencePair <- suppressMessages(emmeans(gaussian.model, pairwise ~ expMatched | sex + selection))
    sexPairEmm <- suppressMessages(emmeans(gaussian.model, pairwise ~ sex | expMatched + selection))
    sexEmm <- suppressMessages(emmeans(gaussian.model, pairwise ~ sex))
    selectionEmm <- suppressMessages(emmeans(gaussian.model, pairwise ~ selection))
    experienceEmm <- suppressMessages(emmeans(gaussian.model, pairwise ~ expMatched))
    sexSelectionPair <- suppressMessages(emmeans(gaussian.model, pairwise ~ sex | selection))
    expSelectionPair <- suppressMessages(emmeans(gaussian.model, pairwise ~ expMatched | selection))
    selectionSexEmm <- suppressMessages(emmeans(gaussian.model, pairwise ~ selection | sex))
    
    
    
    ## one for emmean one for contrast
    
    selectPairEmmean.insert <- selectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectPairEmmean <- cbind(geneID, selectPairEmmean.insert)
    selectionPairEmmean_df_nbinom <- rbind(selectionPairEmmean_df_nbinom, selectPairEmmean)
    
    selectPair.insert <- selectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectPair <- cbind(geneID, selectPair.insert)
    selectionPair_df_nbinom <- rbind(selectionPair_df_nbinom, selectPair)
    
    expPairEmmean.insert <- experiencePair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expPairEmmean <- cbind(geneID, expPairEmmean.insert)
    experiencePairEmmean_df_nbinom <- rbind(experiencePairEmmean_df_nbinom, expPairEmmean)
    
    expPair.insert <- experiencePair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expPair <- cbind(geneID, expPair.insert)
    experiencePair_df_nbinom <- rbind(experiencePair_df_nbinom, expPair)
    
    sexPairEmmean.insert <- sexPairEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexPairEmmean <- cbind(geneID, sexPairEmmean.insert)
    sexPairEmmean_df_nbinom <- rbind(sexPairEmmean_df_nbinom, sexPairEmmean)
    
    sexPair.insert <- sexPairEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexPair <- cbind(geneID, sexPair.insert)
    sexPair_df_nbinom <- rbind(sexPair_df_nbinom, sexPair)
    
    sexContrastEmmean.insert <- sexEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexContrastEmmean <- cbind(geneID, sexContrastEmmean.insert)
    sexContrastEmmean_df_nbinom <- rbind(sexContrastEmmean_df_nbinom, sexContrastEmmean)
    
    sexContrast.insert <- sexEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexContrast <- cbind(geneID, sexContrast.insert)
    sexContrast_df_nbinom <- rbind(sexContrast_df_nbinom, sexContrast)
    
    selectionContrastEmmean.insert <- selectionEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionContrastEmmean <- cbind(geneID, selectionContrastEmmean.insert)
    selectionContrastEmmean_df_nbinom <- rbind(selectionContrastEmmean_df_nbinom, selectionContrastEmmean)
    
    selectionContrast.insert <- selectionEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionContrast <- cbind(geneID, selectionContrast.insert)
    selectionContrast_df_nbinom <- rbind(selectionContrast_df_nbinom, selectionContrast)
    
    experienceContrastEmmean.insert <- experienceEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    experienceContrastEmmean <- cbind(geneID, experienceContrastEmmean.insert)
    experienceContrastEmmean_df_nbinom <- rbind(experienceContrastEmmean_df_nbinom, experienceContrastEmmean)
    
    experienceContrast.insert <- experienceEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    experienceContrast <- cbind(geneID, experienceContrast.insert)
    experienceContrast_df_nbinom <- rbind(experienceContrast_df_nbinom, experienceContrast)
    
    sexSelectionEmmean.insert <- sexSelectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexSelectionEmmean <- cbind(geneID, sexSelectionEmmean.insert)
    sexSelectionPairEmmean_df_nbinom <- rbind(sexSelectionPairEmmean_df_nbinom, sexSelectionEmmean)
    
    sexSelection.insert <- sexSelectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexSelection <- cbind(geneID, sexSelection.insert)
    sexSelectionPair_df_nbinom <- rbind(sexSelectionPair_df_nbinom, sexSelection)
    
    expSelectionEmmean.insert <- expSelectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expSelectionEmmean <- cbind(geneID, expSelectionEmmean.insert)
    expSelectionPairEmmean_df_nbinom <- rbind(expSelectionPairEmmean_df_nbinom, expSelectionEmmean)
    
    expSelection.insert <- expSelectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expSelection <- cbind(geneID, expSelection.insert)
    expSelectionPair_df_nbinom <- rbind(expSelectionPair_df_nbinom, expSelection)
    
    selectionSexPairEmmean.insert <- selectionSexEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionSexPairEmmean <- cbind(geneID, selectionSexPairEmmean.insert)
    selectionSexPairEmmean_df_nbinom <- rbind(selectionSexPairEmmean_df_nbinom, selectionSexPairEmmean)
    
    selectionSexPair.insert <- selectionSexEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionSexPair <- cbind(geneID, selectionSexPair.insert)
    selectionSexPair_df_nbinom <- rbind(selectionSexPair_df_nbinom, selectionSexPair)
  }
}

save(lost_genes_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/lost_genes_nbinom.Rdata")

# Super saving
save(selectionPairEmmean_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/selectionPairEmmean_df_nbinom.Rdata")
save(experiencePairEmmean_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/experiencePairEmmean_df_nbinom.Rdata")
save(sexPairEmmean_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/sexPairEmmean_df_nbinom.Rdata")
save(sexContrastEmmean_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/sexContrastEmmean_df_nbinom.Rdata")
save(selectionContrastEmmean_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/selectionContrastEmmean_df_nbinom.Rdata")
save(experienceContrastEmmean_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/experienceContrastEmmean_df_nbinom.Rdata")
save(sexSelectionPairEmmean_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/sexSelectionPairEmmean_df_nbinom.Rdata")
save(expSelectionPairEmmean_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/expSelectionPairEmmean_df_nbinom.Rdata")
save(selectionSexPairEmmean_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/selectionSexPairEmmean_df_nbinom.Rdata")

save(selectionPair_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/selectionPair_df_nbinom.Rdata")
save(experiencePair_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/experiencePair_df_nbinom.Rdata")
save(sexPair_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/sexPair_df_nbinom.Rdata")
save(sexContrast_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/sexContrast_df_nbinom.Rdata")
save(selectionContrast_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/selectionContrast_df_nbinom.Rdata")
save(experienceContrast_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/experienceContrast_df_nbinom.Rdata")
save(sexSelectionPair_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/sexSelectionPair_df_nbinom.Rdata")
save(expSelectionPair_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/expSelectionPair_df_nbinom.Rdata")
save(selectionSexPair_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/selectionSexPair_df_nbinom.Rdata")

save(sex_anova_df_nbinom, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/sex_anova_df_nbinom.Rdata")
save(exp_anova_df_nbinom, file =  "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/exp_anova_df_nbinom.Rdata")
save(sel_anova_df_nbinom, file =  "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/sel_anova_df_nbinom.Rdata")
save(sexExp_anova_df_nbinom, file =  "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/sexExp_anova_df_nbinom.Rdata")
save(expSel_anova_df_nbinom, file =  "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/expSel_anova_df_nbinom.Rdata")
save(sexSel_anova_df_nbinom, file =  "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/nbinom_downVup/sexSel_anova_df_nbinom.Rdata")


