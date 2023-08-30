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
tx2gene <- AnnotationDbi::select(x = txdb, keys = k, "GENEID", "TXNAME")

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
selectionPairEmmean_df <- data.frame()
experiencePairEmmean_df <- data.frame()
sexPairEmmean_df <- data.frame()
sexContrastEmmean_df <- data.frame()
selectionContrastEmmean_df <- data.frame()
experienceContrastEmmean_df <- data.frame()
sexSelectionPairEmmean_df <- data.frame()
expSelectionPairEmmean_df <- data.frame()
selectionSexPairEmmean_df <- data.frame()

selectionPair_df <- data.frame()
experiencePair_df <- data.frame()
sexPair_df <- data.frame()
sexContrast_df <- data.frame()
selectionContrast_df <- data.frame()
experienceContrast_df <- data.frame()
sexSelectionPair_df <- data.frame()
expSelectionPair_df <- data.frame()
selectionSexPair_df <- data.frame()

sex_anova_df <- data.frame(gene = character(),
                           pval = double(),
                           fullModel = character())

exp_anova_df <- data.frame(gene = character(),
                           pval = double(),
                           fullModel = character())

sel_anova_df <- data.frame(gene = character(),
                           pval = double(),
                           fullModel = character())

sexExp_anova_df <- data.frame(gene = character(),
                              pval = double(),
                              fullModel = character())

expSel_anova_df <- data.frame(gene = character(),
                              pval = double(),
                              fullModel = character())

sexSel_anova_df <- data.frame(gene = character(),
                              pval = double(),
                              fullModel = character())

lost_genes <- list()

## NOTE
# experiencePair_df will double save the contrasts (because emmeans and contrasts are 12 & 6 so will not match)
# so it will save the contrasts twice, remember to remove the duplicates



count <- 1
Sys.time()

for (g in targets) {
  
  count <- count + 1
  if (count%%100 == 0){
    print(count)
    print(Sys.time())
  }
  
  single_cpm <- voomCPM[g,]
  
  
  
  
  gaussian.model <- glmmTMB(single_cpm ~ lane +
                              sex + 
                              expMatched + 
                              selection + 
                              sex:expMatched + 
                              selection:expMatched + 
                              selection:sex + 
                              diag(0 + sex + expMatched | selection:lineage),
                            data = rna.design,
                            control = glmmTMBControl(parallel = nt))

  
 
  
  
  if (anyNA(as.vector(summary(gaussian.model)$coef$cond[,2:4])) == T | 
      anyNA(logLik(gaussian.model)) == T) {
    ## Using reduced model
    gaussian.model.reduced <- glmmTMB(single_cpm ~ lane +
                                        sex + 
                                        expMatched + 
                                        selection + 
                                        sex:expMatched + 
                                        selection:expMatched + 
                                        selection:sex + 
                                        diag(0 + sex | selection:lineage),
                                      data = rna.design,
                                      control = glmmTMBControl(parallel = nt))
   
    
    
    cpm.genes.check[[g]] <- list(summary(gaussian.model.reduced)$coef$cond, "Used Reduced Model")
    
    if (anyNA(as.vector(summary(gaussian.model.reduced)$coef$cond[,2:4])) == T | 
        anyNA(logLik(gaussian.model.reduced)) == T) {
     
      
      
      
      gaussian.model.reduced <- glmmTMB(single_cpm ~ lane +
                                          sex + 
                                          expMatched + 
                                          selection + 
                                          sex:expMatched + 
                                          selection:expMatched + 
                                          selection:sex + 
                                          diag(1 + sex | selection:lineage),
                                        data = rna.design,
                                        control = glmmTMBControl(parallel = nt))
       
      if (anyNA(as.vector(summary(gaussian.model.reduced)$coef$cond[,2:4])) == T | 
          anyNA(logLik(gaussian.model.reduced)) == T) { 
        
        lost_genes <- append(lost_genes, g)
        next;
        }
      
    }
    
   
    
    geneID <- g
    
    anova_out <- car::Anova(gaussian.model.reduced)
    
    sex_anova_res <- anova_out[2,3]
    sex_anova_insert <- data.frame(gene = geneID,
                                   pval = sex_anova_res,
                                   fullModel = "No")
    sex_anova_df <- rbind(sex_anova_df, sex_anova_insert)
    
    
    exp_anova_res <- anova_out[3,3]
    exp_anova_insert <- data.frame(gene = geneID,
                                   pval = exp_anova_res,
                                   fullModel = "No")
    exp_anova_df <- rbind(exp_anova_df, exp_anova_insert)
    
    
    sel_anova_res <- anova_out[4,3]
    sel_anova_insert <- data.frame(gene = geneID,
                                   pval = sel_anova_res,
                                   fullModel = "No")
    sel_anova_df <- rbind(sel_anova_df, sel_anova_insert)
    
    sexExp_anova_res <- anova_out[5,3]
    sexExp_anova_insert <- data.frame(gene = geneID,
                                      pval = sexExp_anova_res,
                                      fullModel = "No")
    sexExp_anova_df <- rbind(sexExp_anova_df, sexExp_anova_insert)
    
    expSel_anova_res <- anova_out[6,3]
    expSel_anova_insert <- data.frame(gene = geneID,
                                      pval = expSel_anova_res,
                                      fullModel = "No")
    expSel_anova_df <- rbind(expSel_anova_df, expSel_anova_insert)
    
    sexSel_anova_res <- anova_out[7,3]
    sexSel_anova_insert <- data.frame(gene = geneID,
                                      pval = sexSel_anova_res,
                                      fullModel = "No")
    sexSel_anova_df <- rbind(sexSel_anova_df, sexSel_anova_insert)
    
    
    
    #emmeans contrasts
    selectionPair <- emmeans(gaussian.model.reduced, pairwise ~ selection | sex + expMatched)
    experiencePair <- emmeans(gaussian.model.reduced, pairwise ~ expMatched | sex + selection)
    sexPairEmm <- emmeans(gaussian.model.reduced, pairwise ~ sex | expMatched + selection)
    sexEmm <- emmeans(gaussian.model.reduced, pairwise ~ sex)
    selectionEmm <- emmeans(gaussian.model.reduced, pairwise ~ selection)
    experienceEmm <- emmeans(gaussian.model.reduced, pairwise ~ expMatched)
    sexSelectionPair <- emmeans(gaussian.model.reduced, pairwise ~ sex | selection)
    expSelectionPair <- emmeans(gaussian.model.reduced, pairwise ~ expMatched | selection)
    selectionSexEmm <- emmeans(gaussian.model.reduced, pairwise ~ selection | sex)
    
    
    
    ## one for emmean one for contrast
    
    selectPairEmmean.insert <- selectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectPairEmmean <- cbind(geneID, selectPairEmmean.insert)
    selectionPairEmmean_df <- rbind(selectionPairEmmean_df, selectPairEmmean)
    
    selectPair.insert <- selectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectPair <- cbind(geneID, selectPair.insert)
    selectionPair_df <- rbind(selectionPair_df, selectPair)
    
    expPairEmmean.insert <- experiencePair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expPairEmmean <- cbind(geneID, expPairEmmean.insert)
    experiencePairEmmean_df <- rbind(experiencePairEmmean_df, expPairEmmean)
    
    expPair.insert <- experiencePair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expPair <- cbind(geneID, expPair.insert)
    experiencePair_df <- rbind(experiencePair_df, expPair)
    
    sexPairEmmean.insert <- sexPairEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexPairEmmean <- cbind(geneID, sexPairEmmean.insert)
    sexPairEmmean_df <- rbind(sexPairEmmean_df, sexPairEmmean)
    
    sexPair.insert <- sexPairEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexPair <- cbind(geneID, sexPair.insert)
    sexPair_df <- rbind(sexPair_df, sexPair)
    
    sexContrastEmmean.insert <- sexEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexContrastEmmean <- cbind(geneID, sexContrastEmmean.insert)
    sexContrastEmmean_df <- rbind(sexContrastEmmean_df, sexContrastEmmean)
    
    sexContrast.insert <- sexEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexContrast <- cbind(geneID, sexContrast.insert)
    sexContrast_df <- rbind(sexContrast_df, sexContrast)
    
    selectionContrastEmmean.insert <- selectionEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionContrastEmmean <- cbind(geneID, selectionContrastEmmean.insert)
    selectionContrastEmmean_df <- rbind(selectionContrastEmmean_df, selectionContrastEmmean)
    
    selectionContrast.insert <- selectionEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionContrast <- cbind(geneID, selectionContrast.insert)
    selectionContrast_df <- rbind(selectionContrast_df, selectionContrast)
    
    experienceContrastEmmean.insert <- experienceEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    experienceContrastEmmean <- cbind(geneID, experienceContrastEmmean.insert)
    experienceContrastEmmean_df <- rbind(experienceContrastEmmean_df, experienceContrastEmmean)
    
    experienceContrast.insert <- experienceEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    experienceContrast <- cbind(geneID, experienceContrast.insert)
    experienceContrast_df <- rbind(experienceContrast_df, experienceContrast)
    
    sexSelectionEmmean.insert <- sexSelectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexSelectionEmmean <- cbind(geneID, sexSelectionEmmean.insert)
    sexSelectionPairEmmean_df <- rbind(sexSelectionPairEmmean_df, sexSelectionEmmean)
    
    sexSelection.insert <- sexSelectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexSelection <- cbind(geneID, sexSelection.insert)
    sexSelectionPair_df <- rbind(sexSelectionPair_df, sexSelection)
    
    expSelectionEmmean.insert <- expSelectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expSelectionEmmean <- cbind(geneID, expSelectionEmmean.insert)
    expSelectionPairEmmean_df <- rbind(expSelectionPairEmmean_df, expSelectionEmmean)
    
    expSelection.insert <- expSelectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expSelection <- cbind(geneID, expSelection.insert)
    expSelectionPair_df <- rbind(expSelectionPair_df, expSelection)
    
    selectionSexPairEmmean.insert <- selectionSexEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionSexPairEmmean <- cbind(geneID, selectionSexPairEmmean.insert)
    selectionSexPairEmmean_df <- rbind(selectionSexPairEmmean_df, selectionSexPairEmmean)
    
    selectionSexPair.insert <- selectionSexEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionSexPair <- cbind(geneID, selectionSexPair.insert)
    selectionSexPair_df <- rbind(selectionSexPair_df, selectionSexPair)
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
    sex_anova_df <- rbind(sex_anova_df, sex_anova_insert)
    
    
    exp_anova_res <- anova_out[3,3]
    exp_anova_insert <- data.frame(gene = geneID,
                                   pval = exp_anova_res,
                                   fullModel = "Yes")
    exp_anova_df <- rbind(exp_anova_df, exp_anova_insert)
    
    
    sel_anova_res <- anova_out[4,3]
    sel_anova_insert <- data.frame(gene = geneID,
                                   pval = sel_anova_res,
                                   fullModel = "Yes")
    sel_anova_df <- rbind(sel_anova_df, sel_anova_insert)
    
    sexExp_anova_res <- anova_out[5,3]
    sexExp_anova_insert <- data.frame(gene = geneID,
                                   pval = sexExp_anova_res,
                                   fullModel = "Yes")
    sexExp_anova_df <- rbind(sexExp_anova_df, sexExp_anova_insert)
    
    expSel_anova_res <- anova_out[6,3]
    expSel_anova_insert <- data.frame(gene = geneID,
                                   pval = expSel_anova_res,
                                   fullModel = "Yes")
    expSel_anova_df <- rbind(expSel_anova_df, expSel_anova_insert)
    
    sexSel_anova_res <- anova_out[7,3]
    sexSel_anova_insert <- data.frame(gene = geneID,
                                   pval = sexSel_anova_res,
                                   fullModel = "Yes")
    sexSel_anova_df <- rbind(sexSel_anova_df, sexSel_anova_insert)
    
    
    
    
    
    #emmeans contrasts
    selectionPair <- emmeans(gaussian.model, pairwise ~ selection | sex + expMatched, type = "response")
    experiencePair <- emmeans(gaussian.model, pairwise ~ expMatched | sex + selection, type = "response")
    sexPairEmm <- emmeans(gaussian.model, pairwise ~ sex | expMatched + selection, type = "response")
    sexEmm <- emmeans(gaussian.model, pairwise ~ sex, type = "response")
    selectionEmm <- emmeans(gaussian.model, pairwise ~ selection, type = "response")
    experienceEmm <- emmeans(gaussian.model, pairwise ~ expMatched, type = "response")
    sexSelectionPair <- emmeans(gaussian.model, pairwise ~ sex | selection, type = "response")
    expSelectionPair <- emmeans(gaussian.model, pairwise ~ expMatched | selection, type = "response")
    selectionSexEmm <- emmeans(gaussian.model, pairwise ~ selection | sex, type = "response")
    
    
    
    ## one for emmean one for contrast
    
    selectPairEmmean.insert <- selectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectPairEmmean <- cbind(geneID, selectPairEmmean.insert)
    selectionPairEmmean_df <- rbind(selectionPairEmmean_df, selectPairEmmean)
    
    selectPair.insert <- selectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectPair <- cbind(geneID, selectPair.insert)
    selectionPair_df <- rbind(selectionPair_df, selectPair)
    
    expPairEmmean.insert <- experiencePair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expPairEmmean <- cbind(geneID, expPairEmmean.insert)
    experiencePairEmmean_df <- rbind(experiencePairEmmean_df, expPairEmmean)
    
    expPair.insert <- experiencePair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expPair <- cbind(geneID, expPair.insert)
    experiencePair_df <- rbind(experiencePair_df, expPair)
    
    sexPairEmmean.insert <- sexPairEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexPairEmmean <- cbind(geneID, sexPairEmmean.insert)
    sexPairEmmean_df <- rbind(sexPairEmmean_df, sexPairEmmean)
    
    sexPair.insert <- sexPairEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexPair <- cbind(geneID, sexPair.insert)
    sexPair_df <- rbind(sexPair_df, sexPair)
    
    sexContrastEmmean.insert <- sexEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexContrastEmmean <- cbind(geneID, sexContrastEmmean.insert)
    sexContrastEmmean_df <- rbind(sexContrastEmmean_df, sexContrastEmmean)
    
    sexContrast.insert <- sexEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexContrast <- cbind(geneID, sexContrast.insert)
    sexContrast_df <- rbind(sexContrast_df, sexContrast)
    
    selectionContrastEmmean.insert <- selectionEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionContrastEmmean <- cbind(geneID, selectionContrastEmmean.insert)
    selectionContrastEmmean_df <- rbind(selectionContrastEmmean_df, selectionContrastEmmean)
    
    selectionContrast.insert <- selectionEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionContrast <- cbind(geneID, selectionContrast.insert)
    selectionContrast_df <- rbind(selectionContrast_df, selectionContrast)
    
    experienceContrastEmmean.insert <- experienceEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    experienceContrastEmmean <- cbind(geneID, experienceContrastEmmean.insert)
    experienceContrastEmmean_df <- rbind(experienceContrastEmmean_df, experienceContrastEmmean)
    
    experienceContrast.insert <- experienceEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    experienceContrast <- cbind(geneID, experienceContrast.insert)
    experienceContrast_df <- rbind(experienceContrast_df, experienceContrast)
    
    sexSelectionEmmean.insert <- sexSelectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexSelectionEmmean <- cbind(geneID, sexSelectionEmmean.insert)
    sexSelectionPairEmmean_df <- rbind(sexSelectionPairEmmean_df, sexSelectionEmmean)
    
    sexSelection.insert <- sexSelectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    sexSelection <- cbind(geneID, sexSelection.insert)
    sexSelectionPair_df <- rbind(sexSelectionPair_df, sexSelection)
    
    expSelectionEmmean.insert <- expSelectionPair$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expSelectionEmmean <- cbind(geneID, expSelectionEmmean.insert)
    expSelectionPairEmmean_df <- rbind(expSelectionPairEmmean_df, expSelectionEmmean)
    
    expSelection.insert <- expSelectionPair$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    expSelection <- cbind(geneID, expSelection.insert)
    expSelectionPair_df <- rbind(expSelectionPair_df, expSelection)
    
    selectionSexPairEmmean.insert <- selectionSexEmm$emmeans %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionSexPairEmmean <- cbind(geneID, selectionSexPairEmmean.insert)
    selectionSexPairEmmean_df <- rbind(selectionSexPairEmmean_df, selectionSexPairEmmean)
    
    selectionSexPair.insert <- selectionSexEmm$contrasts %>%
      summary(infer = TRUE) %>%
      as.data.frame()
    selectionSexPair <- cbind(geneID, selectionSexPair.insert)
    selectionSexPair_df <- rbind(selectionSexPair_df, selectionSexPair)
  }
}

save(lost_genes, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/lost_genes.Rdata")

save(cpm.genes.check, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/cpm.genes.check.Rdata")
save(selectionPairEmmean_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/selectionPairEmmean_df.Rdata")
save(experiencePairEmmean_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/experiencePairEmmean_df.Rdata")
save(sexPairEmmean_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/sexPairEmmean_df.Rdata")
save(sexContrastEmmean_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/sexContrastEmmean_df.Rdata")
save(selectionContrastEmmean_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/selectionContrastEmmean_df.Rdata")
save(experienceContrastEmmean_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/experienceContrastEmmean_df.Rdata")
save(sexSelectionPairEmmean_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/sexSelectionPairEmmean_df.Rdata")
save(expSelectionPairEmmean_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/expSelectionPairEmmean_df.Rdata")
save(selectionSexPairEmmean_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/selectionSexPairEmmean_df.Rdata")

save(selectionPair_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/selectionPair_df.Rdata")
save(experiencePair_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/experiencePair_df.Rdata")
save(sexPair_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/sexPair_df.Rdata")
save(sexContrast_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/sexContrast_df.Rdata")
save(selectionContrast_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/selectionContrast_df.Rdata")
save(experienceContrast_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/experienceContrast_df.Rdata")
save(sexSelectionPair_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/sexSelectionPair_df.Rdata")
save(expSelectionPair_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/expSelectionPair_df.Rdata")
save(selectionSexPair_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/selectionSexPair_df.Rdata")

save(sex_anova_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/sex_anova_df.Rdata")
save(exp_anova_df, file =  "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/exp_anova_df.Rdata")
save(sel_anova_df, file =  "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/sel_anova_df.Rdata")
save(sexExp_anova_df, file =  "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/sexExp_anova_df.Rdata")
save(expSel_anova_df, file =  "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/expSel_anova_df.Rdata")
save(sexSel_anova_df, file =  "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/emmeans/fixed_2023/gauss_downVup/sexSel_anova_df.Rdata")



