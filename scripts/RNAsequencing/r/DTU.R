# DTU with glmmTMB

library(DRIMSeq)
library(DEXSeq)
library(stageR)
library(GenomeInfoDb)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(tximport)
library(BiocParallel)
library(car)
library(emmeans)
library(glmmTMB)
library(stringr)
library(parallel)



#To the directory of counts
#setwd("/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data")
setwd("/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/input")
quant_files <- file.path("quants", list.files("quants"), "quant.sf")
#file.exists(quant_files)
samples <- quant_files
samples <- substring(samples, 8, 26) #grab just sample name from the file name
names(quant_files) <- samples

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

rna.design <- data.frame(sample_id = samples, 
                         #file = quant_files,
                         selection = selection,
                         sex = sex,
                         expMatched = expMatched,
                         lane = lane,
                         lineage = lineage)
rna.design$selection <- relevel(as.factor(rna.design$selection), ref = "down")


txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
k <- keys(txdb, keytype = "TXNAME")
txdb <- AnnotationDbi::select(x = txdb, keys = k, "GENEID", "TXNAME")

txi <- tximport(quant_files, type = "salmon", txOut = TRUE,
                countsFromAbundance = "scaledTPM")

txi2 <- tximport(quant_files, type = "salmon", txOut = TRUE,
                countsFromAbundance = "no")





cts <- txi2$counts

cts2 <- txi2$counts

cts2 <- t(cts2)

save(cts2, file = "/Users/arteen/Desktop/transcriptCounts_notScaled.Rdata")

cts_scaledTPM <- t(cts)

cts <- cts[rowSums(cts) > 0, ]

txdb_match <- txdb[txdb$TXNAME %in% rownames(cts), ]
cts_match <- cts[rownames(cts) %in% txdb_match$TXNAME, ]


txdb <- txdb_match[order(txdb_match$TXNAME), ] #Reorder txdb


cts <- cts_match[order(rownames(cts_match)), ] #Reorder counts

# Confirmation that everything is a-okay
all(rownames(cts) == txdb$TXNAME)


counts <- data.frame(gene_id = txdb$GENEID,
                     feature_id = txdb$TXNAME,
                     cts)


tmp <- t(counts)

save(tmp, file = "/Users/arteen/Desktop/transcriptCounts_notScaled_withGenes.Rdata")

dex_pre_filter <- dmDSdata(counts = counts, samples = rna.design)



dex_filtered_subset <- dmFilter(dex_pre_filter,
                                
                                min_samps_feature_prop = 20,
                                
                                
                                min_feature_prop = 0.05,
                                
                                
                                min_samps_gene_expr = 28,
                                
                                
                                min_gene_expr = 10)

crap1 <- dmFilter(dex_pre_filter)



crap2 <- plotData(dex_filtered_subset)
crap <- plotData(crap1)
crap$data

barplot(table(crap$data), names.arg = 2:10,
        xlab = "transcript number",
        ylab = "frequency")

barplot(table(crap2$data), names.arg = c(2:23,27,28,31,71),
        xlab = "transcript number",
        ylab = "frequency")

crap3 <- as.data.frame(crap2$data)
crap3$tt <- as.factor(crap3$tt) # 4761 genes and 12335 transcripts
postfilter <- ggplot(data = crap3, aes(x = tt)) + geom_bar(fill = "black") + theme_classic() + xlab("Number of transcripts per gene") + ylab("Number of Genes")

crap4 <- as.data.frame(crap$data)
crap4$tt <- as.factor(crap4$tt) # 6559 genes and 21143
prefilter <- ggplot(data = crap4, aes(x = tt)) + geom_bar(fill = "black") + theme_classic() + xlab("Number of transcripts per gene") + ylab("Number of Genes")

plot_grid(prefilter, postfilter, nrow = 2, ncol = 1, labels = c("A)", "B)"), label_size = 11)

sample.data_subset <- DRIMSeq::samples(dex_filtered_subset)

count.data_subset <- round(as.matrix(counts(dex_filtered_subset)[ , -c(1:2)]))

count.data_subset_ids <- counts(dex_filtered_subset)[, c(1:2)]

sample_ids <- colnames( counts(dex_filtered_subset))


dxd_subset <- DEXSeqDataSet(countData = count.data_subset,
                            sampleData = sample.data_subset,
                            design = ~sample + exon + selection:exon,
                            featureID = counts(dex_filtered_subset)$feature_id,
                            groupID = counts(dex_filtered_subset)$gene_id)

dxd_subset  <- DEXSeq::estimateSizeFactors(dxd_subset)
sizeFactors <- colData(dxd_subset)@listData[["sizeFactor"]][1:142] # size factors but not norm factors

sample_table <- as.data.frame(colData(dxd_subset))


gene_id <- unique(count.data_subset_ids$gene_id)

LRT_df <- data.frame(gene = character(),
                           pval = double())

anova_sex_df <- data.frame(gene = character(),
                      pval = double())

anova_selection_df <- data.frame(gene = character(),
                                pval = double())

anova_sexSelInt_df <- data.frame(gene = character(),
                                 pval = double())

emmeans_df <- data.frame(gene = character(),
                         transcript = character(),
                         contrast = character(),
                         estimate = double(),
                         SE = double(),
                         df = double(),
                         t.ratio = double(),
                         p.value = double())

sex_emmeans_df <- data.frame(gene = character(),
                         transcript = character(),
                         contrast = character(),
                         estimate = double(),
                         SE = double(),
                         df = double(),
                         t.ratio = double(),
                         p.value = double())

emmeans_plotting_df <- data.frame(gene = character(),
                                  sex = character(),
                                  transcript = character(),
                                  emmean = double(),
                                  SE = double(),
                                  df = double(),
                                  lower.CL = double(),
                                  upper.CL = double())

sex_emmeans_plotting_df <- data.frame(gene = character(),
                                  sex = character(),
                                  transcript = character(),
                                  emmean = double(),
                                  SE = double(),
                                  df = double(),
                                  lower.CL = double(),
                                  upper.CL = double())

reduced_gene <- c()

nt <- min(parallel::detectCores(), 8)

counting <- 0

for (g in gene_id) {
  
  counting <- counting + 1
  if (counting%%5 == 0){
    print(counting)
    print(Sys.time())
  }
  
  # set up which gene and which transcripts
  index <- which(count.data_subset_ids$gene_id == g)
  count_subset <- count.data_subset[ c(index), ]
  gene_transcript_info <- count.data_subset_ids[ c(index), ]
  transcript_num <- length(gene_transcript_info$feature_id)
  
  empty_counts <- rep(0, 142)
  
  # set up df for first transcript
  transcript1 <- gene_transcript_info[1, 2]
  sample_transcript_df <- sample.data_subset
  gene_toInsert <- rep(g, 142)
  transcript_toInsert <- c(rep(transcript1, 142))
  
  transcript_df <- data.frame(gene_id = gene_toInsert,
                              transcript = transcript_toInsert,
                              counts = empty_counts)
  
  sample_df <- cbind(transcript_df, sample_transcript_df)
  
  # Add in first transcript counts
  for (i in 1:nrow(sample_df)){
    if (sample_df$transcript[i] == transcript1){
      the_sample <- as.character(sample_df$sample_id[i])
      sample_df$counts[i] <- as.integer(count_subset[1, the_sample])
    }
  }
  
  # Do the same but for all other transcripts
  for (i in 2:transcript_num){
    
    next_transcript <- gene_transcript_info[i, 2]
    next_sample_transcript_df <- sample.data_subset
    next_transcript_toInsert <- c(rep(next_transcript, 142))
    
    next_transcript_df <- data.frame(gene_id = gene_toInsert,
                                     transcript = next_transcript_toInsert,
                                     counts = empty_counts)
    
    next_sample_df <- cbind(next_transcript_df, next_sample_transcript_df)
    
    for (j in 1:nrow(next_sample_df)){
      if (next_sample_df$transcript[j] == next_transcript){
        the_sample <- as.character(next_sample_df$sample_id[j])
        next_sample_df$counts[j] <- as.integer(count_subset[i, the_sample])
      }
    }
    
    sample_df <- rbind(sample_df, next_sample_df)
  }
  
  sample_df$transcript <- as.factor(sample_df$transcript) # set transcript as factor
  
  # set up size factors
  the_sizeFactors <- c(rep(sizeFactors, transcript_num))
  names(the_sizeFactors) <- 1:length(the_sizeFactors)
  
  # To deal with complete separation
  countsPlus <- c()
  for (k in 1:nrow(sample_df)) {
    countsPlus <- append(countsPlus, sample_df[k, ]$counts + 1)
  }
  
  sample_df <- cbind(sample_df, countsPlus)
  
  model_data <- sample_df # clean up name
  
  
  ## run our models
  
 # diag(sex + transcript) - 1191 genes fail
 # diag(sex) - 1086 genes fail
  
  model_full <- glmmTMB(countsPlus ~ transcript + transcript:sex +
                          transcript:selection + transcript:sex:selection +
                          diag(sex | selection:lineage) +
                          (1|sample_id),
                        family = nbinom2(),
                        offset = the_sizeFactors,
                        data = model_data,
                        control = glmmTMBControl(parallel = nt))
  
  model_null <- glmmTMB(countsPlus ~ transcript + transcript:sex +
                          diag(sex | selection:lineage) +
                          (1|sample_id) ,
                        family = nbinom2(),
                        offset = the_sizeFactors,
                        data = model_data,
                        control = glmmTMBControl(parallel = nt))
  
  
  
  

  if (anyNA(as.vector(summary(model_full)$coef$cond[,2:4])) == TRUE |
      anyNA(as.vector(summary(model_null)$coef$cond[,2:4])) == TRUE |
    anyNA(logLik(model_full)) == TRUE |
    anyNA(logLik(model_null)) == TRUE ) {
      
    # reduced model because we need to fit
    model_full <- glmmTMB(countsPlus ~ transcript + transcript:sex +
                            transcript:selection + transcript:sex:selection +
                            diag(1 | selection:lineage) +
                            (1|sample_id),
                          family = nbinom2(),
                          offset = the_sizeFactors,
                          data = model_data,
                          control = glmmTMBControl(parallel = nt))
    
    model_null <- glmmTMB(countsPlus ~ transcript + transcript:sex +
                            diag(1 | selection:lineage) +
                            (1|sample_id) ,
                          family = nbinom2(),
                          offset = the_sizeFactors,
                          data = model_data,
                          control = glmmTMBControl(parallel = nt))

    if (anyNA(as.vector(summary(model_full)$coef$cond[,2:4])) == TRUE |
        anyNA(as.vector(summary(model_null)$coef$cond[,2:4])) == TRUE |
        anyNA(logLik(model_full)) == TRUE |
        anyNA(logLik(model_null)) == TRUE ) {
      reduced_gene <- append(reduced_gene, g)
      next;
    }
  }
  
  LRT_out <- anova(model_full, model_null, method = "LRT")
  
  LRT_pval <- LRT_out[2, 8]
  
  LRT_insert <- data.frame(gene = g, pval = LRT_pval)
  LRT_df <- rbind(LRT_df, LRT_insert)
  
  anova_out <- car::Anova(model_full)
  sex_res <- anova_out[2,3]
  selection_res <- anova_out[3, 3]
  interaction_res <- anova_out[4, 3]

  sex_insert <- data.frame(gene = g, pval = sex_res)
  selection_insert <- data.frame(gene = g, pval = selection_res)
  interaction_insert <- data.frame(gene = g, pval = interaction_res)

  anova_sex_df <- rbind(anova_sex_df, sex_insert)
  anova_selection_df <- rbind(anova_selection_df, selection_insert)
  anova_sexSelInt_df <- rbind(anova_sexSelInt_df, interaction_insert)

  # for now keep this just selection, once we get it all working, I'll implement sex and interaction emmeans extraction stuff
  emmeans_insert <- suppressMessages(emmeans(model_full, ~ selection | transcript))
  emmeans_insert_df_tmp <- as.data.frame(pairs(emmeans_insert))
  gene <- rep(g, transcript_num)
  emmeans_insert_df <- cbind(gene, emmeans_insert_df_tmp)
  emmeans_df <- rbind(emmeans_df, emmeans_insert_df)
  
  sex_emmeans_insert <- suppressMessages(emmeans(model_full, ~ sex | transcript))
  sex_emmeans_insert_df_tmp <- as.data.frame(pairs(sex_emmeans_insert))
  gene <- rep(g, transcript_num)
  sex_emmeans_insert_df <- cbind(gene, sex_emmeans_insert_df_tmp)
  sex_emmeans_df <- rbind(sex_emmeans_df, sex_emmeans_insert_df)
  
  sel_transcript_emmean <- suppressMessages(emmeans(model_full, ~  selection, type = "response"))
  plotting_insert <- as.data.frame(sel_transcript_emmean)
  plotting_insert <- as.data.frame(plotting_insert)
  plotting_insert <- cbind(gene, plotting_insert)
  
  emmeans_plotting_df <- rbind(emmeans_plotting_df, plotting_insert)
  
  sex_transcript_emmean <- suppressMessages(emmeans(model_full, ~  sex, type = "response"))
  sex_plotting_insert <- as.data.frame(sex_transcript_emmean)
  sex_plotting_insert <- as.data.frame(sex_plotting_insert)
  sex_plotting_insert <- cbind(gene, sex_plotting_insert)
  
  sex_emmeans_plotting_df <- rbind(sex_emmeans_plotting_df, sex_plotting_insert)
}

save(reduced_gene, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/dtu/fixed_2023/reduced_gene.Rdata")

save(LRT_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/dtu/fixed_2023/LRT_df.Rdata")

save(anova_sex_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/dtu/fixed_2023/anova_sex_df.Rdata")

save(anova_selection_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/dtu/fixed_2023/anova_selection_df.Rdata")

save(anova_sexSelInt_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/dtu/fixed_2023/anova_sexSelInt_df.Rdata")

save(emmeans_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/dtu/fixed_2023/emmeans_df.Rdata")

save(sex_emmeans_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/dtu/fixed_2023/sex_emmeans_df.Rdata")

save(emmeans_plotting_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/dtu/fixed_2023/emmeans_plotting_df.Rdata")

save(sex_emmeans_plotting_df, file = "/home/arteen/projects/def-idworkin/arteen/SociabilityRNA/r/output/dtu/fixed_2023/sex_emmeans_plotting_df.Rdata")





