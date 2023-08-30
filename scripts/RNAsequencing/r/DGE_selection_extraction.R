# DGE Selection Extraction

## Exp dataframe
exp_anova_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/gauss_downVup/exp_anova_df.Rdata"))

## just doing exp right now

exp_anova_df$p.adj <- p.adjust(exp_anova_df$pval, method = "BY")

exp_anova_subset <- subset(exp_anova_df, exp_anova_df$p.adj < 0.05)
dim(exp_anova_subset) # 213 genes

expContrast_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/gauss_downVup/experienceContrast_df.Rdata"))

all_exp_genes <- expContrast_df[expContrast_df$geneID %in% exp_anova_subset$gene, ]

length(intersect(exp_anova_subset$gene, selection_anova_subset$gene)) # 3 overlapping between selection and experience
length(intersect(exp_anova_subset$gene, sex_anova_subset$gene)) # 107 overlapping between sex and experience


## Exp by selection dataframe
expSel_anova_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/gauss_downVup/expSel_anova_df.Rdata"))

## just doing exp right now

expSel_anova_df$p.adj <- p.adjust(expSel_anova_df$pval, method = "BY")

expSel_anova_subset <- subset(expSel_anova_df, expSel_anova_df$p.adj < 0.05)
dim(expSel_anova_subset)


# sex contrast
sexContrast_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/gauss_downVup/sexContrast_df.Rdata"))
sex_anova_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/gauss_downVup/sex_anova_df.Rdata"))

# P adjust and subset
sex_anova_df$p.adj <- p.adjust(sex_anova_df$pval, method = 'BY')
sex_anova_subset <- subset(sex_anova_df, sex_anova_df$p.adj < 0.05)
dim(sex_anova_subset) # 327 genes

# Load anova
sel_anova_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/gauss_downVup/sel_anova_df.Rdata"))

# P adjust and subset
sel_anova_df$p.adj <- p.adjust(sel_anova_df$pval, method = 'BY')
selection_anova_subset <- subset(sel_anova_df, sel_anova_df$p.adj < 0.05)
dim(selection_anova_subset) # 327 genes



selectionContrast_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/gauss_downVup/selectionContrast_df.Rdata"))

all_selection_genes <- selectionContrast_df[selectionContrast_df$geneID %in% selection_anova_subset$gene, ]



# down - up
upVdown <- subset(all_selection_genes, all_selection_genes$contrast == "down - up")

upVdown_subset <- subset(upVdown, upVdown$p.value < 0.05)
dim(upVdown_subset) # 174 genes

upVdownTop <- head(upVdown_subset[order(abs(upVdown_subset$estimate), decreasing = TRUE),] , n = 12)

# down - control
downVcontrol <- subset(all_selection_genes, all_selection_genes$contrast == "down - control")

downVcontrol_subset <- subset(downVcontrol, downVcontrol$p.value < 0.05)
dim(downVcontrol_subset) # 271 genes

# control - up
controlVup <- subset(all_selection_genes, all_selection_genes$contrast == "control - up")

controlVup_subset <- subset(controlVup, controlVup$p.value < 0.05)
dim(controlVup_subset) # 194 genes

## Writing gene List

####
geneDictionary <- read.delim(file = "/Users/arteen/Downloads/FlyGeneDictionary.txt")
geneDictionary$FBgnID <- geneDictionary$validated_id





gene_name <- list()
for (i in 1:nrow(upVdown_subset)) {
  gene <- upVdown_subset$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
LVH_geneList <- cbind(geneName, upVdown_subset)
write.table(LVH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/DGE/LVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####

gene_name <- list()
for (i in 1:nrow(downVcontrol_subset)) {
  gene <- downVcontrol_subset$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
CVL_geneList <- cbind(geneName, downVcontrol_subset)
head(CVL_geneList)
write.table(CVL_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/DGE/CVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####

gene_name <- list()
for (i in 1:nrow(controlVup_subset)) {
  gene <- controlVup_subset$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
CVH_geneList <- cbind(geneName, controlVup_subset)
head(CVH_geneList)
write.table(CVH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/DGE/CVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####
####

gene_name <- list()
for (i in 1:nrow(all_exp_genes)) {
  gene <- all_exp_genes$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
exp_geneList <- cbind(geneName, all_exp_genes)
head(exp_geneList)
write.table(exp_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/DGE/exp_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####






length(Reduce(intersect, list(upVdown_subset$geneID, downVcontrol_subset$geneID, controlVup_subset$geneID)))

length(Reduce(unique, list(upVdown_subset$geneID, downVcontrol_subset$geneID, controlVup_subset$geneID)))


tmp1 <- unique(controlVup_subset$geneID,downVcontrol_subset$geneID )

tmp2 <- unique(upVdown_subset$geneID , tmp1)

tmp <- list("Low Versus High" = upVdown_subset$geneID,
            "Control Versus High" = controlVup_subset$geneID,
            "Low versus Control" = downVcontrol_subset$geneID)

upset(data = fromList(tmp),
      nsets = 3,
      matrix.color = "red",
      main.bar.color = "black",
      mainbar.y.label = "Number of Genes",
      sets.bar.color = "Black",
      sets.x.label = "Total Set Size",
      point.size = 2.2,
      line.size = 0.7,
      order.by = c("freq"),
      decreasing = c(T),
      group.by = "degree")


# reorder all_selection_genes by fold change and grab those genes in order

ordered_genes <- all_selection_genes[order(abs(all_selection_genes$estimate), decreasing = TRUE), ]
head(ordered_genes)

length(unique(ordered_genes$geneID))


upVdown_ordered <- upVdown_subset[order(abs(upVdown_subset$estimate), decreasing = TRUE), ]
head(upVdown_ordered$geneID, n = 10)


