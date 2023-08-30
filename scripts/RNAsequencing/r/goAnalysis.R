## Down V up go analysis


library(topGO)
library(Rgraphviz)

# Big cutie gene to GO mapping
gene_GO <- readMappings("/Users/arteen/Desktop/School/Scripts/SociabilityRNA/fly_to_GO.delim.txt")

# This feels like it should stay constant
gene_filter <- function(allScore){
  return(allScore < 0.05)
}



# Load anova
sel_anova_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/gauss_downVup/sel_anova_df.Rdata"))

# P adjust and subset
sel_anova_df$p.adj <- p.adjust(sel_anova_df$pval, method = 'BY')
selection_anova_subset <- subset(sel_anova_df, sel_anova_df$p.adj < 0.05)
dim(selection_anova_subset) # 327 genes

gene_set <- sel_anova_df$pval
names(gene_set) <- sel_anova_df$gene
gene_set <- gene_set[complete.cases(gene_set)]



# 
# 
# upVdown_subset
# 
# gene_set <- unique(selectionContrast_df$geneID)
# 
# DGE_interesting_genes <- upVdown_subset$geneID
# 
# geneList <- factor(as.integer(gene_set %in% DGE_interesting_genes))
# names(geneList) <- gene_set
# 


## Set up



all_genes <- new("topGOdata",
                 ontology = "BP", 
                 allGenes = gene_set,  
                 geneSel = gene_filter, 
                 nodeSize = 5, 
                 annotationFun = annFUN.gene2GO, 
                 gene2GO = gene_GO) 


## all 
all_fisher <- runTest(all_genes, algorithm = "weight01", statistic = "fisher")

allGO = usedGO(all_genes)

all.table <- GenTable(all_genes, p.value = all_fisher, orderBy = 'weightFisher',
                           topNodes = length(allGO))

#performing BH correction on our p values
p.adj = p.adjust(all.table$p.value, method = 'BH')

# create the file with all the statistics from GO analysis
all_table_final = cbind(all.table, p.adj)

all_table_final = all_table_final[order(all_table_final$p.value), ]

#get list of significant GO before multiple testing correction
results.table.p = all_table_final[which(all_table_final$p.value <= 0.05), ]

#get list of significant GO after multiple testing correction
results.table.bh = all_table_final[which(all_table_final$p.adj <= 0.05), ]

## The actual two lowest p values dont appear in this list because they are in the format of X e-y so it doesnt realize that these two are small values
# so im just going to rbind them on top

results.table.bh <- results.table.bh[order(results.table.bh$p.adj),]

final_total_table <- rbind(results.table.bh, results.table.p)

final_table <- head(final_total_table, n = 15)

write.table(final_total_table, file = "/Users/arteen/Desktop/School/Thesis_Writing/Figures_and_Tables/table_DVU_DGE_goTerms.txt",
            sep = ",", row.names = FALSE, quote = FALSE)

#

showSigOfNodes(all_genes, score(all_fisher),
               useInfo = "def", sigForAll = FALSE,
               firstSigNodes = 2)










