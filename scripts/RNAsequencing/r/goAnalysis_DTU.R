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





anova_selection_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/dtu_output/anova_selection_df.Rdata"))
anova_selection_df$p.adj <- p.adjust(anova_selection_df$pval, method = 'BY')

gene_set <- anova_selection_df$pval
names(gene_set) <- anova_selection_df$gene
gene_set <- gene_set[complete.cases(gene_set)]






# 
# gene_set <- unique(emmeans_df$gene)
# 
# DTU_interesting_genes <- unique(upVdown_dtu_subset$gene)
# 
# geneList <- factor(as.integer(gene_set %in% DTU_interesting_genes))
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

final_table <- head(all_table_final, n = 15)

write.table(results.table.p, file = "/Users/arteen/Desktop/School/Thesis_Writing/Figures_and_Tables/GOtable_DVU_DTU_goTerms.txt",
            sep = ",", row.names = FALSE, quote = FALSE)


View(all_table_final)
#

showSigOfNodes(all_genes, score(all_fisher),
               useInfo = "def", sigForAll = FALSE,
               firstSigNodes = 2)










