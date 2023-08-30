## Down V up go analysis


library(topGO)
library(Rgraphviz)

# Big cutie gene to GO mapping
gene_GO <- readMappings("/Users/arteen/Desktop/School/Scripts/SociabilityRNA/fly_to_GO.delim.txt")

# 103 terms

## This is control vs down
gene_set <- row.names(txi$counts)

interestingGenes <- unique(UVD_popGen_genes)
high_and_moderate_genes <- high_and_moderate

geneList <- factor(as.integer(gene_set %in% interestingGenes))
names(geneList) <- gene_set




## Set up



all_genes <- new("topGOdata",
                 ontology = "BP", 
                 allGenes = geneList,  
                 
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


final_table <- results.table.p

write.table(final_table, file = "/Users/arteen/Desktop/School/Thesis_Writing/Figures_and_Tables/GOTABLE_DVU_popGen.txt",
          sep = ",", row.names = FALSE, quote = FALSE)

#

showSigOfNodes(all_genes, score(all_fisher),
               useInfo = "def", sigForAll = FALSE,
               firstSigNodes = 2)










