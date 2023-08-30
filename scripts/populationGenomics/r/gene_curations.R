library(UpSetR)
`%notin%` <- Negate(`%in%`)



initial_list <- read.csv(  "/Users/arteen/Desktop/School/Projects/SociabilityRNA/SpreadsheetsAndData/geneCuration/geneCuration.csv")


rel_pheno <- subset(initial_list, initial_list$Relavant.phenotype. == "Yes" & initial_list$Real.Not..Yes.No.Maybe. == "Yes")

dim(rel_pheno) ## 33


rel_pheno_emmeans <- upVdown_subset[upVdown_subset$geneID %in% rel_pheno$FBgnID, ]


rel_pheno_emmeans <- rel_pheno_emmeans[order(abs(rel_pheno_emmeans$estimate), decreasing = TRUE),]

rel_pheno_subset_toPlot <- head(rel_pheno_emmeans$geneID, n = 12)




## Want to see overlap of genes between DGE and DTU lists - of total genes

DGE_list <- unique(all_selection_genes$geneID)

DTU_list <- unique(all_selection_genes_dtu$gene)



overlapping_genes <- intersect(DGE_list, DTU_list) 
length(overlapping_genes) # 39 genes 

library(UpSetR)

tmp <- upset(data = fromList(c(DGE_list, DTU_list)), nsets = 2)


compareList <- list(
  "DGE" = DGE_list,
  "DTU" = DTU_list
)




upset(fromList(compareList),
      nsets = 2,
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



# now in the low vs high only



DGE_list_upVdown <- upVdown_subset$geneID

DTU_list_upVdown <- unique(upVdown_dtu_subset$gene)


overlapping_genes_upVdown <- intersect(DGE_list_upVdown, DTU_list_upVdown) 
length(overlapping_genes_upVdown) # 17 genes

# now in the low vs control only



DGE_list_downVcontrol <- downVcontrol_subset$geneID

DTU_list_downVcontrol <- unique(downVcontrol_dtu_subset$gene)


overlapping_genes_downVcontrol <- intersect(DGE_list_downVcontrol, DTU_list_downVcontrol) 
length(overlapping_genes_downVcontrol) # 25 genes

# now in the control vs high only


DGE_list_controlVup <- controlVup_subset$geneID

DTU_list_controlVup <- unique(controlVup_dtu_subset$gene)


overlapping_genes_controlVup <- intersect(DGE_list_controlVup, DTU_list_controlVup)
length(overlapping_genes_controlVup)  # 21 genes

## Also curation of genes from other studies

# Wang et al. 2022

wang_gene_curation <- read.csv(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/SpreadsheetsAndData/geneCuration/WangEtAl_2022_geneComparison.csv", header = TRUE)


wang_gene_curation_list_notRemoved <- wang_gene_curation$FBgnID

wang_gene_curation_list <- subset(wang_gene_curation_list_notRemoved, wang_gene_curation_list_notRemoved != "NotFound")

wangEtAl_subset <- all_selection_genes[all_selection_genes$geneID %in% wang_gene_curation_list, ]

length(unique(wangEtAl_subset$geneID)) # 2

## 2 genes come up


# Shpigler et al. 2017

shpigler_csv <- read.csv(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/SpreadsheetsAndData/geneCuration/ShpiglerEtAl_2017_genes.csv", header = TRUE)

shp_gene_list_notRemoved <- shpigler_csv$FBgnID

shp_gene_list <- na.omit(shp_gene_list_notRemoved)


shp_gene_subset <- all_selection_genes[all_selection_genes$geneID %in% shp_gene_list, ]
length(unique(shp_gene_subset$geneID)) # 14 genes out of 1,057 lol




# Woodard et al. 2011

# WoodardEtAl_2011_genes - dataset is from their test 1


woodard_csv <- read.csv(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/SpreadsheetsAndData/geneCuration/WoodardEtAl_2011_genes.csv", header = TRUE)

woodard_gene_list_notRemoved <- woodard_csv$FBgnID

woodard_gene_list <- subset(woodard_gene_list_notRemoved, woodard_gene_list_notRemoved != "No fly ortholog")


woodard_gene_subset <- all_selection_genes[all_selection_genes$geneID %in% woodard_gene_list, ]
length(unique(woodard_gene_subset$geneID)) # 4 genes out of the 212 here


## Pop Gen

## Vs this listdownVcontrol_subset

pop_gen_CVD_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/controlVLow/snpEff_genes.txt",
                                header = TRUE) #336 genes

length(unique(pop_gen_CVD_genes$GeneId)) # 245 genes 

head(downVcontrol_subset)

popgen_overlap <- downVcontrol_subset[downVcontrol_subset$geneID %in% pop_gen_CVD_genes$GeneId,]

popgen_overlap <- pop_gen_CVD_genes[pop_gen_CVD_genes$GeneId %in% downVcontrol_subset$geneID,]

CVL_DGE_overlap <- data.frame(geneID = unique(popgen_overlap$GeneId),
                              geneName = unique(popgen_overlap$GeneName))

# write.table(CVL_DGE_overlap,
#             file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_DGE/CVL_geneList.csv",
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)

high_and_moderate_index <- c(grep(1, pop_gen_CVD_genes$variants_impact_HIGH),
                             grep(1, pop_gen_CVD_genes$variants_impact_MODERATE),
                             grep(2, pop_gen_CVD_genes$variants_impact_MODERATE),
                             grep(3, pop_gen_CVD_genes$variants_impact_MODERATE))

high_and_moderate <- unique(pop_gen_CVD_genes$GeneId[high_and_moderate_index])

high_and_moderate_overlap <- downVcontrol_subset[downVcontrol_subset$geneID %in% high_and_moderate,]


#### Now for low vs High


pop_gen_UVD_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/lowVHigh/gene_list_UVD_updated.txt",
                                header = TRUE)


# FBgn0262301 - this gene was causing trouble so manually adding it to the list here

## Vs this listdownVcontrol_subset


UVD_popGen_genes <- c("FBgn0262301", pop_gen_UVD_genes$GeneId)
UVD_popGen_geneNames <- c("mir-979", pop_gen_UVD_genes$GeneName)

head(upVdown_subset)

popgen_overlap_UVD <- upVdown_subset[upVdown_subset$geneID %in% UVD_popGen_genes,]

popgen_overlap_UVD <- pop_gen_UVD_genes[pop_gen_UVD_genes$GeneId %in% upVdown_subset$geneID,]

LVH_DGE_overlap <- data.frame(geneID = unique(popgen_overlap_UVD$GeneId),
                              geneName = unique(popgen_overlap_UVD$GeneName))
# 
# write.table(LVH_DGE_overlap,
#             file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_DGE/LVH_geneList.csv",
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)

high_and_moderate_index <- c(grep(1, pop_gen_UVD_genes$variants_impact_MODERATE))

high_and_moderate <- pop_gen_UVD_genes$GeneId[high_and_moderate_index]

high_and_moderate_overlap <- upVdown_subset[upVdown_subset$geneID %in% high_and_moderate,]

#### Now for Control vs High


pop_gen_CVU_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/controlVhigh/snpEff_genes.txt",
                                header = TRUE) 

length(unique(pop_gen_CVU_genes$GeneId)) # 184 genes 




popgen_overlap_CVU <- controlVup_subset[controlVup_subset$geneID %in% pop_gen_CVU_genes$GeneId,]


popgen_overlap_CVU <- pop_gen_CVU_genes[pop_gen_CVU_genes$GeneId %in% controlVup_subset$geneID,]

CVH_DGE_overlap <- data.frame(geneID = unique(popgen_overlap_CVU$GeneId),
                              geneName = unique(popgen_overlap_CVU$GeneName))

# write.table(CVH_DGE_overlap,
#             file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_DGE/CVH_geneList.csv",
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)

high_and_moderate_index <- c(grep(1, pop_gen_CVU_genes$variants_impact_HIGH),
                             grep(1, pop_gen_CVU_genes$variants_impact_MODERATE),
                             grep(2, pop_gen_CVU_genes$variants_impact_MODERATE),
                             grep(3, pop_gen_CVU_genes$variants_impact_MODERATE),
                             grep(5, pop_gen_CVU_genes$variants_impact_MODERATE))

high_and_moderate <- unique(pop_gen_CVU_genes$GeneId[high_and_moderate_index])

high_and_moderate_overlap <- controlVup_subset[controlVup_subset$geneID %in% high_and_moderate,]




## Ancestor vs Low
pop_gen_AVL_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/ancestorVlow/snpEff_genes.txt",
                                header = TRUE) 

length(unique(pop_gen_AVL_genes$GeneId)) # 184 genes 


downVcontrol_subset

popgen_overlap_CVU <- downVcontrol_subset[downVcontrol_subset$geneID %in% pop_gen_AVL_genes$GeneId,]


popgen_overlap_AVL <- pop_gen_AVL_genes[pop_gen_AVL_genes$GeneId %in% downVcontrol_subset$geneID,]

AVL_DGE_overlap <- data.frame(geneID = unique(popgen_overlap_AVL$GeneId),
                              geneName = unique(popgen_overlap_AVL$GeneName))

# write.table(AVL_DGE_overlap,
#             file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_DGE/AVL_geneList.csv",
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)



high_and_moderate_index <- c(grep(1, pop_gen_AVL_genes$variants_impact_HIGH),
                             grep(1, pop_gen_AVL_genes$variants_impact_MODERATE),
                             grep(2, pop_gen_AVL_genes$variants_impact_MODERATE),
                             grep(4, pop_gen_AVL_genes$variants_impact_MODERATE))

high_and_moderate <- unique(pop_gen_AVL_genes$GeneId[high_and_moderate_index])

high_and_moderate_overlap <- controlVup_subset[controlVup_subset$geneID %in% high_and_moderate,]


popgen_overlap_AVL_DTU <- downVcontrol_dtu_subset[downVcontrol_dtu_subset$gene %in% pop_gen_AVL_genes$GeneId,]



## Ancestor vs High
pop_gen_AVH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/ancestorVhigh/snpEff_genes.txt",
                                header = TRUE) 

length(unique(pop_gen_AVH_genes$GeneId)) # 184 genes 


downVcontrol_subset

popgen_overlap_CVU <- controlVup_subset[controlVup_subset$geneID %in% pop_gen_AVH_genes$GeneId,]

popgen_overlap_AVH <- pop_gen_AVH_genes[pop_gen_AVH_genes$GeneId %in% controlVup_subset$geneID,]

AVH_DGE_overlap <- data.frame(geneID = unique(popgen_overlap_AVH$GeneId),
                              geneName = unique(popgen_overlap_AVH$GeneName))

# write.table(AVH_DGE_overlap,
#             file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_DGE/AVH_geneList.csv",
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)




high_and_moderate_index <- c(grep(1, pop_gen_AVH_genes$variants_impact_MODERATE),
                             grep(2, pop_gen_AVH_genes$variants_impact_MODERATE),
                             grep(3, pop_gen_AVH_genes$variants_impact_MODERATE))

high_and_moderate <- unique(pop_gen_AVH_genes$GeneId[high_and_moderate_index])

high_and_moderate_overlap <- controlVup_subset[controlVup_subset$geneID %in% high_and_moderate,]


popgen_overlap_AVH_DTU <- controlVup_dtu_subset[controlVup_dtu_subset$gene %in% pop_gen_AVH_genes$GeneId,]






# pop gen genes overlapping with DTU genes

geneDictionary <- read.delim(file = "/Users/arteen/Downloads/FlyGeneDictionary.txt")
geneDictionary$FBgnID <- geneDictionary$validated_id

DTU_popgen_overlap_UVD <- upVdown_dtu_subset[upVdown_dtu_subset$transcript %in% pop_gen_UVD_genes$TranscriptId,]

DTU_popgen_overlap_CVD <- downVcontrol_dtu_subset[downVcontrol_dtu_subset$transcript %in% pop_gen_CVD_genes$TranscriptId,]

DTU_popgen_overlap_CVU <- controlVup_dtu_subset[controlVup_dtu_subset$transcript %in% pop_gen_CVU_genes$TranscriptId,]

DTU_popgen_overlap_AVL <- downVcontrol_dtu_subset[downVcontrol_dtu_subset$transcript %in% pop_gen_AVL_genes$TranscriptId,]

DTU_popgen_overlap_AVH <- controlVup_dtu_subset[controlVup_dtu_subset$transcript %in% pop_gen_AVH_genes$TranscriptId,]



####
LVH_DTU_overlap <- data.frame(transcript = DTU_popgen_overlap_UVD$transcript,
                              geneID = DTU_popgen_overlap_UVD$gene)
gene_name <- list()
for (i in 1:nrow(LVH_DTU_overlap)) {
  gene <- LVH_DTU_overlap$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
LVH_DTU_overlap <- cbind(LVH_DTU_overlap, geneName)
write.table(LVH_DTU_overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_DTU/LVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####

####
CVL_DTU_overlap <- data.frame(transcript = DTU_popgen_overlap_CVD$transcript,
                              geneID = DTU_popgen_overlap_CVD$gene)
gene_name <- list()
for (i in 1:nrow(CVL_DTU_overlap)) {
  gene <- CVL_DTU_overlap$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
CVL_DTU_overlap <- cbind(CVL_DTU_overlap, geneName)
write.table(CVL_DTU_overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_DTU/CVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####

####
CVH_DTU_overlap <- data.frame(transcript = DTU_popgen_overlap_CVU$transcript,
                              geneID = DTU_popgen_overlap_CVU$gene)
gene_name <- list()
for (i in 1:nrow(CVH_DTU_overlap)) {
  gene <- CVH_DTU_overlap$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
CVH_DTU_overlap <- cbind(CVH_DTU_overlap, geneName)
write.table(CVH_DTU_overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_DTU/CVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####

####
AVL_DTU_overlap <- data.frame(transcript = DTU_popgen_overlap_AVL$transcript,
                              geneID = DTU_popgen_overlap_AVL$gene)
gene_name <- list()
for (i in 1:nrow(AVL_DTU_overlap)) {
  gene <- AVL_DTU_overlap$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
AVL_DTU_overlap <- cbind(AVL_DTU_overlap, geneName)
write.table(AVL_DTU_overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_DTU/AVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####

####
AVH_DTU_overlap <- data.frame(transcript = DTU_popgen_overlap_AVH$transcript,
                              geneID = DTU_popgen_overlap_AVH$gene)
gene_name <- list()
for (i in 1:nrow(AVH_DTU_overlap)) {
  gene <- AVH_DTU_overlap$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
AVH_DTU_overlap <- cbind(AVH_DTU_overlap, geneName)
write.table(AVH_DTU_overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_DTU/AVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####



length(unique(DTU_popgen_overlap_UVD$transcript))
length(unique(DTU_popgen_overlap_CVD$transcript))
length(unique(DTU_popgen_overlap_CVU$transcript))






## Upset plot for the three pop gen lists

to_list <- list("Low versus High" = unique(UVD_popGen_genes),
                "Control versus Low" = (unique(pop_gen_CVD_genes$GeneId)),
                "Control versus High" = unique(pop_gen_CVU_genes$GeneId),
                "Ancestor versus High" = unique(pop_gen_AVH_genes$GeneId),
                "Ancestor versus Low" = unique(pop_gen_AVL_genes$GeneId))



upset(data = fromList(to_list),
      nsets = 5,
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



############# Writing Gene Lists #############
## Now going to make csvs of gene lists
## going to include both gene name and FBgnID


## CMH First

AVL_CMH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/ancestorVlow/cmh/snpEff_genes.txt",
                                header = TRUE) 

AVL_CMH_geneList <- data.frame(geneID = unique(AVL_CMH_genes$GeneId),
                               geneName = unique(AVL_CMH_genes$GeneName))


write.table(AVL_CMH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/cmh_lists/AVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####


AVH_CMH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/ancestorVhigh/cmh/snpEff_genes.txt",
                            header = TRUE) 

AVH_CMH_geneList <- data.frame(geneID = unique(AVH_CMH_genes$GeneId),
                               geneName = unique(AVH_CMH_genes$GeneName))

head(AVH_CMH_geneList)
write.table(AVH_CMH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/cmh_lists/AVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####


LVH_CMH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/lowVHigh/cmh/snpEff_genes.txt",
                            header = TRUE) 

LVH_CMH_geneList <- data.frame(geneID = unique(LVH_CMH_genes$GeneId),
                               geneName = unique(LVH_CMH_genes$GeneName))

head(LVH_CMH_geneList)
write.table(LVH_CMH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/cmh_lists/LVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####


CVL_CMH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/controlVLow/cmh/snpEff_genes.txt",
                            header = TRUE) 

CVL_CMH_geneList <- data.frame(geneID = unique(CVL_CMH_genes$GeneId),
                               geneName = unique(CVL_CMH_genes$GeneName))

head(CVL_CMH_geneList)
write.table(CVL_CMH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/cmh_lists/CVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####


CVH_CMH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/controlVhigh/cmh/snpEff_genes.txt",
                            header = TRUE) 

CVH_CMH_geneList <- data.frame(geneID = unique(CVH_CMH_genes$GeneId),
                               geneName = unique(CVH_CMH_genes$GeneName))

head(CVH_CMH_geneList)
write.table(CVH_CMH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/cmh_lists/CVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####




AVL_FST_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/ancestorVlow/fst/snpEff_genes.txt",
                            header = TRUE) 

AVL_FST_geneList <- data.frame(geneID = unique(AVL_FST_genes$GeneId),
                               geneName = unique(AVL_FST_genes$GeneName))

head(AVL_FST_geneList)
write.table(AVL_FST_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/fst_lists/AVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####




AVH_FST_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/ancestorVhigh/test/snpEff_genes.txt",
                            header = TRUE) 

AVH_FST_geneList <- data.frame(geneID = unique(AVH_FST_genes$GeneId),
                               geneName = unique(AVH_FST_genes$GeneName))

head(AVH_FST_geneList)
write.table(AVH_FST_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/fst_lists/AVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####




LVH_FST_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/lowVHigh/fst/snpEff_genes.txt",
                            header = TRUE) 

LVH_FST_geneList <- data.frame(geneID = unique(LVH_FST_genes$GeneId),
                               geneName = unique(LVH_FST_genes$GeneName))

head(LVH_FST_geneList)
write.table(LVH_FST_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/fst_lists/LVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####




CVL_FST_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/controlVLow/fst/snpEff_genes.txt",
                            header = TRUE) 

CVL_FST_geneList <- data.frame(geneID = unique(CVL_FST_genes$GeneId),
                               geneName = unique(CVL_FST_genes$GeneName))

head(CVL_FST_geneList)
write.table(CVL_FST_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/fst_lists/CVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####




CVH_FST_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/bedtools/controlVhigh/fst/snpEff_genes.txt",
                            header = TRUE) 

CVH_FST_geneList <- data.frame(geneID = unique(CVH_FST_genes$GeneId),
                               geneName = unique(CVH_FST_genes$GeneName))

head(CVH_FST_geneList)
write.table(CVH_FST_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/fst_lists/CVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####


# Now for the overlap between FST and CMH


LVH_CMH_FST_Overlap <- data.frame(geneID = unique(UVD_popGen_genes),
                                  geneName = unique(UVD_popGen_geneNames))
head(LVH_CMH_FST_Overlap)
write.table(LVH_CMH_FST_Overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_cmh_fst_lists/LVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
##
AVL_CMH_FST_Overlap <- data.frame(geneID = unique(pop_gen_AVL_genes$GeneId),
                                  geneName = unique(pop_gen_AVL_genes$GeneName))
head(AVL_CMH_FST_Overlap)
write.table(AVL_CMH_FST_Overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_cmh_fst_lists/AVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####
AVH_CMH_FST_Overlap <- data.frame(geneID = unique(pop_gen_AVH_genes$GeneId),
                                  geneName = unique(pop_gen_AVH_genes$GeneName))
head(AVH_CMH_FST_Overlap)
write.table(AVH_CMH_FST_Overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_cmh_fst_lists/AVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

####
CVL_CMH_FST_Overlap <- data.frame(geneID = unique(pop_gen_CVD_genes$GeneId),
                                  geneName = unique(pop_gen_CVD_genes$GeneName))
head(CVL_CMH_FST_Overlap)
write.table(CVL_CMH_FST_Overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_cmh_fst_lists/CVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####
CVH_CMH_FST_Overlap <- data.frame(geneID = unique(pop_gen_CVU_genes$GeneId),
                                  geneName = unique(pop_gen_CVU_genes$GeneName))
head(CVH_CMH_FST_Overlap)
write.table(CVH_CMH_FST_Overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/overlap_cmh_fst_lists/CVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)



# OLAP w DGE




