## DTU extraction 

library(UpSetR)




lrt_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/dtu_output/lrt_df.Rdata"))
anova_sex_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/dtu_output/anova_sex_df.Rdata"))
anova_selection_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/dtu_output/anova_selection_df.Rdata"))
anova_sexSel_int_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/dtu_output/anova_sexSelInt_df.Rdata"))
emmeans_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/dtu_output/emmeans_df.Rdata"))
emmeans_plotting_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/dtu_output/emmeans_plotting_df.Rdata"))
reduced_gene <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/dtu_output/reduced_gene.Rdata"))
sex_emmeans_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/dtu_output/sex_emmeans_df.Rdata"))
sex_emmeans_plotting_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/dtu_output/sex_emmeans_plotting_df.Rdata"))

length(reduced_gene)

## Sex

anova_sex_df$p.adj <- p.adjust(anova_sex_df$pval, method = "BY")


anova_sex_subset_dtu <- subset(anova_sex_df, anova_sex_df$p.adj < 0.05)
length(unique(anova_sex_subset_dtu$gene)) # 2631

all_sex_genes_dtu <- sex_emmeans_df[sex_emmeans_df$gene %in% anova_sex_subset_dtu$gene, ]

## Selection
anova_selection_df$p.adj <- p.adjust(anova_selection_df$pval, method = "BY")

anova_selection_subset_dtu <- subset(anova_selection_df, anova_selection_df$p.adj < 0.05)
length(unique(anova_selection_subset_dtu$gene))  # 619

all_selection_genes_dtu <- emmeans_df[emmeans_df$gene %in% anova_selection_subset_dtu$gene, ]


lrt_df$p.adj <- p.adjust(lrt_df$pval, method = "BH")

lrt_subset <- subset(lrt_df, lrt_df$p.adj < 0.05)
dim(lrt_subset) # 662 genes 
head(lrt_df)




upVdown_dtu <- subset(all_selection_genes_dtu, all_selection_genes_dtu$contrast == "down - up")

upVdown_dtu_subset <- subset(upVdown_dtu, upVdown_dtu$p.value < 0.05)
length(unique(upVdown_dtu_subset$gene)) # 190 genes


downVcontrol_dtu <- subset(all_selection_genes_dtu, all_selection_genes_dtu$contrast == "down - control")

downVcontrol_dtu_subset <- subset(downVcontrol_dtu, downVcontrol_dtu$p.value < 0.05)
length(unique(downVcontrol_dtu_subset$gene)) # 384 genes


controlVup_dtu <- subset(all_selection_genes_dtu, all_selection_genes_dtu$contrast == "control - up")

controlVup_dtu_subset <- subset(controlVup_dtu, controlVup_dtu$p.value < 0.05)
length(unique(controlVup_dtu_subset$gene)) # 252 genes


geneList1 <- list("Low vs High" = unique(upVdown_dtu_subset$gene),
                 "Low vs Control" = unique(downVcontrol_dtu_subset$gene),
                 "Control vs High" = unique(controlVup_dtu_subset$gene))

##

write.table(upVdown_dtu_subset,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/DTU/LVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
write.table(downVcontrol_dtu_subset,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/DTU/CVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
write.table(controlVup_dtu_subset,
            file = "/Users/arteen/Desktop/pop_gen_out/gene_lists/DTU/CVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####

geneList_upset1 <- upset(fromList(geneList1),
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


anova_sexSel_int_df$p.adj <- p.adjust(anova_sexSel_int_df$pval, method = 'BY')
dtu_interaction_subset <- subset(anova_sexSel_int_df, anova_sexSel_int_df$p.adj < 0.05)
length(unique(dtu_interaction_subset$gene)) # 14 genes here




selectionGeneList <- list(HighVsLow = upVdown_dtu_subset$gene,
                       LowVsControl = downVcontrol_dtu_subset$gene,
                       ControlVsHigh = controlVup_dtu_subset$gene)

selectionGeneList_upset <- upset(fromList(selectionGeneList),
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


selectionTranscriptList <- list(HighVsLow = upVdown_dtu_subset$transcript,
                          LowVsControl = downVcontrol_dtu_subset$transcript,
                          ControlVsHigh = controlVup_dtu_subset$transcript)

selectionTranscriptList_upset <- upset(fromList(selectionTranscriptList),
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


selectionGeneList_upset
selectionTranscriptList_upset



