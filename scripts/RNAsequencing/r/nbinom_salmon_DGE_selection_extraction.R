# DGE Selection Extraction Salmon - Negative Binomial

# Load anova
nbinom_sel_anova_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/nbinom_downVup/sel_anova_df_nbinom.Rdata"))

# P adjust and subset
nbinom_sel_anova_df$p.adj <- p.adjust(nbinom_sel_anova_df$pval, method = 'BY')
nbinom_selection_anova_subset <- subset(nbinom_sel_anova_df, nbinom_sel_anova_df$p.adj < 0.05)
dim(nbinom_selection_anova_subset) # 289 genes



nbinom_selectionContrast_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/nbinom_downVup/selectionContrast_df_nbinom.Rdata"))

nbinom_all_selection_genes <- nbinom_selectionContrast_df[nbinom_selectionContrast_df$gene %in% nbinom_selection_anova_subset$gene, ]



# down - up
nbinom_upVdown <- subset(nbinom_all_selection_genes, nbinom_all_selection_genes$contrast == "down - up")

nbinom_upVdown_subset <- subset(nbinom_upVdown, nbinom_upVdown$p.value < 0.05)
dim(nbinom_upVdown_subset) # 142 genes



# down - control
nbinom_downVcontrol <- subset(nbinom_all_selection_genes, nbinom_all_selection_genes$contrast == "down - control")

nbinom_downVcontrol_subset <- subset(nbinom_downVcontrol, nbinom_downVcontrol$p.value < 0.05)
dim(nbinom_downVcontrol_subset) # 226 genes

# control - up
nbinom_controlVup <- subset(nbinom_all_selection_genes, nbinom_all_selection_genes$contrast == "control - up")

nbinom_controlVup_subset <- subset(nbinom_controlVup, nbinom_controlVup$p.value < 0.05)
dim(nbinom_controlVup_subset) # 176 genes



# reorder nbinom_all_selection_genes by fold change and grab those genes in order

nbinom_ordered_genes <- nbinom_all_selection_genes[order(abs(nbinom_all_selection_genes$estimate), decreasing = TRUE),]
head(nbinom_ordered_genes)

length(unique(nbinom_ordered_genes$geneID))




