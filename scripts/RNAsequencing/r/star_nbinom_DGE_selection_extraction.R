# DGE Selection Extraction STAR - Negative Binomial

# Load anova
nbinom_star_sel_anova_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/star_nbinom_downVup/sel_anova_df_nbinom_star.Rdata"))

# P adjust and subset
nbinom_star_sel_anova_df$p.adj <- p.adjust(nbinom_star_sel_anova_df$pval, method = 'BY')
nbinom_star_selection_anova_subset <- subset(nbinom_star_sel_anova_df, nbinom_star_sel_anova_df$p.adj < 0.05)
dim(nbinom_star_selection_anova_subset) # 319 genes



nbinom_star_selectionContrast_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/star_nbinom_downVup/selectionContrast_df_nbinom_star.Rdata"))

nbinom_star_all_selection_genes <- nbinom_star_selectionContrast_df[nbinom_star_selectionContrast_df$gene %in% nbinom_star_selection_anova_subset$gene, ]



# down - up
nbinom_star_upVdown <- subset(nbinom_star_all_selection_genes, nbinom_star_all_selection_genes$contrast == "down - up")

nbinom_star_upVdown_subset <- subset(nbinom_star_upVdown, nbinom_star_upVdown$p.value < 0.05)
dim(nbinom_star_upVdown_subset) # 154 genes



# down - control
nbinom_star_downVcontrol <- subset(nbinom_star_all_selection_genes, nbinom_star_all_selection_genes$contrast == "down - control")

nbinom_star_downVcontrol_subset <- subset(nbinom_star_downVcontrol, nbinom_star_downVcontrol$p.value < 0.05)
dim(nbinom_star_downVcontrol_subset) # 246 genes

# control - up
nbinom_star_controlVup <- subset(nbinom_star_all_selection_genes, nbinom_star_all_selection_genes$contrast == "control - up")

nbinom_star_controlVup_subset <- subset(nbinom_star_controlVup, nbinom_star_controlVup$p.value < 0.05)
dim(nbinom_star_controlVup_subset) # 199 genes



# reorder nbinom_star_all_selection_genes by fold change and grab those genes in order

nbinom_star_ordered_genes <- nbinom_star_all_selection_genes[order(abs(nbinom_star_all_selection_genes$estimate), decreasing = TRUE),]
head(nbinom_star_ordered_genes)

length(unique(nbinom_star_ordered_genes$geneID))

