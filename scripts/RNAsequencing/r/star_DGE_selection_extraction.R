# DGE Selection Extraction STAR - Gaussian

# Load anova
star_sel_anova_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/star_gauss_downVup/sel_anova_df_star.Rdata"))

# P adjust and subset
star_sel_anova_df$p.adj <- p.adjust(star_sel_anova_df$pval, method = 'BY')
star_selection_anova_subset <- subset(star_sel_anova_df, star_sel_anova_df$p.adj < 0.05)
dim(star_selection_anova_subset) # 345 genes



star_selectionContrast_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/star_gauss_downVup/selectionContrast_df_star.Rdata"))

star_all_selection_genes <- star_selectionContrast_df[star_selectionContrast_df$gene %in% star_selection_anova_subset$gene, ]



# down - up
star_upVdown <- subset(star_all_selection_genes, star_all_selection_genes$contrast == "down - up")

star_upVdown_subset <- subset(star_upVdown, star_upVdown$p.value < 0.05)
dim(star_upVdown_subset) # 183 genes



# down - control
star_downVcontrol <- subset(star_all_selection_genes, star_all_selection_genes$contrast == "down - control")

star_downVcontrol_subset <- subset(star_downVcontrol, star_downVcontrol$p.value < 0.05)
dim(star_downVcontrol_subset) # 282 genes

# control - up
star_controlVup <- subset(star_all_selection_genes, star_all_selection_genes$contrast == "control - up")

star_controlVup_subset <- subset(star_controlVup, star_controlVup$p.value < 0.05)
dim(star_controlVup_subset) # 210 genes



# reorder star_all_selection_genes by fold change and grab those genes in order

star_ordered_genes <- star_all_selection_genes[order(abs(star_all_selection_genes$estimate), decreasing = TRUE),]
head(star_ordered_genes)

length(unique(star_ordered_genes$geneID))








