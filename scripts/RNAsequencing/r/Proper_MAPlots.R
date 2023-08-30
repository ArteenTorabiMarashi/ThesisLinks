## Full plotting of all the contrast MA plots can be found in SocPaperFigures_2023.Rmd
# Just doing experience here as a show of how to get this figure done
#experience


## Exp dataframe
exp_anova_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/gauss_downVup/exp_anova_df.Rdata"))



exp_anova_df$p.adj <- p.adjust(exp_anova_df$pval, method = "BY")

exp_anova_subset <- subset(exp_anova_df, exp_anova_df$p.adj < 0.05)
dim(exp_anova_subset) # 213 genes

expContrast_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/gauss_downVup/experienceContrast_df.Rdata"))

all_exp_genes <- expContrast_df[expContrast_df$geneID %in% exp_anova_subset$gene, ]



head(all_exp_genes) # the 213 genes
head(expContrast_df) # difference
# mean
exp_contrastEmmean_df <- get(load(file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Data/gauss_downVup/experienceContrastEmmean_df.Rdata"))
head(exp_contrastEmmean_df) # to get mean

# Getting mean

arenaEmmean <- subset(exp_contrastEmmean_df, exp_contrastEmmean_df$expMatched == "socArena")
vialEmmean <- subset(exp_contrastEmmean_df, exp_contrastEmmean_df$expMatched == "vial")
meanExperience <- arenaEmmean # just to get column names nice and quick, and we can disregard the expMatched col

meanExperience$emmean <- (arenaEmmean$emmean + vialEmmean$emmean) / 2





mergeFrame <- merge(meanExperience, expContrast_df, by = c("geneID"))
xAxis <- mergeFrame$emmean
yAxis <- mergeFrame$estimate
mergeFrame <- cbind(mergeFrame, yAxis)
mergeFrame <- cbind(mergeFrame, xAxis)

sig_points <- mergeFrame[mergeFrame$geneID %in% all_exp_genes$gene, ]
`%notin%` <- Negate(`%in%`) # creates a little function to make a not in pipe
regular_points <- mergeFrame[mergeFrame$geneID %notin% all_exp_genes$gene, ]

arenaVsVial_MAPlot <-  ggplot() + 
  geom_point(data = regular_points, aes(xAxis, yAxis, color = "blue"), alpha = .1) +
  geom_point(data = sig_points, aes(xAxis, yAxis, color = "red"), alpha = 0.7) +
  # geom_point(aes(color = col, alpha = alpha) ) +
  xlab(expression("Mean Average (" * log[2] * CPM * ")")) +
  ylab(expression("Mean Difference (" * log[2] * CPM * ")")) +
  scale_color_identity() + 
  theme_classic() +
  theme(legend.position="none", plot.title = element_text(face = "bold"))
  
arenaVsVial_MAPlot

