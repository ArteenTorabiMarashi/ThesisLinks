library(ggplot2)
library(cowplot)

fst_insync_Subsetted <- read.csv(file = "/Users/arteen/Desktop/pop_gen_out/fst/sync_subsetted/fst.csv") # Mpileup

# avc = 5, CVD = 8, CVU = 9, DVU = 10

## Mpileup
head(fst_insync_Subsetted, n = 4)
subsetted_AVC_fst_raw <- fst_insync_Subsetted[, c(1:4,7)]

head(subsetted_AVC_fst_raw, n = 4)

subsetted_fst_plotting <- subsetted_AVC_fst_raw[, c(1,2, 5)]


subsetted_fst_plotting_number <- na.omit(subsetted_fst_plotting)
subsetted_fst_plotting_number <- chrNumbering(subsetted_fst_plotting_number)


names(subsetted_fst_plotting_number)[3] <-  "mean"

numbering_chr <- middleChr(subsetted_fst_plotting_number)

DownVsUp <- ggplot(data = subsetted_fst_plotting_number, aes(x = number, y = mean, color = chrom)) + 
  geom_point(size = 0.65, show.legend = F, alpha = 0.2) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab(bquote(F[ST]))+
  scale_x_discrete(limits = c(numbering_chr), labels = c("X", "2L", "2R", '3L', '3R', '4')) +
  geom_smooth(aes(group = as.factor(chrom)), colour = "red", size = 0.5) +
  #ylim(-1, 1) +
  scale_colour_manual(values = c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))


AncestorVsControl
ControlVsDown
ControlVsUp
DownVsUp



sync_sub



plot_grid(mp, mp_sync,
          vcf_sync_scaled, vcf_sync,
          nrow = 2, ncol = 2)


## What comparisons
# AVC 
#CVD UVC DVU
# 4 plots - mp, mp to sync, vcf to sync matched, vcf to sync




## Export top 5 % of Fst
headers <- c("chrom", "chromStart", "chromEnd") # for bed file -ness

export_cleaned <- na.omit(fst_insync_Subsetted)
# AVL
AVL_start <- export_cleaned[ ,c(1:3, 6) ]
AVL_subsetted <- AVL_start[AVL_start$Ancestor.Low > quantile(AVL_start$Ancestor.Low, 0.95), ]
AVL_top5 <- AVL_subsetted[ , c(1:3)]
colnames(AVL_top5) <- headers



write.table(x = AVL_top5, file = "/Users/arteen/Desktop/pop_gen_out/bedtools/ancestorVlow/fst_top.bed",col.names = F, sep = '\t', row.names = FALSE, quote = FALSE)
# AVH
AVH_start <- export_cleaned[ ,c(1:3, 7) ]
AVH_subsetted <- UVD_start[AVH_start$Ancestor.High > quantile(AVH_start$Ancestor.High, 0.95), ]
AVH_top5 <- UVD_subsetted[ , c(1:3)]
colnames(AVH_top5) <- headers



write.table(x = AVH_top5, file = "/Users/arteen/Desktop/pop_gen_out/bedtools/ancestorVhigh/fst_top.bed",col.names = F, sep = '\t', row.names = FALSE, quote = FALSE)





AVC_start <- export_cleaned[ ,c(1:3, 5) ]
AVC_subsetted <- AVC_start[AVC_start$Ancestor.Control > quantile(AVC_start$Ancestor.Control, 0.95), ]
AVC_top5 <- AVC_subsetted[ , c(1:3)]
colnames(AVC_top5) <- headers
AVC_top5$startPos <- paste(AVC_top5$chrom, AVC_top5$chromStart, sep = "_")


CVD_start <- export_cleaned[ ,c(1:3, 8) ]
CVD_subsetted <- CVD_start[CVD_start$Control.Low > quantile(CVD_start$Control.Low, 0.95), ]
CVD_top5 <- CVD_subsetted[ , c(1:3)]
colnames(CVD_top5) <- headers
CVD_top5$startPos <- paste(CVD_top5$chrom, CVD_top5$chromStart, sep = "_")


write.table(x = CVD_top5, file = "/Users/arteen/Desktop/pop_gen_out/bedtools/fst_top.bed",col.names = F, sep = '\t', row.names = FALSE, quote = FALSE)

CVU_start <- export_cleaned[ ,c(1:3, 9) ]
CVU_subsetted <- CVU_start[CVU_start$Control.High > quantile(CVU_start$Control.High, 0.95), ]
CVU_top5 <- CVU_subsetted[ , c(1:3)]
colnames(CVU_top5) <- headers
CVU_top5$startPos <- paste(CVU_top5$chrom, CVU_top5$chromStart, sep = "_")


UVD_start <- export_cleaned[ ,c(1:3, 10) ]
UVD_subsetted <- UVD_start[UVD_start$Low.High > quantile(UVD_start$Low.High, 0.95), ]
UVD_top5 <- UVD_subsetted[ , c(1:3)]
colnames(UVD_top5) <- headers



write.table(x = UVD_top5, file = "/Users/arteen/Desktop/pop_gen_out/bedtools/lowVHigh/fst_top.bed",col.names = F, sep = '\t', row.names = FALSE, quote = FALSE)



 ## subsetting

`%notin%` <- Negate(`%in%`)

CVD_top5_subsetted <- CVD_top5[(CVD_top5$startPos %notin% AVC_top5$startPos), ]
write.table(x = CVD_top5_subsetted, file = "/Users/arteen/Desktop/pop_gen_out/bedtools/CVD_subsetted_fst_top.bed",col.names = F, sep = '\t', row.names = FALSE, quote = FALSE)

CVU_top5_subsetted <- CVU_top5[(CVU_top5$startPos %notin% AVC_top5$startPos), ]
write.table(x = CVU_top5_subsetted, file = "/Users/arteen/Desktop/pop_gen_out/bedtools/controlVhigh/CVU_FST.bed",col.names = F, sep = '\t', row.names = FALSE, quote = FALSE)



subsetted_fst_plotting_number$startPos <- paste(subsetted_fst_plotting_number$chrom, subsetted_fst_plotting_number$start, sep = "_")

to_plot <- subsetted_fst_plotting_number[ subsetted_fst_plotting_number$startPos %notin% AVC_top5$startPos, ]

ggplot(data = to_plot, aes(x = number, y = mean, color = chrom)) + 
  geom_point(size = 0.65, show.legend = F, alpha = 0.2) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab(bquote(F[ST])) +
  ylim(-0.05, 0.35) +
  scale_x_discrete(limits = c(numbering_chr), labels = c("X", "2L", "2R", '3L', '3R', '4')) +
  geom_smooth(aes(group = as.factor(chrom)), colour = "red", size = 0.5) +
  scale_colour_manual(values = c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

# With Drift positions
ggplot(data = subsetted_fst_plotting_number, aes(x = number, y = mean, color = chrom)) + 
  geom_point(size = 0.65, show.legend = F, alpha = 0.2) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab(bquote(F[ST]))+
  ylim(-0.05, 0.35) +
  scale_x_discrete(limits = c(numbering_chr), labels = c("X", "2L", "2R", '3L', '3R', '4')) +
  geom_smooth(aes(group = as.factor(chrom)), colour = "red", size = 0.5) +
  scale_colour_manual(values = c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))




