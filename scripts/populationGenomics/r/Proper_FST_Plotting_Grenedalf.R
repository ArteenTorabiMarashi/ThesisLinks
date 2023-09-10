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


## As revisions - going to split first the X chrom out and pull out top 5 and then for autosomes and then back together


export_cleaned <- na.omit(fst_insync_Subsetted)
# AVL


AVL_start <- export_cleaned[ ,c(1:3, 6) ]


AVL_start_X <- subset(AVL_start, AVL_start$chrom == "X")
AVL_subsetted_X <- AVL_start_X[AVL_start_X$Ancestor.Low > quantile(AVL_start_X$Ancestor.Low, 0.95), ]
AVL_start_autosomes <- subset(AVL_start, AVL_start$chrom != "X")
AVL_subsetted_autosomes <- AVL_start_autosomes[AVL_start_autosomes$Ancestor.Low > quantile(AVL_start_autosomes$Ancestor.Low, 0.95), ]


AVL_subsetted <- rbind(AVL_subsetted_autosomes, AVL_subsetted_X)
AVL_top5 <- AVL_subsetted[ , c(1:3)]
colnames(AVL_top5) <- headers



write.table(x = AVL_top5, file = "/Users/arteen/Desktop/pop_gen_out/revisions/ancestorVlow/fst_top.bed",col.names = F, sep = '\t', row.names = FALSE, quote = FALSE)
# AVH
AVH_start <- export_cleaned[ ,c(1:3, 7) ]


AVH_start_X <- subset(AVH_start, AVH_start$chrom == "X")
AVH_subsetted_X <- AVH_start_X[AVH_start_X$Ancestor.High > quantile(AVH_start_X$Ancestor.High, 0.95), ]
AVH_start_autosomes <- subset(AVH_start, AVH_start$chrom != "X")
AVH_subsetted_autosomes <- AVH_start_autosomes[AVH_start_autosomes$Ancestor.High > quantile(AVH_start_autosomes$Ancestor.High, 0.95), ]


AVH_subsetted <- rbind(AVH_subsetted_autosomes, AVH_subsetted_X)


AVH_top5 <- AVH_subsetted[ , c(1:3)]
colnames(AVH_top5) <- headers



write.table(x = AVH_top5, file = "/Users/arteen/Desktop/pop_gen_out/revisions/ancestorVhigh/fst_top.bed",col.names = F, sep = '\t', row.names = FALSE, quote = FALSE)





AVC_start <- export_cleaned[ ,c(1:3, 5) ]

AVC_start_X <- subset(AVC_start, AVC_start$chrom == "X")
AVC_start_autosomes <- subset(AVC_start, AVC_start$chrom != "X")

AVC_subsetted_X <- AVC_start_X[AVC_start_X$Ancestor.Control > quantile(AVC_start_X$Ancestor.Control, 0.95), ]
AVC_subsetted_autosomes <- AVC_start_autosomes[AVC_start_autosomes$Ancestor.Control > quantile(AVC_start_autosomes$Ancestor.Control, 0.95), ]
AVC_subsetted <- rbind(AVC_subsetted_autosomes, AVC_subsetted_X)
AVC_top5 <- AVC_subsetted[ , c(1:3)]
colnames(AVC_top5) <- headers
AVC_top5$startPos <- paste(AVC_top5$chrom, AVC_top5$chromStart, sep = "_")


CVD_start <- export_cleaned[ ,c(1:3, 8) ]
CVD_start_X <- subset(CVD_start, CVD_start$chrom == "X")
CVD_start_autosomes <- subset(CVD_start, CVD_start$chrom != "X")

CVD_subsetted_X <- CVD_start_X[CVD_start_X$Control.Low > quantile(CVD_start_X$Control.Low, 0.95), ]
CVD_subsetted_autosomes <- CVD_start_autosomes[CVD_start_autosomes$Control.Low > quantile(CVD_start_autosomes$Control.Low, 0.95), ]



CVD_subsetted <- rbind(CVD_subsetted_autosomes, CVD_subsetted_X)
CVD_top5 <- CVD_subsetted[ , c(1:3)]
colnames(CVD_top5) <- headers
CVD_top5$startPos <- paste(CVD_top5$chrom, CVD_top5$chromStart, sep = "_")


CVU_start <- export_cleaned[ ,c(1:3, 9) ]
CVU_start_X <- subset(CVU_start, CVU_start$chrom == "X")
CVU_start_autosomes <- subset(CVU_start, CVU_start$chrom != "X")
CVU_subsetted_X <- CVU_start_X[CVU_start_X$Control.High > quantile(CVU_start_X$Control.High, 0.95), ]
CVU_subsetted_autosomes <- CVU_start_autosomes[CVU_start_autosomes$Control.High > quantile(CVU_start_autosomes$Control.High, 0.95), ]

CVU_subsetted <- rbind(CVU_subsetted_autosomes, CVU_subsetted_X)
CVU_top5 <- CVU_subsetted[ , c(1:3)]
colnames(CVU_top5) <- headers
CVU_top5$startPos <- paste(CVU_top5$chrom, CVU_top5$chromStart, sep = "_")


UVD_start <- export_cleaned[ ,c(1:3, 10) ]
tmp <- UVD_start[UVD_start$Low.High > quantile(UVD_start$Low.High, 0.95), ]

UVD_start_X <- subset(UVD_start, UVD_start$chrom == "X")
UVD_start_autosomes <- subset(UVD_start, UVD_start$chrom != "X")
UVD_subsetted_X <- UVD_start_X[UVD_start_X$Low.High > quantile(UVD_start_X$Low.High, 0.95), ]
UVD_subsetted_autosomes <- UVD_start_autosomes[UVD_start_autosomes$Low.High > quantile(UVD_start_autosomes$Low.High, 0.95), ]

UVD_subsetted <- rbind(UVD_subsetted_autosomes, UVD_subsetted_X)
UVD_top5 <- UVD_subsetted[ , c(1:3)]
colnames(UVD_top5) <- headers



write.table(x = UVD_top5, file = "/Users/arteen/Desktop/pop_gen_out/revisions/lowVhigh/fst_top.bed",col.names = F, sep = '\t', row.names = FALSE, quote = FALSE)



 ## subsetting

`%notin%` <- Negate(`%in%`)

CVD_top5_subsetted <- CVD_top5[(CVD_top5$startPos %notin% AVC_top5$startPos), ]
write.table(x = CVD_top5_subsetted, file = "/Users/arteen/Desktop/pop_gen_out/revisions/lowVcontrol/fst_top.bed",col.names = F, sep = '\t', row.names = FALSE, quote = FALSE)

CVU_top5_subsetted <- CVU_top5[(CVU_top5$startPos %notin% AVC_top5$startPos), ]
write.table(x = CVU_top5_subsetted, file = "/Users/arteen/Desktop/pop_gen_out/revisions/controlVhigh/fst_top.bed",col.names = F, sep = '\t', row.names = FALSE, quote = FALSE)


