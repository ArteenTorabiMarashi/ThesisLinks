library(ggplot2)
library(cowplot)
get_sync <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/actual_sync/unmerged_subsetted.sync", header = FALSE)
headers <- c("chrom","position","ref",
             "Ancestor1","Ancestor2",
             "Ancestor3","Ancestor4",
             "Control1","Control2",
             "Control3","Control4",
             "Down1", "Down2",
             "Down3","Down4",
             "Up1","Up2",
             "Up3","Up4")
colnames(get_sync) <- headers




firstSNP <- subset(get_sync, get_sync$chrom == "2L" & get_sync$position == "2265057")

first_Ancestor <- c(7/83, 7/90, 7/84, 8/93)

ancestor_first <- mean(first_Ancestor)
SE_first_ancestor <- sd(first_Ancestor) / sqrt(length(first_Ancestor)) 

first_control1 <- c(0/102)
first_control2 <- c(39/94)
first_control3 <- c(41/111)
first_control4 <- c(15/95)

first_low1 <- c(0/136)
first_low2 <- c(17/95)
first_low3 <- c(37/110)
first_low4 <- c(0/96)

first_high1 <- c(11/118)
first_high2 <- c(8/91)
first_high3 <- c(13/112)
first_high4 <- c(57/100)

first_alleleF <- c(ancestor_first,
                 first_control1, first_control2, first_control3, first_control4,
                 first_low1, first_low2, first_low3, first_low4,
                 first_high1, first_high2, first_high3, first_high4)



treatment <- c("Ancestor", rep("Control", 4), rep("Low", 4), rep("High", 4))
SE_list <- c(SE_first_ancestor, rep(NA, 12))

first_af_df <- data.frame(alleleFreq = first_alleleF, SE = SE_list, treatment = treatment)

first_af_df$treatment <- factor(first_af_df$treatment, levels = c("Ancestor", "Control", "Low", "High"))



pd = position_dodge(width = 0.3)
first <- ggplot(data = first_af_df, aes(x = treatment, y = alleleFreq)) +
  geom_jitter(aes(color = treatment), position = pd , size = 2.5) +
  geom_errorbar(data = first_af_df, 
                aes(x = treatment, y = alleleFreq, ymin = alleleFreq-SE, ymax = alleleFreq+SE, color = treatment),
                width = 0.15, size = 0.7, position = pd, alpha = 1) +
  theme_classic() +
  theme( plot.title = element_text(), axis.text = element_text(face = "bold")) +
  xlab("Treatment") +
  ylab("Allele Frequency") +
  ggtitle("chr2L:2265057, T>A, Missense: Asn>Lys")
#CG34049|FBgn0054049|transcript|FBtr0302372|protein_coding|1/3|c.255T>A|p.Asn85Lys











secondSNP <- subset(get_sync, get_sync$chrom == "3L" & get_sync$position == "13482719")
stv_Ancestor <- c(78/97, 83/104, 84/103, 78/97)

ancestor_stv <- mean(stv_Ancestor)
SE_stv_ancestor <- sd(stv_Ancestor) / sqrt(length(stv_Ancestor)) 

stv_control1 <- c(93/116)
stv_control2 <- c(56/101)
stv_control3 <- c( 86/124)
stv_control4 <- c(131/131)

stv_low1 <- c(53/114)
stv_low2 <- c(80/116)
stv_low3 <- c(82/134)
stv_low4 <- c(64/120)

stv_high1 <- c(130/144)
stv_high2 <- c(93/130)
stv_high3 <- c(91/103)
stv_high4 <- c( 72/107)


stv_alleleF <- c(ancestor_stv,
                 stv_control1, stv_control2, stv_control3, stv_control4,
                 stv_low1, stv_low2, stv_low3, stv_low4,
                 stv_high1, stv_high2, stv_high3, stv_high4)



treatment <- c("Ancestor", rep("Control", 4), rep("Low", 4), rep("High", 4))
SE_list <- c(SE_stv_ancestor, rep(NA, 12))

stv_af_df <- data.frame(alleleFreq = stv_alleleF, SE = SE_list, treatment = treatment)

stv_af_df$treatment <- factor(stv_af_df$treatment, levels = c("Ancestor", "Control", "Low", "High"))



pd = position_dodge(width = 0.3)
stv <- ggplot(data = stv_af_df, aes(x = treatment, y = alleleFreq)) +
  geom_jitter(aes(color = treatment), position = pd , size = 2.5) +
  geom_errorbar(data = stv_af_df, 
                aes(x = treatment, y = alleleFreq, ymin = alleleFreq-SE, ymax = alleleFreq+SE, color = treatment),
                width = 0.15, size = 0.7, position = pd, alpha = 1) +
  theme_classic() +
  theme(axis.text = element_text(face = "bold")) +
  xlab("Treatment") +
  ylab("Allele Frequency") +
  ggtitle("chr3L:13482719, T>C, Missense: Val>Ala")



plot_grid(first, stv,
          nrow = 1, ncol = 2,
          labels = c('A)','B)'))


