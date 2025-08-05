library(ggplot2)
library(patchwork)
library(tidyverse)

df1 <- read.table("C:\\Users\\ricar\\Downloads\\2021.11.02 CADM2 vs GUSB for R.txt", header=T, sep="\t")
df1 <- df1%>% filter(genotype!="hom_rd1")

df1$genotype <-gsub("wt","WT",df1$genotype)
df1$genotype <-gsub("100kb_het","HET 100kb DEL",df1$genotype)
df1$genotype <-gsub("500bp_del","HET 500bp DEL",df1$genotype)
df1$genotype <-gsub("inv","HET 500bp INV",df1$genotype)
df1$genotype <-gsub("hom_rd2","HOM 500bp DEL",df1$genotype)
df1$genotype = factor(df1$genotype, levels=c("WT","HET 500bp DEL","HET 100kb DEL","HET 500bp INV","HOM 500bp DEL"))

df1$zygosity <- df1$genotype
df1$zygosity <- gsub("HET 500bp INV", "HET",df1$zygosity)
df1$zygosity <- gsub("HET 100kb DEL", "HET",df1$zygosity)
df1$zygosity <- gsub("HET 500bp DEL", "HET",df1$zygosity)
df1$zygosity <- gsub("HOM 500bp DEL", "HOM",df1$zygosity)
df1$zygosity = factor(df1$zygosity, levels=c("WT","HET","HOM"))

pdf("iN_CADM2_vs_GUSB_zyg.pdf", width = 6, height = 6)
ggplot(df1, aes(x=genotype, y=fold_change_vs_wt, fill=zygosity)) +
  geom_boxplot(fatten=NULL, position = 'dodge', outliers = F) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.7, size = 1, linetype = "solid") +
  ylab("Fold Change") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(family="Arial",size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

comparison = list(c("WT","HET 500bp DEL"),c("WT","HET 100kb DEL"),c("WT","HET 500bp INV"),c("WT","HOM 500bp DEL"))
# Mean +/- standard deviation
pdf("iN_CADM2_vs_GUSB_zygosity.pdf", width = 3, height = 6)
ggerrorplot(df1, x = "genotype", y = "fold_change_vs_wt", 
            desc_stat = "mean_sd",
            color = "zygosity",
            size = 1.4,
            legend="center",
            ylab="Fold Change vs WT",
            xlab="") +
  rotate_x_text()
dev.off()
  stat_compare_means(comparisons = comparison)
  


df2 <- read.table("C:\\Users\\ricar\\Downloads\\2021.11.03 AC vs GUSB 45 cyc.txt", header=T, sep="\t")
df2 <- df2%>% filter(genotype!="hom_rd1")

df2$genotype <-gsub("wt","WT",df2$genotype)
df2$genotype <-gsub("100kb_het","HET 100kb DEL",df2$genotype)
df2$genotype <-gsub("500bp_del","HET 500bp DEL",df2$genotype)
df2$genotype <-gsub("inv","HET 500bp INV",df2$genotype)
df2$genotype <-gsub("hom_rd2","HOM 500bp DEL",df2$genotype)
df2$genotype = factor(df2$genotype, levels=c("WT","HET 500bp DEL","HET 100kb DEL","HET 500bp INV","HOM 500bp DEL"))


df2$zygosity <- df2$genotype
df2$zygosity <- gsub("HET 500bp INV", "HET",df2$zygosity)
df2$zygosity <- gsub("HET 100kb DEL", "HET",df2$zygosity)
df2$zygosity <- gsub("HET 500bp DEL", "HET",df2$zygosity)
df2$zygosity <- gsub("HOM 500bp DEL", "HOM",df2$zygosity)
df2$zygosity = factor(df2$zygosity, levels=c("WT","HET","HOM"))

wt <- df2 %>% filter(zygosity=="WT")
het <- df2 %>% filter(zygosity=="HET")
hom <- df2 %>% filter(zygosity=="HOM")

t.test(het$fold_change_vs_wt,wt$fold_change_vs_wt)$p.value
mean(het$fold_change_vs_wt)/mean(wt$fold_change_vs_wt)

mean(hom$fold_change_vs_wt)/mean(wt$fold_change_vs_wt)
t.test(hom$fold_change_vs_wt,wt$fold_change_vs_wt)$p.value

pdf("AC_vs_GUSB.pdf", height = 6, width = 6)
ggplot(df2, aes(x=genotype, y=fold_change_vs_wt, fill=genotype)) +
  geom_boxplot(fatten=NULL, position = 'dodge', outliers = F) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.7, size = 1, linetype = "solid") +
  ylab("Fold Change") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("AC_vs_GUSB_zyg.pdf", height = 6, width = 3)
ggerrorplot(df2, x = "genotype", y = "fold_change_vs_wt", 
            desc_stat = "mean_sd",
            color = "zygosity",
            size = 1.4,
            legend="center",
            ylab="Fold Change vs WT",
            xlab="") +
  rotate_x_text()
dev.off()

df3 <- read.table("C:\\Users\\ricar\\Downloads\\2021.12.18 AC vs GUSB iPSC 45 cyc for R.txt", header=T, sep="\t")
df3 <- df3%>% filter(genotype!="hom_rd1")

df3$genotype <-gsub("wt","WT",df3$genotype)
df3$genotype <-gsub("100kb_het","HET 100kb DEL",df3$genotype)
df3$genotype <-gsub("500bp_del","HET 500bp DEL",df3$genotype)
df3$genotype <-gsub("inv","HET 500bp INV",df3$genotype)
df3$genotype <-gsub("hom_rd2","HOM 500bp DEL",df3$genotype)
df3$genotype = factor(df3$genotype, levels=c("WT","HET 500bp DEL","HET 100kb DEL","HET 500bp INV","HOM 500bp DEL"))


df3$zygosity <- df3$genotype
df3$zygosity <- gsub("HET 500bp INV", "HET",df3$zygosity)
df3$zygosity <- gsub("HET 100kb DEL", "HET",df3$zygosity)
df3$zygosity <- gsub("HET 500bp DEL", "HET",df3$zygosity)
df3$zygosity <- gsub("HOM 500bp DEL", "HOM",df3$zygosity)
df3$zygosity = factor(df3$zygosity, levels=c("WT","HET","HOM"))

pdf("iPSC_AC_vs_GUSB.pdf", height = 6, width = 6)
ggplot(df3, aes(x=genotype, y=fold_change_vs_wt, fill=genotype)) +
  geom_boxplot(fatten=NULL, position = 'dodge', outliers = F) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.7, size = 1, linetype = "solid") +
  ylab("Fold Change") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("iPSC_AC_vs_GUSB_zyg.pdf", height = 6, width = 3)
ggerrorplot(df3, x = "genotype", y = "fold_change_vs_wt", 
            desc_stat = "mean_sd",
            color = "zygosity",
            size = 1.4,
            legend="center",
            ylab="Fold Change vs WT",
            xlab="") +
  rotate_x_text()
dev.off()




df4 <- read.table("C:\\Users\\ricar\\Downloads\\2022.01.04 iPSC CADM2 qFR vs GUSB for R.txt", header=T, sep="\t")
df4 <- df4%>% filter(genotype!="hom_rd1")

df4$genotype <-gsub("wt","WT",df4$genotype)
df4$genotype <-gsub("100kb_het","HET 100kb DEL",df4$genotype)
df4$genotype <-gsub("500bp_del","HET 500bp DEL",df4$genotype)
df4$genotype <-gsub("inv","HET 500bp INV",df4$genotype)
df4$genotype <-gsub("hom_rd2","HOM 500bp DEL",df4$genotype)
df4$genotype = factor(df4$genotype, levels=c("WT","HET 500bp DEL","HET 100kb DEL","HET 500bp INV","HOM 500bp DEL"))

pdf("iPSC_CADM2_qFR_vs_GUSB.pdf", height = 6, width = 6)
ggplot(df4, aes(x=genotype, y=fold_change_vs_wt, fill=genotype)) +
  geom_boxplot(fatten=NULL, position = 'dodge', outliers = F) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.7, size = 1, linetype = "solid") +
  ylab("Fold Change Compared to WT") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



df4$zygosity <- df4$genotype
df4$zygosity <- gsub("HET 500bp INV", "HET",df4$zygosity)
df4$zygosity <- gsub("HET 100kb DEL", "HET",df4$zygosity)
df4$zygosity <- gsub("HET 500bp DEL", "HET",df4$zygosity)
df4$zygosity <- gsub("HOM 500bp DEL", "HOM",df4$zygosity)
df4$zygosity = factor(df4$zygosity, levels=c("WT","HET","HOM"))

pdf("iPSC_CADM2_qFR_vs_GUSB.pdf", height = 6, width = 3)
ggerrorplot(df4, x = "genotype", y = "fold_change_vs_wt", 
            desc_stat = "mean_sd",
            color = "genotype",
            size = 1.4,
            legend="center",
            ylab="Fold Change vs WT",
            xlab="") +
  rotate_x_text()
dev.off()

pdf("iPSC_CADM2_qFR_vs_GUSB_zyg.pdf", height = 6, width = 3)
ggerrorplot(df4, x = "genotype", y = "fold_change_vs_wt", 
            desc_stat = "mean_sd",
            color = "zygosity",
            size = 1.4,
            legend="center",
            ylab="Fold Change vs WT",
            xlab="") +
  rotate_x_text()
dev.off()



############Gene Expression
library(ggplot2)
library(patchwork)
library(tidyverse)

df1 <- read.table("C:\\Users\\ricar\\Downloads\\2021.11.02 CADM2 vs GUSB for R.txt", header=T, sep="\t")
df1 <- df1%>% filter(genotype!="hom_rd1")
df1$genotype <-gsub("wt","WT",df1$genotype)
df1$genotype <-gsub("100kb_het","HET 100kb DEL",df1$genotype)
df1$genotype <-gsub("500bp_del","HET 500bp DEL",df1$genotype)
df1$genotype <-gsub("inv","HET 500bp INV",df1$genotype)
df1$genotype <-gsub("hom_rd2","HOM 500bp DEL",df1$genotype)
df1$genotype = factor(df1$genotype, levels=c("WT","HET 500bp DEL","HET 100kb DEL","HET 500bp INV","HOM 500bp DEL"))

pdf("iN_CADM2_vs_GUSB.pdf", width = 6, height = 6)
ggplot(df1, aes(x=genotype, y=fold_change_vs_wt, fill=genotype)) +
  geom_boxplot(width=0.7, position = 'dodge', outliers = F) +
  ylab("Fold Change") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

df1$zygosity <- df1$genotype
df1$zygosity <- gsub("HET 500bp INV", "HET",df1$zygosity)
df1$zygosity <- gsub("HET 100kb DEL", "HET",df1$zygosity)
df1$zygosity <- gsub("HET 500bp DEL", "HET",df1$zygosity)
df1$zygosity <- gsub("HOM 500bp DEL", "HOM",df1$zygosity)
df1$zygosity = factor(df1$zygosity, levels=c("WT","HET","HOM"))

pdf("iN_CADM2_vs_GUSB_zyg.pdf", height = 6, width = 3)
ggerrorplot(df1, x = "genotype", y = "fold_change_vs_wt", 
            desc_stat = "mean_sd",
            color = "zygosity",
            size = 1.4,
            legend="center",
            ylab="Fold Change vs WT",
            xlab="") +
  rotate_x_text()
dev.off()



