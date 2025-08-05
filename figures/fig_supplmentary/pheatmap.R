library(pheatmap)


setwd("D:\\data\\talkowski\\Samples\\fd_mouse_tissue\\training_MED\\output\\pdf\\no_sex")

counts <- read.csv("DEG_counts_2.csv", header = T)
pvalues <- read.csv("DEG_pvalue_2.csv",header = T)

colnames(counts) <- c("","DRG Brown (778)","DRG Black (596)","DRG Salmon (418)","TG Cyan (370)","TG Red (875)")
colnames(pvalues) <- c("","DRG Brown (778)","DRG Black (596)","DRG Salmon (418)","TG Cyan (370)","TG Red (875)")
rownames(counts) <- c("DEGS in DRG (148)","DEGS in TG (194)","DEGs in PNS Convergence (44)","DRG Brown (778)","DRG Black (596)","DRG Salmon (418)","TG Cyan (370)", "TG Red (875)")
rownames(pvalues) <- c("DEGS in DRG (148)","DEGS in TG (194)","DEGs in PNS Convergence (44)","DRG Brown (778)","DRG Black (596)","DRG Salmon (418)","TG Cyan (370)", "TG Red (875)")
pheatmap(pvalues, display_numbers = counts, color = inferno())

         