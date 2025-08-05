library(tidyverse)
setwd("D:\\data\\talkowski\\Samples\\cadm2\\data\\counts4DESeq\\RNASeq_EditedGTF")
directory <- getwd()

sampleTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\sampleTableCount.txt", header=TRUE, sep = '\t')
#libSizeTable <- read.table("/data/talkowski/Samples/fd_mouse_tissue/training_MED/qc/bamqc/libsz.txt", header=TRUE, sep = '\t')
libSizeTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\data\\Alignment_CADM2_RNAseq_EditedGTF\\libSize.txt", header=TRUE, sep = '\t')#For human elp1 injected mouse reference
conditionTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\conditions.txt", header=TRUE, sep = '\t')
colnames(libSizeTable) <- c("SampleNames","librarySize")
sampleTable$librarySize <- right_join(sampleTable, libSizeTable)$librarySize
sampleTable <- sampleTable[,c("SampleNames","FileNames","Samples","celltype","ID","replicate","zygosity","size","mutation_type","CRISPR_round","librarySize")]
sampleTable$genotype <- paste0(sampleTable$zygosity,"_",sampleTable$size,"_",sampleTable$mutation_type)
sampleTable$Samples <- paste0(sampleTable$celltype,"_",sampleTable$ID,"_",sampleTable$replicate) 
###Add metadata rom Sequencing Manifest
seqTable <- read.table("C:\\Documents and Settings\\ricar\\Downloads\\CADM2\\CADM2_sequencing _metadata.txt", header = T, sep='\t')
sampleTable<- right_join(sampleTable, seqTable)


##

comparison <- c("wt_wt_wt","het_500bp_inv")
sampleTable <- sampleTable %>% 
  filter(celltype=="iPSC") %>%
  filter(genotype %in% comparison) %>%
  filter(SampleNames!="iPSC_rB10_rep2_500bp_inv.count")

#sampleTable$celltype  <- factor(sampleTable$celltype, levels = c("iN","iPSC"))
sampleTable$zygosity <- factor(sampleTable$zygosity,levels=c("wt","het"))
sampleTable$mutation_type <- factor(sampleTable$mutation_type,levels=c("wt","inv"))
sampleTable$replicate <- factor(sampleTable$replicate,levels=c("rep1","rep2", "rep3", "rep4"))
sampleTable$genotype <- factor(sampleTable$genotype,levels=c("wt_wt_wt","het_500bp_inv"))


library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~1)


ddsHTSeq
source("D:\\data\\talkowski\\Samples\\cadm2\\scripts\\countsFilter.R")
map_directory <- "D:\\data\\talkowski\\Samples\\cadm2\\data\\Alignment_CADM2_RNAseq_EditedGTF/"
pathList <- grep("iPSC",list.files(map_directory), value=TRUE)
pathList <- pathList[-grep(c("hom|het_del|iPSC_rB10_rep2_500bp_inv"),pathList)]
#pathList <- pathList[-grep("iPSC_rB10_rep2_500bp_inv",pathList)]#sample had few lo unique reads
gc=counts(ddsHTSeq,normalized=F)
libsz=getTotalReads(pathList, map_directory)


libsz_names <- paste0(rownames(as.data.frame(libsz)), ".count")
gc_names <- colnames(gc)
for (i in 1:length(gc_names)){
  if (is.na(match(gc_names[i],libsz_names[i]))){
    cat(gc_names[i],",",libsz_names[i],"\n")
  }
}
cpm=t(1000000*t(gc)/libsz)
wt <- grep("wt",gc_names)
inv <- grep("inv",gc_names)
exp=apply(cpm,1,FUN=countsFilter,conditionList=list(wt,inv))
dds=ddsHTSeq[exp,]

dds <- estimateSizeFactors(dds)

source("C:\\Users\\ricar\\Downloads\\fd_mouse_tissue\\mouse\\Scripts\\rnaPCA.R")
cell_pca <- rnaPCA(log2(counts(dds, normalized=T) + 1),sampleTable = colData(dds), "genotype", PC1 = 1,PC2 = 2, textLabel = "Samples") 
pdf("..\\..\\output\\uncorrected_iPSC_WT_vs_HET_500bp_INV.pdf")
cell_pca
dev.off() 


####DEGs
library(sva)
dat <- counts(dds, normalized=TRUE)
mod  <- model.matrix(~ genotype, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0)

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]


design(ddssva) <- ~ genotype + SV1
#####Correct and generate PCA Plots
ddssva <- DESeq(ddssva, modelMatrixType="standard")
##########Creating corrected matrix
m=model.matrix(~genotype + SV1, data=colData(ddssva))
m=as.data.frame(m)
colnames(m)
m1=m[,3:ncol(m)]
X=as.matrix(m1)
beta=coef(ddssva)[,3:ncol(m)]
beta=as.matrix(beta)
y=log2(counts(ddssva,normalized=T))
cleany=y-t(as.matrix(X) %*% t(beta))
cleanY=2^cleany
write.table(cleany, "..\\..\\output\\iPSC_WT_vs_HET_500bp_INV_cleany.txt", row.names=TRUE, quote=FALSE, col.names = TRUE)

corrected_PCA_genotype <- rnaPCA(cleany,sampleTable = colData(dds), "genotype", PC1 = 1,PC2 = 2, textLabel = "Samples") 
pdf("..\\..\\output\\corrected_iPSC_WT_vs_HET_500bp_INV_labels.pdf")
corrected_PCA_genotype
dev.off() 


library(Hmisc)
library(corrplot)
source("C:\\Users\\ricar\\Downloads\\fd_mouse_tissue\\mouse\\Scripts\\rnaPCA.R")
test_pc <- rnaPCA(log2(counts(dds, normalized=T) + 1),sampleTable = colData(ddssva), "genotype", PC1 = 1,PC2 = 2, returnData = T)
test1 <- cbind(test_pc[,c("genotype","mutation_type","zygosity", "replicate","librarySize","RIN",
                          "Conc_to_Use",
                          "Input_amount",
                          "Average_Size_bp",
                          "Conc_ng_per_microlitre",
                          "Region_Molarity_nmol_per_l",
                          "qPCR_results_nmol_per_l",
                          "RR_of_Agilent",
                          "Ratio_of_Ratios_TS_vs_qPCR",
                          "Conc_To_use",
                          "sizeFactor")],test_pc[1:7],test_pc[32])


test1$replicate <- as.numeric(test1$replicate)
test1$zygosity <- as.numeric(test1$zygosity)
test1$mutation_type<-as.numeric(test1$mutation_type)
test1$genotype <- as.numeric(test1$genotype)
cormat<- rcorr(as.matrix(test1))
diag(cormat$P) <- 0

pdf("..\\..\\output\\corplot_uncorrected_cleaned_genotype_model_het_500bp_inv_vs_wt.pdf")
corrplot(cormat$r,type="upper", order = "original", p.mat = cormat$P, sig.level = 0.05, insig = "blank")
dev.off()


#####Correlation after correction
library(Hmisc)
library(corrplot)
test_pc <- rnaPCA(cleany,sampleTable = colData(ddssva), "genotype", PC1 = 1,PC2 = 2, returnData = T)
test1 <- cbind(test_pc[,c("genotype","mutation_type","zygosity", "replicate","librarySize","RIN",
                          "Conc_to_Use",
                          "Input_amount",
                          "Average_Size_bp",
                          "Conc_ng_per_microlitre",
                          "Region_Molarity_nmol_per_l",
                          "qPCR_results_nmol_per_l",
                          "RR_of_Agilent",
                          "Ratio_of_Ratios_TS_vs_qPCR",
                          "Conc_To_use",
                          "sizeFactor")],test_pc[1:7],test_pc[32])

test1$replicate <- as.numeric(test1$replicate)
test1$zygosity <- as.numeric(test1$zygosity)
test1$mutation_type<-as.numeric(test1$mutation_type)
test1$genotype <- as.numeric(test1$genotype)
cormat<- rcorr(as.matrix(test1))
diag(cormat$P) <- 0

pdf("..\\..\\output\\corplot_corrected_cleaned_genotype_model_het_500bp_inv_vs_wt.pdf")
corrplot(cormat$r,type="upper", order = "original", p.mat = cormat$P, sig.level = 0.05, insig = "blank")
dev.off()


####Print DEGs FOR hom_del
resultsNames(ddssva)
res <- results(ddssva,contrast=c(0,1,0))
res_genotype <- as.data.frame(res)
diff_genotype <-res_genotype[which(res_genotype$padj < 0.1 & abs(res_genotype$log2FoldChange)>=log2(1)),]
write.table(as.data.frame(rownames(diff_genotype)), "..\\..\\output\\iPSC_DEG_WT_vs_het_500bp_inv_q_0.1_ensembl_names_dropped_sample.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)                                                                                                                                
write.table(res_genotype, "..\\..\\output\\iPSC_DEG_stats_WT_vs_HET_500bp_inv_dropped_sample.txt", row.names=TRUE, quote=FALSE, col.names = TRUE)
write.table(as.data.frame(rownames(dds)),"..\\..\\output\\iPSC_background_ensembl_WT_vs_HET_500bbp_inv_names_dropped_sample.txt",row.names=FALSE, quote=FALSE, col.names = FALSE)

print("Genotype Comparison FDR < 0.1 and FC log2(1.5)")
table(res_genotype$padj < 0.1 & abs(res_genotype$log2FoldChange)>=log2(1.5))

print("Genotype Comparison FDR < 0.1")
table(res_genotype$padj < 0.1 & abs(res_genotype$log2FoldChange)>=log2(1))

print("Genotype Comparison pvalue < 0.05")
table(res_genotype$pvalue < 0.05 & abs(res_genotype$log2FoldChange)>=log2(1))

pdf("..\\..\\output\\iPSC_WT_vs_HET_500_bp_inv_single_comparison_qqplot_dropped_sample.pdf")
my.pvalue.sort = -log10(sort(res$pvalue))
exp.pvalue = runif(length(my.pvalue.sort),min=0,max=1)
exp.pvalue.sort = -log10(sort(exp.pvalue))
qplot(exp.pvalue.sort,my.pvalue.sort) + geom_abline(intercept = 0, slope = 1) + labs(x="-log10 Expected p-value (Uniform)", y="-log10 p-value")
dev.off()

####Get Gene Names
gene_Names<-read.table("C:\\Users\\ricar\\Downloads\\CADM2\\GRCh38.92_ensemblID_geneSymbol_pairs.txt", sep='\t')
ensembl <- as.data.frame(res)
ensembl$V1 <- rownames(ensembl)
gene_ids <-left_join(ensembl, gene_Names)$V2
library(EnhancedVolcano)
pdf("..\\..\\output\\iPSC_WT_vs_HET_500_bp_inv_single_comparison_volcano_dropped_sample.pdf")
EnhancedVolcano(res,
                lab = gene_ids,
                x = 'log2FoldChange',
                y = 'padj',
                #xlim = c(-2,2),
                #ylim=c(0,5),
                pCutoff = 0.1,
                FCcutoff = log2(1.5),
                pointSize = 3,
                labSize = 3)

dev.off()