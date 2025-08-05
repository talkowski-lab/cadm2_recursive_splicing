library(ggplot2)
library(tidyverse)
library(DESeq2)


setwd("D:\\data\\talkowski\\Samples\\cadm2\\data\\counts4DESeq\\RNASeq_EditedGTF")
directory <- getwd()


sampleTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\sampleTableCount.txt", header=TRUE, sep = '\t')
libSizeTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\data\\libSize.txt", header=TRUE, sep = '\t')#For human elp1 injected mouse reference
conditionTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\conditions.txt", header=TRUE, sep = '\t')
colnames(libSizeTable) <- c("SampleNames","librarySize")
sampleTable$librarySize <- right_join(sampleTable, libSizeTable)$librarySize
sampleTable <- sampleTable[,c("SampleNames","FileNames","Samples","celltype","ID","replicate","zygosity","size","mutation_type","CRISPR_round","librarySize")]
sampleTable$genotype <- paste0(sampleTable$zygosity,"_",sampleTable$size,"_",sampleTable$mutation_type)
sampleTable$Samples <- paste0(sampleTable$celltype,"_",sampleTable$ID,"_",sampleTable$replicate) 
###Add metadata rom Sequencing Manifest
#seqTable <- read.table("E:\\Documents and Settings\\ricar\\Downloads\\CADM2\\CADM2_sequencing _metadata.txt", header = T, sep='\t')
#sampleTable<- right_join(sampleTable, seqTable)


##

#comparison <- c("wt_wt_wt","hom_500bp_del")
sampleTable <- sampleTable %>% 
  filter(celltype=="iN") %>%
  filter(ID!="r1B6") %>%
  #filter(genotype %in% comparison)
  #sampleTable <- sampleTable[1:8,]
sampleTable <- sampleTable[1:20,]

sampleTable$zygosity <- factor(sampleTable$zygosity,levels=c("wt","het","hom"))
sampleTable$genotype <- factor(sampleTable$genotype,levels=c("wt_wt_wt","het_500bp_del", "het_500bp_inv","het_100kb_del","hom_500bp_del"))


library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~1)


ddsHTSeq
source("D:\\data\\talkowski\\Samples\\cadm2\\scripts\\countsFilter.R")
map_directory <- "D:\\data\\talkowski\\Samples\\cadm2\\data\\Alignment_CADM2_RNAseq_EditedGTF/"
pathList <- grep("iN",list.files(map_directory), value=TRUE)
pathList <- pathList[-grep(c("r1B6|lsf"),pathList)]
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
wt <- grep("wt|total",gc_names)
het <- grep("het",gc_names)
hom <- grep("hom",gc_names)
exp=apply(cpm,1,FUN=countsFilter,conditionList=list(wt,het, hom))
dds=ddsHTSeq[exp,]

#sizeFactors(ddsTranscript) <- sizeFactors(dds)
dds <- estimateSizeFactors(dds)


source("C:\\Users\\ricar\\Downloads\\fd_mouse_tissue\\mouse\\Scripts\\rnaPCA.R")
cell_pca <- rnaPCA(log2(counts(dds, normalized=T) + 1),sampleTable = colData(dds), "zygosity", PC1 = 1,PC2 = 2, textLabel = "Samples") 
pdf("..\\..\\output\\uncorrected_iN_WT_vs_het_hom_no_sva.pdf")
cell_pca
dev.off() 


####DEGs
library(sva)
dat <- counts(dds, normalized=TRUE)
mod  <- model.matrix(~ zygosity, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0)

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
ddssva$SV4 <- svseq$sv[,4]
ddssva$SV5 <- svseq$sv[,5]

design(ddssva) <- ~ zygosity  + SV1 + SV2+ SV3 + SV4 + SV5
#####Correct and generate PCA Plots
ddssva <- DESeq(ddssva, modelMatrixType="standard")
##########Creating corrected matrix
m=model.matrix(~zygosity + SV1 + SV2+ SV3 + SV4 + SV5, data=colData(ddssva))
m=as.data.frame(m)
colnames(m)
m1=m[,4:ncol(m)]
X=as.matrix(m1)
beta=coef(ddssva)[,4:ncol(m)]
beta=as.matrix(beta)
y=log2(counts(ddssva,normalized=T))
cleany=y-t(as.matrix(X) %*% t(beta))
cleanY=2^cleany
write.table(y, "..\\..\\output\\iN_WT_vs_het_hom_pws_no_sva_cleany.txt", row.names=TRUE, quote=FALSE, col.names = TRUE)

corrected_PCA_genotype <- rnaPCA(counts(ddssva),sampleTable = colData(dds), "genotype", PC1 = 1,PC2 = 2, textLabel = "Samples") 
pdf("..\\..\\output\\corrected_iN_WT_vs_het_hom_sva_labels_genotype.pdf")
corrected_PCA_genotype
dev.off() 


library(Hmisc)
library(corrplot)
source("C:\\Users\\ricar\\Downloads\\fd_mouse_tissue\\mouse\\Scripts\\rnaPCA.R")
test_pc <- rnaPCA(log2(counts(dds, normalized=T) + 1),sampleTable = colData(ddssva), "zygosity", PC1 = 1,PC2 = 2, returnData = T)
test1 <- cbind(test_pc[,c("genotype","mutation_type","zygosity", "replicate","librarySize",
                          "sizeFactor")],test_pc[1:10])


test1$replicate <- as.numeric(test1$replicate)
test1$zygosity <- as.numeric(test1$zygosity)
test1$mutation_type<-as.numeric(test1$mutation_type)
test1$genotype <- as.numeric(test1$genotype)
cormat<- rcorr(as.matrix(test1))
diag(cormat$P) <- 0

pdf("..\\..\\output\\corplot_uncorrected_cleaned_genotype_model_iN_het_hom_vs_wt_no_sva.pdf")
corrplot(cormat$r,type="upper", order = "original", p.mat = cormat$P, sig.level = 0.05, insig = "blank")
dev.off()


#####Correlation after correction
library(Hmisc)
library(corrplot)
test_pc <- rnaPCA(counts(ddssva),sampleTable = colData(ddssva), "zygosity", PC1 = 1,PC2 = 2, returnData = T)
test1 <- cbind(test_pc[,c("genotype","mutation_type","zygosity", "replicate","librarySize",
                          "sizeFactor")],test_pc[1:10])

test1$replicate <- as.numeric(test1$replicate)
test1$zygosity <- as.numeric(test1$zygosity)
test1$mutation_type<-as.numeric(test1$mutation_type)
test1$genotype <- as.numeric(test1$genotype)
cormat<- rcorr(as.matrix(test1))
diag(cormat$P) <- 0

pdf("..\\..\\output\\corplot_corrected_cleaned_genotype_model_iN_het_hom_vs_wt_no_sva.pdf")
corrplot(cormat$r,type="upper", order = "original", p.mat = cormat$P, sig.level = 0.05, insig = "blank")
dev.off()


####Print DEGs FOR hom_del
resultsNames(ddssva)
res <- results(ddssva,contrast=c(0,0,1))
res_genotype <- as.data.frame(res)
diff_genotype <-res_genotype[which(res_genotype$padj < 0.1 & abs(res_genotype$log2FoldChange)>=log2(1)),]
write.table(as.data.frame(rownames(diff_genotype)), "..\\..\\output\\iN_DEG_WT_vs_hom_no_sva_q_0.1_ensembl.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)                                                                                                                                
write.table(res_genotype, "..\\..\\output\\iN_DEG_stats_WT_vs_hom_no_sva.txt", row.names=TRUE, quote=FALSE, col.names = TRUE)
write.table(as.data.frame(rownames(dds)),"..\\..\\output\\iN_background_ensembl_WT_vs_hom_no_sva.txt",row.names=FALSE, quote=FALSE, col.names = FALSE)

print("Genotype Comparison FDR < 0.1 and FC log2(1.5)")
table(res_genotype$padj < 0.1 & abs(res_genotype$log2FoldChange)>=log2(1.5))

print("Genotype Comparison FDR < 0.1")
table(res_genotype$padj < 0.1 & abs(res_genotype$log2FoldChange)>=log2(1))

print("Genotype Comparison pvalue < 0.05")
table(res_genotype$pvalue < 0.05 & abs(res_genotype$log2FoldChange)>=log2(1))

pdf("..\\..\\output\\iN_WT_vs_hom_no_sva_qqplot.pdf")
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
pdf("..\\..\\output\\iN_WT_vs_hom_no_sva_volcano.pdf")
EnhancedVolcano(res,
                lab = gene_ids,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-2.5,2.5),
                #ylim=c(0,5),
                pCutoff = 0.1,
                FCcutoff = log2(1.5),
                pointSize = 3,
                labSize = 3)

dev.off()

