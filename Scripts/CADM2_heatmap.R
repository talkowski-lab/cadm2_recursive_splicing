library(tidyverse)
library(ggplot2)
library(DESeq2)

setwd("D:\\data\\talkowski\\Samples\\cadm2\\data\\counts4DESeq\\RNASeq")
directory <- getwd()


sampleTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\sampleTableCount.txt", header=TRUE, sep = '\t')
libSizeTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\data\\Alignment_CADM2_RNAseq\\cadm2_rnaseq_library_size.txt", header=TRUE, sep = '\t')#For human elp1 injected mouse reference
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

sampleTable <- sampleTable %>% 
  filter(celltype=="iN")

sampleTable <- sampleTable[5:nrow(sampleTable),]

library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~1)


ddsHTSeq
source("D:\\data\\talkowski\\Samples\\cadm2\\scripts\\countsFilter.R")
map_directory <- "D:\\data\\talkowski\\Samples\\cadm2\\data\\Alignment_CADM2_RNAseq/"
pathList <- grep("iN",list.files(map_directory), value=TRUE)
pathList <- pathList[-grep(c(".lsf|r1B6"),pathList)]
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
inv <- grep("het|inv|hom",gc_names)
exp=apply(cpm,1,FUN=countsFilter,conditionList=list(wt,inv))
dds=ddsHTSeq[exp,]

#sizeFactors(ddsTranscript) <- sizeFactors(dds)
dds <- estimateSizeFactors(dds)

count_mat <- counts(dds)


setwd("C:\\Documents and Settings/ricar/Downloads/CADM2")
cell_markers_df <- read.csv(file="cellTypeMarkers.txt",head=T,sep="\t")
#synaptic_markers_df <- read.csv(file="SynpaticMarkers.txt",head=T,sep="\t")

cell_markers_df$gene <- paste(cell_markers_df$Symbol,cell_markers_df$EnsemblID,sep="|")

#exp.design <- metadata[,c(3:55,57,58,60:63)]
#rownames(exp.design) <- colnames(DESeqdata_beforecollapsing)

dds_0 <- DESeq(ddsHTSeq)

nst_0 <- normTransform(dds_0)

norm_counts_0 <- counts(dds_0,normalized=TRUE)

markergenes_norm_counts <- subset(norm_counts_0,(rownames(norm_counts_0) %in% cell_markers_df$EnsemblID))


markergenes_norm_counts_selected <- markergenes_norm_counts[,colnames(markergenes_norm_counts) %in% sampleTable$SampleNames]
markergenes_norm_counts_selected <- as.data.frame(markergenes_norm_counts_selected)
markergenes_norm_counts_selected$EnsemblID <- rownames(markergenes_norm_counts_selected)

markergenes_norm_count_selected_details <- inner_join(cell_markers_df,markergenes_norm_counts_selected,by="EnsemblID")


write.table(markergenes_norm_count_selected_details,file="MarkerGenesNormalizedCounts.txt",sep="\t",quote=F,row.names=F)

markergenes_nst <- subset(assay(nst_0),(rownames(assay(nst_0)) %in% cell_markers_df$EnsemblID))

print(markergenes_nst)

markergenes_nst <- markergenes_nst[match(cell_markers_df$EnsemblID,rownames(markergenes_nst)),]

print(markergenes_nst)

rownames(markergenes_nst) <- cell_markers_df$Symbol

print(markergenes_nst)

annotation_col <- data.frame(
  target="CADM2",
  #                  group=target_group,
  treatment=sampleTable$genotype,
  passage=sampleTable$CRISPR_round)



annotation_row=data.frame(row.names=cell_markers_df$Symbol,class=cell_markers_df$Category)

rownames(annotation_col) <- colnames(markergenes_nst)

annotation_col2 <- data.frame(
  target="CADM2")

rownames(annotation_col2) <- rownames(markergenes_nst)


pdf(file=paste0(prefix,".markergene.heatmap.cluster_col_annotation_row.logtpm.pdf"),height=14,width=8)
pheatmap(markergenes_nst,annotation_row=annotation_row,annotation_col=annotation_col,cluster_rows=FALSE,cluster_cols=TRUE,scale="none",fontsize=9,show_colnames=FALSE)
dev.off()

pdf(file=paste0(prefix,".markergene.heatmap_annotation_row.logtpm.pdf"),height=14,width=8)
pheatmap(markergenes_nst,annotation_row=annotation_row, annotation_col=annotation_col,cluster_rows=FALSE,cluster_cols=TRUE,scale="none",fontsize=9,show_colnames=TRUE,fontsize_col=8)
dev.off()

pdf(file=paste0(prefix,".markergene.heatmap.logtpm.customlegend.pdf"),height=14,width=8)
pheatmap(markergenes_nst,annotation_col=annotation_col2,cluster_rows=FALSE,cluster_cols=FALSE,scale="none",fontsize=9,show_colnames=FALSE,legend_breaks=c(0,3,6,9,12,15,max(markergenes_nst)),
         legend_labels = c("0", "3", "6", "9", "12" , "15", "log2(norm counts)\n"),legend=TRUE)
dev.off()

pdf(file=paste0(prefix,".markergene.heatmap.logtpm.pdf"),height=14,width=8)
pheatmap(markergenes_nst,annotation_row=annotation_row, annotation_col=annotation_col2,cluster_rows=FALSE,cluster_cols=FALSE,scale="none",fontsize=9,show_colnames=FALSE)
dev.off()

pdf(file=paste0(prefix,".markergene.heatmap.logtpm.pdf"),height=14,width=8)
pheatmap(markergenes_nst,annotation_row=annotation_row, annotation_col=annotation_col2,cluster_rows=FALSE,cluster_cols=TRUE,scale="none",fontsize=9,show_colnames=FALSE)
dev.off()