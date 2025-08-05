  cell_markers_df <- read.csv(file="C:\\Users\\ricar\\Downloads\\cellTypeMarkers.txt",head=T,sep="\t")
  
  behaviour_markers <- read.csv(file="C:\\Users\\ricar\\Downloads\\CADM2\\gobp_behaviour_pairs.txt",head=F,sep="\t")
  colnames(behaviour_markers) <-c("EnsemblID","Symbol")
  behaviour_markers$Category <- "GOBP Behaviour"
  
  reg_markers <- read.csv(file="C:\\Users\\ricar\\Downloads\\CADM2\\REG_ACCTTT_UNKNOWN_pairs.txt",head=F,sep="\t")
  colnames(reg_markers) <-c("EnsemblID","Symbol")
  reg_markers$Category <- "REG"
  
  cell_markers_df$Category <- gsub("GABA","iGABA",cell_markers_df$Category)
  cell_markers_df$Category <- gsub("GLUT","iGLUT",cell_markers_df$Category) 
  
  behaviour_markers$gene <- paste(behaviour_markers$Symbol,behaviour_markers$EnsemblID,sep="|")
  
  cell_markers_df <- reg_markers
  #exp.design <- metadata[,c(3:55,57,58,60:63)]
  #rownames(exp.design) <- colnames(DESeqdata_beforecollapsing)
  ########
  
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
  
  
  sampleTable <- sampleTable %>% 
    #filter(celltype=="iN") %>%
    #filter(genotype %in% comparison) %>%
    filter(ID!="r1B6")
  
  sampleTable$genotype <- gsub("wt_wt_wt","wt",sampleTable$genotype)
  sampleTable <- sampleTable[1:40,]
  
  
  library("DESeq2")
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                         directory = directory,
                                         design= ~1)
  
  
  ddsHTSeq
  source("D:\\data\\talkowski\\Samples\\cadm2\\scripts\\countsFilter.R")
  map_directory <- "D:\\data\\talkowski\\Samples\\cadm2\\data\\Alignment_CADM2_RNAseq_EditedGTF/"
  pathList <- grep("iN|iPSC",list.files(map_directory), value=TRUE)
  pathList <- pathList[-grep(c(".lsf|r1B6|total"),pathList)]
  #pathList <- pathList[-grep(c("inv|het_del|.lsf|r1B6"),pathList)]
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
  exp=apply(cpm,1,FUN=countsFilter,conditionList=list(wt))
  dds=ddsHTSeq[exp,]
  
  
  dds_0 <- DESeq(ddsHTSeq)
  
  nst_0 <- normTransform(dds_0)
  
  norm_counts_0 <- counts(dds_0,normalized=TRUE)
  
  
  #############
  
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
    Genotype=sampleTable$genotype,
    Celltype=sampleTable$celltype)
  
  #annotation_col <- data.frame(
  #  target=sampleTable$mutation_type,
  #  #                  group=target_group,
  #  activity=sampleTable$CRISPR_round,
  #  Genotype=sampleTable$genotype,
  #  passage=sampleTable$replicate)
  
  
  
  comparison <- c("Housekeeping gene", "iPSC", "iGABA", "iGLUT", "Neuron","NPC")
  remove_genes <- c("EOMES", "GRIN1", "GRIN2B")
  
  tmp_df <- cell_markers_df %>% filter(Category %in% comparison) %>% filter(Symbol!="EOMES") %>% filter(Symbol!="GRIN1") %>% filter(Symbol!="GRIN2B")
  
  annotation_row=data.frame(row.names=cell_markers_df$Symbol,cell_markers_df$Category)
  
  rownames(annotation_col) <- colnames(markergenes_nst)
  
  annotation_col2 <- data.frame(
    Genotype=sampleTable$genotype,
    Celltype=sampleTable$celltype)
  
  rownames(annotation_col2) <- colnames(markergenes_nst)
  
  ting <- markergenes_nst[rownames(markergenes_nst) %in% cell_markers_df$Symbol,]
  setwd("C:\\Documents and Settings/ricar/Downloads/CADM2/")
  
  newnames <- lapply(
    cell_markers_df$Symbol,
    function(x) bquote(italic(.(x))))
  
  #annot_colors=list(class=c(iGABA="#440154FF",iGLUT="#414487FF",`Housekeeping gene`="#2A788EFF", iPSC="#22A884FF", Neuron="#7AD151FF", NPC="#FDE725FF"))
  annot_colors=list(class=c(Behaviour="#440154FF"))
  
library(pheatmap)
prefix="Behaviour"
pdf(file=paste0(prefix,".markergene.heatmap.cluster_col_annotation_row.logtpm.pdf"),height=11,width=8)
pheatmap(ting,annotation_row=annotation_row,annotation_col=annotation_col,annotation_colors=annot_colors,cluster_rows=FALSE,cluster_cols=TRUE,scale="none",fontsize=9,show_colnames=FALSE,labels_row = as.expression(newnames))
dev.off()

pdf(file=paste0(prefix,".markergene.heatmap_annotation_row.logtpm.pdf"),height=11,width=8)
pheatmap(ting,annotation_row=annotation_row, annotation_col=annotation_col,annotation_colors=annot_colors, cluster_rows=FALSE,cluster_cols=TRUE,scale="none",fontsize=9,show_colnames=TRUE,fontsize_col=8,labels_row = as.expression(newnames))
dev.off()

pdf(file=paste0(prefix,".markergene.heatmap.logtpm.customlegend.pdf"),height=11,width=8)
pheatmap(ting,annotation_col=annotation_col2,annotation_colors=annot_colors, cluster_rows=FALSE,cluster_cols=TRUE,scale="none",fontsize=9,show_colnames=FALSE,legend_breaks=c(0,3,6,9,12,15,max(markergenes_nst)),
         legend_labels = c("0", "3", "6", "9", "12" , "15", "log2(norm counts)\n"),legend=TRUE,   labels_row = as.expression(newnames))
dev.off()

pdf(file=paste0(prefix,".markergene.heatmap.logtpm.pdf"),height=11,width=8)
pheatmap(ting,annotation_row=annotation_row, annotation_col=annotation_col2,annotation_colors=annot_colors,cluster_rows=TRUE,cluster_cols=FALSE,scale="none",fontsize=9,show_colnames=FALSE,labels_row = as.expression(newnames))
dev.off()

pdf(file=paste0(prefix,".markergene.heatmap.logtpm.pdf"),height=11,width=8)
pheatmap(ting,annotation_row=annotation_row, annotation_col=annotation_col2,annotation_colors=annot_colors,cluster_rows=TRUE,cluster_cols=TRUE,scale="none",fontsize=9,show_colnames=FALSE,  labels_row = as.expression(newnames))
dev.off()

