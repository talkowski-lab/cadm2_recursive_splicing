#!/apps/lib-osver/R/3.4.0/bin/Rscript

args<-commandArgs(TRUE)
tissue <- args[1]
c <- args[2]
nameP<-args[3]


.libPaths(c("/PHShome/ry077/R/x86_64-pc-linux-gnu-library/3.4","/apps/lib-osver/R/3.4.0/lib64/R/library")) 
library("scales")
library("DEXSeq")
library("WriteXLS")
library("BiocParallel")
#load("data/RData/annot.RData")
load("D:\\data\\talkowski\\Samples\\mef2c\\data\\RData\\annot.RData")
BPPARAM = MulticoreParam(workers=4)

inDir="data/Counts_1"
flattenedFile = list.files("data/ref", pattern="_Mef2c.gff$", full.names=TRUE)
if(tissue=="NSC"){
    countFiles = list.files(inDir, pattern="^NSC", full.names=TRUE)
    sampleTable=read.table("docs/metaData_NSC.tab",header=T,sep="\t")
    rownames(sampleTable)<-sampleTable$Sample
} else {
    tissue="iN"
    countFiles = list.files(inDir, pattern="^iN", full.names=TRUE)
    sampleTable=read.table("docs/metaData_iNs.tab",header=T,sep="\t")
    rownames(sampleTable)<-sampleTable$Sample
}

eval(parse(text=paste0("sampleTable_0=sampleTable[sampleTable$Geno_DelType == \"WT\" | sampleTable$Geno_DelType == \"",c,"\",]")))
patterns=sampleTable_0$Sample
countFiles_1=countFiles[grepl(paste(patterns, collapse="|"),countFiles)]
sampleTable_1<-sampleTable_0[order(rownames(sampleTable_0)),]
sampleTable_1$Geno_DelType<-droplevels(sampleTable_1$Geno_DelType)
sampleTable_1$Geno_DelType<-relevel(sampleTable_1$Geno_DelType,"WT")
dxd = DEXSeqDataSetFromHTSeq(
countFiles_1,
sampleData=sampleTable_1,
design= ~ exon + Geno_DelType:exon,
flattenedfile=flattenedFile )
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd,BPPARAM=BPPARAM )
dxd = testForDEU( dxd ,fullModel = design(dxd),reducedModel = ~exon,BPPARAM=BPPARAM)
dxd = estimateExonFoldChanges( dxd, fitExpToVar="Geno_DelType")
dxr1 = DEXSeqResults( dxd )
dxr1<-as.data.frame(dxr1)[order(dxr1$padj),]
WriteXLS(dxr1,paste0("tables/GenesIn",c,"_",tissue,"_",nameP,".xlsx"),row.names = FALSE)
gene="ENSG00000175161"
name=annot[gene,7]
pdf(paste0("Figures/",c,"/",name,"_",tissue,"_",nameP,".pdf"))
plotDEXSeq( dxr1,  gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ,fitExpToVar="Geno_DelType")
dev.off()

save(counts(dxd),dxr1,file=paste0("data/RData/DEXSeq_",tissue,"_",c,"_",nameP,".RData"))


##############
library(ggplot2)
library(tidyverse)
library(DESeq2)
library(DEXSeq)
library("scales")
library("DEXSeq")
library("WriteXLS")
library("BiocParallel")
#load("data/RData/annot.RData")
load("D:\\data\\talkowski\\Samples\\mef2c\\data\\RData\\annot.RData")
BPPARAM = MulticoreParam(workers=8)


#setwd("D:\\data\\talkowski\\Samples\\cadm2\\data\\counts4DEXSeq_Exon_Def")
setwd("D:\\data\\talkowski\\Samples\\cadm2\\data\\counts4DEXSeq_Exon_Def\\test2\\")
directory <- getwd()


sampleTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\sampleTableCount.txt", header=TRUE, sep = '\t')
libSizeTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\data\\libSize.txt", header=TRUE, sep = '\t')#For human elp1 injected mouse reference
conditionTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\conditions.txt", header=TRUE, sep = '\t')
colnames(libSizeTable) <- c("SampleNames","librarySize")
sampleTable$librarySize <- right_join(sampleTable, libSizeTable)$librarySize
sampleTable <- sampleTable[,c("SampleNames","FileNames","Samples","celltype","ID","replicate","zygosity","size","mutation_type","CRISPR_round","librarySize")]
sampleTable$genotype <- paste0(sampleTable$zygosity,"_",sampleTable$size,"_",sampleTable$mutation_type)
sampleTable$Samples <- paste0(sampleTable$celltype,"_",sampleTable$ID,"_",sampleTable$replicate) 

sampleTable$Genotype[grepl("hom_del",sampleTable$SampleNames)] <- "HOM 500bp DEL"
sampleTable$Genotype[grepl("500bp_het_del",sampleTable$SampleNames)] <- "HET 500bp DEL"
sampleTable$Genotype[grepl("_100kb_het_del",sampleTable$SampleNames)] <- "HET 100kb DEL"
sampleTable$Genotype[grepl("_500bp_inv",sampleTable$SampleNames)] <- "HET 500bp INV"
sampleTable$Genotype[grepl("_wt",sampleTable$SampleNames)] <- "WT"


sampleTable$Genotype <- factor(sampleTable$Genotype,levels=c("WT","HET 100kb DEL", "HET 500bp DEL", "HET 500bp INV", "HOM 500bp DEL"))
#sampleTable$Genotype <- factor(sampleTable$Genotype,levels=c("HET_DEL_100kb", "HET_DEL_500bp", "HET_INV", "HOM_DEL", "WT"))


sampleTable <- sampleTable %>% 
  filter(celltype=="iN") %>%
  filter(ID!="r1B6")
sampleTable <- sampleTable[1:20,]

#flattenedFile = "D:\\data\\talkowski\\Samples\\cadm2\\data\\ref\\ting.gff"
#system("cat D:\\data\\talkowski\\Samples\\cadm2\\data\\ref\\ting.gff")
#flattenedFile = "C:\\Users\\ricar\\Downloads\\CADM2\\cadm2_R1_R2.gff"
#system("cat C:\\Users\\ricar\\Downloads\\CADM2\\cadm2_R1_R2.gff")
flattenedFile = "C:\\Users\\ricar\\Downloads\\CADM2\\test.gff"
system("cat C:\\Users\\ricar\\Downloads\\CADM2\\test.gff")


dxd = DEXSeqDataSetFromHTSeq(
  sampleTable$FileNames,
  sampleData=sampleTable,
  design= ~ exon + genotype:exon,
  flattenedfile=flattenedFile)


dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions(dxd,BPPARAM=BPPARAM)
dxd = testForDEU( dxd ,fullModel = design(dxd),reducedModel = ~exon)
dxd = estimateExonFoldChanges( dxd, fitExpToVar="Genotype")
dxr1 = DEXSeqResults( dxd )
dxr1=dxr1[order(dxr1$featureID),]
gene="ENSG00000175161"
name=annot[gene,7]
plotDEXSeq( dxr1,  gene, legend=TRUE, cex.axis=1, cex=1, lwd=2 ,fitExpToVar="Genotype")

setwd("C:\\Documents and Settings/ricar/Downloads/CADM2")
pdf("cadm2_dexplot.pdf", width = 20)
plotDEXSeq( dxr1,  gene, legend=TRUE,cex.axis=2, cex=2, lwd=2 ,fitExpToVar="Genotype")
dev.off()

dxr2<-as.data.frame(dxr1)[order(dxr1$padj),]
WriteXLS(dxr2,paste0("dxr1.xlsx"),row.names = FALSE)

colnames(dxr2) <- c( "groupID" ,"featureID","exonBaseMean","dispersion",               
                     "stat","pvalue","padj","WT","HET_DEL_100kb","HET_DEL_500bp","HET_INV",  "HOM_DEL", 
                     "log2fold_HET_DEL_100kb_WT","log2fold_HET_DEL_500bp_WT","log2fold_HET_INV_WT","log2fold_HOM_DEL_WT",      
                     "genomicData.seqnames","genomicData.start","genomicData.end","genomicData.width",        
                     "genomicData.strand","HOM_DEL_count","HOM_DEL_count","HOM_DEL_count",              
                     "HOM_DEL_count","HET_DEL_500bp_count","HET_INV_count","HET_INV_count",              
                     "HET_DEL_500bp_count","HET_DEL_500bp_count","HET_INV_count","HET_INV_count",             
                     "HET_DEL_500bp_count","WT_count","HET_DEL_100kb_count","HET_DEL_100kb_count",             
                     "WT_count","WT_count","HET_DEL_100kb_count","WT_count",           
                     "HET_DEL_100kb_count","transcripts")

ting <- dxr2[22:41]
wt <- ting[c(13,16,17,19)]
hom<- ting[1:4]
het_100 <-ting[c(14,15,18,20)]
het_500 <-ting[c(5,8,9,12)]
inv <-ting[c(6,7,10,11)]

for (exon in 1:nrow(ting)){
  
  #wt_median <- median(as.numeric(wt[exon,]))
  #inv_fc<-mean(as.numeric(inv[exon,])/wt_median)
  #het_100_fc<- mean(as.numeric(het_100[exon,])/wt_median)
  #hom_fc <- mean(as.numeric(hom[exon,])/wt_median)
  #het_500_fc <- mean(as.numeric(het_500[exon,])/wt_median)
  
  #print(paste("FC:", inv_fc, het_500_fc, het_100_fc, hom_fc))
  cat(paste(rownames(ting)[exon],t.test(wt[exon,],inv[exon,])$p.value,
              t.test(wt[exon,],het_500[exon,])$p.value,
              t.test(wt[exon,],het_100[exon,])$p.value,
              t.test(wt[exon,],hom[exon,])$p.value))
}


tmp <- new_table %>% filter(transcript_id==transcript)
tmp_wt <- tmp %>% filter(genotype=="WT")
tmp$FC<- tmp$absolute_counts/median(tmp_wt$absolute_counts)
tmp_wt <- tmp %>% filter(genotype=="WT")
het_inv <- tmp %>% filter(genotype=="HET_INV")
het_500 <- tmp %>% filter(genotype=="HET_DEL_500bp")
het_100 <-tmp %>% filter(genotype=="HET_DEL_100kb")
hom_del <-tmp %>% filter(genotype=="HOM_DEL")
