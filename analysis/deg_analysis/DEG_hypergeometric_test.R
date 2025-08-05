##Make Venn Diagrams of 

library(tidyverse)
library(ggplot2)

setwd("D:\\data\\talkowski\\Samples\\cadm2\\data\\output")
stats_files <- list()

sampleFiles <- grep("iN_DEG_stats",list.files(pattern = "\\.txt"),value=TRUE)
for (sample in sampleFiles){
  temp <- read.table(sample)
  #tissue_name <- sample#paste0(strsplit(sample,"_")[[1]][c(1,4,5,6,7,8)],"_direction")
  counter <- length(stats_files) + 1
  stats_files[[counter]] <- temp %>% 
    filter(pvalue < 0.05) %>%
    mutate(direction = ifelse(log2FoldChange >0,"UpRegulated", ifelse(log2FoldChange <0,"DownRegulated","NoChange"))) %>%
    rownames_to_column("Ensembl") %>%
    select(Ensembl,direction)
  colnames(stats_files[[counter]]) <- c('Ensembl',sample)
  name <- gsub(".txt","",sample)
  write.table(stats_files[[counter]]$Ensembl,file=paste0(name,"_pvalue_ensembl.txt"),sep="\t",quote=F,row.names=F, col.names = F)
}

hom_del <- stats_files[[4]]$Ensembl
het_del_100kb <- stats_files[[1]]$Ensembl

length(hom_del)
length(het_del_100kb)
length(intersect(hom_del, het_del_100kb))

overlapped_degs <-334
pws_degs<-1639
cadm2_degs <-2162
background_genes<- 16271
p <- phyper(overlapped_degs-1, pws_degs, background_genes-pws_degs, cadm2_degs, lower.tail=FALSE)

