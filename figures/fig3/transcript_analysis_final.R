library(tidyverse)
library(ggplot2)
dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/raw_output.tsv", sep='\t', header=T)
colnames(dat) <- c("sample","transcript_id","TPM")
dat$sample <- gsub(" ", "", dat$sample)
dat$transcript_id <- gsub(" ", "", dat$transcript_id)


dat$transcript_id <- gsub("::3:85505218-86067982|::3:84958964-86067982|::3:84958964-86067982|::3:84958958-86067982|::3:85726515-86067982|::3:85802173-85848021|::3:85961626-86065626|::3:85541799-85543286", "", dat$transcript_id)
dat$celltype[grepl("iN_",dat$sample)] <- "iN"
dat$celltype[grepl("iPSC_",dat$sample)] <- "iPSC"

dat$genotype[grepl("hom_del",dat$sample)] <- "HOM_DEL"
dat$genotype[grepl("500bp_het_del",dat$sample)] <- "HET_DEL_500bp"
dat$genotype[grepl("_100kb_het_del",dat$sample)] <- "HET_DEL_100kb"
dat$genotype[grepl("_500bp_inv",dat$sample)] <- "HET_INV"
dat$genotype[grepl("_wt",dat$sample)] <- "WT"
dat$genotype[grepl("GM_Control",dat$sample)] <- "WT"

dat$genotype <- factor(dat$genotype, levels = c("WT","HET_DEL_500bp","HET_DEL_100kb","HET_INV","HOM_DEL"))
dat$transcript_id <- factor(dat$transcript_id, levels=c("CDS_Long+RS1+RS2",
                                                        "CDS_Long+RS2",
                                                        "CDS_Short",
                                                        "AltExon1_RS2_CDS_Short",
                                                        "Novel_1",
                                                        "Novel_2",
                                                        "Novel_3"))


iN_TPM <- dat %>%
  filter(celltype=="iN") %>%
  filter(!grepl("r1B6",sample)) %>%
  filter(transcript_id != "Novel_1", transcript_id!= "Novel_2") %>%
  filter(!grepl("iN_total_GM_Control_1_",sample))

#cleany_dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/cleany_cadm2_no_pws.txt", sep='\t', header = T)
cleany_dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/cadm2_raw_counts.txt", sep='\t', header = T)
cleany_dat$sample <- gsub(".count","",cleany_dat$sample)
df_counts <- cleany_dat %>% group_by(sample) %>% summarise(mean = cleany %>% mean())


new_table <-left_join(iN_TPM,cleany_dat)
#new_table$absolute_counts <- 2^new_table$cleany*((new_table$TPM)/1e6)
new_table$absolute_counts <- new_table$cleany*((new_table$TPM)/1e6)

libSizeTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\data\\libSize.txt", header=TRUE, sep = '\t')#For human elp1 injected mouse reference
colnames(libSizeTable) <- c("sample","librarySize")
libSizeTable$sample <- gsub(".count", "", libSizeTable$sample)
new_table <- left_join(new_table, libSizeTable)
new_table$normalized_counts <- new_table$absolute_counts/new_table$librarySize

wt <- new_table %>% filter(genotype=="WT")
new_table$FC <- new_table$absolute_counts/median(wt$absolute_counts)
hom_del <- new_table %>% filter(genotype =="HOM_DEL")
inv_del <- new_table %>% filter(genotype =="HET_INV")
het_del_100kb <- new_table %>% filter(genotype =="HET_DEL_100kb")
het_del_500bp <- new_table %>% filter(genotype =="HET_DEL_500bp")


df_transcript <-new_table %>% group_by(sample,genotype)  %>% summarise(Normalized_Counts=sum(normalized_counts))


pdf("combined_transcript_genotype_normalized_expression_raw_uncalibrated.pdf")
ggplot(df_transcript, aes(x=genotype, y=Normalized_Counts, fill=genotype)) +
  #(paste0("iN_",transcript)) +
  geom_boxplot(width=0.7, position = 'dodge') +
  #ylim(0,3) +
  xlab("Counts") +
  #geom_hline(yintercept=1) +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #geom_jitter(position = position_jitter(seed = 1)) +
  #geom_text(position = position_jitter(seed = 1), size=3.5) +
  theme(text = element_text(size = 19),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

for (transcript in unique(new_table$transcript_id)){
  tmp <- new_table %>% filter(transcript_id==transcript)
  tmp_wt <- tmp %>% filter(genotype=="WT")
  tmp$FC<- tmp$absolute_counts/median(tmp_wt$absolute_counts)
  tmp_wt <- tmp %>% filter(genotype=="WT")
  het_inv <- tmp %>% filter(genotype=="HET_INV")
  het_500 <- tmp %>% filter(genotype=="HET_DEL_500bp")
  het_100 <-tmp %>% filter(genotype=="HET_DEL_100kb")
  hom_del <-tmp %>% filter(genotype=="HOM_DEL")
  
  tmp_plot <- ggplot(tmp, aes(x=genotype, y=FC, fill=genotype, label=sample)) +
    ggtitle(paste0("iN_",transcript)) +
    geom_boxplot(width=0.7, position = 'dodge') +
    #ylim(0,3) +
    xlab("Counts") +
    geom_hline(yintercept=1) +
    theme_classic() +
    geom_jitter(position=position_jitterdodge(0.5)) +
    #geom_jitter(position = position_jitter(seed = 1)) +
    #geom_text(position = position_jitter(seed = 1), size=3.5) +
    theme(text = element_text(size = 19),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(tmp_plot,file=paste0("iN_no_pws_raw_uncalibrated_",transcript,".pdf"), width=4, height = 4)
  
  #print(paste(transcript, mean(tmp_wt$FC),median(tmp_wt$FC)))
  #print(paste(transcript,t.test(tmp_wt$FC,het_inv$FC)$p.value,
   #               t.test(tmp_wt$FC,het_500$FC)$p.value,
  #               t.test(tmp_wt$FC,het_100$FC)$p.value,
  #                t.test(tmp_wt$FC,hom_del$FC)$p.value))
  #print(paste(transcript,mean(tmp_wt$FC),mean(het_inv$FC),mean(het_500$FC),mean(het_100$FC),mean(hom_del$FC)))
  print(paste(transcript,mean(tmp_wt$absolute_counts),mean(het_inv$absolute_counts),mean(het_500$absolute_counts),mean(het_100$absolute_counts),mean(hom_del$absolute_counts)))
  print(paste(transcript,mean(tmp_wt$TPM),mean(het_inv$TPM),mean(het_500$TPM),mean(het_100$TPM),mean(hom_del$TPM)))
  
  #print(paste(transcript,t.test(tmp_wt$absolute_counts,het_inv$absolute_counts)$p.value,
  #           t.test(tmp_wt$absolute_counts,het_500$absolute_counts)$p.value,
  #            t.test(tmp_wt$absolute_counts,het_100$absolute_counts)$p.value,
  #           t.test(tmp_wt$absolute_counts,hom_del$absolute_counts)$p.value))
  print("-----------------")
}
write.table(new_table, "iN_transcript_stats_no_pws_transcripts_raw_counts.txt", quote = F,row.names = F, col.names = T)



for (transcript in unique(new_table$transcript_id)){
  tmp <- new_table %>% filter(transcript_id==transcript)
  tmp_wt <- tmp %>% filter(genotype=="WT")
  tmp$FC<- tmp$absolute_counts/median(tmp_wt$absolute_counts)
  tmp_wt <- tmp %>% filter(genotype=="WT")
  het_inv <- tmp %>% filter(genotype=="HET_INV")
  het_500 <- tmp %>% filter(genotype=="HET_DEL_500bp")
  het_100 <-tmp %>% filter(genotype=="HET_DEL_100kb")
  hom_del <-tmp %>% filter(genotype=="HOM_DEL")
  
  tmp_plot <- ggplot(tmp, aes(x=genotype, y=FC, fill=genotype, label=sample)) +
    ggtitle(paste0("iN_",transcript)) +
    geom_boxplot(width=0.7, position = 'dodge') +
    #ylim(0,3) +
    xlab("FC") +
    geom_hline(yintercept=1) +
    theme_classic() +
    geom_jitter(position=position_jitterdodge(0.5)) +
    #geom_jitter(position = position_jitter(seed = 1)) +
    #geom_text(position = position_jitter(seed = 1), size=3.5) +
    theme(text = element_text(size = 19),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(tmp_plot,file=paste0("iN_no_pws_raw_uncalibrated_FC_",transcript,".pdf"), width=6, height = 4)
  
  
  df_transcript <- new_table %>% filter(transcript_id==transcript) %>% group_by(genotype) %>% summarise(mean(normalized_counts))
  tmp_plot_1 <-ggplot(df_transcript_RS1_RS2, aes(x=genotype, y=`mean(normalized_counts)`, fill=genotype)) +
    #(paste0("iN_",transcript)) +
    geom_boxplot(width=0.7, position = 'dodge') +
    #ylim(0,3) +
    xlab("") +
    #geom_hline(yintercept=1) +
    theme_classic() +
    geom_jitter(position=position_jitterdodge(0.5)) +
    #geom_jitter(position = position_jitter(seed = 1)) +
    #geom_text(position = position_jitter(seed = 1), size=3.5) +
    theme(text = element_text(size = 19),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(tmp_plot_1,file=paste0("iN_no_pws_raw_uncalibrated_normalized_counts_combined_",transcript,".pdf"), width=6, height = 4)
  #print(paste(transcript, mean(tmp_wt$FC),median(tmp_wt$FC)))
  #print(paste(transcript,t.test(tmp_wt$FC,het_inv$FC)$p.value,
  #                t.test(tmp_wt$FC,het_500$FC)$p.value,
  #               t.test(tmp_wt$FC,het_100$FC)$p.value,
  #                t.test(tmp_wt$FC,hom_del$FC)$p.value))
  #print(paste(transcript,mean(tmp_wt$FC),mean(het_inv$FC),mean(het_500$FC),mean(het_100$FC),mean(hom_del$FC)))
  #print(paste(transcript,mean(tmp_wt$absolute_counts),mean(het_inv$absolute_counts),mean(het_500$absolute_counts),mean(het_100$absolute_counts),mean(hom_del$absolute_counts)))
  #print(paste(transcript,mean(tmp_wt$TPM),mean(het_inv$TPM),mean(het_500$TPM),mean(het_100$TPM),mean(hom_del$TPM)))
  
  print(paste(transcript,t.test(tmp_wt$absolute_counts,het_inv$absolute_counts)$p.value,
              t.test(tmp_wt$absolute_counts,het_500$absolute_counts)$p.value,
              t.test(tmp_wt$absolute_counts,het_100$absolute_counts)$p.value,
              t.test(tmp_wt$absolute_counts,hom_del$absolute_counts)$p.value))
  print("-----------------")
}


#Make cell type plot
cleany_dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/cadm2_raw_counts.txt", sep='\t', header = T)
cleany_dat$sample <- gsub(".count","",cleany_dat$sample)
cleany_dat$celltype[grepl("iPSC",cleany_dat$sample)] <- "iPSC"
cleany_dat$celltype[grepl("iN",cleany_dat$sample)] <- "iN"

libSizeTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\data\\libSize.txt", header=TRUE, sep = '\t')#For human elp1 injected mouse reference
colnames(libSizeTable) <- c("sample","librarySize")
libSizeTable$sample <- gsub(".count", "", libSizeTable$sample)
celltype <- left_join(cleany_dat, libSizeTable)
celltype$normalized_counts <- (celltype$cleany/celltype$librarySize)*1E6

celltype$genotype[grepl("hom_del",celltype$sample)] <- "HOM_DEL"
celltype$genotype[grepl("500bp_het_del",celltype$sample)] <- "HET_DEL_500bp"
celltype$genotype[grepl("_100kb_het_del",celltype$sample)] <- "HET_DEL_100kb"
celltype$genotype[grepl("_500bp_inv",celltype$sample)] <- "HET_INV"
celltype$genotype[grepl("_wt",celltype$sample)] <- "WT"


celltype <- celltype %>%
    filter(genotype=="WT")


pdf("cell_type_expression.pdf")
ggplot(celltype, aes(x=celltype, y=normalized_counts, fill=celltype)) +
  #(paste0("iN_",transcript)) +
  geom_boxplot(width=0.7, position = 'dodge', outliers = NA) +
  #ylim(0,3) +
  xlab("Normalized Expression") +
  #geom_hline(yintercept=1) +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #geom_jitter(position = position_jitter(seed = 1)) +
  #geom_text(position = position_jitter(seed = 1), size=3.5) +
  theme(text = element_text(size = 19),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
