library(tidyverse)
library(ggplot2)
dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/new_RSEM_recalibrated.txt", sep='\t')
colnames(dat) <- c("sample","transcript_id","new_TPM")
dat$sample <- gsub("/", "", dat$sample)


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
dat$transcript_id <- factor(dat$transcript_id, levels=c("CDS_Long",
                                                        "CDS_Long+RS1+RS2",
                                                        "CDS_Long+RS2",
                                                        "CDS_Short",
                                                        "AltExon1_RS2_CDS_Short",
                                                        "Novel_1",
                                                        "Novel_2",
                                                        "Novel_3"))

iPSC_TPM <- dat %>%
  filter(celltype=="iPSC") 

cleany_dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/cleany_iPSC_dat.txt", sep='\t', header = T)
cleany_dat$sample <- gsub(".count","",cleany_dat$sample)

df_list <- list(cleany_dat,iPSC_TPM)

new_table <- df_list %>% reduce(left_join,by='sample')

new_table$absolute_counts <- (2^(new_table$cleany))*((new_table$new_TPM)/1e6)
#new_table[is.na(new_table)] <-0


wt <- new_table %>% filter(genotype=="WT")
hom_del <- new_table %>% filter(genotype =="HOM_DEL")
inv_del <- new_table %>% filter(genotype =="HET_INV")
het_del_100kb <- new_table %>% filter(genotype =="HET_DEL_100kb")
het_del_500bp <- new_table %>% filter(genotype =="HET_DEL_500bp")

#iPSC_counts <- new_table %>% filter(celltype=="iPSC") %>% group_by(genotype, sample) %>% summarise(sum = sum(absolute_counts, na.rm = TRUE))
#iPSC_FC <- new_table %>% filter(celltype=="iPSC") %>% group_by(genotype, sample) %>% summarise(mean = mean(FC, na.rm = TRUE))
write.table(new_table, "iPSC_transcript_stats.txt", quote = F,row.names = F, col.names = T)


thing <- new_table %>% filter(transcript_id!="CDS_Long")
for (transcript in unique(thing$transcript_id)){
  tmp <- thing %>% filter(transcript_id==transcript)
  tmp_wt <- tmp %>% filter(genotype=="WT")
  tmp$FC<- tmp$absolute_counts/median(tmp_wt$absolute_counts)
  tmp_wt <- tmp %>% filter(genotype=="WT")
  het_inv <- tmp %>% filter(genotype=="HET_INV")
  het_500 <- tmp %>% filter(genotype=="HET_DEL_500bp")
  het_100 <-tmp %>% filter(genotype=="HET_DEL_100kb")
  hom_del <-tmp %>% filter(genotype=="HOM_DEL")
  
  tmp_plot <- ggplot(tmp, aes(x=genotype, y=FC, fill=genotype)) +
  ggtitle(paste0("iPSC_",transcript)) +
          xlab("Fold Change") +
    geom_boxplot() +
    geom_hline(yintercept=1) +
    theme_classic() +
    geom_jitter(position=position_jitterdodge(0.5)) +
    theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(tmp_plot,file=paste0("iPSC_",transcript,".pdf"), width=5, height = 5)
  print(paste(transcript,t.test(tmp_wt$FC,het_inv$FC)$p.value,
              t.test(tmp_wt$FC,het_500$FC)$p.value,
              t.test(tmp_wt$FC,het_100$FC)$p.value,
              t.test(tmp_wt$FC,hom_del$FC)$p.value))
  
  #print(paste(transcript,t.test(tmp_wt$absolute_counts,het_inv$absolute_counts)$p.value,
  #           t.test(tmp_wt$absolute_counts,het_500$absolute_counts)$p.value,
  #          t.test(tmp_wt$absolute_counts,het_100$absolute_counts)$p.value,
  #         t.test(tmp_wt$absolute_counts,hom_del$absolute_counts)$p.value))
  
  print("---------------------")
}

tmp_wt <- new_table %>% filter(genotype=="WT")


pdf("iPSC_transcript_genotypes_absolute.pdf")
ggplot(new_table, aes(x=reorder(transcript_id,-absolute_counts), y=absolute_counts, fill=genotype)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Recalibrated Counts")+
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("iPSC_transcripts_TPM.pdf")
ggplot(new_table, aes(x=reorder(transcript_id,-new_TPM), y=new_TPM, fill=genotype)) +
  geom_boxplot(fill="white") +
  ylab("Recalibrated Counts") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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
  
  a <-tmp_wt$new_TPM>0
  b <-het_inv$new_TPM>0
  c <-het_500$new_TPM>0
  d <-het_100$new_TPM>0
  e <-hom_del$new_TPM>0
  print(paste(transcript,length(a[a==TRUE])/length(a),
              length(b[b==TRUE])/length(b),
              length(c[c==TRUE])/length(c),
              length(d[d==TRUE])/length(d),
              length(e[e==TRUE])/length(e)))
  
}

