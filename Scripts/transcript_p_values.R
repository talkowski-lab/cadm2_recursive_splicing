#Making transcripts plots and calculating statistics
###############Building code for fold change matrix
library(tidyverse)
library(ggplot2)
dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/sumtwo.tsv", sep='\t', header=T)
colnames(dat) <- c("sample","transcript_id","TPM")
dat$sample <- gsub(" ", "", dat$sample)
dat$transcript_id <- gsub(" ", "", dat$transcript_id)


dat$transcript_id <- gsub("::3:85505218-86067982|::3:84958964-86067982|::3:84958964-86067982|::3:84958958-86067982|::3:85726515-86067982|::3:85802173-85848021|::3:85961626-86065626|::3:85541799-85543286", "", dat$transcript_id)
dat$celltype[grepl("iN_",dat$sample)] <- "iN"
dat$celltype[grepl("iPSC_",dat$sample)] <- "iPSC"

dat$genotype[grepl("hom_del",dat$sample)] <- "HOM 500bp DEL"
dat$genotype[grepl("500bp_het_del",dat$sample)] <- "HET 500bp DEL"
dat$genotype[grepl("_100kb_het_del",dat$sample)] <- "HET 100kb DEL"
dat$genotype[grepl("_500bp_inv",dat$sample)] <- "HET 500bp INV"
dat$genotype[grepl("_wt",dat$sample)] <- "WT"
dat$genotype[grepl("GM_Control",dat$sample)] <- "WT"

dat$genotype <- factor(dat$genotype, levels = c("WT","HET 500bp DEL","HET 100kb DEL","HET 500bp INV","HOM 500bp DEL"))

dat$transcript_id[grepl("CDS_Long\\+RS1\\+RS2",dat$transcript_id)] <- "CADM2+RS1+RS2"
dat$transcript_id[grepl("CDS_Long\\+RS2",dat$transcript_id)] <- "CAMD2+RS2"
dat$transcript_id[grepl("AltExon1_RS2_CDS_Short",dat$transcript_id)] <- "Alt Exon1+RS2+CADM2 Short"
dat$transcript_id[grepl("Novel_1",dat$transcript_id)] <- "Novel 1"
dat$transcript_id[grepl("Novel_2",dat$transcript_id)] <- "Novel 2"
dat$transcript_id[grepl("Novel_3",dat$transcript_id)] <- "Novel 3"
dat$transcript_id[grepl("CDS_Short",dat$transcript_id)] <- "CADM2 Short"
dat$transcript_id <- factor(dat$transcript_id, levels=c("CADM2+RS1+RS2",
                                                        "CAMD2+RS2",
                                                        "CADM2 Short",
                                                        "Alt Exon1+RS2+CADM2 Short",
                                                        "Novel 1",
                                                        "Novel 2",
                                                        "Novel 3" ))


iN_TPM <- dat %>%
  filter(celltype=="iN") %>%
  filter(!grepl("r1B6",sample)) %>%
  filter(transcript_id!= "Novel 2") %>%
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
hom_del <- new_table %>% filter(genotype =="HOM 500bp DEL")
inv_del <- new_table %>% filter(genotype =="HET 500bp INV")
het_del_100kb <- new_table %>% filter(genotype =="HET 100kb DEL")
het_del_500bp <- new_table %>% filter(genotype =="HET 500bp DEL")


##WT plots
tmp_wt <- read.table("C:\\Users\\ricar\\Downloads\\CADM2\\wt_tpm.tsv", header = F, sep = "\t")
df <- read.table("C:\\Users\\ricar\\OneDrive\\Documents\\cadm2_gene_tpm.txt", header = T,sep='\t')
colnames(tmp_wt)<-c("Sample","transcript_id","TPM_local")
tmp_wt <-left_join(tmp_wt,df)
#tmp_wt <- new_table %>% filter(genotype=="WT") %>% filter(!grepl("iN_total_GM_Control_1*",sample)) %>% filter(transcript_id!="CDS_Long")

tmp_wt$genotype <- "WT"
tmp_wt$transcript_id <- gsub("::3:85505218-86067982|::3:84958964-86067982|::3:84958964-86067982|::3:84958958-86067982|::3:85726515-86067982|::3:85802173-85848021|::3:85961626-86065626|::3:85541799-85543286", "", tmp_wt$transcript_id)
tmp_wt$transcript_id <- gsub("CDS_Long","CADM2",tmp_wt$transcript_id)
tmp_wt$transcript_id <- gsub("RS1","RE1",tmp_wt$transcript_id)
tmp_wt$transcript_id <- gsub("RS2","RE2",tmp_wt$transcript_id)
tmp_wt$transcript_id <- gsub("_"," ",tmp_wt$transcript_id)
tmp_wt$transcript_id <- gsub("AltExon1","Alt Exon 1+",tmp_wt$transcript_id)

tmp_wt$TPM_genomewide <- (tmp_wt$TPM_local/1e6)*tmp_wt$TPM

tmp_wt$transcript_id = factor(tmp_wt$transcript_id, levels=c("Alt Exon 1+ RE2 CDS Short","CADM2+RE2","Novel 3","CADM2+RE1+RE2","CDS Short", "Novel 1", "Novel 2"))

pdf("iN_iPSC_transcripts_counts_WT_no_pws_2.pdf", height = 6,width = 6)
ggerrorplot(tmp_wt, x = "transcript_id", y = "TPM_genomewide", 
            desc_stat = "mean_sd",
            color = "CellType",
            size = 1.4,
            legend="top",
            ylab="Normalized Counts Per Million",
            xlab="") +
  rotate_x_text() +
  facet_wrap(~CellType, scales="free_y")
dev.off()

pdf("iPSC_transcripts_counts_WT_no_pws_4.pdf", height = 6,width = 6)
ggplot(tmp_wt, aes(x=reorder(transcript_id,-TPM_genomewide), y=TPM_genomewide, fill=genotype)) +
  geom_boxplot(fill="white", outliers = F) +
  geom_hline(yintercept=1, linetype="dashed") +
  theme_classic() +
  ylab("TPM") +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 21),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+
#facet_wrap(~CellType, scales="free_y")
dev.off()










###Supporting scripts for analysis
library(tidyverse)
library(ggplot2)

dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/RSEM_results.TPM.tsv", sep='\t', header=T)
#########
#dat[dat$genotype=="PWS_WT",]$genotype <- "WT"
dat$genotype <- factor(dat$genotype, levels = c("WT","HET_DEL_500bp","HET_DEL_100kb","HET_INV","HOM_DEL"))
########

#dat$genotype <- factor(dat$genotype, levels = c("PWS_WT","WT","HET_DEL_500bp","HET_DEL_100kb","HET_INV","HOM_DEL"))
dat$transcript_id <- factor(dat$transcript_id, levels=c("CDS_Long",
                                                        "CDS_Long+RS1+RS2",
                                                        "CDS_Long+RS2",
                                                        "CDS_Short",
                                                        "AltExon1_RS2_CDS_Short",
                                                        "Novel_1",
                                                        "Novel_2",
                                                        "Novel_3"))

iN_transcripts <- dat %>%
  filter(celltype=="iN") %>%
  filter(!grepl("r1B6",sample)) %>%
  filter(genotype %in% c("PWS_WT", "WT"))

iN_transcripts_pws <- dat %>%
                  filter(celltype=="iN") %>%
                  filter(!grepl("_wt",sample)) %>%
                  filter(!grepl("r1B6",sample))

iN_transcripts_wt <- dat %>%
  filter(celltype=="iN") %>%
  filter(!grepl("GM_Control",sample)) %>%
  filter(!grepl("r1B6",sample))
  

iPSC_transcripts <- dat %>%
                    filter(celltype=="iPSC")



pdf("iN_transcripts.pdf", width=11, height=20)
ggplot(iN_transcripts, aes(x=transcript_id, y=TPM, fill=genotype)) +
  geom_boxplot() +
  ylab("TPM") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("iPSC_transcripts.pdf", width=11, height=20)
ggplot(iPSC_transcripts, aes(x=reorder(transcript_id,-TPM), y=TPM, fill=genotype, colour=genotype)) +
  geom_boxplot() +
  ylab("TPM") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


transcripts <- c("CDS_Long",
                 "CDS_Long+RS1+RS2",
                 "CDS_Long+RS2",
                 "CDS_Short",
                 "AltExon1_RS2_CDS_Short",
                 "Novel_1",
                 "Novel_2",
                 "Novel_3")
genotype <- c("HET_INV","HET_DEL_500bp", "HET_DEL_100kb","HOM_DEL")


for(x in transcripts){
  x=
temp_df <- iN_transcripts %>% filter(transcript_id ==x)



#if(t.test(wt$TPM,het_inv$TPM)$p.value < 0.05){print(paste(x,"HET_INV",t.test(wt$TPM,het_inv$TPM)$p.value))}
#if(t.test(wt$TPM,het_500$TPM)$p.value < 0.05){print(paste(x,"HET_500",t.test(wt$TPM,het_500$TPM)$p.value))}
#f(t.test(wt$TPM,het_100$TPM)$p.value < 0.05){print(paste(x,"HET_100",t.test(wt$TPM,het_100$TPM)$p.value))}
#if(t.test(wt$TPM,hom_del$TPM)$p.value < 0.05){print(paste(x,"HOM_DEL",t.test(wt$TPM,hom_del$TPM)$p.value))}
wt <-temp_df %>% filter(genotype == c("WT","PWS_WT"))
het_inv <- temp_df %>% filter(genotype=="HET_INV")
het_500 <- temp_df %>% filter(genotype=="HET_DEL_500bp")
het_100 <-temp_df %>% filter(genotype=="HET_DEL_100kb")
hom_del <-temp_df %>% filter(genotype=="HOM_DEL")

print(paste(x,t.test(wt$TPM,het_inv$TPM)$p.value,
            t.test(wt$TPM,het_500$TPM)$p.value,
            t.test(wt$TPM,het_100$TPM)$p.value,
            t.test(wt$TPM,hom_del$TPM)$p.value))
}


#cleanY_converter <- function(path){
#  cleany <- read.table(path)
#  return(2^t(cleany["ENSG00000175161",]))
}

iN_hom_cleany <-cleanY_converter("E://data/talkowski//Samples//cadm2//data//output//iN_WT_vs_HOM_DEL_cleany.txt")
iN_het_500_cleany <-cleanY_converter("E://data/talkowski//Samples//cadm2//data//output//iN_WT_vs_HET_500bp_DEL_cleany.txt")
iN_het_100_cleany <-cleanY_converter("E://data/talkowski//Samples//cadm2//data//output//iN_WT_vs_HET_100kb_DEL_cleany.txt")
iN_inv_cleany <-cleanY_converter("E://data/talkowski//Samples//cadm2//data//output//iN_WT_vs_HET_500bp_INV_cleany.txt")

df <- as.data.frame(rbind(iN_hom_cleany,iN_het_500_cleany,iN_het_100_cleany,iN_inv_cleany))
df$sample <- gsub(".count","",rownames(df))
new_df <- merge(df,dat)
new_df$counts <- (new_df$ENSG00000175161*new_df$TPM)/1e6

ting <- new_df %>% filter(Genotype != "WT") %>%
  filter(Genotype != "PWS_WT") %>%
  filter(Genotype %in% c("HET_DEL_100kb", "HET_INV"))


ting1 <- new_df %>% filter(Genotype != "WT") %>%
  filter(Genotype != "PWS_WT") %>%
  filter(Genotype %in% c("HET_INV", "HOM_DEL", "HET_DEL_500bp"))

pdf("iN_transcripts.pdf", width=11, height=20)
ggplot(ting, aes(x=reorder(transcript_id,counts), y=counts, fill=Genotype)) +
  geom_boxplot() +
  ylab("Normalized Counts") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("iN_transcripts.pdf", width=11, height=20)
ggplot(ting1, aes(x=reorder(transcript_id,counts), y=counts, fill=Genotype)) +
  geom_boxplot() +
  ylab("Normalized Counts") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

ggplot(new_df, aes(x=reorder(transcript_id,counts), y=counts, fill=Genotype)) +
  geom_boxplot() +
  ylab("Normalized Counts") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##################Figures with Expected counts
library(tidyverse)
library(ggplot2)

dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/new_recalibrated_TPM_expected_counts.txt", sep='\t', header=T)
dat$genotype <- factor(dat$genotype, levels = c("WT","HET_DEL_500bp","HET_DEL_100kb","HET_INV","HOM_DEL"))

dat$transcript_id <- factor(dat$transcript_id, levels=c("CDS_Long",
                                                        "CDS_Long+RS1+RS2",
                                                        "CDS_Long+RS2",
                                                        "CDS_Short",
                                                        "AltExon1_RS2_CDS_Short",
                                                        "Novel_1",
                                                        "Novel_2",
                                                        "Novel_3"))

iN_transcripts <- dat %>%
  filter(celltype=="iN") %>%
  filter(!grepl("r1B6",sample))

pdf("iN_transcripts.pdf", width=11, height=20)
ggplot(iN_transcripts, aes(x=transcript_id, y=expected_count, fill=genotype)) +
  geom_boxplot() +
  ylab("Expected Counts") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


iPSC_transcripts <- dat %>%
  filter(celltype=="iPSC")



pdf("iN_transcripts.pdf", width=11, height=20)
ggplot(iPSC_transcripts, aes(x=reorder(transcript_id,-expected_count), y=expected_count, fill=genotype)) +
  geom_boxplot() +
  ylab("Expected Counts") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#######new TPM and expected counts analysis

library(tidyverse)
library(ggplot2)
dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/new_recalibrated_TPM_expected_counts.txt", sep='\t', header=T)
colnames(dat) <- c("sample","transcript_id","expected_count", "TPM")
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


iN_transcripts <- dat %>%
  filter(celltype=="iN") %>%
  filter(!grepl("r1B6",sample))

cleany_dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/cleany_cadm2.txt", sep='\t', header = T)
cleany_dat$sample <- gsub(".count","",cleany_dat$sample)
cleany_dat$counts <- 2^cleany_dat$cleany
df_counts <- cleany_dat %>% group_by(sample) %>%  summarise(mean = counts %>% mean())

new_table <-left_join(iN_transcripts,df_counts)
new_table$absolute_counts <- (new_table$mean)*((new_table$TPM)/1e6)

new_table <- new_table %>% filter(genotype=="WT")

pdf("iN_transcripts_WT_new_expected_counts.pdf", width=11, height=20)
ggplot(new_table, aes(x=reorder(transcript_id,-expected_count), y=expected_count, fill=genotype)) +
  geom_boxplot() +
  ylab("Recalibrated Expected Counts") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("iN_transcripts_WT_new_TPM.pdf", width=11, height=20)
ggplot(new_table, aes(x=reorder(transcript_id,-TPM), y=TPM, fill=genotype)) +
  geom_boxplot() +
  ylab("Recalibrated TPM") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("iN_transcripts_WT_new_absolute_counts.pdf", width=11, height=20)
ggplot(new_table, aes(x=reorder(transcript_id,-absolute_counts), y=absolute_counts, fill=genotype)) +
  geom_boxplot() +
  ylab("Recalibrated Absolute Counts") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
###TPM
library(tidyverse)
library(ggplot2)
dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/new_RSEM_TPM.tsv", sep='\t')
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

dat[is.na(dat)] <-0

iN_TPM <- dat %>%
  filter(celltype=="iN") %>%
  filter(!grepl("r1B6",sample))

pdf("iN_transcripts_new_expected_counts.pdf", width=11, height=20)
ggplot(iN_TPM, aes(x=transcript_id, y=new_TPM, fill=genotype, label=sample)) +
  geom_boxplot() +
  ylab("TPM") +
  theme_classic() +
  geom_jitter(position = position_jitter(seed = 1)) +
  ggrepel::geom_text_repel(position = position_jitter(seed = 1)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

CDS_Long_new_table <- new_table %>% filter(transcript_id=="CDS_Long")
CDS_Long_  
CDS_Long <- iN_TPM %>% filter(transcript_id=="CDS_Long")
ggplot(CDS_Long, aes(x=genotype, y=new_TPM, fill=genotype, label=sample)) +
      geom_boxplot() +
      ylab("TPM") +
      theme_classic() +
      geom_jitter(position = position_jitter(seed = 1)) +
     geom_text(position = position_jitter(seed = 1)) +
      #theme(legend.position = "none") +
      theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

iPSC_TPM <- dat %>%
  filter(celltype=="iPSC")

pdf("iPSC_transcripts_new_expected_counts.pdf", width=11, height=20)
ggplot(iPSC_TPM, aes(x=transcript_id, y=order(-new_TPM), fill=genotype)) +
  geom_boxplot() +
  ylab("TPM") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#########
aggregate(x$Frequency, by=list(Category=x$Category), FUN=sum)

aggregate(iN_transcripts$new_expected_count, by=list(Category=iN_transcripts$genptype), FUN=sum)
test <- aggregate(new_expected_count ~ genotype, iN_transcripts,sum)

test <- aggregate(. ~ genotype + sample, data=iN_transcripts,sum)

library(dplyr)
new_mtx <- iN_transcripts %>% group_by(genotype,sample) %>% summarize(V=sum(new_expected_count))
ggplot(new_mtx, aes(x=genotype, y=V, fill=genotype)) +
  geom_boxplot()

###############Building code for fold change matrix
library(tidyverse)
library(ggplot2)
dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/sumtwo.tsv", sep='\t', header=T)
colnames(dat) <- c("sample","transcript_id","TPM")
dat$sample <- gsub(" ", "", dat$sample)
dat$transcript_id <- gsub(" ", "", dat$transcript_id)


dat$transcript_id <- gsub("::3:85505218-86067982|::3:84958964-86067982|::3:84958964-86067982|::3:84958958-86067982|::3:85726515-86067982|::3:85802173-85848021|::3:85961626-86065626|::3:85541799-85543286", "", dat$transcript_id)
dat$celltype[grepl("iN_",dat$sample)] <- "iN"
dat$celltype[grepl("iPSC_",dat$sample)] <- "iPSC"

dat$genotype[grepl("hom_del",dat$sample)] <- "HOM 500bp DEL"
dat$genotype[grepl("500bp_het_del",dat$sample)] <- "HET 500bp DEL"
dat$genotype[grepl("_100kb_het_del",dat$sample)] <- "HET 100kb DEL"
dat$genotype[grepl("_500bp_inv",dat$sample)] <- "HET 500bp INV"
dat$genotype[grepl("_wt",dat$sample)] <- "WT"
dat$genotype[grepl("GM_Control",dat$sample)] <- "WT"

dat$genotype <- factor(dat$genotype, levels = c("WT","HET 500bp DEL","HET 100kb DEL","HET 500bp INV","HOM 500bp DEL"))

dat$transcript_id[grepl("CDS_Long\\+RS1\\+RS2",dat$transcript_id)] <- "CADM2+RS1+RS2"
dat$transcript_id[grepl("CDS_Long\\+RS2",dat$transcript_id)] <- "CAMD2+RS2"
dat$transcript_id[grepl("AltExon1_RS2_CDS_Short",dat$transcript_id)] <- "Alt Exon1+RS2+CADM2 Short"
dat$transcript_id[grepl("Novel_1",dat$transcript_id)] <- "Novel 1"
dat$transcript_id[grepl("Novel_2",dat$transcript_id)] <- "Novel 2"
dat$transcript_id[grepl("Novel_3",dat$transcript_id)] <- "Novel 3"
dat$transcript_id[grepl("CDS_Short",dat$transcript_id)] <- "CADM2 Short"
dat$transcript_id <- factor(dat$transcript_id, levels=c("CADM2+RS1+RS2",
                                                        "CAMD2+RS2",
                                                        "CADM2 Short",
                                                        "Alt Exon1+RS2+CADM2 Short",
                                                        "Novel 1",
                                                        "Novel 2",
                                                        "Novel 3" ))


iN_TPM <- dat %>%
  filter(celltype=="iN") %>%
  filter(!grepl("r1B6",sample)) %>%
  filter(transcript_id!= "Novel 2") %>%
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
hom_del <- new_table %>% filter(genotype =="HOM 500bp DEL")
inv_del <- new_table %>% filter(genotype =="HET 500bp INV")
het_del_100kb <- new_table %>% filter(genotype =="HET 100kb DEL")
het_del_500bp <- new_table %>% filter(genotype =="HET 500bp DEL")

iN_counts <- new_table %>% filter(celltype=="iN") %>% group_by(genotype, sample) %>% summarise(sum = sum(absolute_counts, na.rm = TRUE))
iN_FC <- new_table %>% filter(celltype=="iN") %>% group_by(genotype, sample) %>% summarise(mean = mean(FC, na.rm = TRUE))

pdf("iN_counts_transcript_drop_raw_counts.pdf")
#Plot
ggplot(iN_counts, aes(x=genotype, y=sum, fill=genotype)) +
  geom_boxplot() +
  ylab("Corrected Counts") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

iN_FC <- iN_FC %>% filter(genotype!="WT")
pdf("fc_transcript_drop_raw_counts.pdf")
ggplot(iN_FC, aes(x=genotype, y=mean, fill=genotype)) + geom_boxplot() +
  theme_classic() +
  ylab('Fold Change') +
  geom_hline(yintercept =1) +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

new_table$zygosity <- new_table$genotype
new_table$zygosity <- gsub("HET 500bp INV", "HET",new_table$zygosity)
new_table$zygosity <- gsub("HET 100kb DEL", "HET",new_table$zygosity)
new_table$zygosity <- gsub("HET 500bp DEL", "HET",new_table$zygosity)
new_table$zygosity <- gsub("HOM 500bp DEL", "HOM",new_table$zygosity)
new_table$zygosity = factor(new_table$zygosity, levels=c("WT","HET","HOM"))
new_table <- new_table %>% filter(sample!="iN_rC5_rep2_500bp_inv")


for (transcript in unique(new_table$transcript_id)){
  tmp <- new_table %>% filter(transcript_id==transcript)
  tmp_wt <- tmp %>% filter(genotype=="WT")
  tmp$FC<- tmp$absolute_counts/median(tmp_wt$absolute_counts)
  tmp_wt <- tmp %>% filter(genotype=="WT")
  het_inv <- tmp %>% filter(genotype=="HET 500bp INV")
  het_500 <- tmp %>% filter(genotype=="HET 500bp DEL")
  het_100 <-tmp %>% filter(genotype=="HET 100kb DEL")
  hom_del <-tmp %>% filter(genotype=="HOM 500bp DEL")

  
  #tmp_plot <- ggerrorplot(tmp, x = "genotype", y = "FC", 
  #            desc_stat = "mean_sd",
  #            color = "zygosity",
  #            size = 1.4,
  #            legend="center",
  #            ylab="Fold Change",
  #            xlab="", 
  #            ylim=c(0,2)) +
  #  rotate_x_text()
  #tmp_plot <- ggplot(tmp, aes(x=genotype, y=FC, fill=zygosity)) +
  #       ggtitle(paste0("iN ",transcript)) +
  #        geom_boxplot(width=0.7, position = 'dodge', outliers = FALSE) +
  #  ylim(0,3) +
  # xlab("FC") +
  #        geom_hline(yintercept=1) +
  #        theme_classic() +
  #  geom_jitter(position=position_jitterdodge(0.5)) +#

  #        theme(text = element_text(size = 19),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #ggsave(tmp_plot,file=paste0("iN_final_",transcript,"_scale.pdf"), width = 4, height = 4)
  
  #print(paste(transcript, mean(tmp_wt$FC),median(tmp_wt$FC)))
 #print(paste(transcript,t.test(tmp_wt$FC,het_inv$FC)$p.value,
#                t.test(tmp_wt$FC,het_500$FC)$p.value,
#               t.test(tmp_wt$FC,het_100$FC)$p.value,
#               t.test(tmp_wt$FC,hom_del$FC)$p.value))
print(paste(transcript,mean(tmp_wt$FC),mean(het_inv$FC),mean(het_500$FC),mean(het_100$FC),mean(hom_del$FC)))
 #print(paste(transcript,mean(tmp_wt$absolute_counts),mean(het_inv$absolute_counts),mean(het_500$absolute_counts),mean(het_100$absolute_counts),mean(hom_del$absolute_counts)))
 #print(paste(transcript,mean(tmp_wt$TPM),mean(het_inv$TPM),mean(het_500$TPM),mean(het_100$TPM),mean(hom_del$TPM)))
  
   #print(paste(transcript,t.test(tmp_wt$absolute_counts,het_inv$absolute_counts)$p.value,
  #             t.test(tmp_wt$absolute_counts,het_500$absolute_counts)$p.value,
  #             t.test(tmp_wt$absolute_counts,het_100$absolute_counts)$p.value,
  #             t.test(tmp_wt$absolute_counts,hom_del$absolute_counts)$p.value))
  #print("-----------------")
}
write.table(new_table, "iN_transcript_stats_no_pws_drop_transcripts_raw_counts.txt", quote = F,row.names = F, col.names = T)

iN_counts <- iN_counts[!grepl("iN_total_GM_Control_1_", iN_counts$sample),,drop = FALSE]

tmp_wt <-  iN_counts %>% filter(genotype=="WT")
het_inv <- iN_counts %>% filter(genotype=="HET_INV")
het_500 <- iN_counts %>% filter(genotype=="HET_DEL_500bp")
het_100 <-iN_counts %>% filter(genotype=="HET_DEL_100kb")
hom_del <-iN_counts %>% filter(genotype=="HOM_DEL")

print(paste(t.test(tmp_wt$sum,het_inv$sum)$p.value,
t.test(tmp_wt$sum,het_500$sum)$p.value,
t.test(tmp_wt$sum,het_100$sum)$p.value,
t.test(tmp_wt$sum,hom_del$sum)$p.value))
#####
iN_FC <- iN_FC[!grepl("iN_total_GM_Control_1_", iN_FC$sample),,drop = FALSE]
tmp_wt <-  iN_FC %>% filter(genotype=="WT")
het_inv <- iN_FC %>% filter(genotype=="HET_INV")
het_500 <- iN_FC %>% filter(genotype=="HET_DEL_500bp")
het_100 <-iN_FC %>% filter(genotype=="HET_DEL_100kb")
hom_del <-iN_FC %>% filter(genotype=="HOM_DEL")

print(paste(t.test(tmp_wt$mean,het_inv$mean)$p.value,
            t.test(tmp_wt$mean,het_500$mean)$p.value,
            t.test(tmp_wt$mean,het_100$mean)$p.value,
            t.test(tmp_wt$mean,hom_del$mean)$p.value))
########QC to add transcripts together or from samples.
df_transcript <- new_table %>% group_by(genotype) %>% summarise(mean(normalized_counts))
df_transcript$genotype[grepl("hom_del",df_transcript$sample)] <- "HOM_DEL"
df_transcript$genotype[grepl("500bp_het_del",df_transcript$sample)] <- "HET_DEL_500bp"
df_transcript$genotype[grepl("_100kb_het_del",df_transcript$sample)] <- "HET_DEL_100kb"
df_transcript$genotype[grepl("_500bp_inv",df_transcript$sample)] <- "HET_INV"
df_transcript$genotype[grepl("_wt",df_transcript$sample)] <- "WT"

pdf("combined_transcript_genotype_normalized_expression.pdf")
ggplot(df_transcript, aes(x=genotype, y=`mean(normalized_counts)`, fill=genotype)) +
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

##Select the most abundant transcripts
df_transcript <- new_table %>% group_by(genotype) %>% summarise(mean(normalized_counts))


df_transcript_Alt_Exon <- new_table %>% filter(transcript_id=="AltExon1_RS2_CDS_Short") %>% group_by(genotype) %>% summarise(mean(normalized_counts))
df_transcript_RS1_RS2 <- new_table %>% filter(transcript_id=="CDS_Long+RS1+RS2") %>% group_by(genotype) %>% summarise(mean(normalized_counts))
df_transcript_RS2  <- new_table %>% filter(transcript_id=="CDS_Long+RS2") %>% group_by(genotype) %>% summarise(mean(normalized_counts))
df_transcript_short  <- new_table %>% filter(transcript_id=="CDS_Short") %>% group_by(genotype) %>% summarise(mean(normalized_counts))
df_transcript_novel_3  <- new_table %>% filter(transcript_id=="Novel_3") %>% group_by(genotype) %>% summarise(mean(normalized_counts))

pdf("raw_normalized_transcript_Novel_3.pdf")
ggplot(df_transcript_RS1_RS2, aes(x=genotype, y=`mean(normalized_counts)`, fill=genotype)) +
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


##WT plots
tmp_wt <- read.table("C:\\Users\\ricar\\Downloads\\CADM2\\wt_tpm.tsv", header = F, sep = "\t")
df <- read.table("C:\\Users\\ricar\\OneDrive\\Documents\\cadm2_gene_tpm.txt", header = T,sep='\t')
colnames(tmp_wt)<-c("Sample","transcript_id","TPM_local")
tmp_wt <-left_join(tmp_wt,df)
#tmp_wt <- new_table %>% filter(genotype=="WT") %>% filter(!grepl("iN_total_GM_Control_1*",sample)) %>% filter(transcript_id!="CDS_Long")

tmp_wt$genotype <- "WT"
tmp_wt$transcript_id <- gsub("::3:85505218-86067982|::3:84958964-86067982|::3:84958964-86067982|::3:84958958-86067982|::3:85726515-86067982|::3:85802173-85848021|::3:85961626-86065626|::3:85541799-85543286", "", tmp_wt$transcript_id)
tmp_wt$transcript_id <- gsub("CDS_Long","CADM2",tmp_wt$transcript_id)
tmp_wt$transcript_id <- gsub("RS1","RE1",tmp_wt$transcript_id)
tmp_wt$transcript_id <- gsub("RS2","RE2",tmp_wt$transcript_id)
tmp_wt$transcript_id <- gsub("_"," ",tmp_wt$transcript_id)
tmp_wt$transcript_id <- gsub("AltExon1","Alt Exon 1+",tmp_wt$transcript_id)

tmp_wt$TPM_genomewide <- (tmp_wt$TPM_local/1e6)*tmp_wt$TPM

tmp_wt$transcript_id = factor(tmp_wt$transcript_id, levels=c("Alt Exon 1+ RE2 CDS Short","CADM2+RE2","Novel 3","CADM2+RE1+RE2","CDS Short", "Novel 1", "Novel 2"))

pdf("iN_iPSC_transcripts_counts_WT_no_pws_2.pdf", height = 6,width = 6)
ggerrorplot(tmp_wt, x = "transcript_id", y = "TPM_genomewide", 
            desc_stat = "mean_sd",
            color = "CellType",
            size = 1.4,
            legend="top",
            ylab="Normalized Counts Per Million",
            xlab="") +
  rotate_x_text() +
  facet_wrap(~CellType, scales="free_y")
dev.off()

pdf("iPSC_transcripts_counts_WT_no_pws_4.pdf", height = 6,width = 6)
ggplot(tmp_wt, aes(x=reorder(transcript_id,-TPM_genomewide), y=TPM_genomewide, fill=genotype)) +
       geom_boxplot(fill="white", outliers = F) +
       geom_hline(yintercept=1, linetype="dashed") +
       theme_classic() +
       ylab("TPM") +
       geom_jitter(position=position_jitterdodge(0.5)) +
       theme(text = element_text(size = 21),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+
 #facet_wrap(~CellType, scales="free_y")
dev.off()

geom_boxplot(width=0.7, position = 'dodge', outliers = F) +
  ylab("Fold Change") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tmp_wt %>% filter(transcript_id=="Alt Exon 1+ RE2 CDS Short") %>% summarise(mean(TPM_genomewide))
tmp_wt %>% filter(transcript_id=="CADM2+RE1+RE2") %>% summarise(mean(TPM_genomewide))
tmp_wt %>% filter(transcript_id=="CADM2+RE2") %>% summarise(mean(TPM_genomewide))
tmp_wt %>% filter(transcript_id=="Novel 1") %>% summarise(mean(TPM_genomewide))
tmp_wt %>% filter(transcript_id=="Novel 2") %>% summarise(mean(TPM_genomewide))
tmp_wt %>% filter(transcript_id=="Novel 3") %>% summarise(mean(TPM_genomewide))
tmp_wt %>% filter(transcript_id=="CDS Short") %>% summarise(mean(TPM_genomewide))



tmp <- thing %>% filter(transcript_id==transcript)
tmp_wt <- tmp %>% filter(genotype=="WT")
tmp$FC<- tmp$absolute_counts/median(tmp_wt$absolute_counts)
tmp_wt <- tmp %>% filter(genotype=="WT")
het_inv <- tmp %>% filter(genotype=="HET_INV")
het_500 <- tmp %>% filter(genotype=="HET_DEL_500bp")
het_100 <-tmp %>% filter(genotype=="HET_DEL_100kb")
hom_del <-tmp %>% filter(genotype=="HOM_DEL")

tmp_plot <- ggplot(tmp, aes(x=genotype, y=absolute_counts, colour=genotype)) +
  ggtitle(paste0("iN_",transcript)) +
  geom_point() +
  xlab("Fold Change") +
  #geom_hline(yintercept=1) +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(tmp_plot,file=paste0("iN_",transcript, "_point_counts",".pdf"), width=6, height = 5)

##################

CDS_Long_new_table <- new_table %>% filter(transcript_id=="CDS_Long") %>% filter(genotype=="WT")
CDS_Long <- iN_TPM %>% filter(transcript_id=="CDS_Long")


cadm2_expression <- read.table("C://Users//ricar//Downloads//CADM2//cadm2_expression.txt", sep="\t", header=T)

df <- left_join(CDS_Long_new_table, cadm2_expression)

pdf("iN_cleany_combined_absolute_counts_comparison.pdf")
print(ggscatter(df, x = "absolute_counts", y = "cleany_cadm2",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regression line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
))
dev.off()

library(ggcorrplot)
CDS_Long_new_table <- df[,c("cleany","new_TPM","absolute_counts")]
corr <- round(cor(CDS_Long_new_table ), 1)
pdf("cadm2_corrplot.pdf")
ggcorrplot(corr, hc.order = TRUE,
           type = "lower", lab = T)
dev.off()


tmp <- new_table %>% filter(transcript_id=="CDS_Long") %>% filter(genotype)

ggplot(tmp, aes(x=genotype, y=cleany, label=sample)) + geom_boxplot() +
  theme_classic() +
  geom_jitter(position = position_jitter(seed = 1)) +
  ggrepel::geom_text_repel(position = position_jitter(seed = 1)) +
  theme(text = element_text(size = 12),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

###########Make transcript figure with individual genotypes
pdf("iN_transcript_genotypes_counts.pdf")
ggplot(new_table, aes(x=reorder(transcript_id,-absolute_counts), y=absolute_counts, fill=genotype)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Recalibrated Counts") +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


#####calculate percentage of transcripts expressed in each genotype

thing <- new_table %>% filter(transcript_id!="CDS_Long" & transcript_id!="Novel_2")
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
