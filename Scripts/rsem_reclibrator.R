library(tidyverse)

setwd("D://PHShome//rsh38//rsem_output")

rsem_dirs <- dir("D://PHShome//rsh38//rsem_output", pattern="i")
sj_file_dir <- "D://data//talkowski//Samples//cadm2//data//Alignment_CADM2_RNAseq_EditedGTF//"

fileName <- file.path(rsem_dirs, "trinity_DN_rsem01/RSEM.isoforms.results")

for (i in seq(along=fileName)){
  print(fileName[i])
  olddata <- read.table(fileName[i], header=TRUE)
  j1 <- read.table(paste0(sj_file_dir,rsem_dirs[i],'//',rsem_dirs[i],".SJ.out.tab")) %>% filter(V2 ==84959669) %>% filter(V3 ==85238956) %>% select(V7) %>% as.integer()
  j2 <- read.table(paste0(sj_file_dir,rsem_dirs[i],'//',rsem_dirs[i],".SJ.out.tab")) %>% filter(V2 ==84959669) %>% filter(V3 ==85511856) %>% select(V7) %>% as.integer()
  #j3 <- read.table(paste0(sj_file_dir,rsem_dirs[i],'//',rsem_dirs[i],".SJ.out.tab")) %>% filter(V2 ==84959669) %>% filter(V3 ==85726521) %>% select(V7) %>% as.integer()
  if(is.na(j1)){j1<-0}
  if(is.na(j2)){j2<-0}
  #if(is.na(j3)){j3<-0}
  #ratio1 <- j1/(j1+j2+j3)
  #ratio2 <- j2/(j1+j2+j3)
  #ratio3 <- j3/(j1+j2+j3)
  ratio1 <- j1/(j1+j2)
  ratio2 <- j2/(j1+j2)
  olddata$relative_ratio <- 1
  olddata[olddata$transcript_id=="CDS_Long+RS1+RS2::3:84958964-86067982",]$relative_ratio <- ratio1
  olddata[olddata$transcript_id=="CDS_Long+RS2::3:84958964-86067982",]$relative_ratio <- ratio2
  #olddata[olddata$transcript_id=="CDS_Long::3:84958958-86067982",]$relative_ratio <- ratio3
  #if((j1+j2+j3)==0){
  if((j1+j2)==0){
    print(fileName[i])
    olddata$new_expected_count <- olddata$expected_count
  } else {
    olddata$new_expected_count <- olddata$expected_count
    #sum_three <- (olddata[olddata$transcript_id=="CDS_Long+RS1+RS2::3:84958964-86067982","expected_count"] + olddata[olddata$transcript_id=="CDS_Long+RS2::3:84958964-86067982","expected_count"] + olddata[olddata$transcript_id=="CDS_Long::3:84958958-86067982","expected_count"])
    sum_two <- (olddata[olddata$transcript_id=="CDS_Long+RS1+RS2::3:84958964-86067982","expected_count"] + olddata[olddata$transcript_id=="CDS_Long+RS2::3:84958964-86067982","expected_count"])
    #new_expected_counts_tx1 <- sum_three * ratio1
    #new_expected_counts_tx2 <- sum_three * ratio2
    #new_expected_counts_tx3 <- sum_three * ratio3
    
    new_expected_counts_tx1 <- sum_two * ratio1
    new_expected_counts_tx2 <- sum_two * ratio2

    olddata[olddata$transcript_id=="CDS_Long+RS1+RS2::3:84958964-86067982","new_expected_count"] <- new_expected_counts_tx1
    olddata[olddata$transcript_id=="CDS_Long+RS2::3:84958964-86067982","new_expected_count"]<- new_expected_counts_tx2
    #olddata[olddata$transcript_id=="CDS_Long::3:84958958-86067982","new_expected_count"]<- new_expected_counts_tx3

  }
  olddata$New_RPK <- (olddata$new_expected_count*1000)/olddata$effective_length
  new_scaling_factor <- sum(olddata$New_RPK)/1e6
  olddata$new_TPM <- (olddata$New_RPK)/new_scaling_factor
  olddata$Ex1Rs1 <- j1
  olddata$Ex1Rs2 <- j2
  #olddata$Ex1Ex2 <- j3
  new_name <- gsub("RSEM.isoforms.results","sumtwo_RSEM_recalibrated_3.isoforms.results",fileName[i])
  write.table(olddata,new_name,quote = F,row.names = F, sep = "\t")
}

####Plot junctions
dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/cadm2_transcript_junction_counts.txt", sep='\t', header=T)
libSizeTable <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/libSize.txt", header=TRUE, sep = '\t')
colnames(libSizeTable) <- c("sample","librarySize")
libSizeTable$sample <- gsub(".count","", libSizeTable$sample)
dat <- right_join(dat,libSizeTable)
dat$transcript_id <- gsub("::3:85505218-86067982|::3:84958964-86067982|::3:84958964-86067982|::3:84958958-86067982|::3:85726515-86067982|::3:85802173-85848021|::3:85961626-86065626|::3:85541799-85543286", "", dat$transcript_id)
dat$celltype[grepl("iN_",dat$sample)] <- "iN"
dat$celltype[grepl("iPSC_",dat$sample)] <- "iPSC"

dat$genotype[grepl("hom_del",dat$sample)] <- "HOM_DEL"
dat$genotype[grepl("500bp_het_del",dat$sample)] <- "HET_DEL_500bp"
dat$genotype[grepl("_100kb_het_del",dat$sample)] <- "HET_DEL_100kb"
dat$genotype[grepl("_500bp_inv",dat$sample)] <- "HET_INV"
dat$genotype[grepl("_wt",dat$sample)] <- "WT"
dat$genotype[grepl("GM_Control",dat$sample)] <- "WT"
dat$Ex1Rs1 <- dat$Ex1Rs1/dat$librarySize
dat$Ex1Rs2 <- dat$Ex1Rs2/dat$librarySize
dat$Ex1Ex2 <- dat$Ex1Ex2/dat$librarySize

dat$genotype <- factor(dat$genotype, levels = c("WT","HET_DEL_500bp","HET_DEL_100kb","HET_INV","HOM_DEL"))
dat$transcript_id <- factor(dat$transcript_id, levels=c("CDS_Long",
                                                        "CDS_Long+RS1+RS2",
                                                        "CDS_Long+RS2",
                                                        "CDS_Short",
                                                        "AltExon1_RS2_CDS_Short",
                                                        "Novel_1",
                                                        "Novel_2",
                                                        "Novel_3"))

iN_TPM <- dat %>%
  filter(celltype=="iN") %>%
  filter(!grepl("r1B6",sample))

library(reshape2)

iN_dat <- iN_TPM[c("genotype", "Ex1Rs1", "Ex1Rs2", "Ex1Ex2")]
new_dat <- melt(iN_dat,id='genotype')

library(ggplot2)
pdf("C://Documents and Settings/ricar/Downloads/CADM2/junction_genotype.pdf")
  ggplot(new_dat, aes(x=genotype, y=value, fill=variable)) +
  geom_boxplot() +
  ylab("Normalized Counts") +
  scale_fill_discrete(name = "Junction", labels = c("Ex1-RS1", "Ex1-RS2", "Ex1-Ex2")) +
  theme_classic() +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

dat1 <- read.table("C:\\Users\\ricar\\Downloads\\CADM2\\final_junction.txt",header=T, sep = "\t")
library(reshape)
dat_wide <- reshape(dat, idvar = "sample", timevar = "junction", direction = "wide")
dat_table <- dat_wide[,c("sample","unique_map.RS1RS2","unique_map.Ex1Ex2","unique_map.Ex1RS1","unique_map.Ex1RS2","genotype.Ex1RS2" )]
dat_table[is.na(dat_table)] <- 0
write.table(dat_table, "cadm2_junction_counts.tsv", sep = "\t", row.names = F, quote = F)


#############Plot Ex1-RS1
library(tidyverse)
dat1 <- read.table("C:\\Users\\ricar\\Downloads\\CADM2\\final_junction_RNAseq.txt",header=T, sep = "\t")
libSizeTable <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/libSize.txt", header=TRUE, sep = '\t')
colnames(libSizeTable) <- c("sample","librarySize")
libSizeTable$sample <- gsub(".count","", libSizeTable$sample)
dat1$celltype[grepl("iN_",dat1$sample)] <- "iN"
dat1 <- left_join(dat1,libSizeTable)
dat1$genotype <- factor(dat1$genotype, levels=c("WT",
                                                        "INV",
                                                        "HET_500bp",
                                                        "HET_100kb",
                                                        "HOM_DEL"))
dat1 <- dat1 %>% filter(junction=="Ex7Ex8")
dat1 <- dat1[1:20,]
dat1$CPM <- (dat1$unique_map/dat1$librarySize)*1e6

library(ggplot2)
pdf("C://Documents and Settings/ricar/Downloads/CADM2/Ex1_Ex8_junction_counts.pdf")
ggplot(dat1, aes(x=genotype, y=unique_map, fill=genotype)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Counts") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

het_100 <- dat1 %>% filter(genotype=="HET_100kb")
het_500 <- dat1 %>% filter(genotype=="HET_500bp")
hom <-dat1 %>% filter(genotype=="HOM_DEL")
inv <- dat1 %>% filter(genotype=="INV")
wt <- dat1 %>% filter(genotype=="WT")

###% Reduction
((mean(wt$CPM)-mean(hom$CPM))/mean(wt$CPM))*100
((mean(wt$CPM)-mean(het_100$CPM))/mean(wt$CPM))*100
((mean(wt$CPM)-mean(het_500$CPM))/mean(wt$CPM))*100
((mean(wt$CPM)-mean(inv$CPM))/mean(wt$CPM))*100

###

t.test(hom$CPM,wt$CPM)$p.value
t.test(het_100$CPM,wt$CPM)$p.value
t.test(het_500$CPM,wt$CPM)$p.value
t.test(inv$CPM,wt$CPM)$p.value


########Counts
###% Reduction
((mean(wt$unique_map)-mean(hom$unique_map))/mean(wt$unique_map))*100
((mean(wt$unique_map)-mean(het_100$unique_map))/mean(wt$unique_map))*100
((mean(wt$unique_map)-mean(het_500$unique_map))/mean(wt$unique_map))*100
((mean(wt$unique_map)-mean(inv$unique_map))/mean(wt$unique_map))*100

###

t.test(hom$unique_map,wt$unique_map)$p.value
t.test(het_100$unique_map,wt$unique_map)$p.value
t.test(het_500$unique_map,wt$unique_map)$p.value
t.test(inv$unique_map,wt$unique_map)$p.value
