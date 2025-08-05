library(ggplot2)
library(tidyverse)
library(patchwork)

sampleTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\sampleTableCount.txt", header=TRUE, sep = '\t')
libSizeTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\data\\libSize.txt", header=TRUE, sep = '\t')#For human elp1 injected mouse reference
conditionTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\conditions.txt", header=TRUE, sep = '\t')
colnames(libSizeTable) <- c("SampleNames","librarySize")
sampleTable$librarySize <- right_join(sampleTable, libSizeTable)$librarySize
sampleTable <- sampleTable[,c("SampleNames","FileNames","Samples","celltype","ID","replicate","zygosity","size","mutation_type","CRISPR_round","librarySize")]
sampleTable$genotype <- paste0(sampleTable$zygosity,"_",sampleTable$size,"_",sampleTable$mutation_type)
sampleTable$Samples <- paste0(sampleTable$celltype,"_",sampleTable$ID,"_",sampleTable$replicate)
sampleTable <- sampleTable[1:48,]

df <- read.table("C:\\Users\\ricar\\Downloads\\CADM2\\sig_exon.txt",header=F)
colnames(df) <- c("SampleNames", "Exon","Counts")
df1 <- left_join(sampleTable, df)


df1$Exon <- gsub("exon1_long_full","Exon 1",df1$Exon)
df1$Exon <- gsub("RS1_hypothetical_exon_Duff","RE1",df1$Exon)
df1$Exon <- gsub("1st_alt_exon2_and_RS2_exon","RE2",df1$Exon)
df1$Exon <- gsub("1st_alt_exon1","Alt Exon 1",df1$Exon)
df1$Exon <- gsub("exon2","Exon 2",df1$Exon)
df1$Exon <- gsub("exon_3","Exon 3",df1$Exon)
df1$Exon <- gsub("exon_4","Exon 4",df1$Exon)
df1$Exon <- gsub("exon_5","Exon 5",df1$Exon)
df1$Exon <- gsub("exon_6","Exon 6",df1$Exon)
df1$Exon <- gsub("exon_7","Exon 7",df1$Exon)
df1$Exon <- gsub("exon_8","Exon 8",df1$Exon)
df1$Exon <- gsub("exon_9","Exon 9",df1$Exon)
df1$Exon <- gsub("exon_10","Exon 10",df1$Exon)
df1$Exon <- gsub("exon_11_3UTR_long_only","Exon 11 3'UTR",df1$Exon)
df1$Exon <- gsub("exon_11_long_full","Exon 11",df1$Exon)


  
  
  
  df1$Normalized_Counts <- (df1$Counts/df1$librarySize)*1e6
  df2 <- df1 %>%
         filter(Exon %in% c("Exon 1","RE1","RE2", "Alt Exon 1", "Exon 2","Exon 3","Exon 4", "Exon 5", "Exon 6", "Exon 7", "Exon 8", "Exon 9","Exon 10", "Exon 11 3'UTR", "Exon 11"))
  #df2$Exon = factor(df2$Exon, levels =c("exon1_long_full","RS1_hypothetical_exon_Sibley","1st_alt_exon2_and_RS2_exon","exon_10"))
  #df2$genotype = factor(df2$genotype, levels =c("wt_wt_wt","het_500bp_del","het_100kb_del","het_500bp_inv","hom_500bp_del"))



df2$genotype <-gsub("wt_wt_wt","WT",df2$genotype)
df2$genotype <-gsub("het_500bp_del","HET 500bp",df2$genotype)
df2$genotype <-gsub("het_100kb_del","HET 100kb",df2$genotype)
df2$genotype <-gsub("het_500bp_inv","HET INV",df2$genotype)
df2$genotype <-gsub("hom_500bp_del","HOM DEL",df2$genotype)

df2$Exon = factor(df2$Exon, levels=c("Exon 1","RE1","RE2", "Alt Exon 1","Exon 2","Exon 3","Exon 4", "Exon 5", "Exon 6", "Exon 7", "Exon 8","Exon 9","Exon 10", "Exon 11 3'UTR", "Exon 11"))
df2$genotype = factor(df2$genotype, levels=c("WT","HET 500bp","HET 100kb","HET INV","HOM DEL"))

df2$zygosity <- gsub("het","HET",df2$zygosity)
df2$zygosity <- gsub("hom","HOM",df2$zygosity)
df2$zygosity <- gsub("wt","WT",df2$zygosity)

df2$zygosity <- factor(df2$zygosity, levels=c("WT",
                                              "HET",
                                              "HOM"))

iN_df <- df2 %>%
          filter(celltype=="iN") %>%
          filter(ID!="r1B6") %>% filter(SampleNames!="iN_rC5_rep2_500bp_inv.count")

pdf("iN_exons_all.pdf",height=4,width=12)
ggerrorplot(iN_df, x = "genotype", y = "Normalized_Counts", 
            desc_stat = "mean_sd",
            color = "genotype",
            size = 0.5,
            legend="top",
            ylab="Normalized Counts Per Million",
            xlab="") + font("ylab",size=14,face="bold") + font("xy.text",size = 14, face = "bold") +
  rotate_x_text() + facet_wrap(~Exon, scales="free_y")
dev.off()

iPSC_df <- df2 %>%
  filter(celltype=="iPSC") %>%
  filter(ID!="r1B6")
pdf("iPSC_exons_all.pdf",height=3,width=15)
ggerrorplot(iPSC_df, x = "genotype", y = "Normalized_Counts", 
            desc_stat = "mean_sd",
            color = "genotype",
            size = 1.4,
            legend="center",
            ylab="Normalized Counts Per Million",
            xlab="") +
  rotate_x_text() + facet_wrap(~Exon, scales="free_y")
dev.off()


pdf("iN_exons_all_exons_free_y.pdf",height=8,width=15)
 ggplot(iN_df, aes(x=genotype, y=Normalized_Counts, fill=genotype)) +
  geom_boxplot(width=0.7, position = 'dodge', outliers = F) +
  ylab("Normalized Counts Per Million") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   facet_wrap(~ Exon,scales="free_y")
dev.off()


pdf("iPSC_exons_all_free_y.pdf",width=15)
ggplot(iPSC_df, aes(x=genotype, y=Normalized_Counts, fill=genotype)) +
  geom_boxplot(width=0.7, position = 'dodge', outliers = F) +
  ylab("Normalized Counts Per Million") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~ Exon)
dev.off()


 for (exon in unique(iPSC_df$Exon)){
   
   tmp <- iPSC_df %>% filter(Exon==exon)
   wt <- tmp %>% filter(genotype=="WT")
   het_inv <- tmp %>% filter(genotype=="HET INV")
   het_del_100kb <- tmp %>% filter(genotype=="HET 100kb")
   het_del_500bp <- tmp %>% filter(genotype=="HET 500bp")
   hom_del <- tmp %>% filter(genotype=="HOM DEL")
   
  wt_zyg <- tmp %>% filter(zygosity=="WT")
  hom <- tmp %>% filter(zygosity=="HOM")
  het <- tmp %>% filter(zygosity=="HET")
   
   #pdf(paste0("iPSC_",exon,".pdf"),width=15)
   #print(ggplot(tmp, aes(x=genotype, y=Normalized_Counts, fill=genotype)) +
  #   geom_boxplot(width=0.7, position = 'dodge', outliers = F) +
  #   ylab("Normalized Counts Per Million") +
  #   theme_classic() +
  #   geom_jitter(position=position_jitterdodge(0.5)) +
  #   theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  # dev.off()

  print(paste(exon, t.test(wt$Normalized_Counts,het_inv$Normalized_Counts)$p.value,
             t.test(wt$Normalized_Counts,het_del_500bp$Normalized_Counts)$p.value,
                t.test(wt$Normalized_Counts,het_del_100kb$Normalized_Counts)$p.value,
                t.test(wt$Normalized_Counts,hom_del$Normalized_Counts)$p.value,t.test(wt_zyg$Normalized_Counts,het$Normalized_Counts)$p.value))
  #print(paste(exon,mean(het_inv$Normalized_Counts)/mean(wt$Normalized_Counts),
  #           mean(het_del_500bp$Normalized_Counts)/mean(wt$Normalized_Counts),
  #             mean(het_del_100kb$Normalized_Counts)/mean(wt$Normalized_Counts),
  #             mean(hom_del$Normalized_Counts)/mean(wt$Normalized_Counts)))

  #print(paste(exon,t.test(wt_zyg$Normalized_Counts,het$Normalized_Counts)$p.value,t.test(wt_zyg$Normalized_Counts,hom$Normalized_Counts)$p.value))
  #print(paste(exon,mean(het$Normalized_Counts)/mean(wt_zyg$Normalized_Counts),
  #            mean(hom$Normalized_Counts)/mean(wt_zyg$Normalized_Counts)))
  #print("----------")
 }


########Make Junction plots
library(ggplot2)
library(tidyverse)
library(patchwork)

sampleTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\sampleTableCount.txt", header=TRUE, sep = '\t')
libSizeTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\data\\libSize.txt", header=TRUE, sep = '\t')#For human elp1 injected mouse reference
conditionTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\conditions.txt", header=TRUE, sep = '\t')
colnames(libSizeTable) <- c("SampleNames","librarySize")
sampleTable$librarySize <- right_join(sampleTable, libSizeTable)$librarySize
sampleTable <- sampleTable[,c("SampleNames","FileNames","Samples","celltype","ID","replicate","zygosity","size","mutation_type","CRISPR_round","librarySize")]
sampleTable$genotype <- paste0(sampleTable$zygosity,"_",sampleTable$size,"_",sampleTable$mutation_type)
sampleTable$Samples <- paste0(sampleTable$celltype,"_",sampleTable$ID,"_",sampleTable$replicate)
sampleTable <- sampleTable[1:48,]

df <- read.table("C:\\Users\\ricar\\Downloads\\CADM2\\final_junction_RNAseq.txt",header=T)
colnames(df) <- c("SampleNames","unique_map","junction", "Genotype")
df$SampleNames<-paste0(df$SampleNames, ".count")
df$Genotype <-gsub("HET_500bp","HET 500bp DEL",df$Genotype)
df$Genotype <-gsub("HET_100kb","HET 100kb DEL",df$Genotype)
df$Genotype <-gsub("HOM_DEL","HOM 500bp DEL",df$Genotype)
df$Genotype <-gsub("INV","HET 500bp INV",df$Genotype)

df <- left_join(sampleTable, df)
df$Genotype = factor(df$Genotype, levels=c("WT","HET 500bp DEL","HET 100kb DEL","HET 500bp INV","HOM 500bp DEL"))

df$zygosity <- df$Genotype
df$zygosity <- gsub("HET 500bp INV", "HET",df$zygosity)
df$zygosity <- gsub("HET 100kb DEL", "HET",df$zygosity)
df$zygosity <- gsub("HET 500bp DEL", "HET",df$zygosity)
df$zygosity <- gsub("HOM 500bp DEL", "HOM",df$zygosity)
df$zygosity = factor(df$zygosity, levels=c("WT","HET","HOM"))
df$Normalized_Counts <- (df$unique_map/df$librarySize)*1e6

iN_df <- df %>%
  filter(celltype=="iN") %>%
  filter(ID!="r1B6") %>% filter()

iPSC_df <- df %>%
  filter(celltype=="iPSC") %>%
  filter(ID!="r1B6")

pdf("iN_junction.pdf", width = 3, height = 6)
ggerrorplot(iN_df, x = "Genotype", y = "Normalized_Counts", 
            desc_stat = "mean_sd",
            color = "zygosity",
            size = 1.4,
            legend="center",
            ylab="Normalized Counts Per Million",
            xlab="") +
  rotate_x_text() + facet_wrap(~junction)
dev.off()

pdf("iN_junction.pdf",height=8,width=15)
ggplot(df, aes(x=Genotype, y=Normalized_Counts, fill=Genotype)) +
  geom_boxplot(width=0.7, position = 'dodge', outliers = F) +
  ylab("Normalized Counts Per Million") +
  theme_classic() +
  geom_jitter(position=position_jitterdodge(0.5)) +
  theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~ junction, scales="free_y")
dev.off()


for (junctions in unique(iN_df$junction)){
  
  tmp <- iN_df %>% filter(junction==junctions)
  wt <- tmp %>% filter(genotype=="wt_wt_wt")
  het_inv <- tmp %>% filter(genotype=="het_500bp_inv")
  het_del_100kb <- tmp %>% filter(genotype=="het_100kb_del")
  het_del_500bp <- tmp %>% filter(genotype=="het_500bp_del")
  hom_del <- tmp %>% filter(genotype=="hom_500bp_del")
  
  #pdf(paste0("iN_",junctions,"_narrow.pdf"), width = 3)
  #print(ggerrorplot(tmp, x = "Genotype", y = "Normalized_Counts", 
  #                  desc_stat = "mean_sd",
  #                  color = "zygosity",
  #                  size = 1.4,
  #                  legend="center",
  #                  ylab="Normalized Counts Per Million",
  #                  xlab="") +
  #        rotate_x_text())
  #print(ggplot(tmp, aes(x=Genotype, y=Normalized_Counts, fill=Genotype)) +
  #        geom_boxplot(width=0.7, position = 'dodge', outliers = F) +
  #        ylab("Normalized Counts Per Million") +
  #        theme_classic() +
  #        geom_jitter(position=position_jitterdodge(0.5)) +
  #        theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  #dev.off()
  
  print(paste(junctions, t.test(wt$Normalized_Counts,het_inv$Normalized_Counts)$p.value,
              t.test(wt$Normalized_Counts,het_del_500bp$Normalized_Counts)$p.value,
                t.test(wt$Normalized_Counts,het_del_100kb$Normalized_Counts)$p.value,
                t.test(wt$Normalized_Counts,hom_del$Normalized_Counts)$p.value))
  print(paste(junctions,mean(het_inv$Normalized_Counts)/mean(wt$Normalized_Counts),
               mean(het_del_500bp$Normalized_Counts)/mean(wt$Normalized_Counts),
               mean(het_del_100kb$Normalized_Counts)/mean(wt$Normalized_Counts),
               mean(hom_del$Normalized_Counts)/mean(wt$Normalized_Counts)))
}


###########Exon plot from DEXSeq
a <- dxd@assays@data$mu[,1:20]
colnames(a) <- sampleTable$SampleNames

library(reshape2)
b <-melt(a)

colnames(b) <- c("Exon", "SampleNames","mu")

c <- left_join(sampleTable,b)

c$Normalized_Counts1 <- (c$mu/c$librarySize)/1e6 

c$Exon <- gsub("ENSG00000175161:E001","Exon 1",c$Exon)
c$Exon <- gsub("ENSG00000175161:E002","RE1",c$Exon)
c$Exon <- gsub("ENSG00000175161:E004","RE2",c$Exon)
c$Exon <- gsub("ENSG00000175161:E003","Alt Exon 1",c$Exon)
c$Exon <- gsub("ENSG00000175161:E005","Exon 2",c$Exon)
c$Exon <- gsub("ENSG00000175161:E006","Exon 3",c$Exon)
c$Exon <- gsub("ENSG00000175161:E007","Exon 4",c$Exon)
c$Exon <- gsub("ENSG00000175161:E008","Exon 5",c$Exon)
c$Exon <- gsub("ENSG00000175161:E009","Exon 6",c$Exon)
c$Exon <- gsub("ENSG00000175161:E010","Exon 7",c$Exon)
c$Exon <- gsub("ENSG00000175161:E011","Exon 8",c$Exon)
c$Exon <- gsub("ENSG00000175161:E012","Exon 9",c$Exon)
c$Exon <- gsub("ENSG00000175161:E013","Exon 10",c$Exon)
c$Exon <- gsub("ENSG00000175161:E014","Exon 11",c$Exon)

c$Exon = factor(c$Exon, levels=c("Exon 1","RE1","Alt Exon 1","RE2","Exon 2","Exon 3","Exon 4", "Exon 5", "Exon 6", "Exon 7", "Exon 8","Exon 9","Exon 10", "Exon 11"))
c$Genotype = factor(c$Genotype, levels=c("WT","HET 500bp DEL","HET 100kb DEL","HET 500bp INV","HOM 500bp DEL"))
c$zygosity <- c$Genotype
c$zygosity <- gsub("HET 500bp INV", "HET",c$zygosity)
c$zygosity <- gsub("HET 100kb DEL", "HET",c$zygosity)
c$zygosity <- gsub("HET 500bp DEL", "HET",c$zygosity)
c$zygosity <- gsub("HOM 500bp DEL", "HOM",c$zygosity)
c$zygosity = factor(c$zygosity, levels=c("WT","HET","HOM"))

pdf("exon_plot_mu_raw.pdf",width=12)
ggerrorplot(c, x = "Exon", y = "log10(mu)", 
            desc_stat = "mean_sd",
            color = "Genotype",
            size = 1.4,
            legend="center",
            ylab="log10(Corrected Counts)",
            xlab="") +
  rotate_x_text() + scale_y_break(c(-0.5,0.5)) + scale_y_break(c(1.5,1.75)) +  scale_y_break(c(2.25,2.3))
# +  facet_wrap(~ Genotype)
dev.off()

for (exons in unique(c$Exon)){
  
  tmp <- c %>% filter(Exon==exons)
  wt <- tmp %>% filter(genotype=="wt_wt_wt")
  het_inv <- tmp %>% filter(genotype=="het_500bp_inv")
  het_del_100kb <- tmp %>% filter(genotype=="het_100kb_del")
  het_del_500bp <- tmp %>% filter(genotype=="het_500bp_del")
  hom_del <- tmp %>% filter(genotype=="hom_500bp_del")
  
  print(paste(exons, t.test(wt$Normalized_Counts,het_inv$Normalized_Counts)$p.value,
              t.test(wt$Normalized_Counts,het_del_500bp$Normalized_Counts)$p.value,
              t.test(wt$Normalized_Counts,het_del_100kb$Normalized_Counts)$p.value,
              t.test(wt$Normalized_Counts,hom_del$Normalized_Counts)$p.value))
  print(paste(exons,mean(het_inv$Normalized_Counts)/mean(wt$Normalized_Counts),
              mean(het_del_500bp$Normalized_Counts)/mean(wt$Normalized_Counts),
              mean(het_del_100kb$Normalized_Counts)/mean(wt$Normalized_Counts),
              mean(hom_del$Normalized_Counts)/mean(wt$Normalized_Counts)))
  print("--------")
}

