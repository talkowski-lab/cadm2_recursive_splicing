library(tidyverse)
library(ggplot2)
library(patchwork)


extract_slope <- function(table_name, librarySize){
  print(table_name)
  dat <- read.table(table_name, sep="\t", stringsAsFactors = FALSE)
  colnames(dat) <- c("Sequence","Position","Reference Base","Read_Count","Read Results","Quality")
  dat$NormalizedRead_Count <- (dat$Read_Count/librarySize)*1e6
  print(summary(dat$NormalizedRead_Count))
  
  first_segment <- dat[dat$Position > 84959669 & dat$Position < 85238956,]
  middle_segment <- dat[dat$Position > 85238957 & dat$Position < 85511857,]
  end_segment <- dat[dat$Position > 85511858 & dat$Position < 85726144,]
  
  first_segment$Position <-(first_segment$Position - 84959669)/1E6
  middle_segment$Position <- (middle_segment$Position -85238957)/1E6
  end_segment$Position <- (end_segment$Position  - 85511858)/1E6


  #Normalized Count plots for each segment
  p_first <- ggplot(first_segment, aes(x=Position, y=log10(NormalizedRead_Count))) +
    geom_bar(stat="identity") +
    geom_smooth(method=lm , color="blue", fill="red", se=TRUE) +
    ylab("Normalized Read Count") +
    xlab("Relative Position") +
    ylim(0, 3.5e-6) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size=5))
  
  p_middle <- ggplot(middle_segment, aes(x=Position, y=log10(NormalizedRead_Count))) +
    geom_bar(stat="identity") +
    geom_smooth(method=lm , color="blue", fill="red", se=TRUE) +
    ylab("Normalized Read Count") +
    xlab("Relative Position") +
    ylim(0, 3.5e-6) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size=5))
  
  p_end <- ggplot(end_segment, aes(x=Position, y=log10(NormalizedRead_Count))) +
    geom_bar(stat="identity") +
    geom_smooth(method=lm , color="blue", fill="red", se=TRUE) +
    ylab("Normalized Read Count") +
    xlab("Relative Position") +
    ylim(0, 3.5e-6) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size=5))

  pdf(paste0("/data/talkowski/Samples/cadm2/data/output/recursive_plots/",table_name,"_normalized_counts_log.pdf"), width = 3, height=2)
  #pdf(paste0("C:/Documents and Settings/ricar/Downloads/CADM2/saw_tooth_plots/",table_name,"_normalized_counts.pdf"), width = 10)
  print(p_first + p_middle + p_end)
  dev.off()
  
  ## Count plots for each segment
  #p_first <- ggplot(first_segment, aes(x=Position, y=Read_Count)) +
  #  geom_bar(stat="identity") +
  #  geom_smooth(method=lm , color="blue", fill="red", se=TRUE) +
  #  ylab("Read Count") +
  #  xlab("Relative Position") +
  #  ylim(0, 80) +
  #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  #p_middle <- ggplot(middle_segment, aes(x=Position, y=Read_Count)) +
  #  geom_bar(stat="identity") +
  #  geom_smooth(method=lm , color="blue", fill="red", se=TRUE) +
  #  ylab("Read Count") +
  #  xlab("Relative Position") +
  #  ylim(0, 80) +
  #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  #p_end <- ggplot(end_segment, aes(x=Position, y=Read_Count)) +
  #  geom_bar(stat="identity") +
  #  geom_smooth(method=lm , color="blue", fill="red", se=TRUE) +
  #  ylab("Read Count") +
  #  xlab("Relative Position") +
  #  ylim(0, 80) +
  #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  #pdf(paste0("/data/talkowski/Samples/cadm2/data/output/recursive_plots/",table_name,"_counts.pdf"), width = 3, height=2)
  #pdf(paste0("C:/Documents and Settings/ricar/Downloads/CADM2/saw_tooth_plots/",table_name,"_counts.pdf"), width = 10)
  #print(p_first + p_middle + p_end)
  #dev.off()
  
  first_segment_cor <- coef(lm(first_segment$NormalizedRead_Count ~ first_segment$Position))[2]
  middle_segment_cor <- coef(lm(middle_segment$NormalizedRead_Count ~  middle_segment$Position))[2]
  end_segment_cor <- coef(lm(end_segment$NormalizedRead_Count ~ end_segment$Position))[2]
  return(paste(table_name,first_segment_cor,   middle_segment_cor,end_segment_cor,max(dat$Read_Count),min(dat$Read_Count)))
}


df <-NULL
setwd("/data/talkowski/Samples/cadm2/data/IGV_RNAseq_EditedGTF")
directory <- getwd();
libSize <- read.table("/data/talkowski/Samples/cadm2/data/libSize.txt",sep="\t",header = T)
  #libSize <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\data\\libSize.txt",sep="\t",header = T)
libSize$Sample<-gsub(".count",".mpileup",libSize$Sample)

for(file in 1:dim(libSize)[1]){
  print(file)
  df<-rbind(df, extract_slope(libSize[file,1], libSize[file,2]))
}


df <- as.data.frame(df) %>% separate(V1,c("Sample","first","middle","end", "max","min"), sep = "\\s")
write.table(df,"slopes_2.txt", sep="\t", quote=FALSE, row.names = FALSE)
#############
df <- read.table("C:\\Users\\ricar\\Downloads\\CADM2\\slopes_2.txt",header = T,sep = '\t')

df$genotype <-gsub("HET_500bp","HET 500bp DEL",df$genotype)
df$genotype <-gsub("HET_INV","HET 500bp INV",df$genotype)
df$genotype <-gsub("HOM","HOM 500bp DEL",df$genotype)
df$genotype = factor(df$genotype, levels=c("WT","HET 500bp DEL","HET 500bp INV","HOM 500bp DEL"))

df$zygosity <- df$genotype
df$zygosity <- gsub("HET 500bp INV", "HET",df$zygosity)
df$zygosity <- gsub("HET 500bp DEL", "HET",df$zygosity)
df$zygosity <- gsub("HOM 500bp DEL", "HOM",df$zygosity)
df$zygosity = factor(df$zygosity, levels=c("WT","HET","HOM"))

df<- df %>% filter(celltype=="iN") %>% filter(genotype!="HET_100kb")
df <- df[1:16,]

HOM <- df %>%filter(genotype=="HOM")
mean(HOM$first)

HET_500bp <- df %>%filter(genotype=="HET_500bp")
mean(HET_500bp$first)

HET_INV <- df %>%filter(genotype=="HET_INV")
mean(HET_INV$first)

ting_2 <- ting %>% filter(celltype=="iN")

ggplot(ting_2, aes(x=genotype, y=value, fill = variable)) +
  geom_boxplot() +
  #ylim(4e-5,0-7.5e-5) +
  ylab("Read Count") +
  theme_classic() +
  geom_beeswarm()
  #geom_jitter(position = position_jitter(height = .05, width = .05)) +
  #theme(legend.position = "none") +
  theme(text = element_text(size = 25),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(df, aes(x=genotype, y=end, fill=genotype)) +
  geom_violin(aes(alpha=genotype)) +
  scale_fill_grey() +
  #facet_wrap(vars(Tissue), nrow=1, ncol=5) +
  geom_boxplot(width=0.1, fill="white") +
  # ylim(-0.8,0) +
  ylab("Slope") +
  theme_classic() +
  geom_jitter() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_violin_first <- ggplot(df, aes(x=genotype, y=first, fill=genotype)) +
  geom_violin(aes(alpha=genotype)) +
  scale_fill_grey() +
  #facet_wrap(vars(Tissue), nrow=1, ncol=5) +
  geom_boxplot(width=0.1, fill="white") +
  #ylim(-0.8,0) +
  ylab("Slope") +
  theme_classic() +
  geom_jitter() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_violin_middle <- ggplot(df, aes(x=genotype, y=middle, fill=genotype)) +
  geom_violin(aes(alpha=genotype)) +
  scale_fill_grey() +
  #facet_wrap(vars(Tissue), nrow=1, ncol=5) +
  geom_boxplot(width=0.1, fill="white") +
  #ylim(-0.8,0) +
  ylab("Slope") +
  theme_classic() +
  geom_jitter() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_violin_end <- ggplot(df, aes(x=genotype, y=end, fill=genotype)) +
  geom_violin(aes(alpha=genotype)) +
  scale_fill_grey() +
  #facet_wrap(vars(Tissue), nrow=1, ncol=5) +
  geom_boxplot(width=0.1, fill="white") +
 # ylim(-0.8,0) +
  ylab("Slope") +
  theme_classic() +
  geom_jitter() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 21),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("iN_slopes_cpm.pdf", width = 7,height = 7)
p_violin_first + p_violin_middle + p_violin_end
dev.off()

p_box_first <- ggplot(df, aes(x=genotype, y=first, fill = genotype)) +
geom_boxplot() +
  ylim(-2,1) +
  ylab("Slope") +
  theme_classic() +
  geom_jitter(width=0.01) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 25),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_box_middle <- ggplot(df, aes(x=genotype, y=middle, fill=genotype)) +
  geom_boxplot() +
  ylim(-2, 1) +
  ylab("Slope") +
  theme_classic() +
  geom_jitter(width=0.01) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 25),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_box_end <- ggplot(df, aes(x=genotype, y=end, fill=genotype)) +
  geom_boxplot() +
  ylim(-2,1) +
  ylab("Slope") +
  theme_classic() +
  geom_jitter(width=0.01) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 25),axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p_box_first + p_box_middle + p_box_end

pdf("iN_boxplot_slopes_cpm_no_pws.pdf",width =10,height = 4)
p_box_first + p_box_middle + p_box_end
dev.off()

####ggerrorplot
p1 <- ggerrorplot(df, x = "genotype", y = "first", 
            desc_stat = "mean_sd",
            color = "zygosity",
            size = 1.4,
            legend="center",
            ylab="Slope",
            xlab="") +
  ylim(-2,1) +
  rotate_x_text()

p2 <- ggerrorplot(df, x = "genotype", y = "middle", 
            desc_stat = "mean_sd",
            color = "zygosity",
            size = 1.4,
            legend="center",
            ylab="Slope",
            xlab="") +
  ylim(-2,1) +
  rotate_x_text()

p3 <- ggerrorplot(df, x = "genotype", y = "end", 
            desc_stat = "mean_sd",
            color = "zygosity",
            size = 1.4,
            legend="center",
            ylab="Slope",
            xlab="") +
  ylim(-2,1) +
  rotate_x_text()


pdf("iN_boxplot_slopes_cpm_no_pws.pdf",width =10,height = 4)
p1 + p2 + p3
dev.off()

####ggerrorplot
#######


WT_1 <- df %>% filter(genotype=="WT") %>% dplyr::select(first)
INV_1 <-df %>% filter(genotype=="HET_INV") %>% dplyr::select(first)
HOM_1 <- df %>% filter(genotype=="HOM") %>% dplyr::select(first)
DEL_100KB_1 <- df %>% filter(genotype=="HET_100kb") %>% dplyr::select(first)
DEL_500BP_1 <- df %>% filter(genotype=="HET_500bp") %>% dplyr::select(first)

t.test(WT_1,HOM_1)$p.value
t.test(WT_1,INV_1)$p.value
t.test(WT_1,DEL_100KB_1)$p.value
t.test(WT_1,DEL_500BP_1)$p.value



WT_2 <- df %>% filter(genotype=="WT") %>% dplyr::select(middle)
INV_2 <-df %>% filter(genotype=="HET_INV") %>% dplyr::select(middle)
HOM_2 <- df %>% filter(genotype=="HOM") %>% dplyr::select(middle)
DEL_100KB_2 <- df %>% filter(genotype=="HET_100kb") %>% dplyr::select(middle)
DEL_500BP_2 <- df %>% filter(genotype=="HET_500bp") %>% dplyr::select(middle)

t.test(WT_2,HOM_2)$p.value
t.test(WT_2,INV_2)$p.value
t.test(WT_2,DEL_100KB_2)$p.value
t.test(WT_2,DEL_500BP_2)$p.value


WT_3 <- df %>% filter(genotype=="WT") %>% dplyr::select(end)
INV_3 <-df %>% filter(genotype=="HET_INV") %>% dplyr::select(end)
HOM_3 <- df %>% filter(genotype=="HOM") %>% dplyr::select(end)
DEL_100KB_3 <- df %>% filter(genotype=="HET_100kb") %>% dplyr::select(end)
DEL_500BP_3 <- df %>% filter(genotype=="HET_500bp") %>% dplyr::select(end)

t.test(WT_3,HOM_3)$p.value
t.test(WT_3,INV_3)$p.value
t.test(WT_3,DEL_100KB_3)$p.value
t.test(WT_3,DEL_500BP_3)$p.value










ff <- function(number){
  pdf(paste0(number,".pdf"))
  p <- qplot(1:20, 1:20)
  q <- qplot(1:20, 1:20)
  print(p+q)
  dev.off()
  
  return(1)
  

}

for(file in 1:5){
  df <- rbind(ff(file))
}