setwd("C:\\Documents and Settings/ricar/Downloads/CADM2/RNA_seq_EditedGTF_concensus/")
sampleTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\sampleTableCount.txt", header=TRUE, sep = '\t')
#libSizeTable <- read.table("/data/talkowski/Samples/fd_mouse_tissue/training_MED/qc/bamqc/libsz.txt", header=TRUE, sep = '\t')
libSizeTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\data\\Alignment_CADM2_RNAseq_EditedGTF\\libSize.txt", header=TRUE, sep = '\t')#For human elp1 injected mouse reference
conditionTable <- read.table("D:\\data\\talkowski\\Samples\\cadm2\\tables\\conditions.txt", header=TRUE, sep = '\t')
colnames(libSizeTable) <- c("SampleNames","librarySize")
sampleTable$librarySize <- right_join(sampleTable, libSizeTable)$librarySize
sampleTable <- sampleTable[,c("SampleNames","FileNames","Samples","celltype","ID","replicate","zygosity","size","mutation_type","CRISPR_round","librarySize")]
sampleTable$Sample_ID <- paste0(sampleTable$celltype,"_",sampleTable$ID)
sampleTable$celltype  <- factor(sampleTable$celltype, levels = c("iN","iPSC"))
sampleTable$zygosity <- factor(sampleTable$zygosity,levels=c("wt","hom","het"))
sampleTable$mutation_type <- factor(sampleTable$mutation_type,levels=c("del","inv","wt"))
sampleTable$replicate <- factor(sampleTable$replicate,levels=c("rep1","rep2", "rep3", "rep4"))

#Modify this line for cap_seq and RNA-seq files because the filenames are different.
sampleTable$Files <- paste0("cpm","_",sub(".count",".concensus",sampleTable$FileNames))
sampleTable$genotype <- paste0(sampleTable$zygosity,"_",sampleTable$size,"_",sampleTable$mutation_type)


sampleTable_iPSC_hom_500bp_del <- sampleTable %>% 
  filter(celltype=="iPSC") %>%
  filter(genotype=="hom_500bp_del")
  
library(dplyr)
library(readr)
  
  df <- sampleTable$Files %>% 
    lapply(read.table) %>% 
    bind_cols()
  
iPSC_hom_500bp_del <- df %>%
   select(V1...1,V2...2,V3...3,V4...4,V4...8,V4...12,V4...16,V4...20,V4...24,V4...28,V4...32) %>%
  mutate_if(is.character, as.numeric) %>%
   rowwise() %>% mutate(m = mean(c(V4...4,V4...8,V4...12,V4...16,V4...20,V4...24,V4...28,V4...32)))
                        
                        
########################het_500bp_del
sampleTable_iPSC_het_500bp_del <- sampleTable %>% 
  filter(celltype=="iPSC") %>%
  filter(genotype=="het_500bp_del")

library(dplyr)
library(readr)

df <- sampleTable_iPSC_het_500bp_del$Files %>% 
  lapply(read.table) %>% 
  bind_cols()

iPSC_het_500bp_del <- df %>%
  select(V1...1,V2...2,V3...3,V4...4,V4...8,V4...12,V4...16) %>%
  mutate_if(is.character, as.numeric) %>%
  rowwise() %>% mutate(m = mean(c(V4...4,V4...8,V4...12,V4...16)))
######################het_500bp_inv

sampleTable_iPSC_het_500bp_inv <- sampleTable %>% 
  filter(celltype=="iPSC") %>%
  filter(genotype=="het_500bp_inv")

library(dplyr)
library(readr)

df <- sampleTable_iPSC_het_500bp_inv$Files %>% 
  lapply(read.table) %>% 
  bind_cols()

iPSC_het_500bp_inv <- df %>%
  select(V1...1,V2...2,V3...3,V4...4,V4...8,V4...12,V4...16) %>%
  mutate_if(is.character, as.numeric) %>%
  rowwise() %>% mutate(m = mean(c(V4...4,V4...8,V4...12,V4...16)))
######################wt_wt_wt

sampleTable_iPSC_wt_wt_wt <- sampleTable %>% 
  filter(celltype=="iPSC") %>%
  filter(genotype=="wt_wt_wt")

library(dplyr)
library(readr)

df <- sampleTable_iPSC_wt_wt_wt$Files %>% 
  lapply(read.table) %>% 
  bind_cols()

iPSC_wt_wt_wt <- df %>%
  select(V1...1,V2...2,V3...3,V4...4,V4...8,V4...12,V4...16) %>%
  mutate_if(is.character, as.numeric) %>%
  rowwise() %>% mutate(m = mean(c(V4...4,V4...8,V4...12,V4...16)))
######################het_100kb_del

sampleTable_iPSC_het_100kb_del <- sampleTable %>% 
  filter(celltype=="iPSC") %>%
  filter(genotype=="het_100kb_del")

library(dplyr)
library(readr)

df <- sampleTable_iPSC_het_100kb_del$Files %>% 
  lapply(read.table) %>% 
  bind_cols()

iPSC_het_100kb_del <- df %>%
  select(V1...1,V2...2,V3...3,V4...4,V4...8,V4...12,V4...16) %>%
  mutate_if(is.character, as.numeric) %>%
  rowwise() %>% mutate(m = mean(c(V4...4,V4...8,V4...12,V4...16)))

#############
write.table(iPSC_wt_wt_wt[,c(1,2,3,8)], "RNAseq_iPSC_wt_wt_wt.bedGraph",row.names=FALSE, quote=FALSE, col.names = FALSE,sep = "\t")
write.table(iPSC_hom_500bp_del[,c(1,2,3,12)], "RNAseq_iPSC_hom_500bp_del.bedGraph",row.names=FALSE, quote=FALSE, col.names = FALSE,sep = "\t")
write.table(iPSC_het_500bp_del[,c(1,2,3,8)], "RNAseq_iPSC_het_500bp_del.bedGraph",row.names=FALSE, quote=FALSE, col.names = FALSE,sep = "\t")
write.table(iPSC_het_100kb_del[,c(1,2,3,8)], "RNAseq_iPSC_het_100kb_del.bedGraph",row.names=FALSE, quote=FALSE, col.names = FALSE,sep = "\t")
write.table(iPSC_het_500bp_inv[,c(1,2,3,8)], "RNAseq_iPSC_het_500bp_inv.bedGraph",row.names=FALSE, quote=FALSE, col.names = FALSE,sep = "\t")


#############
library(data.table)
library(Gviz)
setwd("C:\\Documents and Settings/ricar/Downloads/CADM2/RNA_seq_EditedGTF_concensus/")


bedgraph_iPSC_wt_wt_wt <- fread(
  'RNAseq_iPSC_wt_wt_wt.bedGraph',
  col.names = c('chromosome', 'start', 'end', 'value')
)

bedgraph_iPSC_hom_500bp_del <- fread(
  'RNAseq_iPSC_hom_500bp_del.bedGraph',
  col.names = c('chromosome', 'start', 'end', 'value')
)


bedgraph_iPSC_het_500bp_del <- fread(
  'RNAseq_iPSC_het_500bp_del.bedGraph',
  col.names = c('chromosome', 'start', 'end', 'value')
)

bedgraph_iPSC_het_100kb_del <- fread(
  'RNAseq_iPSC_het_100kb_del.bedGraph',
  col.names = c('chromosome', 'start', 'end', 'value')
)

bedgraph_iPSC_het_500bp_inv <- fread(
  'RNAseq_iPSC_het_500bp_inv.bedGraph',
  col.names = c('chromosome', 'start', 'end', 'value')
)


bedgraph_iPSC_wt_wt_wt$value <- log10(bedgraph_iPSC_wt_wt_wt$value + 1)
bedgraph_iPSC_het_500bp_del$value <- log10(bedgraph_iPSC_het_500bp_del$value + 1)
bedgraph_iPSC_het_500bp_inv$value <- log10(bedgraph_iPSC_het_500bp_inv$value + 1)
bedgraph_iPSC_het_100kb_del$value <- log10(bedgraph_iPSC_het_100kb_del$value + 1)
bedgraph_iPSC_hom_500bp_del$value <- log10(bedgraph_iPSC_hom_500bp_del$value + 1)

summary(bedgraph_iPSC_wt_wt_wt$value)
summary(bedgraph_iPSC_hom_500bp_del$value)
summary(bedgraph_iPSC_het_500bp_del$value)
summary(bedgraph_iPSC_het_500bp_inv$value)
summary(bedgraph_iPSC_het_100kb_del$value)

# Specifiy the range to plot
thechr <- "3"
#st <- 84959546
#en <- 86066835
st_in1 <- 84953953
en_in1 <- 85732440

bedgraph_dt_one_chr <- bedgraph_iPSC_wt_wt_wt[chromosome == thechr]
dtrack_iPSC_wt_wt_wt <- DataTrack(
  range = bedgraph_dt_one_chr,
  type = "a",
  genome = 'hg38',
  name = "iPSC_wt",
  ylim=c(0,3),
  #smooth
)


bedgraph_dt_one_chr <- bedgraph_iPSC_hom_500bp_del[chromosome == thechr]
dtrack_iPSC_hom_500bp_del <- DataTrack(
  range = bedgraph_dt_one_chr,
  type = "a",
  genome = 'hg38',
  name = "hom_del",
  ylim=c(0,3),
  #smooth
)


bedgraph_dt_one_chr <- bedgraph_iPSC_het_500bp_del[chromosome == thechr]
dtrack_iPSC_het_500bp_del <- DataTrack(
  range = bedgraph_dt_one_chr,
  type = "a",
  genome = 'hg38',
  name = "het_500bp_del",
  ylim=c(0,3),
  #smooth
)


bedgraph_dt_one_chr <- bedgraph_iPSC_het_100kb_del[chromosome == thechr]
dtrack_iPSC_het_100kb_del <- DataTrack(
  range = bedgraph_dt_one_chr,
  type = "a",
  genome = 'hg38',
  name = "het_100kb_del",
  ylim=c(0,3),
  #smooth
)

bedgraph_dt_one_chr <- bedgraph_iPSC_het_500bp_inv[chromosome == thechr]
dtrack_iPSC_het_500bp_inv <- DataTrack(
  range = bedgraph_dt_one_chr,
  type = "a",
  genome = 'hg38',
  name = "het_inv",
  ylim=c(0,3),
  #smooth
)

####Gene structure

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

ucscGenes <- UcscTrack(genome="hg38", table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome="chr3", rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)

mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z
plotTracks(list(ucscGenes, ucscGenes2), chromosome = "chr3", from = st_in1, to = en_in1, transcriptAnnotation = "symbol")


###
library(rtracklayer)
RS1 <- import.bed('C:\\Users\\ricar\\Downloads\\CADM2\\CADM2_ExonDef_RS1.bed')

RS1track <- GeneRegionTrack(
  RS1,
  genome = "hg38",
  chromosome = thechr, start = st_in1, end = en_in1,
  showId = TRUE,
  name = "RS1",
  group = c("RS1-Sibley",
            "RS1-Duff",
            "Alt-Exon1")
)

RS1_annoT <- AnnotationTrack(RS1,
                         group = c("RS1-Sibley",
                                   "RS1-Duff",
                                   "Alt-Exon1"))



CRISPR_Edits <- import.bed('C:\\Users\\ricar\\Downloads\\CADM2\\CADM2_CRISPR_del.bed')

CRISPRtrack <- GeneRegionTrack(
  CRISPR_Edits ,
  genome = "hg38",
  chromosome = thechr, start = st_in1, end = en_in1,
  showId = TRUE,
  name = "RS1",
  group = c("500bp_outer_del",
            "100kb_outer_del",
            "100kb_inner_del")
)

CRISPR_annoT <- AnnotationTrack(CRISPR_Edits,
                         group = c("500bp_outer_del",
                                   "100kb_outer_del",
                                   "100kb_inner_del"))



geneTrack <- BiomartGeneRegionTrack(genome="hg38",
                                    chromosome=thechr, start=st_in1, end=en_in1,
                                    transcriptAnnotation="symbol", name="Genes")
####

itrack <- IdeogramTrack(
  genome = "hg38", chromosome = thechr
)

grtrack <- GeneRegionTrack(
  txdb,
  chromosome = thechr, start = st_in1, end = en_in1,
  showId = TRUE,
  name = "Gene Annotation"
)

  pdf("RNA_iPSC_intron1_avg_gviz.pdf", height = 12)
    plotTracks(
        list(itrack,ucscGenes2,RS1_annoT,CRISPR_annoT,dtrack_iPSC_wt_wt_wt,dtrack_iPSC_het_500bp_del,dtrack_iPSC_het_500bp_inv,dtrack_iPSC_het_100kb_del,dtrack_iPSC_hom_500bp_del),
      from = st_in1, to = en_in1,
      transcriptAnnotation = "symbol",groupAnnotation="group", type=c('histogram','l')
    )
  dev.off()