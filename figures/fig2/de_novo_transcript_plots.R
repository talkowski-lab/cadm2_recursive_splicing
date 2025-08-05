library(rtracklayer)
library(ggplot2)
library(tidyverse)
library(ggtranscript)

bed_file <- read.table("C:\\Users\\ricar\\Downloads\\CADM2\\test.bed", header = FALSE,sep='\t')
colnames(bed_file) <-c("seqnames","start","end","strand","gene_name","transcript_name")
bed_file$transcript_name = factor(bed_file$transcript_name, levels=c("Novel 3","Novel 2","Novel 1","CADM2 Short","Alt Exon 1+RE2+CADM2 Short","CADM2+RE2","CADM2+RE1+RE2"))
bed <- bed_file %>% dplyr::as_tibble()


 tplot <- bed %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = "exon")
  ) +
  geom_intron(
    data = to_intron(bed, "transcript_name"),
    aes(strand = strand),
    arrow.min.intron.length = 500,
  ) +
  ylab("Transcript Name") +
  theme_classic() +
  theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(vjust = 0.5, hjust=1)) +
   scale_y_discrete(expand = c(0,0))
  
 pdf("transcripts.pdf", width=11, height = 3)
 tplot
 dev.off()
 

 bed_file <- read.table("C:\\Users\\ricar\\Downloads\\CADM2\\novel3.bed", header = FALSE,sep='\t')
 colnames(bed_file) <-c("seqnames","start","end","strand","gene_name","transcript_name")
 bed <- bed_file %>% dplyr::as_tibble()
 
 novel3 <- bed %>%
   ggplot(aes(
     xstart = start,
     xend = end,
     y = transcript_name
   )) +
   geom_range(
     aes(fill = "exon",height=0.25)
   ) +
   geom_intron(
     data = to_intron(bed, "transcript_name"),
     aes(strand = strand),
     arrow.min.intron.length = 500,
   ) +
   ylab("Transcript Name") +
   theme_classic() +
   theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(vjust = 0.5, hjust=1))
 
 pdf("novel3_zoom.pdf", width=11, height = 3)
 novel3
 dev.off()
 
 ####Canonical CADM2
 
############From gtf for canonical camd2
 gtf_path <- "C:\\Users\\ricar\\Downloads\\Homo_sapiens.GRCh38.105.chr.gtf.gz"
 gtf <- rtracklayer::import(gtf_path)
 
 class(gtf)
 
 gtf <- gtf %>% dplyr::as_tibble()
  
 gene_of_interest <- "CADM2"
 
 cadm2_annotation_from_gtf <- gtf %>% 
   dplyr::filter(
     !is.na(gene_name), 
     gene_name == gene_of_interest
   ) 
 
 # extract the required annotation columns
 cadm2_annotation_from_gtf <- cadm2_annotation_from_gtf %>% 
   dplyr::select(
     seqnames,
     start,
     end,
     strand,
     type,
     gene_name,
     transcript_name,
     transcript_biotype
   )
 
 cadm2_annotation_from_gtf %>% head() 
 
 # extract exons
 cadm2_exons <- cadm2_annotation_from_gtf %>% dplyr::filter(type == "exon")
 cadm2_exons <- cadm2_exons %>% select(seqnames,start,end,strand,gene_name,transcript_name)
 
 cadm2_exons$transcript_name  <- gsub("CADM2-201","CADM2-ENST00000383699",cadm2_exons$transcript_name)
 cadm2_exons$transcript_name  <- gsub("CADM2-203","CADM2-ENST00000407528",cadm2_exons$transcript_name)
 cadm2_exons$transcript_name  <- gsub("CADM2-202","CADM2-ENST00000405615",cadm2_exons$transcript_name)
 cadm2_exons$transcript_name  <- gsub("CADM2-204","CADM2-ENST00000473523",cadm2_exons$transcript_name)
 cadm2_exons$transcript_name  <- gsub("CADM2-205","CADM2-ENST00000485126",cadm2_exons$transcript_name)

 cadm2_transcripts <- rbind(bed, cadm2_exons)
 cadm2_transcripts$transcript_name = factor(cadm2_transcripts$transcript_name, levels=c("Novel 3","Novel 2","Novel 1","CADM2 Short","Alt Exon 1+RE2+CADM2 Short","CADM2+RE2","CADM2+RE1+RE2",
                                                                               "CADM2-ENST00000473523",
                                                                               "CADM2-ENST00000485126","CADM2-ENST00000405615",
                                                                               "CADM2-ENST00000407528","CADM2-ENST00000383699"))
 
 tplot <- cadm2_transcripts %>%
   ggplot(aes(
     xstart = start,
     xend = end,
     y = transcript_name
   )) +
   geom_range(
     aes(fill = "exon")
   ) +
   geom_intron(
     data = to_intron(cadm2_transcripts, "transcript_name"),
     aes(strand = strand),
     arrow.min.intron.length = 500,
   ) +
   ylab("Transcript Name") +
   theme_classic() +
   theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(vjust = 0.5, hjust=1)) +
   scale_y_discrete(expand = c(0,0))
 
 pdf("transcripts.pdf", width=11, height = 3)
 tplot
 dev.off()
 
 
 
tmp <- cadm2_transcripts %>% filter(transcript_name %in% c("CADM2-ENST00000485126","CADM2-ENST00000473523"))
 
cadm2_small <- tmp %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = "exon")
  ) +
  geom_intron(
    data = to_intron(tmp, "transcript_name"),
    aes(strand = strand),
    arrow.min.intron.length = 500,
  ) +
  ylab("Transcript Name") +
  theme_classic() +
  theme(text = element_text(size = 19),legend.position="none",axis.title.x=element_blank(),axis.text.x = element_text(vjust = 0.5, hjust=1)) +
  scale_y_discrete(expand = c(0,0))

pdf("cadm2_small_zoom.pdf", width=11, height = 3)
cadm2_small
dev.off()
 