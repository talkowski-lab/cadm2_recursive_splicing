library(data.table)
library(Gviz)
setwd("C:\\Documents and Settings\\ricar\\Downloads/")


bedgraph_iPSC_wt <- fread(
  'iPSC_wt.bedGraph',
  col.names = c('chromosome', 'start', 'end', 'value')
)

bedgraph_iPSC_del <- fread(
  'iPSC_hom_del.bedGraph',
  col.names = c('chromosome', 'start', 'end', 'value')
)


bedgraph_iN_wt <- fread(
  'iN_wt.bedGraph',
  col.names = c('chromosome', 'start', 'end', 'value')
)


bedgraph_iN_del <- fread(
  'iN_hom_del.bedGraph',
  col.names = c('chromosome', 'start', 'end', 'value')
)


# Specifiy the range to plot
thechr <- "3"
st <- 84959546
en <- 86066835


bedgraph_dt_one_chr <- bedgraph_iPSC_wt[chromosome == thechr]
dtrack_iPSC_WT <- DataTrack(
  range = bedgraph_dt_one_chr,
  type = "a",
  genome = 'hg38',
  name = "Seq. Depth"
)


bedgraph_dt_one_chr <- bedgraph_iPSC_del [chromosome == thechr]
dtrack_iPSC_del <- DataTrack(
  range = bedgraph_dt_one_chr,
  type = "a",
  genome = 'hg38',
  name = "Seq. Depth"
)


bedgraph_dt_one_chr <- bedgraph_iN_wt[chromosome == thechr]
dtrack_iN_WT <- DataTrack(
  range = bedgraph_dt_one_chr,
  type = "a",
  genome = 'hg38',
  name = "Seq. Depth"
)


bedgraph_dt_one_chr <- bedgraph_iN_del[chromosome == thechr]

dtrack_iN_del <- DataTrack(
  range = bedgraph_dt_one_chr,
  type = "a",
  genome = 'hg38',
  name = "Seq. Depth"
)


plotTracks(
  list(dtrack_iN_del),
  from = st, to = en
)

itrack <- IdeogramTrack(
  genome = "hg38", chromosome = thechr
)
gtrack <- GenomeAxisTrack()

plotTracks(
  list(itrack, gtrack, dtrack),
  from = st, to = en
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
plotTracks(list(ucscGenes, ucscGenes2), chromosome = "chr3", from = st, to = en, transcriptAnnotation = "symbol")

  
grtrack <- GeneRegionTrack(
 txdb,
    chromosome = thechr, start = st, end = en,
    showId = TRUE,
    name = "Gene Annotation"
  )
 
pdf("test_gviz.pdf") 
 plotTracks(
    list(itrack,ucscGenes2,dtrack_iPSC_WT,dtrack_iPSC_del,dtrack_iN_WT,dtrack_iN_del),
    from = st, to = en,
    transcriptAnnotation = "symbol"
  )
 dev.off()

#####Annotate with different gene names


library(data.table)
library(Gviz)
setwd("C:\\Documents and Settings\\ricar\\Downloads/CADM2")

directory <- getwd();

sampleFiles <- grep("i",list.files(directory, pattern = "\\.bedGraph"),value=TRUE)

myfiles = lapply(sampleFiles, fread, col.names = c('chromosome', 'start', 'end', 'value'))

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
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
plotTracks(list(ucscGenes, ucscGenes2), chromosome = "chr3", from = st, to = en, transcriptAnnotation = "symbol")


grtrack <- GeneRegionTrack(
  txdb,
  chromosome = thechr, start = st, end = en,
  showId = TRUE,
  name = "Gene Annotation"
)

for (stats_file in sampleFiles){
  
  temp <- read.table(stats_file)
  bedgraph_iPSC_wt <- fread(
    'iPSC_wt.bedGraph',
    col.names = c('chromosome', 'start', 'end', 'value')
  )
}


d_plots <- d_nested %>%
  mutate(plt=map(data, function(x){
    browser()
    ggplot(aes(Position, value), data=x) +
      geom_bar(stat="identity") +
      stat_smooth(aes(y=value), method = "gam",se=F,formula=y~s(x,k=7))
  }))

