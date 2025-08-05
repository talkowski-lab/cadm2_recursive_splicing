library(tidyverse)
library(UpSetR)
library(ComplexUpset)
setwd("D:/data/talkowski/Samples/fd_mouse_tissue/training_MED/output/final_results")
################## Make lists of up and down regulated DEGs

directory <- "."
sampleFiles <- grep("_DEG_stats",list.files(directory, pattern = "\\0.1.txt"),value=TRUE)

DEG_filter_file <- function(ExpTableName) {
  ExpTable <- read.table(ExpTableName)
  final_table <- ExpTable %>% 
    filter(pvalue < 0.1 ) %>%  
    #filter(log2FoldChange > log2(1.2) | log2FoldChange < log2(0.8))
  return(final_table)
}


cortex_stats <- DEG_filter_file(sampleFiles[1])
drg_stats <- DEG_filter_file(sampleFiles[2])
med_stats <- DEG_filter_file(sampleFiles[3])
sc_stats <- DEG_filter_file(sampleFiles[4])
tg_stats <- DEG_filter_file(sampleFiles[5])

cortex_DEGs=rownames(cortex_stats)
drg_DEGs=rownames(drg_stats)
med_DEGs=rownames(med_stats)
sc_DEGs=rownames(sc_stats)
tg_DEGs=rownames(tg_stats)

#library(ComplexUpset)

test <- list(DRG=rownames(drg_stats),
             TG=rownames(tg_stats),
             MED=rownames(med_stats),
             Cortex=rownames(cortex_stats),
             SC=rownames(sc_stats)
             )

library(ComplexHeatmap)
m1 <- make_comb_mat(test,mode = "intersect")
comb_name(m1)
m2=m1[comb_name(m1) %in% c("10000","01000","00100","00010","00001", "11000","00110","00101", "00011", "00111")]
m2=m2[,c("10000","01000","00100","00010","00001", "11000","00110","00101", "00011", "00111")]
pdf("all_tissues_DEGs_upsetR_q_0.1_LFC_20_percent.pdf")
UpSet(m2,set_order = c(1:5), comb_order = 1:10)
dev.off()

pdf("all_tissues_DEGs_upsetR_q_0.1_LFC_20_percent.pdf")
ComplexUpset::upset(fromList(test),c("Cortex","DRG","MED","SC","TG"),mode='inclusive_intersection', intersections='all', width_ratio = 0.00001, n_intersections=37,themes=upset_default_themes(text=element_text(size=25)))
dev.off()

all_tissue_genes <- intersect(drg_DEGs,intersect(tg_DEGs,intersect(med_DEGs,intersect(sc_DEGs,cortex_DEGs))))
write.table(all_tissue_genes, "all_tissue_genes_q_0.1_20_percent.txt", row.names=TRUE, quote=FALSE, col.names = TRUE)

cns <- list(Cortex=cortex_DEGs,
             MED=med_DEGs,
             SC=sc_DEGs)
pdf("CNS_tissues_DEGs_upsetR_q_0.1_20_percent.pdf")
cns_upset<-ComplexUpset::upset(fromList(cns),c("Cortex","MED","SC"),mode='inclusive_intersection', intersections='all', width_ratio = 0.00001, n_intersections=12,themes=upset_default_themes(text=element_text(size=25)))
dev.off()

cns_genes <-intersect(cortex_DEGs,intersect(med_DEGs,sc_DEGs))
write.table(cns_genes, "cns_genes_q_0.1_20_percent.txt", row.names=TRUE, quote=FALSE, col.names = TRUE)

pns <- list(DRG=drg_DEGs,
             TG=tg_DEGs)
pdf("PNS_tissues_DEGs_upsetR_q_0.1_20_percent.pdf")
pns_upset <- ComplexUpset::upset(fromList(pns),c("DRG","TG"),mode='inclusive_intersection', intersections='all', width_ratio = 0.00001, n_intersections=37,themes=upset_default_themes(text=element_text(size=25)))
dev.off()
pns_genes <- intersect(drg_DEGs,tg_DEGs)
write.table(pns_genes, "PNS_genes_q_0.1_20_percent.txt", row.names=TRUE, quote=FALSE, col.names = TRUE)

#######
library(UpSetR)


# example list
gene_set_list <- list(Cortex=cortex_DEGs,
                      DRG=drg_DEGs,
                      MED=med_DEGs,
                      SC=sc_DEGs,
                      TG=tg_DEGs
)
# make df from list 
all_genes <- Reduce(union, gene_set_list)
df <- sapply(all_genes,
             function(gene){
               memberships <- sapply(gene_set_list,
                                     function(gene_set){
                                       ifelse(gene %in% gene_set,
                                              1,
                                              0
                                       )
                                     }
               )
             }
) %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var='gene_name') %>%
  as_tibble() %>%
  rename_with(~ paste0('Category_', .),)
category_colnames <- df %>% select(starts_with('Category_')) %>% colnames()

################
category_colnames <- c("Gene","Cortex","DRG","MED","SC","TG")
colnames(df) <- c("Gene","Cortex","DRG","MED","SC","TG")

plt <- upset(df,
      category_colnames,
      mode='inclusive_intersection',  
      min_degree=2,  # ignore intersections with < 2 groups
      height_ratio=1, # make matrix same size as bar graph proportionally
      set_sizes=(  # include bar plot of group sizes (genes in each tissue)
        # left side of intersection matrix
        upset_set_size(position='left') +  
          # grey bar graph with white text labels
          geom_text(aes(label=..count..),
                    col='white',
                    size=5.5,
                    vjust=0.2,
                    hjust=-0.05,
                    stat='count'
          ) +
          # Vertical X axis labels
          theme(axis.text.x=element_text(size=25,
                                         angle=90
          ),
          axis.title.x=element_text(size=30)
          )
      ),
      matrix=(  # increase size of dots in intersection matrix
        intersection_matrix(
          geom=geom_point(size=6),
          segment=geom_segment(linetype='solid')
        )
      ),
      themes=upset_modify_themes(
        list(
          'intersections_matrix'=theme(
            axis.text=element_text(size=30, face='bold'),
            axis.title.x=element_text(size=0),
            axis.title=element_text(size=30),
          )
        )
      )
)

pdf("UpSetR.pdf", width = 12)
plt
dev.off()
