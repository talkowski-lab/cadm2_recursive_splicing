library(ggplot2)
library(tidyverse)
library(ggpubr)

dat <- read.table("C:\\Documents and Settings/ricar/Downloads/CADM2/plotting_cadm2_go.txt",sep='\t', header = TRUE)

dat_gsea_hom_del <- dat %>% filter(DB=="GSEA" & Genotype=="HOM_DEL")
dat_syngo_hom_del <- dat %>% filter(DB=="SynGO" & Genotype=="HOM_DEL")
dat_gsea_het_del <- dat %>% filter(DB=="GSEA" & Genotype=="HET_100KB")
dat_syngo_het_del <- dat %>% filter(DB=="SynGO" & Genotype=="HET_100KB")
dat_gsea_convergence <- dat %>% filter(DB=="GSEA" & Genotype=="CONVERGENCE")
dat_syngo_convergence <-dat %>% filter(DB=="SynGO" & Genotype=="CONVERGENCE")

dat_GSEA <- dat %>% filter(DB=="GSEA") %>% filter(Term %in% c("GOBP_NEUROGENESIS",
                                                              "GOBP_NEURON_DIFFERENTIATION",
                                                              "GOBP_NEURON_DEVELOPMENT",
                                                              "GOBP_AXON_DEVELOPMENT",
                                                              "GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION",
                                                              "GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_NEURON_DIFFERENTIATION",
                                                              "GOBP_REGULATION_OF_MULTICELLULAR_ORGANISMAL_DEVELOPMENT",
                                                              "GOBP_REGULATION_OF_NERVOUS_SYSTEM_DEVELOPMENT",
                                                              "GOBP_CELL_PART_MORPHOGENESIS"))
dat_GSEA$Term <- gsub("_"," ",dat_GSEA$Term)
dat_GSEA$Term <- gsub("GOBP","",dat_GSEA$Term)

p <- ggplot(dat_GSEA, aes(x =reorder( Term,-log10(FDR)), y = -log10(FDR), fill=Genotype)) +
  geom_bar(stat="identity", width = 0.7, position = "dodge2") +
  coord_flip() + 
  xlab("Term") + ylab("-log10(FDR)") +
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1)) + 
  theme_classic() +
  geom_hline(yintercept = -log10(0.1))
  theme(text = element_text(size =30),legend.position = "bottom") 

# Horizontal bar plot
pdf("GSEA_GO_all.pdf",width=9, height=8)
p
dev.off()


