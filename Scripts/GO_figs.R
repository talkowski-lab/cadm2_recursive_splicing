library(tidyverse)
library(ggplot2)
dat_table <- read.table("C:\\Users\\ricar\\OneDrive\\Documents\\SynGOBP.csv", header = T,sep=",")

pns <- dat_table %>%
  filter(Tissue %in% c("DRG","TG")) %>%
  mutate(Description=factor(Description, levels=unique(Description))) %>%
  mutate(Tissue=factor(Tissue, levels=c("DRG", "TG")))
pns$Description=as.vector(pns$Description)
pns$Description[9] <- 'Regulation of translation at presynapse, modulating synaptic transmission1'
pns$Description=factor(pns$Description,levels = rev(as.vector(pns$Description)))

cns <- dat_table %>%
  filter(Tissue %in% c("Cortex","MED", "SC")) %>%
  mutate(Description=factor(Description, levels=unique(Description))) %>%
  mutate(Tissue=factor(Tissue, levels=c("Cortex", "MED", "SC")))
cns$Description=as.vector(cns$Description)
cns$Description[4] <- 'Process in the presynapse1'
cns$Description[5] <- 'Process in the synapse1'
cns$Description[7] <- 'Regulation of presynaptic cytosolic calcium levels1'
cns$Description[9] <- 'Synaptic vesicle exocytosis1'
cns$Description=factor(cns$Description,levels = rev(as.vector(cns$Description)))

convergence <- dat_table %>%
  filter(Tissue == "PNS") %>%
  mutate(Description=factor(Description, levels=unique(Description))) %>%
  mutate(Tissue=factor(Tissue, levels=c("PNS")))
convergence$Description=as.vector(convergence$Description)
convergence$Description=factor(convergence$Description,levels = rev(as.vector(convergence$Description)))


pdf("cns_syngo_BP.pdf")  
ggplot(cns, aes(x=-log10(FDR.q.value), y=Description,fill=Tissue)) +
  geom_bar(position = "dodge",width=0.5,stat="identity") +
  geom_vline(xintercept=-log10(0.05), linetype="dashed",color="black") +
  theme(text = element_text(size = 20))   +
  facet_wrap(vars(Tissue), scales="free_y", nrow=3, ncol=1) +
  theme(legend.position="none") +
  xlab("-log10(p-value)") +
  ylab("Biological Process")
dev.off()


##############Module plotting##############
library(tidyverse)
library(ggplot2)
dat_table <- read.table("C:\\Users\\ricar\\OneDrive\\Documents\\GO_high_cor_mod.csv", header = T,sep=",")

syngo_annotation <- dat_table %>%
  filter(GO == "CC") %>%
  mutate(Description=factor(Description, levels=unique(Description))) %>%
  mutate(Tissue=factor(Tissue, levels=c("DRG Black", "DRG Salmon","DRG Brown","TG Cyan", "TG Red","MED Turquoise", "MED Antiquewhite4")))

######CC
syngo_annotation$Description=as.vector(syngo_annotation$Description)
syngo_annotation$Description[1] <- 'Integral component of postsynaptic membrane1'
syngo_annotation$Description[9] <- 'Integral component of postsynaptic membrane2'
syngo_annotation$Description[6] <- 'Integral component of presynaptic active zone membrane1'
syngo_annotation$Description[23] <- 'Integral component of presynaptic membrane1'
syngo_annotation$Description[2] <- 'Postsynapse1'
syngo_annotation$Description[13] <- 'Postsynapse2'
syngo_annotation$Description[14] <- 'Postsynaptic density1'
syngo_annotation$Description[16] <- 'Postsynaptic density2'
syngo_annotation$Description[3] <- 'Postsynaptic membrane1'
syngo_annotation$Description[18] <- 'Postsynaptic specialization1'
syngo_annotation$Description[12] <- 'Postsynaptic specialization2'
syngo_annotation$Description[15] <- 'Presynapse1'
syngo_annotation$Description[8] <- 'Presynaptic active zone1'
syngo_annotation$Description[4] <- 'Presynaptic membrane1'
syngo_annotation$Description[24] <- 'Presynaptic membrane2'
syngo_annotation$Description[2] <- 'Synapse1'
syngo_annotation$Description[11] <- 'Synapse2'
syngo_annotation$Description[31] <- 'Synapse3'
syngo_annotation$Description[5] <- 'Synapse4'

#####BP
syngo_annotation$Description=as.vector(syngo_annotation$Description)
syngo_annotation$Description[4] <- 'Chemical synaptic transmission1'
syngo_annotation$Description[9] <- 'Postsynapse organization1'
syngo_annotation$Description[7] <- 'Presynaptic process involved in chemical synaptic transmission1'
syngo_annotation$Description[30] <- 'Presynaptic process involved in chemical synaptic transmission2'
syngo_annotation$Description[23] <- 'Process in the postsynapse1'
syngo_annotation$Description[24] <- 'Process in the synapse1'
syngo_annotation$Description[12] <- 'Process in the synapse2'
syngo_annotation$Description[29] <- 'Synaptic signaling1'
syngo_annotation$Description[1] <- 'Trans-synaptic signaling1'


syngo_annotation$Description=factor(syngo_annotation$Description,levels = rev(as.vector(syngo_annotation$Description)))


pdf("modules_high_cor_syngo_cc.pdf")  
ggplot(syngo_annotation, aes(x=-log10(FDR.q.value), y=Description,fill=Tissue)) +
  geom_bar(position = "dodge",width=0.5,stat="identity") +
  geom_vline(xintercept=-log10(0.05), linetype="dashed",color="black") +
  theme(text = element_text(size = 10))   +
  facet_wrap(vars(Tissue), scales="free_y", nrow=8, ncol=1) +
  theme(legend.position="none") +
  theme_classic() +
  xlab("-log10(p-value)") +
  ylab("Biological Process")
dev.off()


##########make multi tissue GO plot
library(tidyverse)
library(tidytext)
library(ggplot2)
dat_table <- read.table("C:\\Users\\ricar\\OneDrive\\Documents\\condensed_GO.csv", header = T,sep=",")
############
tissue_go <- dat_table %>%
  filter(Keep =="0") %>%
  mutate(Description=factor(Description, levels=unique(Description)))
tissue_go$Description=as.vector(tissue_go$Description)
##############Trying to reorder_within

tissue_go <- dat_table %>%
  filter(Keep =="0") %>%
  mutate(Description=factor(Description), Description = reorder_within(Description, FDR.q.value, Tissue))
tissue_go$Description=as.vector(tissue_go$Description)
##############

tissue_go$Description=as.vector(tissue_go$Description)
tissue_go$Tissue <- factor(tissue_go$Tissue,      # Reordering group factor levels
                           levels = c("DRG Black", "DRG Brown", "DRG Salmon", "TG Cyan", "TG Red","MED Turquoise", "MED Antiquewhite4"))
pdf("module_high_cor_go.pdf", width = 25)  
ggplot(tissue_go, aes(x=-log10(FDR.q.value), y=Description,fill=Database)) +
  geom_bar(position = "dodge",width=0.5,stat="identity") +
  #geom_vline(xintercept=-log10(0.05), linetype="dashed",color="black") +
  theme_classic() +
  theme(text = element_text(size = 9))   +
  facet_wrap(vars(Tissue), nrow = 1, ncol = 7) +
  theme(panel.grid.major.y =  element_line(color = "gray",
                                           size = 0.5,
                                           linetype = 2)) +
  xlab("-log10(p-value)") +
  ylab("Biological Process")
dev.off()




library(tidyverse)
library(ggplot2)
dat_table <- read.table("C:\\Users\\ricar\\OneDrive\\Documents\\combined_mods_go.csv", header = T,sep=",")
tissue_go <- dat_table %>%
  filter(Keep =="0") %>%
  mutate(Description=factor(Description, levels=unique(Description)))
tissue_go$Description=as.vector(tissue_go$Description)
tissue_go$Tissue <- factor(tissue_go$Tissue,      # Reordering group factor levels
                           levels = c("DRG", "TG", "MED"))
pdf("module_combined_high_cor_go.pdf", width = 25)  
ggplot(tissue_go, aes(x=sort(-log10(FDR.q.value)), y=Description,fill=Database)) +
  geom_bar(position = "dodge",width=0.5,stat="identity") +
  #geom_vline(xintercept=-log10(0.05), linetype="dashed",color="black") +
  theme_classic() +
  theme(text = element_text(size = 9))   +
  facet_wrap(vars(Tissue), nrow = 1, ncol = 7) +
  theme(panel.grid.major.y =  element_line(color = "gray",
                                           size = 0.5,
                                           linetype = 2)) +
  xlab("-log10(p-value)") +
  ylab("Biological Process")
dev.off()


#################
dat_table <- read.table("C:\\Users\\ricar\\OneDrive\\Documents\\tissue_GO.csv", header = T,sep=",")

tissue_go <- dat_table %>%
  filter(Keep==0) %>%
  filter(Tissue %in% c("SC")) %>%
  filter(Database == "GSEA") %>%
  arrange(desc(FDR.q.value))%>%
  mutate(Description=factor(Description, levels=unique(Description)))
tissue_go$Description=as.vector(tissue_go$Description)
tissue_go$Tissue <- factor(tissue_go$Tissue,      # Reordering group factor levels
                           levels = c("SC"))

pdf("tissue_go.pdf", width = 25)  
c <- ggplot(tissue_go, aes(x=-log10(FDR.q.value), y=sort(Description))) +
  geom_bar(position = "dodge",width=0.5,stat="identity") +
  geom_vline(xintercept=-log10(0.05), linetype="dashed",color="black") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25)) +
  theme(text = element_text(size = 16), legend.position = "none")   +
  facet_grid(~Tissue) +
  theme(panel.grid.major.y =  element_line(color = "gray",
                                        size = 0.5,
                                        linetype = 2)) +
  xlim(0,6) +
  xlab("-log10(p-value)") +
  ylab("Biological Process") +
  theme(aspect.ratio = 0.3) +
  theme_classic()
dev.off()

pdf("tissue_go_cns_narrow.pdf", width = 18, height = 14)
grid.arrange(a,b,c,ncol=3)
dev.off()

dat_table <- read.table("C:\\Users\\ricar\\OneDrive\\Documents\\pns_GO.csv", header = T,sep=",")
tissue_go <- dat_table %>%
  filter(Database=="GSEA") %>%
  mutate(Description=factor(Description, levels=unique(Description)))
tissue_go$Description=as.vector(tissue_go$Description)

tissue_go$Description <- factor(tissue_go$Description,levels = rev(as.vector(tissue_go$Description)))

pdf("pns_go_gsea.pdf", width = 10)  
ggplot(tissue_go, aes(x=-log10(FDR.q.value), y=Description,fill=Database)) +
  geom_bar(position = "dodge",width=0.5,stat="identity") +
  geom_vline(xintercept=-log10(0.05), linetype="dashed",color="black") +
  theme_classic() +
  theme(text = element_text(size = 20), legend.position = "None")   +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25)) +
  facet_wrap(vars(Tissue), scales='free_x', nrow = 1, ncol = 5) +
  xlab("-log10(p-value)") +
  ylab("Biological Process")
dev.off()

#####################
library(tidyverse)
library(ggplot2)
dat_table <- read.table("C:\\Users\\ricar\\OneDrive\\Documents\\combined_dose_response_mods_go.csv", header = T,sep=",")
tissue_go <- dat_table %>%
  filter(Keep =="0") %>%
  filter(Database=="GSEA") %>%
  arrange(desc(FDR.q.value))%>%
  mutate(Description=factor(Description, levels=unique(Description)))
tissue_go$Description=as.vector(tissue_go$Description)
tissue_go$Tissue <- factor(tissue_go$Tissue,      # Reordering group factor levels
                           levels = c("DRG", "TG", "MED"))
pdf("elp1_dose_response1.pdf", width = 11, height = 7)  
ggplot(tissue_go, aes(x=-log10(FDR.q.value), y=sort(Description),fill=Database)) +
  geom_bar(position = "dodge",width=0.85,stat="identity") +
  geom_vline(xintercept=-log10(0.05), linetype="dashed",color="black") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  theme_classic() +
  theme(text = element_text(size = 16), legend.position = "none")   +
  facet_wrap(vars(Tissue), nrow = 1, ncol = 7) +
  theme(panel.grid.major.y =  element_line(color = "gray",
                                           size = 0.5,
                                           linetype = 2)) +
  xlab("-log10(p-value)") +
  ylab("Biological Process") +
  theme_classic()
dev.off()
