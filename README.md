# CADM2 Transcriptomic Analysis

This repository contains scripts, data, and analyses related to transcriptomic profiling of CRISPR edits for the first recursive splice site in intron 1 of `CADM2` gene across hiPSCs and hiPSC-derived neurons. The work includes isoform-specific quantification, recursive splicing analysis, exon usage, and functional enrichment due to heterozygous and homozygous deletions and inversion events.

---

## Directory Structure

### `data/`
Contains differential expression results for various genomic perturbations:
- `iN_DEG_stats_WT_vs_*.txt`: Differential expression stats for induced neurons (iNs) comparing wild-type to different zygosity and structural variation conditions.

---

### `Scripts/`
Core R and shell scripts for alignment, transcript quantification, plotting, and statistical testing:
- `celltype_markers.R`, `CADM2_heatmap.R`: Marker gene expression and heatmap generation.
- `pheatmap.R`, `plot_GO_CADM2.R`, `GO_figs.R`: Gene ontology and visualization utilities.
- `iPSC_transcript_pvalues+plots.R`, `transcript_analysis_final.R`: Isoform differential analysis and final transcript summaries.
- `Alignment_trim.sh`, `dexseq_count.sh`: Shell scripts for alignment preprocessing and exon counting.
- `ComplexUpset_plotting.R`: Visualize overlaps across DEG sets.

---

### `analysis/deg_analysis/`
Differential expression (DEG) analyses across multiple comparisons:
- `iN_DEG_WT_vs_HET_zygosity_single_comparison.R`: Scripts stratifying zygosity effects.
- `iN_HET_*`, `iN_HOM_*`: Subgroup analyses for HET and HOM comparisons.
- `iPSC_DEG_WT_vs_*`: DEG analyses in iPSC-derived progenitor cells under different structural variants.

---

### `analysis/de_novo_transcript/`
Pipeline scripts for novel transcript discovery and quantification:
- `trinity_*`: Trinity-based transcriptome assembly for various cell types and genotypes.
- `rsem_index.py`, `rsem_quant_p3.py`: Scripts for transcript quantification using RSEM.
- `validate_transcripts_iN_WT.py`: Cross-validation of novel transcripts.

---

### `analysis/exon_analysis/`
DEXSeq-based analysis of alternative exon usage:
- `DEXSeq_perType.R`, `Exon_plots.R`: R scripts for differential exon usage analysis.
- `dexseq_count.sh`: Shell script for DEXSeq-compatible exon count generation.

---

### `analysis/recursive_splicing_quantification/`
Scripts for identifying and quantifying recursive splicing events:
- `RS_slope_calculator.R`: Calculates RS scores.
- `RNA_seq_*_concensus_to_bedGraph.R`: Generates bedGraph signal files from consensus splice signals.

---

### `figures/`
Organized scripts for each main and supplementary figure in the manuscript:
- `fig1/`, `fig2/`, `fig3/`, `fig4/`: R scripts and plotting utilities for specific figure panels.
- `fig_supplmentary/`: Supplementary figure generation including heatmaps and cell type marker plots.

---

### `refs/genomic_ref/`
Genomic reference files used throughout the project:
- `gencode.v47.transcripts.fa`, `GRCh38.primary_assembly.genome.fa.gz`: Reference genome and transcriptome.
- `cadm2_exons.bed`, `CADm2_canonical.bed`: CADM2-specific annotations and BED files.

---

### `Scripts/data_for_plots/`
Intermediate and final output tables for figure plotting:
- Read count matrices, normalized expression tables, exon statistics, and significance tables.

---

##  Requirements
- R (= 4.0)
- Python (= 3.6)
- Trinity, RSEM, STAR, HTSeq, DEXSeq, ComplexUpset
- Recommended packages: `ggplot2`, `pheatmap`, `DESeq2`, `Bioconductor`, `dplyr`, `tidyr`

---

##  License
*Add your license here (e.g., MIT, GPLv3, etc.)*

---

##  Contact
*Add your contact information or lab page link here*

---




