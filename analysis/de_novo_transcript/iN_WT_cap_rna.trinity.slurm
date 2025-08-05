#!/bin/bash
#SBATCH --job-name=iN_WT_cap_rna_merged.rmdup.bam
#SBATCH --partition=normal
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --output=log%J.out
#SBATCH --error=log%J.err
cd /data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna
conda activate ML
module load java samtools salmon jellyfish bowtie2 trinity
Trinity --SS_lib_type RF --min_contig_length 100 --genome_guided_max_intron 100000 --max_memory 64G --CPU 8 --genome_guided_bam iN_WT_cap_rna_merged.rmdup.bam --output /data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/iN_WT_cap_rna_merged_trinity
