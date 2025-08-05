#!/bin/bash
#SBATCH --job-name=iN_rsem_index
#SBATCH --partition=normal
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --output=log%J.out
#SBATCH --error=log%J.err

module load perl/5.30.3 bowtie2/2.4.2 samtools/1.11 rsem/1.3.3
cd /data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/iPSC_rsem_reference

rsem-prepare-reference --bowtie2 --transcript-to-gene-map /data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/rsem_index/rsem_index_Uniform/gene_trans.map /data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/rsem_index/cadm2_iN_all_transcripts.fa /data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/rsem_index/cadm2_iN_all_transcripts.fa
