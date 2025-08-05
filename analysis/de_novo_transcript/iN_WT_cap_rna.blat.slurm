#!/bin/bash
#SBATCH --job-name=iN_WT_cap_rna_merged_trinity/
#SBATCH --partition=normal
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --output=log%J.out
#SBATCH --error=log%J.err
cd /data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/iN_WT_cap_rna_merged_trinity/
module load perl/5.30.3 GCCcore/7.3.0 zlib/1.2.11-GCCcore-7.3.0 libpng/1.6.34-GCCcore-7.3.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/PHShome/ry077/lib
/PHShome/ry077/bin/blat /data/talkowski/Samples/cadm2/data/de_novo/trinity_bam/ensembl_chr3.fa Trinity_1line.fasta -fine raw.psl
echo iN_WT_cap_rna_merged_trinity Completed
