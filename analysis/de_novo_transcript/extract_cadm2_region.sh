#!/bin/bash 
#SBATCH --job-name=single_cpu_example 
#SBATCH --partition=short 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=2G 
#SBATCH --output=log%J.out 
#SBATCH --error=log%J.err 

module load samtools

cd /data/talkowski/Samples/cadm2/data/de_novo/cadm2_capseq_bam
samtools view -b /data/talkowski/Samples/cadm2/data/Alignment_CADM2_capseq/iN_rA6/iN_rA6.Aligned.sortedByCoord.out.bam 3:84958546-86067835 > iN_rA6_capseq_trimmed.bam

cd /data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna
