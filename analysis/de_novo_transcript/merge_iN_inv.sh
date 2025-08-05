#!/bin/bash 
#SBATCH --job-name=single_cpu_example 
#SBATCH --partition=short 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=2G 
#SBATCH --output=log%J.out 
#SBATCH --error=log%J.err 
module load samtools java

cd /data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna
samtools merge /data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/iN_inv_cap_rna_merged.bam /data/talkowski/Samples/cadm2/data/de_novo/cadm2_capseq_bam/iN_rB10*_capseq_trimmed.bam /data/talkowski/Samples/cadm2/data/de_novo/cadm2_capseq_bam/iN_rC5*_capseq_trimmed.bam /data/talkowski/Samples/cadm2/data/de_novo/cadm2_rnaseq_bam/iN_*inv*trimmed.bam
samtools index iN_inv_cap_rna_merged.bam
java -jar -Xmx3g /apps/lib/picard/2.7.1/picard.jar MarkDuplicates INPUT=iN_inv_cap_rna_merged.bam OUTPUT=iN_inv_cap_rna_merged.rmdup.bam METRICS_FILE=iN_inv_cap_rna_merged.dup_metrics CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
