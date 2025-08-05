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
samtools merge /data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/iN_het_del_100_cap_rna_merged.bam /data/talkowski/Samples/cadm2/data/de_novo/cadm2_capseq_bam/iN_siB6*_capseq_trimmed.bam /data/talkowski/Samples/cadm2/data/de_novo/cadm2_capseq_bam/iN_siC1*_capseq_trimmed.bam /data/talkowski/Samples/cadm2/data/de_novo/cadm2_capseq_bam/iN_soD8*_capseq_trimmed.bam /data/talkowski/Samples/cadm2/data/de_novo/cadm2_capseq_bam/iN_soF6*_capseq_trimmed.bam /data/talkowski/Samples/cadm2/data/de_novo/cadm2_rnaseq_bam/iN_*100kb_het_del*trimmed.bam
samtools index iN_het_del_100_cap_rna_merged.bam
java -jar -Xmx3g /apps/lib/picard/2.7.1/picard.jar MarkDuplicates INPUT=iN_het_del_100_cap_rna_merged.bam OUTPUT=iN_het_del_100_cap_rna_merged.rmdup.bam METRICS_FILE=iN_het_del_100_cap_rna_merged.bam.dup_metrics CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
