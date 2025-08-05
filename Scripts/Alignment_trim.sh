#!/bin/bash

###############################################
#
# Aligning bams and reprocessing them
# Rachita Yadav (ryadav1@mgh.harvard.edu)
# Talkowski Lab
#
###############################################
fdata=$1
prefix=$2
DIR=`pwd`
DATADIR=$DIR/$fdata
WORKDIR=$DIR/data/Alignment_${prefix}
#DIRS=`find $DATADIR -type "d" -name "SRR*"`
mkdir -p $DIR/Jobs
mkdir -p $WORKDIR
module load star/2.7.3

for file in `ls -1 $DATADIR/*.R1.fastq.trim.gz`
do

ID=`basename "$file"`
ID=${ID/.R1.fastq.trim.gz/}
#ID=${ID/_R1.fastq.gz/}
mkdir -p $WORKDIR/$ID

echo $DIR
echo $ID

bsub -sla miket_sc -q big-multi -n 8 -M 50000 -J ${ID} \
      -o Jobs/${ID}.out -e Jobs/${ID}.err \
      "STAR --runThreadN 8 \
      --genomeDir /data/talkowski/Samples/cadm2/data/ref/star_index \
      --twopassMode Basic \
      --outSAMunmapped Within \
      --outFilterMultimapNmax 1 \
      --outFilterMismatchNoverLmax 0.1 \
      --outSAMtype BAM SortedByCoordinate \
      --readFilesCommand zcat \
      --alignEndsType EndToEnd \
      --alignIntronMin 21 \
      --alignIntronMax 0 \
      --quantMode GeneCounts \
      --outFileNamePrefix $WORKDIR/${ID}/${ID}. \
      --readFilesIn $DATADIR/${ID}.R1.fastq.trim.gz $DATADIR/${ID}.R2.fastq.trim.gz"
       ######     --readFilesIn $DATADIR/${ID}_R1.fastq.gz $DATADIR/${ID}_R2.fastq.gz"
    
done 
