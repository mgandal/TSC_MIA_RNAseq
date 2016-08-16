#!/bin/bash

SAMTOOLS_call=/home/mgandal/bin/samtools-1.3/samtools
STAR_call=/home/mgandal/bin/STAR-2.5.2a/bin/Linux_x86_64/STAR
jav=/share/apps/jre1.8.0_92/bin/java
pic=/home/mgandal/bin/picard.jar
genomeFA=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/GRCm38.p4.genome.fa
genomeDir=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/STAR_index
gtfFile=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.gtf
rootdir=/hp_shares/mgandal/projects/TSC_MIA_Silva
refFlat=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.refFlat.txt.gz
refBed=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.bed


#qrsh -l h_rt=24:00:00,h_data=80G -pe shared 8
$STAR_call \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/STAR_index \
  --genomeFastaFiles /hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/GRCm38.p4.genome.fa \
  --sjdbGTFfile /hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.gtf \
  --sjdbOverhang 49

##Check Basic Directory Strcuture
#/data/fastq_merged
#/data/STAR_bam
#/data/QC_picard
