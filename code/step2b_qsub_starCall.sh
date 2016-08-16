#!/bin/bash

SAMTOOLS_call=/home/mgandal/bin/samtools-1.3/samtools
STAR_call=/home/mgandal/bin/STAR-2.5.2a/bin/Linux_x86_64/STAR
jav=/share/apps/jre1.8.0_92/bin/java
pic=/home/mgandal/bin/picard.jar
genomeFA=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/GRCm38.p4.genome.fa
genomeDir=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/STAR_index
gtfFile=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.gtf
rootdir=/hp_shares/mgandal/projects/TSC_MIA_Silva


name=$1

if [ ! -s $rootdir/data/STAR_bam/$name ]; then mkdir $rootdir/data/STAR_bam/$name; fi

#STAR Alignment to Genome
$STAR_call \
--runThreadN 2 \
--genomeDir ${genomeDir} \
--outFileNamePrefix ${rootdir}/data/STAR_bam/$name/${name}. \
--readFilesCommand gunzip -c \
--readFilesIn ${rootdir}/data/fastq_merged/${name}_R1_001.fastq.gz ${rootdir}/data/fastq_merged/${name}_R2_001.fastq.gz \
--outSAMtype BAM Unsorted SortedByCoordinate

#Samtools sorting and indexing of BAM file
#$SAMTOOLS_call sort ${rootdir}/data/STAR_sam/$name -n -o ${rootdir}/data/STAR_bam/${name}.nameSorted.bam
$SAMTOOLS_call index ${rootdir}/data/STAR_bam/$name/${name}.Aligned.sortedByCoord.out.bam
$SAMTOOLS_call idxstats ${rootdir}/data/STAR_bam/$name/${name}.Aligned.sortedByCoord.out.bam

#PicardTools RNA alignment & quality metrics
$jav -Xmx4g -jar $pic ReorderSam \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.Aligned.sortedByCoord.out.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  REFERENCE=$genomeFA ## Reorder the .bam file according to the reference at hand
