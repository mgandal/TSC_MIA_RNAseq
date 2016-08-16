#!/bin/bash

STAR_call=/home/mgandal/bin/STAR-2.5.2a/bin/Linux_x86_64/STAR
genomeFA=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/GRCm38.p4.genome.fa
genomeDir=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/STAR_index
gtfFile=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.gtf
rootdir=/hp_shares/mgandal/projects/TSC_MIA_Silva

name=$1

$STARcall \
--runThreadN 2 \
--genomeDir ${genomeDir} \
--outFileNamePrefix ${rootdir}/data/STAR_sam/$name \
--readFilesCommand gunzip -c \
--readFilesIn ${rootdir}/data/fastq_merged/${name}_R1_001.fastq.gz.fastq.gz ${rootdir}/data/fastq_merged/${name}_R2_001.fastq.gz.fastq.gz \
--outSAMtype BAM Unsorted SortedByCoordinate
