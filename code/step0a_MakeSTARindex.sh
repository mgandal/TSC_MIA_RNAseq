#!/bin/bash

STAR_call=/home/mgandal/bin/STAR-2.5.2a/bin/Linux_x86_64/STAR

qrsh -l h_rt=24:00:00,h_data=80G -pe shared 8
$STAR_call
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/STAR_index \
  --genomeFastaFiles /hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/GRCm38.p4.genome.fa \
  --sjdbGTFfile /hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.gtf \
  --sjdbOverhang 49
