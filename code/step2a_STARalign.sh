#!/bin/bash


STAR_call=/home/mgandal/bin/STAR-2.5.2a/bin/Linux_x86_64/STAR
genomeFA=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/GRCm38.p4.genome.fa
genomeDir=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/STAR_index
gtfFile=/hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.gtf
rootdir=/hp_shares/mgandal/projects/TSC_MIA_Silva

cd ${rootdir}/data/fastq_merged

for file in *_R1_001.fastq.gz; do

name=`basename $file _R1_001.fastq.gz`
if ! [ -s ${rootdir}/STAR_sam/${name}*.sam ]; then
   echo ${name}_R1_001.fastq.gz ${name}_R1_002.fastq.gz
   qsub -o ${rootdir}/code/log -e ${rootdir}/code/log -l h_rt=8:00:00,h_data=25G,highmem -pe shared 4  ${rootdir}/code/step2b_qsub_starCall.sh ${name}
fi
done
