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

#PicardTools RNA alignment & quality metrics
$jav -Xmx4g -jar $pic ReorderSam \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.Aligned.sortedByCoord.out.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  REFERENCE=$genomeFA ## Reorder the .bam file according to the reference at hand

$jav -Xmx4g -jar $pic CollectAlignmentSummaryMetrics \
  REFERENCE_SEQUENCE=${genomeFA} \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/alignment_stats.txt \
  ASSUME_SORTED=false \
  ADAPTER_SEQUENCE=null ## Collect alignment metrics if the file is not present

$jav -Xmx4g -jar $pic CollectRnaSeqMetrics \
  REFERENCE_SEQUENCE=${genomeFA} \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/rnaseq_stats.txt \
  STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
  REF_FLAT=$refFlat \
  ASSUME_SORTED=false ## Collect sequencing metrics if the file is not present

$jav -Xmx4g -jar $pic CollectGcBiasMetrics \
  REFERENCE_SEQUENCE=${genomeFA} \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/gcbias_stats.txt \
  ASSUME_SORTED=false \
  CHART_OUTPUT=${rootdir}/data/STAR_bam/$name/gcbias_chart.pdf \
  SUMMARY_OUTPUT=${rootdir}/data/STAR_bam/$name/gcbias_summary.txt

$jav -Xmx4g -jar $pic CollectInsertSizeMetrics \
  REFERENCE_SEQUENCE=${genomeFA} \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/insert_size_metrics.txt \
  ASSUME_SORTED=false \
  HISTOGRAM_FILE=${rootdir}/data/STAR_bam/$name/insert_size_histogram.pdf

$jav -Xmx4g -jar $pic MarkDuplicates \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  METRICS_FILE=${rootdir}/data/STAR_bam/$name/duplication_stats.txt \
  ASSUME_SORTED=false \
  OUTPUT=${rootdir}/data/STAR_bam/$name/reordered_duplication_marked_reads.bam \
  REMOVE_DUPLICATES=TRUE

$jav -Xmx4g -jar $pic SortSam \
  INPUT=${rootdir}/data/STAR_bam/$name/reordered_duplication_marked_reads.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/reordered_duplication_marked_reads_sorted.bam \
  SORT_ORDER=coordinate ## Sort the de-duplicated file

$jav -Xmx4g -jar $pic BuildBamIndex \
  INPUT=${rootdir}/data/STAR_bam/$name/reordered_duplication_marked_reads_sorted.bam

rm ${rootdir}/data/STAR_bam/$name/reordered_duplication_marked_reads.bam
rm ${rootdir}/data/STAR_bam/$name/reordered_reads.bam



bam_stat.py -i ${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam > ${rootdir}/data/STAR_bam/$name/rseqc_bamstat.txt
clipping_profile.py \
  -i ${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  -s "PE" -o ${rootdir}/data/STAR_bam/$name/clipping_profile

geneBody_coverage.py \
-i ${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
-r $refBed -o ${rootdir}/data/STAR_bam/$name/gene_body

inner_distance.py \
-i ${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
-o ${rootdir}/data/STAR_bam/$name/inner_distance \
-r $refBed

read_distribution.py
-i ${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam
-r

#read_duplication.py -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -o output

#read_GC.py -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -o output

#tin.py -i -r
