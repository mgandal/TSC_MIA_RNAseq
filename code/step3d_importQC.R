picard.align <- read.table("../data/QC/RNAseqAlign.txt",skip=1)
rownames(picard.align) <- picard.align[,1]
picard.align <- picard.align[,-c(1)]

colnames(picard.align) <- c("CATEGORY","TOTAL_READS","PF_READS","PCT_PF_READS","PF_NOISE_READS","PF_READS_ALIGNED","PCT_PF_READS_ALIGNED","PF_ALIGNED_BASES","PF_HQ_ALIGNED_READS","PF_HQ_ALIGNED_BASES","PF_ HQ_ALIGNED_Q20_BASES","PF_HQ_MEDIAN_MISMATCHES","PF_MISMATCH_RATE","PF_HQ_ERROR_RATE","PF_INDEL_RATE","MEAN_READ_LENGTH","READS_ALIGNED_IN_PAIRS","PCT_READS_ALIGNED_IN_PAIRS","BAD_CYCLES","STRAND_BALANCE","PCT_CHIMERAS","PCT_ADAPTER")

#summary(picard.align)
perc.hq <- picard.align[,"PF_HQ_ALIGNED_BASES"]/picard.align[,"PF_ALIGNED_BASES"]


## Transcript level QC
picard.rnaseq.qc <- read.table("../data/QC/RNAseqQC.txt",skip=1)
colnames(picard.rnaseq.qc) <- c("FileName","PF_BASES","PF_ALIGNED_BASES","CODING_BASES","UTR_BASES","INTRONIC_BASES","INTERGENIC_BASES","IGNORED_READS","CORRECT_STRAND_READS","INCORRECT_STRAND_READS",
                                "PCT_CODING_BASES","PCT_UTR_BASES", "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES", "PCT_USABLE_BASES", "PCT_CORRECT_STRAND_READS", "MEDIAN_CV_COVERAGE",
                                "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS", "MEDIAN_5PRIME_TO_3PRIME_BIAS")
                                
                                
                                #1st Qu.:1.020e+09   1st Qu.:7.851e+08   1st Qu.: 74684310
                                # Median :1.137e+09   Median :9.493e+08   Median :105254783
                                # Mean   :1.117e+09   Mean   :9.070e+08   Mean   :124770007
                                # 3rd Qu.:1.306e+09   3rd Qu.:1.056e+09   3rd Qu.:154292913
                                # Max.   :1.646e+09   Max.   :1.299e+09   Max.   :430650412
                                #
                                # INTERGENIC_BASES    IGNORED_READS CORRECT_STRAND_READS INCORRECT_STRAND_READS
                                # Min.   :        0   Min.   :0     Min.   :0            Min.   :0
                                # 1st Qu.: 47638852   1st Qu.:0     1st Qu.:0            1st Qu.:0         BASES PCT_UTR_BASES    PCT_INTRONIC_BASES PCT_INTERGENIC_BASES
                                # Min.   :0.3640   Min.   :0.2805   Min.   :0.00000    Min.   :0.00000
                                # 1st Qu.:0.4791   1st Qu.:0.3887   1st Qu.:0.03381    1st Qu.:0.02270
                                # Median :0.5072   Median :0.4051   Median :0.04643    Median :0.02840
                                # Mean   :0.5037   Mean   :0.4104   Mean   :0.05563    Mean   :0.03026
                                # 3rd Qu.:0.5355   3rd Qu.:0.4331   3rd Qu.:0.06650    3rd Qu.:0.03389
                                # Max.   :0.6158   Max.   :0.5222   Max.   :0.16635    Max.   :0.16281
                                #
                                # PCT_MRNA_BASES   PCT_USABLE_BASES PCT_CORRECT_STRAND_READS MEDIAN_CV_COVERAGE
                                # Min.   :0.6708   Min.   :0.6708   Min.   :0                Min.   :0.0000
                                # 1st Qu.:0.9002   1st Qu.:0.9001   1st Qu.:0                1st Qu.:0.4666
                                # Median :0.9255   Median :0.9254   Median :0                Median :0.5073
                                # Mean   :0.9141   Mean   :0.9140   Mean   :0                Mean   :0.5349
                                # 3rd Qu.:0.9427   3rd Qu.:0.9425   3rd Qu.:0                3rd Qu.:0.5732
                                # Max.   :1.0000   Max.   :1.0000   Max.   :0                Max.   :0.9846
                                #
                                # MEDIAN_5PRIME_BIAS MEDIAN_3PRIME_BIAS MEDIAN_5PRIME_TO_3PRIME_BIAS
                                # Min.   :0.0000     Min.   :0.0000     Min.   :0.0000
                                # 1st Qu.:0.3602     1st Qu.:0.6674     1st Qu.:0.4938
                                # Median :0.4208     Median :0.7137     Median :0.6048
                                # Mean   :0.4304     Mean   :0.7224     Mean   :0.5914
                                # 3rd Qu.:0.5055     3rd Qu.:0.7481     3rd Qu.:0.7122
                                # Max.   :0.6967     Max.   :1.0695     Max.   :0.9045
                               
#quantile(picard.rnaseq.qc[,"MEDIAN_5PRIME_TO_3PRIME_BIAS"],c(seq(0,1,0.05)))
#quantile(picard.rnaseq.qc[,"MEDIAN_5PRIME_TO_3PRIME_BIAS"],c(0.025,0.975))
#quantile(picard.rnaseq.qc[,"PF_BASES"]/100) ## gives the number of fragments aligned - 38.9 million is the median
#quantile(picard.rnaseq.qc[,"CODING_BASES"]/100) ## gives the number of fragments aligned to coding regions - 8.5 million on average
#quantile(picard.rnaseq.qc[,"INTRONIC_BASES"]/100) ## gives the number of fragments aligned to intronic regions - 15 million on average
#quantile(picard.rnaseq.qc[,"INTERGENIC_BASES"]/100) ## gives the number of fragments aligned to intronic regions - 2.7 million on average
           
                     
tx.coverage <- read.table("../data/QC/TranscriptCoverage.txt") ## Need to go into this file and manually remove blank lines: HSB155.V1C_SE.bam was blank
keep <- match(unique(tx.coverage[,1]),tx.coverage[,1])
tx.coverage <- tx.coverage[keep,]
rownames(tx.coverage) <- tx.coverage[,1]
tx.coverage.cols <- as.character(tx.coverage[1,seq(4,204,by=2)])
tx.coverage <- tx.coverage[,seq(5,205,by=2)]
colnames(tx.coverage) <- tx.coverage.cols
tx.coverage.quant <- apply(tx.coverage,2,quantile,c(0.025,0.5,0.975)) #error bar


par(mfrow=c(1,1))
boxplot(cbind(perc.hq,picard.rnaseq.qc[,c(15,11:14)]),ylab="Percent of bases",main="RNA-seq QC metrics across samples",names=c("High\nQuality","mRNA\nBases","Protein\nCoding","Untranslated\nRegion","Intronic\nRegion","Intergenic\nRegion"),ylim=c(0,1))

plot(x=as.numeric(colnames(tx.coverage.quant)),y=tx.coverage.quant[2,],xlab="Percentile of gene body (5' -> 3')",ylab="Coverage relative to whole transcript",pch=19,ylim=c(0,1.7),main="Relative transcript coverage - median with 95% CIs across samples")
lines(x=as.numeric(colnames(tx.coverage.quant)),y=tx.coverage.quant[2,],col="black")
lines(x=as.numeric(colnames(tx.coverage.quant)),y=tx.coverage.quant[1,],col="grey")
lines(x=as.numeric(colnames(tx.coverage.quant)),y=tx.coverage.quant[3,],col="grey")
                               
                                
#GC content
# Columns c("GC","WINDOWS","READ_STARTS","MEAN_BASE_QUALITY","NORMALIZED_COVERAGE","ERROR_BAR_WIDTH")
# tmp <- read.table("../data/QC/RNAseqGC.txt",skip=1)
# picard.rnaseq.gc <- matrix(NA,nrow=nrow(tmp),ncol=100)
# colnames(picard.rnaseq.gc) <- as.character(tmp[1,seq(4,801,by=8)])
# rownames(picard.rnaseq.gc) <- tmp[,1]
# picard.rnaseq.gc <- tmp[,seq(8,801,by=8)]
# gc.coverage.quant <- apply(picard.rnaseq.gc,2,quantile,c(0.025,0.5,0.975))
# ###colnames(gc.coverage.quant) <- as.character(tmp[1,seq(2,601,by=6)])

#GC quality
# picard.rnaseq.gcqual <- matrix(NA,nrow=nrow(tmp),ncol=100)
# colnames(picard.rnaseq.gcqual) <- tmp[1,seq(2,601,by=6)]
# rownames(picard.rnaseq.gcqual) <- tmp[,1]
# picard.rnaseq.gcqual <- tmp[,seq(5,601,by=6)]
# gc.qual.quant <- apply(picard.rnaseq.gcqual,2,quantile,c(0.025,0.5,0.975))
# picard.rnaseq.gcbin <- matrix(NA,nrow=nrow(tmp),ncol=100)
# colnames(picard.rnaseq.gcbin) <- tmp[1,seq(2,601,by=6)]
# rownames(picard.rnaseq.gcbin) <- tmp[,1]
# picard.rnaseq.gcbin <- tmp[,seq(3,601,by=6)]
# gc.bins.quant <- apply(picard.rnaseq.gcbin,2,quantile,c(0.5))
# gc.bins.quant <- gc.bins.quant / sum(gc.bins.quant)
# par(mfrow = c(3,1))
# plot(x = as.numeric(colnames(gc.coverage.quant)),y = gc.coverage.quant[2,],xlab =
#     "% GC content",ylab = "Normalized coverage",pch = 19,ylim = c(0,15),main =
#     "Coverage by GC percentage - median with 95% CIs across samples")
# lines(x = as.numeric(colnames(gc.coverage.quant)),y = gc.coverage.quant[2,],col = "black")
# lines(x = as.numeric(colnames(gc.coverage.quant)),y = gc.coverage.quant[1,],col = "grey")
# lines(x = as.numeric(colnames(gc.coverage.quant)),y = gc.coverage.quant[3,],col = "grey")
# plot(x = as.numeric(colnames(gc.coverage.quant)),y = gc.qual.quant[2,],xlab =
#     "%GC content",ylab = "Read quality score",pch = 19,ylim = c(0,36),main =
#     "Read quality by GC percent - median with 95% CIs across samples")
# lines(x = as.numeric(colnames(gc.coverage.quant)),y = gc.qual.quant[2,],col = "black")
# lines(x = as.numeric(colnames(gc.coverage.quant)),y = gc.qual.quant[1,],col = "grey")
# lines(x = as.numeric(colnames(gc.coverage.quant)),y = gc.qual.quant[3,],col = "grey")
# plot.density(x = as.numeric(colnames(gc.coverage.quant)),y = gc.bins.quant,xlab = "%GC content",ylab =
#     "Proportion of 100bp bins",pch = 19,ylim = c(0,max(gc.bins.quant) + 0.01),main =
#     "Proportion of bins corresponding to each GC percentile")


## GC metrics
gc.summary <- read.table("../data/QC/GCsummary.txt", skip=1)


## Compile most important QC metrics
PQCdat <- cbind(picard.align[,c(2,6,9)],picard.rnaseq.qc[,c(3:7,18:21)],gc.summary[,c(7,8)]) #gc.summary index=7 ~ AT_DROPOUT, index=8  ~ GC_DROPOUT

colnames(PQCdat)[13:14] = c("AT_DROPOUT", "GC_DROPOUT") 
write.csv(PQCdat,"../data/QC/PicardToolsQC.csv")
