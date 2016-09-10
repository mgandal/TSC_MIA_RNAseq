
options(stringsAsFactors = F)
rm(list=ls()) #Clear workspace

setwd("/Users/mgandal/Documents/Github/TSC_MIA_RNAseq")

#Install Necessary Packages
#-------------------------
#install.packages()   #Only for package on CRAN repository
#source("http://bioconductor.org/biocLite.R") #Packages on Bioconductor
#biocLite(c("DESeq2", "CQN"))

#Load Libraries
library(DESeq2); library(cqn)



#Import Raw Data
load("./data/HTseqCounts.RData")
#datExpr

datMeta$Genotype = "Het"
datMeta$Genotype[datMeta$Subject %in% c(442,466,469)] = "WT"
datMeta$Treatment = "PolyIC"
datMeta$Treatment[datMeta$Subject %in% c(420, 455, 447)] = "Saline"
datMeta$Group = as.factor(paste(datMeta$Genotype, "_", datMeta$Treatment,sep=""))

table(datMeta$Group)









###ADAPT MONKEY CODE FROM BELOW


## NHP_MIA-Step1-ImportData.R
rm(list=ls())
options(stringsAsFactors = T)
setwd("~/Dropbox/GeschwindLab/Projects/CONTE-SCZ-MIA/")

#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(DESeq2); library(biomaRt); library(ggplot2); library(reshape); library(cqn); library(limma); library(edgeR)
library(gProfileR); library(gplots); library(venneuler); library(nlme)
options(stringsAsFactors = F)


##----Load Expression Data
datExpr = read.csv("./raw_data/SeqBatch1+2/datExpr.htseq.20160619.csv")
rownames(datExpr) = datExpr$X; datExpr = datExpr[,-1]
colnames(datExpr) = gsub("X", "", colnames(datExpr))
datSeq = as.data.frame(t(datExpr[(nrow(datExpr)-4):(nrow(datExpr)-3),]))
datSeq$Depth = apply(datExpr,2,sum)
datExpr = datExpr[1:30246, ]


##----Load MetaData
datMeta= read.csv("./raw_data/SeqBatch1+2/datMeta_20160601.csv")
rownames(datMeta) = datMeta$Sample
idx = match(colnames(datExpr), rownames(datMeta))
datMeta = datMeta[idx,]
datMeta$MIA = "MIA"; datMeta$MIA[datMeta$Group=="con"] = "CON"
datMeta$MIA = factor(datMeta$MIA, levels=c("CON", "MIA"))
datMeta$Group = factor(datMeta$Group, levels=c("con","poly1", "poly2"))
datMeta$SeqBatch = as.factor(datMeta$SeqBatch)
datMeta$Region = as.factor(datMeta$Region)

##-----Load Sequencing Statistics from Picard
datSeq2 = read.csv("./raw_data/SeqBatch1+2/qc_picard/PicardToolsQC.csv")
rownames(datSeq2) = gsub("_cat.coordSorted.bam", "", datSeq2$X);
datSeq2=datSeq2[,-1]
idx = match(rownames(datSeq),rownames(datSeq2)); datSeq = cbind(datSeq,datSeq2[idx,])
rm(datSeq2)
datSeq[,c(1:2, 5:14)] = (datSeq[,c(1:2, 5:14)] )
PC.datSeq = prcomp(t(scale(datSeq,scale=T)),center=F)
varexp <- (PC.datSeq$sdev)^2 / sum(PC.datSeq$sdev^2)
corrplot::corrplot(cor(cbind(PC.datSeq$rotation[,1:5], datSeq)))
datMeta$seqPC1=PC.datSeq$rotation[,1]
datMeta$seqPC2=PC.datSeq$rotation[,2]
datMeta$seqPC3=PC.datSeq$rotation[,3]
datMeta$seqPC4=PC.datSeq$rotation[,4]
datMeta$SeqDepth = apply(datExpr,2,sum)

## Annotate Probes
if(FALSE) {
  bm = useMart("ensembl", "mmulatta_gene_ensembl")
  a = listAttributes(bm); f= listFilters(bm)
  bm1 = getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene", "chromosome_name", "start_position", "end_position", "percentage_gc_content"),
              filters = "ensembl_gene_id", values=rownames(datExpr.htsc),mart=bm)
  
  idx = match(rownames(datExpr.htsc),bm1$ensembl_gene_id)
  datProbes = bm1[idx,]
  save(datProbes, file="./raw_data/SeqBatch1+2/datProbes_20160411.rda")
}
load("./raw_data/SeqBatch1+2/datProbes_20160411.rda")
idx = match(rownames(datExpr), datProbes$ensembl_gene_id)
datProbes = datProbes[idx,]

##--Filter Genes
to_keep = apply(datExpr>10, 1, sum)
to_keep = to_keep >= 0.5*ncol(datExpr) #At least 10 counts in 50% of samples
table(to_keep)
datExpr = datExpr[to_keep,]
datProbes = datProbes[to_keep,]


##---CPM for QC
datExpr.cpm = voom(calcNormFactors(DGEList(datExpr)), data=datMeta)$E

##----------------QC, Normalization, Outlier Removal ----------------
#pdf("./figures/Fig1-QC-Prenorm.pdf")
par(mfrow=c(1,1))
boxplot(datExpr.cpm, col=as.numeric(datMeta$MIA), ylab ="log2 CPM"); legend("topright", c("CTL", "MIA"), col=c("black","red"), pch=19, cex = 0.7)
i = 1; plot(density((datExpr.cpm[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   
for(i in 2:dim(datExpr.cpm)[2]) {     lines(density((datExpr.cpm[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))

mds = cmdscale(dist(t(datExpr.cpm)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points[datMeta$MIA=="MIA",1], mds$points[datMeta$MIA=="MIA",2], col="grey60", pch=16,main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
for(i in 1:30) points(mds$points, col=as.numeric(as.factor(datMeta$Region)),cex=1.1);   
legend("topleft", levels(datMeta$Region), col=c(1:length(levels(datMeta$Region))), pch=16, cex=0.8)
legend("bottomleft", levels(datMeta$MIA), col=c("white", "grey60"), pch=16, cex=0.8)


tree = hclust(dist(t(datExpr.cpm)), method = "average");   
plot(tree)
rin_col = numbers2colors(datMeta$X260.280, blueWhiteRed(100), signed=FALSE, centered = FALSE, lim=c(min(datMeta$X260.280),max(datMeta$X260.280)))
seqdepth_col = numbers2colors(datMeta$SeqDepth)
plotDendroAndColors(tree,colors = cbind(as.numeric(datMeta$MIA), as.numeric(datMeta$Region), as.numeric(as.factor(datMeta$RNAisoBatch)), as.numeric(datMeta$SeqBatch),rin_col, seqdepth_col), groupLabels = c("Group", "Region", "RNAisoBatch", "SeqBatch", "RNA:260/280", "Seq Depth"))


##Outlier Removal
sdout <- 3; normadj <- (0.5+0.5*bicor(datExpr.cpm, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); 
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(Z.K, col = as.numeric(datMeta$Group), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-3, lty=2)
datExpr.cpm = datExpr.cpm[,!outliers]
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]


##-----------------DESeq2 -- Full Model
dds = DESeqDataSetFromMatrix(datExpr, datMeta, ~MIA + Region + SeqBatch + X260.230 + seqPC1 + seqPC2)
dds = estimateSizeFactors(dds); dds = DESeq(dds); 
res.all.mia = results(dds,contrast=c("MIA", "MIA", "CON"))
table(res.all.mia$padj<.1)
res.all.mia$hsapiens_homolog=datProbes$hsapiens_homolog_ensembl_gene; res.all.mia$gene =datProbes$external_gene_name

r = rlog(dds,blind = F); datExpr.rlog = assay(r)
datExpr.vst = assay(varianceStabilizingTransformation(dds,blind=F))
datMeta$sizeFactor = colData(dds)$sizeFactor

resOrdered = res.all.mia[order(res.all.mia$padj),]
head(resOrdered)
table(res.all.mia$padj<.1)


dds = DESeqDataSetFromMatrix(datExpr, datMeta, ~Group + Region + X260.230 + seqPC1 + seqPC2)
dds = estimateSizeFactors(dds); dds = DESeq(dds); 
res.all.poly1 = results(dds,contrast=c("Group", "poly1", "con"))
res.all.poly2 = results(dds,contrast=c("Group", "poly2", "con"))



#pdf("./figures/Fig4-Volcano.pdf")
par(mfrow=c(2,1),mar=c(4,4,2,2))
DESeq2::plotMA(res.all.mia,ylim=c(-1,1))
c= rgb(t(col2rgb(as.numeric(1))),alpha=100,maxColorValue = 255)
plot(res.all.mia$log2FoldChange, -log10(res.all.mia$padj),xlab="Log2 Fold Change", ylab="log10(P.adj)",pch=19,col=c,cex=.5, ylim=c(0,20))
idx = which(res.all.mia$padj<0.05)
points(res.all.mia$log2FoldChange[idx], -log10(res.all.mia$padj)[idx],col="red",cex=0.6)
idx = which(res.all.mia$padj<0.005)
text(res.all.mia$log2FoldChange[idx], -log10(res.all.mia$padj)[idx], labels = datProbes$external_gene_name[idx],cex=.5, pos = 3)



dge = as.data.frame(res.all.mia[res.all.mia$padj<.1,c("gene", "log2FoldChange", "pvalue", "padj","hsapiens_homolog")])
dge[,2:4]=apply(dge[,2:4],2,signif,2)
dge.up = dge[dge$log2FoldChange>0,]; 
dge.up = dge.up[order(dge.up$log2FoldChange, decreasing = T),]
dge.down = dge[dge$log2FoldChange<0,]; 
dge.down = dge.down[order(dge.down$log2FoldChange),]

plot.new();grid.table(dge.up[1:33,1:4])
plot.new();grid.table(dge.down[1:30,1:4])

query = rownames(dge.down)[order(dge.down$log2FoldChange, decreasing = F)]
go.mmu = gprofiler(query, organism="mmulatta", custom_bg = datProbes$ensembl_gene_id, 
                   correction_method = "fdr",hier_filtering = "moderate", ordered_query = T, significant = T, exclude_iea = F,
                   region_query = F,max_p_value = 0.05, max_set_size=1000, numeric_ns = "",
                   include_graph = F,src_filter = c("GO", "KEGG", "REACTOME"))

go = go.mmu[order(go.mmu$p.value),]
ttl = "HC - Upregulate"
par(oma=c(0,15,0,0));
bp = barplot(-log10(as.numeric(na.omit(go$p.value[10:1]))), main=ttl, horiz=T, yaxt='n', col="blue", xlab='-log10(p)',cex.main=0.7, cex.axis = .7)
axis(2,at=bp,labels=na.omit(go$term.name[10:1]),tick=FALSE,las=2,cex.axis=.7);
abline(v=-log10(0.05), col="red", lwd=2,lty=2)



go.hs = gprofiler(dge.down$hsapiens_homolog, custom_bg = datProbes$hsapiens_homolog_ensembl_gene,
                  correction_method = "fdr",hier_filtering = "none", ordered_query = T, significant = T, exclude_iea = F,
                  region_query = F,max_p_value = 0.05, max_set_size=1000, numeric_ns = "",
                  include_graph = F,src_filter = c("GO", "KEGG", "REACTOME"))
go = go.hs[order(go.hs$p.value),]
