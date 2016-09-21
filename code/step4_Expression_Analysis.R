options(stringsAsFactors = F)

rm(list=ls()) #Clear workspace

setwd("/Users/sepid/Documents/Geschwind Lab/TSC_MIA_RNAseq/code")



#Install Necessary Packages
#install.packages()   #Only for package on CRAN repository
#source("http://bioconductor.org/biocLite.R") #Packages on Bioconductor
#biocLite(c("DESeq2", "cqn"))
#biocLite(c("WGCNA", "biomaRt", "ggplot2", "reshape", "limma", "edgeR", "gProfileR", "gplots", "venneuler", "nlme"))

#Load Libraries
library(WGCNA); library(DESeq2); library(biomaRt); library(ggplot2); library(reshape); library(cqn); library(limma); library(edgeR)
library(gProfileR); library(gplots); library(venneuler); library(nlme)
library(limma)
library(gridExtra)
library(gProfileR)


load("../data/HTseqCounts.RData")
RIN_data <- read.csv("../data/Processed Silva_TSC_MIA_MsRNAseq20160510.csv")


#Import Raw Data
#Meta Data
datMeta$Genotype = "Het"
datMeta$Genotype[datMeta$Subject %in% c(442,466,469)] = "WT"
datMeta$Treatment = "PolyIC"
datMeta$Treatment[datMeta$Subject %in% c(420, 455, 447)] = "Saline"
datMeta$Group = as.factor(paste(datMeta$Genotype, "_", datMeta$Treatment,sep=""))

idx = match(datMeta$Sample, gsub("_","-",RIN_data$Sample))
datMeta$RIN = RIN_data$RNA.RIN[idx]
datMeta$RNA_concentration = RIN_data$Conc...ng.ul..from.NanoChip[idx]

#table(datMeta$Group) #tabulates number of subjects in each group

#Expression Data
gene_ens <- rownames(datExpr) 
seq_depth <- apply(datExpr,2,sum)
datSeq <- data.frame(seq_depth)
colnames(datSeq) = c("SeqDepth")
datMeta$SeqDepth = seq_depth

#removes numbers following decimal from ensembl ID
gene_ens_truncated = character(length=length(gene_ens))
for (i in 1:length(gene_ens)){
  s <- gene_ens[i]
  truncated <- substr(s, 1, 18)
  gene_ens_truncated[i] <- truncated
}
rownames(datExpr) = gene_ens_truncated
idx = grep("ENSMUS", rownames(datExpr))
datExpr = datExpr[idx,]


#Sequencing Statistics from Picard
#######datSeq
datSeq2 = read.delim("../data/QC/RNAseqQC.txt", sep = "")
datSeq2$SeqDepth = seq_depth
colnames(datSeq)[1] = "Sample"
datSeq2$Sample = datMeta$Sample



#Annotate Probes
annotate = FALSE

if(annotate){
  bm = useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl")
  marts = listMarts(useMart("ensembl"))
  datasets = listDatasets(useMart("ENSEMBL_MART_ENSEMBL"))

  bm = useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl")

  a = listAttributes(bm); f= listFilters(bm)
  bm1 = getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene", "chromosome_name", "start_position", "end_position", "percentage_gc_content"),
       filters = "ensembl_gene_id",
       values=gene_ens_truncated,mart=bm)
  idx = match(gene_ens_truncated,bm1$ensembl_gene_id)
  datProbes = bm1[idx,]
}

#Load Annotated Probes
load("../data/datProbes.rda")
idx = match(gene_ens_truncated, datProbes$ensembl_gene_id)
datProbes = datProbes[idx,] 


#Filter Genes
to_keep = apply(datExpr>10, 1, sum)
to_keep = to_keep >= 0.5*ncol(datExpr) #At least 10 counts in 50% of samples
keep_ind = which(to_keep)
table(to_keep)
datExpr = datExpr[to_keep,]
datProbes = datProbes[keep_ind,]



#CPM for QC
datExpr.cpm = voom(calcNormFactors(DGEList(datExpr)), data=datMeta)$E


#QC, Normalization, Outlier Removal
#pdf("./figures/Fig1-QC-Prenorm.pdf")
par(mfrow=c(1,1))
boxplot(datExpr.cpm, col=as.numeric(datMeta$Group), ylab ="log2 CPM"); legend("topright", c("Het_PolyIC", "Het_Saline", "WT_PolyIC"), col=c("black","red","green"), pch=19, cex = 0.7)
i = 1; plot(density((datExpr.cpm[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   
for(i in 2:dim(datExpr.cpm)[2]) {     lines(density((datExpr.cpm[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))

mds = cmdscale(dist(t(datExpr.cpm)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col="grey60", pch=16,main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
plot(mds$points, col=as.factor(datMeta$Region), pch=16,main="MDS Plot by Region", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
plot(mds$points, col=as.factor(datMeta$Genotype), pch=16,main="MDS Plot by Genotype", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
plot(mds$points, col=as.factor(datMeta$Treatment), pch=16,main="MDS Plot by Treatment", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))


tree = hclust(dist(t(datExpr.cpm)), method = "average");   
plot(tree)
rin_col = numbers2colors(datMeta$RIN, blueWhiteRed(100), signed=TRUE, centered = TRUE, lim=c(0,10))
seqdepth_col = numbers2colors(datMeta$SeqDepth)
plotDendroAndColors(tree,colors = cbind(as.factor(datMeta$Region), as.factor(datMeta$Genotype), as.factor(datMeta$Treatment), as.factor(datMeta$Hemisphere), seqdepth_col, rin_col), groupLabels = c("Region", "Genotype", "Treatment", "Hemisphere","SeqDepth", "RIN"))


#Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr.cpm, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); 
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(Z.K, col = as.numeric(datMeta$Group), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-2, lty=2)
datExpr.cpm = datExpr.cpm[,!outliers]
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]


#Run DESeq by region
regions = c("all", "cbl", "hc", "pfc")

dds.global = DESeqDataSetFromMatrix(datExpr, datMeta, ~Genotype + Treatment + Region + Hemisphere + RIN) # + seqPC1 + seqPC2
dds.global = estimateSizeFactors(dds.global); dds.global = DESeq(dds.global); 



datExpr.cbl = datExpr[,grep("CB", colnames(datExpr))]
datMeta.cbl = datMeta[which(datMeta$Region == "CBL"),]
dds.cbl = DESeqDataSetFromMatrix(datExpr.cbl, datMeta.cbl, ~Genotype + Treatment + Hemisphere + RIN) # + seqPC1 + seqPC2
dds.cbl = estimateSizeFactors(dds.cbl); dds.cbl = DESeq(dds.cbl);


datExpr.hc = datExpr[,grep("HP", colnames(datExpr))]
datMeta.hc = datMeta[which(datMeta$Region == "HC"),]
dds.hc = DESeqDataSetFromMatrix(datExpr.hc, datMeta.hc, ~Genotype + Treatment + Hemisphere + RIN) # + seqPC1 + seqPC2
dds.hc = estimateSizeFactors(dds.hc); dds.hc = DESeq(dds.hc);


datExpr.pfc = datExpr[,grep("PFC", colnames(datExpr))]
datMeta.pfc = datMeta[which(datMeta$Region == "PFC"),]
dds.pfc = DESeqDataSetFromMatrix(datExpr.pfc, datMeta.pfc, ~Genotype + Treatment + Hemisphere + RIN) # + seqPC1 + seqPC2
dds.pfc = estimateSizeFactors(dds.pfc); dds.pfc = DESeq(dds.pfc);


dds <- list(dds.global, dds.cbl, dds.hc, dds.pfc)


list_contrasts = c("genotype", "treatment")


#initialize data frames to store tallies  
tally_genotype <- data.frame(matrix(0, ncol = 2, nrow = 4))
rownames(tally_genotype) = regions; colnames(tally_genotype) <- c("up", "down")
tally_treatment <- data.frame(matrix(0, ncol = 2, nrow = 4))
rownames(tally_treatment) = regions; colnames(tally_treatment) <- c("up", "down")

plot_volcano <- TRUE

#Find significantly regulated genes and GO in each region
for (reg_ind in 1:length(regions)){
  
  curr_dds <- dds[[reg_ind]]
  curr_reg <- regions[reg_ind]
  
  #Output of DEseq2 for the region query
  res.genotype <- results(dds[[reg_ind]],contrast=c("Genotype", "Het", "WT"))
  res.treatment <- results(dds[[reg_ind]],contrast=c("Treatment", "PolyIC", "Saline"))
  
  res.genotype$hsapiens_homolog=datProbes$hsapiens_homolog_ensembl_gene; res.genotype$gene =datProbes$external_gene_name
  res.treatment$hsapiens_homolog=datProbes$hsapiens_homolog_ensembl_gene; res.treatment$gene =datProbes$external_gene_name
  
  
  res = list(res.genotype, res.treatment)
  
  #make plots, identify signifcantly regulated genes, and GO
  for (i in 1:length(list_contrasts)){
    
    curr_res = res[[i]]
    contrast = list_contrasts[i]
    
    #plots - MA and volcano
    if (plot_volcano){
    
    plot.new()
    
    par(mfrow=c(2,1),mar=c(5,4,2,2))
    DESeq2::plotMA(curr_res, main=paste(curr_reg, contrast), ylim=c(-2,2))
    
    
    c= rgb(t(col2rgb(as.numeric(1))),alpha=100,maxColorValue = 255)
    minp = min(curr_res$padj[curr_res$padj>0])
    curr_res$padj[curr_res$padj==0] = minp
    maxp = 1.1*max(-log10(curr_res$padj));
    
    plot(curr_res$log2FoldChange, -log10(curr_res$padj), xlab="Log2 Fold Change", ylab="log10(P.adj)",pch=19,col=c,cex=.5, xlim = c(-2,2), ylim=c(0,maxp))
    idx = which(curr_res$padj<0.05)
    points(curr_res$log2FoldChange[idx], -log10(curr_res$padj)[idx],col="red",cex=0.6)
    idx = which(curr_res$padj<0.005)
    text(curr_res$log2FoldChange[idx], -log10(curr_res$padj)[idx], labels = datProbes$external_gene_name[idx],cex=.5, pos = 3)
    
    
    #significant genes
    ind_sig_genotype = which(curr_res$padj<.1)
    dge = as.data.frame(curr_res[ind_sig_genotype,c("gene", "log2FoldChange", "pvalue", "padj","hsapiens_homolog")])
    dge[,2:4]=apply(dge[,2:4],2,signif,2)
    dge.up = dge[dge$log2FoldChange>0,]; 
    dge.up = dge.up[order(dge.up$log2FoldChange, decreasing = T),]
    dge.down = dge[dge$log2FoldChange<0,]; 
    dge.down = dge.down[order(dge.down$log2FoldChange),]
    
    
    
    #downregulated GO
    query = rownames(dge.down)
    
    go.mus = gprofiler(query, organism="mmusculus", custom_bg = datProbes$ensembl_gene_id, 
                       correction_method = "fdr",hier_filtering = "moderate", ordered_query = T, significant = T, exclude_iea = F,
                       region_query = F, max_p_value = 0.05, max_set_size=1000, numeric_ns = "",
                       include_graph = F,src_filter = c("GO", "KEGG", "REACTOME"))
    go = go.mus[order(go.mus$p.value),]
    
    
    par(oma=c(0,15,0,0))
    ttl = paste("GO downregulated", contrast, curr_reg)
    n_go_show = min(10, dim(dge.down)[1])
    bp = barplot(-log10(as.numeric(na.omit(go$p.value[n_go_show:1]))), main = ttl, horiz=T, yaxt='n', col="blue", xlab='-log10(p)',cex.main=0.7, cex.axis = .7)
    axis(2,at=bp,labels=na.omit(go$term.name[n_go_show:1]),tick=FALSE,las=2,cex.axis=.7);
    abline(v=-log10(0.05), col="red", lwd=2,lty=2)
    
    
    #upregulated GO
    query = rownames(dge.up)[order(dge.up$log2FoldChange, decreasing = T)]

    go.mus = gprofiler(query, organism="mmusculus", custom_bg = datProbes$ensembl_gene_id,
                       correction_method = "fdr",hier_filtering = "moderate", ordered_query = T, significant = T, exclude_iea = F,
                       region_query = F, max_p_value = 0.05, max_set_size=1000, numeric_ns = "",
                       include_graph = F,src_filter = c("GO", "KEGG", "REACTOME"))
    
    go = go.mus[order(go.mus$p.value),]
    
    
    if (nrow(go) > 0){
    
    
    #tally up total upregulated and downregulated
    n_up = dim(dge.up)[1]
    n_down = dim(dge.down)[1]
    
   if (contrast == "genotype") {
      tally_genotype[reg_ind,] <- c(n_up, n_down)
      tally_treatment[reg_ind,] <- c(n_up, n_down)
   }
  }
}


dev.off()



##Network Analysis
datMeta = datMeta
datExpr.vst = varianceStabilizingTransformation(dds.global)
datExpr.vst = assay(datExpr.vst)


multiExpr = vector(mode="list", length=1)
multiExpr[[1]] = list(data=as.data.frame(t(datExpr.vst)), meta=datMeta)
<<<<<<< HEAD
#multiExpr[[2]] = list(data=as.data.frame(t(dE[,!datMeta$Region=="HC"])), meta = datMeta.cortex)
#multiExpr[[3]] = list(data=as.data.frame(t(dE[,datMeta$Region=="HC"])), meta= datMeta.hc)
#multiExpr[[4]] = list(data=as.data.frame(t(dE[,datMeta$Region=="DLPFC"])), meta=datMeta.pfc)
#names(multiExpr) = c("ALL", "Cortical", "HC", "DLPFC")


##Step 1 Choose Soft Threshold Power 

if(FALSE)
{
  bsize=6000
  powers=seq(2,30,by=2)
  for(n in 1:1) { 
    multiExpr[[n]]$softThresh = pickSoftThreshold(data= multiExpr[[n]]$data, networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers,blockSize = bsize)
    
    sft = multiExpr[[n]]$softThresh
    par(mfrow=c(1,2))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Thresh Power", ylab="Scale free R^2",type="n", main=names(multiExpr)[n])
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2",main=names(multiExpr)[n])
    abline(h=0.8, col="black")
    plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
  }
}



wgcna_parameters = list(powers =  c(18))
wgcna_parameters$minModSize = 100
wgcna_parameters$minHeight = 0.1
wgcna_parameters$bsize = 18000
wgcna_parameters$ds = 2
wgcna_parameters$networkType = "signed"
wgcna_parameters$corFnc = "bicor"
wgcna_parameters$pamStage = TRUE

if(TRUE) {
  for (n in 1:4) {
    ##Calculate TOM, save to file
    rootdir= getwd()
    filenm <- paste(rootdir, "/processed_data/WGCNA/network_signed_exprSet_cqn.noregress_", as.character(n),sep="")
    multiExpr[[n]]$netData = blockwiseModules(datExpr=multiExpr[[n]]$data, maxBlockSize=wgcna_parameters$bsize, networkType=wgcna_parameters$networkType, corType = wgcna_parameters$corFnc ,  power = wgcna_parameters$powers[n], mergeCutHeight= wgcna_parameters$minHeight, nThreads=23, 
                                              saveTOMFileBase=filenm, saveTOMs=TRUE, minModuleSize= wgcna_parameters$minModSize, pamStage=wgcna_parameters$pamStage, reassignThreshold=1e-6, verbose = 3, deepSplit=wgcna_parameters$ds)
  }
}


set.idx = 1
load("./processed_data/WGCNA/network_signed_exprSet_cqn_3-block.1.RData")

datMeta = multiExpr[[set.idx]]$meta
geneTree = hclust(as.dist(1-TOM), method = "average"); 


colors = vector(mode="list"); labels = vector(mode="list"); labels=""
pam=F; minModSize=100; ds=2; dthresh=0.1
tree = cutreeHybrid(dendro = geneTree, minClusterSize= minModSize, pamStage=pam, cutHeight = 0.999, deepSplit=ds, distM=as.matrix(1-TOM))
merged = mergeCloseModules(exprData= multiExpr[[1]]$data, colors = tree$labels, cutHeight=dthresh)
colors = labels2colors(merged$colors)
dev.off(); plotDendroAndColors(geneTree, colors, groupLabels=labels, addGuide= TRUE, dendroLabels=FALSE, main="Dendrogram", cex.colorLabels=0.5)
table(colors)

MEs = moduleEigengenes(expr = (multiExpr[[set.idx]]$data), colors = labels2colors(merged$colors), softPower = wgcna_parameters$powers[set.idx])
kME = signedKME(multiExpr[[set.idx]]$data, MEs$eigengenes,corFnc = "bicor")
