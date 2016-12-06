
#--------------Initialize workspace-------------
options(stringsAsFactors = F)
rm(list=ls()) #Clear workspace


#-----------------Install necessary packages---------------
#install.packages()   #Only for package on CRAN repository
#source("http://bioconductor.org/biocLite.R") #Packages on Bioconductor
#biocLite(c("DESeq2", "cqn"))
#biocLite(c("WGCNA", "biomaRt", "ggplot2", "reshape", "limma", "edgeR", "gProfileR", "gplots", "venneuler", "nlme"))


#-------------Load necessary libraries--------------
library(WGCNA); library(DESeq2); library(biomaRt); library(ggplot2); library(reshape); library(cqn); library(limma); library(edgeR)
library(gProfileR); library(gplots); library(venneuler); library(nlme)
library(limma)
library(gridExtra)
library(gProfileR)


#-----------Designate PDF file to save figures---------
pdf(file="../figures/TSC-graphs-Het.pdf")


#--------------Load expression counts, meta data, Picard stats, and RIN data------------
load("../data/HTseqCounts.RData")
datSeq <- read.csv("../data/QC/PicardToolsQC.csv")
RIN_data <- read.csv("../data/Processed Silva_TSC_MIA_MsRNAseq20160510.csv")


#-------------Include only Het subjects, exclude WT------------
datMeta$Genotype = "Het"
datMeta$Genotype[datMeta$Subject %in% c(442,466,469)] = "WT"
datMeta$Genotype <- as.factor(datMeta$Genotype)
#datMeta$Genotype <- relevel(datMeta$Genotype, "WT")
Het_inds = which(datMeta$Genotype == "Het")
datMeta = datMeta[Het_inds,]
rownames(datMeta) = c(1:dim(datMeta)[1])

datExpr = datExpr[,Het_inds]

datSeq = datSeq[Het_inds,]
rownames(datSeq) = c(1:dim(datSeq)[1])



#------------Organize remaining meta data--------------
datMeta$Treatment = "PolyIC"
datMeta$Treatment[datMeta$Subject %in% c(420, 455, 447)] = "Saline"
datMeta$Treatment <- as.factor(datMeta$Treatment)
datMeta$Treatment <- relevel(datMeta$Treatment, "Saline")

datMeta$Group = as.factor(paste(datMeta$Genotype, "_", datMeta$Treatment,sep=""))

datMeta$Region <- as.factor(datMeta$Region)

idx = match(datMeta$Sample, gsub("_","-",RIN_data$Sample))
datMeta$RIN = RIN_data$RNA.RIN[idx]
datMeta$RNA_concentration = RIN_data$Conc...ng.ul..from.NanoChip[idx]



#------------------------Organize expression data----------------------------
#--------------Remove numbers following decimal from ensembl ID--------------
gene_ens <- rownames(datExpr)
gene_ens_truncated = character(length=length(gene_ens))
for (i in 1:length(gene_ens)){
  s <- gene_ens[i]
  truncated <- substr(s, 1, 18)
  gene_ens_truncated[i] <- truncated
}
rownames(datExpr) = gene_ens_truncated


#---------------Look for rows that actually represent genes----------------
idx = grep("ENSMUS", rownames(datExpr))
datExpr = datExpr[idx,]


#------------Organize Picard sequencing stats------------
seq_depth <- apply(datExpr,2,sum)
datSeq$SeqDepth = seq_depth
datMeta$SeqDepth = seq_depth
datSeq$Sample = datMeta$Sample

datSeq_unlab = datSeq[,2:(length(datSeq)-1)] #removes sample ID column to allow PCA


#----------Compute principal components of Picard QC data---------
PC.datSeq = prcomp(t(scale(datSeq_unlab,scale=T)),center=F)
varexp <- (PC.datSeq$sdev)^2 / sum(PC.datSeq$sdev^2)
corrplot::corrplot(cor(cbind(PC.datSeq$rotation[,1:5], datSeq_unlab)))
datMeta$seqPC1=PC.datSeq$rotation[,1]
datMeta$seqPC2=PC.datSeq$rotation[,2]
datMeta$seqPC3=PC.datSeq$rotation[,3]
datMeta$seqPC4=PC.datSeq$rotation[,4]


#--------------Annotate probes if not done already-----------------
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
  save(datProbes, file="../data/datProbes.rda")
}



#--------------------Load annotated probes-------------------
load("../data/datProbes.rda")
idx = match(gene_ens_truncated, datProbes$ensembl_gene_id)
datProbes = datProbes[idx,] 


#------------Filter genes to include only those with adequate expression--------
to_keep = apply(datExpr>10, 1, sum) #compute n samples with at least 10 counts
to_keep = to_keep >= 0.5*ncol(datExpr) #>=50% of samples have >= 10 counts of the gene
keep_ind = which(to_keep)
datExpr = datExpr[keep_ind,]
datProbes = datProbes[keep_ind,]



#-------------CPM (counts per million) for QC--------------
datExpr.cpm = voom(calcNormFactors(DGEList(datExpr)), data=datMeta)$E



#-------------------QC visualization---------------
#------------------log2CPM boxplot-----------------
par(mfrow=c(1,1))
boxplot(datExpr.cpm, col=as.numeric(datMeta$Group), ylab ="log2 CPM"); legend("topright", c("Het_PolyIC", "Het_Saline", "WT_PolyIC"), col=c("black","red","green"), pch=19, cex = 0.7)

#------------------log2CPM density-----------------
i = 1; plot(density((datExpr.cpm[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   
for(i in 2:dim(datExpr.cpm)[2]) {     lines(density((datExpr.cpm[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))

#--------------MDS (multidimensional scaling) plots---------
mds = cmdscale(dist(t(datExpr.cpm)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col="grey60", pch=16,main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""))
plot(mds$points, col=as.factor(datMeta$Region), pch=16,main="MDS Plot by Region", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""), cex.lab=1.3) ; legend("topleft", levels(datMeta$Region),cex=1.5, text.col=1:length(levels(datMeta$Region)))
#plot(mds$points, col=as.factor(datMeta$Genotype), pch=16,main="MDS Plot by Genotype", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""))
plot(mds$points, col=as.factor(datMeta$Treatment), pch=16,main="MDS Plot by Treatment", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); legend("topleft", levels(datMeta$Treatment),cex=1, text.col=1:length(levels(datMeta$Treatment)))


#------------------Dendrogram and trait correlations-----------
tree = hclust(dist(t(datExpr.cpm)), method = "average");   
rin_col = numbers2colors(datMeta$RIN, colors=blueWhiteRed(100), signed=TRUE, centered = TRUE)
seqdepth_col = numbers2colors(datMeta$SeqDepth, colors=blueWhiteRed(100), signed=TRUE, centered = TRUE)
seqPC1_col = numbers2colors(datMeta$seqPC1, colors=blueWhiteRed(100), signed=TRUE, centered = TRUE)
seqPC2_col = numbers2colors(datMeta$seqPC2, colors=blueWhiteRed(100), signed=TRUE, centered = TRUE)
plotDendroAndColors(tree,colors = cbind(as.factor(datMeta$Region), as.factor(datMeta$Treatment), as.factor(datMeta$Hemisphere), seqdepth_col, rin_col, seqPC1_col, seqPC2_col), groupLabels = c("Region", "Treatment", "Hemisphere","SeqDepth", "RIN", "SeqPC1", "SeqPC2"), cex.lab = 1.3, cex.colorLabels = 1.1)



#----------------Outlier removal---------------------
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




#-----------EXPRESSION ANALYSIS-----------
#-------Choose which kind of differential expression analysis will be used--------
DE_method = "DESeq"
DE_method = "limma"
DE_method = "edgeR"



#-----------Run differential analysis by region using limma VOOM----------
if (DE_method == "limma"){
  datExprvoom = voom(calcNormFactors(DGEList(datExpr)), data=datMeta)
  design = model.matrix(~Treatment + Region + Hemisphere + RIN, data=datMeta)
  fit <- lmFit(datExprvoom, design)
  fit <- eBayes(fit, trend=TRUE)
  
  n_genes = dim(datExpr)[1]
  limma_all = topTable(fit, coef=2,genelist = datProbes$ensembl_gene_id, number = n_genes)
  limma_all = limma_all[rownames(datExpr),]
  FC_limma_all = limma_all$logFC
  
  limma_beta = as.data.frame(fit$coefficients)$TreatmentPolyIC
  
  reg_genes_limma=topTable(fit, coef=2,genelist = datProbes$external_gene_name,number = 15)
  reg_genes_limma=reg_genes_limma[reg_genes_limma$adj.P.Val<0.1,]
}


#--------Run differential analysis by region using edgeR--------
#let saline = 1 and polyIC = 2
group <- rep(2,dim(datExpr)[2]) #subjects 420, 455, 447 treated w/ saline
group[grep("420|455|447", colnames(datExpr))] = 1
as.factor(group)
y <- DGEList(counts=datExpr,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
et <- exactTest(y)

glm <- glmFit(y, design)
edgeR_genes = rownames(glm$coefficients)
edgeR_beta = as.data.frame(glm$coefficients)$group

FC_edgeR_all = as.data.frame(et$table)$logFC

top_edgeR = topTags(et)
top_edgeR = as.data.frame(top_edgeR)
top_edgeR$Symbol = datProbes$external_gene_name[match(rownames(top_edgeR),datProbes$ensembl_gene_id)]
top_edgeR = top_edgeR[,c(5,1:4)] #put gene symbol as first column
top_edgeR = top_edgeR[top_edgeR$FDR<0.1,]




#-----------Run differential analysis by region using DESeq--------
if (DE_method == "DESeq"){
regions = c("all", "cbl", "hc", "pfc")

dds.global = DESeqDataSetFromMatrix(datExpr, datMeta, ~Treatment + Region + Hemisphere + RIN + seqPC1 + seqPC2)
dds.global = estimateSizeFactors(dds.global); dds.global = DESeq(dds.global); 
dds.global$Treatment <- relevel(dds.global$Treatment, ref="Saline")
dds.cbl$Genotype <- relevel(dds.cbl$Genotype, ref="WT")

datExpr.cbl = datExpr[,grep("CB", colnames(datExpr))]
datMeta.cbl = datMeta[which(datMeta$Region == "CBL"),]
dds.cbl = DESeqDataSetFromMatrix(datExpr.cbl, datMeta.cbl, ~Treatment + Hemisphere + RIN + seqPC1 + seqPC2)
dds.cbl = estimateSizeFactors(dds.cbl); dds.cbl = DESeq(dds.cbl);
dds.cbl$Treatment <- relevel(dds.cbl$Treatment, ref="Saline")
#dds.cbl$Genotype <- relevel(dds.cbl$Genotype, ref="WT")

datExpr.hc = datExpr[,grep("HP", colnames(datExpr))]
datMeta.hc = datMeta[which(datMeta$Region == "HC"),]
dds.hc = DESeqDataSetFromMatrix(datExpr.hc, datMeta.hc, ~Treatment + Hemisphere + RIN + seqPC1 + seqPC2)
dds.hc = estimateSizeFactors(dds.hc); dds.hc = DESeq(dds.hc);
dds.hc$Treatment <- relevel(dds.hc$Treatment, ref="Saline")
#dds.hc$Genotype <- relevel(dds.hc$Genotype, ref="WT")


datExpr.pfc = datExpr[,grep("PFC", colnames(datExpr))]
datMeta.pfc = datMeta[which(datMeta$Region == "PFC"),]
dds.pfc = DESeqDataSetFromMatrix(datExpr.pfc, datMeta.pfc, ~Treatment + Hemisphere + RIN + seqPC1 + seqPC2)
dds.pfc = estimateSizeFactors(dds.pfc); dds.pfc = DESeq(dds.pfc);
dds.pfc$Treatment <- relevel(dds.pfc$Treatment, ref="Saline")
#dds.pfc$Genotype <- relevel(dds.pfc$Genotype, ref="WT")

dds <- list(dds.global, dds.cbl, dds.hc, dds.pfc)




#--------Choose what contrasts to use on DEseq data
#list_contrasts = c("genotype", "treatment")
list_contrasts = c("treatment")



#------------Initialize data frames to store tallies-----------  
#tally_genotype <- data.frame(matrix(0, ncol = 2, nrow = length(regions)))
#rownames(tally_genotype) = regions; colnames(tally_genotype) <- c("up", "down")
tally_treatment <- data.frame(matrix(0, ncol = 2, nrow = length(regions)))
rownames(tally_treatment) = regions; colnames(tally_treatment) <- c("up", "down")


#----Indicate if volcano & MA plots should be produced and if GO should be evaluated---
plot_volcano <- TRUE
eval_GO <- FALSE


#----Find significantly regulated genes and perform GO (gene ontology) in each region----
for (reg_ind in 1:length(regions)){
  
  curr_dds <- dds[[reg_ind]]
  curr_reg <- regions[reg_ind]
  
  #-----------Output of DEseq2 for the region query-----------
  #res.genotype <- results(dds[[reg_ind]],contrast=c("Genotype", "Het", "WT"))
  res.treatment <- results(dds[[reg_ind]],contrast=c("Treatment", "PolyIC", "Saline"))
  
  #res.genotype$hsapiens_homolog=datProbes$hsapiens_homolog_ensembl_gene; res.genotype$gene=datProbes$external_gene_name
  res.treatment$hsapiens_homolog=datProbes$hsapiens_homolog_ensembl_gene; res.treatment$gene=datProbes$external_gene_name
  
  #write.csv(res.genotype, file = paste("../data/regulated genes lists QC/res.genotype", curr_reg, ".csv"))
  #write.csv(res.treatment, file = paste("../data/regulated genes lists QC/res.treatment", curr_reg, ".csv"))
  
  #res = list(res.genotype, res.treatment)
  res = list(res.treatment)
  
  
  #-----for each contrast: fold change plots, ID sig regulated genes, perform/plot GO----
  for (i in 1:length(list_contrasts)){
    
    curr_res = res[[i]]
    contrast = list_contrasts[i]
    
    #----------Make MA and volcano plots-------
    if (plot_volcano){
      
      plot.new()
      
      par(mfrow=c(2,1),mar=c(5,4,2,2))
      DESeq2::plotMA(curr_res, main=paste(curr_reg, contrast), ylim=c(-2,2))
      
      c= rgb(t(col2rgb(as.numeric(1))),alpha=100,maxColorValue = 255)
      minp = min(curr_res$padj[curr_res$padj>0], na.rm = TRUE) #ignore indices where padj == 0 when finding min
      curr_res$padj[curr_res$padj==0] = minp
      maxp = 1.1*max(-log10(curr_res$padj), na.rm = TRUE);
      
      plot(curr_res$log2FoldChange, -log10(curr_res$padj), xlab="Log2 Fold Change", ylab="log10(P.adj)",pch=19,col=c,cex=.5, xlim = c(-0.3,0.3), ylim=c(0,4))
      idx = which(curr_res$padj<0.1)
      points(curr_res$log2FoldChange[idx], -log10(curr_res$padj)[idx],col="red",cex=0.6)
      idx = which(curr_res$padj<0.1 & curr_res$log2 < 0)
      if(length(idx) > 0){
      text(curr_res$log2FoldChange[idx], -log10(curr_res$padj)[idx], labels = datProbes$external_gene_name[idx],cex=.5, pos = 2)
      }
      idx = which(curr_res$padj<0.1 & curr_res$log2 > 0)
      if(length(idx) > 0){
        text(curr_res$log2FoldChange[idx], -log10(curr_res$padj)[idx], labels = datProbes$external_gene_name[idx],cex=.5, pos = 4)
      }
      
    }
    
    
    #--------Find which genes significantly up/downregulated-----
    DESeq_beta = as.data.frame(coef(dds.global))$TreatmentPolyIC
    FC_DESeq_all = curr_res$log2FoldChange
    
    ind_sig_genotype = which(curr_res$padj<.1)
    dge = as.data.frame(curr_res[ind_sig_genotype,c("gene", "log2FoldChange", "pvalue", "padj","hsapiens_homolog")])
    dge[,2:4]=apply(dge[,2:4],2,signif,2)
    dge.up = dge[dge$log2FoldChange>0,]; 
    dge.up = dge.up[order(dge.up$log2FoldChange, decreasing = T),]
    dge.down = dge[dge$log2FoldChange<0,]; 
    dge.down = dge.down[order(dge.down$log2FoldChange),]
    
    write.csv(dge.up, file = paste("../data/regulated genes lists QC/dge.up-Het", curr_reg, contrast, ".csv"))
    write.csv(dge.down, file = paste("../data/regulated genes lists QC/dge.down-Het", curr_reg, contrast, ".csv"))
    
    
    #------GO for downregulated genes-----
    #----Access gprofiler for GO if not already computed------
    if(eval_GO){
      query = rownames(dge.down)
      
      go.mus = gprofiler(query, organism="mmusculus", custom_bg = datProbes$ensembl_gene_id, 
                         correction_method = "fdr",hier_filtering = "moderate", ordered_query = T, significant = T, exclude_iea = F,
                         region_query = F, max_p_value = 0.05, max_set_size=1000, numeric_ns = "",
                         include_graph = F,src_filter = c("GO", "KEGG", "REACTOME"))
      go = go.mus[order(go.mus$p.value),]
      
      saveRDS(go, file = paste("../data/regulated genes lists QC/go-down-Het", curr_reg, contrast, ".RDS", sep="-"))
    }
  
    #----------Read in saved GO---------
    go <- readRDS(paste("../data/regulated genes lists QC/go-down-Het", curr_reg, contrast, ".RDS", sep="-"))
  
    
    #-------Plot gene ontologies if any were returned------
    if (nrow(go) > 0){
      par(oma=c(0,15,0,0))
      ttl = paste("GO downregulated", contrast, curr_reg)
      n_go_show = min(10, dim(dge.down)[1])
      bp = barplot(-log10(as.numeric(na.omit(go$p.value[n_go_show:1]))), main = ttl, horiz=T, yaxt='n', col="blue", xlab='-log10(p)',cex.main=0.7, cex.axis = .7)
      axis(2,at=bp,labels=na.omit(go$term.name[n_go_show:1]),tick=FALSE,las=2,cex.axis=.7);
      abline(v=-log10(0.05), col="red", lwd=2,lty=2)
    }
    
    
    #-------GO for upregulated genes-----
    #----Access gprofiler for GO if not already computed------
    if (eval_GO){
    query = rownames(dge.up)[order(dge.up$log2FoldChange, decreasing = T)]
    
    go.mus = gprofiler(query, organism="mmusculus", custom_bg = datProbes$ensembl_gene_id,
                       correction_method = "fdr",hier_filtering = "moderate", ordered_query = T, significant = T, exclude_iea = F,
                       region_query = F, max_p_value = 0.05, max_set_size=1000, numeric_ns = "",
                       include_graph = F,src_filter = c("GO", "KEGG", "REACTOME"))
    
    go = go.mus[order(go.mus$p.value),]

    saveRDS(go, file = paste("../data/regulated genes lists QC/go-up-Het", curr_reg, contrast, ".RDS",sep="-"))
    }
    
    #------Read in saved GO-------
    go <- readRDS(paste("../data/regulated genes lists QC/go-up-Het", curr_reg, contrast, ".RDS",sep="-"))
    
    #-------Plot gene ontologies if any were returned------
    if (nrow(go) > 0){
      par(oma=c(0,15,0,0));
      ttl = paste("GO upregulated", contrast, curr_reg)
      n_go_show = min(10, dim(dge.up)[1])
      bp = barplot(-log10(as.numeric(na.omit(go$p.value[n_go_show:1]))), main=ttl, horiz=T, yaxt='n', col="blue", xlab='-log10(p)',cex.main=0.7, cex.axis = .7)
      axis(2,at=bp,labels=na.omit(go$term.name[n_go_show:1]),tick=FALSE,las=2,cex.axis=.7);
      abline(v=-log10(0.05), col="red", lwd=2,lty=2)
    }
    
    
    #tally up total upregulated and downregulated
    n_up = dim(dge.up)[1]
    n_down = dim(dge.down)[1]
    
    if (contrast == "genotype") {
      tally_genotype[reg_ind,] <- c(n_up, n_down)
    }else if (contrast == "treatment"){
      tally_treatment[reg_ind,] <- c(n_up, n_down)
    }
  }
}

dev.off()
}



#-------------make expression boxplots for significant genes--------------

datMeta$Group = relevel(datMeta$Group, "Het_Saline")

pdf(file="../figures/Gene-boxplots.pdf")
for (i in 1:length(regions)){
  dge.up = read.csv(paste("../data/regulated genes lists QC/dge.up-Het", regions[i], "treatment .csv")) #load significantly upregulated genes
  dge.down = read.csv(paste("../data/regulated genes lists QC/dge.down-Het", regions[i], "treatment .csv")) #load significantly downregulated genes
  
  if (regions[i]=="all"){
    #------make plots for upregulated genes------
    if(nrow(dge.up)>0){
      for(gene_ind in 1:dim(dge.up)[1]){
        curr_gene = dge.up$X[gene_ind]
        expr_vec = datExpr[grep(curr_gene,rownames(datExpr)),]
        rownames(expr_vec) = c("Expression")
        dat = cbind(t(expr_vec),datMeta)
        print(ggplot(dat, aes(x=Group, y=Expression)) + geom_boxplot(aes(fill=Group)) + theme(legend.position="none") 
              + geom_point(position=position_jitter(.2),size=2) + facet_wrap(~Region, scales="free", nrow=3, ncol=1)
              + ggtitle(paste(dge.up$gene[gene_ind], "upregulated, padj = ", dge.up$padj[gene_ind], sep=)))
      }
    }
    
    #------make plots for downregulated genes------
    if(nrow(dge.down)>0){
      for (gene_ind in 1:dim(dge.down)[1]){
        curr_gene = dge.down$X[gene_ind]
        expr_vec = datExpr[grep(curr_gene,rownames(datExpr)),]
        rownames(expr_vec) = c("Expression")
        dat = cbind(t(expr_vec),datMeta)
        print(ggplot(dat, aes(x=Group, y=Expression)) + geom_boxplot(aes(fill=Group)) + theme(legend.position="none") 
              + geom_point(position=position_jitter(.2),size=2) + facet_wrap(~Region, scales="free", nrow=3, ncol=1)
              + ggtitle(paste(dge.down$gene[gene_ind], "downregulated, padj = ", dge.down$padj[gene_ind], sep=)))
      }
    }
  }
  
  else{
    #------make plots for upregulated genes------
    if(nrow(dge.up)>0){
      for (gene_ind in 1:dim(dge.up)[1]){
        curr_gene = dge.up$X[gene_ind]
        expr_vec = datExpr[grep(curr_gene,rownames(datExpr)),]
        rownames(expr_vec) = c("Expression")
        sample_inds = grep(toupper(regions[i]),datMeta$Region)
        dat = cbind(t(expr_vec[sample_inds]),datMeta[sample_inds,])
        print(ggplot(dat, aes(x=Group, y=Expression)) + geom_boxplot(aes(fill=Group)) + theme(legend.position="none")
              + geom_point(position=position_jitter(.2),size=2) + ggtitle(paste(dge.up$gene[gene_ind], regions[i], "upregulated, padj = ", dge.up$padj[gene_ind], sep=)))
      }
    }
    
    #------make plots for downregulated genes------
    if(nrow(dge.down)>0){
      for (gene_ind in 1:dim(dge.down)[1]){
        curr_gene = dge.down$X[gene_ind]
        expr_vec = datExpr[grep(curr_gene,rownames(datExpr)),]
        rownames(expr_vec) = c("Expression")
        sample_inds = grep(toupper(regions[i]),datMeta$Region)
        dat = cbind(t(expr_vec[sample_inds]),datMeta[sample_inds,])
        print(ggplot(dat, aes(x=Group, y=Expression)) + geom_boxplot(aes(fill=Group)) + theme(legend.position="none")
              + geom_point(position=position_jitter(.2),size=2) + ggtitle(paste(dge.down$gene[gene_ind], regions[i], "downregulated, padj = ", dge.down$padj[gene_ind], sep=)))
      }
    }
  }
}

dev.off()


#------Plot results from the 3 differential analysis methods against one another----
reg_genes_DESeq = dge
reg_genes_edgeR = top_edgeR
genes = rownames(reg_genes_DESeq)[order(reg_genes_DESeq$log2FoldChange)]  #sort by increasing FC
genes = intersect(genes, rownames(reg_genes_edgeR)) #take genes that were significant by all methods
reg_genes_DESeq = reg_genes_DESeq[genes,]
reg_genes_limma = reg_genes_limma[genes,]
gene_symbols = reg_genes_limma$ID
reg_genes_edgeR = reg_genes_edgeR[genes,]
reg_genes_edgeR$Symbol = gene_symbols
FC_DESeq = reg_genes_DESeq$log2FoldChange
FC_limma = reg_genes_limma$logFC
FC_edgeR = reg_genes_edgeR$logFC


pdf(file="../figures/DE-comparisons.pdf")

plot(FC_DESeq, FC_limma, main = paste("Sig only FC, r = ", signif(cor(FC_DESeq, FC_limma),digits=4)))
text(FC_DESeq, FC_limma, gene_symbols, pos = 3)

plot(FC_DESeq, FC_edgeR, main = paste("Sig only FC, r = ", signif(cor(FC_DESeq, FC_edgeR),digits=4)))
text(FC_DESeq, FC_edgeR, gene_symbols, pos = 3)

plot(FC_edgeR, FC_limma, main = paste("Sig only FC, r = ", signif(cor(FC_edgeR, FC_limma),digits=4)))
text(FC_edgeR, FC_limma, gene_symbols, pos = 3)

# ----Compare differential expression (FC) across all genes----

# edgeR_beta eliminates some genes when running glm -> use only genes returned by all 3 methods
keep_inds = match(edgeR_genes,rownames(datExpr))

plot(FC_DESeq_all[keep_inds], FC_limma_all[keep_inds], main = paste("FC across all genes, r = ", signif(cor(FC_DESeq_all[keep_inds], FC_limma_all[keep_inds]),digits=4)))
plot(FC_DESeq_all[keep_inds], FC_edgeR_all, main = paste("FC across all genes, r = ", signif(cor(FC_DESeq_all[keep_inds], FC_edgeR_all),digits=4)))
plot(FC_edgeR_all, FC_limma_all[keep_inds], main = paste("FC across all genes, r = ", signif(cor(FC_edgeR_all, FC_limma_all[keep_inds]),digits=4)))

# -----Compare effect sizes across all genes-----
plot(DESeq_beta[keep_inds], limma_beta[keep_inds], main = paste("Effect sizes across all genes, r = ", signif(cor(DESeq_beta[keep_inds], limma_beta[keep_inds]),digits=4)))
plot(DESeq_beta[keep_inds], edgeR_beta, main = paste("Effect sizes across all genes, r = ", signif(cor(DESeq_beta[keep_inds], edgeR_beta),digits=4)))
plot(edgeR_beta, limma_beta[keep_inds], main = paste("Effect sizes across all genes, r = ", signif(cor(edgeR_beta, limma_beta[keep_inds]),digits=4)))


dev.off()


