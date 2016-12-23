options(stringsAsFactors = F)
rm(list=ls()) #Clear workspace
setwd("/Users/mgandal/Documents/Github/TSC_MIA_RNAseq/code")

#Load Libraries
library(WGCNA); library(DESeq2); library(biomaRt); library(ggplot2);
library(reshape); library(cqn); library(limma)
library(gridExtra)
library(gProfileR); library(fdrtool)
library(corrplot)
library(venneuler)
library(pSI)

if(FALSE) {
  load("../data/HTseqCounts.RData")
  datExpr.hts = datExpr; rm(datExpr)
  rownames(datExpr.hts) = substr(rownames(datExpr.hts),0,18)
  RIN_data <- read.csv("../data/Processed Silva_TSC_MIA_MsRNAseq20160510.csv")

  load("../data/datExpr.kallisto.txi.RData")
  datExpr.kalliso = txi$counts

  #Organize Meta Data
  datMeta$Genotype = "Het"
  datMeta$Genotype[datMeta$Subject %in% c(442,466,469)] = "WT"
  datMeta$Genotype <- as.factor(datMeta$Genotype)
  datMeta$Genotype <- relevel(datMeta$Genotype, "WT")
  datMeta$Treatment = "PolyIC"
  datMeta$Treatment[datMeta$Subject %in% c(420, 455, 447)] = "Saline"
  datMeta$Treatment <- as.factor(datMeta$Treatment)
  datMeta$Treatment <- relevel(datMeta$Treatment, "Saline")

  datMeta$Group = as.factor(paste(datMeta$Genotype, "_", datMeta$Treatment,sep=""))
  datMeta$Region <- as.factor(datMeta$Region)

  idx = match(datMeta$Sample, gsub("_","-",RIN_data$Sample))
  datMeta$RIN = RIN_data$RNA.RIN[idx]
  datMeta$RNA_concentration = RIN_data$Conc...ng.ul..from.NanoChip[idx]


  #Load Sequencing Statistics from Picard
  datSeq <- read.csv("../data/QC/PicardToolsQC.csv",row.names=1)
  rownames(datSeq) = gsub("/u/home/g/gandalm/project-geschwind/TSC_MIA_Silva/data/STAR_bam/", "", rownames(datSeq))
  datSeq = datSeq[match(datMeta$Sample,rownames(datSeq)),]


  #Compute principle components of QC
  PC.datSeq = prcomp(t(scale(datSeq,scale=T)),center=F)
  varexp <- (PC.datSeq$sdev)^2 / sum(PC.datSeq$sdev^2)
  corrplot::corrplot(cor(cbind(PC.datSeq$rotation[,1:5], datSeq)))
  datMeta$seqPC1=PC.datSeq$rotation[,1]
  datMeta$seqPC2=PC.datSeq$rotation[,2]
  datMeta$seqPC3=PC.datSeq$rotation[,3]
  datMeta$seqPC4=PC.datSeq$rotation[,4]

  #Annotate Probes
  bm = useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl", host="mar2015.archive.ensembl.org")
  a = listAttributes(bm); f= listFilters(bm)
  bm1 = getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene", "chromosome_name", "start_position", "end_position", "percentage_gc_content"),
              mart=bm)
  datProbes = bm1[match(rownames(datExpr.kalliso), bm1$ensembl_gene_id),]
  datProbes$length = datProbes$end_position - datProbes$start_position

  save(file="../data/kallisoForDGE.Rdata", datExpr.kalliso, datMeta,datProbes, txi)
}


load("../data/kallisoForDGE.Rdata")

dds = DESeqDataSetFromTximport(txi, datMeta,design = ~Treatment + Region + Hemisphere + RIN + seqPC1 + seqPC2)

#Filter Genes: 10+ counts in half samples
genes_to_keep = apply(counts(dds)>=10,1,sum) >= round(0.5 * ncol(dds))
table(genes_to_keep)
dds = dds[genes_to_keep,]
datProbes = datProbes[genes_to_keep,]


#QC, Normalization
datExpr = assay(vst(dds))
par(mfrow=c(2,1))
mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.factor(datMeta$Region), pch=16,main="MDS Plot by Region", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
legend("topleft",legend = levels(as.factor(datMeta$Region)), col=1:3,pch=19)

i = 1; plot(density((datExpr[,i]+1), na.rm=T), col = as.numeric(datMeta$Region[i]), main="Expression Histogram", xlab = "log2 Normalized Counts");
for(i in 2:dim(datExpr)[2]) {     lines(density(((datExpr[,i]+1)), na.rm=T), col = as.numeric(datMeta$Region[i]),) }

#Outlier Removal --> no outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(Z.K, col = as.numeric(colData(dds)$Region), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-2, lty=2)



#CQN
cqn.dat <- cqn(counts(dds),lengths = as.numeric(datProbes$length), x = as.numeric(datProbes$percentage_gc_content),
                lengthMethod=c("smooth"),sqn=FALSE, verbose = TRUE) ## Run cqn with specified depths and with no quantile normalization
normalizationFactors(dds) <- exp(cqn.dat$glm.offset)


##Remove WT from DGE analysis
no_wt = datMeta$Genotype=="Het"
dds = dds[,no_wt]

##DESeq2: Calculate DGE
dds <- DESeq(dds,betaPrior = T)
res.all <- results(dds, contrast = c("Treatment", "PolyIC", "Saline"))
DESeq2::plotMA(res.all, ylim=c(-2,2))

#Correct pvalue
hist(res.all$pvalue, main="DESeq2 P value histogram", xlab="P value",col="grey",breaks = 50,border = "grey60")
res.all <- as.data.frame(res.all[ !is.na(res.all$pvalue), ])
corrected <- fdrtool(res.all$stat, statistic= "normal", plot = F)
res.all$qval <- corrected$qval
res.all$gene = datProbes$external_gene_name
table(res.all$qval<.1)

#Volcano Plot
plot(res.all$log2FoldChange, -log10(res.all$qval), xlab="log2 Fold Change", ylab="-log10(qval)",cex=.7, main="Treatment DGE (across all regions)", ylim=c(0,13))
points(res.all$log2FoldChange[res.all$qval<.1], -log10(res.all$qval[res.all$qval<.1]),pch=19,col="red", cex=.7)
text(res.all$log2FoldChange[res.all$qval<.1], -log10(res.all$qval[res.all$qval<.1]), labels = datProbes$external_gene_name[res.all$qval<.1],cex=.6, pos = 3)
abline(h=-log10(.1), lty=2)


#Plot Top Individual Genes
gene_plot = list()
idx = which(res.all$qval<.1)
idx = idx[order(res.all$qval[idx])]

for(i in 1:length(idx)) {
  gene = rownames(res.all)[idx[i]]
  d= plotCounts(dds, gene, intgroup=c("Group", "Region"),returnData = T, transform = TRUE)
  d$Group = relevel(d$Group, "Het_Saline")
  gene_plot[[i]]=
    ggplot(d,aes(x=Group, y=count, fill=Group))  +
    geom_boxplot(outlier.shape = NA, aes(fill=Group)) + geom_point(position=position_jitterdodge(0.8), alpha=0.5,size=1) +
    xlab("") + ylab("log2 normalized counts") + facet_grid(~Region) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust=.5),
          axis.text.x = element_blank(), axis.ticks.x =element_blank(),axis.title.x=element_blank(),
          axis.title.y = element_text(size=7), axis.text.y = element_text(size=6),
          legend.title = element_text(size=7), legend.text = element_text(size=6),
          legend.key.width=unit(0.3,"cm"))+
    ggtitle(datProbes$external_gene_name[match(gene,datProbes$ensembl_gene_id)],
            subtitle = paste("log2FC=", signif(res.all[gene, "log2FoldChange"],2), ", q=", signif(res.all[gene, "qval"],2),sep=""))
  if(i!=1) {
    gene_plot[[i]]=gene_plot[[i]] + theme(legend.position="none")
  }

}

pdf("../figures/Fig3-TopGenes.pdf",width=12,height=20)
grid.arrange(grobs=gene_plot[1:length(gene_plot)], ncol=3)
dev.off()




#Calculate DGE within each region separately
res.region = vector(mode="list",length = 4)
names(res.region)=  c("AllRegion", "PFC", "HC", "CBL")
res.region$AllRegion = res.all
par(mfrow=c(1,2))

for(region in c("HC", "PFC", "CBL")) {
  idx = (colData(dds)$Region == region)
  dds_region = dds[,idx]
  design(dds_region) = ~Treatment + RIN + seqPC1 + seqPC2

  dds_region <- DESeq(dds_region)

  res.region[[region]] = as.data.frame(results(dds_region, contrast = c("Treatment", "PolyIC", "Saline")))
  to_keep = !is.na(res.region[[region]]$stat)
  res.region[[region]] = res.region[[region]][to_keep,]
  corrected <- fdrtool(res.region[[region]]$stat, statistic= "normal", plot = F)
  res.region[[region]]$qval <- corrected$qval
  res.region[[region]]$gene = datProbes$external_gene_name[to_keep]

  hist(res.region[[region]]$pvalue, main=paste("DESeq2 P value histogram\n",region,sep=""), xlab="P value",col="grey",breaks = 50,border = "grey60")

  x = res.region[[region]]$log2FoldChange
  y = -log10(res.region[[region]]$qval)
  sig = which(res.region[[region]]$qval < .1)
  plot(x,y, xlab="log2 Fold Change", ylab="-log10(qval)",cex=.7, main=paste("Treatment DGE (", region, ")",sep=""))
  points(x[sig], y[sig],pch=19,col="red", cex=.5)
  sig = sig[order(res.region[[region]]$qval[sig])][1:20]
  text(x[sig], y[sig], labels = datProbes$external_gene_name[sig],cex=.6, pos = 3)
  abline(h=-log10(.1), lty=2)
}



#Compare DGE log2FC signature across all regions
genes = intersect(rownames(res.region[[1]]), rownames(res.region[[2]]))
for(i in 3:4) genes = intersect(genes, rownames(res.region[[i]]))

for(i in 1:4) res.region[[i]] = res.region[[i]][match(genes, rownames(res.region[[i]])),]

corrplot.mixed(cor(data.frame(AllRegion=res.region$AllRegion$log2FoldChange,
                              PFC = res.region$PFC$log2FoldChange,
                              HC = res.region$HC$log2FoldChange,
                              CBL = res.region$CBL$log2FoldChange),method="spearman"),main="log2FC Correlation")


#Compare DGE qval < 0.1 genes across all regions
dge_genes = matrix(0, nrow=length(genes), ncol=4)
colnames(dge_genes) = names(res.region)
rownames(dge_genes) = genes
for(i in 1:4) {
  idx = res.region[[i]]$qval<.1
  dge_genes[idx,i] = 1
}
vd = venneuler(dge_genes)
plot(vd, main="DGE genes (fdr < 0.1)")


#Calculate gene ontology for each region
go_results = data.frame()

for(i in 1:4) {
  region = names(res.region)[[i]]
  query.up = rownames(res.region[[i]])[res.region[[i]]$qval<.1 & res.region[[i]]$log2FoldChange>0]
  query.up = query.up[order(res.region[[i]][query.up, "log2FoldChange"], decreasing = T)]
  query.down = rownames(res.region[[i]])[res.region[[i]]$qval<.1 & res.region[[i]]$log2FoldChange<0]
  query.down = query.down[order(res.region[[i]][query.down, "log2FoldChange"], decreasing = F)]

  go.up = gprofiler(query.up, organism = "mmusculus", custom_bg = rownames(res.region[[i]]),
                    src_filter=c("KEGG", "GO"), correction_method = "fdr", ordered_query = T,
                    significant = T, exclude_iea = F,  region_query = F, max_p_value = 0.05,
                    numeric_ns = "", max_set_size = 1000,hier_filtering = "strong")
  go.up = go.up[order(go.up$p.value)[1:min(5,nrow(go.up))],]

  go.down = gprofiler(query.down,  organism = "mmusculus", custom_bg = rownames(res.region[[i]]),
                      src_filter=c("KEGG", "GO"), correction_method = "fdr", ordered_query = T,
                      significant = T, exclude_iea = F,  region_query = F, max_p_value = 0.05,
                      numeric_ns = "", max_set_size = 1000,hier_filtering = "strong")
  go.down = go.down[order(go.down$p.value)[1:min(5,nrow(go.down))],]

  go_results = rbind(go_results, data.frame(Region = region, Change="Up", go.up))
  go_results = rbind(go_results, data.frame(Region = region, Change="Down", go.down))
}


go_results = go_results[order(-log10(go_results$p.value),decreasing = F),]
go_results$Region = factor(go_results$Region,levels=names(res.region))

plot.down = ggplot(go_results[go_results$Change=="Down",],aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity") + coord_flip() +
  facet_grid(Region ~ Change, scales="free_y") +  xlab("") + geom_abline(slope = 0, intercept = -log10(.05), lty=2, lwd=.5) + theme(axis.text.y = element_text(size=6))
plot.up = ggplot(go_results[go_results$Change=="Up",],aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity", aes(x=reorder(term.name, -log10(p.value)))) + coord_flip() +
  facet_grid(Region ~ Change, scales="free_y") +  xlab("") + geom_abline(slope = 0, intercept = -log10(.05), lty=2, lwd=.5)+ theme(axis.text.y = element_text(size=6))

grid.arrange(grobs=list(plot.down, plot.up), ncol=2)



##Cell-Type Enrichment Goes Here
zhang = read.csv("../data/barres_RNAseq.csv")
bm = useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl", host="mar2015.archive.ensembl.org")
bm1 = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),filters="external_gene_name", values=zhang$Gene.symbol,mart=bm)
bm1 = bm1[match(zhang$Gene.symbol, bm1$external_gene_name),]
zhang$ensg = bm1$ensembl_gene_id
zhang = zhang[!is.na(zhang$ensg),]
rownames(zhang) = zhang$ensg
zhang = zhang[,c(3:9)]

pSI.zhang= specificity.index(zhang)

celltype_results = data.frame()
for(i in 1:length(res.region)) {
  region = names(res.region)[[i]]

  enrich.up = fisher.iteration(pSI.zhang, candidate.genes = rownames(res.region[[i]])[res.region[[i]]$qval <.1 & res.region[[i]]$log2FoldChange > 0])
  enrich.down = fisher.iteration(pSI.zhang, candidate.genes = rownames(res.region[[i]])[res.region[[i]]$qval <.1 & res.region[[i]]$log2FoldChange < 0])

  celltype_results = rbind(celltype_results,
                           data.frame(Region=region, Change="Up", CellType = rownames(enrich.up), P.value=enrich.up[,1]))
  celltype_results = rbind(celltype_results,
                           data.frame(Region=region, Change="Down", CellType = rownames(enrich.down), P.value=enrich.down[,1]))
}

plot.down = ggplot(celltype_results[celltype_results$Change=="Down",],aes(x=CellType,  y=log10(P.value), fill=CellType)) + geom_bar(stat="identity") + coord_flip() +
  facet_grid(Region ~ Change, scales="free_y") +  xlab("") + geom_abline(slope = 0, intercept = log10(.05), lty=2, lwd=.5) + theme(axis.text.y = element_text(size=6)) + theme(legend.position="none")

plot.up = ggplot(celltype_results[celltype_results$Change=="Up",],aes(x=CellType,  y=-log10(P.value), fill=CellType)) + geom_bar(stat="identity") + coord_flip() +
  facet_grid(Region ~ Change, scales="free_y") +  xlab("") + geom_abline(slope = 0, intercept = -log10(.05), lty=2, lwd=.5) + theme(axis.text.y = element_text(size=6))

grid.arrange(grobs=list(plot.down, plot.up), layout_matrix=rbind(c(1,1,2,2,2),c(1,1,2,2,2)), ncol=2)



go_down = go_results[go_results$Change=="Down",]
go_down$p.value = as.numeric(go_down$p.value)
go_down$logPval = -log10(go_down$p.value)
ggplot(go_results,aes(x=term.name, y=-log10(p.value))) + geom_bar(stat="identity") + coord_flip()


+ coord_flip() +
    xlab("") + geom_abline(slope = 0, intercept = -log10(.05), lty=2, lwd=.5) + theme(axis.text.y = element_text(size=6))





    query = rownames(datExpr)[colors=="midnightblue"]

    go.mus = gprofiler(query, organism="mmusculus") #, custom_bg = datProbes$ensembl_gene_id,
                       correction_method = "fdr",hier_filtering = "moderate", ordered_query = T, significant = T, exclude_iea = F,
                       region_query = F, max_p_value = 0.05, max_set_size=1000, numeric_ns = "",
                       include_graph = F,src_filter = c("GO", "KEGG", "REACTOME"))
    go = go.mus[order(go.mus$p.value),]

    saveRDS(go, file = paste("../data/regulated genes lists QC/go down", curr_reg, contrast, ".RDS"))


    #read in saved GO
    go <- readRDS(paste("../data/regulated genes lists QC/go down", curr_reg, contrast, ".RDS"))

    if (nrow(go) > 0){
      par(oma=c(0,15,0,0))
      ttl = paste("GO downregulated", contrast, curr_reg)
      n_go_show = min(10, dim(dge.down)[1])
      bp = barplot(-log10(as.numeric(na.omit(go$p.value[n_go_show:1]))), main = ttl, horiz=T, yaxt='n', col="blue", xlab='-log10(p)',cex.main=0.7, cex.axis = .7)
      axis(2,at=bp,labels=na.omit(go$term.name[n_go_show:1]),tick=FALSE,las=2,cex.axis=.7);
      abline(v=-log10(0.05), col="red", lwd=2,lty=2)
    }


    #upregulated GO

    #access gprofiler for GO
    query = rownames(dge.up)[order(dge.up$log2FoldChange, decreasing = T)]

    go.mus = gprofiler(query, organism="mmusculus", custom_bg = datProbes$ensembl_gene_id,
                       correction_method = "fdr",hier_filtering = "moderate", ordered_query = T, significant = T, exclude_iea = F,
                       region_query = F, max_p_value = 0.05, max_set_size=1000, numeric_ns = "",
                       include_graph = F,src_filter = c("GO", "KEGG", "REACTOME"))

    go = go.mus[order(go.mus$p.value),]

    saveRDS(go, file = paste("../data/regulated genes lists QC/go up", curr_reg, contrast, ".RDS"))



    #read in saved GO
    go <- readRDS(paste("../data/regulated genes lists QC/go up", curr_reg, contrast, ".RDS"))


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

write.csv(tally_genotype, file = "../data/regulated genes lists QC/tally_genotype.csv"); write.csv(tally_treatment, file = "../data/regulated genes lists QC/tally_treatment.csv")

dev.off()
