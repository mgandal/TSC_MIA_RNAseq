library(WGCNA)
rm(list=ls())

load("../data/kallisoForDGE.Rdata")
dds = DESeqDataSetFromTximport(txi, datMeta,design = ~Treatment + Region + Hemisphere + RIN + seqPC1 + seqPC2)


#Filter Genes: 10+ counts in half samples
genes_to_keep = apply(counts(dds)>=10,1,sum) >= round(0.5 * ncol(dds))
table(genes_to_keep)
dds = dds[genes_to_keep,]
datProbes = datProbes[genes_to_keep,]

##Network Analysis
cqn.dat <- cqn(counts(dds),lengths = as.numeric(datProbes$length), x = as.numeric(datProbes$percentage_gc_content),
               lengthMethod=c("smooth"),sqn=FALSE) ## Run cqn with specified depths and with no quantile normalization
normalizationFactors(dds) <- exp(cqn.dat$glm.offset)

datExpr.vst = assay(vst(dds,blind=FALSE))

multiExpr = vector(mode="list", length= 3)
multiExpr[[1]]= list(data=as.data.frame(t(datExpr.vst[,datMeta$Region == "PFC"])))
multiExpr[[2]]= list(data=as.data.frame(t(datExpr.vst[,datMeta$Region == "HC"])))
multiExpr[[3]]= list(data=as.data.frame(t(datExpr.vst[,datMeta$Region == "CBL"])))

gsg = goodSamplesGenesMS(multiExpr)
for(i in 1:length(multiExpr)) multiExpr[[i]]$data = multiExpr[[i]]$data[,gsg$goodGenes]

bsize = 6000; powers = c(seq(1,9,by=1),seq(10,30,by=2))

for(i in 1:length(multiExpr)) {
  sft=multiExpr[[i]]$sft=pickSoftThreshold(data= multiExpr[[i]]$data, networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers,blockSize = bsize)
  
  par(mfrow=c(1,2))
  plot(sft$fitIndices[,1], sft$fitIndices[,1], ylim=c(-1,1),xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,4], labels = powers, cex = 0.7, col="blue",  xlab="Soft Thresh Power", ylab="Scale free R^2")
  text(0,1,labels = "SFT.R.sq", col="red",cex=.7,pos = 4)
  text(0,.95,labels = "truncated.R.sq", col="blue",cex=.7,pos = 4)
  abline(h=0.8, col="black")
  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
}

softPower=16; cquant=0.5
net = blockwiseConsensusModules(multiExpr,maxBlockSize=25000,power=softPower,mergeCutHeight=0.25,nThreads=10,
                                networkCalibration  = 'single quantile',
                                calibrationQuantile = 0.8,consensusQuantile = cquant,
                                corType="bicor", useMean = T,
                                minModuleSize=20,pamStage=FALSE,
                                reassignThresholdPS = 1e-10,verbose=Inf,networkType="signed",
                                saveConsensusTOMs = TRUE, saveIndividualTOMs = FALSE,
                                consensusTOMFileNames = "../data/WGCNA_consensusTOM%b.RData")





load(net$TOMFiles)
# 
# 
# net = blockwiseModules(datExpr=multiExpr[[1]]$data, maxBlockSize=46339, networkType="signed", corType = "bicor" ,  
#                        power = 7, 
#                        mergeCutHeight= 0.1, saveTOMFileBase=paste("../data/WGCNA_network_signed_exprSet_noCBL"), saveTOMs=TRUE, 
#                        minModuleSize= 50, pamStage=FALSE, reassignThreshold=1e-6, verbose = Inf, deepSplit=2)
# load("../data/WGCNA_network_signed_exprSet_noCBL-block.1.RData")
# 

geneTree= hclust(1-consTomDS, method="average")
tree =cutreeHybrid(dendro = geneTree, minClusterSize= 100, pamStage=FALSE, cutHeight = 0.999, deepSplit=2, distM=as.matrix(1-consTomDS))
plotDendroAndColors(geneTree, labels2colors(tree$labels), "Modules", dendroLabels = FALSE)

# merged = mergeCloseModules(exprData= t(datExpr.vst[gsg$goodGenes,]), colors = tree$labels, cutHeight=0.1)
# plotDendroAndColors(geneTree, labels2colors(merged$colors), "Modules", dendroLabels = FALSE)



if(FALSE) {
  geneCov = matrix(NA, nrow=nrow(datExpr.vst[gsg$goodGenes,]),ncol=9)
  rownames(geneCov) = rownames(datExpr.vst)[gsg$goodGenes]
  colnames(geneCov) = c("PFC", "HC", "CBL", "Genotype", "Treatment", "RIN", "RNA_concentration", "seqPC1", "seqPC2")
  for(i in 1:nrow(geneCov)) {
    if(i%%1000==0) print(i)
    geneCov[i,"PFC"] = bicor(datExpr.vst[i,], (datMeta$Region=="PFC"))^2
    geneCov[i,"HC"] = bicor(datExpr.vst[i,], (datMeta$Region=="HC"))^2
    geneCov[i,"CBL"] = bicor(datExpr.vst[i,], (datMeta$Region=="CBL"))^2
    for(j in 4:ncol(geneCov)) {
      var = colnames(geneCov)[j]
      geneCov[i,j] = bicor(datExpr.vst[i,], as.numeric(datMeta[,var]))^2
    }
  }
  geneCov.color = numbers2colors(geneCov,signed=F,lim = c(0,1))
}

colors=labels2colors(tree$labels) 
plotDendroAndColors(geneTree, cbind(colors,geneCov.color), groupLabels = c("Modules", colnames(geneCov)),dendroLabels=FALSE)



MEs = moduleEigengenes(t(datExpr.vst[gsg$goodGenes,]), colors,excludeGrey = T)
tg_only = datMeta$Genotype=="Het"
i=0
i=i+1; summary(lm(MEs$eigengenes[tg_only,i] ~ Treatment + Region + Hemisphere + RIN + seqPC1 + seqPC2, data=datMeta[tg_only,]))


boxplot(MEs$eigengenes$MEcyan~ datMeta$Group)

pairwise.t.test(MEs$eigengenes[tg_only,i], datMeta$Treatment[tg_only] : datMeta$Region[tg_only])
pairwise.t.test(MEs$eigengenes[,i], datMeta$Group : datMeta$Region)


kME = signedKME(t(datExpr.vst[gsg$goodGenes,]),MEs$eigengenes)


