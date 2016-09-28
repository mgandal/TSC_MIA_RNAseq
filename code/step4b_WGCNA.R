#WGCNA

##Network Analysis
datMeta = datMeta
datExpr.vst = varianceStabilizingTransformation(dds.global)
datExpr.vst = assay(datExpr.vst)
dE = datExpr.vst


#set up expression MS
multiExpr = vector(mode="list", length=1)
multiExpr[[1]] = list(data=as.data.frame(t(datExpr.vst)), meta=datMeta)
multiExpr[[2]] = list(data=as.data.frame(t(dE[,datMeta$Region=="CBL"])), meta = datMeta.cbl)
multiExpr[[3]] = list(data=as.data.frame(t(dE[,datMeta$Region=="HC"])), meta= datMeta.hc)
multiExpr[[4]] = list(data=as.data.frame(t(dE[,datMeta$Region=="PFC"])), meta=datMeta.pfc)
names(multiExpr) = c("ALL", "Cortical", "HC", "DLPFC")


#Subselect genes to filter out ones with poor expression
g <- goodSamplesGenesMS(multiExpr)
good_genes <- which(g$goodGenes == TRUE)
dE = dE[good_genes,]

multiExpr[[1]] = list(data=as.data.frame(t(datExpr.vst)), meta=datMeta)
multiExpr[[2]] = list(data=as.data.frame(t(dE[,datMeta$Region=="CBL"])), meta = datMeta.cbl)
multiExpr[[3]] = list(data=as.data.frame(t(dE[,datMeta$Region=="HC"])), meta= datMeta.hc)
multiExpr[[4]] = list(data=as.data.frame(t(dE[,datMeta$Region=="PFC"])), meta=datMeta.pfc)

write.csv(multiExpr, "../data/multiExpr.csv")



##Step 1 Choose Soft Threshold Power 

if(TRUE)
{
  pdf(file="../figures/Soft threshold power graphs.pdf")
  bsize=6000
  powers=seq(2,30,by=2)
  for(n in 1:length(multiExpr)) { 
    #multiExpr[[n]]$softThresh = pickSoftThreshold(data= multiExpr[[n]]$data, networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers,blockSize = bsize)
    
    sft = multiExpr[[n]]$softThresh
    par(mfrow=c(1,2))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Thresh Power", ylab="Scale free R^2",type="n", main=names(multiExpr)[n])
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2",main=names(multiExpr)[n])
    abline(h=0.8, col="black")
    plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
  }
  graphics.off()
}



wgcna_parameters = list(powers =  c(18,14,14,10)) #confirm these with Mike
wgcna_parameters$minModSize = 100
wgcna_parameters$minHeight = 0.1
wgcna_parameters$bsize = 18000
wgcna_parameters$ds = 2
wgcna_parameters$networkType = "signed"
wgcna_parameters$corFnc = "bicor"
wgcna_parameters$pamStage = TRUE

if(TRUE) {
  
  for (n in 1:length(multiExpr)) {
    t_start <- Sys.time()
    ##Calculate TOM, save to file
    rootdir = getwd()
    filenm <- paste(rootdir, "/processed_data/WGCNA/network_signed_exprSet_cqn.noregress_", as.character(n),sep="")
    multiExpr[[n]]$netData = blockwiseModules(datExpr=multiExpr[[n]]$data, maxBlockSize=wgcna_parameters$bsize, networkType=wgcna_parameters$networkType, corType = wgcna_parameters$corFnc ,  power = wgcna_parameters$powers[n], mergeCutHeight= wgcna_parameters$minHeight, nThreads=23, 
                                              saveTOMFileBase=filenm, saveTOMs=TRUE, minModuleSize= wgcna_parameters$minModSize, pamStage=wgcna_parameters$pamStage, reassignThreshold=1e-6, verbose = 3, deepSplit=wgcna_parameters$ds)
    t_end <- Sys.time()
    t = t_end-t_start
    print(t)
  }
}



#merge modules and produce dendrograms for each region
pdf(file="../figures/Modules figures.pdf")

for (set.idx in 1:length(multiExpr)){
  
  load(paste("./processed_data/WGCNA/network_signed_exprSet_cqn.noregress_", as.character(set.idx),sep="", "-block.1.RData"))
  
  geneTree = hclust(as.dist(1-TOM), method = "average"); 
  colors = vector(mode="list"); labels = vector(mode="list"); labels=""
  pam=F; minModSize=100; ds=2; dthresh=0.1
  tree = cutreeHybrid(dendro = geneTree, minClusterSize= minModSize, pamStage=pam, cutHeight = 0.999, deepSplit=ds, distM=as.matrix(1-TOM))
  merged = mergeCloseModules(exprData= multiExpr[[set.idx]]$data, colors = tree$labels, cutHeight=dthresh)
  colors = labels2colors(merged$colors)
  #plotDendroAndColors(geneTree, colors, groupLabels=labels, addGuide= TRUE, dendroLabels=FALSE, main="Dendrogram", cex.colorLabels=0.5)
  #table(colors)
  
  MEs = moduleEigengenes(expr = (multiExpr[[set.idx]]$data), colors = labels2colors(merged$colors), softPower = wgcna_parameters$powers[set.idx])
  kME = signedKME(multiExpr[[set.idx]]$data, MEs$eigengenes,corFnc = "bicor")

  ##Gene-level WGCNA
  
  #for aggregate of regions
  if (regions[set.idx] == "all"){
    traitmat = as.matrix(model.matrix(~0+datMeta$Region + datMeta$Genotype + datMeta$Treatment + datMeta$Hemisphere + datMeta$RIN))
    rownames(traitmat) = rownames(datMeta)
    
    geneSigs = matrix(NA, nrow= ncol(traitmat), ncol = ncol(multiExpr[[set.idx]]$data))
    
    
    #loop through genes to find correlation of that gene's expression with the various traits
    for (i in 1:ncol(geneSigs)) {
      exprvec = as.numeric(multiExpr[[set.idx]]$data[,i])
      
      RegionCBLr = bicor(traitmat[,"datMeta$RegionCBL"], exprvec, use="pairwise.complete.obs")
      RegionHCr = bicor(traitmat[,"datMeta$RegionHC"], exprvec, use="pairwise.complete.obs")
      RegionPFCr = bicor(traitmat[,"datMeta$RegionPFC"], exprvec, use="pairwise.complete.obs")
      Genotyper = bicor(traitmat[,"datMeta$GenotypeWT"], exprvec, use="pairwise.complete.obs")
      Treatmentr = bicor(traitmat[,"datMeta$TreatmentSaline"], exprvec, use="pairwise.complete.obs")
      Hemispherer = bicor(traitmat[,"datMeta$HemisphereR"], exprvec, use="pairwise.complete.obs")
      rinr = bicor(traitmat[,"datMeta$RIN"], exprvec, use="pairwise.complete.obs")
      geneSigs[,i] = c(RegionCBLr, RegionHCr, RegionPFCr, Genotyper, Treatmentr, Hemispherer, rinr)
    }
    
    
    #for individual regions
    #generate traitmat for that region
  } else{
      if(regions[set.idx] == "cbl"){
      traitmat = as.matrix(model.matrix(~0+datMeta.cbl$Genotype + datMeta.cbl$Treatment + datMeta.cbl$Hemisphere + datMeta.cbl$RIN))
      rownames(traitmat) = rownames(datMeta.cbl)
      } else if(regions[set.idx] == "hc") {
        traitmat = as.matrix(model.matrix(~0+datMeta.hc$Genotype + datMeta.hc$Treatment + datMeta.hc$Hemisphere + datMeta.hc$RIN))
        rownames(traitmat) = rownames(datMeta.hc)
      } else if(regions[set.idx] == "pfc") {
        traitmat = as.matrix(model.matrix(~0+datMeta.pfc$Genotype + datMeta.pfc$Treatment + datMeta.pfc$Hemisphere + datMeta.pfc$RIN))
        rownames(traitmat) = rownames(datMeta.pfc)
      }
    geneSigs = matrix(NA, nrow= ncol(traitmat), ncol = ncol(multiExpr[[set.idx]]$data))
    
    
    #loop through genes to find correlation of that gene's expression with the various traits
    for (i in 1:ncol(geneSigs)) {
      exprvec = as.numeric(multiExpr[[set.idx]]$data[,i])
      
      Genotyper = bicor(traitmat[,"datMeta$GenotypeWT"], exprvec, use="pairwise.complete.obs")
      Treatmentr = bicor(traitmat[,"datMeta$TreatmentSaline"], exprvec, use="pairwise.complete.obs")
      Hemispherer = bicor(traitmat[,"datMeta$HemisphereR"], exprvec, use="pairwise.complete.obs")
      rinr = bicor(traitmat[,"datMeta$RIN"], exprvec, use="pairwise.complete.obs")
      geneSigs[,i] = c(Genotyper, Treatmentr, Hemispherer, rinr)
    }
  }
  
  multiExpr[[set.idx]]$geneSigs = geneSigs
  
  
  geneSigsColors <- matrix(0, dim(geneSigs)[1], dim(geneSigs)[2])
  for (i in 1:ncol(traitmat)) {
    geneSigsColors[i,] = numbers2colors(as.numeric(geneSigs[i,])^2, blueWhiteRed(100), signed=FALSE, centered = TRUE, lim=c(0,1))
  }
  multiExpr[[set.idx]]$geneColors = geneSigsColors
  
  colors = cbind(colors, t(geneSigsColors))
  labels = c(multiExpr[[set.idx]]$netData$cutParameters, colnames(traitmat))
  
  plotDendroAndColors(multiExpr[[set.idx]]$netData$dendrograms[[1]], colors=colors, groupLabels=c("Modules", labels), dendroLabels=FALSE, main =regions[set.idx])
  
}
dev.off()


