#WGCNA

##Network Analysis
datMeta = datMeta
datExpr.vst = varianceStabilizingTransformation(dds.global)
datExpr.vst = assay(datExpr.vst)
dE = datExpr.vst

X = model.matrix(~RIN + Hemisphere + seqPC1 + seqPC2, data=datMeta) 
beta  = (solve(t(X) %*% X) %*%t(X)) %*% t(Y)

dE.regressed = dE - t(X[,2:5] %*% beta[2:5,])


#set up expression MS
multiExpr = vector(mode="list", length=1)
multiExpr[[1]] = list(data=as.data.frame(t(dE.regressed)), meta=datMeta)
multiExpr[[2]] = list(data=as.data.frame(t(dE.regressed[,datMeta$Region=="CBL"])), meta = datMeta.cbl)
multiExpr[[3]] = list(data=as.data.frame(t(dE.regressed[,datMeta$Region=="HC"])), meta= datMeta.hc)
multiExpr[[4]] = list(data=as.data.frame(t(dE.regressed[,datMeta$Region=="PFC"])), meta=datMeta.pfc)
names(multiExpr) = c("ALL", "Cortical", "HC", "PFC")


#Subselect genes to filter out ones with poor expression
g <- goodSamplesGenesMS(multiExpr)
good_genes <- which(g$goodGenes == TRUE)
dE.regressed = dE.regressed[good_genes,]

multiExpr[[1]] = list(data=as.data.frame(t(dE.regressed)), meta=datMeta)
multiExpr[[2]] = list(data=as.data.frame(t(dE.regressed[,datMeta$Region=="CBL"])), meta = datMeta.cbl)
multiExpr[[3]] = list(data=as.data.frame(t(dE.regressed[,datMeta$Region=="HC"])), meta= datMeta.hc)
multiExpr[[4]] = list(data=as.data.frame(t(dE.regressed[,datMeta$Region=="PFC"])), meta=datMeta.pfc)

write.csv(multiExpr, "../data/multiExpr regr.csv")

##Step 1 Choose Soft Threshold Power 

if(TRUE)
{
  pdf(file="../figures/Soft threshold power graphs regr.pdf")
  bsize=6000
  powers=seq(2,30,by=2)
  for(n in 1:length(multiExpr)) { 
    multiExpr[[n]]$softThresh = pickSoftThreshold(data= multiExpr[[n]]$data, networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers,blockSize = bsize)
    
    sft = multiExpr[[n]]$softThresh
    par(mfrow=c(1,2))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Thresh Power", ylab="Scale free R^2",type="n", main=names(multiExpr)[n])
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2",main=names(multiExpr)[n])
    abline(h=0.8, col="black")
    plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
  }
  graphics.off()
  dev.off()
}


#chosen_powers = c(18,14,14,10) #non-QC model
#chosen_powers = c(18, 8, 10, 10) #QC model
chosen_powers = c(8, 8, 8, 8) #regressed model
wgcna_parameters = list(powers = chosen_powers)
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
    filenm <- paste(rootdir, "/processed_data/WGCNA QC/network_signed_exprSet_cqn.regress_", as.character(n),sep="")
    multiExpr[[n]]$netData = blockwiseModules(datExpr=multiExpr[[n]]$data, maxBlockSize=wgcna_parameters$bsize, networkType=wgcna_parameters$networkType, corType = wgcna_parameters$corFnc ,  power = wgcna_parameters$powers[n], mergeCutHeight= wgcna_parameters$minHeight, nThreads=23, 
                                              saveTOMFileBase=filenm, saveTOMs=TRUE, minModuleSize= wgcna_parameters$minModSize, pamStage=wgcna_parameters$pamStage, reassignThreshold=1e-6, verbose = 3, deepSplit=wgcna_parameters$ds)
    t_end <- Sys.time()
    t = t_end-t_start
    print(t)
    save(multiExpr, file = "multiExpr regr")
  }
}


load("multiExpr QC")

#merge modules and produce dendrograms for each region
pdf(file="../figures/Consensus modules figures regr.pdf")

pdf(file="../figures/Consensus significant modules regr.pdf")

calc_MEs = FALSE; stats = TRUE;

for (set.idx in 1:length(multiExpr)){
  
  #load(paste("./processed_data/WGCNA QC/network_signed_exprSet_cqn.regress_", as.character(set.idx),sep="", "-block.1.RData"))
  load(paste("./processed_data/WGCNA QC/tsc_consensus_regr-small-block1.Rdata"))
  
  if(calc_MEs){
  geneTree = hclust(as.dist(1-TOM), method = "average"); 
  colors = vector(mode="list"); labels = vector(mode="list"); labels=""
  pam=F; minModSize=100; ds=2; dthresh=0.1
  tree = cutreeHybrid(dendro = geneTree, minClusterSize= minModSize, pamStage=pam, cutHeight = 0.999, deepSplit=ds, distM=as.matrix(1-TOM))
  merged = mergeCloseModules(exprData= multiExpr[[set.idx]]$data, colors = tree$labels, cutHeight=dthresh)
  colors = labels2colors(merged$colors)
  #plotDendroAndColors(geneTree, colors, groupLabels=labels, addGuide= TRUE, dendroLabels=FALSE, main="Dendrogram", cex.colorLabels=0.5)
  #table(colors)
  
  MEs = moduleEigengenes(expr = (multiExpr[[set.idx]]$data), colors = labels2colors(merged$colors), softPower = wgcna_parameters$powers[set.idx])
  save(MEs, file = paste("../data/MEs/MEs cons regr", regions[set.idx], sep = "-"))
  
  kME = signedKME(multiExpr[[set.idx]]$data, MEs$eigengenes,corFnc = "bicor")
  save(kME, file = paste("../data/MEs/kMEs cons regr", regions[set.idx], sep = "-"))
  }
  
  
  #subselect correct samples from datMeta
  datMeta_curr = datMeta
  if(regions[set.idx]!="all"){
    datMeta_curr = datMeta[which(datMeta$Region == toupper(regions[set.idx])),]
  }
  
  #load(paste("../data/MEs/MEs cons regr", regions[set.idx], sep = "-"))
  #load(paste("../data/MEs/kMEs cons regr", regions[set.idx], sep = "-"))

  
  #Loop through each eigengene and assess statistical significance with group
  if(stats){
  for(i in 1:ncol(MEs)) { #MEs$eigengenes
    m = colnames(MEs)[i]; #MEs$eigengenes
    #c = substr(m,3,nchar(m))
    c = substr(m,8,nchar(m))
    
    dat= cbind(data.frame(ME = unlist(MEs[m])), datMeta_curr) #MEs$eigengenes
    #ggplot(dat, aes(x=MIA, y=ME)) + geom_boxplot() + geom_point(aes(color=Group), position=position_jitter(.1),size=3) + facet_wrap(~Region) + ggtitle(m)
    
    expr = MEs[,i] #MEs$eigengenes
    
    if(regions[set.idx]!="all"){
      a = summary(lm(expr ~ Genotype + Treatment, data=datMeta_curr)) 
      } else{
        a = summary(lm(expr ~ Region + Genotype + Treatment, data=datMeta_curr)) 
      }
    
    
    gen_p = a$coefficients["GenotypeHet","Pr(>|t|)"]
    treat_p = a$coefficients["TreatmentPolyIC","Pr(>|t|)"]
    
    
    if(gen_p < 0.05 || treat_p < 0.05) {
      title <- paste(regions[set.idx],m)
      if(gen_p < 0.05){
        title <- paste(title, "Gen", gen_p)
        boxplot(expr ~ datMeta_curr$Genotype, main = paste(regions[set.idx],m,"Gen"))
      } 
      if(treat_p < 0.05){
        title <- paste(title, "Treat", treat_p)
        boxplot(expr ~ datMeta_curr$Treatment, main = paste(regions[set.idx],m,"Treat"))
      }
      
      
      #plot hub genes
      hub_genes=colnames(multiExpr[[set.idx]]$data)[order(kME[,i],decreasing = T)[1:25]]
      hub_gene.symbol = datProbes$external_gene_name[order(kME[,i],decreasing = T)[1:25]]
      hub_gene.symbol[hub_gene.symbol==""] = hub_genes[hub_gene.symbol==""]
      gene_idx = match(hub_genes,colnames(multiExpr[[set.idx]]$data))
      adjMat = adjacency((multiExpr[[set.idx]]$data[,gene_idx]),type = "signed",corFnc = "bicor")
      adjMat[adjMat<=0]=0
      g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=FALSE)
      
      plot.igraph(g1, vertex.label = hub_gene.symbol,
                  vertex.label.dist=0.3, 
                  vertex.size=4,
                  vertex.label.color="black",
                  vertex.color = c,
                  vertex.label.cex=0.7,
                  vertex.frame.color="black",
                  layout=layout.fruchterman.reingold(g1),
                  edge.color="green",
                  main = title)
      
      
      c = substr(m,3,nchar(m))
      go.mus = gprofiler(datProbes$ensembl_gene_id[colors==c], organism="mmusculus", custom_bg = datProbes$ensembl_gene_id, 
                         correction_method = "fdr",hier_filtering = "moderate", ordered_query = F, significant = T, exclude_iea = F,
                         region_query = F,max_p_value = 0.05, max_set_size=1000, numeric_ns = "",
                         include_graph = F,src_filter = c("GO", "KEGG", "REACTOME"))
      go = go.mus[order(go.mus$p.value),]
      
      #par(oma=c(0,15,0,0));
      if(nrow(go)>0) {
        bp = barplot(-log10(as.numeric(na.omit(go$p.value[10:1]))), main=paste(c, "Module"), horiz=T, yaxt='n', col="blue", xlab='-log10(p)',cex.main=0.7, cex.axis = .7)
        axis(2,at=bp,labels=na.omit(go$term.name[10:1]),tick=FALSE,las=2,cex.axis=.7);
        abline(v=-log10(0.05), col="red", lwd=2,lty=2)
      }
    }
  }
  }
#}
#dev.off()
  

  ##Gene-level WGCNA
  
  #for aggregate of regions
  if (regions[set.idx] == "all"){
    traitmat = as.matrix(model.matrix(~0+datMeta$Region + datMeta$Genotype + datMeta$Treatment))
    rownames(traitmat) = rownames(datMeta)
    
    geneSigs = matrix(NA, nrow= ncol(traitmat), ncol = ncol(multiExpr[[set.idx]]$data))
    
    
    #loop through genes to find correlation of that gene's expression with the various traits
    for (i in 1:ncol(geneSigs)) {
      exprvec = as.numeric(multiExpr[[set.idx]]$data[,i])
      
      RegionCBLr = bicor(traitmat[,"datMeta$RegionCBL"], exprvec, use="pairwise.complete.obs")
      RegionHCr = bicor(traitmat[,"datMeta$RegionHC"], exprvec, use="pairwise.complete.obs")
      RegionPFCr = bicor(traitmat[,"datMeta$RegionPFC"], exprvec, use="pairwise.complete.obs")
      Genotyper = bicor(traitmat[,"datMeta$GenotypeHet"], exprvec, use="pairwise.complete.obs")
      Treatmentr = bicor(traitmat[,"datMeta$TreatmentPolyIC"], exprvec, use="pairwise.complete.obs")
      #Hemispherer = bicor(traitmat[,"datMeta$HemisphereR"], exprvec, use="pairwise.complete.obs")
      #rinr = bicor(traitmat[,"datMeta$RIN"], exprvec, use="pairwise.complete.obs")
      #seqPC1r = bicor(traitmat[,"datMeta$seqPC1"], exprvec, use="pairwise.complete.obs")
      #seqPC2r = bicor(traitmat[,"datMeta$seqPC2"], exprvec, use="pairwise.complete.obs")
      geneSigs[,i] = c(RegionCBLr, RegionHCr, RegionPFCr, Genotyper, Treatmentr)
   }
    
  
    
    #for individual regions
    #generate traitmat for that region
  } else{
      if(regions[set.idx] == "cbl"){
      traitmat = as.matrix(model.matrix(~datMeta.cbl$Genotype + datMeta.cbl$Treatment)) # + datMeta.cbl$Hemisphere + datMeta.cbl$RIN + datMeta.cbl$seqPC1 + datMeta.cbl$seqPC2))
      traitmat = traitmat[,-1]
      rownames(traitmat) = rownames(datMeta.cbl)
      } else if(regions[set.idx] == "hc") {
        traitmat = as.matrix(model.matrix(~datMeta.hc$Genotype + datMeta.hc$Treatment)) # + datMeta.hc$Hemisphere + datMeta.hc$RIN + datMeta.hc$seqPC1 + datMeta.hc$seqPC2))
        traitmat = traitmat[,-1]
        rownames(traitmat) = rownames(datMeta.hc)
      } else if(regions[set.idx] == "pfc") {
        traitmat = as.matrix(model.matrix(~datMeta.pfc$Genotype + datMeta.pfc$Treatment)) # + datMeta.pfc$Hemisphere + datMeta.pfc$RIN + datMeta.pfc$seqPC1 + datMeta.pfc$seqPC2))
        traitmat = traitmat[,-1]
        rownames(traitmat) = rownames(datMeta.pfc)
      }
    colnames(traitmat) = c("datMeta$GenotypeHet", "datMeta$TreatmentPolyIC") # "datMeta$HemisphereR", "datMeta$RIN", "datMeta$seqPC1", "datMeta$seqPC2")
    geneSigs = matrix(NA, nrow= ncol(traitmat), ncol = ncol(multiExpr[[set.idx]]$data))
    
    
    #loop through genes to find correlation of that gene's expression with the various traits
    for (i in 1:ncol(geneSigs)) {
      exprvec = as.numeric(multiExpr[[set.idx]]$data[,i])
      
      Genotyper = bicor(traitmat[,"datMeta$GenotypeHet"], exprvec, use="pairwise.complete.obs")
      Treatmentr = bicor(traitmat[,"datMeta$TreatmentPolyIC"], exprvec, use="pairwise.complete.obs")
      #Hemispherer = bicor(traitmat[,"datMeta$HemisphereR"], exprvec, use="pairwise.complete.obs")
      #rinr = bicor(traitmat[,"datMeta$RIN"], exprvec, use="pairwise.complete.obs")
      #seqPC1r = bicor(traitmat[,"datMeta$seqPC1"], exprvec, use="pairwise.complete.obs")
      #seqPC2r = bicor(traitmat[,"datMeta$seqPC2"], exprvec, use="pairwise.complete.obs")
      
      geneSigs[,i] = c(Genotyper, Treatmentr) #, Hemispherer, rinr, seqPC1r, seqPC2r)
    }
  }
  
  multiExpr[[set.idx]]$geneSigs = geneSigs
  
  #come up with colors to represent geneSigs
  geneSigsColors <- matrix(0, dim(geneSigs)[1], dim(geneSigs)[2])
  for (i in 1:ncol(traitmat)) {
    #geneSigsColors[i,] = numbers2colors(as.numeric(geneSigs[i,])^2, blueWhiteRed(100), signed=FALSE, centered = TRUE, lim=c(0,1))
    geneSigsColors[i,] = numbers2colors(as.numeric(geneSigs[i,]), blueWhiteRed(100), signed=TRUE, centered = TRUE, lim=c(-1,1))
    if (i >= ncol(traitmat)-1){
      geneSigsColors[i,] = numbers2colors(as.numeric(geneSigs[i,]), blueWhiteRed(100), signed=TRUE, centered = TRUE, lim=c(-1,1))
    }
  }
  multiExpr[[set.idx]]$geneColors = geneSigsColors
  
  colors = cbind(labels2colors(merged$colors), t(geneSigsColors))
  trait_labels = substr(colnames(traitmat), 9, nchar(colnames(traitmat))) #remove "datMeta$" from variable names
  
  plotDendroAndColors(multiExpr[[set.idx]]$netData$dendrograms[[1]], colors=colors, groupLabels=c("Modules", trait_labels), dendroLabels=FALSE, main =regions[set.idx])
  
}
dev.off()

save(multiExpr, file = "multiExpr regr")



#build consensus network
#----------------------
multiExpr_reg = multiExpr
multiExpr_reg[[1]] <- NULL

net = blockwiseConsensusModules(multiExpr_reg, maxBlockSize = 20000, corType = "bicor", power = chosen_powers[2:4],
                                networkType = "signed", saveIndividualTOMs = T, saveConsensusTOMs = T,
                                consensusTOMFileNames = "./processed_data/WGCNA QC/tsc_consensus_regr-small-block%b.Rdata",
                                consensusQuantile = 0.2, deepSplit = 2, pamStage = FALSE, 
                                detectCutHeight = 0.995, mergeCutHeight = 0.2, 
                                minModuleSize = 100, 
                                reassignThresholdPS = 1e-10,verbose=3)

load("./processed_data/WGCNA QC/tsc_consensus_regr-small-block1.Rdata")

matTOM = as.matrix(consTomDS)

moduleColors = net$colors
consTree = net$dendrograms[[1]]

geneSigsColors = multiExpr[[1]]$geneColors

colors = cbind(moduleColors, t(geneSigsColors))
traitmat = as.matrix(model.matrix(~0+datMeta$Region + datMeta$Genotype + datMeta$Treatment + datMeta$Hemisphere + datMeta$RIN + datMeta$seqPC1 + datMeta$seqPC2))
trait_labels = substr(colnames(traitmat), 9, nchar(colnames(traitmat))) #remove "datMeta$" from variable names



pdf(file="../figures/Consensus modules unsigned.pdf")
plotDendroAndColors(consTree, colors=colors,
                    groupLabels=c("Modules", trait_labels),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram")
dev.off()


#calculate consensus MEs and kMEs
#-------------------------------
consMEs = merge(net$multiMEs[[1]], net$multiMEs[[2]], all=TRUE)
consMEs = merge(consMEs, net$multiMEs[[3]], all=TRUE)
kME = signedKME(multiExpr[[1]]$data, consMEs, corFnc = "bicor")
MEs = consMEs


if(FALSE){
#Make web plot for consensus network
#---------------------------------
pdf("./figures/ConsensusModulesNetworks.pdf", height=8, width=8)

unique_colors = unique(moduleColors)

for(mod in unique_colors) {
  
  module_genes = colnames(multiExpr[[1]]$data)[moduleColors == mod]
  mod_name = paste("kMEta.ME", mod, sep="")
  MM = order(kME[module_genes,mod_name],decreasing = T) #sort genes in order of highest MM
  module = kME[module_genes[MM], mod_name]
  
  numgenesingraph = 100
  numconnections2keep = 500
  MMtoDisp = MM[1:numgenesingraph]
  
  idx = which(colnames(multiExpr[[1]]$data) %in% module_genes)
  #----------------
  modTOM = matTOM[idx,idx]
  dimnames(modTOM) = list(genes, genes)
  
  reducedTOM = modTOM
  orderind = order(reducedTOM,decreasing=TRUE);
  connections2keep = orderind[(1+nrow(modTOM)):(numconnections2keep+nrow(modTOM))];
  reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
  reducedTOM[connections2keep] = 1;
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[1:5,1:5]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMata <- layout.circle(g0)
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[6:29,6:29]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatb <- layout.circle(g0)
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[30:ncol(reducedTOM),30:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatc <- layout.circle(g0)
  
  g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.5, layoutMatc)
  plot(g1,edge.color="grey",vertex.color=mod,vertex.label=as.character(genes),vertex.label.cex=0.7,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=module[1:numgenesingraph]^3*4,main=paste(mod,"module"))
}

dev.off()



adjMat = covAll^1
adjMat2 = adjMat; adjMat2[adjMat2 <0.5] = 0
adjMat[(adjMat)<0.2] = 0

genes_nogrey=rownames(modColors)[modColors!="grey"]

keepgenes = rownames(cons_kme)[gene_idx]
adjMat = adjMat[gene_idx,gene_idx]
adjMat2 = adjMat2[gene_idx,gene_idx]
numcors =Inf
topcors = sort(as.numeric(adjMat), decreasing = T)[numcors]
adjMat[adjMat<=topcors]=0

geneSymbols = asd_datProbes[keepgenes, "external_gene_id"]
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=FALSE)
g2 = graph.adjacency(as.matrix(adjMat2),mode="undirected",weighted=T,diag=FALSE)
layoutFR <- layout.fruchterman.reingold(g2,dim=2)



edgecolors = numbers2colors(E(g1)$weight, colors = redWhiteGreen(100, gamma=2), signed=T, centered=T, lim=c(-0.93,0.93))
plot.igraph(g1, vertex.label = geneSymbols,
            vertex.label.dist=0.3, 
            vertex.size=4,
            vertex.label.color="black",
            vertex.color = modColors[keepgenes,"color"],
            vertex.label.cex=0.6,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=edgecolors,
            edge.width = 2*E(g1)$weight)

edgecolors = numbers2colors(E(g2)$weight, colors = redWhiteGreen(100, gamma =2), signed=T, centered=T, lim=c(-1,1))
plot.igraph(g2, vertex.label = geneSymbols, add=T,
            vertex.label.dist=0.3, 
            vertex.size=4,
            vertex.label.color="black",
            vertex.color = modColors[keepgenes,"color"],
            vertex.label.cex=0.6,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=edgecolors,
            edge.width = 4*E(g1)$weight^4)


##PPI (protein-protein interaction)
ppi3 = ppiMat
ppi3[ppiMat>0] = 1
deg3 = apply(ppi3,2,sum)
idx = match(geneSymbols, rownames(ppi3))
ppiSmall = ppi3[idx,idx]

ppiSmall[is.na(ppiSmall)]=0


g2 = graph.adjacency(as.matrix(ppiSmall), mode="undirected", weighted=T, diag=F)
plot.igraph(g2, vertex.label = "", 
            vertex.label.dist=0.3, add=T,
            vertex.size=0,
            vertex.shape="none",
            vertex.label.cex=0.6,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color="gray",
            edge.width=3)
}




# Y = datExpr.vst
# X = model.matrix(~RIN + Hemisphere + seqPC1 + seqPC2, data=datMeta) 
# beta  = (solve(t(X) %*% X) %*%t(X)) %*% t(Y)
# 
# Y.regressed = Y - t(X[,2:5] %*% beta[2:5,])
# 
# 
# expr = datExpr.vst[1,]
# beta2 = lm(expr ~ RIN + Hemisphere + seqPC1 + seqPC2, data=datMeta)

