##Gene-level WGCNA

## Calculate Gene-Wise Statistics
## ------------------------------
scz_traitmat = cbind(as.factor(scz_datPheno$"characteristics: Group"), as.factor(scz_datPheno$"characteristics: Brain Region"), as.factor(scz_datPheno$"characteristics: Sex"), as.numeric(scz_datPheno$"characteristics: Age"), as.numeric(scz_datPheno$"characteristics: PMI(min)"))
colnames(scz_traitmat) = c("Group", "Region", "Sex", "Age", "PMI")
rownames(scz_traitmat) = rownames(scz_datPheno)
scz_geneSigs = matrix(NA, nrow= ncol(scz_traitmat), ncol = ncol(scz_multiExpr[[1]]$data))

asd_traitmat = cbind(as.factor(asd_datPheno$A.C), as.factor(asd_datPheno$Brain.area), as.factor(asd_datPheno$Sex), as.numeric(asd_datPheno$Age), as.numeric(asd_datPheno$PMI), as.numeric(asd_datPheno$RIN))
colnames(asd_traitmat) = c("Group", "Region", "Sex", "Age", "PMI", "RIN")
rownames(asd_traitmat) = rownames(asd_datPheno)
asd_geneSigs = matrix(NA, nrow= ncol(asd_traitmat), ncol = ncol(asd_multiExpr[[1]]$data))


setIndex = 1; #COMBINED DATASET

for (i in 1:ncol(scz_geneSigs)) {
  exprvec = as.numeric(scz_multiExpr[[setIndex]]$data[,i])
  
  groupr = sqrt(max(summary(lm(exprvec ~ as.factor(scz_datPheno$"characteristics: Group")))$adj.r.squared,0))
  regionr = sqrt(max(summary(lm(exprvec ~ as.factor(scz_datPheno$"characteristics: Brain Region")))$adj.r.squared,0))
  sexr = sqrt(max(summary(lm(exprvec ~ as.factor(scz_datPheno$"characteristics: Sex")))$adj.r.squared,0))
  ager = bicor(scz_traitmat[,"Age"], exprvec)
  pmir = bicor(SCZmat[,"PMI"], exprvec)
  
  scz_geneSigs[,i] = c(groupr, regionr, sexr, ager, pmir)
}

scz_multiExpr[[setIndex]]$geneSigs = scz_geneSigs
scz_geneSigs[1,] = numbers2colors(as.numeric(scz_geneSigs[1,], blueWhiteRed(100), signed=FALSE, centered = FALSE, lim=c(0,1)))
scz_geneSigs[2,] = numbers2colors(as.numeric(scz_geneSigs[2,], blueWhiteRed(100), signed=FALSE, centered = FALSE, lim=c(0,1)))
scz_geneSigs[3,] = numbers2colors(as.numeric(scz_geneSigs[3,], blueWhiteRed(100), signed=FALSE, centered = FALSE, lim=c(0,1)))
scz_geneSigs[4,] = numbers2colors(as.numeric(scz_geneSigs[4,], blueWhiteRed(100), signed=TRUE, centered = TRUE, lim=c(-1,1)))
scz_geneSigs[5,] = numbers2colors(as.numeric(scz_geneSigs[5,], blueWhiteRed(100), signed=TRUE, centered = TRUE, lim=c(-1,1)))
scz_multiExpr[[setIndex]]$geneColors = scz_geneSigs

scz_colors = cbind(labels2colors(scz_merged$colors), t(scz_geneSigs))
scz_labels = c(scz_multiExpr[[setIndex]]$netData$cutParameters, colnames(scz_traitmat))
plotDendroAndColors(scz_multiExpr[[setIndex]]$netData$dendrograms[[1]], colors=scz_colors, groupLabels=scz_labels, dendroLabels=FALSE)

for (i in 1:ncol(asd_geneSigs)) {
  exprvec = as.numeric(asd_multiExpr[[setIndex]]$data[,i])
  
  groupr = sqrt(max(summary(lm(exprvec ~ as.factor(asd_datPheno$A.C)))$adj.r.squared,0))
  regionr = sqrt(max(summary(lm(exprvec ~ as.factor(asd_datPheno$Brain.area)))$adj.r.squared,0))
  sexr = sqrt(max(summary(lm(exprvec ~ as.factor(asd_datPheno$Sex)))$adj.r.squared,0))
  ager = bicor(asd_traitmat[,"Age"], exprvec, use="pairwise.complete.obs")
  pmir = bicor(asd_traitmat[,"PMI"], exprvec, use="pairwise.complete.obs")
  rinr = bicor(asd_traitmat[,"RIN"], exprvec, use="pairwise.complete.obs")
  asd_geneSigs[,i] = c(groupr, regionr, sexr, ager, pmir, rinr)
}
asd_multiExpr[[setIndex]]$geneSigs = asd_geneSigs
asd_geneSigs[1,] = numbers2colors(as.numeric(asd_geneSigs[1,], blueWhiteRed(100), signed=FALSE, centered = FALSE, lim=c(0,1)))
asd_geneSigs[2,] = numbers2colors(as.numeric(asd_geneSigs[2,], blueWhiteRed(100), signed=FALSE, centered = FALSE, lim=c(0,1)))
asd_geneSigs[3,] = numbers2colors(as.numeric(asd_geneSigs[3,], blueWhiteRed(100), signed=FALSE, centered = FALSE, lim=c(0,1)))
asd_geneSigs[4,] = numbers2colors(as.numeric(asd_geneSigs[4,], blueWhiteRed(100), signed=TRUE, centered = TRUE, lim=c(-1,1)))
asd_geneSigs[5,] = numbers2colors(as.numeric(asd_geneSigs[5,], blueWhiteRed(100), signed=TRUE, centered = TRUE, lim=c(-1,1)))
asd_geneSigs[6,] = numbers2colors(as.numeric(asd_geneSigs[6,], blueWhiteRed(100), signed=TRUE, centered = TRUE, lim=c(-1,1)))
asd_multiExpr[[setIndex]]$geneColors = asd_geneSigs

asd_colors = cbind(labels2colors(asd_merged$colors), t(asd_geneSigs))
asd_labels = c(asd_multiExpr[[setIndex]]$netData$cutParameters, colnames(asd_traitmat))

pdf("./figures/ModuleGeneRelationships_collapseRows_overlap_8851genes.pdf")
plotDendroAndColors(asd_multiExpr[[setIndex]]$netData$dendrograms[[1]], colors=asd_colors, groupLabels=asd_labels, dendroLabels=FALSE, main ="ASD Dendrogram")
plotDendroAndColors(scz_multiExpr[[setIndex]]$netData$TOMdendrogram, colors=scz_colors, groupLabels=scz_labels, dendroLabels=FALSE)
dev.off()
