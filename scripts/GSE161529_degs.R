.libPaths(c("/home/aurora.savino/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
library(GEOquery)
gse <- getGEO("GSE161529",GSEMatrix=F)

gsmids <- Meta(gse)$sample_id
gsmlist <- sapply(gsmids,getGEO)
for(i in 1:length(gsmlist)){
getGEOSuppFiles(gsmids[i])
}

library(edgeR)
mtx<-list.files(pattern="mtx", recursive = T)
barcodes<-list.files(pattern="tsv", recursive = T)
barcodes<-barcodes[-1]
data<-list()
for(i in 1:length(mtx)){ 
  data[[i]]<-edgeR::read10X(mtx=mtx[i], genes="/single cell/GSE161529_features.tsv", barcodes=gsub("-matrix.mtx.gz", "-barcodes.tsv.gz", mtx[i]))
}

for(i in 1:length(mtx)){ 
  colnames(data[[i]])<-paste(mtx[i], 1:ncol(data[[i]]), sep="-")
}


names(data)<-mtx
cancers<-data[c(38:64,66,68,70,72,74)]

rm(data)

cancersdata<-list()
for(i in 1:length(cancers)){ 
  cancersdata[[i]]<-cancers[[i]]$counts
}

rm(cancers)

genes<-read.csv("/home/aurora.savino/crossWGCNA/single cell/GSE161529_features.tsv", sep="\t", header=F)
for(i in 1:length(cancersdata)){ 
rownames(cancersdata[[i]])<-genes[match(rownames(cancersdata[[i]]), genes[,1]),2]
}

library(Seurat)

TNmerge<-cancersdata[[1]]
for(i in 2:8){
  TNmerge<-cbind(TNmerge, cancersdata[[i]])
}
save(TNmerge, file="TNmerge.RData")

HER2merge<-cancersdata[[9]]
for(i in 10:14){
  HER2merge<-cbind(HER2merge, cancersdata[[i]])
}

save(HER2merge, file="HER2merge.RData")

ERmerge<-cancersdata[[15]]
for(i in 16:32){
  ERmerge<-cbind(ERmerge, cancersdata[[i]])
}

save(ERmerge, file="ERmerge.RData")

save(cancersdata, file="cancersdata.RData")


mergedTNBC<-CreateSeuratObject(TNmerge)
mergedTNBC <- NormalizeData(mergedTNBC, verbose = FALSE)
mergedTNBC <- FindVariableFeatures(mergedTNBC, selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
mergedTNBC <- ScaleData(mergedTNBC, verbose = FALSE)
mergedTNBC <- RunPCA(mergedTNBC, npcs = 30, verbose = FALSE)
mergedTNBC <- RunUMAP(mergedTNBC, reduction = "pca", dims = 1:30)
  
FeaturePlot(mergedTNBC, features = c("CD34",  "COL8A1", "NOTCH3", "C1QC", "CD8A",  "KRT17"))
FeaturePlot(mergedTNBC, features = c("PDGFRA", "COL10A1","COL8A1", "KRT17", "KRT5"))
FeaturePlot(mergedTNBC, features = c("CD4", "CD8A","FOXP3", "PDGFRA", "MCAM", "KRT17"))
FeaturePlot(mergedTNBC, features = c("PDGFRA", "KRT5"))

FeaturePlot(mergedTNBC, features = c("NKG7", "CD19","CD68", "CD14", "VWF",  "CD34"))
FeaturePlot(mergedTNBC, features = c("MCF2L", "GREB1","PIP", "TFAP2B", "PTPN4",  "TFF3"))
##epithelial
pdf("TNBC_epi1.pdf", 15,20)
FeaturePlot(mergedTNBC, features = c("EPCAM", "EGFR","CDH1", "KRT14", "ITGA6",  "KRT5"))
dev.off()
pdf("TNBC_epi2.pdf", 15,20)
FeaturePlot(mergedTNBC, features = c("TP63", "KRT17","MME", "KRT8", "KRT18",  "KRT19"))
dev.off()
pdf("TNBC_epi3.pdf", 15,20)
FeaturePlot(mergedTNBC, features = c("FOXA1", "GATA3","MUC1", "CD24", "KIT",  "GABRP"))
dev.off()
##stromal
pdf("TNBC_str1.pdf", 15,20)
FeaturePlot(mergedTNBC, features = c("FAP", "COL1A1","COL3A1", "COL5A1", "ACTA2",  "TAGLN"))
dev.off()
pdf("TNBC_str2.pdf", 15,20)
FeaturePlot(mergedTNBC, features = c("LUM", "FBLN1","COL6A3", "COL1A2", "COL6A1",  "COL6A2"))
dev.off()
##endothelial
pdf("TNBC_endo.pdf", 15,20)
FeaturePlot(mergedTNBC, features = c("PECAM1", "VWF","CDH5", "SELE"))
dev.off()
#immune
pdf("TNBC_imm.pdf", 15,20)
FeaturePlot(mergedTNBC, features = c("LAPTM5", "IL2RG"))
dev.off()
#TCell
pdf("TNBC_Tcell.pdf", 15,20)
FeaturePlot(mergedTNBC, features = c("CD2", "CD3D","CD3E", "CD3G", "CD8A",  "CD8B"))
dev.off()
#BCell
pdf("TNBC_Bcell.pdf", 15,20)
FeaturePlot(mergedTNBC, features = c("MS4A1", "CD79A","CD79B", "BLNK"))
dev.off()
#Macrophage
pdf("TNBC_macro.pdf", 15,20)
FeaturePlot(mergedTNBC, features = c("CD14", "CD68","CD163", "CSF1R"))
dev.off()

mergedTNBC <- FindNeighbors(mergedTNBC, dims = 1:30)
mergedTNBC <- FindClusters(mergedTNBC, resolution = 0.1)

DoHeatmap(mergedTNBC, features = c("PDGFRA", "COL10A1","COL8A1", "KRT17", "KRT5",
                                   "CD4", "CD8A","FOXP3", "PDGFRA", "MCAM","NKG7", "CD19","CD68", "CD14", "VWF",  "CD34"),size = 4,
          angle = 90) 
DimPlot(mergedTNBC, reduction = "umap")

markers <- FindMarkers(object = mergedTNBC, ident.1 = 13)

save(mergedTNBC, file="mergedTNBC.RData")

dataTNBC<-mergedTNBC@assays$RNA

plot<-DimPlot(mergedTNBC, reduction = "umap")

pdf("TNBC_clusts.pdf", 10,10)
LabelClusters(plot=plot, id="ident")
dev.off()

#epi
TNBCepi<-WhichCells(mergedTNBC, ident = c(0,5,3,13,9,4))
#caf
TNBCstr <-  WhichCells(mergedTNBC, ident = 8)
#Tcell
TNBCTcell<-WhichCells(mergedTNBC, ident = c(1,10))
#Bcell
TNBCBcell<-WhichCells(mergedTNBC, ident = 7)
#Endo
TNBCendo<-WhichCells(mergedTNBC, ident = 11)
#Macrophage
TNBCmacro<-WhichCells(mergedTNBC, ident = 2)
#Plasma
TNBCplasma<-WhichCells(mergedTNBC, ident = 6)


save(TNBCstr, file="TNBCstr.RData")
save(TNBCepi, file="TNBCepi.RData")
save(TNBCTcell, file="TNBCTcell.RData")
save(TNBCBcell, file="TNBCBcell.RData")
save(TNBCendo, file="TNBCendo.RData")
save(TNBCmacro, file="TNBCmacro.RData")
save(TNBCplasma, file="TNBCplasma.RData")

mergedHER2<-CreateSeuratObject(HER2merge)

mergedHER2 <- NormalizeData(mergedHER2, verbose = FALSE)
mergedHER2 <- FindVariableFeatures(mergedHER2, selection.method = "vst", 
                                   nfeatures = 2000, verbose = FALSE)
mergedHER2 <- ScaleData(mergedHER2, verbose = FALSE)
mergedHER2 <- RunPCA(mergedHER2, npcs = 30, verbose = FALSE)
mergedHER2 <- RunUMAP(mergedHER2, reduction = "pca", dims = 1:30)

FeaturePlot(mergedHER2, features = c("CD34",  "COL8A1", "NOTCH3", "C1QC", "CD8A",  "KRT17"))
FeaturePlot(mergedHER2, features = c("PDGFRA", "COL10A1","COL8A1", "KRT17", "KRT5"))
FeaturePlot(mergedHER2, features = c("CD4", "CD8A","FOXP3", "PDGFRA", "MCAM", "KRT17"))
FeaturePlot(mergedHER2, features = c("PDGFRA", "KRT5"))

FeaturePlot(mergedHER2, features = c("NKG7", "CD19","CD68", "CD14", "VWF",  "CD34"))
FeaturePlot(mergedHER2, features = c("MCF2L", "GREB1","PIP", "TFAP2B", "PTPN4",  "TFF3"))

mergedHER2 <- FindNeighbors(mergedHER2, dims = 1:30)
mergedHER2 <- FindClusters(mergedHER2, resolution = 0.1)

##epithelial
pdf("HER2_epi1.pdf", 15,20)
FeaturePlot(mergedHER2, features = c("EPCAM", "EGFR","CDH1", "KRT14", "ITGA6",  "KRT5"))
dev.off()
pdf("HER2_epi2.pdf", 15,20)
FeaturePlot(mergedHER2, features = c("TP63", "KRT17","MME", "KRT8", "KRT18",  "KRT19"))
dev.off()
pdf("HER2_epi3.pdf", 15,20)
FeaturePlot(mergedHER2, features = c("FOXA1", "GATA3","MUC1", "CD24", "KIT",  "GABRP"))
dev.off()
##stromal
pdf("HER2_str1.pdf", 15,20)
FeaturePlot(mergedHER2, features = c("FAP", "COL1A1","COL3A1", "COL5A1", "ACTA2",  "TAGLN"))
dev.off()
pdf("HER2_str2.pdf", 15,20)
FeaturePlot(mergedHER2, features = c("LUM", "FBLN1","COL6A3", "COL1A2", "COL6A1",  "COL6A2"))
dev.off()
##endothelial
pdf("HER2_endo.pdf", 15,20)
FeaturePlot(mergedHER2, features = c("PECAM1", "VWF","CDH5", "SELE"))
dev.off()
#immune
pdf("HER2_imm.pdf", 15,20)
FeaturePlot(mergedHER2, features = c("LAPTM5", "IL2RG"))
dev.off()
#TCell
pdf("HER2_Tcell.pdf", 15,20)
FeaturePlot(mergedHER2, features = c("CD2", "CD3D","CD3E", "CD3G", "CD8A",  "CD8B"))
dev.off()
#BCell
pdf("HER2_Bcell.pdf", 15,20)
FeaturePlot(mergedHER2, features = c("MS4A1", "CD79A","CD79B", "BLNK"))
dev.off()
#Macrophage
pdf("HER2_macro.pdf", 15,20)
FeaturePlot(mergedHER2, features = c("CD14", "CD68","CD163", "CSF1R"))
dev.off()

DoHeatmap(mergedHER2, features = c("PDGFRA", "COL10A1","COL8A1", "KRT17", "KRT5",
                                   "CD4", "CD8A","FOXP3", "PDGFRA", "MCAM","NKG7", "CD19","CD68", "CD14", "VWF",  "CD34"),size = 4,
          angle = 90) 
plot<-DimPlot(mergedHER2, reduction = "umap")
pdf("HER2_ctust.pdf",10,10)
LabelClusters(plot=plot, id="ident")
dev.off()

save(mergedHER2, file="mergedHER2.RData")

#epi
HER2epi<-WhichCells(mergedHER2, ident = c(9,5,10,1,3))
#caf
HER2str<-WhichCells(mergedHER2, ident = c(7,13))
#Tcell
HER2Tcell<-WhichCells(mergedHER2, ident = c(0,12))
#Bcell
HER2Bcell<-WhichCells(mergedHER2, ident = c(8))
#Endo
HER2Endo<-WhichCells(mergedHER2, ident = c(11))
#Macrophage
HER2macro<-WhichCells(mergedHER2, ident = c(4))
#Plasma
HER2plasma<-WhichCells(mergedHER2, ident = c(6))


save(HER2epi, file="HER2epi.RData")
save(HER2str, file="HER2str.RData")
save(HER2Tcell, file="HER2Tcell.RData")
save(HER2Bcell, file="HER2Bcell.RData")
save(HER2Endo, file="HER2Endo.RData")
save(HER2macro, file="HER2macro.RData")
save(HER2plasma, file="HER2plasma.RData")

mergedER<-CreateSeuratObject(ERmerge)
mergedER <- NormalizeData(mergedER, verbose = FALSE)
mergedER <- FindVariableFeatures(mergedER, selection.method = "vst", 
                                   nfeatures = 2000, verbose = FALSE)
mergedER <- ScaleData(mergedER, verbose = FALSE)
mergedER <- RunPCA(mergedER, npcs = 30, verbose = FALSE)
mergedER <- RunUMAP(mergedER, reduction = "pca", dims = 1:30)

FeaturePlot(mergedER, features = c("CD34",  "COL8A1", "NOTCH3", "C1QC", "CD8A",  "KRT17"))
FeaturePlot(mergedER, features = c("PDGFRA", "COL10A1","COL8A1", "KRT17", "KRT5"))
FeaturePlot(mergedER, features = c("CD4", "CD8A","FOXP3", "PDGFRA", "MCAM", "KRT17"))
FeaturePlot(mergedER, features = c("PDGFRA", "KRT5"))

FeaturePlot(mergedER, features = c("NKG7", "CD19","CD68", "CD14", "VWF",  "CD34"))

##epithelial
pdf("ER_epi1.pdf", 15,20)
FeaturePlot(mergedER, features = c("EPCAM", "EGFR","CDH1", "KRT14", "ITGA6",  "KRT5"))
dev.off()
pdf("ER_epi2.pdf", 15,20)
FeaturePlot(mergedER, features = c("TP63", "KRT17","MME", "KRT8", "KRT18",  "KRT19"))
dev.off()
pdf("ER_epi3.pdf", 15,20)
FeaturePlot(mergedER, features = c("FOXA1", "GATA3","MUC1", "CD24", "KIT",  "GABRP"))
dev.off()
##stromal
pdf("ER_str1.pdf", 15,20)
FeaturePlot(mergedER, features = c("FAP", "COL1A1","COL3A1", "COL5A1", "ACTA2",  "TAGLN"))
dev.off()
pdf("ER_str2.pdf", 15,20)
FeaturePlot(mergedER, features = c("LUM", "FBLN1","COL6A3", "COL1A2", "COL6A1",  "COL6A2"))
dev.off()
##endothelial
pdf("ER_endo.pdf", 15,20)
FeaturePlot(mergedER, features = c("PECAM1", "VWF","CDH5", "SELE"))
dev.off()
#immune
pdf("ER_imm.pdf", 15,20)
FeaturePlot(mergedER, features = c("LAPTM5", "IL2RG"))
dev.off()
#TCell
pdf("ER_Tcell.pdf", 15,20)
FeaturePlot(mergedER, features = c("CD2", "CD3D","CD3E", "CD3G", "CD8A",  "CD8B"))
dev.off()
#BCell
pdf("ER_Bcell.pdf", 15,20)
FeaturePlot(mergedER, features = c("MS4A1", "CD79A","CD79B", "BLNK"))
dev.off()
#Macrophage
pdf("ER_macro.pdf", 15,20)
FeaturePlot(mergedER, features = c("CD14", "CD68","CD163", "CSF1R"))
dev.off()

mergedER <- FindNeighbors(mergedER, dims = 1:30)
mergedER <- FindClusters(mergedER, resolution = 0.1)



DoHeatmap(mergedER, features = c("PDGFRA", "COL10A1","COL8A1", "KRT17", "KRT5",
                                   "CD4", "CD8A","FOXP3", "PDGFRA", "MCAM","NKG7", "CD19","CD68", "CD14", "VWF",  "CD34"),size = 4,
          angle = 90) 
plot<-DimPlot(mergedER, reduction = "umap")
pdf("ER_clust.pdf", 10,10)
LabelClusters(plot=plot, id="ident")
dev.off()

#epi
ERepi<-WhichCells(mergedER, ident = c(0:20)[-c(1,2,4,17,18,15,16)])
#caf
ERstr <- WhichCells(mergedER, ident = 4)
#Tcell
ERTcell <- WhichCells(mergedER, ident = 1)
#Bcell
ERBcell <- WhichCells(mergedER, ident = 18)
#Endo
ERendo <- WhichCells(mergedER, ident = 17)
#Macrophage
ERmacro <- WhichCells(mergedER, ident = 2)

save(ERepi, file="ERepi.RData")
save(ERstr, file="ERstr.RData")
save(ERTcell, file="ERTcell.RData")
save(ERBcell, file="ERBcell.RData")
save(ERmacro, file="ERmacro.RData")
save(ERendo, file="ERendo.RData")

save(mergedER, file="mergedER.RData")


#######################################
#####Versione con pi? cluster epi e stroma
###############################################


###selezione delle sole cellule epiteliali e CAF
load("/home/aurora.savino/crossWGCNA/cancersdata.RData")
load("/home/aurora.savino/crossWGCNA/TNBCepi.RData")
load("/home/aurora.savino/crossWGCNA/HER2epi.RData")
load("/home/aurora.savino/crossWGCNA/ERepi.RData")
load("/home/aurora.savino/crossWGCNA/TNBCstr.RData")
load("/home/aurora.savino/crossWGCNA/HER2str.RData")
load("/home/aurora.savino/crossWGCNA/ERstr.RData")


epithelial<-c(TNBCepi, HER2epi, ERepi)
cafs<-c(TNBCstr, HER2str, ERstr)

load("/home/aurora.savino/crossWGCNA/mergedTNBC.RData")

dataTNBCepi <- as.matrix(GetAssayData(mergedTNBC[,TNBCepi], slot = "data", assay="RNA" ))
dataTNBCcafs <- as.matrix(GetAssayData(mergedTNBC[,TNBCstr], slot = "data", assay="RNA" ))

averagedEpiTNBC<-matrix(nrow=nrow(dataTNBCepi), ncol=7)
averagedStromaTNBC<-matrix(nrow=nrow(dataTNBCcafs), ncol=7)
j<-1
for(i in 2:8){
  averagedEpiTNBC[,j]<-rowMeans(dataTNBCepi[,intersect(colnames(cancersdata[[i]]), TNBCepi)], na.rm=T)
  averagedStromaTNBC[,j]<-rowMeans(dataTNBCcafs[,intersect(colnames(cancersdata[[i]]), TNBCstr)], na.rm=T)
  j<-1+j
}

save(averagedEpiTNBC,  file="averagedEpiTNBC.RData")
save(averagedStromaTNBC,  file="averagedStromaTNBC.RData")

load("/home/aurora.savino/crossWGCNA/mergedHER2.RData")

dataHER2epi <- as.matrix(GetAssayData(mergedHER2[,HER2epi], slot = "data", assay="RNA" ))
dataHER2cafs <- as.matrix(GetAssayData(mergedHER2[,HER2str], slot = "data", assay="RNA" ))

averagedEpiHER2<-matrix(nrow=nrow(dataHER2epi), ncol=5)
averagedStromaHER2<-matrix(nrow=nrow(dataHER2cafs), ncol=5)
j<-1
for(i in 10:14){
  averagedEpiHER2[,j]<-rowMeans(dataHER2epi[,intersect(colnames(cancersdata[[i]]), HER2epi)], na.rm=T)
  averagedStromaHER2[,j]<-rowMeans(dataHER2cafs[,intersect(colnames(cancersdata[[i]]), HER2str)], na.rm=T)
j<-j+1
  }

save(averagedEpiHER2,  file="averagedEpiHER2.RData")
save(averagedStromaHER2,  file="averagedStromaHER2.RData")

load("/home/aurora.savino/crossWGCNA/mergedER.RData")

dataERepi1 <- as.matrix(GetAssayData(mergedER[,ERepi[1:30000]], slot = "data", assay="RNA" ))
dataERepi2 <- as.matrix(GetAssayData(mergedER[,ERepi[-c(1:30000)]], slot = "data", assay="RNA" ))
dataERcafs <- as.matrix(GetAssayData(mergedER[,ERstr], slot = "data", assay="RNA" ))
dataERepi<-cbind(dataERepi1, dataERepi2)
rm(dataERepi1)
rm(dataERepi2)
gc()

averagedEpiER<-matrix(nrow=nrow(dataERepi), ncol=17)
averagedStromaER<-matrix(nrow=nrow(dataERcafs), ncol=17)
j<-1
for(i in 16:32){
  averagedEpiER[,j]<-rowMeans(dataERepi[,intersect(colnames(cancersdata[[i]]), ERepi)], na.rm=T)
  averagedStromaER[,j]<-rowMeans(dataERcafs[,intersect(colnames(cancersdata[[i]]), ERstr)], na.rm=T)
  j<-j+1
}

save(averagedEpiER,  file="averagedEpiER.RData")
save(averagedStromaER,  file="averagedStromaER.RData")


################################
######Double WGCNA epi stroma tutte le cellule
##############################
load("/home/aurora.savino/crossWGCNA/averagedEpiTNBC.RData")
load("/home/aurora.savino/crossWGCNA/averagedEpiHER2.RData")
load("/home/aurora.savino/crossWGCNA/averagedEpiER.RData")
load("/home/aurora.savino/crossWGCNA/averagedStromaTNBC.RData")
load("/home/aurora.savino/crossWGCNA/averagedStromaHER2.RData")
load("/home/aurora.savino/crossWGCNA/averagedStromaER.RData")
load("~/crossWGCNA/TNmerge.RData")

library(WGCNA)
library(xlsx)
library(limma)
library(pheatmap)
library(Hmisc)
averagedEpi<-cbind(averagedEpiTNBC, averagedEpiHER2, averagedEpiER)
averagedStroma<-cbind(averagedStromaTNBC, averagedStromaHER2, averagedStromaER)

rownames(averagedEpi)<-rownames(TNmerge)
rownames(averagedStroma)<-rownames(TNmerge)

###filtering data
ind<-which(colSums(is.na(averagedEpi))==0 & colSums(is.na(averagedStroma))==0)

averagedEpi<-averagedEpi[, ind]
averagedStroma<-averagedStroma[, ind]

genes<-unique(c(which(apply(averagedStroma,1,var)>=quantile(apply(averagedStroma,1,var),0.5)),
                which(apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5))))
stroma<-averagedStroma[genes,]
epi<-averagedEpi[genes,]


#merging stroma and epi
rownames(stroma)<-paste(rownames(stroma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(stroma)

data_merged<-rbind(stroma, epi)
var_merged<-apply(data_merged, 1, var)
data_merged<-data_merged[which(var_merged!=0),]

#Functions da A_L.R
degsc<-network(data=data_merged, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(degsc, file="degsc.RData")


################################
######Double WGCNA epi - Tcell tutte le cellule
##############################
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedTcellTNBC.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedTcellHER2.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedTcellER.RData")
load("/home/aurora.savino/crossWGCNA/genenames.RData")

library(WGCNA)
library(xlsx)
library(limma)
library(pheatmap)
library(Hmisc)
averagedEpi<-cbind(averagedEpiTNBC, averagedEpiHER2, averagedEpiER)
averagedTcell<-cbind(averagedTcellTNBC, averagedTcellHER2, averagedTcellER)

rownames(averagedEpi)<-genenames
rownames(averagedTcell)<-genenames

###filtering data
ind<-which(colSums(is.na(averagedEpi))==0 & colSums(is.na(averagedTcell))==0)

averagedEpi<-averagedEpi[, ind]
averagedTcell<-averagedTcell[, ind]

genes<-unique(c(which(apply(averagedTcell,1,var)>=quantile(apply(averagedTcell,1,var),0.5)),
                which(apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5))))
Tcell<-averagedTcell[genes,]
epi<-averagedEpi[genes,]

#merging stroma and epi
rownames(Tcell)<-paste(rownames(Tcell), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(Tcell)

data_merged<-rbind(Tcell, epi)

#Functions da A_L.R
degsc_ET<-network(data=data_merged, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(degsc_ET, file="degsc_ET.RData")




################################
######Double WGCNA epi - Bcell tutte le cellule
##############################
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedBcellTNBC.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedBcellHER2.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedBcellER.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/genenames.RData")

library(WGCNA)
library(xlsx)
library(limma)
library(pheatmap)
library(Hmisc)
averagedEpi<-cbind(averagedEpiTNBC, averagedEpiHER2, averagedEpiER)
averagedBcell<-cbind(averagedBcellTNBC, averagedBcellHER2, averagedBcellER)

rownames(averagedEpi)<-genenames
rownames(averagedBcell)<-genenames

###filtering data
ind<-which(colSums(is.na(averagedEpi))==0 & colSums(is.na(averagedBcell))==0)

averagedEpi<-averagedEpi[, ind]
averagedBcell<-averagedBcell[, ind]

genes<-unique(c(which(apply(averagedBcell,1,var)>=quantile(apply(averagedBcell,1,var),0.5)),
                which(apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5))))
Bcell<-averagedBcell[genes,]
epi<-averagedEpi[genes,]

#merging stroma and epi
rownames(Bcell)<-paste(rownames(Bcell), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(Bcell)

data_mergedB<-rbind(Bcell, epi)

#Functions da A_L.R
degsc_EB<-network(data=data_mergedB, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(degsc_EB, file="degsc_EB.RData")

################################
######Double WGCNA epi - Macrophage tutte le cellule
##############################
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedmacroTNBC.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedmacroHER2.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedmacroER.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/genenames.RData")

library(WGCNA)
library(xlsx)
library(limma)
library(pheatmap)
library(Hmisc)
averagedEpi<-cbind(averagedEpiTNBC, averagedEpiHER2, averagedEpiER)
averagedmacro<-cbind(averagedmacroTNBC, averagedmacroHER2, averagedmacroER)

rownames(averagedEpi)<-genenames
rownames(averagedmacro)<-genenames

###filtering data
ind<-which(colSums(is.na(averagedEpi))==0 & colSums(is.na(averagedmacro))==0)

averagedEpi<-averagedEpi[, ind]
averagedmacro<-averagedmacro[, ind]

genes<-unique(c(which(apply(averagedmacro,1,var)>=quantile(apply(averagedmacro,1,var),0.5)),
                which(apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5))))
macro<-averagedmacro[genes,]
epi<-averagedEpi[genes,]


#merging stroma and epi
rownames(macro)<-paste(rownames(macro), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(macro)

data_mergedM<-rbind(macro, epi)

#Functions da A_L.R
degsc_EM<-network(data=data_mergedM, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(degsc_EM, file="degsc_EM.RData")

################################
######Double WGCNA epi - Macrophage tutte le cellule
##############################
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedendoTNBC.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedendoHER2.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedendoER.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/genenames.RData")

library(WGCNA)
library(xlsx)
library(limma)
library(pheatmap)
library(Hmisc)
averagedEpi<-cbind(averagedEpiTNBC, averagedEpiHER2, averagedEpiER)
averagedendo<-cbind(averagedendoTNBC, averagedendoHER2, averagedendoER)

rownames(averagedEpi)<-genenames
rownames(averagedendo)<-genenames

###filtering data
ind<-which(colSums(is.na(averagedEpi))==0 & colSums(is.na(averagedendo))==0)

averagedEpi<-averagedEpi[, ind]
averagedendo<-averagedendo[, ind]

genes<-unique(c(which(apply(averagedendo,1,var)>=quantile(apply(averagedendo,1,var),0.5)),
                which(apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5))))
endo<-averagedendo[genes,]
epi<-averagedEpi[genes,]

#merging stroma and epi
rownames(endo)<-paste(rownames(endo), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(endo)

data_mergedE<-rbind(endo, epi)

#Functions da A_L.R
degsc_EE<-network(data=data_mergedE, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(degsc_EE, file="degsc_EE.RData")

################################
######Double WGCNA epi - Macrophage tutte le cellule
##############################
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedplasmaTNBC.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedplasmaHER2.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedplasmaER.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/genenames.RData")

library(WGCNA)
library(xlsx)
library(limma)
library(pheatmap)
library(Hmisc)
averagedEpi<-cbind(averagedEpiTNBC, averagedEpiHER2, averagedEpiER)
averagedplasma<-cbind(averagedplasmaTNBC, averagedplasmaHER2, averagedplasmaER)

rownames(averagedEpi)<-genenames
rownames(averagedplasma)<-genenames

###filtering data
ind<-which(colSums(is.na(averagedEpi))==0 & colSums(is.na(averagedplasma))==0)

averagedEpi<-averagedEpi[, ind]
averagedplasma<-averagedplasma[, ind]

genes<-unique(c(which(apply(averagedplasma,1,var)>=quantile(apply(averagedplasma,1,var),0.5)),
which(apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5))))
plasma<-averagedplasma[genes,]
epi<-averagedEpi[genes,]

#merging stroma and epi
rownames(plasma)<-paste(rownames(plasma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(plasma)

data_mergedP<-rbind(plasma, epi)

#Functions da A_L.R
degsc_EP<-network(data=data_mergedP, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(degsc_EP, file="degsc_EP.RData")

