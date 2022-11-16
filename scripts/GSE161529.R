library(GEOquery)
gse <- getGEO("GSE161529",GSEMatrix=F)

gsmids <- Meta(gse)$sample_id
gsmlist <- sapply(gsmids[1:5],getGEO)
getGEOSuppFiles(gsmids[1])

library(edgeR)
mtx<-list.files(pattern="mtx")
barcodes<-list.files(pattern="tsv")
barcodes<-barcodes[-1]
data<-list()
for(i in 1:length(mtx)){ 
  data[[i]]<-edgeR::read10X(mtx=mtx[i], genes="features.tsv", barcodes=gsub("-matrix.mtx.gz", "-barcodes.tsv.gz", mtx[i]))
  
}

for(i in 1:length(mtx)){ 
  colnames(data[[i]])<-paste(mtx[i], 1:ncol(data[[i]]), sep="-")
  
}


names(data)<-mtx
cancers<-data[c(29:55,57,59,61,63,65,68)]

rm(data)

cancersdata<-list()
for(i in 1:length(cancers)){ 
  cancersdata[[i]]<-cancers[[i]]$counts
}

rm(cancers)

genes<-read.csv("features.tsv", sep="\t", header=F)
for(i in 1:length(cancersdata)){ 
rownames(cancersdata[[i]])<-genes[match(rownames(cancersdata[[i]]), genes[,1]),2]
}

library(Seurat)
seuratobj<-list()
for(i in 1:length(cancersdata)){ 
seuratobj[[i]] <- CreateSeuratObject(cancersdata[[i]]) 
}

for (i in 1:length(seuratobj)) {
  seuratobj[[i]] <- NormalizeData(seuratobj[[i]], verbose = FALSE)
  seuratobj[[i]] <- FindVariableFeatures(seuratobj[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

for (i in 1:length(seuratobj)) {
  seuratobj[[i]] <- ScaleData(seuratobj[[i]], verbose = FALSE)
  seuratobj[[i]] <- RunPCA(seuratobj[[i]], npcs = 30, verbose = FALSE)
  seuratobj[[i]] <- RunUMAP(seuratobj[[i]], reduction = "pca", dims = 1:30)
}

FeaturePlot(seuratobj[[2]], features = c("CD34",  "COL8A1", "NOTCH3", "C1QC", "CD8A",  "KRT17"))
FeaturePlot(seuratobj[[4]], features = c("PDGFRA", "COL10A1","COL8A1", "KRT17"))
FeaturePlot(seuratobj[[4]], features = c("PDGFRA", "COL10A1","COL8A1", "KRT17", "KRT5"))


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
for(i in 16:33){
  ERmerge<-cbind(ERmerge, cancersdata[[i]])
}

save(ERmerge, file="ERmerge.RData")

save(seuratobj, file="seuratobj.RData")
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
FeaturePlot(mergedTNBC, features = c("EPCAM", "EGFR","CDH1", "KRT14", "ITGA6",  "KRT5"))
FeaturePlot(mergedTNBC, features = c("TP63", "KRT17","MME", "KRT8", "KRT18",  "KRT19"))
FeaturePlot(mergedTNBC, features = c("FOXA1", "GATA3","MUC1", "CD24", "KIT",  "GABRP"))
##stromal
FeaturePlot(mergedTNBC, features = c("FAP", "COL1A1","COL3A1", "COL5A1", "ACTA2",  "TAGLN"))
FeaturePlot(mergedTNBC, features = c("LUM", "FBLN1","COL6A3", "COL1A2", "COL6A1",  "COL6A2"))
##endothelial
FeaturePlot(mergedTNBC, features = c("PECAM1", "VWF","CDH5", "SELE"))
#immune
FeaturePlot(mergedTNBC, features = c("LAPTM5", "IL2RG"))
#TCell
FeaturePlot(mergedTNBC, features = c("CD2", "CD3D","CD3E", "CD3G", "CD8A",  "CD8B"))
#BCell
FeaturePlot(mergedTNBC, features = c("MS4A1", "CD79A","CD79B", "BLNK"))
#Macrophage
FeaturePlot(mergedTNBC, features = c("CD14", "CD68","CD163", "CSF1R"))


mergedTNBC <- FindNeighbors(mergedTNBC, dims = 1:30)
mergedTNBC <- FindClusters(mergedTNBC, resolution = 0.1)

DoHeatmap(mergedTNBC, features = c("PDGFRA", "COL10A1","COL8A1", "KRT17", "KRT5",
                                   "CD4", "CD8A","FOXP3", "PDGFRA", "MCAM","NKG7", "CD19","CD68", "CD14", "VWF",  "CD34"),size = 4,
          angle = 90) 
DimPlot(mergedTNBC, reduction = "umap")

markers <- FindMarkers(object = mergedTNBC, ident.1 = 13)

save(mergedTNBC, file="mergedTNBC.RData")

dataTNBC<-mergedTNBC@assays$RNA$
mergedTNBC@meta.data[which(mergedTNBC@meta.data$seurat_clusters==0),1]

plot<-DimPlot(mergedTNBC, reduction = "umap")
LabelClusters(plot=plot, id="ident")

#epi
TNBC0 <- WhichCells(mergedTNBC, ident = 0)
TNBC5 <- WhichCells(mergedTNBC, ident = 5)
TNBCepiall<-WhichCells(mergedTNBC, ident = c(0,5,3,13,9,4))
#caf
TNBC8 <-  WhichCells(mergedTNBC, ident = 8)
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


save(TNBC0, file="TNBC0.RData")
save(TNBC5, file="TNBC5.RData")
save(TNBC8, file="TNBC8.RData")

save(TNBCepiall, file="TNBCepiall.RData")
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
FeaturePlot(mergedHER2, features = c("EPCAM", "EGFR","CDH1", "KRT14", "ITGA6",  "KRT5"))
FeaturePlot(mergedHER2, features = c("TP63", "KRT17","MME", "KRT8", "KRT18",  "KRT19"))
FeaturePlot(mergedHER2, features = c("FOXA1", "GATA3","MUC1", "CD24", "KIT",  "GABRP"))
##stromal
FeaturePlot(mergedHER2, features = c("FAP", "COL1A1","COL3A1", "COL5A1", "ACTA2",  "TAGLN"))
FeaturePlot(mergedHER2, features = c("LUM", "FBLN1","COL6A3", "COL1A2", "COL6A1",  "COL6A2"))
##endothelial
FeaturePlot(mergedHER2, features = c("PECAM1", "VWF","CDH5", "SELE"))
#immune
FeaturePlot(mergedHER2, features = c("LAPTM5", "IL2RG"))
#TCell
FeaturePlot(mergedHER2, features = c("CD2", "CD3D","CD3E", "CD3G", "CD8A",  "CD8B"))
#BCell
FeaturePlot(mergedHER2, features = c("MS4A1", "CD79A","CD79B", "BLNK"))
#Macrophage
FeaturePlot(mergedHER2, features = c("CD14", "CD68","CD163", "CSF1R"))


DoHeatmap(mergedHER2, features = c("PDGFRA", "COL10A1","COL8A1", "KRT17", "KRT5",
                                   "CD4", "CD8A","FOXP3", "PDGFRA", "MCAM","NKG7", "CD19","CD68", "CD14", "VWF",  "CD34"),size = 4,
          angle = 90) 
plot<-DimPlot(mergedHER2, reduction = "umap")
LabelClusters(plot=plot, id="ident")

save(mergedHER2, file="mergedHER2.RData")

#epi
dataHER210 <- as.matrix(GetAssayData(mergedHER2, slot = "scale.data" )[, WhichCells(mergedHER2, ident = 10)])
#caf
dataHER27 <- as.matrix(GetAssayData(mergedHER2, slot = "scale.data" )[, WhichCells(mergedHER2, ident = 7)])

#epi
HER210 <- WhichCells(mergedHER2, ident = 10)
HER2epiall<-WhichCells(mergedHER2, ident = c(9,5,10,1,3))
#caf
HER27 <- WhichCells(mergedHER2, ident = 7)
HER2CAFall<-WhichCells(mergedHER2, ident = c(7,13))
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


save(HER210, file="HER210.RData")
save(HER27, file="HER27.RData")

save(HER2epiall, file="HER2epiall.RData")
save(HER2CAFall, file="HER2CAFall.RData")
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
FeaturePlot(mergedER, features = c("EPCAM", "EGFR","CDH1", "KRT14", "ITGA6",  "KRT5"))
FeaturePlot(mergedER, features = c("TP63", "KRT17","MME", "KRT8", "KRT18",  "KRT19"))
FeaturePlot(mergedER, features = c("FOXA1", "GATA3","MUC1", "CD24", "KIT",  "GABRP"))
##stromal
FeaturePlot(mergedER, features = c("FAP", "COL1A1","COL3A1", "COL5A1", "ACTA2",  "TAGLN"))
FeaturePlot(mergedER, features = c("LUM", "FBLN1","COL6A3", "COL1A2", "COL6A1",  "COL6A2"))
##endothelial
FeaturePlot(mergedER, features = c("PECAM1", "VWF","CDH5", "SELE"))
#immune
FeaturePlot(mergedER, features = c("LAPTM5", "IL2RG"))
#TCell
FeaturePlot(mergedER, features = c("CD2", "CD3D","CD3E", "CD3G", "CD8A",  "CD8B"))
#BCell
FeaturePlot(mergedER, features = c("MS4A1", "CD79A","CD79B", "BLNK"))
#Macrophage
FeaturePlot(mergedER, features = c("CD14", "CD68","CD163", "CSF1R"))

mergedER <- FindNeighbors(mergedER, dims = 1:30)
mergedER <- FindClusters(mergedER, resolution = 0.1)



DoHeatmap(mergedER, features = c("PDGFRA", "COL10A1","COL8A1", "KRT17", "KRT5",
                                   "CD4", "CD8A","FOXP3", "PDGFRA", "MCAM","NKG7", "CD19","CD68", "CD14", "VWF",  "CD34"),size = 4,
          angle = 90) 
plot<-DimPlot(mergedER, reduction = "umap")
LabelClusters(plot=plot, id="ident")

#epi
ER10 <- WhichCells(mergedER, ident = 10)
EREpiall<-WhichCells(mergedER, ident = c(1:21)[-c(2,19,17,16,4,18)])
#caf
ER4 <- WhichCells(mergedER, ident = 4)
ERCAFall<-WhichCells(mergedER, ident = c(4,16))
#Tcell
ER0 <- WhichCells(mergedER, ident = 0)
#Bcell
ER19 <- WhichCells(mergedER, ident = 19)
#Endo
ER18 <- WhichCells(mergedER, ident = 18)
#Macrophage
ER2 <- WhichCells(mergedER, ident = 2)
#Plasma
ER17 <- WhichCells(mergedER, ident = 17)

save(ER10, file="ER10.RData")
save(ER4, file="ER4.RData")
save(EREpiall, file="EREpiall.RData")
save(ERCAFall, file="ERCAFall.RData")
save(ER0, file="ER0.RData")
save(ER19, file="ER19.RData")
save(ER18, file="ER18.RData")
save(ER2, file="ER2.RData")
save(ER17, file="ER17.RData")

save(mergedER, file="mergedER.RData")

###primissima versione
###selezione delle sole cellule epiteliali e CAF
epithelial<-c(TNBC0, TNBC5, HER210, ER10)
cafs<-c(TNBC8, HER27, ER4)

cancersdataEpi<-list()
for(i in 1:length(cancersdata)){ 
  cancersdataEpi[[i]]<-cancersdata[[i]][,which(colnames(cancersdata[[i]]) %in% c(epithelial))]
}

cancersdataStroma<-list()
for(i in 1:length(cancersdata)){ 
  cancersdataStroma[[i]]<-cancersdata[[i]][,which(colnames(cancersdata[[i]]) %in% c(cafs))]
}


datasel<-cbind(cancersdataEpi[[1]], cancersdataStroma[[1]])
for(i in 2:length(cancersdata)){ 
  datasel<-cbind(datasel, cancersdataEpi[[i]], cancersdataStroma[[i]])
}

#normalizzazione con  Seurat (tutte insieme)
mergedAll<-CreateSeuratObject(datasel)
mergedAll <- NormalizeData(mergedAll, verbose = FALSE)
#mergedAll <- FindVariableFeatures(mergedAll, selection.method = "vst", 
#                                 nfeatures = 2000, verbose = FALSE)
mergedAll <- ScaleData(mergedAll, verbose = FALSE)
#mergedAll <- RunPCA(mergedAll, npcs = 30, verbose = FALSE)
#mergedAll <- RunUMAP(mergedAll, reduction = "pca", dims = 1:30)

FeaturePlot(mergedAll, features = c("PDGFRA", "KRT5"))

#di ogni paziente, valore medio
datasel <- as.matrix(GetAssayData(mergedAll, slot = "counts" ))

averagedEpi<-matrix(nrow=nrow(datasel), ncol=33)
averagedStroma<-matrix(nrow=nrow(datasel), ncol=33)

for(i in 1:length(cancersdata)){
  averagedEpi[,i]<-rowMeans(datasel[,colnames(cancersdataEpi[[i]])], na.rm=T)
  averagedStroma[,i]<-rowMeans(datasel[,colnames(cancersdataStroma[[i]])], na.rm=T)
  }

save(averagedEpi,  file="averagedEpi.RData")
save(averagedStroma,  file="averagedStroma.RData")

#######################################
#####Versione con più cluster epi e stroma
###############################################


###selezione delle sole cellule epiteliali e CAF
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/cancersdata.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/TNBCepiall.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/HER2epiall.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/EREpiall.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/TNBC8.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/HER2CAFall.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/ERCAFall.RData")


epithelial<-c(TNBCepiall, HER2epiall, EREpiall)
cafs<-c(TNBC8, HER2CAFall, ERCAFall)

load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/mergedTNBC.RData")

dataTNBCepi <- as.matrix(GetAssayData(mergedTNBC[,TNBCepiall], slot = "data", assay="RNA" ))
dataTNBCcafs <- as.matrix(GetAssayData(mergedTNBC[,TNBC8], slot = "data", assay="RNA" ))

averagedEpiTNBC<-matrix(nrow=nrow(dataTNBCepi), ncol=8)
averagedStromaTNBC<-matrix(nrow=nrow(dataTNBCcafs), ncol=8)
for(i in 1:8){
  averagedEpiTNBC[,i]<-rowMeans(dataTNBCepi[,intersect(colnames(cancersdata[[i]]), TNBCepiall)], na.rm=T)
  averagedStromaTNBC[,i]<-rowMeans(dataTNBCcafs[,intersect(colnames(cancersdata[[i]]), TNBC8)], na.rm=T)
}

save(averagedEpiTNBC,  file="averagedEpiTNBC.RData")
save(averagedStromaTNBC,  file="averagedStromaTNBC.RData")

load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/mergedHER2.RData")

dataHER2epi <- as.matrix(GetAssayData(mergedHER2[,HER2epiall], slot = "data", assay="RNA" ))
dataHER2cafs <- as.matrix(GetAssayData(mergedHER2[,HER2CAFall], slot = "data", assay="RNA" ))

averagedEpiHER2<-matrix(nrow=nrow(dataHER2epi), ncol=6)
averagedStromaHER2<-matrix(nrow=nrow(dataHER2cafs), ncol=6)
j<-1
for(i in 9:14){
  averagedEpiHER2[,j]<-rowMeans(dataHER2epi[,intersect(colnames(cancersdata[[i]]), HER2epiall)], na.rm=T)
  averagedStromaHER2[,j]<-rowMeans(dataHER2cafs[,intersect(colnames(cancersdata[[i]]), HER2CAFall)], na.rm=T)
j<-j+1
  }

save(averagedEpiHER2,  file="averagedEpiHER2.RData")
save(averagedStromaHER2,  file="averagedStromaHER2.RData")

load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/mergedER.RData")

dataERepi1 <- as.matrix(GetAssayData(mergedER[,EREpiall[1:30000]], slot = "data", assay="RNA" ))
dataERepi2 <- as.matrix(GetAssayData(mergedER[,EREpiall[-c(1:30000)]], slot = "data", assay="RNA" ))
dataERcafs <- as.matrix(GetAssayData(mergedER[,ERCAFall], slot = "data", assay="RNA" ))
dataERepi<-cbind(dataERepi1, dataERepi2)
rm(dataERepi1)
rm(dataERepi2)
gc()

averagedEpiER<-matrix(nrow=nrow(dataERepi), ncol=19)
averagedStromaER<-matrix(nrow=nrow(dataERcafs), ncol=19)
j<-1
for(i in 15:33){
  averagedEpiER[,j]<-rowMeans(dataERepi[,intersect(colnames(cancersdata[[i]]), EREpiall)], na.rm=T)
  averagedStromaER[,j]<-rowMeans(dataERcafs[,intersect(colnames(cancersdata[[i]]), ERCAFall)], na.rm=T)
  j<-j+1
}

save(averagedEpiER,  file="averagedEpiER.RData")
save(averagedStromaER,  file="averagedStromaER.RData")


###selezione delle sole cellule T, B, macrophages, endothelial, plasma
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/cancersdata.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/TNBCTcell.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/TNBCBcell.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/TNBCendo.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/TNBCmacro.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/TNBCplasma.RData")


load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/mergedTNBC.RData")

dataTNBCTcell <- as.matrix(GetAssayData(mergedTNBC[,TNBCTcell], slot = "data", assay="RNA" ))
dataTNBCBcell <- as.matrix(GetAssayData(mergedTNBC[,TNBCBcell], slot = "data", assay="RNA" ))
dataTNBCendo <- as.matrix(GetAssayData(mergedTNBC[,TNBCendo], slot = "data", assay="RNA" ))
dataTNBCmacro <- as.matrix(GetAssayData(mergedTNBC[,TNBCmacro], slot = "data", assay="RNA" ))
dataTNBCplasma <- as.matrix(GetAssayData(mergedTNBC[,TNBCplasma], slot = "data", assay="RNA" ))

averagedTcellTNBC<-matrix(nrow=nrow(dataTNBCTcell), ncol=8)
averagedBcellTNBC<-matrix(nrow=nrow(dataTNBCBcell), ncol=8)
averagedendoTNBC<-matrix(nrow=nrow(dataTNBCendo), ncol=8)
averagedmacroTNBC<-matrix(nrow=nrow(dataTNBCmacro), ncol=8)
averagedplasmaTNBC<-matrix(nrow=nrow(dataTNBCplasma), ncol=8)

for(i in 1:8){
  averagedTcellTNBC[,i]<-rowMeans(dataTNBCTcell[,intersect(colnames(cancersdata[[i]]), TNBCTcell)], na.rm=T)
  averagedBcellTNBC[,i]<-rowMeans(dataTNBCBcell[,intersect(colnames(cancersdata[[i]]), TNBCBcell)], na.rm=T)
  averagedendoTNBC[,i]<-rowMeans(dataTNBCendo[,intersect(colnames(cancersdata[[i]]), TNBCendo)], na.rm=T)
  averagedmacroTNBC[,i]<-rowMeans(dataTNBCmacro[,intersect(colnames(cancersdata[[i]]), TNBCmacro)], na.rm=T)
  if(length(intersect(colnames(cancersdata[[i]]), TNBCplasma))>1){
  averagedplasmaTNBC[,i]<-rowMeans(dataTNBCplasma[,intersect(colnames(cancersdata[[i]]), TNBCplasma)], na.rm=T)
}
  }

save(averagedTcellTNBC,  file="averagedTcellTNBC.RData")
save(averagedBcellTNBC,  file="averagedBcellTNBC.RData")
save(averagedendoTNBC,  file="averagedendoTNBC.RData")
save(averagedmacroTNBC,  file="averagedmacroTNBC.RData")
save(averagedplasmaTNBC,  file="averagedplasmaTNBC.RData")

##HER2

load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/mergedHER2.RData")

load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/HER2Tcell.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/HER2Bcell.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/HER2Endo.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/HER2macro.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/HER2plasma.RData")


dataHER2Tcell <- as.matrix(GetAssayData(mergedHER2[,HER2Tcell], slot = "data", assay="RNA" ))
dataHER2Bcell <- as.matrix(GetAssayData(mergedHER2[,HER2Bcell], slot = "data", assay="RNA" ))
dataHER2endo <- as.matrix(GetAssayData(mergedHER2[,HER2Endo], slot = "data", assay="RNA" ))
dataHER2macro <- as.matrix(GetAssayData(mergedHER2[,HER2macro], slot = "data", assay="RNA" ))
dataHER2plasma <- as.matrix(GetAssayData(mergedHER2[,HER2plasma], slot = "data", assay="RNA" ))

averagedTcellHER2<-matrix(nrow=nrow(dataHER2Tcell), ncol=6)
averagedBcellHER2<-matrix(nrow=nrow(dataHER2Bcell), ncol=6)
averagedendoHER2<-matrix(nrow=nrow(dataHER2endo), ncol=6)
averagedmacroHER2<-matrix(nrow=nrow(dataHER2macro), ncol=6)
averagedplasmaHER2<-matrix(nrow=nrow(dataHER2plasma), ncol=6)
j<-1
for(i in 9:14){
  averagedTcellHER2[,j]<-rowMeans(dataHER2Tcell[,intersect(colnames(cancersdata[[i]]), HER2Tcell)], na.rm=T)
  averagedBcellHER2[,j]<-rowMeans(dataHER2Bcell[,intersect(colnames(cancersdata[[i]]), HER2Bcell)], na.rm=T)
  averagedendoHER2[,j]<-rowMeans(dataHER2endo[,intersect(colnames(cancersdata[[i]]), HER2Endo)], na.rm=T)
  averagedmacroHER2[,j]<-rowMeans(dataHER2macro[,intersect(colnames(cancersdata[[i]]), HER2macro)], na.rm=T)
  averagedplasmaHER2[,j]<-rowMeans(dataHER2plasma[,intersect(colnames(cancersdata[[i]]), HER2plasma)], na.rm=T)
  j<-j+1
  }

save(averagedTcellHER2,  file="averagedTcellHER2.RData")
save(averagedBcellHER2,  file="averagedBcellHER2.RData")
save(averagedendoHER2,  file="averagedendoHER2.RData")
save(averagedmacroHER2,  file="averagedmacroHER2.RData")
save(averagedplasmaHER2,  file="averagedplasmaHER2.RData")

####ER

load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/mergedER.RData")

load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/ER0.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/ER19.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/ER18.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/ER2.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/ER17.RData")

load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/mergedTNBC.RData")

dataERTcell <- as.matrix(GetAssayData(mergedER[,ER0], slot = "data", assay="RNA" ))
dataERBcell <- as.matrix(GetAssayData(mergedER[,ER19], slot = "data", assay="RNA" ))
dataERendo <- as.matrix(GetAssayData(mergedER[,ER18], slot = "data", assay="RNA" ))
dataERmacro <- as.matrix(GetAssayData(mergedER[,ER2], slot = "data", assay="RNA" ))
dataERplasma <- as.matrix(GetAssayData(mergedER[,ER17], slot = "data", assay="RNA" ))

averagedTcellER<-matrix(nrow=nrow(dataERTcell), ncol=19)
averagedBcellER<-matrix(nrow=nrow(dataERBcell), ncol=19)
averagedendoER<-matrix(nrow=nrow(dataERendo), ncol=19)
averagedmacroER<-matrix(nrow=nrow(dataERmacro), ncol=19)
averagedplasmaER<-matrix(nrow=nrow(dataERplasma), ncol=19)

j<-1
for(i in 15:33){
  averagedTcellER[,j]<-rowMeans(dataERTcell[,intersect(colnames(cancersdata[[i]]), ER0)], na.rm=T)
  averagedBcellER[,j]<-rowMeans(dataERBcell[,intersect(colnames(cancersdata[[i]]), ER19)], na.rm=T)
  averagedendoER[,j]<-rowMeans(dataERendo[,intersect(colnames(cancersdata[[i]]), ER18)], na.rm=T)
  averagedmacroER[,j]<-rowMeans(dataERmacro[,intersect(colnames(cancersdata[[i]]), ER2)], na.rm=T)
  if(length(intersect(colnames(cancersdata[[i]]), ER17))>1){
  averagedplasmaER[,j]<-rowMeans(dataERplasma[,intersect(colnames(cancersdata[[i]]), ER17)], na.rm=T)
  }
  j<-j+1
  }

save(averagedTcellER,  file="averagedTcellER.RData")
save(averagedBcellER,  file="averagedBcellER.RData")
save(averagedendoER,  file="averagedendoER.RData")
save(averagedmacroER,  file="averagedmacroER.RData")
save(averagedplasmaER,  file="averagedplasmaER.RData")


#doubleWGCNA #primaversione

library(WGCNA)
library(xlsx)
library(limma)
library(pheatmap)
library(Hmisc)

load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedEpi.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedStroma.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/genenames.RData")
rownames(averagedEpi)<-genenames
rownames(averagedStroma)<-genenames

###filtering data
ind<-which(colSums(is.na(averagedEpi))==0 & colSums(is.na(averagedStroma))==0)

averagedEpi<-averagedEpi[, ind]
averagedStroma<-averagedStroma[, ind]


stroma<-averagedStroma[apply(averagedStroma,1,var)>=quantile(apply(averagedStroma,1,var),0.5),]
epi<-averagedEpi[apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5),]

#merging stroma and epi
rownames(stroma)<-paste(rownames(stroma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(stroma)

data_merged_GSE5847<-rbind(stroma, epi)



################################
######Double WGCNA epi stroma tutte le cellule
##############################

library(WGCNA)
library(xlsx)
library(limma)
library(pheatmap)
library(Hmisc)
averagedEpi<-cbind(averagedEpiTNBC, averagedEpiHER2, averagedEpiER)
averagedStroma<-cbind(averagedStromaTNBC, averagedStromaHER2, averagedStromaER)

rownames(averagedEpi)<-genenames
rownames(averagedStroma)<-genenames

###filtering data
ind<-which(colSums(is.na(averagedEpi))==0 & colSums(is.na(averagedStroma))==0)

averagedEpi<-averagedEpi[, ind]
averagedStroma<-averagedStroma[, ind]


stroma<-averagedStroma[apply(averagedStroma,1,var)>=quantile(apply(averagedStroma,1,var),0.5),]
epi<-averagedEpi[apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5),]

#merging stroma and epi
rownames(stroma)<-paste(rownames(stroma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(stroma)

data_merged<-rbind(stroma, epi)
save(data_merged, file="data_merged_EpiStroma.RData")

#Functions da A_L.R
netsc<-network(data=data_merged, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(netsc, file="netsc.RData")


################################
######Double WGCNA epi - Tcell tutte le cellule
##############################
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedTcellTNBC.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedTcellHER2.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/averagedTcellER.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/genenames.RData")

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


Tcell<-averagedTcell[apply(averagedTcell,1,var)>=quantile(apply(averagedTcell,1,var),0.5),]
epi<-averagedEpi[apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5),]

#merging stroma and epi
rownames(Tcell)<-paste(rownames(Tcell), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(Tcell)

data_merged<-rbind(Tcell, epi)

#Functions da A_L.R
netsc_ET<-network(data=data_merged, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(netsc_ET, file="netsc_ET.RData")




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


Bcell<-averagedBcell[apply(averagedBcell,1,var)>=quantile(apply(averagedBcell,1,var),0.5),]
epi<-averagedEpi[apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5),]

#merging stroma and epi
rownames(Bcell)<-paste(rownames(Bcell), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(Bcell)

data_mergedB<-rbind(Bcell, epi)

#Functions da A_L.R
netsc_EB<-network(data=data_mergedB, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(netsc_EB, file="netsc_EB.RData")

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


macro<-averagedmacro[apply(averagedmacro,1,var)>=quantile(apply(averagedmacro,1,var),0.5),]
epi<-averagedEpi[apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5),]

#merging stroma and epi
rownames(macro)<-paste(rownames(macro), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(macro)

data_mergedM<-rbind(macro, epi)

#Functions da A_L.R
netsc_EM<-network(data=data_mergedM, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(netsc_EM, file="netsc_EM.RData")

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


endo<-averagedendo[apply(averagedendo,1,var)>=quantile(apply(averagedendo,1,var),0.5),]
epi<-averagedEpi[apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5),]

#merging stroma and epi
rownames(endo)<-paste(rownames(endo), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(endo)

data_mergedE<-rbind(endo, epi)

#Functions da A_L.R
netsc_EE<-network(data=data_mergedE, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(netsc_EE, file="netsc_EE.RData")

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


plasma<-averagedplasma[apply(averagedplasma,1,var)>=quantile(apply(averagedplasma,1,var),0.5),]
epi<-averagedEpi[apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5),]

#merging stroma and epi
rownames(plasma)<-paste(rownames(plasma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(plasma)

data_mergedP<-rbind(plasma, epi)

#Functions da A_L.R
netsc_EP<-network(data=data_mergedP, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
save(netsc_EP, file="netsc_EP.RData")


########
##Osservazioni

##Le CAF in tutti i sottotipi comprendono una popolazione con KRT e CD24
##E una popolazione con alte TGLN e ACTA2
#In Epi si trova una simile popolazione TAGLN e ACTA2

#se si prendono le T cells in ER e si ripete il clustering, una metà ha più alta espressione di CD24 e anche di KRT19
