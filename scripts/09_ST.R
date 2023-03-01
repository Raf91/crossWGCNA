##10x dataset
library(Seurat)
library(ggplot2)

BC <- Load10X_Spatial(data.dir = "data/Spatial transcriptomics/10x/",
                      filename = "Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5")


VlnPlot(BC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(BC, features = "nCount_Spatial") + theme(legend.position = "right")

BC <- SCTransform(BC, assay = "Spatial", verbose = FALSE)

pdf("results/featureplot_epi.pdf", 10,10)
SpatialFeaturePlot(BC, features = c("EPCAM", "EGFR", "CDH1","KRT5"), alpha = c(0.5, 1))
dev.off()
pdf("results/featureplot_stroma.pdf", 10,10)
SpatialFeaturePlot(BC, features = c("FAP", "COL1A1", "COL3A1","ACTA2"), alpha = c(0.1, 1))
dev.off()

BC <- RunPCA(BC, assay = "SCT", verbose = FALSE)
BC <- FindNeighbors(BC, reduction = "pca", dims = 1:30)
BC <- FindClusters(BC, verbose = FALSE)
BC <- RunUMAP(BC, reduction = "pca", dims = 1:30)

p1 <- DimPlot(BC, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(BC, label = TRUE, label.size = 3)
pdf("results/ST_clusters.pdf", 12,6)
p1 + p2
dev.off()

pdf("HE_BC.pdf", 6,6)
SpatialFeaturePlot(BC, features = "nCount_Spatial", alpha=0) + theme(legend.position = "right")
dev.off()

###markers to define epithelial and stromal clusters
pdf("results/vlplot_stroma1.pdf", 12,6)
VlnPlot(BC, features = c("FAP", "COL1A1", "COL3A1", "COL5A1", "ACTA2", "TAGLN"), pt.size = 0.1, group.by="seurat_clusters")
dev.off()
pdf("results/vlplot_stroma2.pdf", 12,6)
VlnPlot(BC, features = c("LUM", "FBLN1", "COL6A3", "COL1A2", "COL6A1", "COL6A2"), pt.size = 0.1, group.by="seurat_clusters")
dev.off()

pdf("results/vlplot_epi1.pdf", 12,6)
VlnPlot(BC, features = c("EPCAM", "EGFR", "CDH1", "KRT14", "ITGA6", "KRT5"), pt.size = 0.1, group.by="seurat_clusters")
dev.off()

pdf("results/vlplot_epi2.pdf", 12,6)
VlnPlot(BC, features = c("TP63", "KRT17", "MME", "KRT8", "KRT18", "KRT19"), pt.size = 0.1, group.by="seurat_clusters")
dev.off()

pdf("results/vlplot_epi3.pdf", 12,6)
VlnPlot(BC, features = c("FOXA1", "GATA3", "MUC1", "CD24", "KIT", "GABRP"), pt.size = 0.1, group.by="seurat_clusters")
dev.off()


###define spots coordinates
y<-GetTissueCoordinates(BC)[,1]
p<-hist(y, breaks=1000)#1000 isway higher than the number of peaks. If using another array, change this number
#peaks separated by 0
breaks<-p$breaks[p$counts==0]
pos_break<-which(p$counts==0)
breaks<-breaks[-which(pos_break %in% (pos_break+1))]
peaks<-cut(y, breaks = c(min(y),breaks, max(y)))
levels(peaks)<-c(1:length(unique(peaks)))
y_bin<-peaks


x<-GetTissueCoordinates(BC)[,2]
p<-hist(x, breaks=1000)
breaks<-p$breaks[p$counts==0]
pos_break<-which(p$counts==0)
breaks<-breaks[-which(pos_break %in% (pos_break+1))]
peaks<-cut(x, breaks = c(min(x),breaks, max(x)))
levels(peaks)<-c(1:length(unique(peaks)))
x_bin<-peaks

x_bin<-as.numeric(x_bin)
y_bin<-as.numeric(y_bin)
#transform ycoordinates so that they represent the height of an equilateral triangle
y_bin<-(y_bin*sqrt(3))
#compute the distances between spots
coords<-cbind(x_bin, y_bin)
eucl_dist<-dist(coords)

##############define compartments and get expression data
epi_spots<-which(BC$seurat_clusters %in% c(3, 7, 9:15))
stroma_spots<-which(BC$seurat_clusters %in% c(0, 1:2,4:6,8,16,17))
expr_data<-BC@assays$SCT@data
#####################################

class<-rep(NA, ncol(expr_data))
class[epi_spots]<-"Epi"
class[stroma_spots]<-"Stroma"

coords<-data.frame(x=x_bin, y=y_bin, class=class)

###smooths gene expression using the weighted mean of neighbouring spots in the same compartment
averaged_expr_all<-matrix(ncol=ncol(expr_data), nrow=nrow(expr_data))

for(es in 1:ncol(expr_data)){
  class_es<-class[es]
  if(class_es %in% c("Epi", "Stroma")){
    dist_es<-as.matrix(eucl_dist)[,es]
    sel_spots<-which(dist_es<5 & class==class_es)
    weights<-1/(dist_es[sel_spots]+1)
    if(length(sel_spots)>1){
      averaged_expr_all[,es]<-apply((exp(expr_data[,sel_spots])-1), 1, function(x){log(weighted.mean(x, weights)+1)})
    } else {
      averaged_expr_all[,es]<-expr_data[,sel_spots]
    }
  }
}
rownames(averaged_expr_all)<-rownames(expr_data)
save(averaged_expr_all, file="results/averaged_expr_all3.RData")

###coordinates of epi spots
pdf("results/spots_class.pdf",5,6)
plot(as.numeric(x_bin[epi_spots]),-as.numeric(y_bin[epi_spots]), pch=19, cex=0.5, xlab="x", ylab="y")
points(as.numeric(x_bin[stroma_spots]),-as.numeric(y_bin[stroma_spots]), pch=19, col="red", cex=0.5)
dev.off()

##remove the most intermixed clusters as likely too difficult to distinguish epi and stroma
epi_spots<-which(BC$seurat_clusters %in% c(7, 9:12))
stroma_spots<-which(BC$seurat_clusters %in% c(0, 2,4,6,8,16,17))

pdf("results/spots_class_filt.pdf",5,6)
plot(as.numeric(x_bin[epi_spots]),-as.numeric(y_bin[epi_spots]), pch=19, cex=0.5, xlab="x", ylab="y")
points(as.numeric(x_bin[stroma_spots]),-as.numeric(y_bin[stroma_spots]), pch=19, col="red", cex=0.5)
dev.off()

###selects epi and stroma spots not isolated
##epi with at least 1 neighbouring epi spot
##stroma with at least 1 neighbouring stroma spot
##epi with at least 1 selected sroma spot
included_es_spots<-c()
for(es in epi_spots){
  which_x_coord<-which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
  which_y_coord<-which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
  sel_spots<-intersect(which_x_coord, which_y_coord)
  neighbour_spots_epi<-intersect(epi_spots, sel_spots)
  if(length(neighbour_spots_epi)>1){
    included_es_spots<-c(included_es_spots, es)
  }
}


included_ss_spots<-c()
for(ss in stroma_spots){
  which_x_coord<-which(x_bin>=x_bin[ss]-2 & x_bin<=x_bin[ss]+2)
  which_y_coord<-which(y_bin>=y_bin[ss]-sqrt(3)-0.1 & y_bin<=y_bin[ss]+sqrt(3)+0.1)
  sel_spots<-intersect(which_x_coord, which_y_coord)
  neighbour_spots_stroma<-intersect(stroma_spots, sel_spots)
  if(length(neighbour_spots_stroma)>1){
    included_ss_spots<-c(included_ss_spots, ss)
  }
}

pdf("results/spots_class_filt2.pdf",5,6)
plot(as.numeric(x_bin[included_es_spots]),-as.numeric(y_bin[included_es_spots]), cex=0.5, pch=19, xlab="x", ylab="y")
points(as.numeric(x_bin[included_ss_spots]),-as.numeric(y_bin[included_ss_spots]), cex=0.5, pch=19, col="red")
dev.off()

sel_es_spots<-c()
for(es in included_es_spots){
  which_x_coord<-which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
  which_y_coord<-which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
  sel_spots<-intersect(which_x_coord, which_y_coord)
  neighbour_spots_stroma<-intersect(included_ss_spots, sel_spots)
  if(length(neighbour_spots_stroma)>0){
    sel_es_spots<-c(sel_es_spots, es)
  }
}


pdf("results/spots_class_filt3.pdf",5,6)
plot(as.numeric(x_bin[sel_es_spots]),-as.numeric(y_bin[sel_es_spots]), cex=0.5, pch=19, xlab="x", ylab="y")
points(as.numeric(x_bin[included_ss_spots]),-as.numeric(y_bin[included_ss_spots]), cex=0.5, pch=19, col="red")
dev.off()

#############
###takes the stromal spots next to each epi spot and
#creates a matched epi and stroma expression matrices
#averages the expression of stroma spots in contact with a specific epithelial spot

included_spots<-c()
epi_expr_all<-matrix(ncol=1, nrow=nrow(averaged_expr_all))
stroma_expr_all<-matrix(ncol=1, nrow=nrow(averaged_expr_all))

for(es in sel_es_spots){
  which_x_coord<-which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
  which_y_coord<-which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
  sel_spots<-intersect(which_x_coord, which_y_coord)
  sel_spots_stroma<-intersect(stroma_spots, sel_spots)
  neighbour_spots_epi<-intersect(sel_es_spots, sel_spots)

  if(length(sel_spots_stroma)>0 & length(neighbour_spots_epi)>1){
    if(length(sel_spots_stroma)==1){
      stroma_expr<-averaged_expr_all[,sel_spots_stroma]
    } else if(length(sel_spots_stroma)>1 ){
      stroma_expr<-log(rowMeans(exp(averaged_expr_all[,sel_spots_stroma])-1)+1)
    }

    #epi_expr<-log(rowMeans(exp(averaged_expr_all[,neighbour_spots_epi])-1)+1)
    epi_expr<-averaged_expr_all[,es]
    epi_expr_all<-cbind(epi_expr_all, epi_expr)
    stroma_expr_all<-cbind(stroma_expr_all, stroma_expr)

    included_spots<-c(included_spots, es)
  }

}
rownames(epi_expr_all)<-rownames(averaged_expr_all)
rownames(stroma_expr_all)<-rownames(averaged_expr_all)
epi_expr_all<-epi_expr_all[,-1]
stroma_expr_all<-stroma_expr_all[,-1]

###################################
included_spots_epi<-c()
included_spots_stroma<-c()
for(es in included_spots){
  which_x_coord<-which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
  which_y_coord<-which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
  sel_spots<-intersect(which_x_coord, which_y_coord)
  sel_spots_stroma<-intersect(stroma_spots, sel_spots)
  included_spots_stroma<-c(included_spots_stroma, sel_spots_stroma)
  included_spots_epi<-c(included_spots_epi, rep(es, length(sel_spots_stroma)))
}

##points in contact
pdf("results/spots_class_filt4.pdf",5,6)
plot(as.numeric(x_bin[included_spots_epi]),-as.numeric(y_bin[included_spots_epi]), cex=0.5, pch=19, xlab="x", ylab="y")
points(as.numeric(x_bin[included_spots_stroma]),-as.numeric(y_bin[included_spots_stroma]), cex=0.5, pch=19, col="red")
dev.off()

##midpoints coordinates
df<-data.frame(x_coord= c(x_bin[included_spots_epi],x_bin[included_spots_stroma]),y_coord= c(-(y_bin[included_spots_epi]), -(y_bin[included_spots_stroma])))

midpoints_x<-(df$x_coord[c(1:length(included_spots_epi))]+df$x_coord[-c(1:length(included_spots_epi))])/2
midpoints_y<-(df$y_coord[c(1:length(included_spots_epi))]+df$y_coord[-c(1:length(included_spots_epi))])/2

df<-data.frame(x_coord= c(x_bin[included_spots_epi],x_bin[included_spots_stroma], midpoints_x),y_coord= c(-(y_bin[included_spots_epi]), -(y_bin[included_spots_stroma]), midpoints_y),
               compartment=c(rep("epi", length(included_spots_epi)),rep("stroma", length(included_spots_stroma)), rep("midpoint", length(midpoints_x))))

ggplot(data=df, aes(x=x_coord, y=y_coord, colour=compartment))+geom_point(size=1)
###################
###trim isolated spots


###filters the most variable genes in both.
#Many with var==0. They have to be removed from both matrices
#otherwise correlations in crossWGCNA are NA

var_stroma<-apply(stroma_expr_all,1,var)
var_epi<-apply(epi_expr_all,1,var)

genes<-intersect(which(var_stroma>quantile(var_stroma, 0.75, na.rm=T)),
                 which(var_epi>quantile(var_epi, 0.75, na.rm=T)))
genes<-rownames(averaged_expr_all)[genes]

stroma<-stroma_expr_all[genes,]
epi<-epi_expr_all[genes,]

#merge stroma and epi
rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")
colnames(epi)<-colnames(stroma)

data_merged<-rbind(stroma, epi)

library(WGCNA)
source("scripts/crossWGCNA_functions_netdiff.R")

adj<-Adjacency(data=data_merged, Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
save(adj, file="results/ST_adj_weighted_netdiff3_s.RData")

net<-network(data=data_merged, Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
save(net, file="results/ST_net_weighted_netdiff3_s.RData")

mods<-crossWGCNA(data=data_merged, Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
#functional enrichment of the most communicative genes
#######################
library(msigdbr)
library(fgsea)

m_df = msigdbr(species = "Homo sapiens", category = "C5")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

coef_gsea1 <- net$kExt1 / net$kInt1
names(coef_gsea1) <- gsub("_tis1", "", names(coef_gsea1))

fgseaRes1 <- fgseaMultilevel(m_list, coef_gsea1)
write.xlsx(data.frame(fgseaRes1), file="results/ST_GSEA_stroma_ratio2_s.xlsx")


coef_gsea2 <- net$kExt2 / net$kInt2
names(coef_gsea2) <- gsub("_tis2", "", names(coef_gsea2))

fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea2)
write.xlsx(data.frame(fgseaRes2), file="results/ST_GSEA_epi_ratio2_s.xlsx")

sort(coef_gsea1, decreasing = T)[1:10]
sort(coef_gsea2, decreasing = T)[1:10]

plot(coef_gsea1, coef_gsea2)

sort(adj["IGLC2_tis1", grep("tis2", colnames(adj))], decreasing=T)[1:10]
sort(A["GSTP1_tis1", grep("tis2", colnames(A))], decreasing=T)[1:10]

sort(adj["MT-CO2_tis2", grep("tis1", colnames(adj))], decreasing=T)[1:10]
sort(A["VIM_tis2", grep("tis1", colnames(A))], decreasing=T)[1:10]

#MALAT1 UGCG
# SFRP2 MORF4L2

m<-max(adj[grep("tis1", rownames(adj)), grep("tis2", colnames(adj))])
m<-sort(adj[grep("tis1", rownames(adj)), grep("tis2", colnames(adj))], decreasing=T)[3]
which(adj==m, arr.ind=T)
#COL6A2, EFNA1
gene1<-"LGALS9"
gene2<-"ACTG1"

df<-data.frame(x_coord= c(x_bin[epi_spots],x_bin[stroma_spots], x_bin[included_spots]),y_coord= c(-(y_bin[epi_spots]), -(y_bin[stroma_spots]), -y_bin[included_spots]),
               compartment=c(rep("epi", length(epi_spots)),rep("stroma", length(stroma_spots)), rep("edge", length(included_spots))),
               gene1=c(averaged_expr_all[gene1,epi_spots], averaged_expr_all[gene1,stroma_spots], averaged_expr_all[gene1,included_spots]),
               gene2=c(averaged_expr_all[gene2,epi_spots], averaged_expr_all[gene2,stroma_spots], averaged_expr_all[gene2,included_spots]))


df$gene1.scale <- with(df, (gene1-min(gene1, na.rm=T))/diff(range(gene1, na.rm=T)))
df$gene2.scale <- with(df, (gene2-min(gene2, na.rm=T))/diff(range(gene2, na.rm=T)))
df<-df[!is.na(df$gene1.scale) & !is.na(df$gene2.scale),]

#ggplot(data=df, aes(x=x_coord, y=y_coord))+geom_tile(size=2, colour=rgb(red=df$gene1.scale,green=df$gene2.scale,blue=0))+
#  geom_point(data=df, aes(x=x_coord, y=y_coord, size=compartment), color="black") +
#  scale_size_manual(values = c("epi" = 0, "stroma"=0, "edge"=0.8))

pdf(paste("results/", gene1, gene2, "scatter.pdf", sep="_"),6,6)
plot(stroma[paste(gene1, "tis1", sep="_"),], epi[paste(gene2, "tis2", sep="_"),], pch=19, xlab=paste(gene1, "stroma"), ylab=paste(gene2, "epi"))
dev.off()


##########plot communication

##color midboint by the amount of communication
library(ggnewscale)
df<-data.frame(x_coord= c(x_bin[included_spots_epi],x_bin[included_spots_stroma], midpoints_x),
               y_coord= c(-(y_bin[included_spots_epi]), -(y_bin[included_spots_stroma]), midpoints_y),
               compartment=c(rep("epi", length(included_spots_epi)),rep("stroma", length(included_spots_stroma)), rep("midpoint", length(midpoints_x))),
               gene1=c(averaged_expr_all[gene1,included_spots_epi], averaged_expr_all[gene1,included_spots_stroma], rep(NA,  length(midpoints_x))),
               gene2=c(averaged_expr_all[gene2,included_spots_epi], averaged_expr_all[gene2,included_spots_stroma], rep(NA,  length(midpoints_x))),
               comm_score=c(rep(NA, length(c(included_spots_epi, included_spots_stroma))), averaged_expr_all[gene1,included_spots_stroma]*averaged_expr_all[gene2,included_spots_epi]))

ggplot(data=df[df$compartment!="midpoint",], aes(x=x_coord, y=y_coord, colour=gene1 ))+geom_point()+scale_color_continuous(type = "viridis")+
  new_scale_colour()+
  geom_point(data=df[df$compartment=="midpoint",], size=1, aes(colour=comm_score))+scale_colour_gradientn(colours = c("purple", "orange"))



df<-data.frame(x_coord= c(x_bin[epi_spots],x_bin[stroma_spots], x_bin[included_spots]),y_coord= c(-(y_bin[epi_spots]), -(y_bin[stroma_spots]), -y_bin[included_spots]),
               compartment=c(rep("epi", length(epi_spots)),rep("stroma", length(stroma_spots)), rep("edge", length(included_spots))),
               gene1=c(averaged_expr_all[gene1,epi_spots], averaged_expr_all[gene1,stroma_spots], averaged_expr_all[gene1,included_spots]),
               gene2=c(averaged_expr_all[gene2,epi_spots], averaged_expr_all[gene2,stroma_spots], averaged_expr_all[gene2,included_spots]))


df$gene1.scale <- with(df, (gene1-min(gene1, na.rm=T))/diff(range(gene1, na.rm=T)))
df$gene2.scale <- with(df, (gene2-min(gene2, na.rm=T))/diff(range(gene2, na.rm=T)))
df<-df[!is.na(df$gene1.scale) & !is.na(df$gene2.scale),]

df_midpoint<-data.frame(x_coord= midpoints_x,
                        y_coord= midpoints_y,
                        comm_score=averaged_expr_all[gene1,included_spots_stroma]*averaged_expr_all[gene2,included_spots_epi])

ggplot(data=df, aes(x=x_coord, y=y_coord))+geom_point(size=1, colour=rgb(red=df$gene1.scale,green=df$gene2.scale,blue=0))+
  new_scale_colour()+geom_point(data=df_midpoint, size=1, aes(colour=comm_score))+scale_color_continuous(type = "viridis")


pdf(paste("results/", gene1, gene2, "comm.pdf", sep="_"),7,6)
ggplot(data=df, aes(x=x_coord, y=y_coord, colour=compartment))+
  geom_point(data=df_midpoint, aes(x=x_coord, y=y_coord, colour=comm_score))+scale_colour_gradientn(colours = c("purple", "orange"))+
  theme_classic()
dev.off()


########
gene<-gene1
df<-data.frame(x_coord= c(x_bin[epi_spots],x_bin[stroma_spots], x_bin[included_spots]),y_coord= c(-(y_bin[epi_spots]), -(y_bin[stroma_spots]), -y_bin[included_spots]),
               compartment=c(rep("epi", length(epi_spots)),rep("stroma", length(stroma_spots)), rep("edge", length(included_spots))),
               gene=c(averaged_expr_all[gene,epi_spots], averaged_expr_all[gene,stroma_spots], averaged_expr_all[gene,included_spots]))
df_midpoint<-data.frame(x_coord= midpoints_x,
                        y_coord= midpoints_y,
                        gene=averaged_expr_all[gene,included_spots_stroma])

pdf(paste("results/", gene1, "expr.pdf", sep="_"),7,6)
ggplot(data=df, aes(x=x_coord, y=y_coord, colour=gene))+geom_point()+
  scale_color_continuous(type = "viridis")+
  new_scale_colour()+geom_point(data=df_midpoint, size=1, colour="black")+theme_classic()
dev.off()

gene<-gene2
df<-data.frame(x_coord= c(x_bin[epi_spots],x_bin[stroma_spots], x_bin[included_spots]),y_coord= c(-(y_bin[epi_spots]), -(y_bin[stroma_spots]), -y_bin[included_spots]),
               compartment=c(rep("epi", length(epi_spots)),rep("stroma", length(stroma_spots)), rep("edge", length(included_spots))),
               gene=c(averaged_expr_all[gene,epi_spots], averaged_expr_all[gene,stroma_spots], averaged_expr_all[gene,included_spots]))
df_midpoint<-data.frame(x_coord= midpoints_x,
                        y_coord= midpoints_y,
                        gene=averaged_expr_all[gene,included_spots_stroma])

pdf(paste("results/", gene2, "expr.pdf", sep="_"),7,6)
ggplot(data=df, aes(x=x_coord, y=y_coord, colour=gene))+geom_point()+
  scale_color_continuous(type = "viridis")+
  new_scale_colour()+geom_point(data=df_midpoint, size=1, colour="black")+theme_classic()
dev.off()


###########
##compare the expression before and after the smoothing
###########

gene<-gene1
df<-data.frame(x_coord= c(x_bin[c(stroma_spots, epi_spots)], x_bin[c(stroma_spots, epi_spots)]),
               y_coord= c(-(y_bin[c(stroma_spots, epi_spots)]), -(y_bin[c(stroma_spots, epi_spots)])),
               gene=c(expr_data[gene,c(stroma_spots, epi_spots)], averaged_expr_all[gene,c(stroma_spots, epi_spots)]),
               type=c(rep("original", length(c(stroma_spots, epi_spots))), rep("averaged", length(c(stroma_spots, epi_spots)))))
df$type<-factor(df$type, levels=c("original", "averaged"))
df_midpoint<-data.frame(x_coord= midpoints_x,
                        y_coord= midpoints_y,
                        gene=averaged_expr_all[gene,included_spots_stroma])

pdf(paste("results/", gene1, "origVSave.pdf", sep="_"),12,6)
ggplot(data=df, aes(x=x_coord, y=y_coord, colour=gene))+geom_point()+
  scale_color_continuous(type = "viridis")+facet_grid(~type)+
  new_scale_colour()+geom_point(data=df_midpoint, size=1, colour="black")+theme_classic()
dev.off()


gene<-gene2
df<-data.frame(x_coord= c(x_bin[c(stroma_spots, epi_spots)], x_bin[c(stroma_spots, epi_spots)]),
               y_coord= c(-(y_bin[c(stroma_spots, epi_spots)]), -(y_bin[c(stroma_spots, epi_spots)])),
               gene=c(expr_data[gene,c(stroma_spots, epi_spots)], averaged_expr_all[gene,c(stroma_spots, epi_spots)]),
               type=c(rep("original", length(c(stroma_spots, epi_spots))), rep("averaged", length(c(stroma_spots, epi_spots)))))
df$type<-factor(df$type, levels=c("original", "averaged"))
df_midpoint<-data.frame(x_coord= midpoints_x,
                        y_coord= midpoints_y,
                        gene=averaged_expr_all[gene,included_spots_stroma])

pdf(paste("results/", gene2, "origVSave.pdf", sep="_"),12,6)
ggplot(data=df, aes(x=x_coord, y=y_coord, colour=gene))+geom_point()+
  scale_color_continuous(type = "viridis")+facet_grid(~type)+
  new_scale_colour()+geom_point(data=df_midpoint, size=1, colour="black")+theme_classic()
dev.off()

##########
##scheme of averaging
##########

  class_es<-class[es]
     dist_es<-as.matrix(eucl_dist)[,es]
    sel_spots<-which(dist_es<5 & class==class_es)
    weights<-1/(dist_es[sel_spots]+1)
 all_weights<-rep(0,ncol(expr_data))
 all_weights[sel_spots]<-weights
#all_weights<-dist_es
 df<-data.frame(x_coord= c(x_bin[c(stroma_spots, epi_spots)]),
                   y_coord= c(-(y_bin[c(stroma_spots, epi_spots)])),
                   weights=all_weights[c(stroma_spots, epi_spots)])
    df_midpoint<-data.frame(x_coord= midpoints_x,
                            y_coord= midpoints_y,
                            gene=averaged_expr_all[gene,included_spots_stroma])

    ggplot(data=df, aes(x=x_coord, y=y_coord, colour=weights))+geom_point()+scale_color_gradient(low = "white", high = "red")



    cor_inspect<-function(data=data_merged_GSE88715, gene1, gene2, comp1="_tis1", comp2="_tis2"){

      df<-data.frame(gene1=c(data[paste(gene1, comp1, sep=""),], data[paste(gene1, comp2, sep=""),], data[paste(gene1, comp1, sep=""),], data[paste(gene1, comp2, sep=""),]),
                     gene2=c(data[paste(gene2, comp1, sep=""),], data[paste(gene2, comp2, sep=""),], data[paste(gene2, comp2, sep=""),], data[paste(gene2, comp2, sep=""),]),
                     compartment=c(rep(c("comp1 vs comp1","comp2 vs comp2","comp1 vs comp2", "comp2 vs comp1"), each=ncol(data)) ))
      p<-ggplot(df, aes(x=gene1, y=gene2))+geom_point()+facet_wrap(.~compartment)+geom_smooth(method = "lm")+theme_classic()+labs(x=gene1, y=gene2)
      return(p)
    }


    pdf("results/ST_COL6A2_EFNA1_cor_inspect.pdf", 6, 6)
    cor_inspect(data_merged, gene1="COL6A2", gene2="EFNA1")
    dev.off()

    pdf("results/ST_HLA-DRB1_ACTG1_cor_inspect.pdf", 6, 6)
    cor_inspect(data_merged, gene1="HLA-DRB1", gene2="ACTG1")
    dev.off()

    pdf("results/ST_LGALS9_ACTG1_cor_inspect.pdf", 6, 6)
    cor_inspect(data_merged, gene1="LGALS9", gene2="ACTG1")
    dev.off()


    df<-cytoscape_net(adjacency = adj, dataset=as.matrix(data_merged), gene="COL6A2", comp1="tis1", comp2="tis2", num=10)
    df$Target<-factor(df$Target, levels=df$Target[df$Edge_type=="inter"][order(df$Weight[df$Edge_type=="inter"], decreasing = T)])
    pdf("results/ST_COL6A2_top10_corrs.pdf", 8,7 )
    ggplot(df, aes(x=Target, y=Weight))+geom_col()+facet_grid(Edge_type~.)+theme_classic()
    dev.off()


    ################
    #### GO for modules
    #################
    library(clusterProfiler)
    modules<-mods[[2]]$colors

    #GO for stroma
    background_1<-gsub("_tis1", "", names(modules)[grep("_tis1",names(modules))])
    background_2<-gsub("_tis2", "", names(modules)[grep("_tis2",names(modules))])
    mod_sel<-unique(modules)
    ego_1<-list()
    for(i in 1:length(mod_sel)){
      module_1<-gsub("_tis1", "", names(modules)[modules==mod_sel[i]][grep("_tis1",names(modules)[modules==mod_sel[i]])])
      if(length(module_1)>0 & mod_sel[i]!=0){
        ego_1[[i]]<- enrichGO(gene = module_1,
                              keyType="SYMBOL",
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 1,
                              universe=background_1, OrgDb="org.Hs.eg.db")

        if(nrow(summary(ego_1[[i]]))>0){
          write.xlsx( summary(ego_1[[i]]) , paste("Stroma", "ST", mod_sel[i], "GO_s.xlsx", sep="_"))
        }
      }
    }


    ##GO for epi
    ego_2<-list()
    for(i in 1:length(mod_sel)){
      module_2<-gsub("_tis2", "", names(modules)[modules==mod_sel[i]][grep("_tis2",names(modules)[modules==mod_sel[i]])])
      if(length(module_2)>0){
        ego_2[[i]]<- enrichGO(gene = module_2,
                              keyType="SYMBOL",
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 1,
                              universe=background_2, OrgDb="org.Hs.eg.db")

        if(nrow(summary(ego_2[[i]]))>0){
          write.xlsx( summary(ego_2[[i]]) , paste("Epi", "ST", mod_sel[i], "GO_s.xlsx", sep="_"))
        }

      }
    }


    kwithin<-degrees_mod(data=data_merged, modules=mods[[2]]$colors, Adj_type = "signed",
                         cortype = "spearman",
                         pval = "none",
                         thr = 0.05,
                         beta = 6,
                         comp1 = "_tis1",
                         comp2 = "_tis2")

    #weighted module expression

    genes<-unique(gsub(comp2, "", gsub(comp1, "", mod)))
    c(kwithin[[1]]$kExt1, kwithin[[1]]$kExt2)[mod]
    c(kwithin[[1]]$kExt1, kwithin[[1]]$kExt2)[mod]

    kwithin[[1]]$kExt2[mod[grep("_tis2", mod)]]

    ########
    i<-3
    mod<-names(modules)[which(modules==unique(modules)[i])]

     weights<-kwithin[[i]]$kExt1[mod[grep("_tis1", mod)]]
    epi_wm1<-apply((averaged_expr_all[gsub("_tis1", "", mod[grep("_tis1", mod)]),epi_spots]), 2, function(x){weighted.mean(x, weights)})
    stroma_wm1<-apply((averaged_expr_all[gsub("_tis1", "", mod[grep("_tis1", mod)]),stroma_spots]), 2, function(x){weighted.mean(x, weights)})
    edge_wm1<-apply((averaged_expr_all[gsub("_tis1", "", mod[grep("_tis1", mod)]),included_spots]), 2, function(x){weighted.mean(x, weights)})

    #epi_wm<-colMeans(averaged_expr_all[gsub("_tis1", "", mod[grep("_tis1", mod)]),epi_spots])
    #stroma_wm<-colMeans(averaged_expr_all[gsub("_tis1", "", mod[grep("_tis1", mod)]),stroma_spots])
    #edge_wm<-colMeans(averaged_expr_all[gsub("_tis1", "", mod[grep("_tis1", mod)]),included_spots])

    df<-data.frame(x_coord= c(x_bin[epi_spots],x_bin[stroma_spots], x_bin[included_spots]),y_coord= c(-(y_bin[epi_spots]), -(y_bin[stroma_spots]), -y_bin[included_spots]),
                   compartment=c(rep("epi", length(epi_spots)),rep("stroma", length(stroma_spots)), rep("edge", length(included_spots))),
                   gene=c(epi_wm1, stroma_wm1, edge_wm1))
    df_midpoint<-data.frame(x_coord= midpoints_x,
                            y_coord= midpoints_y,
                            gene=averaged_expr_all[gene,included_spots_stroma])

    p1<-ggplot(data=df, aes(x=x_coord, y=y_coord, colour=gene))+geom_point()+
      scale_color_continuous(type = "viridis")+
      geom_point(data=df_midpoint, size=1, colour="black")+theme_classic()

    weights<-kwithin[[i]]$kExt2[mod[grep("_tis2", mod)]]
    epi_wm2<-apply((averaged_expr_all[gsub("_tis2", "", mod[grep("_tis2", mod)]),epi_spots]), 2, function(x){weighted.mean(x, weights)})
    stroma_wm2<-apply((averaged_expr_all[gsub("_tis2", "", mod[grep("_tis2", mod)]),stroma_spots]), 2, function(x){weighted.mean(x, weights)})
    edge_wm2<-apply((averaged_expr_all[gsub("_tis2", "", mod[grep("_tis2", mod)]),included_spots]), 2, function(x){weighted.mean(x, weights)})

    #epi_wm<-colMeans(averaged_expr_all[gsub("_tis2", "", mod[grep("_tis2", mod)]),epi_spots])
    #stroma_wm<-colMeans(averaged_expr_all[gsub("_tis2", "", mod[grep("_tis2", mod)]),stroma_spots])
    #edge_wm<-colMeans(averaged_expr_all[gsub("_tis2", "", mod[grep("_tis2", mod)]),included_spots])

    df<-data.frame(x_coord= c(x_bin[epi_spots],x_bin[stroma_spots], x_bin[included_spots]),y_coord= c(-(y_bin[epi_spots]), -(y_bin[stroma_spots]), -y_bin[included_spots]),
                   compartment=c(rep("epi", length(epi_spots)),rep("stroma", length(stroma_spots)), rep("edge", length(included_spots))),
                   gene=c(epi_wm2, stroma_wm2, edge_wm2))
    df_midpoint<-data.frame(x_coord= midpoints_x,
                            y_coord= midpoints_y,
                            gene=averaged_expr_all[gene,included_spots_stroma])

    p2<-ggplot(data=df, aes(x=x_coord, y=y_coord, colour=gene))+geom_point()+
      scale_color_continuous(type = "viridis")+
      geom_point(data=df_midpoint, size=1, colour="black")+theme_classic()


    ##midpoints coordinates
    df<-data.frame(x_coord= c(x_bin[included_spots],x_bin[included_spots]),y_coord= c(-(y_bin[included_spots]), -(y_bin[included_spots])))

    midpoints_x1<-(df$x_coord[c(1:length(included_spots))]+df$x_coord[-c(1:length(included_spots))])/2
    midpoints_y1<-(df$y_coord[c(1:length(included_spots))]+df$y_coord[-c(1:length(included_spots))])/2

    df_midpoint<-data.frame(x_coord= midpoints_x1,
                            y_coord= midpoints_y1,
                            comm_score=edge_wm1*edge_wm2)

    ggplot(data=df, aes(x=x_coord, y=y_coord))+geom_point(data=df_midpoint, size=1, aes(colour=comm_score))+scale_color_continuous(type = "viridis")


    p1
    p2


    mod<-names(modules)[which(modules==unique(modules)[i])]
    weights<-kwithin[[i]]$kExt1[mod[grep("_tis1", mod)]]
    wm1<-apply((data_merged[mod[grep("_tis1", mod)],]), 2, function(x){weighted.mean(x, weights)})

    weights<-kwithin[[i]]$kExt2[mod[grep("_tis2", mod)]]
    wm2<-apply((data_merged[mod[grep("_tis2", mod)],]), 2, function(x){weighted.mean(x, weights)})

    plot(wm1, wm2)

pdf("results/ST_mod1_tis1.pdf", 5,6)
p1
dev.off()

pdf("results/ST_mod1_tis2.pdf", 5,6)
p2
dev.off()

pdf("results/ST_mod1_corr.pdf", 5,5)
plot(wm1, wm2,pch=19, xlab="weighted expression comp1", ylab="weighted expression comp2")
dev.off()

pdf("results/ST_mod1_comm.pdf", 5,6)
ggplot(data=df, aes(x=x_coord, y=y_coord))+geom_point(data=df_midpoint, size=1, aes(colour=comm_score))+scale_color_continuous(type = "viridis")
dev.off()



