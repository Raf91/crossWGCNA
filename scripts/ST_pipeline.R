##10x dataset
source("scripts/crossWGCNA_functions_netdiff.R")
source("scripts/ST_utilities.R")
library(WGCNA)
library(Seurat)
library(ggplot2)
library(openxlsx)

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

coords<-spots_coords(BC)

epi_spots<-which(BC$seurat_clusters %in% c(3, 7, 9:15))
stroma_spots<-which(BC$seurat_clusters %in% c(0, 1:2,4:6,8,16,17))
expr_data<-BC@assays$SCT@data
#####################################

class<-rep(NA, ncol(expr_data))
class[epi_spots]<-"Epi"
class[stroma_spots]<-"Stroma"

averaged_expr_all_fun<-expr_smooth(expr_data=expr_data, coords=coords, max_dist=5, spots_class=class, sel_class=c("Epi", "Stroma"))

epi_spots<-which(BC$seurat_clusters %in% c(7, 9:12))
stroma_spots<-which(BC$seurat_clusters %in% c(0, 2,4,6,8,16,17))


sf<-spots_filt(coords, tis1_spots=epi_spots, tis2_spots=stroma_spots)

md<-merged_dataset(sel_spots=sf, coords=coords, averaged_expr_all=averaged_expr_all_fun, var_thr=0.75, comp1="_tis1", comp2="_tis2")


adj<-Adjacency(data=md[[1]], Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
save(adj, file="results/ST_adj_pipeline.RData")

net<-network(data=md[[1]], Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
save(net, file="results/ST_net_pipeline.RData")

mods<-crossWGCNA(data=md[[1]], Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
save(mods, file="results/ST_mods_pipeline.RData")


#functional enrichment of the most communicative genes
#######################
library(msigdbr)
library(fgsea)

m_df = msigdbr(species = "Homo sapiens", category = "C5")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

coef_gsea1 <- net$kExt1 / net$kInt1
names(coef_gsea1) <- gsub("_tis1", "", names(coef_gsea1))

fgseaRes1 <- fgseaMultilevel(m_list, coef_gsea1)
write.xlsx(data.frame(fgseaRes1), file="results/ST_GSEA_epi_ratio2_s.xlsx")


coef_gsea2 <- net$kExt2 / net$kInt2
names(coef_gsea2) <- gsub("_tis2", "", names(coef_gsea2))

fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea2)
write.xlsx(data.frame(fgseaRes2), file="results/ST_GSEA_stroma_ratio2_s.xlsx")

m<-max(adj[grep("tis2", rownames(adj)), grep("tis1", colnames(adj))])
m<-sort(adj[grep("tis1", rownames(adj)), grep("tis2", colnames(adj))], decreasing=T)[2]
which(adj==m, arr.ind=T)

#COL6A2, EFNA1
gene1<-"EFNA1"
gene2<-"COL6A2"
x_bin<-coords[,1]
y_bin<-coords[,2]

df<-data.frame(x_coord= c(x_bin[epi_spots],x_bin[stroma_spots], x_bin[md[[2]]]),y_coord= c(-(y_bin[epi_spots]), -(y_bin[stroma_spots]), -y_bin[md[[2]]]),
               compartment=c(rep("epi", length(epi_spots)),rep("stroma", length(stroma_spots)), rep("edge", length(md[[2]]))),
               gene1=c(averaged_expr_all_fun[gene1,epi_spots], averaged_expr_all_fun[gene1,stroma_spots], averaged_expr_all_fun[gene1,md[[2]]]),
               gene2=c(averaged_expr_all_fun[gene2,epi_spots], averaged_expr_all_fun[gene2,stroma_spots], averaged_expr_all_fun[gene2,md[[2]]]))


df$gene1.scale <- with(df, (gene1-min(gene1, na.rm=T))/diff(range(gene1, na.rm=T)))
df$gene2.scale <- with(df, (gene2-min(gene2, na.rm=T))/diff(range(gene2, na.rm=T)))
df<-df[!is.na(df$gene1.scale) & !is.na(df$gene2.scale),]


pdf(paste("results/", gene1, gene2, "scatter.pdf", sep="_"),5,5)
plot(md[[1]][paste(gene2, "tis2", sep="_"),], md[[1]][paste(gene1, "tis1", sep="_"),], pch=19, ylab=paste(gene1, "epi"), xlab=paste(gene2, "stroma"))
dev.off()


##########plot communication
sp<-boundary_spots(coords, included_spots=md[[2]], tis2_spots = stroma_spots)
##color midboint by the amount of communication
library(ggnewscale)
midpoints<-midpoints_def(coords, sp)

p<-plot_comm(gene1, gene2, averaged_expr_all=averaged_expr_all_fun, coords, included_spots=md[[2]], sel_spots=sp, tis1_spots=epi_spots, tis2_spots=stroma_spots, midpoints)

pdf(paste("results/", gene1, gene2, "comm.pdf", sep="_"),6,5)
print(p)
dev.off()


########
p1<-plot_expr(gene1, averaged_expr_all=averaged_expr_all_fun, coords, included_spots=md[[2]],
          tis1_spots=stroma_spots, tis2_spots=epi_spots, midpoints)
p2<-plot_expr(gene2, averaged_expr_all=averaged_expr_all_fun, coords, included_spots=md[[2]],
              tis1_spots=stroma_spots, tis2_spots=epi_spots, midpoints )

pdf(paste("results/", gene1, "expr.pdf", sep="_"),6,5)
print(p1)
dev.off()

pdf(paste("results/", gene2, "expr.pdf", sep="_"),6,5)
print(p2)
dev.off()


###########
##compare the expression before and after the smoothing
###########

gene<-gene1
df<-data.frame(x_coord= c(x_bin[c(stroma_spots, epi_spots)], x_bin[c(stroma_spots, epi_spots)]),
               y_coord= c(-(y_bin[c(stroma_spots, epi_spots)]), -(y_bin[c(stroma_spots, epi_spots)])),
               gene=c(expr_data[gene,c(stroma_spots, epi_spots)], averaged_expr_all_fun[gene,c(stroma_spots, epi_spots)]),
               type=c(rep("original", length(c(stroma_spots, epi_spots))), rep("averaged", length(c(stroma_spots, epi_spots)))))
df$type<-factor(df$type, levels=c("original", "averaged"))
df_midpoint<-data.frame(x_coord= midpoints[[1]],
                        y_coord= midpoints[[2]])

pdf(paste("results/", gene1, "origVSave.pdf", sep="_"),9,5)
ggplot(data=df, aes(x=x_coord, y=y_coord, colour=gene))+geom_point()+
        scale_color_continuous(type = "viridis")+facet_grid(~type)+
        geom_point(data=df_midpoint, size=1, colour="black")+theme_classic()
dev.off()


gene<-gene2
df<-data.frame(x_coord= c(x_bin[c(stroma_spots, epi_spots)], x_bin[c(stroma_spots, epi_spots)]),
               y_coord= c(-(y_bin[c(stroma_spots, epi_spots)]), -(y_bin[c(stroma_spots, epi_spots)])),
               gene=c(expr_data[gene,c(stroma_spots, epi_spots)], averaged_expr_all_fun[gene,c(stroma_spots, epi_spots)]),
               type=c(rep("original", length(c(stroma_spots, epi_spots))), rep("averaged", length(c(stroma_spots, epi_spots)))))
df$type<-factor(df$type, levels=c("original", "averaged"))

pdf(paste("results/", gene2, "origVSave.pdf", sep="_"),9,5)
ggplot(data=df, aes(x=x_coord, y=y_coord, colour=gene))+geom_point()+
        scale_color_continuous(type = "viridis")+facet_grid(~type)+geom_point(data=df_midpoint, size=1, colour="black")+theme_classic()
dev.off()


pdf("results/ST_COL6A2_EFNA1_cor_inspect_s.pdf", 6, 6)
cor_inspect(md[[1]], gene1, gene2)
dev.off()

pdf("results/ST_HLA-DRB1_ACTG1_cor_inspect_s.pdf", 6, 6)
cor_inspect(md[[1]], gene1="ACTG1", gene2="HLA-DRB1")
dev.off()

pdf("results/ST_LGALS9_ACTG1_cor_inspect_s.pdf", 6, 6)
cor_inspect(md[[1]], gene1="ACTG1", gene2="LGALS9")
dev.off()


#comp1 should be set as the component where the gene was measured
df<-cytoscape_net(adjacency = adj, dataset=as.matrix(md[[1]]), gene="COL6A2", comp1="tis2", comp2="tis1", num=10)
df$Target<-factor(df$Target, levels=df$Target[df$Edge_type=="inter"][order(df$Weight[df$Edge_type=="inter"], decreasing = T)])
pdf("results/ST_COL6A2_top10_corrs_s.pdf", 8,7 )
ggplot(df, aes(x=Target, y=Weight))+geom_col()+facet_grid(Edge_type~.)+theme_classic()
dev.off()


################
#### GO for modules
#################
library(clusterProfiler)
modules<-mods[[2]]$colors

#GO for epi
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
                        write.xlsx( summary(ego_1[[i]]) , paste("Epi", "ST", mod_sel[i], "GO_s.xlsx", sep="_"))
                }
        }
}


##GO for stroma
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
                        write.xlsx( summary(ego_2[[i]]) , paste("Stroma", "ST", mod_sel[i], "GO_s.xlsx", sep="_"))
                }

        }
}


kwithin<-degrees_mod(data=md[[1]], modules=mods[[2]]$colors, Adj_type = "signed",
                     cortype = "spearman",
                     pval = "none",
                     thr = 0.05,
                     beta = 6,
                     comp1 = "_tis1",
                     comp2 = "_tis2")

#weighted module expression
modules<-mods[[2]]$colors

wm<-weighted_mod(modules=mods[[2]]$colors, kwithin, mod_sel=1, averaged_expr_all=averaged_expr_all_fun, comp1="_tis1", comp2="_tis2")


df<-data.frame(x_coord= c(x_bin[epi_spots],x_bin[stroma_spots]),y_coord= c(-(y_bin[epi_spots]), -(y_bin[stroma_spots])),
               compartment=c(rep("epi", length(epi_spots)),rep("stroma", length(stroma_spots))),
               gene=c(wm[[1]][epi_spots], wm[[1]][stroma_spots]))
df_midpoint<-data.frame(x_coord= midpoints[[1]],
                        y_coord= midpoints[[2]])

p1<-ggplot(data=df, aes(x=x_coord, y=y_coord, colour=gene))+geom_point()+
        scale_color_continuous(type = "viridis")+
        geom_point(data=df_midpoint, size=1, colour="black")+theme_classic()


df<-data.frame(x_coord= c(x_bin[epi_spots],x_bin[stroma_spots]),y_coord= c(-(y_bin[epi_spots]), -(y_bin[stroma_spots])),
               compartment=c(rep("epi", length(epi_spots)),rep("stroma", length(stroma_spots))),
               gene=c(wm[[2]][epi_spots], wm[[2]][stroma_spots]))
df_midpoint<-data.frame(x_coord= midpoints[[1]],
                        y_coord= midpoints[[2]])

p2<-ggplot(data=df, aes(x=x_coord, y=y_coord, colour=gene))+geom_point()+
        scale_color_continuous(type = "viridis")+
        geom_point(data=df_midpoint, size=1, colour="black")+theme_classic()



df_midpoint<-data.frame(x_coord= midpoints[[1]],
                        y_coord= midpoints[[2]],
                        comm_score=wm[[1]][sp[[1]]]*wm[[2]][sp[[2]]])


pdf("results/ST_mod1_tis1.pdf", 6,5)
p1
dev.off()

pdf("results/ST_mod1_tis2.pdf", 6,5)
p2
dev.off()

pdf("results/ST_mod1_corr.pdf", 5,5)
plot(wm[[1]][sp[[1]]], wm[[2]][sp[[2]]],pch=19, xlab="weighted expression comp1", ylab="weighted expression comp2")
dev.off()

pdf("results/ST_mod1_comm.pdf", 6,5)
ggplot(data=df, aes(x=x_coord, y=y_coord))+geom_point(data=df_midpoint, size=1, aes(colour=comm_score))+scale_color_continuous(type = "viridis")+theme_classic()
dev.off()




