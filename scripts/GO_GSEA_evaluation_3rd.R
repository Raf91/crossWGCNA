
load("nets_3rdSep.RData")
library(clusterProfiler)
library(xlsx)

names(nets)<-c("GSE5847", "GSE10797", "GSE14548", "GSE83591",
               "GSE68744", "GSE88715")
nets<-nets[c(1,2,4,5,6)]#il 3 eliminato perché troppo piccolo

####Calculate networks for each dataset
for(j in 1:length(nets)){
  modules<-nets[[j]]$merged$colors

background_1<-gsub("_1", "", names(modules)[grep("_1",names(modules))])
background_2<-gsub("_2", "", names(modules)[grep("_2",names(modules))])
mod_sel<-unique(modules)
  ##GO for liver
  ego_1<-list()
  for(i in 1:length(mod_sel)){
    module_1<-gsub("_1", "", names(modules)[modules==mod_sel[i]][grep("_1",names(modules)[modules==mod_sel[i]])])
    if(length(module_1)>0 & mod_sel[i]!=0){
    ego_1[[i]]<- enrichGO(gene = module_1,
                          keyType="SYMBOL",
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          universe=background_1, OrgDb="org.Hs.eg.db")
    
    if(nrow(summary(ego_1[[i]]))>0){
      write.xlsx( summary(ego_1[[i]]) , paste("Stroma", names(nets)[j], mod_sel[i], "GO.xlsx", sep="_"))
    }
    }
  }
  
  
  ##GO for adipose
  ego_2<-list()
  for(i in 1:length(mod_sel)){
    module_2<-gsub("_2", "", names(modules)[modules==mod_sel[i]][grep("_2",names(modules)[modules==mod_sel[i]])])
    if(length(module_2)>0){
    ego_2[[i]]<- enrichGO(gene = module_2,
                          keyType="SYMBOL",
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          universe=background_2, OrgDb="org.Hs.eg.db")
    
    if(nrow(summary(ego_2[[i]]))>0){
      write.xlsx( summary(ego_2[[i]]) , paste("Epi", names(nets)[j], mod_sel[i], "GO.xlsx", sep="_"))
    }
    
  }
  
}
}



  ###GSEA
  library(msigdbr)
  library(fgsea)
  library(data.table)
library(org.Hs.eg.db)
  
  m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
  m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  for(j in 1:length(nets)){
 coef_gsea<-nets[[j]][[1]]$kExt1/nets[[j]][[1]]$kInt1
  names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
  
  fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
  if(nrow(fgseaRes)>0){
  write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]), paste("fgsea_ratio_stroma_", names(nets)[j], ".xlsx", sep=""))
  }
  
  coef_gsea<-nets[[j]][[1]]$kExt2/nets[[j]][[1]]$kInt2
  names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
  
  fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes2<-fgseaRes2[fgseaRes2$padj<0.05,]
  if(nrow(fgseaRes2)>0){
  write.xlsx(as.data.frame(fgseaRes2[,c(1,2,3,6,7)]), paste("fgsea_ratio_epi_", names(nets)[j], ".xlsx", sep=""))
  }
  
  }
  
  png("GSE10797_gsea_stroma.png", res=300, 2200, 1200)
  ggplot(fgseaRes, aes(colour=NES, y=reorder(pathway, -log10(padj)), x=-log10(padj)))+geom_point(size=4)+theme_bw()+ scale_color_gradient(low="blue", high="red")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  dev.off()
  
  fgseaRes2<-fgseaRes2[order(fgseaRes2$padj),]
  fgseaRes2<-fgseaRes2[1:10,]
  
  png("GSE10797_gsea_epi.png", res=300, 2200, 1200)
  ggplot(fgseaRes2, aes(colour=NES, y=reorder(pathway, -log10(padj)), x=-log10(padj)))+geom_point(size=4)+theme_bw()+ scale_color_gradient(low="blue", high="red")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  dev.off()
  
  ###############
  ###REACTOME # non ha finito di girare
  ##########
  
  library(ReactomePA)
  
    head(y)
  
  for(j in 1:length(nets)){
    modules<-nets[[j]][[2]]$merged$colors
    
    background_1<-gsub("_1", "", names(modules)[grep("_1",names(modules))])
    background_2<-gsub("_2", "", names(modules)[grep("_2",names(modules))])
    background_1<-bitr(background_1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    background_2<-bitr(background_2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    mod_sel<-unique(modules)
    
    
    ##GO for liver
    ego_1<-list()
    for(i in 1:length(mod_sel)){
      module_1<-gsub("_1", "", names(modules)[modules==mod_sel[i]][grep("_1",names(modules)[modules==mod_sel[i]])])
      module_1<-bitr(module_1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      if(nrow(module_1)>0){
        ego_1[[i]]<-enrichPathway(gene=module_1[,2], pvalueCutoff = 0.05, readable=TRUE, universe=background_1[,2])
        if(nrow(summary(ego_1[[i]]))>0){
          write.xlsx( summary(ego_1[[i]]) , paste("Stroma", names(nets)[j], mod_sel[i], "reactome.xlsx", sep="_"))
        }
      }
    }
    
    
    ##GO for adipose
    ego_2<-list()
    for(i in 1:length(mod_sel)){
      module_2<-gsub("_2", "", names(modules)[modules==mod_sel[i]][grep("_2",names(modules)[modules==mod_sel[i]])])
      module_2<-bitr(module_2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      if(nrow(module_2)>0){
        ego_2[[i]]<-enrichPathway(gene=module_2[,2], pvalueCutoff = 0.05, readable=TRUE, universe=background_2[,2])
        
        if(nrow(summary(ego_2[[i]]))>0){
          write.xlsx( summary(ego_2[[i]]) , paste("Epi", names(nets)[j], mod_sel[i], "reactome.xlsx", sep="_"))
        }
        
      }
      
    }
  }
  
  ###GSEA
  library(msigdbr)
  library(fgsea)
  library(data.table)
  
   
  for(j in 1:length(nets)){
    coef_gsea<-nets[[j]][[1]]$kExt1/nets[[j]][[1]]$kInt1
    names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
    names(coef_gsea)<-bitr( names(coef_gsea), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)[,2]
    coef_gsea<-sort(coef_gsea, decreasing=T)
    fgseaRes <- gsePathway(coef_gsea, 
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH", 
                    verbose = FALSE)
    
    
    if(nrow(fgseaRes)>0){
      write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]), paste("fgsea_ratio_stroma_", names(nets)[j], "reactome.xlsx", sep=""))
    }
    
    coef_gsea<-nets[[j]][[1]]$kExt2/nets[[j]][[1]]$kInt2
    names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
    names(coef_gsea)<-bitr( names(coef_gsea), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)[,2]
    coef_gsea<-sort(coef_gsea, decreasing=T)
    fgseaRes2 <- gsePathway(coef_gsea, 
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH", 
                           verbose = FALSE)
    
    
    if(nrow(fgseaRes2)>0){
      write.xlsx(as.data.frame(fgseaRes2[,c(1,2,3,6,7)]), paste("fgsea_ratio_epi_", names(nets)[j], "reactome.xlsx", sep=""))
    }
    
  }
  
    head(coef_gsea)
  #######################
  ###Wikipathways
  #######################
  
  data(geneList, package="DOSE")
  gene <- names(geneList)[abs(geneList) > 2]
  
  enrichWP(gene, organism = "Homo sapiens") 
  
  gseWP(geneList, organism = "Homo sapiens")

  
  ###similitudine tra i geni dei moduli
  jac<-list()
  for(j in 1:length(nets)){
    jac[[j]]<-vector()
    modules<-nets[[j]][[2]]$merged$colors
    mod_sel<-sort(unique(modules))[-1]
    for(i in 1:length(mod_sel)){
      module_1<-gsub("_1", "", names(modules)[modules==mod_sel[i]][grep("_1",names(modules)[modules==mod_sel[i]])])
      module_2<-gsub("_2", "", names(modules)[modules==mod_sel[i]][grep("_2",names(modules)[modules==mod_sel[i]])])
      jac[[j]]<-c(jac[[j]], length(intersect(module_1, module_2))/length(union(module_1, module_2)))
    }
  }
  ##il primo dataset ha un'altissima intersezione, mentre gli altri no
  
  signaling<-read.csv("D:/Dropbox (MBC)/Aurora/R analyses/GO_term_summary_20210508_165918.txt", sep="\t", row.names=NULL)
  
  
  ##net 2 CD69, TRIM10
  #net3 VEGFA
  #net4 TNFRSF11A, ID2 (sempre NF-kB)
  ###per selezionare le più importanti, usare il rank product (ma in comune sono pochissimi)
  #un rank più basso (1) rappresenta i geni più importanti per la comunicazione
  
  ranked<-list()
  for(j in 1:length(nets)){
  coef_gsea<-nets[[j]][[1]]$kExt1/nets[[j]][[1]]$kInt1
  names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
   ranked[[j]]<-rank(-coef_gsea)
  }
   
  incommon<-Reduce(intersect, lapply(ranked, names))
  
  rankings<-matrix(nrow=length(incommon), ncol=length(nets))
  rownames(rankings)<-incommon
  for(j in 1:length(nets)){
    rankings[,j]<-ranked[[j]][incommon]
  }
  rankprod<-apply(rankings,1,prod)
  names(rankprod)<-incommon
  
  write.xlsx(rankprod, file="rankprod_stroma.xlsx")
  
  
  #################################################
  ########### Confronti tra network
  ###############################################
  names(nets[[4]][[1]]$kInt1)<-gsub("_1_1", "_1", names(nets[[4]][[1]]$kInt1))
  names(nets[[4]][[1]]$kInt2)<-gsub("_2_2", "_2", names(nets[[4]][[1]]$kInt2))
  names(nets[[4]][[1]]$kExt1)<-gsub("_1_1", "_1", names(nets[[4]][[1]]$kExt1))
  names(nets[[4]][[1]]$kExt2)<-gsub("_2_2", "_2", names(nets[[4]][[1]]$kExt2))
  
  
  ####confronti tra network
 
  corInt1<-matrix(nrow=5, ncol=5)
  for(i in 1:6){
    for(j in 1:6){
      incommon<-intersect(names(nets[[i]][[1]]$kInt1), names(nets[[j]][[1]]$kInt1))
      corInt1[i,j]<-cor(nets[[i]][[1]]$kInt1[incommon], nets[[j]][[1]]$kInt1[incommon], use='pairwise.complete.obs', method="s")
    }
  }
  hist(corInt1)
  pheatmap(corInt1, cluster_rows = F, cluster_cols = F)
  
  corInt2<-matrix(nrow=5, ncol=5)
  for(i in 1:6){
    for(j in 1:6){
      incommon<-intersect(names(nets[[i]][[1]]$kInt2), names(nets[[j]][[1]]$kInt2))
      corInt2[i,j]<-cor(nets[[i]][[1]]$kInt2[incommon], nets[[j]][[1]]$kInt2[incommon], use='pairwise.complete.obs', method="s")
    }
  }
  hist(corInt2)
  
  corExt1<-matrix(nrow=5, ncol=5)
  for(i in 1:6){
    for(j in 1:6){
      incommon<-intersect(names(nets[[i]][[1]]$kExt1), names(nets[[j]][[1]]$kExt1))
      corExt1[i,j]<-cor(nets[[i]][[1]]$kExt1[incommon], nets[[j]][[1]]$kExt1[incommon], use='pairwise.complete.obs', method="s")
    }
  }
  hist(corExt1, breaks=5)
  
  pheatmap(corExt1, cluster_rows = F, cluster_cols = F)
  
  corExt2<-matrix(nrow=6, ncol=6)
  for(i in 1:6){
    for(j in 1:6){
      incommon<-intersect(names(nets[[i]][[1]]$kExt2), names(nets[[j]][[1]]$kExt2))
      corExt2[i,j]<-cor(nets[[i]][[1]]$kExt2[incommon], nets[[j]][[1]]$kExt2[incommon], use='pairwise.complete.obs', method="s")
    }
  }
  hist(corExt2, breaks=20)
  
  names(nets[[4]][[1]]$kInt1)<-gsub("_1", "", names(nets[[4]][[1]]$kInt1))
  names(nets[[4]][[1]]$kInt1)<-paste(names(nets[[4]][[1]]$kInt1), "_1", sep="")
  
  names(nets[[4]][[1]]$kExt1)<-gsub("_1", "", names(nets[[4]][[1]]$kExt1))
  names(nets[[4]][[1]]$kExt1)<-paste(names(nets[[4]][[1]]$kExt1), "_1", sep="")
  
  names(nets[[4]][[1]]$kTot1)<-gsub("_1", "", names(nets[[4]][[1]]$kTot1))
  names(nets[[4]][[1]]$kTot1)<-paste(names(nets[[4]][[1]]$kTot1), "_1", sep="")
  
  names(nets[[4]][[1]]$kInt2)<-gsub("_2", "", names(nets[[4]][[1]]$kInt2))
  names(nets[[4]][[1]]$kInt2)<-paste(names(nets[[4]][[1]]$kInt2), "_2", sep="")
  
  names(nets[[4]][[1]]$kExt2)<-gsub("_2", "", names(nets[[4]][[1]]$kExt2))
  names(nets[[4]][[1]]$kExt2)<-paste(names(nets[[4]][[1]]$kExt2), "_2", sep="")
  
  names(nets[[4]][[1]]$kTot2)<-gsub("_2", "", names(nets[[4]][[1]]$kTot2))
  names(nets[[4]][[1]]$kTot2)<-paste(names(nets[[4]][[1]]$kTot2), "_2", sep="")
  
  
  corRatio1<-matrix(nrow=6, ncol=6)
  corRatio1P<-matrix(nrow=6, ncol=6)
  for(i in 1:6){
    for(j in 1:6){
      incommon<-intersect(names(nets[[i]][[1]]$kExt1), names(nets[[j]][[1]]$kExt1))
      corRatio1[i,j]<-cor(nets[[i]][[1]]$kExt1[incommon]/nets[[i]][[1]]$kInt1[incommon], nets[[j]][[1]]$kExt1[incommon]/nets[[j]][[1]]$kInt1[incommon], use='pairwise.complete.obs', method="s")
      corRatio1P[i,j]<-cor.test(nets[[i]][[1]]$kExt1[incommon]/nets[[i]][[1]]$kInt1[incommon], nets[[j]][[1]]$kExt1[incommon]/nets[[j]][[1]]$kInt1[incommon], use='pairwise.complete.obs', method="s")[[3]]
      
      }
  }
  hist(corRatio1, breaks=10)
  
  corRatio1[corRatio1P>0.05]<-NA
  corRatio1[upper.tri(corRatio1)]
  
  pheatmap(corRatio1, cluster_rows = F, cluster_cols = F)
  
  corRatio2<-matrix(nrow=6, ncol=6)
  corRatio2P<-matrix(nrow=6, ncol=6)
  
  for(i in 1:6){
    for(j in 1:6){
      incommon<-intersect(names(nets[[i]][[1]]$kExt2), names(nets[[j]][[1]]$kExt2))
      corRatio2[i,j]<-cor(nets[[i]][[1]]$kExt2[incommon]/nets[[i]][[1]]$kInt2[incommon], nets[[j]][[1]]$kExt2[incommon]/nets[[j]][[1]]$kInt2[incommon], use='pairwise.complete.obs', method="s")
      corRatio2P[i,j]<-cor.test(nets[[i]][[1]]$kExt2[incommon]/nets[[i]][[1]]$kInt2[incommon], nets[[j]][[1]]$kExt2[incommon]/nets[[j]][[1]]$kInt2[incommon], use='pairwise.complete.obs', method="s")[[3]]
    }
  }
  hist(corRatio2, breaks=10)
  
  corRatio2[corRatio2P>0.05]<-NA
  corRatio2[upper.tri(corRatio2)]
  
  pheatmap(corRatio2, cluster_rows = F, cluster_cols = F)
  
  
  
  library(infotheo)
  all_modules<-matrix(nrow=length(common_all), ncol=6)
  all_modules[,1]<-net_GSE5847[[2]][[1]][common_all]
  all_modules[,2]<-net_GSE10797[[2]][[1]][common_all]
  all_modules[,3]<-net_GSE14548[[2]][[1]][common_all]
  all_modules[,4]<-net_GSE83591[[2]][[1]][common_all]
  all_modules[,5]<-net_GSE68744[[2]][[1]][common_all]
  all_modules[,6]<-net_GSE88715[[2]][[1]][common_all]
  all_modules<-as.data.frame(all_modules)
  MI<-mutinformation(X=(all_modules))
  
  
  
  library(aricode)
  Nmut<-matrix(ncol=6, nrow=6)
  for(i in 1:6){
    for(j in 1:6){
      incommon<-c(intersect(names(nets[[i]][[1]]$kExt1), names(nets[[j]][[1]]$kExt1)), intersect(names(nets[[i]][[1]]$kExt2), names(nets[[j]][[1]]$kExt2)))
      
      Nmut[i,j]<-NMI(nets[[i]][[2]]$merged$colors[incommon],nets[[j]][[2]]$merged$colors[incommon])
    }
  }
  
  png("NMI.png")
  pheatmap(Nmut, cluster_rows = F, cluster_cols = F)
  dev.off()
  
  MI<-mutinformation(X=(all_modules))
  
  
  #####Confronto con single cell dataset
  ranked<-list()
  for(j in 1:length(nets)){
    coef_gsea<-nets[[j]][[1]]$kExt1/nets[[j]][[1]]$kInt1
    names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
    ranked[[j]]<-rank(-coef_gsea)
  }
  
  incommon<-Reduce(intersect, lapply(ranked, names))
  
  rankings<-matrix(nrow=length(incommon), ncol=length(nets))
  rownames(rankings)<-incommon
  for(j in 1:length(nets)){
    rankings[,j]<-ranked[[j]][incommon]
  }
  rankprod<-apply(rankings,1,prod)
  names(rankprod)<-incommon
  
  write.xlsx(rankprod, file="rankprod_stroma.xlsx")
  
  
  load("F:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/netsc.RData")
  singlecell<-netsc[[1]]$kExt1/netsc[[1]]$kInt1
  for(j in 1:length(nets)){
    coef_gsea<-nets[[j]][[1]]$kExt1/nets[[j]][[1]]$kInt1
    plot(coef_gsea, singlecell[match(names(coef_gsea), names(singlecell))])
    
  }
  
  plot(log2(rankprod), singlecell[match(names(rankprod), gsub("_1", "", names(singlecell)))])
  
  which(log2(rankprod)<60 & singlecell[match(names(rankprod), gsub("_1", "", names(singlecell)))] >1)
  
  df<-data.frame(logFC1=log2(rankprod), logFC2=singlecell[match(names(rankprod), gsub("_1", "", names(singlecell)))])
  
  png("Singlecell_scatter_stroma.png", res=300, 2000, 1500)
  ggplot(df, aes(x=logFC1, y=logFC2))+geom_hex(bins=70)+xlab("Rank product")+ylab("Comm. score single cell")+theme_bw()+geom_smooth(method = "lm", col="black")
  dev.off()
  cor.test(df$logFC1, df$logFC2)
  
  
  top<-names(sort(rankprod)[1:100])
  bottom<-names(sort(rankprod, decreasing=T)[1:100])
  nontop<-names(sort(rankprod)[-c(1:100)])
  
  boxplot(singlecell[match(top, gsub("_1", "", names(singlecell)))], singlecell[match(bottom, gsub("_1", "", names(singlecell)))])
  boxplot(singlecell[match(top, gsub("_1", "", names(singlecell)))], singlecell[match(nontop, gsub("_1", "", names(singlecell)))])
  
  df<-data.frame(score=c(singlecell[match(top, gsub("_1", "", names(singlecell)))], singlecell[match(bottom, gsub("_1", "", names(singlecell)))]), class=rep(c("Top", "Bottom"), each=100))
  df$class<-factor(df$class, levels=c("Top", "Bottom"))
  
  png("Singlecell_boxplot_stroma.png", res=300, 700, 1500)
  ggplot(df, aes(x=class, y=score))+geom_boxplot()+geom_signif(comparisons = list(c("Top", "Bottom")), map_signif_level=TRUE) + theme_bw()+
    theme(axis.text.x = element_text(angle = 90))+ylab("Comm. score single cell")+xlab("")
  dev.off()
  
  ############################
  ### Epi score
  
  ranked<-list()
  for(j in 1:length(nets)){
    coef_gsea<-nets[[j]][[1]]$kExt2/nets[[j]][[1]]$kInt2
    names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
    ranked[[j]]<-rank(-coef_gsea)
  }
  
  incommon<-Reduce(intersect, lapply(ranked, names))
  
  rankings<-matrix(nrow=length(incommon), ncol=length(nets))
  rownames(rankings)<-incommon
  for(j in 1:length(nets)){
    rankings[,j]<-ranked[[j]][incommon]
  }
  rankprod<-apply(rankings,1,prod)
  names(rankprod)<-incommon
  
  write.xlsx(rankprod, file="rankprod_epi.xlsx")
  
  
  singlecell<-netsc[[1]]$kExt2/netsc[[1]]$kInt2
  for(j in 1:length(nets)){
    coef_gsea<-nets[[j]][[1]]$kExt2/nets[[j]][[1]]$kInt2
    plot(coef_gsea, singlecell[match(names(coef_gsea), names(singlecell))])
    
  }
  
  plot(log2(rankprod), singlecell[match(names(rankprod), gsub("_2", "", names(singlecell)))])
  
  df<-data.frame(logFC1=log2(rankprod), logFC2=singlecell[match(names(rankprod), gsub("_2", "", names(singlecell)))])
  png("Singlecell_scatter_epi.png", res=300, 2000, 1500)
  ggplot(df, aes(x=logFC1, y=logFC2))+geom_hex(bins=70)+xlab("Rank product")+ylab("Comm. score single cell")+theme_bw()+geom_smooth(method = "lm", col="black")
  dev.off()
  cor.test(df$logFC1, df$logFC2)
  
  
  top<-names(sort(rankprod)[1:100])
  bottom<-names(sort(rankprod, decreasing=T)[1:100])
  nontop<-names(sort(rankprod)[-c(1:100)])
  
  df<-data.frame(score=c(singlecell[match(top, gsub("_2", "", names(singlecell)))], singlecell[match(bottom, gsub("_2", "", names(singlecell)))]), class=rep(c("Top", "Bottom"), each=100))
  df$class<-factor(df$class, levels=c("Top", "Bottom"))
  
  png("Singlecell_boxplot_epi.png", res=300, 700, 1500)
  ggplot(df, aes(x=class, y=score))+geom_boxplot()+geom_signif(comparisons = list(c("Top", "Bottom")), map_signif_level=TRUE) + theme_bw()+
    theme(axis.text.x = element_text(angle = 90))+ylab("Comm. score single cell")+xlab("")
  dev.off()
  
  
  
  #######################
  ######## Geni con cui CP è correlato
  ######################
  
  adj_GSE5847<-Adjacency(data=data_merged_GSE5847, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6)
  CP_GSE5847<-adj_GSE5847["CP_1",]
  adj_GSE10797<-Adjacency(data=data_merged_GSE10797, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6)
  CP_GSE10797<-adj_GSE10797["CP_1",]
  
  adj_GSE14548<-Adjacency(data=data_merged_GSE14548, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6)
  CP_GSE14548<-adj_GSE14548["CP_1",]
  adj_GSE83591<-Adjacency(data=data_merged_GSE83591, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6)
  CP_GSE83591<-adj_GSE83591["CP_1",]
  
  adj_GSE68744<-Adjacency(data=data_merged_GSE68744, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6)
  adj_GSE88715<-Adjacency(data=data_merged_GSE88715, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6)
  
  incommon<-Reduce(intersect, list(names(CP_GSE5847), names(CP_GSE10797), names(CP_GSE14548),
                 names(CP_GSE83591), names(CP_GSE68744), names(CP_GSE88715)))
  
  incommon<-incommon[grep("_2", incommon)]
  
  CP_corr<-cbind(CP_GSE5847[incommon], CP_GSE10797[incommon], CP_GSE14548[incommon],
  CP_GSE83591[incommon], CP_GSE68744[incommon], CP_GSE88715[incommon])
 rownames(CP_corr)<-gsub("_2", "", incommon)
 CP_corr<-cbind(CP_corr, rowMeans(CP_corr))
 colnames(CP_corr)<-c("GSE5847", "GSE10797", "GSE14548", "GSE83591", "GSE68744", "GSE88715", "average")
 
 write.csv(CP_corr, file="CP_corr_LCM.csv")
 
 library(msigdbr)
 library(fgsea)
 library(data.table)
 library(org.Hs.eg.db)
 
 m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
 m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
 
 
   coef_gsea<-CP_corr[,7]
   names(coef_gsea)<-rownames(CP_corr)
   
   fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
   fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
   write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]),"fgsea_CP_corr_epi.xlsx")
   
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/data_merged_EpiStroma.RData")
adj_sc<-Adjacency(data=data_merged, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6)
CP_sc<-adj_sc["CP_1",grep("_2", colnames(adj_sc))]
names(CP_sc)<-gsub("_2", "", names(CP_sc))
CP_sc<-data.frame(gene=names(CP_sc), Adj=CP_sc)
CP_sc<-CP_sc[order(CP_sc[,2], decreasing = T),]
write.xlsx(CP_sc, file="CP_corr_sc.xlsx", row.names = F)
CP_corr<-read.xlsx("CP_corr_LCM.xlsx", 1)

rankings<-matrix(nrow=nrow(CP_corr), ncol=ncol(CP_corr)-2)
rownames(rankings)<-CP_corr[,1]
for(j in 2:ncol(CP_corr)-1){
  rankings[,j-1]<-rank(CP_corr[,j])
}
rankprod<-apply(rankings,1,prod)
names(rankprod)<-rownames(rankings)

rownames(CP_sc)<-CP_sc[,1]
inboth<-intersect(CP_sc[,1], CP_corr[,1])
plot(CP_sc[match(inboth, CP_sc[,1]),2], rankprod[inboth])


m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)


coef_gsea<-rankprod
names(coef_gsea)<-names(rankprod)

fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]),"fgsea_CP_corr_epi_rankprod.xlsx")

coef_gsea<-CP_sc[,2]
names(coef_gsea)<-CP_sc[,1]

fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]),"fgsea_CP_corr_epi_sc.xlsx")
