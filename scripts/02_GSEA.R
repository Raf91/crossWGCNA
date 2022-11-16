
load("results/degs_3rd.RData")
library(clusterProfiler)
library(openxlsx)
names(degs)<-c("GSE5847", "GSE10797", "GSE14548", "GSE83591",
               "GSE68744", "GSE88715")


   ###GSEA
  library(msigdbr)
  library(fgsea)
  library(data.table)
library(org.Hs.eg.db)

  m_df = msigdbr(species = "Homo sapiens", category = "C5")
  m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

  for(j in 1:length(degs)){
 coef_gsea<-degs[[j]]$kExt1/degs[[j]]$kInt1
  names(coef_gsea)<-gsub("_tis1", "", names(coef_gsea))

  fgseaRes <- fgseaMultilevel(m_list, coef_gsea)
  fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
  if(nrow(fgseaRes)>0){
  write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]), paste("fgsea_ratio_stroma_", names(degs)[j], ".xlsx", sep=""))
  }

  coef_gsea<-degs[[j]]$kExt2/degs[[j]]$kInt2
  names(coef_gsea)<-gsub("_tis2", "", names(coef_gsea))

  fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea)
  fgseaRes2<-fgseaRes2[fgseaRes2$padj<0.05,]
  if(nrow(fgseaRes2)>0){
  write.xlsx(as.data.frame(fgseaRes2[,c(1,2,3,6,7)]), paste("fgsea_ratio_epi_", names(degs)[j], ".xlsx", sep=""))
  }
  }



  m_df = msigdbr(species = "Homo sapiens", category = "C5")
  m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

  fgseaRes<-list()
  fgseaRes2<-list()
  for(j in 1:length(degs)){
    coef_gsea<-degs[[j]]$kExt1/degs[[j]]$kInt1
    names(coef_gsea)<-gsub("_tis1", "", names(coef_gsea))

    fgseaRes[[j]] <- fgseaMultilevel(m_list, coef_gsea)

    coef_gsea<-degs[[j]]$kExt2/degs[[j]]$kInt2
    names(coef_gsea)<-gsub("_tis2", "", names(coef_gsea))

    fgseaRes2[[j]] <- fgseaMultilevel(m_list, coef_gsea)

  }




  #############################
  ####### stroma
  ###############################

  allpaths<-c()
  for(i in 1:length(fgseaRes)){
    allpaths<-union(allpaths, fgseaRes[[i]]$pathway[fgseaRes[[i]]$pval<0.05])
  }
  allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgseaRes))
  rownames(allpaths_mat)<-allpaths
  for(i in 1:length(fgseaRes)){
    allpaths_mat[fgseaRes[[i]]$pathway[which(fgseaRes[[i]]$pval<0.05)],i]<-1
  }


  library(pheatmap)
  shared_paths<-allpaths_mat[rowSums(allpaths_mat)>3, ]
  allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgseaRes))
  rownames(allpaths_mat)<-allpaths
  for(i in 1:length(fgseaRes)){
    allpaths_mat[fgseaRes[[i]]$pathway[which(fgseaRes[[i]]$pval<0.05)],i]<-fgseaRes[[i]]$NES[which(fgseaRes[[i]]$pval<0.05)]
  }

   hist(rowSums(sign(allpaths_mat)))

  toplot<-allpaths_mat[rowSums(sign(allpaths_mat))>2,]
  #toplot<-allpaths_mat[rowSums(sign(allpaths_mat))<(-2),]
  #rownames(toplot)<-ego_up$Description[match(rownames(toplot), ego_up$ID)]
  paletteLength <- 50
  myColor <- colorRampPalette(c( "white", "red"))(paletteLength)
  myBreaks <- c(seq(0, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength)))
  length(myBreaks) == length(paletteLength) + 1
  pheatmap(toplot,   cellwidth=8, cellheight=8, breaks=myBreaks, color = myColor, keep.dendro=T)

  pdf("results/summary_GSEA_stroma.pdf", 10, 15)
  pheatmap(toplot,   cellwidth=10, cellheight=10, breaks=myBreaks, color = myColor, keep.dendro=T)
  dev.off()

  #############################
  ####### epi
  ###############################

  allpaths<-c()
  for(i in 1:length(fgseaRes2)){
    allpaths<-union(allpaths, fgseaRes2[[i]]$pathway[fgseaRes2[[i]]$pval<0.05])
  }
  allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgseaRes2))
  rownames(allpaths_mat)<-allpaths
  for(i in 1:length(fgseaRes2)){
    allpaths_mat[fgseaRes2[[i]]$pathway[which(fgseaRes2[[i]]$pval<0.05)],i]<-1
  }


  library(pheatmap)
  shared_paths<-allpaths_mat[rowSums(allpaths_mat)>3, ]
  allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgseaRes2))
  rownames(allpaths_mat)<-allpaths
  for(i in 1:length(fgseaRes2)){
    allpaths_mat[fgseaRes2[[i]]$pathway[which(fgseaRes2[[i]]$pval<0.05)],i]<-fgseaRes2[[i]]$NES[which(fgseaRes2[[i]]$pval<0.05)]
  }

  hist(rowSums(sign(allpaths_mat)))

  toplot<-allpaths_mat[rowSums(sign(allpaths_mat))>3,]
  #toplot<-allpaths_mat[rowSums(sign(allpaths_mat))<(-2),]
  #rownames(toplot)<-ego_up$Description[match(rownames(toplot), ego_up$ID)]
  paletteLength <- 50
  myColor <- colorRampPalette(c( "white", "red"))(paletteLength)
  myBreaks <- c(seq(0, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength)))
  length(myBreaks) == length(paletteLength) + 1
  pheatmap(toplot,   cellwidth=5, cellheight=5, breaks=myBreaks, color = myColor, keep.dendro=T)

  pdf("results/summary_GSEA_epi.pdf", 10, 15)
  pheatmap(toplot,   cellwidth=10, cellheight=10, breaks=myBreaks, color = myColor, keep.dendro=T)
  dev.off()

