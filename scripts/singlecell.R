load("F:/Dropbox (MBC)/Aurora/R analyses/Double WGCNA/DoubleWGCNA stroma-epi/Systematic evaluation/thirdquartile bis/nets_3rd.RData")
load("F:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/netsc.RData")

###rank product dei 6 network
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

#write.xlsx(rankprod, file="rankprod_stroma.xlsx")

coef_gsea_sc<-netsc[[1]]$kExt1/netsc[[1]]$kInt1
names(coef_gsea_sc)<-gsub("_1", "", names(coef_gsea_sc))

#confronto con rankprod
inboth<-intersect(names(rankprod), names(coef_gsea_sc))
plot(coef_gsea_sc[inboth], log2(rankprod[inboth]))
cor.test(coef_gsea_sc[inboth], log2(rankprod[inboth]))

which(coef_gsea_sc[inboth]>1&(log2(rankprod[inboth])<55))

#confronto con singoli network
j<-5
coef_gsea<-nets[[j]][[1]]$kExt1/nets[[j]][[1]]$kInt1
names(coef_gsea)<-gsub("_1", "", names(coef_gsea))

inboth<-intersect(names(coef_gsea), names(coef_gsea_sc))
plot(coef_gsea_sc[inboth], (coef_gsea[inboth]))
cor.test(coef_gsea_sc[inboth], (coef_gsea[inboth]))

top6<-names(which(coef_gsea_sc[inboth]>1&(coef_gsea[inboth]>1)))

sort(table(c(top1,top2,top3,top4,top5,top6)), decreasing=T)


quantile(rankprod, c(0.1, 0.9))

rankprodinboth<-rankprod[inboth]

top<-names(rankprodinboth)[which(rankprodinboth<quantile(rankprodinboth, c(0.1)))]
bottom<-names(rankprodinboth)[which(rankprodinboth>quantile(rankprodinboth, c(0.9)))]

boxplot(coef_gsea_sc[top], coef_gsea_sc[bottom])


###GSEA
library(msigdbr)
library(fgsea)
library(data.table)
library(openxlsx)

m_df = msigdbr(species = "Homo sapiens", category = "C5")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)


  coef_gsea<-netsc[[1]]$kExt1/netsc[[1]]$kInt1
  names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
  
  fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
  if(nrow(fgseaRes)>0){
    write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]), paste("fgsea_ratio_stroma_sc.xlsx", sep=""))
  }
  
  coef_gsea<-netsc[[1]]$kExt2/netsc[[1]]$kInt2
  names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
  
  fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes2<-fgseaRes2[fgseaRes2$padj<0.05,]
  if(nrow(fgseaRes2)>0){
    write.xlsx(as.data.frame(fgseaRes2[,c(1,2,3,6,7)]), paste("fgsea_ratio_epi_sc.xlsx", sep=""))
  
  
}

  
  load("F:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/netsc_ET.RData")
  
  m_df = msigdbr(species = "Homo sapiens", category = "C5")
  m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  
  coef_gsea<-netsc_ET[[1]]$kExt1/netsc_ET[[1]]$kInt1
  names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
  
  fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
  if(nrow(fgseaRes)>0){
    write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]), paste("fgsea_ratio_stroma_scET.xlsx", sep=""))
  }
  
  coef_gsea<-netsc_ET[[1]]$kExt2/netsc_ET[[1]]$kInt2
  names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
  
  fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes2<-fgseaRes2[fgseaRes2$padj<0.05,]
  if(nrow(fgseaRes2)>0){
    write.xlsx(as.data.frame(fgseaRes2[,c(1,2,3,6,7)]), paste("fgsea_ratio_epi_scET.xlsx", sep=""))
    
    
  }
  
  
  
  load("F:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/netsc_EE.RData")
  
  m_df = msigdbr(species = "Homo sapiens", category = "C5")
  m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  
  coef_gsea<-netsc_EE[[1]]$kExt1/netsc_EE[[1]]$kInt1
  names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
  
  fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
  if(nrow(fgseaRes)>0){
    write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]), paste("fgsea_ratio_stroma_scEE.xlsx", sep=""))
  }
  
  coef_gsea<-netsc_EE[[1]]$kExt2/netsc_EE[[1]]$kInt2
  names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
  
  fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes2<-fgseaRes2[fgseaRes2$padj<0.05,]
  if(nrow(fgseaRes2)>0){
    write.xlsx(as.data.frame(fgseaRes2[,c(1,2,3,6,7)]), paste("fgsea_ratio_epi_scEE.xlsx", sep=""))
    
    
  }
  
  
  
  load("F:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/netsc_EB.RData")
  
  m_df = msigdbr(species = "Homo sapiens", category = "C5")
  m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  
  coef_gsea<-netsc_EB[[1]]$kExt1/netsc_EB[[1]]$kInt1
  names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
  
  fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
  if(nrow(fgseaRes)>0){
    write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]), paste("fgsea_ratio_stroma_scEB.xlsx", sep=""))
  }
  
  coef_gsea<-netsc_EB[[1]]$kExt2/netsc_EB[[1]]$kInt2
  names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
  
  fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes2<-fgseaRes2[fgseaRes2$padj<0.05,]
  if(nrow(fgseaRes2)>0){
    write.xlsx(as.data.frame(fgseaRes2[,c(1,2,3,6,7)]), paste("fgsea_ratio_epi_scEB.xlsx", sep=""))
    
    
  }
  
  
  load("F:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/netsc_EM.RData")
  
  m_df = msigdbr(species = "Homo sapiens", category = "C5")
  m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  
  coef_gsea<-netsc_EM[[1]]$kExt1/netsc_EM[[1]]$kInt1
  names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
  
  fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
  if(nrow(fgseaRes)>0){
    write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]), paste("fgsea_ratio_stroma_scEM.xlsx", sep=""))
  }
  
  coef_gsea<-netsc_EM[[1]]$kExt2/netsc_EM[[1]]$kInt2
  names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
  
  fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes2<-fgseaRes2[fgseaRes2$padj<0.05,]
  if(nrow(fgseaRes2)>0){
    write.xlsx(as.data.frame(fgseaRes2[,c(1,2,3,6,7)]), paste("fgsea_ratio_epi_scEM.xlsx", sep=""))
    
    
  }
  
  
  load("F:/Dropbox (MBC)/Aurora/R analyses/Data/Breast cancer datasets/single cell/GSE161529/netsc_EP.RData")
  
  m_df = msigdbr(species = "Homo sapiens", category = "C5")
  m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  
  coef_gsea<-netsc_EP[[1]]$kExt1/netsc_EP[[1]]$kInt1
  names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
  
  fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
  if(nrow(fgseaRes)>0){
    write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]), paste("fgsea_ratio_stroma_scEP.xlsx", sep=""))
  }
  
  coef_gsea<-netsc_EP[[1]]$kExt2/netsc_EP[[1]]$kInt2
  names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
  
  fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes2<-fgseaRes2[fgseaRes2$padj<0.05,]
  if(nrow(fgseaRes2)>0){
    write.xlsx(as.data.frame(fgseaRes2[,c(1,2,3,6,7)]), paste("fgsea_ratio_epi_scEP.xlsx", sep=""))
    
    
  }