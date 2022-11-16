rm(list=ls())
library(openxlsx)
library(ggplot2)
library(ggpubr)
load("results/degs_3rd.RData")
names(degs)<-c("GSE5847", "GSE10797", "GSE14548", "GSE83591",
               "GSE68744", "GSE88715")
load("degsc.RData")

library(biomaRt)
ensembl.human<- useEnsembl(biomart = 'genes', 
                           dataset = 'hsapiens_gene_ensembl')
#extracellular region
GO<-getBM(attributes = c('hgnc_symbol'),
          filters = 'go', 
          values = "GO:0005576", 
          mart = ensembl.human)

extracellular<-GO[,1]


########STROMA

###rank product dei 6 network
ranked<-list()
for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt1/degs[[j]]$kInt1
  names(coef_gsea)<-gsub("_tis1", "", names(coef_gsea))
  ranked[[j]]<-rank(-coef_gsea)
}

incommon<-Reduce(intersect, lapply(ranked, names))

rankings<-matrix(nrow=length(incommon), ncol=length(degs))
rownames(rankings)<-incommon
for(j in 1:length(degs)){
  rankings[,j]<-ranked[[j]][incommon]
}
rankprod<-apply(rankings,1,prod)
names(rankprod)<-incommon

write.xlsx(rankprod, file="results/rankprod_stroma.xlsx")

sort(rankprod[intersect(names(rankprod), extracellular)], decreasing=F)[1:10]

coef_gsea_sc<-degsc$kExt1/degsc$kInt1
names(coef_gsea_sc)<-gsub("_1", "", names(coef_gsea_sc))

#confronto con rankprod
inboth<-intersect(names(rankprod), names(coef_gsea_sc))
cor.test(coef_gsea_sc[inboth], log2(rankprod[inboth]))

df<-data.frame(rankproduct=log2(rankprod[inboth]), sc=coef_gsea_sc[inboth])
pdf("results/singlecell_scatter_stroma.pdf",5,4)
ggplot(df, aes(x=rankproduct,y=sc))+geom_hex(bins=60)+ylab(label = "Comm. score single cell")+xlab(label = "Rank product")+geom_smooth(method="lm")+theme_classic()
dev.off()

library(ggpubr)
rankprodinboth<-rankprod[inboth]
top<-names(rankprodinboth)[which(rankprodinboth<quantile(rankprodinboth, c(0.1)))]
bottom<-names(rankprodinboth)[which(rankprodinboth>quantile(rankprodinboth, c(0.9)))]
df<-data.frame(sc=c(coef_gsea_sc[top], coef_gsea_sc[bottom]), dir=factor(c(rep("Top", length(top)), rep("Bottom", length(bottom)))))

my_comparisons <- list( c("Top", "Bottom"))
pdf("results/singlecell_boxplot_stroma.pdf",3,5)
ggboxplot(df, x="dir", y="sc")+stat_compare_means(comparisons=my_comparisons,label = "p.signif")+ylab("Comm. score single cell")+xlab("")
dev.off()

top50LCM<-names(sort(rankprod, decreasing=F))[1:100]
top50sc<-names(sort(coef_gsea_sc, decreasing=T))[1:100]
intersect(top50LCM, top50sc)


####################EPI
###rank product dei 6 network
ranked<-list()
for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt2/degs[[j]]$kInt2
  names(coef_gsea)<-gsub("_tis2", "", names(coef_gsea))
  ranked[[j]]<-rank(-coef_gsea)
}

incommon<-Reduce(intersect, lapply(ranked, names))

rankings<-matrix(nrow=length(incommon), ncol=length(degs))
rownames(rankings)<-incommon
for(j in 1:length(degs)){
  rankings[,j]<-ranked[[j]][incommon]
}
rankprod<-apply(rankings,1,prod)
names(rankprod)<-incommon

write.xlsx(rankprod, file="rankprod_epi.xlsx")
sort(rankprod[intersect(names(rankprod), extracellular)], decreasing=F)[1:10]

coef_gsea_sc<-degsc$kExt2/degsc$kInt2
names(coef_gsea_sc)<-gsub("_2", "", names(coef_gsea_sc))

#confronto con rankprod
inboth<-intersect(names(rankprod), names(coef_gsea_sc))
cor.test(coef_gsea_sc[inboth], log2(rankprod[inboth]))

df<-data.frame(rankproduct=log2(rankprod[inboth]), sc=coef_gsea_sc[inboth])
pdf("results/singlecell_scatter_epi.pdf",5,4)
ggplot(df, aes(x=rankproduct,y=sc))+geom_hex(bins=60)+ylab(label = "Comm. score single cell")+xlab(label = "Rank product")+geom_smooth(method="lm")+theme_classic()
dev.off()

library(ggpubr)
rankprodinboth<-rankprod[inboth]
top<-names(rankprodinboth)[which(rankprodinboth<quantile(rankprodinboth, c(0.1)))]
bottom<-names(rankprodinboth)[which(rankprodinboth>quantile(rankprodinboth, c(0.9)))]
df<-data.frame(sc=c(coef_gsea_sc[top], coef_gsea_sc[bottom]), dir=factor(c(rep("Top", length(top)), rep("Bottom", length(bottom)))))

my_comparisons <- list( c("Top", "Bottom"))
pdf("results/singlecell_boxplot_epi.pdf",3,5)
ggboxplot(df, x="dir", y="sc")+stat_compare_means(comparisons=my_comparisons,label = "p.signif")+ylab("Comm. score single cell")+xlab("")
dev.off()


top50LCM<-names(sort(rankprod, decreasing=F))[1:500]
top50sc<-names(sort(coef_gsea_sc, decreasing=T))[1:500]
intersect(top50LCM, top50sc)

###GSEA
library(msigdbr)
library(fgsea)
library(data.table)
library(openxlsx)

m_df = msigdbr(species = "Homo sapiens", category = "C5")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)


  coef_gsea<-degsc$kExt1/degsc$kInt1
  names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
  
  fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
  if(nrow(fgseaRes)>0){
    write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]), paste("results/fgsea_ratio_stroma_sc.xlsx", sep=""))
  }
  
  coef_gsea<-degsc$kExt2/degsc$kInt2
  names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
  
  fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
  fgseaRes2<-fgseaRes2[fgseaRes2$padj<0.05,]
  if(nrow(fgseaRes2)>0){
    write.xlsx(as.data.frame(fgseaRes2[,c(1,2,3,6,7)]), paste("results/fgsea_ratio_epi_sc.xlsx", sep=""))
}

  
