load("F:/Dropbox (MBC)/Aurora/R analyses/Metanalyses/cocultures/cocultures_CAFsonly.RData")
load("F:/Dropbox (MBC)/Aurora/R analyses/Double WGCNA/DoubleWGCNA stroma-epi/Systematic evaluation/thirdquartile bis/nets_3rdSep.RData")
library(clusterProfiler)
library(xlsx)

names(nets)<-c("GSE5847", "GSE10797", "GSE14548", "GSE83591",
               "GSE68744", "GSE88715")




#Creare una lista DE con tutti i geni differenziali nei vari confronti


enrich_coculture<-function(background_data, modules, DEGs, comp="_1"){
  
  ngenes_stroma<-length(grep(comp, background_data))
  modules_stroma<-modules[grep(comp, background_data)]
  names(modules_stroma)<-gsub(comp, "",names(modules_stroma))
  
  pval_DEGsstroma<-c()
  for(i in 1:max(modules)){
    a<-length(intersect(names(modules_stroma)[modules_stroma==i],DEGs))
    b<-length(intersect(names(modules_stroma)[modules_stroma!=i&modules_stroma!=0],DEGs))
    c<-length(setdiff(names(modules_stroma)[modules_stroma==i], DEGs))
    d<-length(setdiff(names(modules_stroma)[modules_stroma!=i&modules_stroma!=0],DEGs))
    pval_DEGsstroma<-c(pval_DEGsstroma, fisher.test(matrix(c(a,b,c,d), byrow = T, nrow=2), alternative = "greater")$p.value)
  }
  
  return(pval_DEGsstroma)
}

enrich_epi<-list()
for(i in 1:length(DE_epi)){
  enrich_epi[[i]]<-list()
  for(n in 1:length(nets)){
    enrich_epi[[i]][[n]]<-enrich_coculture(background_data=names(nets[[n]][[2]]$merged$colors), modules=nets[[n]][[2]]$merged$colors, DEGs=DE_epi[[i]], comp="_2")
  }
}


enrich_stroma<-list()
for(i in 1:length(DE_stroma)){
  enrich_stroma[[i]]<-list()
  for(n in 1:length(nets)){
    enrich_stroma[[i]][[n]]<-enrich_coculture(background_data=names(nets[[n]][[2]]$merged$colors), modules=nets[[n]][[2]]$merged$colors, DEGs=DE_stroma[[i]], comp="_1")
  }
}



#ci sono alcuni specifici moduli sempre arricchiti?
enrich_mat1<-enrich_epi[[1]][[1]]
for(i in 2:length(enrich_epi)){
  enrich_mat1<-cbind(enrich_mat1, enrich_epi[[i]][[1]])
}

enrich_mat2<-enrich_epi[[1]][[2]]
for(i in 2:length(enrich_epi)){
  enrich_mat2<-cbind(enrich_mat2, enrich_epi[[i]][[2]])
}

enrich_mat3<-enrich_epi[[1]][[3]]
for(i in 2:length(enrich_epi)){
  enrich_mat3<-cbind(enrich_mat3, enrich_epi[[i]][[3]])
}

enrich_mat4<-enrich_epi[[1]][[4]]
for(i in 2:length(enrich_epi)){
  enrich_mat4<-cbind(enrich_mat4, enrich_epi[[i]][[4]])
}

enrich_mat5<-enrich_epi[[1]][[5]]
for(i in 2:length(enrich_epi)){
  enrich_mat5<-cbind(enrich_mat5, enrich_epi[[i]][[5]])
}

enrich_mat6<-enrich_epi[[1]][[6]]
for(i in 2:length(enrich_epi)){
  enrich_mat6<-cbind(enrich_mat6, enrich_epi[[i]][[6]])
}

rowSums(enrich_mat1<0.05)
rowSums(enrich_mat2<0.05)
rowSums(enrich_mat3<0.05)
rowSums(enrich_mat4<0.05)
rowSums(enrich_mat5<0.05)
rowSums(enrich_mat6<0.05)

dim(enrich_mat6<0.05)
sum(enrich_mat6<0.05)

#ci sono alcuni specifici moduli sempre arricchiti?
enrich_mat1s<-enrich_stroma[[1]][[1]]
for(i in 2:length(enrich_stroma)){
  enrich_mat1s<-cbind(enrich_mat1s, enrich_stroma[[i]][[1]])
}

enrich_mat2s<-enrich_stroma[[1]][[2]]
for(i in 2:length(enrich_stroma)){
  enrich_mat2s<-cbind(enrich_mat2s, enrich_stroma[[i]][[2]])
}

enrich_mat3s<-enrich_stroma[[1]][[3]]
for(i in 2:length(enrich_stroma)){
  enrich_mat3s<-cbind(enrich_mat3s, enrich_stroma[[i]][[3]])
}

enrich_mat4s<-enrich_stroma[[1]][[4]]
for(i in 2:length(enrich_stroma)){
  enrich_mat4s<-cbind(enrich_mat4s, enrich_stroma[[i]][[4]])
}

enrich_mat5s<-enrich_stroma[[1]][[5]]
for(i in 2:length(enrich_stroma)){
  enrich_mat5s<-cbind(enrich_mat5s, enrich_stroma[[i]][[5]])
}

enrich_mat6s<-enrich_stroma[[1]][[6]]
for(i in 2:length(enrich_stroma)){
  enrich_mat6s<-cbind(enrich_mat6s, enrich_stroma[[i]][[6]])
}

plot(rowSums(enrich_mat1<0.05), rowSums(enrich_mat1s<0.05))
plot(rowSums(enrich_mat2<0.05), rowSums(enrich_mat2s<0.05))
forphe<-rbind(rowSums(enrich_mat2<0.05)/ncol(enrich_mat2), rowSums(enrich_mat2s<0.05)/ncol(enrich_mat2s))
colnames(forphe)<-1:ncol(forphe)
rownames(forphe)<-c("Epi", "Stroma")

png("GSE10797_numberenrich_coculturebis.png", res=300, 2200, 300)
pheatmap(forphe, cluster_cols = F, cluster_rows = F)
dev.off()

png("GSE10797_numberenrich_coculturebis2.png", res=300, 2200, 1000)
pheatmap(forphe, cluster_cols = F, cluster_rows = F)
dev.off()

plot(rowSums(enrich_mat3<0.05), rowSums(enrich_mat3s<0.05))
plot(rowSums(enrich_mat4<0.05), rowSums(enrich_mat4s<0.05))
plot(rowSums(enrich_mat5<0.05), rowSums(enrich_mat5s<0.05))
plot(rowSums(enrich_mat6<0.05), rowSums(enrich_mat6s<0.05))


rowSums(enrich_mat2<0.05)


###############
##Quanti test complessivamente significativi?
##############
sum(unlist(enrich_epi)<0.05)/length(unlist(enrich_epi))
sum(unlist(enrich_stroma)<0.05)/length(unlist(enrich_stroma))

#i geni DE sono tra i top interacting? (fare un gsea)

###GSEA
library(msigdbr)
library(fgsea)
library(data.table)

names(DE_stroma)<-c(1:20)
gsea_stroma<-matrix(nrow=20, ncol=6)

for(j in 1:length(nets)){
  coef_gsea<-nets[[j]][[1]]$kExt1/nets[[j]][[1]]$kInt1
  names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
  
  fgseaRes <- fgseaMultilevel(DE_stroma, coef_gsea, scoreType = "pos")
  gsea_stroma[,j]<-unlist(fgseaRes[,2])
  
}

names(DE_epi)<-c(1:35)
gsea_epi<-matrix(nrow=35, ncol=6)

for(j in 1:length(nets)){
  coef_gsea<-nets[[j]][[1]]$kExt2/nets[[j]][[1]]$kInt2
  names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
  
  fgseaRes <- fgseaMultilevel(DE_epi, coef_gsea, scoreType = "pos")
  gsea_epi[,j]<-unlist(fgseaRes[,2])
  
}

random_repeated_stroma<-c()
for(i in 1:100){
names(DE_stroma)<-c(1:20)
gsea_stroma_rand<-matrix(nrow=20, ncol=6)

for(j in 1:length(nets)){
  coef_gsea<-nets[[j]][[1]]$kExt1/nets[[j]][[1]]$kInt1
  names(coef_gsea)<-gsub("_1", "", names(coef_gsea))
  names(coef_gsea)<-sample(names(coef_gsea), length(names(coef_gsea)), replace = F)
  
  fgseaRes <- fgseaMultilevel(DE_stroma, coef_gsea, scoreType = "pos")
  gsea_stroma_rand[,j]<-unlist(fgseaRes[,2])
  
}

random_repeated_stroma<-c(random_repeated_stroma, sum(gsea_stroma_rand<0.05)) #1
}
#quantile(unlist(gsea_stroma_rand),0.05)
#5% 
#0.081
data <- data.frame(
  name=c("observed", "random"),
  value=c(sum(gsea_stroma<0.05), median(random_repeated_stroma)),
  sd=c(0, sd(random_repeated_stroma))
)

png("stroma_gsea_signif.png", res=300, 700, 1000)
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="darkgreen") +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, size=1.3) + theme_bw()
dev.off()


random_repeated_epi<-c()
for(i in 1:100){
names(DE_epi)<-c(1:35)
gsea_epi_rand<-matrix(nrow=35, ncol=6)

for(j in 1:length(nets)){
  coef_gsea<-nets[[j]][[1]]$kExt2/nets[[j]][[1]]$kInt2
  names(coef_gsea)<-gsub("_2", "", names(coef_gsea))
  names(coef_gsea)<-sample(names(coef_gsea), length(names(coef_gsea)), replace = F)
  
  fgseaRes <- fgseaMultilevel(DE_epi, coef_gsea, scoreType = "pos")
  gsea_epi_rand[,j]<-unlist(fgseaRes[,2])
  
}

random_repeated_epi<-c(random_repeated_epi, sum(gsea_epi_rand<0.05))
}#5

data <- data.frame(
  name=c("observed", "random"),
  value=c(sum(gsea_epi<0.05)/length(gsea_epi), median(random_repeated_epi)),
  sd=c(0, sd(random_repeated_epi))
)

png("epi_gsea_signif.png", res=300, 700, 1000)
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="brown2") +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, size=1.3) + theme_bw()
dev.off()
#quantile(unlist(gsea_epi_rand),0.05)
#5% 
#0.065