load("data/cocultures_CAFsonly.RData")
load("results/nets_3rdSep.RData")
library(clusterProfiler)
library(openxlsx)

names(nets)<-c("GSE5847", "GSE10797", "GSE14548", "GSE83591",
               "GSE68744", "GSE88715")


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

set.seed(46956305)
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

png("stroma_gsea_signif_sep.png", res=300, 700, 1000)
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="darkgreen") +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, size=1.3) + theme_bw()
dev.off()

set.seed(46956305)
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
  value=c(sum(gsea_epi<0.05), median(random_repeated_epi)),
  sd=c(0, sd(random_repeated_epi))
)

png("epi_gsea_signif_sep.png", res=300, 700, 1000)
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="brown2") +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, size=1.3) + theme_bw()
dev.off()
#quantile(unlist(gsea_epi_rand),0.05)
#5% 
#0.065