load("data/CAFs_cocultures.RData")
load("results/degs_3rd_netdiff.RData")
library(clusterProfiler)
library(openxlsx)
library(ggplot2)

names(degs)<-c("GSE5847", "GSE10797", "GSE14548", "GSE83591",
               "GSE68744", "GSE88715")


#i geni DE sono tra i top interacting? (fare un gsea)

###GSEA
library(msigdbr)
library(fgsea)
library(data.table)

names(DE_stroma)<-c(1:20)
gsea_stroma<-matrix(nrow=20, ncol=6)

for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt1/degs[[j]]$kInt1
  names(coef_gsea)<-gsub("_tis1", "", names(coef_gsea))

  fgseaRes <- fgseaMultilevel(DE_stroma, coef_gsea)
  gsea_stroma[,j]<-unlist(fgseaRes[,2])

}

names(DE_epi)<-c(1:35)
gsea_epi<-matrix(nrow=35, ncol=6)

for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt2/degs[[j]]$kInt2
  names(coef_gsea)<-gsub("_tis2", "", names(coef_gsea))

  fgseaRes <- fgseaMultilevel(DE_epi, coef_gsea)
  gsea_epi[,j]<-unlist(fgseaRes[,2])

}

set.seed(46956305)
random_repeated_stroma<-c()
for(i in 1:100){
names(DE_stroma)<-c(1:20)
gsea_stroma_rand<-matrix(nrow=20, ncol=6)

for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt1/degs[[j]]$kInt1
  names(coef_gsea)<-gsub("_tis1", "", names(coef_gsea))
  names(coef_gsea)<-sample(names(coef_gsea), length(names(coef_gsea)), replace = F)

  fgseaRes <- fgseaMultilevel(DE_stroma, coef_gsea)
  gsea_stroma_rand[,j]<-unlist(fgseaRes[,2])

}

random_repeated_stroma<-c(random_repeated_stroma, sum(gsea_stroma_rand<0.05, na.rm = T)) #1
}
#quantile(unlist(gsea_stroma_rand),0.05)
#5%
#0.081
data <- data.frame(
  name=c("observed", "random"),
  value=c(sum(gsea_stroma<0.05, na.rm=T), median(random_repeated_stroma, na.rm=T)),
  sd=c(0, sd(random_repeated_stroma))
)

pdf("results/stroma_gsea_signif_netdiff.pdf", 4,5)
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="darkgreen") +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, size=1.3) + theme_bw()
dev.off()

save(random_repeated_stroma, file="results/random_repeated_stroma_netdiff.RData")
load(file="results/random_repeated_stroma_netdiff.RData")

set.seed(46956305)
random_repeated_epi<-c()
for(i in 1:100){
names(DE_epi)<-c(1:35)
gsea_epi_rand<-matrix(nrow=35, ncol=6)

for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt2/degs[[j]]$kInt2
  names(coef_gsea)<-gsub("_tis2", "", names(coef_gsea))
  names(coef_gsea)<-sample(names(coef_gsea), length(names(coef_gsea)), replace = F)

  fgseaRes <- fgseaMultilevel(DE_epi, coef_gsea)
  gsea_epi_rand[,j]<-unlist(fgseaRes[,2])

}

random_repeated_epi<-c(random_repeated_epi, sum(gsea_epi_rand<0.05, na.rm = T))
}#5

save(random_repeated_epi, file="results/random_repeated_epi_netdiff.RData")
load(file="results/random_repeated_epi_netdiff.RData")

data <- data.frame(
  name=c("observed", "random"),
  value=c(sum(gsea_epi<0.05, na.rm=T), median(random_repeated_epi, na.rm=T)),
  sd=c(0, sd(random_repeated_epi))
)

pdf("results/epi_gsea_signif_netdiff.pdf", 4,5)
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="brown2") +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, size=1.3) + theme_bw()
dev.off()
#quantile(unlist(gsea_epi_rand),0.05)
#5%
#0.065
