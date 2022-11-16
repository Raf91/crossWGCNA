library(WGCNA)
library(xlsx)
library(limma)
library(pheatmap)
library(Hmisc)
load("D:/Dropbox (MBC)/Aurora/R analyses/Metanalyses/Breast tumor stroma/Shared/Aurora/15Apr2021/metadataAll.RData")
load("D:/Dropbox (MBC)/Aurora/R analyses/Metanalyses/Breast tumor stroma/Shared/Aurora/15Apr2021/datasetsAll.RData")
changenames<-function(data, anno){
  annotation_sel=anno[match( rownames(data), anno[,1]),2]
  
  if(length(which(annotation_sel==""))>0){
    data<-data[-which(annotation_sel==""),]
    annotation_sel<-annotation_sel[-which(annotation_sel=="")]
  }
  
  a<-which(duplicated(annotation_sel))
  while(length(a)>0){
    for(i in 1:length(unique(annotation_sel))){
      if(length(which(annotation_sel==unique(annotation_sel)[i]))>1){
        m=which.max(rowMeans(data[which(annotation_sel==unique(annotation_sel)[i]),], na.rm=T))
        data=data[-which(annotation_sel==unique(annotation_sel)[i])[-m],]
        annotation_sel=annotation_sel[-which(annotation_sel==unique(annotation_sel)[i])[-m]]
      }
    }
    
    data=data[which(is.na(annotation_sel)==F),]
    annotation_sel=na.omit(annotation_sel)
    a<-which(duplicated(annotation_sel))
  }
  
  rownames(data)=annotation_sel
  return(data)
}

#Datasets
###dataset GSE5847

stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

stromaMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

inboth<-intersect(stromaMatch, epiMatch)

stromaID<-stromaID[match(inboth, stromaMatch)]
epiID<-epiID[match(inboth, epiMatch)]

stroma<-datasetsAll[["GSE5847"]][,stromaID]
epi<-datasetsAll[["GSE5847"]][,epiID]

###filtering data
stroma<-stroma[apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25),]
epi<-epi[apply(epi,1,var)>=quantile(apply(epi,1,var),0.25),]

#merging stroma and epi
rownames(stroma)<-paste(rownames(stroma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(stroma)

data_merged_GSE5847<-rbind(stroma, epi)


############ GSE10797
stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE10797" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE10797" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

stromaMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE10797" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE10797" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

inboth<-intersect(stromaMatch, epiMatch)

stromaID<-stromaID[match(inboth, stromaMatch)]
epiID<-epiID[match(inboth, epiMatch)]

stroma<-datasetsAll[["GSE10797"]][,stromaID]
epi<-datasetsAll[["GSE10797"]][,epiID]

###filtering data
stroma<-stroma[apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25),]
epi<-epi[apply(epi,1,var)>=quantile(apply(epi,1,var),0.25),]
##############All genes double WGCNA
###use WGCNA with twice the genes
rownames(stroma)<-paste(rownames(stroma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE10797<-rbind(stroma, epi)

##############GSE14548
stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE14548" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE14548" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

stromaMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE14548" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE14548" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

inboth<-intersect(stromaMatch, epiMatch)

stromaID<-stromaID[match(inboth, stromaMatch)]
epiID<-epiID[match(inboth, epiMatch)]

stroma<-datasetsAll[["GSE14548"]][,stromaID]
epi<-datasetsAll[["GSE14548"]][,epiID]

rownames(stroma)<-paste(rownames(stroma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")

###filtering data
stroma<-stroma[apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25),]
epi<-epi[apply(epi,1,var)>=quantile(apply(epi,1,var),0.25),]

colnames(epi)<-colnames(stroma)
data_merged_GSE14548<-rbind(stroma, epi)


#############GSE83591
load("D:/Dropbox (MBC)/Aurora/R analyses/Data/Stroma datasets/Processed/GSE83591.RData")
GSE83591_meta<-read.csv("D:/Dropbox (MBC)/Aurora/R analyses/Data/Stroma datasets/Metadata/GSE83591_metadata.txt", sep="\t")
GSE83591_meta<-t(GSE83591_meta)
colnames(GSE83591_meta)<-GSE83591_meta[1,]
GSE83591_meta<-GSE83591_meta[-1,]

GSE83591_meta2<-read.xlsx("D:/Dropbox (MBC)/Aurora/R analyses/Data/Stroma datasets/Metadata/GSE83591_Correspondence_LCM.xlsx",1)
GSE83591_meta2$ID<-unlist(strsplit(GSE83591_meta2$GSE83591, "_"))[seq(1,109*9,9)]

#reorder based on annotation file
GSE83591<-GSE83591[,match(GSE83591_meta2$ID, colnames(GSE83591))]
GSE83591_meta<-GSE83591_meta[match(GSE83591_meta2$ID, GSE83591_meta[,1]),]

GSE83591<-GSE83591[,which(GSE83591_meta2$Cy3 %in%c("TE", "TS"))]
GSE83591_meta<-GSE83591_meta[which(GSE83591_meta2$Cy3 %in%c("TE", "TS")),]
GSE83591_meta2<-GSE83591_meta2[which(GSE83591_meta2$Cy3 %in%c("TE", "TS")),]

GSE83591<-GSE83591[,-c(15,34,67)]
GSE83591_meta<-GSE83591_meta[-c(15,34,67),]
GSE83591_meta2<-GSE83591_meta2[-c(15,34,67),]

stroma<-GSE83591[,which(GSE83591_meta2$Cy3=="TS")]
epi<-GSE83591[,which(GSE83591_meta2$Cy3=="TE")]

###filtering data
stroma<-stroma[apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25),]
epi<-epi[apply(epi,1,var)>=quantile(apply(epi,1,var),0.25),]

rownames(stroma)<-paste(rownames(stroma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE83591<-rbind(stroma, epi)

#######GSE68744
stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE68744" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE68744" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]
#sono matched, controllato

stroma<-datasetsAll[["GSE68744"]][,stromaID]
epi<-datasetsAll[["GSE68744"]][,epiID]

###filtering data
stroma<-stroma[apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25),]
epi<-epi[apply(epi,1,var)>=quantile(apply(epi,1,var),0.25),]

rownames(stroma)<-paste(rownames(stroma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE68744<-rbind(stroma, epi)

#######GSE88715
stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE88715" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE88715" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

stromaMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE88715" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE88715" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

inboth<-intersect(stromaMatch, epiMatch)

stromaID<-stromaID[match(inboth, stromaMatch)]
epiID<-epiID[match(inboth, epiMatch)]

stroma<-datasetsAll[["GSE88715"]][,stromaID]
epi<-datasetsAll[["GSE88715"]][,epiID]

###filtering data
stroma<-stroma[apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25),]
epi<-epi[apply(epi,1,var)>=quantile(apply(epi,1,var),0.25),]

rownames(stroma)<-paste(rownames(stroma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE88715<-rbind(stroma, epi)


#####################
########################
############# FUNCTIONS
#######################
##########################

#cortype: pearson, spearman, bicor
#Adj_type: signed, unsigned
#pval :threshold, weight, none
#thr: default 0.05
#beta: default 6 (per avere adjacency=none mettere beta=1)
###kInt/kExt
Adjacency<-function(data, Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6 ){
  if(pval=="none"){
    if(cortype=="bicor"){
      A<-bicor(t(data))
    } else {
      A<-cor(t(data), method=cortype)
    }
  } else {
    mat<-rcorr(t(data), type=cortype)
    A<-mat[[1]]
    p.val<-mat[[3]]
    rm(mat)
    if(pval=="threhsold"){
      A[which(p.val>thr)]<-NA
    } else if (pval=="weight"){
      A<-A*(1-p.val)
    }
  }
  
  if(Adj_type=="signed"){
    A<-(0.5 * (1+A) )^beta
  } else if(Adj_type=="unsigned"){
    A<-(abs(A))^beta  
  } 
  rownames(A)<-rownames(A)
  colnames(A)<-rownames(A)
  return(A)
}

degrees<-function(A, comp1="_1", comp2="_2"){
  #remove correlations with the same gene
  
  genes1<-gsub(comp1, "", rownames(A)[grep(comp1, rownames(A))])
  genes2<-gsub(comp2, "", rownames(A)[grep(comp2, rownames(A))])
  Idx<-cbind(grep(comp1, rownames(A)), grep(comp2, rownames(A))[match(genes1, genes2)])
  A[Idx]<-0
  diag(A)<-0
  
  kTot<-rowSums(A)
  kInt_1<-rowSums(A[grep(comp1, rownames(A)),grep(comp1, colnames(A))])
  kInt_2<-rowSums(A[grep(comp2, rownames(A)),grep(comp2, colnames(A))])
  kExt_1<-rowSums(A[grep(comp1, rownames(A)),grep(comp2, colnames(A))])
  kExt_2<-rowSums(A[grep(comp2, rownames(A)),grep(comp1, colnames(A))])
  k<-list(kInt1=kInt_1, kInt2=kInt_2, kExt1=kExt_1, kExt2=kExt_2, kTot1=kTot[grep(comp1, colnames(A))], kTot2=kTot[grep(comp2, colnames(A))])
  return(k)
}


clusteringWGCNA<-function(A, data, comp1="_1", comp2="_2", TOM=T, ds=1, crossOnly=T ){
  if(TOM==T){
    similarity<-TOMsimilarity(A, TOMType="signed")
  } else {
    similarity<-A
  }
  #remove correlations with the same gene
  genes1<-gsub(comp1, "", rownames(A)[grep(comp1, rownames(A))])
  genes2<-gsub(comp2, "", rownames(A)[grep(comp2, rownames(A))])
  Idx<-cbind(grep(comp1, rownames(A)), grep(comp2, rownames(A))[match(genes1, genes2)])
  similarity[Idx]<-0
  
  similarity[grep(comp1, rownames(data)), grep(comp2, rownames(data))[match(genes1, genes2)]]<-0
  
  if(crossOnly==T){
    similarity[grep(comp1, rownames(data)), grep(comp1, rownames(data))]<-0
    similarity[grep(comp2, rownames(data)), grep(comp2, rownames(data))]<-0
  }
  
  conTree=hclust(as.dist(1-similarity), method="average")
  unmergedLabels=cutreeDynamic(dendro=conTree, dist=1-similarity, deepSplit=ds)
  names(unmergedLabels)<-rownames(data)
  merged=mergeCloseModules(t(data), unmergedLabels, cutHeight=0.25, verbose=3)
  names(merged$colors)<-rownames(data)
  
  clusters<-list(unmerged=unmergedLabels, merged=merged)
  return(clusters)
}


network<-function(data, Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=1, crossOnly=T){
  Adj<-Adjacency(data=data, Adj_type=Adj_type, cortype=cortype,pval=pval, thr=thr, beta=beta )
  k<-degrees(A=Adj, comp1=comp1, comp2=comp2)
  clusters<-clusteringWGCNA(A=Adj, data=data, comp1=comp1, comp2=comp2, TOM=doTOM, ds=ds, crossOnly=crossOnly )
  net_out<-list(k,clusters)
  return(net_out)
}


####enrichment for coculture degs
#modules
DEGs_GSE41678<-function(data, ind1, ind2, ind3, ind4){
  data_BC<-data[,c(ind1, ind2)]
  design <- model.matrix(~ 0+factor(c(rep("Mono",length(ind1)), rep("Co",length(ind2)))))
  colnames(design) <- c("Co", "Mono")
  
  cont_matrix <- makeContrasts(CovsMono = Co-Mono, levels=design)
  
  fit <- lmFit(data_BC, design)
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  fit_contrast <- eBayes(fit_contrast)
  top_genes <- topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")
  result <- decideTests(fit_contrast)
  
  downstroma<-rownames(top_genes)[top_genes$adj.P.Val<0.05&top_genes$logFC<(-1)]
  upstroma<-rownames(top_genes)[top_genes$adj.P.Val<0.05&top_genes$logFC>1]
  DEGs_stroma<-c( upstroma)
  
  data_BC<-data[,c(ind3, ind4)]
  design <- model.matrix(~ 0+factor(c(rep("Mono", length(ind3)), rep("Co",length(ind4)))))
  colnames(design) <- c("Co", "Mono")
  
  cont_matrix <- makeContrasts(CovsMono = Co-Mono, levels=design)
  
  fit <- lmFit(data_BC, design)
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  fit_contrast <- eBayes(fit_contrast)
  top_genes <- topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")
  result <- decideTests(fit_contrast)
  
  downepi<-rownames(top_genes)[top_genes$adj.P.Val<0.05&top_genes$logFC<(-1)]
  upepi<-rownames(top_genes)[top_genes$adj.P.Val<0.05&top_genes$logFC>1]
  DEGs_epi<-c(upepi)
  DEGs<-list(DEGs_stroma, DEGs_epi)
  return(DEGs)
}

enrich_coculture<-function(data, modules, DEGs_stroma, DEGs_epi, comp1="_1", comp2="_2"){
  
  ngenes_stroma<-length(grep(comp1, rownames(data)))
  ngenes_epi<-length(grep(comp2, rownames(data)))
  modules_stroma<-modules[1:ngenes_stroma]
  names(modules_stroma)<-gsub(comp1, "",names(modules_stroma))
  modules_epi<-modules[ngenes_stroma+1:length(modules)]
  names(modules_epi)<-gsub(comp2, "",names(modules_epi))
  
pval_DEGsstroma<-c()
for(i in 1:max(modules)){
  a<-length(intersect(names(modules_stroma)[modules_stroma==i],DEGs_stroma))
  b<-length(intersect(names(modules_stroma)[modules_stroma!=i&modules_stroma!=0],DEGs_stroma))
  c<-length(setdiff(names(modules_stroma)[modules_stroma==i], DEGs_stroma))
  d<-length(setdiff(names(modules_stroma)[modules_stroma!=i&modules_stroma!=0],DEGs_stroma))
  pval_DEGsstroma<-c(pval_DEGsstroma, fisher.test(matrix(c(a,b,c,d), byrow = T, nrow=2), alternative = "greater")$p.value)
}


pval_DEGsepi<-c()
for(i in 1:max(modules)){
  a<-length(intersect(names(modules_epi)[modules_epi==i],DEGs_epi))
  b<-length(intersect(names(modules_epi)[modules_epi!=i&modules_epi!=0],DEGs_epi))
  c<-length(setdiff(names(modules_epi)[modules_epi==i], DEGs_epi))
  d<-length(setdiff(names(modules_epi)[modules_epi!=i&modules_epi!=0],DEGs_epi))
  pval_DEGsepi<-c(pval_DEGsstroma, fisher.test(matrix(c(a,b,c,d), byrow = T, nrow=2), alternative = "greater")$p.value)
}

pvals<-list(pval_DEGsstroma, pval_DEGsepi)
return(pvals)
}

whole_script<-function(todir, Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=1, crossOnly=T){
setwd(dir)
dir.create(todir)
setwd(todir)

argg <- c(as.list(environment()))
write.csv(argg, "parameters.csv")
####Calculate networks for each dataset

net_GSE5847<-network(data=data_merged_GSE5847, Adj_type=Adj_type, cortype=cortype, pval=pval, thr=thr, beta=beta, comp1=comp1, comp2=comp2, doTOM=doTOM, ds=ds, crossOnly=crossOnly)
net_GSE10797<-network(data=data_merged_GSE10797, Adj_type=Adj_type, cortype=cortype, pval=pval, thr=thr, beta=beta, comp1=comp1, comp2=comp2, doTOM=doTOM, ds=ds, crossOnly=crossOnly)
net_GSE14548<-network(data=data_merged_GSE14548, Adj_type=Adj_type, cortype=cortype, pval=pval, thr=thr, beta=beta, comp1=comp1, comp2=comp2, doTOM=doTOM, ds=ds, crossOnly=crossOnly)
net_GSE83591<-network(data=data_merged_GSE83591, Adj_type=Adj_type, cortype=cortype, pval=pval, thr=thr, beta=beta, comp1=comp1, comp2=comp2, doTOM=doTOM, ds=ds, crossOnly=crossOnly)
net_GSE68744<-network(data=data_merged_GSE68744, Adj_type=Adj_type, cortype=cortype, pval=pval, thr=thr, beta=beta, comp1=comp1, comp2=comp2, doTOM=doTOM, ds=ds, crossOnly=crossOnly)
net_GSE88715<-network(data=data_merged_GSE88715, Adj_type=Adj_type, cortype=cortype, pval=pval, thr=thr, beta=beta, comp1=comp1, comp2=comp2, doTOM=doTOM, ds=ds, crossOnly=crossOnly)


nets<-list( GSE5847=net_GSE5847, GSE10797=net_GSE10797, GSE14548=net_GSE14548, GSE83591=net_GSE83591, GSE68744=net_GSE68744, GSE88715=net_GSE88715)
save(nets, file="nets.RData")
}


net_GSE5847<-network(data=data_merged_GSE5847, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
net_GSE10797<-network(data=data_merged_GSE10797, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
net_GSE14548<-network(data=data_merged_GSE14548, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
net_GSE83591<-network(data=data_merged_GSE83591, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
net_GSE68744<-network(data=data_merged_GSE68744, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)
net_GSE88715<-network(data=data_merged_GSE88715, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=4, crossOnly=T)

nets<-list(net_GSE5847, net_GSE10797, net_GSE14548, 
           net_GSE83591, net_GSE68744, net_GSE88715)

save(nets, file="nets_3rd.RData")