library(WGCNA)
library(openxlsx)
library(limma)
library(pheatmap)
library(Hmisc)
load("data/metadataAll.RData")
load("data/datasetsAll.RData")
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
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

#merging stroma and epi
rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")
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
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]
##############All genes double WGCNA
###use WGCNA with twice the genes
rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

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

###filtering data
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE14548<-rbind(stroma, epi)


#############GSE83591
load("data/GSE83591.RData")
GSE83591_meta<-read.csv("data/GSE83591_metadata.txt", sep="\t")
GSE83591_meta<-t(GSE83591_meta)
colnames(GSE83591_meta)<-GSE83591_meta[1,]
GSE83591_meta<-GSE83591_meta[-1,]

GSE83591_meta2<-read.xlsx("data/GSE83591_Correspondence_LCM.xlsx",1)
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
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE83591<-rbind(stroma, epi)

#######GSE68744
stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE68744" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE68744" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]
#sono matched, controllato

stroma<-datasetsAll[["GSE68744"]][,stromaID]
epi<-datasetsAll[["GSE68744"]][,epiID]

###filtering data
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

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
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE88715<-rbind(stroma, epi)


#######################RUN
source("scripts/crossWGCNA_functions_netdiff.R")
net_GSE5847<-crossWGCNA(data=data_merged_GSE5847, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
net_GSE10797<-crossWGCNA(data=data_merged_GSE10797, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
net_GSE14548<-crossWGCNA(data=data_merged_GSE14548, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
net_GSE83591<-crossWGCNA(data=data_merged_GSE83591, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
net_GSE68744<-crossWGCNA(data=data_merged_GSE68744, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
net_GSE88715<-crossWGCNA(data=data_merged_GSE88715, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")

source("scripts/crossWGCNA_functions_selfloops.R")
net_GSE5847_old<-crossWGCNA(data=data_merged_GSE5847, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
net_GSE10797_old<-crossWGCNA(data=data_merged_GSE10797, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
net_GSE14548_old<-crossWGCNA(data=data_merged_GSE14548, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
net_GSE83591_old<-crossWGCNA(data=data_merged_GSE83591, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
net_GSE68744_old<-crossWGCNA(data=data_merged_GSE68744, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
net_GSE88715_old<-crossWGCNA(data=data_merged_GSE88715, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")

save.image(file="net_compare_env.RData")


prop_same_genes<-function(net){
prop_same_all<-c()
for(i in 1:length(unique(net[[2]][[1]]))){
  clust<-unique(net[[2]][[1]])[i]
  clust_genes<-names(which(net[[2]][[1]]==clust))
  prop_same<-sum(table(gsub("_tis1|_tis2", "", clust_genes))==2)/length(clust_genes)
  prop_same_all<-c(prop_same_all, prop_same)
}
return(prop_same_all)
}


df_new<-data.frame(proportion=c(prop_same_genes(net_GSE5847),
prop_same_genes(net_GSE10797),
prop_same_genes(net_GSE14548),
prop_same_genes(net_GSE83591),
prop_same_genes(net_GSE68744),
prop_same_genes(net_GSE88715)),
dataset=c(rep("GSE5847", length(prop_same_genes(net_GSE5847))),
rep("GSE10797", length(prop_same_genes(net_GSE10797))),
rep("GSE14548", length(prop_same_genes(net_GSE14548))),
rep("GSE83591", length(prop_same_genes(net_GSE83591))),
rep("GSE68744", length(prop_same_genes(net_GSE68744))),
rep("GSE88715", length(prop_same_genes(net_GSE88715)))
),
algorithm=rep("indirect loops", length(c(prop_same_genes(net_GSE5847),
prop_same_genes(net_GSE10797),
prop_same_genes(net_GSE14548),
prop_same_genes(net_GSE83591),
prop_same_genes(net_GSE68744),
prop_same_genes(net_GSE88715)))))

df_old<-data.frame(proportion=c(prop_same_genes(net_GSE5847_old),
prop_same_genes(net_GSE10797_old),
prop_same_genes(net_GSE14548_old),
prop_same_genes(net_GSE83591_old),
prop_same_genes(net_GSE68744_old),
prop_same_genes(net_GSE88715_old)),
dataset=c(rep("GSE5847", length(prop_same_genes(net_GSE5847_old))),
rep("GSE10797", length(prop_same_genes(net_GSE10797_old))),
rep("GSE14548", length(prop_same_genes(net_GSE14548_old))),
rep("GSE83591", length(prop_same_genes(net_GSE83591_old))),
rep("GSE68744", length(prop_same_genes(net_GSE68744_old))),
rep("GSE88715", length(prop_same_genes(net_GSE88715_old)))
),
algorithm=rep("self loops", length(c(prop_same_genes(net_GSE5847_old),
prop_same_genes(net_GSE10797_old),
prop_same_genes(net_GSE14548_old),
prop_same_genes(net_GSE83591_old),
prop_same_genes(net_GSE68744_old),
prop_same_genes(net_GSE88715_old)))))

df_tot<-rbind.data.frame(df_old, df_new)

pdf("results/algorithm_compare.pdf", 10, 5)
ggplot(df_tot, aes(x=dataset, y=proportion, fill=algorithm))+geom_boxplot()+theme_classic()
dev.off()


#######
##Confronto adjacencies
######
source("scripts/crossWGCNA_functions_netdiff.R")
adj_GSE5847<-Adjacency(data=data_merged_GSE5847, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")
source("scripts/crossWGCNA_functions_selfloops.R")
adj_GSE5847_old<-Adjacency(data=data_merged_GSE5847, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2")

##la correlazione di un gene con glo altri all'interno e all'esterno del tessuto è simile perché in entrambi i tessuti è coinvolto negli stessi pathway
i<-10
corr_tis1<-adj_GSE5847_old[i,grep("tis1", colnames(adj_GSE5847_old))][-i]
corr_tis2<-adj_GSE5847_old[i,grep("tis2", colnames(adj_GSE5847_old))][-i]
cor(corr_tis1, corr_tis2)
pdf("Example_cor_selfloops.pdf",5,5)
plot(corr_tis1, corr_tis2, pch=19)
dev.off()

cor_tot<-c()
for(i in 1:(nrow(adj_GSE5847_old)/2)){
  corr_tis1<-adj_GSE5847_old[i,grep("tis1", colnames(adj_GSE5847_old))][-i]
corr_tis2<-adj_GSE5847_old[i,grep("tis2", colnames(adj_GSE5847_old))][-i]
cor_tot<-c(cor_tot, cor(corr_tis1, corr_tis2))
}

pdf("All_cor_selfloops.pdf",5,5)
hist(cor_tot)
dev.off()



i<-10
corr_tis1<-adj_GSE5847[i,grep("tis1", colnames(adj_GSE5847))][-i]
corr_tis2<-adj_GSE5847[i,grep("tis2", colnames(adj_GSE5847))][-i]
cor(corr_tis1, corr_tis2)
pdf("Example_cor_netdiff.pdf",5,5)
plot(corr_tis1, corr_tis2, pch=19)
dev.off()

cor_tot_alt<-c()
for(i in 1:(nrow(adj_GSE5847)/2)){
  corr_tis1<-adj_GSE5847[i,grep("tis1", colnames(adj_GSE5847))][-i]
corr_tis2<-adj_GSE5847[i,grep("tis2", colnames(adj_GSE5847))][-i]
cor_tot_alt<-c(cor_tot_alt, cor(corr_tis1, corr_tis2))
}

pdf("All_cor_netdiff.pdf",5,5)
hist(cor_tot_alt)
dev.off()
