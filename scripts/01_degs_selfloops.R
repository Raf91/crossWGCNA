library(WGCNA)
library(openxlsx)
library(limma)
library(pheatmap)
library(Hmisc)
load("data/metadataAll.RData")
load("data/datasetsAll.RData")
source("scripts/01_degs_selfloops.R")

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

degrees_GSE5847<-network(data=data_merged_GSE5847, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2", doTOM=T, ds=4, crossOnly=T)
degrees_GSE10797<-network(data=data_merged_GSE10797, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2", doTOM=T, ds=4, crossOnly=T)
degrees_GSE14548<-network(data=data_merged_GSE14548, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2", doTOM=T, ds=4, crossOnly=T)
degrees_GSE83591<-network(data=data_merged_GSE83591, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2", doTOM=T, ds=4, crossOnly=T)
degrees_GSE68744<-network(data=data_merged_GSE68744, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2", doTOM=T, ds=4, crossOnly=T)
degrees_GSE88715<-network(data=data_merged_GSE88715, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_tis1", comp2="_tis2", doTOM=T, ds=4, crossOnly=T)

degs<-list(degrees_GSE5847, degrees_GSE10797, degrees_GSE14548,
           degrees_GSE83591, degrees_GSE68744, degrees_GSE88715)

save(degs, file="results/degs_3rd_selfloops.RData")
