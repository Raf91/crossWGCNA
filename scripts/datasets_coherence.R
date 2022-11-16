load("results/degs_3rd_alt.RData")
source("scripts/crossWGCNA_alt2.R")
degs<-list(net_GSE5847, net_GSE10797, net_GSE14548, net_GSE83591, net_GSE68744, net_GSE88715)
library(clusterProfiler)
library(openxlsx)
names(degs)<-c("GSE5847", "GSE10797", "GSE14548", "GSE83591",
               "GSE68744", "GSE88715")

corrs<-matrix(nrow=length(degs), ncol=length(degs))
for(i in 1:length(degs)){
    for(j in 1:length(degs)){
        incommon<-intersect(names(degs[[i]][[1]])[grep("tis1", names(degs[[i]][[1]]))], 
        names(degs[[j]][[1]])[grep("tis1", names(degs[[j]][[1]]))])
corrs[i,j]<-cor(degs[[i]][[1]][incommon], degs[[j]][[1]][incommon])
    }
}

degs_old<-list(net_GSE5847_old, net_GSE10797_old, net_GSE14548_old, net_GSE83591_old, net_GSE68744_old, net_GSE88715_old)
save(degs_old, file="degs_old.RData")
load("/Volumes/iorio/Aurora/crossWGCNA/degs_old.RData")

corrs_old<-matrix(nrow=length(degs_old), ncol=length(degs_old))
for(i in 1:length(degs_old)){
    for(j in 1:length(degs_old)){
        incommon<-intersect(names(degs_old[[i]][[1]][[1]])[grep("tis1", names(degs_old[[i]][[1]][[1]]))], 
        names(degs_old[[j]][[1]][[1]])[grep("tis1", names(degs_old[[j]][[1]][[1]]))])
corrs_old[i,j]<-cor(degs_old[[i]][[1]][[1]][incommon], degs_old[[j]][[1]][[1]][incommon])
    }
}
save(degs, file="degs.RData")

hist(corrs_old, breaks=20)
hist(corrs, breaks=20)
