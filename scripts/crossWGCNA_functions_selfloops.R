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
Adjacency<-function(data, comp1="_1",comp2="_2", Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6, sign_list=1, which.sign="none" ){
  comp1<-paste(comp1, "$", sep="")
  comp2<-paste(comp2, "$", sep="")

  if(pval=="none"){
    if(cortype=="bicor"){
      A<-bicor(t(data))
    } else {
      A<-cor(t(data), method=cortype)
    }
  } else {
    if(cortype=="bicor"){
      paste("Can't use pval with bicor")
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
}
  ####se si ha un network con segno in input
  if(which.sign=="comp2"){
    A[grep(comp2, rownames(A)), grep(comp1, rownames(A))]<-A[grep(comp2, rownames(A)), grep(comp1, rownames(A))]*sign_list
    A[grep(comp1, rownames(A)), grep(comp2, rownames(A))]<-t(A[grep(comp2, rownames(A)), grep(comp1, rownames(A))]*sign_list)
  } else if(which.sign=="comp1"){
    A[grep(comp1, rownames(A)), grep(comp2, rownames(A))]<-A[grep(comp1, rownames(A)), grep(comp2, rownames(A))]*sign_list
    A[grep(comp2, rownames(A)), grep(comp1, rownames(A))]<-t(A[grep(comp1, rownames(A)), grep(comp2, rownames(A))]*sign_list)
  }

  ####


  if(Adj_type=="signed"){
    A<-(0.5 * (1+A) )^beta
  } else if(Adj_type=="unsigned"){
    A<-(abs(A))^beta
  } else if(Adj_type=="keep sign"){
    A<-((abs(A))^beta )*sign(A)
  }
  return(A)
}

degrees<-function(A, comp1="_1", comp2="_2"){
  #remove correlations with the same gene
  comp1<-paste(comp1, "$", sep="")
  comp2<-paste(comp2, "$", sep="")

  genes1<-gsub(comp1, "", rownames(A)[grep(comp1, rownames(A))])
  genes2<-gsub(comp2, "", rownames(A)[grep(comp2, rownames(A))])
  Idx1<-cbind(grep(comp1, rownames(A)), grep(comp2, rownames(A))[match(genes1, genes2)])
  A[Idx1]<-0
  Idx2<-cbind(grep(comp2, rownames(A)), grep(comp1, rownames(A))[match(genes1, genes2)])
  A[Idx2]<-0

  kTot<-rowSums(A)
  kInt_1<-rowSums(A[grep(comp1, rownames(A)),grep(comp1, colnames(A))])
  kInt_2<-rowSums(A[grep(comp2, rownames(A)),grep(comp2, colnames(A))])
  kExt_1<-rowSums(A[grep(comp1, rownames(A)),grep(comp2, colnames(A))])
  kExt_2<-rowSums(A[grep(comp2, rownames(A)),grep(comp1, colnames(A))])
  k<-list(kInt1=kInt_1, kInt2=kInt_2, kExt1=kExt_1, kExt2=kExt_2, kTot1=kTot[grep(comp1, colnames(A))], kTot2=kTot[grep(comp2, colnames(A))])
  return(k)
}


clusteringWGCNA<-function(A, data, comp1="_1", comp2="_2", TOM=T, ds=1, crossOnly=T){
  comp1<-paste(comp1, "$", sep="")
  comp2<-paste(comp2, "$", sep="")


  #remove correlations with the same gene

  if(crossOnly==T){
    A[grep(comp1, rownames(A)), grep(comp1, rownames(A))]<-0
    A[grep(comp2, rownames(A)), grep(comp2, rownames(A))]<-0
  }

  genes1<-gsub(comp1, "", rownames(A)[grep(comp1, rownames(A))])
  genes2<-gsub(comp2, "", rownames(A)[grep(comp2, rownames(A))])
  Idx1<-cbind(grep(comp1, rownames(A)), grep(comp2, rownames(A))[match(genes1, genes2)])
  A[Idx1]<-0
  Idx2<-cbind(grep(comp2, rownames(A)), grep(comp1, rownames(A))[match(genes1, genes2)])
  A[Idx2]<-0

  if(TOM==T){
    similarity<-TOMsimilarity(A, TOMType="signed")
  } else {
    similarity<-A
  }
  rownames(similarity)<-rownames(A)
  colnames(similarity)<-colnames(A)


  if(crossOnly==T){
    similarity[grep(comp1, rownames(similarity)), grep(comp1, rownames(similarity))]<-0
    similarity[grep(comp2, rownames(similarity)), grep(comp2, rownames(similarity))]<-0
  }
  #similarity[grep(comp1, rownames(data)), grep(comp2, rownames(data))[match(genes1, genes2)]]<-0

  conTree=hclust(as.dist(1-similarity), method="average")
  unmergedLabels=cutreeDynamic(dendro=conTree, dist=1-similarity, deepSplit=ds)
  names(unmergedLabels)<-rownames(similarity)
  merged=mergeCloseModules(t(data), unmergedLabels, cutHeight=0.25, verbose=3)
  names(merged$colors)<-rownames(similarity)

  return(merged)
}


crossWGCNA<-function(data, Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=1, crossOnly=T, sign_list=1, which.sign="none"){
  comp1<-paste(comp1, "$", sep="")
  comp2<-paste(comp2, "$", sep="")
  Adj<-Adjacency(data=data, Adj_type=Adj_type, cortype=cortype,pval=pval, thr=thr, beta=beta, sign_list=sign_list, which.sign=which.sign )
  k<-degrees(A=Adj, comp1=comp1, comp2=comp2)
  clusters<-clusteringWGCNA(A=Adj, data=data, comp1=comp1, comp2=comp2, TOM=doTOM, ds=ds, crossOnly=crossOnly )
  net_out<-list(k,clusters)
  return(net_out)
}

network<-function(data, Adj_type="signed", cortype="spearman", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", doTOM=T, ds=1, crossOnly=T, sign_list=1, which.sign="none"){
  comp1<-paste(comp1, "$", sep="")
  comp2<-paste(comp2, "$", sep="")
  Adj<-Adjacency(data=data, Adj_type=Adj_type, cortype=cortype,pval=pval, thr=thr, beta=beta, sign_list=sign_list, which.sign=which.sign )
  k<-degrees(A=Adj, comp1=comp1, comp2=comp2)
  return(k)
}
