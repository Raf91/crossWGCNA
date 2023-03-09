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
  Idx2<-cbind(grep(comp2, rownames(A)), grep(comp1, rownames(A))[match(genes2, genes1)])
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
  Idx2<-cbind(grep(comp2, rownames(A)), grep(comp1, rownames(A))[match(genes2, genes1)])
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


degrees_mod<-function(data,
                      modules,
                      Adj_type = "signed",
                      cortype = "spearman",
                      pval = "none",
                      thr = 0.05,
                      beta = 6,
                      comp1 = "_1",
                      comp2 = "_2") {


  k<-list()
  for(i in 1:length(unique(modules))){
    mod<-names(modules)[which(modules==unique(modules)[i])]
    genes<-unique(gsub(comp2, "", gsub(comp1, "", mod)))
    Adj <-
      Adjacency(
        data = data[c(paste(genes, comp1, sep=""), paste(genes, comp2, sep="")),],
        Adj_type = Adj_type,
        cortype = cortype,
        pval = pval,
        thr = thr,
        beta = beta,
        comp1 = comp1,
        comp2 = comp2
      )
    Adj<-Adj[mod, mod]
    k[[i]] <- degrees(A = Adj, comp1 = comp1, comp2 = comp2)

  }
  return(k)
}


cytoscape_net<-function(adjacency=adj_GSE10797, dataset=data_merged_GSE10797, gene, comp1, comp2, num, corr="spearman"){
  interactors<-c(names(sort(adjacency[paste(gene, comp1, sep="_"),grep(comp2, colnames(adjacency))], decreasing=T))[1:num])

  inter_edges<-cor(dataset[paste(gene, comp1, sep="_"),], t(dataset[intersect(interactors, rownames(dataset)),]), method=corr)
  intra_1<-cor(dataset[paste(gene, comp1, sep="_"),], t(dataset[gsub(comp2, comp1, intersect(interactors, rownames(dataset))),]), method=corr)
  intra_2<-cor(dataset[paste(gene, comp2, sep="_"),], t(dataset[intersect(interactors, rownames(dataset)),]), method=corr)

  df<-data.frame(Source=c(rep(gene, length(inter_edges)),
                          rep(gene, length(intra_1)),
                          rep(gene, length(intra_2))),
                 Target=c(gsub(paste("_", comp2, sep=""), "", colnames(inter_edges)),
                          gsub(paste("_", comp1, sep=""), "", colnames(intra_1)),
                          gsub(paste("_", comp2, sep=""), "", colnames(intra_2))
                 ),
                 Weight=c(inter_edges,
                          intra_1,
                          intra_2
                 ),
                 Edge_type=c(rep("inter", length(inter_edges)),
                             rep("intra1", length(intra_1)),
                             rep("intra2", length(intra_2))),
                 Source_type=c(rep(comp1, length(inter_edges)),
                               rep(comp1, length(intra_1)),
                               rep(comp2,  length(intra_2))),
                 Target_type=c(rep(comp2, length(inter_edges)),
                               rep(comp1, length(intra_1)),
                               rep(comp2,  length(intra_2)))
  )

  write.csv(df, file=paste("Cytoscape", gene, comp1, "egdes.csv", sep="_"))

  return(df)
}

cor_inspect<-function(data=data_merged_GSE88715, gene1, gene2, comp1="_tis1", comp2="_tis2"){

  df<-data.frame(gene1=c(data[paste(gene1, comp1, sep=""),], data[paste(gene1, comp2, sep=""),], data[paste(gene1, comp1, sep=""),], data[paste(gene1, comp2, sep=""),]),
                 gene2=c(data[paste(gene2, comp1, sep=""),], data[paste(gene2, comp2, sep=""),], data[paste(gene2, comp2, sep=""),], data[paste(gene2, comp2, sep=""),]),
                 compartment=c(rep(c("comp1 vs comp1","comp2 vs comp2","comp1 vs comp2", "comp2 vs comp1"), each=ncol(data)) ))
  p<-ggplot(df, aes(x=gene1, y=gene2))+geom_point()+facet_wrap(.~compartment)+geom_smooth(method = "lm")+theme_classic()+labs(x=gene1, y=gene2)
  return(p)
}
