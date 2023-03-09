
Adjacency <-
  function(data,
           comp1 = "_1",
           comp2 = "_2",
           Adj_type = "signed",
           cortype = "spearman",
           pval = "none",
           thr = 0.05,
           beta = 6,
           sign_list = 1,
           compartment_sel = "none" ,
           selgenes = NA) {

    comp1 <- paste(comp1, "$", sep = "")
    comp2 <- paste(comp2, "$", sep = "")

    ###Creates the correlation matrix A
    if (pval == "none") {
      if (cortype == "bicor") {
        A <- bicor(t(data))
      } else {
        A <- cor(t(data), method = cortype)
      }
    } else {
      if (cortype == "bicor") {
        paste("Can't use pval with bicor")
      } else {
        mat <- rcorr(t(data), type = cortype)
        A <- mat[[1]]
        p.val <- mat[[3]]
        rm(mat)
        if (pval == "threhsold") {
          A[which(p.val > thr)] <- NA
        } else if (pval == "weight") {
          A <- A * (1 - p.val)
        }
      }
    }
    A_orig<-A
    #computes average conserved interactions
    #assumes that genes are in the same order in the two compartments
    avgpath <- matrix(ncol = nrow(A) / 2, nrow = nrow(A) / 2)
    for (x in 1:(nrow(A) / 2)) {
      avgpath[x, ] <-
        (A[grep(comp1, rownames(A))[x], grep(comp1, rownames(A))] + A[grep(comp2, rownames(A))[x], grep(comp2, rownames(A))]) /
        2
    }

    #removes average conserved interactions
    A[grep(comp1, rownames(A)), grep(comp2, rownames(A))] <-
      A[grep(comp1, rownames(A)), grep(comp2, rownames(A))] - avgpath
    A[grep(comp2, rownames(A)), grep(comp1, rownames(A))] <-
      A[grep(comp2, rownames(A)), grep(comp1, rownames(A))] - avgpath

    ##take the lowest absolute value
    diff<-abs(A_orig[grep(comp1, rownames(A_orig)), grep(comp2, rownames(A_orig))])-abs(A[grep(comp1, rownames(A)), grep(comp2, rownames(A))])
    A[grep(comp1, rownames(A)), grep(comp2, rownames(A))][diff<0]<-A_orig[grep(comp1, rownames(A_orig)), grep(comp2, rownames(A_orig))][diff<0]

    diff<-abs(A_orig[grep(comp2, rownames(A_orig)), grep(comp1, rownames(A_orig))])-abs(A[grep(comp2, rownames(A)), grep(comp1, rownames(A))])
    A[grep(comp2, rownames(A)), grep(comp1, rownames(A))][diff<0]<-A_orig[grep(comp2, rownames(A_orig)), grep(comp1, rownames(A_orig))][diff<0]


    ###make A ranges between -1 and 1
    A <- A / 2

    ####when there is a network with sign in the input
    if (compartment_sel == "comp2") {
      if (is.na(selgenes)) {
        stop("No selected genes provided")
      }

      sign_list <- sign_list[which(!is.na(sign_list))]

      names(sign_list) <-
        paste(names(sign_list), gsub("\\$", "", comp2), sep = "")

      selgenes <-
        intersect(paste(selgenes, gsub("\\$", "", comp2), sep = ""), rownames(A))

      A <-
        A[c(grep(comp1, rownames(A)), which(rownames(A) %in% selgenes)), c(grep(comp1, rownames(A)), which(rownames(A) %in% selgenes))]

      A[grep(comp2, rownames(A)), grep(comp1, rownames(A))] <-
        A[grep(comp2, rownames(A)), grep(comp1, rownames(A))] * (sign_list[rownames(A)[grep(comp2, rownames(A))]])

      A[grep(comp1, rownames(A)), grep(comp2, rownames(A))] <-
        t(A[grep(comp2, rownames(A)), grep(comp1, rownames(A))] * (sign_list[rownames(A)[grep(comp2, rownames(A))]]))

      } else if (compartment_sel == "comp1") {

        if (is.na(selgenes)) {
        stop("No selected genes provided")
        }

      sign_list <- sign_list[which(!is.na(sign_list))]

      names(sign_list) <-
        paste(names(sign_list), gsub("\\$", "", comp1), sep = "")

      selgenes <-
        intersect(paste(selgenes, gsub("\\$", "", comp1), sep = ""), rownames(A))

      A <-
        A[c(which(rownames(A) %in% selgenes), grep(comp2, rownames(A))), c(which(rownames(A) %in% selgenes), grep(comp2, rownames(A)))]

      A[grep(comp1, rownames(A)), grep(comp2, rownames(A))] <-
        A[grep(comp1, rownames(A)), grep(comp2, rownames(A))] * (sign_list[rownames(A)[grep(comp1, rownames(A))]])

      A[grep(comp2, rownames(A)), grep(comp1, rownames(A))] <-
        t(A[grep(comp1, rownames(A)), grep(comp2, rownames(A))] * (sign_list[rownames(A)[grep(comp1, rownames(A))]]))
    }

    ####

    if (Adj_type == "signed") {
      A <- (0.5 * (1 + A)) ^ beta
    } else if (Adj_type == "unsigned") {
      A <- (abs(A)) ^ beta
    } else if (Adj_type == "keep sign") {
      A <- ((abs(A)) ^ beta) * sign(A)
    }
    return(A)
  }

##clustering with WGCNA functions on pre-computed Adjacency
clusteringWGCNA <-
  function(A,
           data,
           comp1 = "_1",
           comp2 = "_2",
           TOM = T,
           ds = 1,
           crossOnly = T) {

    comp1 <- paste(comp1, "$", sep = "")
    comp2 <- paste(comp2, "$", sep = "")

    ##crossOnly == T sets intra-connectivities to 0
    if (crossOnly == T) {
      A[grep(comp1, rownames(A)), grep(comp1, rownames(A))] <- 0
      A[grep(comp2, rownames(A)), grep(comp2, rownames(A))] <- 0
    }

    #remove self loops
    genes1 <- gsub(comp1, "", rownames(A)[grep(comp1, rownames(A))])
    genes2 <- gsub(comp2, "", rownames(A)[grep(comp2, rownames(A))])

    Idx1 <-
      cbind(grep(comp1, rownames(A)), grep(comp2, rownames(A))[match(genes1, genes2)])
    A[Idx1] <- 0

    Idx2 <-
      cbind(grep(comp2, rownames(A)), grep(comp1, rownames(A))[match(genes2, genes1)])
    A[Idx2] <- 0

    if (TOM == T) {
      similarity <- TOMsimilarity(A, TOMType = "signed")
      rownames(similarity) <- rownames(A)
      colnames(similarity) <- colnames(A)
      A <- similarity
      rm(similarity)
    }

    ##crossOnly == T sets intra-connectivities to 0
    if (crossOnly == T) {
      A[grep(comp1, rownames(A)), grep(comp1, rownames(A))] <- 0
      A[grep(comp2, rownames(A)), grep(comp2, rownames(A))] <- 0
    }

    conTree = hclust(as.dist(1 - A), method = "average")
    unmergedLabels = cutreeDynamic(dendro = conTree,
                                   dist = 1 - A,
                                   deepSplit = ds)
    names(unmergedLabels) <- rownames(A)
    merged = mergeCloseModules(t(data),
                               unmergedLabels,
                               cutHeight = 0.25,
                               verbose = 3)
    names(merged$colors) <- rownames(A)

    return(merged)
  }

###computes intra- and inter-tissue connectivities
degrees <- function(A, comp1 = "_1", comp2 = "_2") {
  comp1 <- paste(comp1, "$", sep = "")
  comp2 <- paste(comp2, "$", sep = "")

  #remove self loops
  genes1 <- gsub(comp1, "", rownames(A)[grep(comp1, rownames(A))])
  genes2 <- gsub(comp2, "", rownames(A)[grep(comp2, rownames(A))])

  Idx1 <-
    cbind(grep(comp1, rownames(A)), grep(comp2, rownames(A))[match(genes1, genes2)])
  A[Idx1] <- 0

  Idx2 <-
    cbind(grep(comp2, rownames(A)), grep(comp1, rownames(A))[match(genes2, genes1)])
  A[Idx2] <- 0

  kTot <- rowSums(A)
  kInt_1 <-
    rowSums(A[grep(comp1, rownames(A)), grep(comp1, colnames(A))])
  kInt_2 <-
    rowSums(A[grep(comp2, rownames(A)), grep(comp2, colnames(A))])
  kExt_1 <-
    rowSums(A[grep(comp1, rownames(A)), grep(comp2, colnames(A))])
  kExt_2 <-
    rowSums(A[grep(comp2, rownames(A)), grep(comp1, colnames(A))])

  k <-
    list(
      kInt1 = kInt_1,
      kInt2 = kInt_2,
      kExt1 = kExt_1,
      kExt2 = kExt_2,
      kTot1 = kTot[grep(comp1, colnames(A))],
      kTot2 = kTot[grep(comp2, colnames(A))]
    )

  return(k)
}


##whole crossWGCNA pipeline
crossWGCNA <-
  function(data,
           Adj_type = "signed",
           cortype = "spearman",
           pval = "none",
           thr = 0.05,
           beta = 6,
           comp1 = "_1",
           comp2 = "_2",
           doTOM = T,
           ds = 1,
           crossOnly = T,
           sign_list = 1,
           compartment_sel = "none",
           selgenes = NA) {

    comp1 <- paste(comp1, "$", sep = "")
    comp2 <- paste(comp2, "$", sep = "")

    Adj <-
      Adjacency(
        data = data,
        Adj_type = Adj_type,
        cortype = cortype,
        pval = pval,
        thr = thr,
        beta = beta,
        comp1 = comp1,
        comp2 = comp2,
        sign_list = sign_list,
        compartment_sel = compartment_sel,
        selgenes = selgenes
      )

    k <- degrees(A = Adj, comp1 = comp1, comp2 = comp2)

    clusters <-
      clusteringWGCNA(
        A = Adj,
        data = data,
        comp1 = comp1,
        comp2 = comp2,
        TOM = doTOM,
        ds = ds,
        crossOnly = crossOnly
      )

    net_out <- list(k, clusters)
    return(net_out)
  }

##crossWGCNA pipeline avoiding clustering (output: kInt, kExt, kTot)
network <-
  function(data,
           Adj_type = "signed",
           cortype = "spearman",
           pval = "none",
           thr = 0.05,
           beta = 6,
           comp1 = "_1",
           comp2 = "_2",
           sign_list = 1,
           compartment_sel = "none",
           selgenes = NA) {

    comp1 <- paste(comp1, "$", sep = "")
    comp2 <- paste(comp2, "$", sep = "")

    Adj <-
      Adjacency(
        data = data,
        Adj_type = Adj_type,
        cortype = cortype,
        pval = pval,
        thr = thr,
        beta = beta,
        comp1 = comp1,
        comp2 = comp2,
        sign_list = sign_list,
        compartment_sel = compartment_sel,
        selgenes = selgenes
      )

    k <- degrees(A = Adj, comp1 = comp1, comp2 = comp2)
    return(k)
  }


changenames <- function(data, anno) {
  annotation_sel = anno[match(rownames(data), anno[, 1]), 2]

  if (length(which(annotation_sel == "")) > 0) {
    data <- data[-which(annotation_sel == ""), ]
    annotation_sel <- annotation_sel[-which(annotation_sel == "")]
  }

  a <- which(duplicated(annotation_sel))
  while (length(a) > 0) {
    for (i in 1:length(unique(annotation_sel))) {
      if (length(which(annotation_sel == unique(annotation_sel)[i])) > 1) {
        m = which.max(rowMeans(data[which(annotation_sel == unique(annotation_sel)[i]), ], na.rm =
                                 T))
        data = data[-which(annotation_sel == unique(annotation_sel)[i])[-m], ]
        annotation_sel = annotation_sel[-which(annotation_sel == unique(annotation_sel)[i])[-m]]
      }
    }

    data = data[which(is.na(annotation_sel) == F), ]
    annotation_sel = na.omit(annotation_sel)
    a <- which(duplicated(annotation_sel))
  }

  rownames(data) = annotation_sel
  return(data)
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
