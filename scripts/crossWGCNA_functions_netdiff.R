
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
      cbind(grep(comp2, rownames(A)), grep(comp1, rownames(A))[match(genes1, genes2)])
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
    cbind(grep(comp2, rownames(A)), grep(comp1, rownames(A))[match(genes1, genes2)])
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

