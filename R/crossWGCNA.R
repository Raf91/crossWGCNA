rm_selfloop <- function(A,comp1=comp1,comp2=comp2,verbose=TRUE)
{
  if(verbose) cat("Removing self-loops only...\n")
  genes_comp1 <- grep(comp1, rownames(A))
  genes_comp2 <- grep(comp2, rownames(A))
  genes1 <- gsub(comp1,"",rownames(A)[genes_comp1])
  genes2 <- gsub(comp2,"",rownames(A)[genes_comp2])
  Idx1 <- cbind(genes_comp1,genes_comp2[match(genes1,genes2)])
  Idx2 <- cbind(genes_comp2,genes_comp1[match(genes2,genes1)])
  A[Idx1] <- 0
  A[Idx2] <- 0
  if(verbose) cat("..Done!\n")
  return(A)
}

rm_netdiff <- function(A,comp1=comp1,comp2=comp2,verbose=TRUE)
{
  if(verbose) cat("Computing average conserved interactions...\n")
  genes_comp1 <- grep(comp1, rownames(A))
  genes_comp2 <- grep(comp2, rownames(A))
  A_orig <- A
  genes_comp1_orig <- grep(comp1, rownames(A_orig))
  genes_comp2_orig <- grep(comp2, rownames(A_orig))
  avgpath <- matrix(ncol=nrow(A)/2, nrow=nrow(A)/2)
  for(x in 1:(nrow(A)/2)){
    avgpath[x, ] <- (
      A[genes_comp1[x], genes_comp1]+
      A[genes_comp2[x], genes_comp2])/2
  }
  if(verbose) cat("Removing average conserved interactions...\n")
  A[genes_comp1, genes_comp2] <- A[genes_comp1, genes_comp2]-avgpath
  A[genes_comp2, genes_comp1] <- A[genes_comp2, genes_comp1]-avgpath
  diff <- abs(A_orig[genes_comp1_orig, genes_comp2_orig])-abs(A[genes_comp1, genes_comp2])
  A[genes_comp1, genes_comp2][diff<0] <- A_orig[genes_comp1_orig, genes_comp2_orig][diff<0]
  diff <- abs(A_orig[genes_comp2_orig, genes_comp1_orig])-abs(A[genes_comp2, genes_comp1])
  A[genes_comp2, genes_comp1][diff<0] <- A_orig[genes_comp2_orig, genes_comp1_orig][diff<0]
  A <- A/2
  if(verbose) cat("Removing self-loops...\n")
  A <- rm_selfloop(A,comp1=comp1,comp2=comp2,verbose=FALSE)
  if(verbose) cat("..Done!\n")
  return(A)
}

Adjacency <- function(
  data,
  method="selfloop",
  comp1="_1",
  comp2="_2",
  Adj_type="signed",
  cortype="spearman",
  pval="none",
  thr=0.05,
  beta=6,
  verbose=FALSE)
  {
    if(!(method %in% c("netdiff","selfloop"))){
      stop("Please select a valid method. Should be 'netdiff' or 'selfloop'.")
    }
    if(!(Adj_type %in% c("signed","unsigned","keep sign"))){
      stop("\n'Adj_type' argument is different from all the admitted values.\nPlease refer to manual for further details.")
    }
    if(!(cortype %in% c("pearson","spearman","bicor"))){
      stop("'cortype' argument is different from all the admitted values.\nPlease refer to manual for further details.")
    }
    if(!(pval %in% c("none","threshold","weight"))){
      stop("'pval' argument is different from all the admitted values.\nPlease refer to manual for further details.")
    }
    if(!is.numeric(thr)){
      stop("'thr' argument should be numeric. Default is 0.05.")
    }
    if(!is.numeric(beta)){
      stop("'beta' argument should be an integer. Default is 6.")
    }

    if(verbose){
      cat("Computing correlation matrix...\n")
    }

    if (pval=="none") {
      if (cortype=="bicor") {
        A <- bicor(t(data))
      } else {
        A <- cor(t(data), method = cortype)
      }
    } else {
      if (cortype=="bicor") {
        cat("Can't use pval with bicor\n")
      } else {
        mat <- rcorr(t(data), type=cortype)
        A <- mat[[1]]
        p.val <- mat[[3]]
        rm(mat)
        if (pval=="threhsold") {
          A[which(p.val>thr)] <- NA
        } else if (pval=="weight") {
          A <- A * (1-p.val)
        }
      }
    }
    if(verbose){
      cat("..Done!\n")
    }
    
    comp1 <- paste(comp1,"$",sep="")
    comp2 <- paste(comp2,"$",sep="")

    if(method=="selfloop"){
      A <- rm_selfloop(A,comp1=comp1,comp2=comp2)
    }
    
    if(method=="netdiff"){
      A <- rm_netdiff(A,comp1=comp1,comp2=comp2)
    }

    if(verbose){
      cat("Computing adjacency matrix...\n")
    }

    if (Adj_type=="signed"){
      A <- (0.5 * (1+A))^beta
    } else if (Adj_type=="unsigned"){
      A <- (abs(A))^beta
    } else if (Adj_type=="keep sign"){
      A <- ((abs(A))^beta)*sign(A)
    }

    if(verbose){
      cat("..Done!\n")
    }

    return(A)
  }


clusteringWGCNA <- function(A,data,comp1="_1",comp2="_2",TOM=TRUE,ds=1,crossOnly=TRUE)
{
  comp1 <- paste(comp1, "$", sep="")
  comp2 <- paste(comp2, "$", sep="")
  genes_comp1 <- grep(comp1, rownames(A))
  genes_comp2 <- grep(comp2, rownames(A))
  
  if (crossOnly){
    A[genes_comp1, genes_comp1] <- 0
    A[genes_comp2, genes_comp2] <- 0
  }

  A <- rm_selfloop(A,comp1=comp1,comp2=comp2,verbose=F)
    
  if(TOM){
    similarity <- TOMsimilarity(A, TOMType="signed")
    rownames(similarity) <- rownames(A)
    colnames(similarity) <- colnames(A)
    A <- similarity
    rm(similarity)
    if(crossOnly){
      A[genes_comp1, genes_comp1] <- 0
      A[genes_comp2, genes_comp2] <- 0
    }
  }

  conTree <- hclust(as.dist(1-A), method="average")
  unmergedLabels <- cutreeDynamic(dendro=conTree,dist=1-A,deepSplit=ds)
  names(unmergedLabels) <- rownames(A)
  merged <- mergeCloseModules(t(data),unmergedLabels,cutHeight=0.25,verbose=3)
  names(merged$colors) <- rownames(A)
  return(merged)
}


degrees <- function(A, comp1="_1", comp2="_2")
{
  comp1 <- paste(comp1, "$", sep="")
  comp2 <- paste(comp2, "$", sep="")
  kTot <- rowSums(A)
  kInt_1 <- rowSums(A[grep(comp1, rownames(A)), grep(comp1, colnames(A))])
  kInt_2 <- rowSums(A[grep(comp2, rownames(A)), grep(comp2, colnames(A))])
  kExt_1 <- rowSums(A[grep(comp1, rownames(A)), grep(comp2, colnames(A))])
  kExt_2 <- rowSums(A[grep(comp2, rownames(A)), grep(comp1, colnames(A))])
  k <- list(
    kInt1=kInt_1,
    kInt2=kInt_2,
    kExt1=kExt_1,
    kExt2=kExt_2,
    kTot1=kTot[grep(comp1, colnames(A))],
    kTot2=kTot[grep(comp2, colnames(A))]
  )
  return(k)
}

## whole crossWGCNA pipeline
crossWGCNA <- function(
  data,
  method="netdiff",
  Adj_type="signed",
  cortype="spearman",
  pval="none",
  thr=0.05,
  beta=6,
  comp1="_1",
  comp2="_2",
  doClusters=TRUE,
  doTOM=TRUE,
  ds=1,
  crossOnly=TRUE,
  verbose=TRUE)
  {
    comp1 <- paste(comp1, "$", sep = "")
    comp2 <- paste(comp2, "$", sep = "")
    
    Adj <- Adjacency(
      data=data,
      method=method,
      Adj_type=Adj_type,
      cortype=cortype,
      pval=pval,
      thr=thr,
      beta=beta,
      comp1=comp1,
      comp2=comp2,
      verbose=TRUE
    )

    if(verbose) cat("Computing intra- and inter-tissue connectivities...\n")
    k <- degrees(A=Adj,comp1=comp1,comp2=comp2)

    if(doClusters){
      if(verbose) cat("Computing clusters...\n")
      clusters <- clusteringWGCNA(
        A=Adj,
        data=data,
        comp1=comp1,
        comp2=comp2,
        TOM=doTOM,
        ds=ds,
        crossOnly=crossOnly
      )
      out <- list(k, clusters)
    }

    if(verbose) cat("..Done!\n")
    if(doClusters){
      return(out)
    } else {
      return(k)
    }
  }

changenames <- function(data, anno)
{
  annotation_sel <- anno[match(rownames(data), anno[, 1]), 2]
  if (length(which(annotation_sel==""))>0){
    data <- data[-which(annotation_sel==""), ]
    annotation_sel <- annotation_sel[-which(annotation_sel=="")]
  }
  a <- which(duplicated(annotation_sel))
  while(length(a)>0) 
  {
    for(i in 1:length(unique(annotation_sel))){
      if(length(which(annotation_sel==unique(annotation_sel)[i])) > 1) {
        m <- which.max(rowMeans(data[which(annotation_sel==unique(annotation_sel)[i]), ], na.rm =T))
        data <- data[-which(annotation_sel==unique(annotation_sel)[i])[-m], ]
        annotation_sel <- annotation_sel[-which(annotation_sel==unique(annotation_sel)[i])[-m]]
      }
    }
    data <- data[which(is.na(annotation_sel)==F), ]
    annotation_sel <- na.omit(annotation_sel)
    a <- which(duplicated(annotation_sel))
  }
  rownames(data)-annotation_sel
  return(data)
}

degrees_mod <- function(
  data,
  method="netdiff",
  modules,
  Adj_type="signed",
  cortype="spearman",
  pval="none",
  thr=0.05,
  beta=6,
  comp1="_1",
  comp2="_2") 
  {
    k <- list()
    for (i in 1:length(unique(modules))){
      cat(sprintf("Computing node degrees and k-ratio for module %i\n",i))
      mod <- names(modules)[which(modules==unique(modules)[i])]
      genes <- unique(gsub(comp2, "", gsub(comp1, "", mod)))
      Adj <- Adjacency(
        data=data[c(paste(genes, comp1, sep=""), paste(genes, comp2, sep="")),],
        method=method,
        Adj_type=Adj_type,
        cortype=cortype,
        pval=pval,
        thr=thr,
        beta=beta,
        comp1=comp1,
        comp2=comp2,
      )
      Adj <- Adj[mod, mod]
      k[[i]] <- degrees(A=Adj, comp1=comp1, comp2=comp2)
    }

   return(k)
  }

cytoscape_net <- function(A, data, gene, comp1, comp2, num, corr="spearman")
{
  interactors <- c(names(sort(A[paste(gene, comp1, sep="_"),grep(comp2, colnames(A))], decreasing=T))[1:num])
  inter_edges <- cor(data[paste(gene, comp1, sep="_"),],t(data[intersect(interactors, rownames(data)),]), method=corr)
  intra_1 <- cor(data[paste(gene, comp1, sep="_"),],t(data[gsub(comp2, comp1, intersect(interactors, rownames(data))),]),method=corr)
  intra_2 <- cor(data[paste(gene, comp2, sep="_"),],t(data[intersect(interactors, rownames(data)),]), method=corr)

  df <- data.frame(
    Source=c(
      rep(gene, length(inter_edges)),
      rep(gene, length(intra_1)),
      rep(gene, length(intra_2))),
    Target=c(
      gsub(paste("_", comp2, sep=""), "", colnames(inter_edges)),
      gsub(paste("_", comp1, sep=""), "", colnames(intra_1)),
      gsub(paste("_", comp2, sep=""), "", colnames(intra_2))),
    Weight=c(
      inter_edges,
      intra_1,
      intra_2),
    Edge_type=c(
      rep("inter", length(inter_edges)),
      rep("intra1", length(intra_1)),
      rep("intra2", length(intra_2))),
    Source_type=c(
      rep(comp1, length(inter_edges)),
      rep(comp1, length(intra_1)),
      rep(comp2,  length(intra_2))),
    Target_type=c(
      rep(comp2, length(inter_edges)),
      rep(comp1, length(intra_1)),
      rep(comp2,  length(intra_2)))
  )

  write.csv(df, file=paste("Cytoscape", gene, comp1, "egdes.csv", sep="_"))
  return(df)
}

cor_inspect <- function(data,gene1,gene2,comp1="_tis1",comp2="_tis2")
{
  df <- data.frame(
    gene1=c(
      data[paste(gene1, comp1, sep=""),], 
      data[paste(gene1, comp2, sep=""),], 
      data[paste(gene1, comp1, sep=""),], 
      data[paste(gene1, comp2, sep=""),]),
    gene2=c(
      data[paste(gene2, comp1, sep=""),], 
      data[paste(gene2, comp2, sep=""),], 
      data[paste(gene2, comp2, sep=""),], 
      data[paste(gene2, comp2, sep=""),]),
    compartment=c(
      rep(c(
        "comp1 vs comp1",
        "comp2 vs comp2",
        "comp1 vs comp2", 
        "comp2 vs comp1"), 
        each=ncol(data))))
  
  p <- ggplot(df, aes(x=gene1, y=gene2))+
    geom_point()+
    facet_wrap(.~compartment)+
    geom_smooth(method = "lm")+
    theme_classic()+
    labs(x=gene1, y=gene2)
  
  return(p)
}

###define spots coordinates
#data is the Seurat object

ST_spots_coords <- function(data, br=1000)
{
  y <- GetTissueCoordinates(data)[,1]
  p <- hist(y, breaks=br)
  #1000 is way higher than the number of peaks. If using another array, change this number
  #peaks separated by 0
  breaks <- p$breaks[p$counts==0]
  pos_break <- which(p$counts==0)
  breaks<-breaks[-which(pos_break %in% (pos_break+1))]
  peaks <- cut(y, breaks = c(min(y),breaks, max(y)))
  levels(peaks) <- c(1:length(unique(peaks)))
  y_bin <- peaks
  x <- GetTissueCoordinates(data)[,2]
  p <- hist(x, breaks=br)
  breaks <- p$breaks[p$counts==0]
  pos_break <- which(p$counts==0)
  breaks <- breaks[-which(pos_break %in% (pos_break+1))]
  peaks <- cut(x, breaks = c(min(x),breaks, max(x)))
  levels(peaks) <- c(1:length(unique(peaks)))
  x_bin <- peaks
  x_bin <- as.numeric(x_bin)
  y_bin <- as.numeric(y_bin)
  #transform ycoordinates so that they represent the height of an equilateral triangle
  y_bin <- (y_bin*sqrt(3))
  coords <- cbind(x_bin, y_bin)
  return(coords)
}

###smooths gene expression using the weighted mean of neighbouring spots in the same compartment
#spots_class class of spots represented in expr_data, same order
#coords, output of spots_coords
#spots_dist

ST_expr_smooth <- function(expr_data, coords, max_dist=5, spots_class, sel_class=c("Epi", "Stroma"))
{
  spots_dist <- dist(coords)
  averaged_expr_all <- matrix(ncol=ncol(expr_data), nrow=nrow(expr_data))
  for(es in 1:ncol(expr_data)){
    class_es <- spots_class[es]
    if(class_es %in% sel_class){
      dist_es <- as.matrix(spots_dist)[,es]
      sel_spots <- which(dist_es<max_dist & spots_class==class_es)
      weights <- 1/(dist_es[sel_spots]+1)
      if(length(sel_spots)>1){
        averaged_expr_all[,es] <- apply(
          (exp(expr_data[,sel_spots])-1), 1, function(x){log(weighted.mean(x, weights)+1)})
      } else {
        averaged_expr_all[,es] <- expr_data[,sel_spots]
      }
    }
  }
  rownames(averaged_expr_all) <- rownames(expr_data)
  return(averaged_expr_all)
}

###selects epi and stroma spots not isolated
##epi with at least 1 neighbouring epi spot
##stroma with at least 1 neighbouring stroma spot
##epi with at least 1 selected sroma spot

ST_spots_filt <- function(coords, tis1_spots, tis2_spots)
{  
  x_bin <- coords[,1]
  y_bin <- coords[,2]
  
  included_es_spots <- c()
  
  for(es in tis1_spots){
    which_x_coord <- which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
    which_y_coord <- which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
    sel_spots <- intersect(which_x_coord, which_y_coord)
    neighbour_spots_epi <- intersect(tis1_spots, sel_spots)
    if(length(neighbour_spots_epi)>1){
      included_es_spots <- c(included_es_spots, es)
    }
  }

  included_ss_spots <- c()
  
  for(ss in tis2_spots){
    which_x_coord <- which(x_bin>=x_bin[ss]-2 & x_bin<=x_bin[ss]+2)
    which_y_coord <- which(y_bin>=y_bin[ss]-sqrt(3)-0.1 & y_bin<=y_bin[ss]+sqrt(3)+0.1)
    sel_spots <- intersect(which_x_coord, which_y_coord)
    neighbour_spots_stroma <- intersect(tis2_spots, sel_spots)
    if(length(neighbour_spots_stroma)>1){
      included_ss_spots <- c(included_ss_spots, ss)
    }
  }
  return(list(included_es_spots, included_ss_spots))
}

#############
###takes the stromal spots next to each epi spot and
#creates a matched epi and stroma expression matrices
#averages the expression of stroma spots in contact with a specific epithelial spot
#averaged_expr_all output of expr_smooth
#coords output of spots_coords
#sel_spots output of spots_filt
#stroma spots not the same as sel_ss_spots, is the output identical?
#tis1 here is the epithelium, so cor consistency with the rest of the scripts comp1 and comp2 should be inverted

ST_merged_dataset <- function(sel_spots, coords, averaged_expr_all, var_thr=0.75, comp1="_tis1", comp2="_tis2")
{
  included_es_spots <- sel_spots[[1]]
  included_ss_spots <- sel_spots[[2]]
  x_bin <- coords[,1]
  y_bin <- coords[,2]
  sel_es_spots <- c()
  for (es in included_es_spots){
    which_x_coord <- which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
    which_y_coord <- which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
    sel_spots <- intersect(which_x_coord, which_y_coord)
    neighbour_spots <- intersect(included_ss_spots, sel_spots)
    if(length(neighbour_spots)>0){
      sel_es_spots <- c(sel_es_spots, es)
    }
  }

  included_spots <- c()
  tis1_expr_all <- matrix(ncol=1, nrow=nrow(averaged_expr_all))
  tis2_expr_all <- matrix(ncol=1, nrow=nrow(averaged_expr_all))

  for(es in sel_es_spots){
    which_x_coord <- which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
    which_y_coord <- which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
    sel_spots <- intersect(which_x_coord, which_y_coord)
    sel_spots_tis2 <- intersect(included_ss_spots, sel_spots)
    neighbour_spots_tis1 <- intersect(sel_es_spots, sel_spots)

    if (length(sel_spots_tis2)>0 & length(neighbour_spots_tis1)>1){
      if (length(sel_spots_tis2)==1){
        tis2_expr <- averaged_expr_all[,sel_spots_tis2]
      } else if (length(sel_spots_tis2)>1 ){
        tis2_expr <- log(rowMeans(exp(averaged_expr_all[,sel_spots_tis2])-1)+1)
      }

      tis1_expr <- averaged_expr_all[,es]
      tis1_expr_all <- cbind(tis1_expr_all, tis1_expr)
      tis2_expr_all <- cbind(tis2_expr_all, tis2_expr)
      included_spots <- c(included_spots, es)
    }
  }
  
  rownames(tis1_expr_all) <- rownames(averaged_expr_all)
  rownames(tis2_expr_all) <- rownames(averaged_expr_all)
  
  tis1_expr_all <- tis1_expr_all[,-1]
  tis2_expr_all <- tis2_expr_all[,-1]

  var_tis2 <- apply(tis2_expr_all,1,var)
  var_tis1 <- apply(tis1_expr_all,1,var)

  genes <- intersect(
    which(var_tis2>quantile(var_tis2, var_thr, na.rm=T)), 
    which(var_tis1>quantile(var_tis1, var_thr, na.rm=T)))

  genes <- rownames(averaged_expr_all)[genes]
  tis2 <- tis2_expr_all[genes,]
  tis1 <- tis1_expr_all[genes,]

  #merge tis2 and tis1
  rownames(tis2) <- paste(rownames(tis2), comp2, sep="")
  rownames(tis1) <- paste(rownames(tis1), comp1, sep="")
  colnames(tis1) <- colnames(tis2)
  data_merged <- rbind(tis2, tis1)
  return(list(data_merged, included_spots))
}

##finds spots included in the boundaries
#included_spots output of merged_dataset [[2]]

ST_boundary_spots <- function(included_spots, coords, tis2_spots)
{  
  x_bin <- coords[,1]
  y_bin <- coords[,2]
  
  included_spots_tis1 <- c()
  included_spots_tis2 <- c()
  
  for(es in included_spots){
    which_x_coord <- which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
    which_y_coord <- which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
    sel_spots <- intersect(which_x_coord, which_y_coord)
    sel_spots_tis2 <- intersect(tis2_spots, sel_spots)
    included_spots_tis2 <- c(included_spots_tis2, sel_spots_tis2)
    included_spots_tis1 <- c(included_spots_tis1, rep(es, length(sel_spots_tis2)))
  }
  return(list(included_spots_tis1, included_spots_tis2))
}

##midpoints coordinates
#output of boundary_spots
ST_midpoints_def <- function(coords, sel_spots)
{  
  x_bin <- coords[,1]
  y_bin <- coords[,2]
  
  included_spots_tis1 <- sel_spots[[1]]
  included_spots_tis2 <- sel_spots[[2]]
  
  df <- data.frame(
    x_coord=c(
      x_bin[included_spots_tis1],
      x_bin[included_spots_tis2]),
    y_coord=c(
      -(y_bin[included_spots_tis1]),
      -(y_bin[included_spots_tis2])))

  midpoints_x <- (df$x_coord[c(1:length(included_spots_tis1))]+df$x_coord[-c(1:length(included_spots_tis1))])/2
  midpoints_y <- (df$y_coord[c(1:length(included_spots_tis1))]+df$y_coord[-c(1:length(included_spots_tis1))])/2
  return(list(midpoints_x, midpoints_y))
}

###visualize filtered spots
#plot(as.numeric(x_bin[included_es_spots]),-as.numeric(y_bin[included_es_spots]), cex=0.5, pch=19, xlab="x", ylab="y")
#points(as.numeric(x_bin[included_ss_spots]),-as.numeric(y_bin[included_ss_spots]), cex=0.5, pch=19, col="red")

###visualize gene expression in space
#midpoints from midpoints_def
ST_plot_expr <- function(
  gene, 
  averaged_expr_all, 
  coords, 
  included_spots, 
  tis1_spots, 
  tis2_spots, 
  midpoints)
  {
    midpoints_x <- midpoints[[1]]
    midpoints_y <- midpoints[[2]]

    df <- data.frame(
      x_coord= c(
        x_bin[tis1_spots],
        x_bin[tis2_spots], 
        x_bin[included_spots]),
      y_coord=c(
        -(y_bin[tis1_spots]), 
        -(y_bin[tis2_spots]), 
        -y_bin[included_spots]),
      compartment=c(
        rep("tis1", length(tis1_spots)),
        rep("tis2", length(tis2_spots)), 
        rep("edge", length(included_spots))),
      gene=c(
        averaged_expr_all[gene,tis1_spots], 
        averaged_expr_all[gene,tis2_spots], 
        averaged_expr_all[gene,included_spots]))
  
  df_midpoint <- data.frame(x_coord= midpoints_x,y_coord=midpoints_y)

  p <- ggplot(
    data=df, aes(x=x_coord, y=y_coord, colour=gene))+
    geom_point()+scale_color_continuous(type = "viridis")+
    geom_point(data=df_midpoint, size=1, colour="black")+
    theme_classic()

  return(p)
  }

##visualize gene communication in space
#gene1 measured in tis1
#averaged_expr_all output of expr_smooth
#coords output of spots_coords
#midpoints from midpoints_def

ST_plot_comm <- function(
  gene1, 
  gene2, 
  averaged_expr_all, 
  coords, 
  included_spots, 
  sel_spots, 
  tis1_spots, 
  tis2_spots, 
  midpoints)
  {
    midpoints_x <- midpoints[[1]]
    midpoints_y <- midpoints[[2]]

    x_bin <- coords[,1]
    y_bin <- coords[,2]

    included_spots_tis1 <- sel_spots[[1]]
    included_spots_tis2 <- sel_spots[[2]]
  
    df <- data.frame(
      x_coord=c(
        x_bin[tis1_spots],
        x_bin[tis2_spots], 
        x_bin[included_spots]),
      y_coord=c(
        -(y_bin[tis1_spots]), 
        -(y_bin[tis2_spots]), 
        -y_bin[included_spots]),
      compartment=c(
        rep("tis1", length(tis1_spots)),
        rep("tis2", length(tis2_spots)), 
        rep("edge", length(included_spots))),
      gene1=c(
        averaged_expr_all[gene1,tis1_spots], 
        averaged_expr_all[gene1,tis2_spots], 
        averaged_expr_all[gene1,included_spots]),
      gene2=c(
        averaged_expr_all[gene2,tis1_spots], 
        averaged_expr_all[gene2,tis2_spots], 
        averaged_expr_all[gene2,included_spots])
    )

    df$gene1.scale <- with(df, (gene1-min(gene1, na.rm=T))/diff(range(gene1, na.rm=T)))
    df$gene2.scale <- with(df, (gene2-min(gene2, na.rm=T))/diff(range(gene2, na.rm=T)))
    df <- df[!is.na(df$gene1.scale) & !is.na(df$gene2.scale),]

    df_midpoint <- data.frame(
      x_coord=midpoints_x,
      y_coord=midpoints_y,
      comm_score=averaged_expr_all[gene1,included_spots_tis1]*averaged_expr_all[gene2,included_spots_tis2])
        
    p <- ggplot(
      data=df, 
      aes(x=x_coord, y=y_coord))+
      geom_point(data=df_midpoint, size=1, aes(colour=comm_score))+
      scale_color_continuous(type = "viridis")+
      theme_classic()
      
    return(p)
  }

## compute weighted average of a module
# modules vector with module labels assignments
ST_weighted_mod <- function(modules,kw,mod_sel,averaged_expr_all,comp1="_tis1",comp2="_tis2")
{
  mod <- names(modules)[which(modules==mod_sel)]
  weights <- kw[[which(unique(modules)==mod_sel)]]$kExt1[mod[grep(comp1, mod)]]
  wm1 <- apply((averaged_expr_all[gsub(comp1, "", mod[grep(comp1, mod)]),]), 2, function(x){weighted.mean(x, weights)})
  weights <- kw[[which(unique(modules)==mod_sel)]]$kExt2[mod[grep(comp2, mod)]]
  wm2 <- apply((averaged_expr_all[gsub(comp2, "", mod[grep(comp2, mod)]),]), 2, function(x){weighted.mean(x, weights)})
  return(list(wm1, wm2))
}