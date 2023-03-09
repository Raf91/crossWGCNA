###define spots coordinates
#data is the Seurat object
spots_coords<-function(data, br=1000){
  require(Seurat)
  y<-GetTissueCoordinates(data)[,1]
  p<-hist(y, breaks=br)
  #1000 is way higher than the number of peaks. If using another array, change this number
  #peaks separated by 0
  breaks<-p$breaks[p$counts==0]
  pos_break<-which(p$counts==0)
  breaks<-breaks[-which(pos_break %in% (pos_break+1))]
  peaks<-cut(y, breaks = c(min(y),breaks, max(y)))
  levels(peaks)<-c(1:length(unique(peaks)))
  y_bin<-peaks


  x<-GetTissueCoordinates(data)[,2]
  p<-hist(x, breaks=br)
  breaks<-p$breaks[p$counts==0]
  pos_break<-which(p$counts==0)
  breaks<-breaks[-which(pos_break %in% (pos_break+1))]
  peaks<-cut(x, breaks = c(min(x),breaks, max(x)))
  levels(peaks)<-c(1:length(unique(peaks)))
  x_bin<-peaks

  x_bin<-as.numeric(x_bin)
  y_bin<-as.numeric(y_bin)
  #transform ycoordinates so that they represent the height of an equilateral triangle
  y_bin<-(y_bin*sqrt(3))
  coords<-cbind(x_bin, y_bin)
  return(coords)
}

###smooths gene expression using the weighted mean of neighbouring spots in the same compartment
#spots_class class of spots represented in expr_data, same order
#coords, output of spots_coords
#spots_dist
#
expr_smooth<-function(expr_data, coords, max_dist=5, spots_class, sel_class=c("Epi", "Stroma")){
  spots_dist<-dist(coords)
  averaged_expr_all<-matrix(ncol=ncol(expr_data), nrow=nrow(expr_data))

for(es in 1:ncol(expr_data)){
  class_es<-spots_class[es]
  if(class_es %in% sel_class){
    dist_es<-as.matrix(spots_dist)[,es]
    sel_spots<-which(dist_es<max_dist & class==class_es)
    weights<-1/(dist_es[sel_spots]+1)
    if(length(sel_spots)>1){
      averaged_expr_all[,es]<-apply((exp(expr_data[,sel_spots])-1), 1, function(x){log(weighted.mean(x, weights)+1)})
    } else {
      averaged_expr_all[,es]<-expr_data[,sel_spots]
    }
  }
}
rownames(averaged_expr_all)<-rownames(expr_data)
return(averaged_expr_all)
}

###selects epi and stroma spots not isolated
##epi with at least 1 neighbouring epi spot
##stroma with at least 1 neighbouring stroma spot
##epi with at least 1 selected sroma spot
spots_filt<-function(coords, tis1_spots, tis2_spots){
  x_bin<-coords[,1]
  y_bin<-coords[,2]

included_es_spots<-c()
for(es in tis1_spots){
  which_x_coord<-which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
  which_y_coord<-which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
  sel_spots<-intersect(which_x_coord, which_y_coord)
  neighbour_spots_epi<-intersect(tis1_spots, sel_spots)
  if(length(neighbour_spots_epi)>1){
    included_es_spots<-c(included_es_spots, es)
  }
}


included_ss_spots<-c()
for(ss in tis2_spots){
  which_x_coord<-which(x_bin>=x_bin[ss]-2 & x_bin<=x_bin[ss]+2)
  which_y_coord<-which(y_bin>=y_bin[ss]-sqrt(3)-0.1 & y_bin<=y_bin[ss]+sqrt(3)+0.1)
  sel_spots<-intersect(which_x_coord, which_y_coord)
  neighbour_spots_stroma<-intersect(tis2_spots, sel_spots)
  if(length(neighbour_spots_stroma)>1){
    included_ss_spots<-c(included_ss_spots, ss)
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

merged_dataset<-function(sel_spots, coords, averaged_expr_all, var_thr=0.75, comp1="_tis1", comp2="_tis2"){
  included_es_spots<-sel_spots[[1]]
  included_ss_spots<-sel_spots[[2]]

  x_bin<-coords[,1]
  y_bin<-coords[,2]

  sel_es_spots<-c()
  for(es in included_es_spots){
    which_x_coord<-which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
    which_y_coord<-which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
    sel_spots<-intersect(which_x_coord, which_y_coord)
    neighbour_spots<-intersect(included_ss_spots, sel_spots)
    if(length(neighbour_spots)>0){
      sel_es_spots<-c(sel_es_spots, es)
    }
  }

included_spots<-c()
tis1_expr_all<-matrix(ncol=1, nrow=nrow(averaged_expr_all))
tis2_expr_all<-matrix(ncol=1, nrow=nrow(averaged_expr_all))

for(es in sel_es_spots){
  which_x_coord<-which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
  which_y_coord<-which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
  sel_spots<-intersect(which_x_coord, which_y_coord)
  sel_spots_tis2<-intersect(included_ss_spots, sel_spots)
  neighbour_spots_tis1<-intersect(sel_es_spots, sel_spots)

  if(length(sel_spots_tis2)>0 & length(neighbour_spots_tis1)>1){
    if(length(sel_spots_tis2)==1){
      tis2_expr<-averaged_expr_all[,sel_spots_tis2]
    } else if(length(sel_spots_tis2)>1 ){
      tis2_expr<-log(rowMeans(exp(averaged_expr_all[,sel_spots_tis2])-1)+1)
    }

    tis1_expr<-averaged_expr_all[,es]
    tis1_expr_all<-cbind(tis1_expr_all, tis1_expr)
    tis2_expr_all<-cbind(tis2_expr_all, tis2_expr)

    included_spots<-c(included_spots, es)
  }

}
rownames(tis1_expr_all)<-rownames(averaged_expr_all)
rownames(tis2_expr_all)<-rownames(averaged_expr_all)
tis1_expr_all<-tis1_expr_all[,-1]
tis2_expr_all<-tis2_expr_all[,-1]

var_tis2<-apply(tis2_expr_all,1,var)
var_tis1<-apply(tis1_expr_all,1,var)

genes<-intersect(which(var_tis2>quantile(var_tis2, var_thr, na.rm=T)),
                 which(var_tis1>quantile(var_tis1, var_thr, na.rm=T)))
genes<-rownames(averaged_expr_all)[genes]

tis2<-tis2_expr_all[genes,]
tis1<-tis1_expr_all[genes,]

#merge tis2 and tis1
rownames(tis2)<-paste(rownames(tis2), comp2, sep="")
rownames(tis1)<-paste(rownames(tis1), comp1, sep="")
colnames(tis1)<-colnames(tis2)

data_merged<-rbind(tis2, tis1)

return(list(data_merged, included_spots))
}


##finds spots included in the boundaries
#included_spots output of merged_dataset [[2]]
boundary_spots<-function(included_spots, coords, tis2_spots){
  x_bin<-coords[,1]
  y_bin<-coords[,2]

  included_spots_tis1<-c()
  included_spots_tis2<-c()
  for(es in included_spots){
    which_x_coord<-which(x_bin>=x_bin[es]-2 & x_bin<=x_bin[es]+2)
    which_y_coord<-which(y_bin>=y_bin[es]-sqrt(3)-0.1 & y_bin<=y_bin[es]+sqrt(3)+0.1)
    sel_spots<-intersect(which_x_coord, which_y_coord)
    sel_spots_tis2<-intersect(tis2_spots, sel_spots)
    included_spots_tis2<-c(included_spots_tis2, sel_spots_tis2)
    included_spots_tis1<-c(included_spots_tis1, rep(es, length(sel_spots_tis2)))
  }

  return(list(included_spots_tis1, included_spots_tis2))
}


##midpoints coordinates
#output of boundary_spots
midpoints_def<-function(coords, sel_spots){
  x_bin<-coords[,1]
  y_bin<-coords[,2]

  included_spots_tis1<-sel_spots[[1]]
  included_spots_tis2<-sel_spots[[2]]

df<-data.frame(x_coord= c(x_bin[included_spots_tis1],x_bin[included_spots_tis2]),y_coord= c(-(y_bin[included_spots_tis1]), -(y_bin[included_spots_tis2])))

midpoints_x<-(df$x_coord[c(1:length(included_spots_tis1))]+df$x_coord[-c(1:length(included_spots_tis1))])/2
midpoints_y<-(df$y_coord[c(1:length(included_spots_tis1))]+df$y_coord[-c(1:length(included_spots_tis1))])/2
return(list(midpoints_x, midpoints_y))
}

###visualize filtered spots
#plot(as.numeric(x_bin[included_es_spots]),-as.numeric(y_bin[included_es_spots]), cex=0.5, pch=19, xlab="x", ylab="y")
#points(as.numeric(x_bin[included_ss_spots]),-as.numeric(y_bin[included_ss_spots]), cex=0.5, pch=19, col="red")

###visualize gene expression in space
#midpoints from midpoints_def
plot_expr<-function(gene, averaged_expr_all, coords, included_spots, tis1_spots, tis2_spots, midpoints){
  require(ggplot2)

  midpoints_x<-midpoints[[1]]
  midpoints_y<-midpoints[[2]]

df<-data.frame(x_coord= c(x_bin[tis1_spots],x_bin[tis2_spots], x_bin[included_spots]),y_coord= c(-(y_bin[tis1_spots]), -(y_bin[tis2_spots]), -y_bin[included_spots]),
               compartment=c(rep("tis1", length(tis1_spots)),rep("tis2", length(tis2_spots)), rep("edge", length(included_spots))),
               gene=c(averaged_expr_all[gene,tis1_spots], averaged_expr_all[gene,tis2_spots], averaged_expr_all[gene,included_spots]))
df_midpoint<-data.frame(x_coord= midpoints_x,
                        y_coord= midpoints_y)

p<-ggplot(data=df, aes(x=x_coord, y=y_coord, colour=gene))+geom_point()+
  scale_color_continuous(type = "viridis")+
  geom_point(data=df_midpoint, size=1, colour="black")+theme_classic()

return(p)
}

##visualize gene communication in space
#gene1 measured in tis1
#averaged_expr_all output of expr_smooth
#coords output of spots_coords
#midpoints from midpoints_def

plot_comm<-function(gene1, gene2, averaged_expr_all, coords, included_spots, sel_spots, tis1_spots, tis2_spots, midpoints){
  require(ggplot2)
  require(ggnewscale)

  midpoints_x<-midpoints[[1]]
  midpoints_y<-midpoints[[2]]

  x_bin<-coords[,1]
  y_bin<-coords[,2]

  included_spots_tis1<-sel_spots[[1]]
  included_spots_tis2<-sel_spots[[2]]



  df<-data.frame(x_coord= c(x_bin[tis1_spots],x_bin[tis2_spots], x_bin[included_spots]),y_coord= c(-(y_bin[tis1_spots]), -(y_bin[tis2_spots]), -y_bin[included_spots]),
                 compartment=c(rep("tis1", length(tis1_spots)),rep("tis2", length(tis2_spots)), rep("edge", length(included_spots))),
                 gene1=c(averaged_expr_all[gene1,tis1_spots], averaged_expr_all[gene1,tis2_spots], averaged_expr_all[gene1,included_spots]),
                 gene2=c(averaged_expr_all[gene2,tis1_spots], averaged_expr_all[gene2,tis2_spots], averaged_expr_all[gene2,included_spots]))


  df$gene1.scale <- with(df, (gene1-min(gene1, na.rm=T))/diff(range(gene1, na.rm=T)))
  df$gene2.scale <- with(df, (gene2-min(gene2, na.rm=T))/diff(range(gene2, na.rm=T)))
  df<-df[!is.na(df$gene1.scale) & !is.na(df$gene2.scale),]

  df_midpoint<-data.frame(x_coord= midpoints_x,
                          y_coord= midpoints_y,
                          comm_score=averaged_expr_all[gene1,included_spots_tis1]*averaged_expr_all[gene2,included_spots_tis2])

 p <- ggplot(data=df, aes(x=x_coord, y=y_coord))+geom_point(data=df_midpoint, size=1, aes(colour=comm_score))+scale_color_continuous(type = "viridis")+theme_classic()

 return(p)
}



##compute weighted average of a module
#modules vector with module labels assignments
weighted_mod<-function(modules, kwithin, mod_sel, averaged_expr_all, comp1="_tis1", comp2="_tis2"){
  mod<-names(modules)[which(modules==mod_sel)]
  weights<-kwithin[[which(unique(modules)==mod_sel)]]$kExt1[mod[grep(comp1, mod)]]
  wm1<-apply((averaged_expr_all[gsub(comp1, "", mod[grep(comp1, mod)]),]), 2, function(x){weighted.mean(x, weights)})

  weights<-kwithin[[which(unique(modules)==mod_sel)]]$kExt2[mod[grep(comp2, mod)]]
  wm2<-apply((averaged_expr_all[gsub(comp2, "", mod[grep(comp2, mod)]),]), 2, function(x){weighted.mean(x, weights)})

  return(list(wm1, wm2))
}

