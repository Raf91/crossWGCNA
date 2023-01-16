library(pheatmap)
library(ggpubr)
save_pheatmap_pdf <- function(x, filename, width=5, height=5) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}


load("results/degs_3rd_netdiff.RData")
degs_netdiff<-degs
load("results/degs_3rd_selfloops.RData")


names(degs) <- c("GSE5847",
                 "GSE10797",
                 "GSE14548",
                 "GSE83591",
                 "GSE68744",
                 "GSE88715")
names(degs_netdiff) <- c("GSE5847",
                 "GSE10797",
                 "GSE14548",
                 "GSE83591",
                 "GSE68744",
                 "GSE88715")

k_type<-c("kInt_stroma", "kInt_epi", "kExt_stroma", "kExt_epi", "kTot_stroma", "kTot_epi")

for(k in 1:6){

corrs <- matrix(nrow = length(degs_netdiff), ncol = length(degs_netdiff))
for (i in 1:length(degs_netdiff)) {
    for (j in 1:length(degs_netdiff)) {
        incommon <-
            intersect(names(degs_netdiff[[i]][[k]]),
                      names(degs_netdiff[[j]][[k]]))
        corrs[i, j] <- cor(degs_netdiff[[i]][[k]][incommon], degs_netdiff[[j]][[k]][incommon])
    }
}


corrs_sl <- matrix(nrow = length(degs), ncol = length(degs))
for (i in 1:length(degs)) {
    for (j in 1:length(degs)) {
        incommon <-
            intersect(names(degs[[i]][[k]]),
                      names(degs[[j]][[k]]))
        corrs_sl[i, j] <-
            cor(degs[[i]][[k]][incommon], degs[[j]][[k]][incommon])
    }
}

colnames(corrs_sl)<-names(degs)
rownames(corrs_sl)<-names(degs)
colnames(corrs)<-names(degs)
rownames(corrs)<-names(degs)

df<-data.frame(R=c(corrs_sl[upper.tri(corrs_sl)],corrs[upper.tri(corrs)]),
               method=c(rep("selfloops", length=length(corrs_sl[upper.tri(corrs_sl)])),
                                                         rep("netdiff", length=length(corrs[upper.tri(corrs)]))))

pdf(paste("results/Compare_corr_",k_type[k],".pdf",sep=""), 4, 6)
my_comparisons<-list(c("selfloops", "netdiff"))
p<-ggviolin(df, x = "method", y = "R",
            palette = "jco", add="boxplot") + stat_compare_means(comparisons=my_comparisons) +
    rotate_x_text(angle = 45)
print(p)
dev.off()


paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)
myBreaks <- c(seq(0,max(unlist(corrs_sl), unlist(corrs), na.rm=T),
              length.out=floor(paletteLength)))
length(myBreaks) == length(paletteLength) + 1

p<-pheatmap(corrs_sl,   cellwidth=10, cellheight=10, breaks=myBreaks, color = myColor)
save_pheatmap_pdf(p, paste("results/corr_",k_type[k],"_selfloops.pdf", sep=""))

p<-pheatmap(corrs,   cellwidth=10, cellheight=10, breaks=myBreaks, color = myColor)
save_pheatmap_pdf(p, paste("results/corr_",k_type[k],"_netdiff.pdf", sep=""))
}
