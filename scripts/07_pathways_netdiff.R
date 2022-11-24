rm(list = ls())
source("scripts/crossWGCNA_functions_netdiff.R")
####################################
library(openxlsx)
cellcycle <- read.xlsx("data/Pathways.xlsx", 1)[-1, 1]
hippo <- read.xlsx("data/Pathways.xlsx", 2)[-1, 1]
myc <- read.xlsx("data/Pathways.xlsx", 3)[-1, 1]
notch <- read.xlsx("data/Pathways.xlsx", 4)[-1, 1]
nrf2 <- read.xlsx("data/Pathways.xlsx", 5)[-1, 1]
pi3k <- read.xlsx("data/Pathways.xlsx", 6)[-1, 1]
tgfb <- read.xlsx("data/Pathways.xlsx", 7)[-1, 1]
rtkras <- read.xlsx("data/Pathways.xlsx", 8)[-1, 1]
p53 <- read.xlsx("data/Pathways.xlsx", 9)[-1, 1]
wnt <- read.xlsx("data/Pathways.xlsx", 10)[-1, 1]



library(openxlsx)
sign_cellcycle <- read.xlsx("data/Pathways.xlsx", 1)[-1, 3]
sign_hippo <- read.xlsx("data/Pathways.xlsx", 2)[-1, 3]
sign_myc <- read.xlsx("data/Pathways.xlsx", 3)[-1, 3]
sign_notch <- read.xlsx("data/Pathways.xlsx", 4)[-1, 3]
sign_nrf2 <- read.xlsx("data/Pathways.xlsx", 5)[-1, 3]
sign_pi3k <- read.xlsx("data/Pathways.xlsx", 6)[-1, 3]
sign_tgfb <- read.xlsx("data/Pathways.xlsx", 7)[-1, 3]
sign_rtkras <- read.xlsx("data/Pathways.xlsx", 8)[-1, 3]
sign_p53 <- read.xlsx("data/Pathways.xlsx", 9)[-1, 3]
sign_wnt <- read.xlsx("data/Pathways.xlsx", 10)[-1, 3]

sign_cellcycle[which(sign_cellcycle == "OG")] <- 1
sign_cellcycle[which(sign_cellcycle == "TSG")] <- (-1)
sign_cellcycle[-which(sign_cellcycle %in% c("1", "-1"))] <- NA
cellcycle <- cellcycle[which(sign_cellcycle %in% c("1", "-1"))]
sign_cellcycle <-
  sign_cellcycle[which(sign_cellcycle %in% c("1", "-1"))]



sign_hippo[which(sign_hippo == "OG")] <- 1
sign_hippo[which(sign_hippo == "TSG")] <- (-1)
sign_hippo[-which(sign_hippo %in% c("1", "-1"))] <- NA
hippo <- hippo[which(sign_hippo %in% c("1", "-1"))]
sign_hippo <- sign_hippo[which(sign_hippo %in% c("1", "-1"))]


sign_myc[which(sign_myc == "OG")] <- 1
sign_myc[which(sign_myc == "TSG")] <- (-1)
sign_myc[-which(sign_myc %in% c("1", "-1"))] <- NA
myc <- myc[which(sign_myc %in% c("1", "-1"))]
sign_myc <- sign_myc[which(sign_myc %in% c("1", "-1"))]


sign_notch[which(sign_notch == "OG")] <- 1
sign_notch[which(sign_notch == "TSG")] <- (-1)
sign_notch[-which(sign_notch %in% c("1", "-1"))] <- NA
notch <- notch[which(sign_notch %in% c("1", "-1"))]
sign_notch <- sign_notch[which(sign_notch %in% c("1", "-1"))]


sign_nrf2[which(sign_nrf2 == "OG")] <- 1
sign_nrf2[which(sign_nrf2 == "TSG")] <- (-1)
sign_nrf2[-which(sign_nrf2 %in% c("1", "-1"))] <- NA
nrf2 <- nrf2[which(sign_nrf2 %in% c("1", "-1"))]
sign_nrf2 <- sign_nrf2[which(sign_nrf2 %in% c("1", "-1"))]


sign_pi3k[which(sign_pi3k == "OG")] <- 1
sign_pi3k[which(sign_pi3k == "TSG")] <- (-1)
sign_pi3k[-which(sign_pi3k %in% c("1", "-1"))] <- NA
pi3k <- pi3k[which(sign_pi3k %in% c("1", "-1"))]
sign_pi3k <- sign_pi3k[which(sign_pi3k %in% c("1", "-1"))]


sign_tgfb[which(sign_tgfb == "OG")] <- 1
sign_tgfb[which(sign_tgfb == "TSG")] <- (-1)
sign_tgfb[-which(sign_tgfb %in% c("1", "-1"))] <- NA
tgfb <- tgfb[which(sign_tgfb %in% c("1", "-1"))]
sign_tgfb <- sign_tgfb[which(sign_tgfb %in% c("1", "-1"))]

sign_rtkras[which(sign_rtkras == "OG")] <- 1
sign_rtkras[which(sign_rtkras == "TSG")] <- (-1)
sign_rtkras[-which(sign_rtkras %in% c("1", "-1"))] <- NA
rtkras <- rtkras[which(sign_rtkras %in% c("1", "-1"))]
sign_rtkras <- sign_rtkras[which(sign_rtkras %in% c("1", "-1"))]


sign_p53[which(sign_p53 == "OG")] <- 1
sign_p53[which(sign_p53 == "TSG")] <- (-1)
sign_p53[-which(sign_p53 %in% c("1", "-1"))] <- NA
p53 <- p53[which(sign_p53 %in% c("1", "-1"))]
sign_p53 <- sign_p53[which(sign_p53 %in% c("1", "-1"))]


sign_wnt[which(sign_wnt == "OG")] <- 1
sign_wnt[which(sign_wnt == "TSG")] <- (-1)
sign_wnt[-which(sign_wnt %in% c("1", "-1"))] <- NA
wnt <- wnt[which(sign_wnt %in% c("1", "-1"))]
sign_wnt <- sign_wnt[which(sign_wnt %in% c("1", "-1"))]



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


path <-
  c("cellcycle",
    "hippo",
    "myc",
    "notch",
    "nrf2",
    "pi3k",
    "tgfb",
    "rtkras",
    "p53",
    "wnt")
for (i in 7:10) {
  selgenes <- get(path[i])
  sign <- as.numeric(get(paste("sign_", path[i], sep = "")))
  names(sign) <- selgenes

  #Datasets
  ###dataset GSE5847

  load("data/metadataAll.RData")
  load("data/datasetsAll.RData")


  #Datasets
  ###dataset GSE5847

  stromaID <-
    metadataAll$id[which(
      metadataAll$dataset == "GSE5847" &
        metadataAll$compartment == "Stroma" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]
  epiID <-
    metadataAll$id[which(
      metadataAll$dataset == "GSE5847" &
        metadataAll$compartment == "Epi" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]

  stromaMatch <-
    metadataAll$matching[which(
      metadataAll$dataset == "GSE5847" &
        metadataAll$compartment == "Stroma" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]
  epiMatch <-
    metadataAll$matching[which(
      metadataAll$dataset == "GSE5847" &
        metadataAll$compartment == "Epi" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]

  inboth <- intersect(stromaMatch, epiMatch)

  stromaID <- stromaID[match(inboth, stromaMatch)]
  epiID <- epiID[match(inboth, epiMatch)]

  stroma <- datasetsAll[["GSE5847"]][, stromaID]
  epi <- datasetsAll[["GSE5847"]][, epiID]

  ###filtering data
  genes <-
    unique(c(which(
      apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
    ),
    which(
      apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
    )))
  stroma <- stroma[genes, ]
  epi <- epi[genes, ]

  #merging stroma and epi
  rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
  rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")
  colnames(epi) <- colnames(stroma)

  data_merged_GSE5847 <- rbind(stroma, epi)


  ############ GSE10797
  stromaID <-
    metadataAll$id[which(
      metadataAll$dataset == "GSE10797" &
        metadataAll$compartment == "Stroma" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]
  epiID <-
    metadataAll$id[which(
      metadataAll$dataset == "GSE10797" &
        metadataAll$compartment == "Epi" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]

  stromaMatch <-
    metadataAll$matching[which(
      metadataAll$dataset == "GSE10797" &
        metadataAll$compartment == "Stroma" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]
  epiMatch <-
    metadataAll$matching[which(
      metadataAll$dataset == "GSE10797" &
        metadataAll$compartment == "Epi" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]

  inboth <- intersect(stromaMatch, epiMatch)

  stromaID <- stromaID[match(inboth, stromaMatch)]
  epiID <- epiID[match(inboth, epiMatch)]

  stroma <- datasetsAll[["GSE10797"]][, stromaID]
  epi <- datasetsAll[["GSE10797"]][, epiID]

  ###filtering data
  genes <-
    unique(c(which(
      apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
    ),
    which(
      apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
    )))
  stroma <- stroma[genes, ]
  epi <- epi[genes, ]
  ##############All genes double WGCNA
  ###use WGCNA with twice the genes
  rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
  rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")

  colnames(epi) <- colnames(stroma)
  data_merged_GSE10797 <- rbind(stroma, epi)

  ##############GSE14548
  stromaID <-
    metadataAll$id[which(
      metadataAll$dataset == "GSE14548" &
        metadataAll$compartment == "Stroma" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]
  epiID <-
    metadataAll$id[which(
      metadataAll$dataset == "GSE14548" &
        metadataAll$compartment == "Epi" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]

  stromaMatch <-
    metadataAll$matching[which(
      metadataAll$dataset == "GSE14548" &
        metadataAll$compartment == "Stroma" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]
  epiMatch <-
    metadataAll$matching[which(
      metadataAll$dataset == "GSE14548" &
        metadataAll$compartment == "Epi" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]

  inboth <- intersect(stromaMatch, epiMatch)

  stromaID <- stromaID[match(inboth, stromaMatch)]
  epiID <- epiID[match(inboth, epiMatch)]

  stroma <- datasetsAll[["GSE14548"]][, stromaID]
  epi <- datasetsAll[["GSE14548"]][, epiID]

  ###filtering data
  genes <-
    unique(c(which(
      apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
    ),
    which(
      apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
    )))
  stroma <- stroma[genes, ]
  epi <- epi[genes, ]

  rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
  rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")

  colnames(epi) <- colnames(stroma)
  data_merged_GSE14548 <- rbind(stroma, epi)


  #############GSE83591
  load("data/GSE83591.RData")
  GSE83591_meta <- read.csv("data/GSE83591_metadata.txt", sep = "\t")
  GSE83591_meta <- t(GSE83591_meta)
  colnames(GSE83591_meta) <- GSE83591_meta[1, ]
  GSE83591_meta <- GSE83591_meta[-1, ]

  GSE83591_meta2 <- read.xlsx("data/GSE83591_Correspondence_LCM.xlsx", 1)
  GSE83591_meta2$ID <-
    unlist(strsplit(GSE83591_meta2$GSE83591, "_"))[seq(1, 109 * 9, 9)]

  #reorder based on annotation file
  GSE83591 <- GSE83591[, match(GSE83591_meta2$ID, colnames(GSE83591))]
  GSE83591_meta <-
    GSE83591_meta[match(GSE83591_meta2$ID, GSE83591_meta[, 1]), ]

  GSE83591 <- GSE83591[, which(GSE83591_meta2$Cy3 %in% c("TE", "TS"))]
  GSE83591_meta <-
    GSE83591_meta[which(GSE83591_meta2$Cy3 %in% c("TE", "TS")), ]
  GSE83591_meta2 <-
    GSE83591_meta2[which(GSE83591_meta2$Cy3 %in% c("TE", "TS")), ]

  GSE83591 <- GSE83591[, -c(15, 34, 67)]
  GSE83591_meta <- GSE83591_meta[-c(15, 34, 67), ]
  GSE83591_meta2 <- GSE83591_meta2[-c(15, 34, 67), ]

  stroma <- GSE83591[, which(GSE83591_meta2$Cy3 == "TS")]
  epi <- GSE83591[, which(GSE83591_meta2$Cy3 == "TE")]

  ###filtering data
  genes <-
    unique(c(which(
      apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
    ),
    which(
      apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
    )))
  stroma <- stroma[genes, ]
  epi <- epi[genes, ]

  rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
  rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")

  colnames(epi) <- colnames(stroma)
  data_merged_GSE83591 <- rbind(stroma, epi)

  #######GSE68744
  stromaID <-
    metadataAll$id[which(
      metadataAll$dataset == "GSE68744" &
        metadataAll$compartment == "Stroma" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]
  epiID <-
    metadataAll$id[which(
      metadataAll$dataset == "GSE68744" &
        metadataAll$compartment == "Epi" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]
  #sono matched, controllato

  stroma <- datasetsAll[["GSE68744"]][, stromaID]
  epi <- datasetsAll[["GSE68744"]][, epiID]

  ###filtering data
  genes <-
    unique(c(which(
      apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
    ),
    which(
      apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
    )))
  stroma <- stroma[genes, ]
  epi <- epi[genes, ]

  rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
  rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")

  colnames(epi) <- colnames(stroma)
  data_merged_GSE68744 <- rbind(stroma, epi)

  #######GSE88715
  stromaID <-
    metadataAll$id[which(
      metadataAll$dataset == "GSE88715" &
        metadataAll$compartment == "Stroma" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]
  epiID <-
    metadataAll$id[which(
      metadataAll$dataset == "GSE88715" &
        metadataAll$compartment == "Epi" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]

  stromaMatch <-
    metadataAll$matching[which(
      metadataAll$dataset == "GSE88715" &
        metadataAll$compartment == "Stroma" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]
  epiMatch <-
    metadataAll$matching[which(
      metadataAll$dataset == "GSE88715" &
        metadataAll$compartment == "Epi" &
        metadataAll$diseaseStatus == "InvasiveBC"
    )]

  inboth <- intersect(stromaMatch, epiMatch)

  stromaID <- stromaID[match(inboth, stromaMatch)]
  epiID <- epiID[match(inboth, epiMatch)]

  stroma <- datasetsAll[["GSE88715"]][, stromaID]
  epi <- datasetsAll[["GSE88715"]][, epiID]

  ###filtering data
  genes <-
    unique(c(which(
      apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
    ),
    which(
      apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
    )))
  stroma <- stroma[genes, ]
  epi <- epi[genes, ]

  rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
  rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")


  colnames(epi) <- colnames(stroma)
  data_merged_GSE88715 <- rbind(stroma, epi)

  #ext sono NA, ma non se faccio girare le funzioni da dentro. Il problema è in Adj
  degs_GSE5847 <-
    network(
      selgenes = selgenes,
      data = data_merged_GSE5847,
      Adj_type = "keep sign",
      cortype = "pearson",
      pval = "none",
      thr = 0.05,
      beta = 6,
      comp1 = "_tis1$",
      comp2 = "_tis2$",
      sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE5847)[grep("_tis2$", rownames(data_merged_GSE5847))])],
      compartment_sel = "comp2"
    )
  degs_GSE10797 <-
    network(
      selgenes = selgenes,
      data = data_merged_GSE10797,
      Adj_type = "keep sign",
      cortype = "pearson",
      pval = "none",
      thr = 0.05,
      beta = 6,
      comp1 = "_tis1$",
      comp2 = "_tis2$",
      sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE10797)[grep("_tis2$", rownames(data_merged_GSE10797))])],
      compartment_sel = "comp2"
    )
  degs_GSE14548 <-
    network(
      selgenes = selgenes,
      data = data_merged_GSE14548,
      Adj_type = "keep sign",
      cortype = "pearson",
      pval = "none",
      thr = 0.05,
      beta = 6,
      comp1 = "_tis1$",
      comp2 = "_tis2$",
      sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE14548)[grep("_tis2$", rownames(data_merged_GSE14548))])],
      compartment_sel = "comp2"
    )
  degs_GSE83591 <-
    network(
      selgenes = selgenes,
      data = data_merged_GSE83591,
      Adj_type = "keep sign",
      cortype = "pearson",
      pval = "none",
      thr = 0.05,
      beta = 6,
      comp1 = "_tis1$",
      comp2 = "_tis2$",
      sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE83591)[grep("_tis2$", rownames(data_merged_GSE83591))])],
      compartment_sel = "comp2"
    )
  degs_GSE68744 <-
    network(
      selgenes = selgenes,
      data = data_merged_GSE68744,
      Adj_type = "keep sign",
      cortype = "pearson",
      pval = "none",
      thr = 0.05,
      beta = 6,
      comp1 = "_tis1$",
      comp2 = "_tis2$",
      sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE68744)[grep("_tis2$", rownames(data_merged_GSE68744))])],
      compartment_sel = "comp2"
    )
  degs_GSE88715 <-
    network(
      selgenes = selgenes,
      data = data_merged_GSE88715,
      Adj_type = "keep sign",
      cortype = "pearson",
      pval = "none",
      thr = 0.05,
      beta = 6,
      comp1 = "_tis1$",
      comp2 = "_tis2$",
      sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE88715)[grep("_tis2$", rownames(data_merged_GSE88715))])],
      compartment_sel = "comp2"
    )

  degs_path_sign <- list(
    degs_GSE5847,
    degs_GSE10797,
    degs_GSE14548,
    degs_GSE83591,
    degs_GSE68744,
    degs_GSE88715
  )

  save(degs_path_sign,
       file = paste("results/degs_", path[i], "_KS.RData", sep = ""))
}

load("~/crossWGCNA/results/degs_cellcycle_KS.RData")
incommon <- Reduce(intersect, list(
  names(degs_path_sign[[1]]$kExt1),
  names(degs_path_sign[[2]]$kExt1),
  names(degs_path_sign[[3]]$kExt1),
  names(degs_path_sign[[4]]$kExt1),
  names(degs_path_sign[[5]]$kExt1),
  names(degs_path_sign[[6]]$kExt1)
))


library(biomaRt)
ensembl.human <- useEnsembl(biomart = 'genes',
                            dataset = 'hsapiens_gene_ensembl')

#extracellular region
GO <- getBM(
  attributes = c('hgnc_symbol'),
  filters = 'go',
  values = "GO:0005576",
  mart = ensembl.human
)

extracellular <- intersect(incommon, paste(GO[, 1], "_tis1", sep = ""))

i <- 1
load(file = paste("results/degs_", path[i], "_KS.RData", sep = ""))
all_sign <-
  (rank(-((degs_path_sign[[1]]$kExt1 / degs_path_sign[[1]]$kInt1)[extracellular]
  )) *
    rank(-((degs_path_sign[[2]]$kExt1 / degs_path_sign[[2]]$kInt1)[extracellular]
    )) *
    rank(-((degs_path_sign[[3]]$kExt1 / degs_path_sign[[3]]$kInt1)[extracellular]
    )) *
    rank(-((degs_path_sign[[4]]$kExt1 / degs_path_sign[[4]]$kInt1)[extracellular]
    )) *
    rank(-((degs_path_sign[[5]]$kExt1 / degs_path_sign[[5]]$kInt1)[extracellular]
    )) *
    rank(-((degs_path_sign[[6]]$kExt1 / degs_path_sign[[6]]$kInt1)[extracellular]
    )))
for (i in 2:10) {
  load(file = paste("results/degs_", path[i], "_KS.RData", sep = ""))
  rankprod <-
    (rank(-((degs_path_sign[[1]]$kExt1 / degs_path_sign[[1]]$kInt1)[extracellular]
    )) *
      rank(-((degs_path_sign[[2]]$kExt1 / degs_path_sign[[2]]$kInt1)[extracellular]
      )) *
      rank(-((degs_path_sign[[3]]$kExt1 / degs_path_sign[[3]]$kInt1)[extracellular]
      )) *
      rank(-((degs_path_sign[[4]]$kExt1 / degs_path_sign[[4]]$kInt1)[extracellular]
      )) *
      rank(-((degs_path_sign[[5]]$kExt1 / degs_path_sign[[5]]$kInt1)[extracellular]
      )) *
      rank(-((degs_path_sign[[6]]$kExt1 / degs_path_sign[[6]]$kInt1)[extracellular]
      )))

  all_sign <- cbind(all_sign, rankprod)

}


colnames(all_sign) <- path
rownames(all_sign) <- gsub("_tis1", "", rownames(all_sign))


#all contiene tutti i rankproduct nell'ordine dei geni extracellular
save(all_sign, file = "results/allRankprodExt_KS.RData")

genes <-
  which(
    rank(all_sign[, 1]) <= 5 |
      rank(all_sign[, 2]) <= 5 |
      rank(all_sign[, 3]) <= 5 |
      rank(all_sign[, 4]) <= 5 |
      rank(all_sign[, 5]) <= 5 |
      rank(all_sign[, 6]) <= 5 |
      rank(all_sign[, 7]) <= 5 |
      rank(all_sign[, 8]) <= 5 |
      rank(all_sign[, 9]) <= 5 | rank(all_sign[, 10]) <= 5
  )

png("Pathways_topInteractorsKS.png", res = 300, 3000, 1500)
pheatmap(t(log10(all_sign[genes, ])), cellwidth = 15, cellheight = 15)
dev.off()

tabella <-
  cbind(
    names(sort(rank(all_sign[, 1]))[1:5]),
    names(sort(rank(all_sign[, 2]))[1:5]),
    names(sort(rank(all_sign[, 3]))[1:5]),
    names(sort(rank(all_sign[, 4]))[1:5]),
    names(sort(rank(all_sign[, 5]))[1:5]),
    names(sort(rank(all_sign[, 6]))[1:5]),
    names(sort(rank(all_sign[, 7]))[1:5]),
    names(sort(rank(all_sign[, 8]))[1:5]),
    names(sort(rank(all_sign[, 9]))[1:5]),
    names(sort(rank(all_sign[, 10]))[1:5])
  )

colnames(tabella) <- path

write.xlsx(data.frame(tabella), "top5paths.xlsx")
