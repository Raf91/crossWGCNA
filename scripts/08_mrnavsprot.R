library(openxlsx)

meta <-
  read.xlsx(
    "/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Data/CPTAC/BRCA/Mertins/CPTAC_BC_SupplementaryTable01.xlsx",
    1
  )
rownames(meta) <- meta[, 2]

CPTAC_BRCA <-
  read.csv(
    "/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Data/CPTAC/BRCA/Mertins/CPTAC_BC_SupplementaryTable03.csv",
    sep = ";",
    dec = ","
  )
rownames(CPTAC_BRCA) <- CPTAC_BRCA[, 1]
##remove replicates
CPTAC_BRCA <- CPTAC_BRCA[, -c(93:95)]
#change colnames
colnames(CPTAC_BRCA)[13:92] <-
  gsub("[.]", "-", colnames(CPTAC_BRCA)[13:92])
colnames(CPTAC_BRCA)[13:92] <-
  gsub("TCGA", "", colnames(CPTAC_BRCA)[13:92])
colnames(CPTAC_BRCA)[13:92] <-
  gsub('([^-]*)-([^-]*)-.*', '\\1-\\2', colnames(CPTAC_BRCA)[13:92])

###TCGA
library(TCGAbiolinks)
load(
  "/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Lanciano/TCGA_BRCA_tumor.RData"
)
colnames(TCGA_tumor) <-
  gsub('([^-]*)-([^-]*)-([^-]*)-.*',
       '\\1-\\2-\\3',
       colnames(TCGA_tumor))
subtype <- TCGAquery_subtype(tumor = "brca")
colnames(TCGA_tumor) <- gsub("TCGA-", "", colnames(TCGA_tumor))

##tenuti i geni espressi ad alti livelli in almeno 50 pazienti
TCGA_tumor <- TCGA_tumor[which(rowSums(TCGA_tumor < 1) <= 1040), ]
TCGA_tumor <- log2(TCGA_tumor + 1)

inboth <- intersect(colnames(CPTAC_BRCA), colnames(TCGA_tumor))
##only CPTAC samples in TCGA ID sel
CPTAC_BRCA_sel <-
  cbind.data.frame(CPTAC_BRCA[, 11], CPTAC_BRCA[, inboth])
TCGA_sel <- TCGA_tumor[, inboth]

source("scripts/crossWGCNA_functions_netdiff.R")
CPTAC_BRCA_sel <-
  changenames(data = CPTAC_BRCA_sel[, -1],
              anno = cbind(rownames(CPTAC_BRCA_sel), CPTAC_BRCA_sel[, 1]))

CPTAC_BRCA_sel <-
  CPTAC_BRCA_sel[-which(rowSums(is.na(CPTAC_BRCA_sel)) > 0), ]

genes <-
  unique(intersect(rownames(CPTAC_BRCA_sel)[which(apply(CPTAC_BRCA_sel, 1, var) >=
                                                    quantile(apply(CPTAC_BRCA_sel, 1, var), 0.25, na.rm = T))],
                   rownames(TCGA_sel)[which(apply(TCGA_sel, 1, var) >= quantile(apply(TCGA_sel, 1, var), 0.25, na.rm =
                                                                                  T))]))

CPTAC_BRCA_sel <- CPTAC_BRCA_sel[genes, ]
TCGA_sel <- TCGA_sel[genes, ]

rownames(TCGA_sel) <- paste(rownames(TCGA_sel), "mRNA", sep = "_")
rownames(CPTAC_BRCA_sel) <-
  paste(rownames(CPTAC_BRCA_sel), "prot", sep = "_")
data_merged <- rbind(TCGA_sel, CPTAC_BRCA_sel)

net <-
  network(
    data = data_merged,
    Adj_type = "signed",
    cortype = "pearson",
    pval = "none",
    thr = 0.05,
    beta = 6,
    comp1 = "_mRNA",
    comp2 = "_prot"
  )
save(net, file = "results/net_mRNAvsprot.RData")


library(msigdbr)
library(fgsea)

m_df = msigdbr(species = "Homo sapiens", category = "C5")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

coef_gsea <- net$kExt1 / net$kInt1
names(coef_gsea) <- gsub("_mRNA", "", names(coef_gsea))

fgseaRes <- fgseaMultilevel(m_list, coef_gsea)
write.xlsx(data.frame(fgseaRes), file = "results/GSEA_mRNA_ratio.xlsx")



coef_gsea <- net$kExt2 / net$kInt2
names(coef_gsea) <- gsub("_prot", "", names(coef_gsea))

fgseaRes <- fgseaMultilevel(m_list, coef_gsea)
write.xlsx(data.frame(fgseaRes), file = "results/GSEA_prot_ratio.xlsx")
