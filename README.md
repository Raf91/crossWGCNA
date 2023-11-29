# Cross-tissue gene expression interactions from bulk, single cell and spatial transcriptomics with crossWGCNA

CrossWGCNA is an open-source R package for the unsupervised identification of inter-tissue interactions from bulk, single-cell and spatial transcriptomics data. 
We employed crossWGCNA for the analysis of bulk RNA-seq data from laser-capture microdissected stroma and epithelium, single-cell RNA-seq data and spatial transcriptomics, for studying the tumour-stroma cross-talk in breast cancer.

![alt text](https://github.com/Raf91/crossWGCNA/blob/main/Graphical_Abstract.png)

The crossWGCNA pipeline involves several key steps for analyzing subject-matched transcriptomic data from two tissues. Here is a summary of each step:

####  1. Data Pre-processing:
- Input: Subject-matched transcriptomic data for two tissues.
- Gene Filtering: Genes are filtered based on expression variance across samples to reduce memory requirements.
- Expression Matrices: Two expression matrices (genes x samples) are generated, one for each tissue. They are combined row-wise, with added labels to distinguish tissues.

#### 2. Adjacency Calculation:
- Spearman Correlations: Calculated between each pair of rows in the pre-processed data.
- Corrections: Two correction methods are applied to remove external factors: i) setting self-loops to zero, or ii) setting self-loops to zero and subtracting the average correlation within each tissue.
- Soft Thresholding: A soft thresholding parameter is used in the adjacency calculation to preserve co-expression information.

#### 3. Adjacency Calculation Methods:
- "Signed" or "Unsigned" methods are applied, with the "signed" method employed in this analysis.
- Co-expression similarity (sij) is calculated based on gene expression profiles.

#### 4. Degrees Calculation:
- Degree Vectors: Three types of degree vectors (kInt, kExt, kTot) are assembled for each tissue.
- kInt: Number of connections within the same tissue.
- kExt: Number of connections to genes in the other tissue.
- kTot: Total number of connections involving the gene (kInt + kExt).
- kRatio: Measure of a gene's relevance for inter-tissue communication (kExt/kInt).

#### 4. Topological Overlap:
- Topological Overlap Matrix (TOM): Obtained using the WGCNA package.
- TOM combines the adjacency of two genes with their adjacencies to shared first neighbors in the network.
- TOMij formula is used for calculation.
- Hierarchical clustering is performed on the TOM after setting self-loops and intra-tissue TOMs to zero.

#### 5. Clustering:
- Dynamic Tree Cut: A dynamic tree cut algorithm is applied using the cutreeDynamic function with a distance metric of 1-TOMsimilarity and the method set to "average."
- Modules are identified using the mergeCloseModules function with a cutHeight of 0.25.

### Built with 

R (>= 4.2)

### Installation and usage

You can install the released version of crossWGCNA with:

```
# install.packages("devtools")
devtools::github_install(Raf91/crossWGCNA)
```

For details on how to apply crossWGCNA to matched stroma-epithelium RNAseq bulk data or Spatial Transcriptomics data please refer to vignettes provided within the package ![here](https://github.com/Raf91/crossWGCNA/blob/main/vignettes/). 
These are fully computable and reproducible R.markdown files guiding the user through each step of the pipeline analysis. For further details please refer to manual ![here](https://github.com/Raf91/crossWGCNA/blob/main/man/crossWGCNA_manual.pdf)

## References

1. Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 559.

3. Shannon P et al. Cytoscape: a software environment for integrated models of biomolecular interaction networks.
Genome Research 2003 Nov

4. Hao Y et al. (2023). “Dictionary learning for integrative, multimodal and scalable single-cell analysis.” Nature Biotechnology. 

5. Hao Y et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048.

6. Stuart T et al. (2019). “Comprehensive Integration of Single-Cell Data.” Cell, 177, 1888-1902. doi:10.1016/j.cell.2019.05.031.
    
7. Butler A et al. (2018). “Integrating single-cell transcriptomic data across different conditions, technologies, and species.” Nature Biotechnology, 36, 411-420.

8. Satija R et al. (2015). “Spatial reconstruction of single-cell gene expression data.” Nature Biotechnology, 33, 495-502. 

