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

### Usage

For details on how to apply crossWGCNA to matched stroma-epithelium RNAseq bulk data or Spatial Transcriptomics data please refer to vignettes provided within the package ![here](https://github.com/Raf91/crossWGCNA/blob/main/vignettes/). These are fully computable R.markdown guiding the user through each step of the pipeline analyses. 
