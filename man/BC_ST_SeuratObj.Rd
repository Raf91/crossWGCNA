\name{Breast Cancer ST Seurat Object}
\alias{Breast Cancer ST Seurat Object}

\title{ST Seurat Object}

\description{
    Pre-computed clusters from a Breast Cancer spatial transcriptomics dataset provided by Illumina (Visium platform).
    A common Seurat (version \code{4.1.3}) analysis pipeline has been applied on the \code{"h5"} file to finally obtain clusters and umaps. 
}

\details{
    See Seurat package for details.
    The object is available for download at zenodo web page: code{\url{https://zenodo.org/record/8268805}} . 
}

\examples{
    download.file(url="https://zenodo.org/record/8268805/files/BC.RData?download=1", destfile = "~/BC.RData")
    load("~/BC.RData")
}

