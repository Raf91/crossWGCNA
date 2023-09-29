\name{degrees_mod}
\alias{degrees_mod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  degrees_mod
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
degrees_mod(data, method = "netdiff", modules, Adj_type="signed", 
            cortype="spearman", pval="none", thr=0.05, 
            beta=6, comp1="_tis1", comp2="_tis2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Input data matrix, where rows are genes and columns are samples. Gene expression measures in the two tissues are combined by row.
}
  \item{method}{
    String character identifying the method to use to compute cross-tissue adjacencies. Can be \code{"self-loop"} or \code{"netdiff"}.
}
  \item{modules}{
    List of co-expression modules. Output of \code{\link{clusteringWGCNA}}
}
  \item{Adj_type}{
    String character identifying the type of adjacency to compute. Can be \code{"signed"}, \code{"unsigned"} or \code{"keep sign"}. Default=\code{"signed"}. 
    \code{"keep sign"} is equivalent to \code{"unsigned"} in absolute values, but with sign corresponding to the sign of correlation.
    See WGCNA for details.
}
  \item{cortype}{
    String character identifying the type of correlation. Can be \code{"pearson"}, \code{"spearman"} or \code{"bicor"}. Default=\code{"spearman"}. 
    See WGCNA for details.
}
  \item{pval}{
    String character identifying whether to use the correlations' p-value in edge weight.
    Possible values are: \code{"none"} (does not include p-values), \code{"threshold"} (discard edges with a p-value higher than a threshold),
    \code{"weight"} (multiply edges' adjacencies by 1-pvalue). Default=\code{"none"}.
    Can be used with \code{"pearson"} and \code{"spearman"} correlations.
}
  \item{thr}{
    Threshold for cutting edges based on p-values. Default=\code{0.05}.
}
  \item{beta}{
    Power for adjacency calculation. Default=\code{6}.
}
  \item{comp1}{
    Suffix (string character) to identify genes of tissue 1. Default=\code{"_tis1"}.
}
  \item{comp2}{
    Suffix (string character) to identify genes of tissue 2. Default=\code{"_tis2"}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
For each module returns a list of degrees' vectors:
  \item{kInt1}{Intra-tissue degree for tissue 1}
  \item{kInt2}{Intra-tissue degree for tissue 2}
  \item{kExt1}{Inter-tissue degree for tissue 1 genes with tissue 2}
  \item{kExt2}{Inter-tissue degree for tissue 2 genes with tissue 1}
  \item{kTot1}{Overall connectivity of genes in tissue 1 with genes in both tissue 1 and 2}
  \item{kTot2}{Overall connectivity of genes in tissue 2 with genes in both tissue 1 and 2}
}
\references{
  Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 559.
}
\author{
Aurora Savino, Raffaele M. Iannuzzi.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{degrees}},
  \code{\link{Adjacency}},
  \code{\link{clusteringWGCNA}},
  \code{\link{crossWGCNA}}
}
\examples{
See package vignette.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
