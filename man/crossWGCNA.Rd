\name{crossWGCNA}
\alias{crossWGCNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  crossWGCNA
}
\description{
  One-step function to obtain connectivities of each gene with all genes of the same or of the alternate tissue/organ, 
  and gene clusters based on inter-tissue connectivities.
}
\usage{
crossWGCNA(data,method="netdiff",Adj_type="signed",cortype="spearman",
           pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2", 
           doClusters=TRUE, doTOM=TRUE, ds=1, crossOnly=TRUE,verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Input data matrix, where rows are genes and columns are samples. Gene expression measures in the two tissues are combined by row.
}
  \item{method}{
    String character identifying the method to use to compute cross-tissue adjacencies. Can be \code{"self-loop"} or \code{"netdiff"}.
}
  \item{Adj_type}{
    String character identifying the type of adjacency to compute. Can be \code{"signed"}, \code{"unsigned"} or \code{"keep sign"}. Default=\code{"signed"}. 
    \code{"keep sign"} is equivalent to \code{"unsigned"} in absolute values, but with sign corresponding to the sign of correlation.
    See WGCNA package for details.
}
  \item{cortype}{
    String character identifying the type of correlation. Can be \code{"pearson"}, \code{"spearman"} or \code{"bicor"}. Default=\code{"spearman"}. 
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
    Suffix (string character) to identify genes of tissue 1. Default=\code{"_1"}.
}
  \item{comp2}{
    Suffix (string character) to identify genes of tissue 2. Default=\code{"_2"}.
}
  \item{doClusters}{
    Boolean value. If \code{TRUE} compute clusters. Default=\code{TRUE}.
}
  \item{doTOM}{
    Boolean value. If \code{TRUE} compute the Topological Overlap Matrix prior to clustering. Default=\code{TRUE}.
}
  \item{ds}{
    Parameter influencing the final number of clusters. It sets the \code{deepSplit} parameter in the \code{cutreeDynamic} function from the WGCNA package.
}
  \item{crossOnly}{
    Boolean value. If \code{TRUE} remove intra-tissue adjacencies prior to clustering (set to 0). Default=\code{TRUE}.
}
  \item{verbose}{
    Boolean value. If \code{TRUE} prints logging output. Default=\code{TRUE}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If \code{doClusters} is set to \code{FALSE}, returns a list of connectivities otherwise returns both list of connectivities and gene modules.
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
\code{\link{Adjacency}}, 
\code{\link{rm_selfloop}},
\code{\link{rm_netdiff}},
\code{\link{degrees}},
\code{\link{clusteringWGCNA}}.
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
