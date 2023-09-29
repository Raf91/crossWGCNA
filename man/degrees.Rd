\name{degrees}
\alias{degrees}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  degrees
}
\description{
  Function to compute the weighted degree of a gene with all genes of the same or of the alternate tissue/organ.
}
\usage{
degrees(A, comp1="_tis1", comp2="_tis2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
    Matrix of adjacencies, output of the function \code{Adjacency}.
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
A list of degrees' vectors:
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{kInt1}{Intra-tissue degree for tissue 1}
  \item{kInt2}{Intra-tissue degree for tissue 2}
  \item{kExt1}{Inter-tissue degree for tissue 1 genes with tissue 2}
  \item{kExt2}{Inter-tissue degree for tissue 2 genes with tissue 1}
  \item{kTot1}{Overall connectivity of genes in tissue 1 with genes in both tissue 1 and 2}
  \item{kTot2}{Overall connectivity of genes in tissue 2 with genes in both tissue 1 and 2}
}
\references{
  Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation network analysis.Returns a list of connectivities  BMC Bioinformatics, 559.
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
\code{\link{clusteringWGCNA}},
\code{\link{crossWGCNA}},
\code{\link{degrees_mod}}.
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
