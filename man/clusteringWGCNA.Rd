\name{clusteringWGCNA}
\alias{clusteringWGCNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  clusteringWGCNA
}
\description{
  Function to perform WGCNA clustering based on an adjacency matrix.
}
\usage{
clusteringWGCNA(A, data, comp1, comp2, TOM, ds, crossOnly)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
    Adjacency matrix. Should be the output of the function \code{\link{Adjacency}}.
    Must be a symmetric matrix.
}
  \item{data}{
    Input data matrix, where rows are genes and columns are samples. Gene expression measures in the two tissues are combined by row.
}
  \item{comp1}{
    Suffix (string character) to identify genes of tissue 1. Default=\code{"_1"}.
}
  \item{comp2}{
    Suffix (string character) to identify genes of tissue 2. Default=\code{"_2"}.
}
  \item{TOM}{
    Boolean value. If TRUE compute the Topological Overlap Matrix prior to clustering. Default=\code{TRUE}.
}
  \item{ds}{
    Parameter influencing the final number of clusters. It sets the \code{deepSplit} parameter in the \code{cutreeDynamic} function from the WGCNA package.
}
  \item{crossOnly}{
    Boolean value. If \code{TRUE} remove intra-tissue adjacencies prior to clustering (set to 0). Default=\code{TRUE}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of connectivities output of \code{\link{degrees}} and of gene clusters
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Aurora Savina, Raffaele M. Iannuzzi
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Adjacency}}, \code{\link{degrees}}.
}
\examples{
See package vignette
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
