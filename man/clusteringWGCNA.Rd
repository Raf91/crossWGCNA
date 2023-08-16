\name{clusteringWGCNA}
\alias{clusteringWGCNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  clusteringWGCNA
}
\description{

}
\usage{
clusteringWGCNA(A, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
    Adjacency matrix. Must be a symmetric matrix. 
}
  \item{data}{
    Input data matrix, where rows are genes and columns are samples. Gene expression measures in the two tissues are combined by row.
}
  \item{comp1}{
    Suffix (string character) to identify genes of tissue 1. Default="_1".
}
  \item{comp2}{
    Suffix (string character) to identify genes of tissue 2. Default="_2".
}
  \item{TOM}{
    Boolean value. If TRUE, it computes the TOM similarity matrix. 
    See WGCNA for details
}
  \item{ds}{
%%     ~~Describe \code{ds} here~~
}
  \item{crossOnly}{
%%     ~~Describe \code{crossOnly} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(crossWGCNA)
library(WGCNA)
inboth <- intersect(rownames(stroma), rownames(epi))
stroma <- stroma[inboth,]
epi <- epi[inboth,]
rownames(stroma) <- paste(rownames(stroma), "S", sep="_")
colnames(stroma) <- paste("Sample", c(1:34), sep="_")
rownames(epi) <- paste(rownames(epi), "E", sep="_")
colnames(epi) <- paste("Sample", c(1:34), sep="_")
BC_merged <- rbind(stroma, epi)
Adj <- Adjacency(data=BC_merged, method="netdiff", Adj_type="signed", comp1="_S", comp2="_E")
cl <- clusteringWGCNA(A=Adj,data=BC_merged,ds=2,comp1="_S",comp2="_E")
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
