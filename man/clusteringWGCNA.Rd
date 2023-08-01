\name{clusteringWGCNA}
\alias{clusteringWGCNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  clusteringWGCNA
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clusteringWGCNA(A, data, comp1 = "_1", comp2 = "_2", TOM = T, ds = 1, crossOnly = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{comp1}{
%%     ~~Describe \code{comp1} here~~
}
  \item{comp2}{
%%     ~~Describe \code{comp2} here~~
}
  \item{TOM}{
%%     ~~Describe \code{TOM} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (A, data, comp1 = "_1", comp2 = "_2", TOM = T, ds = 1, 
    crossOnly = T) 
{
    comp1 <- paste(comp1, "$", sep = "")
    comp2 <- paste(comp2, "$", sep = "")
    if (crossOnly == T) {
        A[grep(comp1, rownames(A)), grep(comp1, rownames(A))] <- 0
        A[grep(comp2, rownames(A)), grep(comp2, rownames(A))] <- 0
    }
    genes1 <- gsub(comp1, "", rownames(A)[grep(comp1, rownames(A))])
    genes2 <- gsub(comp2, "", rownames(A)[grep(comp2, rownames(A))])
    Idx1 <- cbind(grep(comp1, rownames(A)), grep(comp2, rownames(A))[match(genes1, 
        genes2)])
    A[Idx1] <- 0
    Idx2 <- cbind(grep(comp2, rownames(A)), grep(comp1, rownames(A))[match(genes2, 
        genes1)])
    A[Idx2] <- 0
    if (TOM == T) {
        similarity <- TOMsimilarity(A, TOMType = "signed")
        rownames(similarity) <- rownames(A)
        colnames(similarity) <- colnames(A)
        A <- similarity
        rm(similarity)
    }
    if (crossOnly == T) {
        A[grep(comp1, rownames(A)), grep(comp1, rownames(A))] <- 0
        A[grep(comp2, rownames(A)), grep(comp2, rownames(A))] <- 0
    }
    conTree = hclust(as.dist(1 - A), method = "average")
    unmergedLabels = cutreeDynamic(dendro = conTree, dist = 1 - 
        A, deepSplit = ds)
    names(unmergedLabels) <- rownames(A)
    merged = mergeCloseModules(t(data), unmergedLabels, cutHeight = 0.25, 
        verbose = 3)
    names(merged$colors) <- rownames(A)
    return(merged)
  }
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
