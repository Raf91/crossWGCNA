\name{Adjacency}
\alias{Adjacency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Adjacency
}
\description{
  Wrapper function to compute cross-tissue gene expression adjacencies. The adjacencies are computed accordingly to the method selected (i.e., "self-loop" or "netdiff").
}
\usage{
  Adjacency(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Input data matrix, where rows are genes and columns are samples. Gene expression measures in the two tissues are combined by row.
}
\item{method}{
Character string identifying the method to use to compute cross-tissue adjacencies. Can be "self-loop" or "netdiff".
}
\item{comp1}{
Suffix to identify genes of tissue 1. Default="_1".
}
\item{comp2}{
Suffix to identify genes of tissue 2. Default="_2".
}
\item{Adj_type}{
Type of adjacency. Can be "signed", "unsigned" or "keep sign". Default="signed". See the package WGCNA for details. "keep sign" is equivalent to "unsigned" in absolute values, but with sign corresponding to the sign of correlation.
}
\item{cortype}{
Type of correlation. Can be "pearson", "spearman" or "bicor". Default="spearman". See the package WGCNA for details.
}
\item{pval}{
Whether to use the correlations' p-value in edge weight.
Possible values are: "none"" (does not include p-values), "threshold"" (discard edges with a p-value higher than a threshold),
"weight" (multiply edges' adjacencies by 1-pvalue). Default="none".
Can be used with "pearson" and "spearman" correlations.
}
\item{thr}{
Threshold for cutting edges based on p-values. Default="0.05".
}
\item{beta}{
Power for adjacency calculation. Default=6.
}
\item{sign_list}{
When considering genes associated with a pathway with sign (activation: +1, inactivation: -1), a vector with a value for each gene in the pathway. Default=1, meaning that each gene has the same positive weight.
}
\item{compartment_sel}{
In which tissue the signed pathway's expression is considered. Can be "none", "comp1" or "comp2".
}
\item{selgenes}{
Vector of genes associated with the pathway of interest.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a matrix with pairwise gene adjacencies.
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
See package vignette
}

