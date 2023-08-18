\name{cytoscape_net}
\alias{cytoscape_net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  cytoscape_net
}
\description{
Finds the top cross-tissue interactors of a gene and saves them as a Cytoscape input.
}
\usage{
cytoscape_net(A, data, gene, comp1, comp2, num, corr="spearman")
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
  \item{gene}{
    Gene to be analyzed (string character).
}
  \item{comp1}{
    Suffix (string character) to identify genes of tissue 1. Default=\code{"_1"}.
}
  \item{comp2}{
    Suffix (string character) to identify genes of tissue 2. Default=\code{"_2"}.
}
  \item{num}{
    Number (numeric) of first closest neigbours of the selected genes to consider.
}
  \item{corr}{
    String character identifying the type of correlation. Can be \code{"pearson"}, \code{"spearman"} or \code{"bicor"}. Default=\code{"spearman"}. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Dataframe with the top interactors of the selected gene and their pairwise adjacencies within and between tissues. Columns of the dataframe are:
\item{Source}{
  source gene
}
\item{Target}{
  connected gene
}
\item{Weight}{
  connection strenght: adjacency
}
\item{Edge_type}{
  type of interaction: within compartment 1 (intra1), compartment 2 (intra2) or between compartments (inter)
}
\item{Source_type}{
  label indicating whether the source gene has been measured in compartment 1 or compartment 2
}
\item{Target_type}{
  label indicating whether the target gene has been measured in compartment 1 or compartment 2
}
}
\references{
  Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T.
  Cytoscape: a software environment for integrated models of biomolecular interaction networks.
  Genome Research 2003 Nov; 13(11):2498-504 
}
\author{
Aurora Savino, Raffaele M. Iannuzzi.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

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
