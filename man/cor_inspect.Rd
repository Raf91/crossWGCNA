\name{cor_inspect}
\alias{cor_inspect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  cor_inspect
}
\description{
Plot the correlation between two genes in the same or in different compartments/tissues.
}
\usage{
cor_inspect(data, gene1, gene2, comp1="_tis1", comp2="_tis2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Input data matrix, where rows are genes and columns are samples. Gene expression measures in the two tissues are combined by row.
}
  \item{gene1}{
    First gene to plot
}
  \item{gene2}{
    Second gene to plot
}
  \item{comp1}{
    Suffix (string character) to identify genes of tissue 1. Default=\code{"_1"}.
}
  \item{comp2}{
    Suffix (string character) to identify genes of tissue 2. Default=\code{"_2"}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A plot with four panels displaying the scatterplots of gene1 and gene2 relationships within and between tissues/compartments.
}
\references{
}
\author{
Aurora Savino, Raffaele M. Iannzzi.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
