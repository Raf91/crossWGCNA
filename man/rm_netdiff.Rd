\name{rm_netdiff}
\alias{rm_netdiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  rm_netdiff
}
\description{
Removes the average connectivity between each gene pair from both tissues’ networks. 
}
\usage{
rm_netdiff(A, comp1="_tis1", comp2="_tis2", verbose=TRUE)
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
  \item{verbose}{
    Boolean value. If \code{TRUE} prints logging output. Default=\code{FALSE}.
}
}
\details{
See the manuscript for further details.
}
\value{
Returns the adjusted Adjacency matrix. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Aurora Savino, Raffaele M. Iannuzzi.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rm_netdiff}}
}
\examples{

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
