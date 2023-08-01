\name{network}
\alias{network}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  network
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
network(data, method = "netdiff", Adj_type = "signed", cortype = "spearman", pval = "none", thr = 0.05, beta = 6, comp1 = "_1", comp2 = "_2", sign_list = 1, compartment_sel = "none", selgenes = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{Adj_type}{
%%     ~~Describe \code{Adj_type} here~~
}
  \item{cortype}{
%%     ~~Describe \code{cortype} here~~
}
  \item{pval}{
%%     ~~Describe \code{pval} here~~
}
  \item{thr}{
%%     ~~Describe \code{thr} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{comp1}{
%%     ~~Describe \code{comp1} here~~
}
  \item{comp2}{
%%     ~~Describe \code{comp2} here~~
}
  \item{sign_list}{
%%     ~~Describe \code{sign_list} here~~
}
  \item{compartment_sel}{
%%     ~~Describe \code{compartment_sel} here~~
}
  \item{selgenes}{
%%     ~~Describe \code{selgenes} here~~
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
function (data, method = "netdiff", Adj_type = "signed", cortype = "spearman", 
    pval = "none", thr = 0.05, beta = 6, comp1 = "_1", comp2 = "_2", 
    sign_list = 1, compartment_sel = "none", selgenes = NA) 
{
    comp1 <- paste(comp1, "$", sep = "")
    comp2 <- paste(comp2, "$", sep = "")
    Adj <- Adjacency(data = data, method = method, Adj_type = Adj_type, 
        cortype = cortype, pval = pval, thr = thr, beta = beta, 
        comp1 = comp1, comp2 = comp2, sign_list = sign_list, 
        compartment_sel = compartment_sel, selgenes = selgenes)
    k <- degrees(A = Adj, comp1 = comp1, comp2 = comp2)
    return(k)
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
