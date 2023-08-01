\name{Adjacency_sl}
\alias{Adjacency_sl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Adjacency_sl (selfloop)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Adjacency_sl(data, comp1 = "_1", comp2 = "_2", Adj_type = "signed", cortype = "spearman", pval = "none", thr = 0.05, beta = 6, sign_list = 1, which.sign = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{comp1}{
%%     ~~Describe \code{comp1} here~~
}
  \item{comp2}{
%%     ~~Describe \code{comp2} here~~
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
  \item{sign_list}{
%%     ~~Describe \code{sign_list} here~~
}
  \item{which.sign}{
%%     ~~Describe \code{which.sign} here~~
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
function (data, comp1 = "_1", comp2 = "_2", Adj_type = "signed", 
    cortype = "spearman", pval = "none", thr = 0.05, beta = 6, 
    sign_list = 1, which.sign = "none") 
{
    comp1 <- paste(comp1, "$", sep = "")
    comp2 <- paste(comp2, "$", sep = "")
    if (pval == "none") {
        if (cortype == "bicor") {
            A <- bicor(t(data))
        }
        else {
            A <- cor(t(data), method = cortype)
        }
    }
    else {
        if (cortype == "bicor") {
            paste("Can't use pval with bicor")
        }
        else {
            mat <- rcorr(t(data), type = cortype)
            A <- mat[[1]]
            p.val <- mat[[3]]
            rm(mat)
            if (pval == "threhsold") {
                A[which(p.val > thr)] <- NA
            }
            else if (pval == "weight") {
                A <- A * (1 - p.val)
            }
        }
    }
    if (which.sign == "comp2") {
        A[grep(comp2, rownames(A)), grep(comp1, rownames(A))] <- A[grep(comp2, 
            rownames(A)), grep(comp1, rownames(A))] * sign_list
        A[grep(comp1, rownames(A)), grep(comp2, rownames(A))] <- t(A[grep(comp2, 
            rownames(A)), grep(comp1, rownames(A))] * sign_list)
    }
    else if (which.sign == "comp1") {
        A[grep(comp1, rownames(A)), grep(comp2, rownames(A))] <- A[grep(comp1, 
            rownames(A)), grep(comp2, rownames(A))] * sign_list
        A[grep(comp2, rownames(A)), grep(comp1, rownames(A))] <- t(A[grep(comp1, 
            rownames(A)), grep(comp2, rownames(A))] * sign_list)
    }
    if (Adj_type == "signed") {
        A <- (0.5 * (1 + A))^beta
    }
    else if (Adj_type == "unsigned") {
        A <- (abs(A))^beta
    }
    else if (Adj_type == "keep sign") {
        A <- ((abs(A))^beta) * sign(A)
    }
    return(A)
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
