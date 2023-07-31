\name{degrees}
\alias{degrees}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
degrees(A, comp1 = "_1", comp2 = "_2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{comp1}{
%%     ~~Describe \code{comp1} here~~
}
  \item{comp2}{
%%     ~~Describe \code{comp2} here~~
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
function (A, comp1 = "_1", comp2 = "_2") 
{
    comp1 <- paste(comp1, "$", sep = "")
    comp2 <- paste(comp2, "$", sep = "")
    genes1 <- gsub(comp1, "", rownames(A)[grep(comp1, rownames(A))])
    genes2 <- gsub(comp2, "", rownames(A)[grep(comp2, rownames(A))])
    Idx1 <- cbind(grep(comp1, rownames(A)), grep(comp2, rownames(A))[match(genes1, 
        genes2)])
    A[Idx1] <- 0
    Idx2 <- cbind(grep(comp2, rownames(A)), grep(comp1, rownames(A))[match(genes2, 
        genes1)])
    A[Idx2] <- 0
    kTot <- rowSums(A)
    kInt_1 <- rowSums(A[grep(comp1, rownames(A)), grep(comp1, 
        colnames(A))])
    kInt_2 <- rowSums(A[grep(comp2, rownames(A)), grep(comp2, 
        colnames(A))])
    kExt_1 <- rowSums(A[grep(comp1, rownames(A)), grep(comp2, 
        colnames(A))])
    kExt_2 <- rowSums(A[grep(comp2, rownames(A)), grep(comp1, 
        colnames(A))])
    k <- list(kInt1 = kInt_1, kInt2 = kInt_2, kExt1 = kExt_1, 
        kExt2 = kExt_2, kTot1 = kTot[grep(comp1, colnames(A))], 
        kTot2 = kTot[grep(comp2, colnames(A))])
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
