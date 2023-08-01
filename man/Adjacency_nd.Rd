\name{Adjacency_nd}
\alias{Adjacency_nd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Adjacency_nd (netdiff)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
    Adjacency_nd(data, comp1 = "_1", comp2 = "_2", Adj_type = "signed", cortype = "spearman", pval = "none", thr = 0.05, beta = 6, sign_list = 1, compartment_sel = "none", selgenes = NA)
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
function (data, comp1 = "_1", comp2 = "_2", Adj_type = "signed", 
    cortype = "spearman", pval = "none", thr = 0.05, beta = 6, 
    sign_list = 1, compartment_sel = "none", selgenes = NA) 
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
    A_orig <- A
    avgpath <- matrix(ncol = nrow(A)/2, nrow = nrow(A)/2)
    for (x in 1:(nrow(A)/2)) {
        avgpath[x, ] <- (A[grep(comp1, rownames(A))[x], grep(comp1, 
            rownames(A))] + A[grep(comp2, rownames(A))[x], grep(comp2, 
            rownames(A))])/2
    }
    A[grep(comp1, rownames(A)), grep(comp2, rownames(A))] <- A[grep(comp1, 
        rownames(A)), grep(comp2, rownames(A))] - avgpath
    A[grep(comp2, rownames(A)), grep(comp1, rownames(A))] <- A[grep(comp2, 
        rownames(A)), grep(comp1, rownames(A))] - avgpath
    diff <- abs(A_orig[grep(comp1, rownames(A_orig)), grep(comp2, 
        rownames(A_orig))]) - abs(A[grep(comp1, rownames(A)), 
        grep(comp2, rownames(A))])
    A[grep(comp1, rownames(A)), grep(comp2, rownames(A))][diff < 
        0] <- A_orig[grep(comp1, rownames(A_orig)), grep(comp2, 
        rownames(A_orig))][diff < 0]
    diff <- abs(A_orig[grep(comp2, rownames(A_orig)), grep(comp1, 
        rownames(A_orig))]) - abs(A[grep(comp2, rownames(A)), 
        grep(comp1, rownames(A))])
    A[grep(comp2, rownames(A)), grep(comp1, rownames(A))][diff < 
        0] <- A_orig[grep(comp2, rownames(A_orig)), grep(comp1, 
        rownames(A_orig))][diff < 0]
    A <- A/2
    if (compartment_sel == "comp2") {
        if (is.na(selgenes)) {
            stop("No selected genes provided")
        }
        sign_list <- sign_list[which(!is.na(sign_list))]
        names(sign_list) <- paste(names(sign_list), gsub("\\$", 
            "", comp2), sep = "")
        selgenes <- intersect(paste(selgenes, gsub("\\$", "", 
            comp2), sep = ""), rownames(A))
        A <- A[c(grep(comp1, rownames(A)), which(rownames(A) \%in\% 
            selgenes)), c(grep(comp1, rownames(A)), which(rownames(A) \%in\% 
            selgenes))]
        A[grep(comp2, rownames(A)), grep(comp1, rownames(A))] <- A[grep(comp2, 
            rownames(A)), grep(comp1, rownames(A))] * (sign_list[rownames(A)[grep(comp2, 
            rownames(A))]])
        A[grep(comp1, rownames(A)), grep(comp2, rownames(A))] <- t(A[grep(comp2, 
            rownames(A)), grep(comp1, rownames(A))] * (sign_list[rownames(A)[grep(comp2, 
            rownames(A))]]))
    }
    else if (compartment_sel == "comp1") {
        if (is.na(selgenes)) {
            stop("No selected genes provided")
        }
        sign_list <- sign_list[which(!is.na(sign_list))]
        names(sign_list) <- paste(names(sign_list), gsub("\\$", 
            "", comp1), sep = "")
        selgenes <- intersect(paste(selgenes, gsub("\\$", "", 
            comp1), sep = ""), rownames(A))
        A <- A[c(which(rownames(A) \%in\% selgenes), grep(comp2, 
            rownames(A))), c(which(rownames(A) \%in\% selgenes), 
            grep(comp2, rownames(A)))]
        A[grep(comp1, rownames(A)), grep(comp2, rownames(A))] <- A[grep(comp1, 
            rownames(A)), grep(comp2, rownames(A))] * (sign_list[rownames(A)[grep(comp1, 
            rownames(A))]])
        A[grep(comp2, rownames(A)), grep(comp1, rownames(A))] <- t(A[grep(comp1, 
            rownames(A)), grep(comp2, rownames(A))] * (sign_list[rownames(A)[grep(comp1, 
            rownames(A))]]))
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
