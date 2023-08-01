\name{cor_inspect}
\alias{cor_inspect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  cor_inspect
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cor_inspect(data = data_merged_GSE88715, gene1, gene2, comp1 = "_tis1", comp2 = "_tis2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{gene1}{
%%     ~~Describe \code{gene1} here~~
}
  \item{gene2}{
%%     ~~Describe \code{gene2} here~~
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
function (data = data_merged_GSE88715, gene1, gene2, comp1 = "_tis1", 
    comp2 = "_tis2") 
{
    df <- data.frame(gene1 = c(data[paste(gene1, comp1, sep = ""), 
        ], data[paste(gene1, comp2, sep = ""), ], data[paste(gene1, 
        comp1, sep = ""), ], data[paste(gene1, comp2, sep = ""), 
        ]), gene2 = c(data[paste(gene2, comp1, sep = ""), ], 
        data[paste(gene2, comp2, sep = ""), ], data[paste(gene2, 
            comp2, sep = ""), ], data[paste(gene2, comp2, sep = ""), 
            ]), compartment = c(rep(c("comp1 vs comp1", "comp2 vs comp2", 
        "comp1 vs comp2", "comp2 vs comp1"), each = ncol(data))))
    p <- ggplot(df, aes(x = gene1, y = gene2)) + geom_point() + 
        facet_wrap(. ~ compartment) + geom_smooth(method = "lm") + 
        theme_classic() + labs(x = gene1, y = gene2)
    return(p)
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
