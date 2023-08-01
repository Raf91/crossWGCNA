\name{changenames}
\alias{changenames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  changenames
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
changenames(data, anno)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{anno}{
%%     ~~Describe \code{anno} here~~
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
function (data, anno) 
{
    annotation_sel = anno[match(rownames(data), anno[, 1]), 2]
    if (length(which(annotation_sel == "")) > 0) {
        data <- data[-which(annotation_sel == ""), ]
        annotation_sel <- annotation_sel[-which(annotation_sel == 
            "")]
    }
    a <- which(duplicated(annotation_sel))
    while (length(a) > 0) {
        for (i in 1:length(unique(annotation_sel))) {
            if (length(which(annotation_sel == unique(annotation_sel)[i])) > 
                1) {
                m = which.max(rowMeans(data[which(annotation_sel == 
                  unique(annotation_sel)[i]), ], na.rm = T))
                data = data[-which(annotation_sel == unique(annotation_sel)[i])[-m], 
                  ]
                annotation_sel = annotation_sel[-which(annotation_sel == 
                  unique(annotation_sel)[i])[-m]]
            }
        }
        data = data[which(is.na(annotation_sel) == F), ]
        annotation_sel = na.omit(annotation_sel)
        a <- which(duplicated(annotation_sel))
    }
    rownames(data) = annotation_sel
    return(data)
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
