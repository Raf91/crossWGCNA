\name{ST_merged_dataset}
\alias{ST_merged_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Creates a merged dataset with gene expression from the two selected compartments
}
\description{
  For each selected spot in compartment 1, takes the neighbouring spots in compartment 2 and creates a gene expression matrix with matched gene expression from the two selected compartments.
}
\usage{
ST_merged_dataset(sel_spots, coords, averaged_expr_all, var_thr=0.75, comp1="_tis1", comp2="_tis2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sel_spots}{
    List of two vectors, one for each compartment, with indexes of spots to be included. Output of the \code{\link{ST_spots_filt function}}.
}
  \item{coords}{
     X and y coordinates for each spot, obtained from the \code{\link{ST_spots_coords}} function.
  }
  \item{averaged_expr_all}{
    Gene expression data matrix where rows are genes and columns are spots. Output of the \code{\link{ST_expr_smooth}}function.
  }
  \item{var_thr}{
    Threshold (numeric) with the lowest gene expression variance percentile to keep a gene in the final data matrix.
    Default=\code{0.75}.
  }
  \item{comp1}{
    Suffix (string character) to identify genes of tissue 1. Default=\code{"_tis1"}.
  }
  \item{comp2}{
    Suffix (string character) to identify genes of tissue 2. Default=\code{"_tis2"}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List of two values:
  \item{data_merged}{
    Matrix with gene expression from the two selected compartments, where columns represent positions in space (of selected compartment 1 spots), and rows are genes measured either in one or in the other compartment (denoted by different suffixes).
}
  \item{included_spots}{
    Vector of indexes for the compartment 1 spots included in the final data matrix.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Aurora Savino, Raffaele M. Iannuzzi
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ST_spots_coords}}, \code{\link{ST_expr_smooth}}, 
\code{\link{ST_spots_filt}}
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
