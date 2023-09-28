\name{ST_expr_smooth}
\alias{ST_expr_smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ST_expr_smooth
}
\description{
  Smooths gene expression using the weighted mean of neighbouring spots in the same tissue compartment for spatial transcriptomics data..
}
\usage{
  ST_expr_smooth(expr_data, coords, max_dist, spots_class, sel_class)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr_data}{
    Input data matrix, where rows are genes and columns are spots.
}
  \item{coords}{
    Dataframe with columns specyfing x and y coordinates for each spot, obtained from the \code{\link{ST_spots_coords}} function.
}
  \item{max_dist}{
    Integer specifying maximum spots distance to be included in the weighed mean. Default=\code{5}.
}
  \item{spots_class}{
    Charater vector with labels indicating the compartment to which each spot belongs.
}
  \item{sel_class}{
    Character vector with the labels of the two compartments to analyze. Can be \code{"Epi"} and \code{"Stroma"}.
}
}
\details{
  Input data is obtained from \code{@assay$SCT@data} slot of the Seurat object.
}
\value{
  Smoothed data matrix, where rows are genes and columns are spots.
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
\code{\link{ST_boundary_spots}}, 
\code{\link{ST_merged_dataset}},
\code{\link{ST_midpoints_def}},
\code{\link{ST_plot_comm}},
\code{\link{ST_plot_expr}},
\code{\link{ST_spots_coords}},
\code{\link{ST_spots_filt}}, 
\code{\link{ST_weighted_mod}}.
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
