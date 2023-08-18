\name{ST_boundary_spots}
\alias{ST_boundary_spots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ST_boundary_spots
}
\description{
  Generates matched vectors of indexes for compartment1 and compartment2 spots.
}
\usage{
  ST_boundary_spots(included_spots, coords)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{included_spots}{
    Spots of compartment1 along the boundaries of compartment 1 and 2.
}
  \item{coords}{
    X and y coordinates for each spot, obtained from the \code{\link{ST_spots_coords}} function.
  }
}
\details{

}
\value{
List of two vectors of indexes (matched) for compartment 1 and compartment 2 spots. Wherever the same spot is neighbouring multiple spots of the alternate compartment, the spot index is repeated to represent both.
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
\code{\link{ST_expr_smooth}}, 
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
