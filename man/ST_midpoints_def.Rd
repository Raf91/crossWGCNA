\name{ST_midpoints_def}
\alias{ST_midpoints_def}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ST_midpoints_def
}
\description{
  Defines boundaries between two tissues/compartments as the points equally distant from two neighbouring compartment1 and compartment2 spots.
}
\usage{
ST_midpoints_def(coords, sel_spots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coords}{
    X and y coordinates for each spot, obtained from the \code{\link{ST_spots_coords}} function.
}
  \item{sel_spots}{
    List of two vectors of indexes (matched) for compartment 1 and compartment 2 spots. Output of the \code{\link{ST_boundary_spots}} function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List of two vectors with x and y coordinates of the midpoints.
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
\code{\link{ST_plot_comm}},
\code{\link{ST_plot_expr}},
\code{\link{ST_spots_coords}},
\code{\link{ST_spots_filt}}, 
\code{\link{ST_boundary_spots}}, 
\code{\link{ST_weighted_mod}}
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
