\name{ST_spots_filt}
\alias{ST_spots_filt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Select non-isolated spots.
}
\description{
  Selects spots with at least one neighbouring spot belonging to the same tissue compartment for spatial transcriptomics data.
}
\usage{
  ST_spots_filt(coords, tis1_spots, tis2_spots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coords}{
    X and y coordinates for each spot, obtained from the \code{\link{ST_spots_coords}} function.
}
  \item{tis1_spots}{
    Vector of indexes of spots belonging to tissue 1
}
  \item{tis2_spots}{
    Vector of indexes of spots belonging to tissue 2
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List of two vectors: indexes of selected spots for tissue 1; indexes of selected spots for tissue 2.
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
\code{\link{ST_plot_expr}},
\code{\link{ST_plot_comm}},
\code{\link{ST_spots_coords}},
\code{\link{ST_spots_filt}}, 
\code{\link{ST_boundary_spots}}, 
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
