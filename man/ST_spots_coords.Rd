\name{ST_spots_coords}
\alias{ST_spots_coords}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ST_spots_coords
}
\description{
Defines the x and y coordinates of each spot in the array, assigning the value of 1 to the leftmost spot and considering an inter-spot x distance of 1.
}
\usage{
ST_spots_coords(data,br=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Seurat object with the spatial transcriptomics data.
}
\item{br}{
    Number of breaks (numeric) for finding the coordinate peaks. Needs to be higher than the number of actual x or y coordinates of spots in the array.
    Default=\code{1000}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a matrix with two columns for x and y coordinates of each spot.
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
\code{\link{ST_plot_comm}}
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
