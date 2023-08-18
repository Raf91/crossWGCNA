\name{ST_plot_expr}
\alias{ST_plot_expr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot gene expression in space
}
\description{
Generates a plot showing the expression of a gene in space and the boundaries of the two compartments.
}
\usage{
ST_plot_expr(gene, averaged_expr_all, coords, included_spots, tis1_spots, tis2_spots, midpoints)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene}{
    Gene to be plotted.
}
  \item{averaged_expr_all}{
    Gene expression data matrix where rows are genes and columns are spots. Output of the \code{\link{ST_expr_smooth}} function.
}
  \item{coords}{
    X and y coordinates for each spot, obtained from the \code{\link{ST_spots_coords}} function.
}
  \item{included_spots}{
    Spots of compartment 1 along the boundaries of compartment 1 and 2.
}
  \item{tis1_spots}{
    Vector of indexes of spots belonging to tissue 1
}
  \item{tis2_spots}{
    Vector of indexes of spots belonging to tissue 2
}
  \item{midpoints}{
    List of two vectors with x and y coordinates of the midpoints between neighbouring compartment 1 and compartment 2 spots. 
    Output of the \code{\link{ST_midpoints_df}} function.
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
\code{\link{ST_spots_coords}},
\code{\link{ST_spots_filt}}, 
\code{\link{ST_boundary_spots}}, 
\code{\link{ST_weighted_mod}}.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
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
