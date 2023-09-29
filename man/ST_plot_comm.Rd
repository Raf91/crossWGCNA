\name{ST_plot_comm}
\alias{ST_plot_comm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    ST_plot_comm
}
\description{
Generates a plot showing the communication score of two genes along the boundries between the two tissue compartments for spatial transcriptomics data.
}
\usage{
ST_plot_comm(gene1, gene2, averaged_expr_all, coords, included_spots, 
             sel_spots, tis1_spots, tis2_spots, midpoints)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene1}{
    First gene, measured in compartment 1.
}
 \item{gene2}{
    Second gene, measured in compartment 2.
}
 \item{averaged_expr_all}{
    Gene expression data matrix where rows are genes and columns are spots. Output of the \code{ST_expr_smooth} function.
}
 \item{coords}{
    X and y coordinates for each spot, obtained from the \code{ST_spots_coords} function.
}
 \item{included_spots}{
    Spots of compartment1 along the boundaries of compartment 1 and 2.
}
 \item{sel_spots}{
    List of two vectors of indexes (matched) for compartment 1 and compartment 2 spots. 
    Output of the \code{ST_boundary_spots} function.
}
 \item{tis1_spots}{
    Indexes of compartment 1 spots.
}
 \item{tis2_spots}{
    Indexes of compartment 2 spots.
}
 \item{midpoints}{
    List of two vectors with x and y coordinates of the midpoints between 
    neighbouring compartment 1 and compartment 2 spots. Output of the \code{ST_midpoints_df} function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a plot.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Aurosa Savino, Raffaele M. Iannuzzi.
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
