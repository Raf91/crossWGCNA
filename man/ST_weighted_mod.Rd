\name{ST_weighted_mod}
\alias{ST_weighted_mod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ST_weighted_mod
}
\description{
Computes a weighted summmary of a module's expression for spatial transcriptomics data.
}
\usage{
ST_weighted_mod(modules, kw, mod_sel, averaged_expr_all, 
                comp1="_tis1", comp2="_tis2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modules}{
    Vector with module labels for each gene.
}
  \item{kw}{
    Connectivities computed with the function kwithin.
}
  \item{mod_sel}{
    Label of the selected module.
}
  \item{averaged_expr_all}{
    Gene expression data matrix where rows are genes and columns are spots. Output of the \code{\link{ST_expr_smooth}} function.
}
  \item{comp1}{
    Suffix (string character) to identify genes of tissue 1. Default=\code{"_tis1"}.
}
  \item{comp2}{
    Suffix (string character) to identify genes of tissue 2. Default=\code{"_tis2"}.
}
}
\details{
The computation is based on the kExt of the genes measured in each compartment separately.
}
\value{
List of two vectors with the weighted mean of the expression of module's genes belonging to compartment 1 or compartment 2.
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
\code{\link{ST_boundary_spots}}.
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
