\name{cytoscape_net}
\alias{cytoscape_net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  cytoscape_net
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cytoscape_net(adjacency = adj_GSE10797, dataset = data_merged_GSE10797, gene, comp1, comp2, num, corr = "spearman")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjacency}{
%%     ~~Describe \code{adjacency} here~~
}
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{gene}{
%%     ~~Describe \code{gene} here~~
}
  \item{comp1}{
%%     ~~Describe \code{comp1} here~~
}
  \item{comp2}{
%%     ~~Describe \code{comp2} here~~
}
  \item{num}{
%%     ~~Describe \code{num} here~~
}
  \item{corr}{
%%     ~~Describe \code{corr} here~~
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
function (adjacency = adj_GSE10797, dataset = data_merged_GSE10797, 
    gene, comp1, comp2, num, corr = "spearman") 
{
    interactors <- c(names(sort(adjacency[paste(gene, comp1, 
        sep = "_"), grep(comp2, colnames(adjacency))], decreasing = T))[1:num])
    inter_edges <- cor(dataset[paste(gene, comp1, sep = "_"), 
        ], t(dataset[intersect(interactors, rownames(dataset)), 
        ]), method = corr)
    intra_1 <- cor(dataset[paste(gene, comp1, sep = "_"), ], 
        t(dataset[gsub(comp2, comp1, intersect(interactors, rownames(dataset))), 
            ]), method = corr)
    intra_2 <- cor(dataset[paste(gene, comp2, sep = "_"), ], 
        t(dataset[intersect(interactors, rownames(dataset)), 
            ]), method = corr)
    df <- data.frame(Source = c(rep(gene, length(inter_edges)), 
        rep(gene, length(intra_1)), rep(gene, length(intra_2))), 
        Target = c(gsub(paste("_", comp2, sep = ""), "", colnames(inter_edges)), 
            gsub(paste("_", comp1, sep = ""), "", colnames(intra_1)), 
            gsub(paste("_", comp2, sep = ""), "", colnames(intra_2))), 
        Weight = c(inter_edges, intra_1, intra_2), Edge_type = c(rep("inter", 
            length(inter_edges)), rep("intra1", length(intra_1)), 
            rep("intra2", length(intra_2))), Source_type = c(rep(comp1, 
            length(inter_edges)), rep(comp1, length(intra_1)), 
            rep(comp2, length(intra_2))), Target_type = c(rep(comp2, 
            length(inter_edges)), rep(comp1, length(intra_1)), 
            rep(comp2, length(intra_2))))
    write.csv(df, file = paste("Cytoscape", gene, comp1, "egdes.csv", 
        sep = "_"))
    return(df)
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
