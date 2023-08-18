\name{Adjacency}
\alias{Adjacency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Adjacency
}
\description{
  Function to compute cross-tissue gene expression adjacencies. 
}
\usage{
  Adjacency(data, method="selfloop",comp1="_1",comp2="_2",
            Adj_type="signed", cortype="spearman", pval="none", 
            thr=0.05, beta=6, sign_list=1, compartment_sel="none", 
            selgenes=NA, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Input data matrix, where rows are genes and columns are samples. Gene expression measures in the two tissues are combined by row.
  }
  \item{method}{
    String character identifying the method to use to compute cross-tissue adjacencies. Can be \code{"self-loop"} or \code{"netdiff"}.
  }
  \item{comp1}{
    Suffix (string character) to identify genes of tissue 1. Default=\code{"_1"}.
  }
  \item{comp2}{
    Suffix (string character) to identify genes of tissue 2. Default=\code{"_2"}.
  }
  \item{Adj_type}{
    String character identifying the type of adjacency to compute. Can be \code{"signed"}, \code{"unsigned"} or \code{"keep sign"}. Default=\code{"signed"}. 
    \code{"keep sign"} is equivalent to \code{"unsigned"} in absolute values, but with sign corresponding to the sign of correlation.
    See WGCNA package for details.
  }
  \item{cortype}{
    String character identifying the type of correlation. Can be \code{"pearson"}, \code{"spearman"} or \code{"bicor"}. Default=\code{"spearman"}. 
    See WGCNA package for details.
  }
  \item{pval}{
    String character identifying whether to use the correlations' p-value in edge weight.
    Possible values are: \code{"none"} (does not include p-values), \code{"threshold"} (discard edges with a p-value higher than a threshold),
    \code{"weight"} (multiply edges' adjacencies by 1-pvalue). Default=\code{"none"}.
    Can be used with \code{"pearson"} and \code{"spearman"} correlations.
  }
  \item{thr}{
    Threshold for cutting edges based on p-values. Default=\code{0.05}.
  }
  \item{beta}{
    Power for adjacency calculation. Default=\code{6}.
  }
  \item{sign_list}{
    When considering genes associated with a pathway with sign (activation: +1, inactivation: -1), a vector with a value for each gene in the pathway. 
    Default=\code{1}, meaning that each gene has the same positive weight.
  }
  \item{compartment_sel}{
    String character identifying the tissue to be considered when providing a pathway. Can be \code{"none"}, \code{"comp1"} or \code{"comp2"}.
  }
  \item{selgenes}{
    Vector of genes associated with the pathway of interest.
  }
  \item{verbose}{
    Boolean value. If \code{TRUE} prints logging output. Default=\code{FALSE}.
  }
}
\details{
  The adjacencies are computed accordingly to the method selected (i.e., \code{"self-loop"} or \code{"netdiff"}).
}
\value{
Returns a matrix with pair-wise gene adjacencies.
}
\references{
Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 559.
}
\author{
Aurora Savino, Raffaele M. Iannuzzi.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{rm_selfloop}},
  \code{\link{rm_netdiff}}
}
\examples{
See package vignette.
}


