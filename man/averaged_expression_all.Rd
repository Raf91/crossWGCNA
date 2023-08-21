\name{ST_expr_smooth_output}
\alias{ST_expr_smooth_output}
\docType{data}
\title{
Averaged gene expression of \code{BC@assays$SCT@data}
}
\description{
Pre-computed smooths spatial transcriptomics gene expression using the weighted mean of neighbouring spots in one compartment.
}
\usage{load("ST_expr_smooth_output.RData")}
\format{
  A data frame 20227 genes x 4898 spots
}
\details{
  Output of ST_expr_smooth function. 
  The object is available for download at zenodo web page: code{\url{https://zenodo.org/record/8268805}} . 
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
  options(timeout=100)
  download.file(url="https://zenodo.org/record/8268805/files/ST_expr_smooth_out.RData?download=1",destfile="~/ST_expr_smooth_out.RData")
  load("~/ST_expr_smooth_output.RData")
}
\keyword{datasets}
