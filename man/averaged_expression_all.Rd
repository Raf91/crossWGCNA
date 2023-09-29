\name{ST_expr_smooth_output}
\alias{ST_expr_smooth_output}
\docType{data}
\title{
Averaged gene expression of \code{BC@assays$SCT@data}
}
\description{
Pre-computed smooths spatial transcriptomics gene expression profiles using the weighted mean of neighbouring spots in one compartment.
Data were derived from a Breast Cancer Visium dataset provided by the 10x genomics.
}
\usage{load("ST_expr_smooth_out.RData")}
\format{
  A data frame \code{20227} genes \code{x} \code{4898} spots. 
}
\details{
  Output of \code{ST_expr_smooth} function. 
  The object is available for download at zenodo web page: code{\url{https://zenodo.org/record/8268805}} . 
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
  options(timeout=200)
  download.file(
    url="https://zenodo.org/record/8268805/files/ST_expr_smooth_out.RData?download=1",
    destfile="~/ST_expr_smooth_out.RData")
  load("~/ST_expr_smooth_out.RData")
}
\seealso{
\code{\link{ST_expr_smooth}}
}
\keyword{datasets}
