\docType{data}
\name{yangLiver}
\alias{yangLiver}
\title{scRNA-Seq data from Yang et al}
\format{An expressionset dataset object containing data from Yang et al defined over
  22 liver-specific transcription factors and 447 cells from 7 embryonic
  developmental timepoints.
}

\usage{
data(yangLiver)
}
\description{
This expression dataset object contains differentiation activity (TFA) values
for 22 liver-specific regulatory factors and 447 cells, drawn from 7
embryonic developmental timepoints (E10.5, E11.5, E12.5, E13.5, E14.5,
E15.5, E17.5). 
}

\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{yangLiver}: The expression dataset object. Cell-type info
  can be accessed as usual with \code{phenoData}. TFA-values can be
  accessed \code{exprs} function.
}
}
\references{
Yang et al. A single-cell transcriptomic analysis reveals precise
pathways and regulatory mechanisms underlying hepatoblast
differentiation, Hepatology 2017.
}
\keyword{datasets}
