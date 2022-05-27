\docType{data}
\name{ranzoniBranchPr}
\alias{ranzoniBranchPr}
\title{Cell-fate probability data towards MEMP and GP lineages}
\description{
A matrix of differentiation probabilities towards MEMPs and GPs for 1200 (multipotent
progenitors) MPPs, 177 megakaryote-erythrocyte-mast progenitors (MEMPs) and 192 granulocyte progenitors (GPs), obtained by using Palantir, a Python module, on the scRNA-Seq data from Ranzoni et al. 
}
\usage{data(ranzoniBranchPr)}
\format{
  A matrix of differentiation probabilities, rows matched with cells from \code{ranzoni.m} and columns specifying the lineage.
}

\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{ranzoni.branch.pr.m}: matrix of differentiation probabilities towards MEMPs and granulocyte progenitors GPs for 1200 MPPs, 192 MEMPs and 177 GPs.
}
}
\keyword{datasets}
