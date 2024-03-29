% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindCandidateFactors.R
\name{FindCandidateFactors}
\alias{FindCandidateFactors}
\title{Identifies all candidate priming regulatory factors}
\usage{
FindCandidateFactors(data.m, rf.v, celltypes.v, mpp, linCT.v, sigth =
NULL, type1fth="Mean")
}
\arguments{
\item{data.m}{The input data matrix, representing the normalized gene expression data matrix
over genes (rows) and single cells (columns).}

\item{rf.v}{A vector of regulatory factor identifiers. Must be the same identifier as used for
the rownames of \code{data.m}.}

\item{celltypes.v}{A vector of cell-type identifiers. Matches columns of \code{data.m}.}

\item{mpp}{Character specifying the cell-type that defines the "multi-potent"
cell population, i.e. the population where priming is to be assessed.}

\item{linCT.v}{Character vector specifying the competing downstream lineages, for which
priming is to be assessed.}

\item{sigth}{Significance threshold for Benjamini-Hochberg adjusted P-values. The default value NULL makes a Bonferroni correction at 0.05 level.}


\item{type1fth}{The threshold on the frequency of expression among multipotent cells to filter RFs. Only RFs
 with a frequency lower than the threshold are considered for downstream analyses. For the default value "Mean", the mean expression frequency over all differentially expressed RFs is taken. If NULL, no filtering is performed.}
}

\value{
\item{candRF.lv}{A list of vectors with regulatory factor identifiers
  and their frequencies of expression in MPPs for candidate priming RF.
  Each list entry corresponds to one competing lineage.}
}

\description{
This function performs differential overexpression analysis in order to 
identify candidate regulatory factors that may play a role in priming of 
various lineages. Optionally, it also filters candidate regulatory factors
by the frequency of expression across multipotent cells.
}
\examples{


}
\references{
Qi Luo, Alok K Maity, Andrew E Teschendorff
\emph{Distance Covariance Entropy reveals primed states and bifurcation dynamics in single-cell RNA-Seq data}
Submitted.
}
