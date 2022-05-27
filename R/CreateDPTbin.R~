#' @title 
#' Computes the  Ordinary Covariance Entropy (CE)
#' 
#' @aliases CompCE
#'  
#' @description 
#' This is the main user function for computing the Covariance
#' Entropy (CE). It takes as input the normalized gene expression data matrix
#' defined over the regulatory factors of interest (rows) and single cells (columns).
#' 
#' @param data.m
#' The input data matrix, representing the normalized gene expression data matrix
#' over regulatory factors (rows) and single cells (columns).
#' 
#' @param useSpear
#' A logical to specify whether to use Spearman correlations when computing the correlations between regulatory factors/genes. By default this is FALSE
#'
#' @param std
#' A logical to specify whether the expression profiles should be standardized to unit variance.
#' 
#' @return A list with three elements
#' 
#' @return ce
#' The CE value with constant term removed
#'
#' @return ceN
#' The CE value with constant term (i.e. taking the dimensionality of the dataset into account)
#' 
#' @return evals
#' The eigenvalues of the covariance matrix
#' 
#' #' 
#' @references 
#' Qi Luo, Alok K Maity, Andrew E Teschendorff
#' \emph{Distance Covariance Entropy reveals primed states and bifurcation dynamics in single-cell RNA-Seq data}
#' Submitted.
#' 
#' 
#' @examples 
#'
#' 
#' @import corpcor
#' 
#' @export
#'
CompCE <- function(data.m, useSpear=FALSE, std=TRUE){

  n <- ncol(data.m);
  p <- nrow(data.m);
  data.m <- data.m - rowMeans(data.m);
  ### standardise rows to unit variance (mean is already 0)
  if(std){
   sd.v <- apply(data.m,1,sd);
   data.m <- data.m/sd.v; ### now rows have unit variance
  }
  if(nrow(data.m) <= 0.5*n*(n-1)){
    if(useSpear==FALSE){
        cov.m <- cov(t(data.m));
    }
    else if (useSpear){
        cov.m <- cor(t(data.m),method="spear");
    }
  }
  else {
    print("Using shrinkage");
    if(useSpear==FALSE){
       cov.m <- cov.shrink(t(data.m));
    }
    else if (useSpear){
        print("The Spearman option is not possible if the number of regulatory factors is greater than the number of distance elements, i.e. 0.5*ncells*(ncells-1). Please set useSpear to FALSE!");
    }
  }
  eigen.o <- eigen(cov.m,symmetric=TRUE,only.values=TRUE);
  ce <- 0.5*sum(log(eigen.o$values));
  ceN <- ce + 0.5*p*(1+log(2*pi));
  return(list(ce=ce,ceN=ceN,evals=eigen.o$values));
}




