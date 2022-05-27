#' @title 
#' Computes the  Distance Covariance Entropy (DICE)
#' 
#' @aliases CompDICE
#'  
#' @description 
#' This is the main user function for computing the Distance Covariance
#' Entropy (DICE). It takes as input the normalized gene expression data matrix
#' defined over the regulatory factors of interest (rows) and single cells (columns).
#' 
#' @param data.m
#' The input data matrix, representing the normalized gene expression data matrix
#' over regulatory factors (rows) and single cells (columns).
#' 
#' @param useSpear
#' A logical to specify whether to use Spearman correlations when computing the correlation between the distance vectors. By default this is FALSE
#'
#' @param std
#' A logical to specify whether the distance vectors should be standardized to unit variance.
#' 
#' @return A list with three elements
#' 
#' @return dice
#' The DICE value with constant term removed
#'
#' @return diceN
#' The DICE value with constant term (i.e. taking the dimensionality of the dataset into account)
#' 
#' @return evals
#' The eigenvalues of the distance covariance matrix
#' 
#' 
#' 
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
CompDICE <- function(data.m, useSpear=FALSE, std=TRUE){

  n <- ncol(data.m);
  p <- nrow(data.m);
  ### for each regulatory factor (gene), compute pairwise distances across cells
  dist.m <- t(apply(data.m,1,CompDistElem));
  dist.m <- dist.m - rowMeans(dist.m);
  ### standardise rows to unit variance (mean is already 0)
  if(std){
   sd.v <- apply(dist.m,1,sd);
   dist.m <- dist.m/sd.v; ### now rows have unit variance
  }
  if(nrow(data.m) <= 0.5*n*(n-1)){
    if(useSpear==FALSE){
        cov.m <- cov(t(dist.m));
    }
    else if (useSpear){
        cov.m <- cor(t(dist.m),method="spear");
    }
  }
  else {
    print("Using shrinkage");
    if(useSpear==FALSE){
       cov.m <- cov.shrink(t(dist.m));
    }
    else if (useSpear){
        print("The Spearman option is not possible if the number of regulatory factors is greater than the number of distance elements, i.e. 0.5*ncells*(ncells-1). Please set useSpear to FALSE!");
    }
  }
  eigen.o <- eigen(cov.m,symmetric=TRUE,only.values=TRUE);
  dice <- 0.5*sum(log(eigen.o$values));
  diceN <- dice + 0.5*p*(1+log(2*pi));
  return(list(dice=dice,diceN=diceN,evals=eigen.o$values));
}

#### Auxilliary functions
#' 
#' 
CompDistElem <- function(v){

    dist0.m <- as.matrix(dist(v,diag=TRUE,upper=TRUE));
    ### now double center
    dist.m <- dist0.m - matrix(rowMeans(dist0.m),nrow=nrow(dist0.m),ncol=ncol(dist0.m),byrow=TRUE) - matrix(colMeans(dist0.m),nrow=nrow(dist0.m),ncol=ncol(dist0.m),byrow=FALSE) + mean(dist0.m);
    dist.v <- as.vector(dist.m[upper.tri(dist.m)]);
    return(dist.v);
}

