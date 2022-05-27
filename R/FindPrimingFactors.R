#' @title 
#' Identifies putative priming regulatory factor pairs
#' 
#' @aliases FindPrimingFactors
#'  
#' @description 
#' If priming is present (as evaluated using \code{CompPrimingPval} function), then
#' one wishes to find the regulatory factors driving such priming. This function ranks
#' regulatory factors and regulatory factor pairs according to their likelihood of being
#' involved in priming.
#' 
#' 
#' @param data.m
#' The input data matrix, representing the normalized gene expression data matrix
#' over candidate regulatory factors (rows) and single cells (columns).
#' 
#' 
#' @return A list of four elements:
#' 
#' @return rankedLOO
#' A vector of regulatory factors ranked in decreasing order of priming potential
#'
#' @return rankedPairs
#' A table ranking the most likely pairs of regulatory factors involved in priming,
#' with 1st column labeling the Spearman Rank Correlation and 2nd column labeling the
#' corresponding P-value.
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
#' 
#' @export
#'
FindPrimingFactors <- function(data.m){

  diceObs <- CompDICE(data.m)$dice;

  diceRmO.v <- vector();
  print("Doing perturbation LOO DICE analysis");
  for(rf in 1:nrow(data.m)){
    tmp.idx <- setdiff(1:nrow(data.m),rf);
    diceRmO.v[rf] <- CompDICE(data.m[tmp.idx,])$dice;
  }
  names(diceRmO.v) <- rownames(data.m);
  dDice.v <- diceRmO.v - diceObs;
  tmp.s <- sort(dDice.v,decreasing=TRUE,index.return=TRUE);
  

  print("Doing Spearman analysis");
  npairs <- 0.5*nrow(data.m)*(nrow(data.m)-1);
  tab.m <- matrix(nrow=npairs,ncol=2);
  colnames(tab.m) <- c("SCC","P");
  tmp.v <- vector();
  ii <- 1;    
  for(tf1 in 1:(nrow(data.m)-1)){
   for(tf2 in (tf1+1):nrow(data.m)){
       ct.o <- suppressWarnings(cor.test(data.m[tf1,],data.m[tf2,],method="spearman"));
       tab.m[ii,1:2] <- c(ct.o$est,ct.o$p.value);
       tmp.v[ii] <- paste(rownames(data.m)[tf1],":",rownames(data.m)[tf2],sep="");
       ii <- ii+1;
   }
  }
  rownames(tab.m) <- tmp.v;
  tab.s <- sort(tab.m[,1],decreasing=TRUE,index.return=TRUE);

  return(list(rankedLOO=names(dDice.v)[tmp.s$ix],rankedPairs=tab.m[tab.s$ix,]));

}

