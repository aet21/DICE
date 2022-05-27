#' @title 
#' Bins pseudotime into bins of equal cell numbers if possible.
#' 
#' @aliases CreateDPTbin
#'  
#' @description 
#' Given as input a vector of pseudotime values, and the desired number of bins,
#' it places cells into bins of increased pseudotime values.
#' 
#' @param dpt.v
#' The input vector of pseudotime values.
#' 
#' @param minBins
#' The number of bins required. The actual number will be at least as big as this one.
#'
#' 
#' @return A list with two elements
#' 
#' @return binList
#' A list of cellular indices, one for each bin, in increasing order of pseudotime.
#'
#' @return idx
#' A cellular index vector specifying the bin each cell falls into.
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
#' 
#' @export
#'
CreateDPTbin <- function(dpt.v,minBins=7){

  avNcells <- floor(length(dpt.v)/minBins);
  dpt.s <- sort(dpt.v,decreasing=FALSE,index.return=TRUE);
  i <- 1;
  f <- avNcells;
  bin.li <- list();
  bin.idx <- vector();
  for(bin in 1:(minBins-1)){
        bin.idx[dpt.s$ix[i:f]] <- bin;
        bin.li[[bin]] <- dpt.s$ix[i:f];
        i <- f+1;
        f <- i+avNcells-1
  }
  bin.idx[dpt.s$ix[i:length(dpt.v)]] <- minBins;
  bin.li[[minBins]] <- dpt.s$ix[i:length(dpt.v)];
  names(bin.li) <- 1:minBins;
  return(list(binList=bin.li,idx=bin.idx));
}




