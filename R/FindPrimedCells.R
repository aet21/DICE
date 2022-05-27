#' @title 
#' Identifies the putative cells primed to a given lineage
#' 
#' @aliases FindPrimedCells
#'  
#' @description 
#' If priming is present (as evaluated using \code{CompPrimingPval} function), and
#' having identified the putative priming factors, this function helps identify the
#' corresponding primed cells
#' 
#' 
#' @param data.m
#' The input data matrix, representing the normalized gene expression data matrix
#' over regulatory factors (rows) and single cells representing MPPs (columns). 
#' 
#' 
#' @param prf.lv
#' A list of vectors specifying the putative regulatory factors involved in priming of each lineage.
#' Each entry of the list corresponds to one lineage and should have a name.
#' 
#' @param thE
#' A scalar threshold on expression for calling a RF expressed in a given cell. By default 0.
#' 
#' @return primedCells.li
#' A list of vectors specifying the indices for the putative primed cells of each lineage.
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
#' 
#' 
#' @export
#'
FindPrimedCells <- function(data.m, prf.lv, thE=0){
 primed.li <- list();
 for(lin in 1:length(prf.lv)){
  tmp.m <- data.m[match(prf.lv[[lin]],rownames(data.m)),];
  primed.idx <- 1:ncol(tmp.m);
  for(r in 1:nrow(tmp.m)){
    primed.idx <- intersect(primed.idx,which(tmp.m[r,]>thE));
  }
  primed.li[[lin]] <- primed.idx;
 }
 ## now remove common primed cells
 primedCells.li <- list();
 for(lin in 1:length(prf.lv)){
    otherLIN <- setdiff(1:length(prf.lv),lin);
    primedCells.li[[lin]] <- setdiff(primed.li[[lin]],unlist(primed.li[otherLIN]));
 }
 names(primedCells.li) <- names(prf.lv);
 return(primedCells.li);
}

