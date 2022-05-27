#' @title 
#' Identifies all candidate priming regulatory factors
#' 
#' @aliases FindCandidateFactors.R
#'  
#' @description 
#' This function performs differential overexpression analysis in order to 
#' identify candidate regulatory factors that may play a role in priming of 
#' various lineages. Optionally, it also filters candidate regulatory factors
#' by the frequency of expression across multipotent cells.
#' 
#' @param data.m
#' The input data matrix, representing the normalized gene expression data matrix
#' over genes (rows) and single cells (columns).
#' 
#' @param rf.v
#' A vector of regulatory factor identifiers. Must be the same identifier as used for
#' the rownames of \code{data.m}.
#'
#' @param celltypes.v
#' A vector of cell-type identifiers. Matches columns of \code{data.m}.
#'
#' @param mpp
#' Character specifying the cell-type that defines the "multi-potent"
#' cell population, i.e. the population where priming is to be assessed.
#'
#' @param linCT.v
#' Character vector specifying the competing downstream lineages, for which
#' priming is to be assessed.
#' 
#' @param sigth
#' Significance threshold for Benjamini-Hochberg adjusted P-values. The default value NULL
#' makes a Bonferroni correction at 0.05 level. 
#'
#' @param type1fth
#' The threshold on the frequency of expression among multipotent cells to filter RFs. Only RFs
#' with a frequency lower than the threshold are considered for downstream analyses. For the default
#' value "Mean", the mean expression frequency over all differentially expressed RFs is taken.
#' If NULL, no filtering is performed.
#' 
#' 
#' @return candRF.lv
#' A list of vectors with regulatory factor identifiers and their
#' frequencies of expression in MPPs for candidate priming RF.
#' Each list entry corresponds to one competing lineage.
#' 
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
#' @export
#'
FindCandidateFactors <- function(data.m,rf.v,celltypes.v,mpp,linCT.v,sigth=NULL,type1fth="Mean"){

  rfdata.m <- data.m[match(intersect(rf.v,rownames(data.m)),rownames(data.m)),];
  mpp.idx <- which(celltypes.v==mpp);

  fExpr.v <- apply(rfdata.m[,mpp.idx],1,fExpr);
    
  candRF.lv <- list();
  for(i in 1:length(linCT.v)){
    lin.idx <- which(celltypes.v==linCT.v[i]);
    pval1.v <- apply(rfdata.m,1,doWTfP,mpp.idx,lin.idx);
    fdr1.v <- p.adjust(pval1.v,method="BH");
    oth <- setdiff(1:length(linCT.v),i);
    othLIN.idx <- which(celltypes.v %in% linCT.v[oth]);
    pval2.v <- apply(rfdata.m,1,doWTfP,othLIN.idx,lin.idx);
    fdr2.v <- p.adjust(pval2.v,method="BH");
    if(is.null(sigth)){
      sigth <- 0.05/nrow(rfdata.m);      
      sig.idx <- intersect(which(pval1.v < sigth),which(pval2.v < sigth)); 
    }
    else {
      sig.idx <- intersect(which(fdr1.v < sigth),which(fdr2.v < sigth));
    }
    print(paste("Found ",length(sig.idx)," candidate RFs for lineage ",linCT.v[i],sep=""));
    if(length(sig.idx)<=1){ ### demand at least 2
        print("Try relaxing significance threshold!");
    }
    fExprLin.v <- fExpr.v[sig.idx];
    if(is.null(type1fth)){
      candRF.lv[[i]] <- fExprLin.v;
    }
    else if (type1fth=="Mean"){
       print("Filtering candidate RFs by frequency of expression in MPPs");
       sel.idx <- which(fExprLin.v < mean(fExprLin.v));
       print(paste("Left with ",length(sel.idx)," candidate RFs for lineage ",linCT.v[i],sep=""));
       candRF.lv[[i]] <- fExprLin.v[sel.idx];
    }
    else {
       print("Filtering candidate RFs by frequency of expression in MPPs");
       sel.idx <- which(fExprLin.v < type1fth);
       print(paste("Left with ",length(sel.idx)," candidate RFs for lineage ",linCT.v[i],sep=""));
       candRF.lv[[i]] <- fExprLin.v[sel.idx]; 
    }
   }
  names(candRF.lv) <- linCT.v;      
  return(candRF.lv);
}

#### Auxilliary functions
#' 
#' 
doWTfP <- function(v,mpp.idx,lin.idx){
    pval <- wilcox.test(v[mpp.idx],v[lin.idx],alt="less")$p.value
    return(pval);
}

fExpr <- function(v,th=0){
    f <- length(which(v > th))/length(v);
    return(f);
}

