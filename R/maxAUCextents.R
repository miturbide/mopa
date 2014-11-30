

#' @title Index of the minimum extents corresponding to the test statistics 
#' in the specified percentile
#' @description Given a matrix as returned by function loadTextValues
#' extracts an index of the extents at which the AUC, kappa or tss statistics obtain
#' the maximum score
#' 
#' @param testmat Object returned by function loadTestValues (matrix with extents 
#' in columns and species in rows)
#' @param qprob Percentile. Value between 0 and 1. Default is 0.9
#'   
#' @return Named integer. Index of the extents at which the AUC, kappa or tss 
#' statistics obtain the maximum score
#' 
#' @details detail.
#' 
#' 
#' 
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples
#' data(presaus)
#' data(biostack)
#' ##modeling
#' modirs <-allModeling(data = presaus, varstack = biostack, k = 10, "mars") 
#' ##loading#'  
#' auc_mars <-loadTestValues(data = presaus, "auc", "mars") 
#' 
#' ind<-indextent(auc_mars)
#' 
#' @export


indextent<-function(testmat, qprob = .9){
  e<-list()
  for (i in 1:nrow(testmat)){
    q<-quantile(na.omit(testmat[i,]), qprob)
    qs <- which(na.omit(testmat[i, ]) > q)  
    e[[i]]<-qs[which(qs == min(qs))]
  }
  return(unlist(e))
}
