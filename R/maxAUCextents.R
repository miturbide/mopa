

#' @title Index of the extents with the maximum scores of the test statistics
#' @description Given a matrix as returned by function loadTextValues
#' extracts an index of the extents at which the AUC, kappa or tss statistics obtain
#' the maximum score
#' 
#' @param testmat Object returned by function loadTestValues (matrix with extents 
#' in columns and species in rows)
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
#' ##delimit study area
#' data(Oak_phylo2)
#' data(sp_grid)
#' data(presaus)
#' data(biostack)
#' oak.extension<-boundingCoords(Oak_phylo2)
#' box.grid<-delimit(oak.extension, sp_grid, names(Oak_phylo2))
#' ##modelling
#' modirs <-allModeling(data = presaus, varstack = biostack, k = 10, "mars") 
#' ##loading#'  
#' auc_mars <-loadTestValues(data = presaus, "auc", "mars") 
#' 
#' ind<-indextent(auc_mars)
#' 
#' @export





indextent<-function(testmat){
  a<-list()
  for ( i in 1:nrow(testmat)){
    o<-which((na.omit(testmat[i,])==max(na.omit(testmat[i,])))==T)
    if (length(o) > 1){
    a[[i]]<-o[which(o == min(o))]
  } else {
    a[[i]]<-o
  }
  }
  return(unlist(a))
}
