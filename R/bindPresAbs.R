
#' @title Bind presences and absences  
#' @description Binds presence and absence data for each background extension 
#' 
#' @param presences Data frame or list of data frames with coordinates for presence data 
#' (each row is a point)
#' @param absences Object returned by function \code{\link[mopa]{PseudoAbsences}}. 
#' List/s of data frames with coordinates for absence data (each row is a point)
#' 
#' @return  List/s of matrixes with xy coordinates for presence/pseudo-absence data.
#' Each matrix correspond to a different background extent
#' 
#' 
#'
#' 
#' 
#' 
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' ##delimit study area
#' data(Oak_phylo2)
#' data(sp_grid)
#' oak.extension<-boundingCoords(Oak_phylo2)
#' box.grid<-delimit(oak.extension, sp_grid, names(Oak_phylo2))
#' ## environmental profiling
#' data(biostack)
#' unsuitable.bg <-OCSVMprofiling(xy = Oak_phylo2, varstack = biostack, 
#' bbs.grid = box.grid$bbs.grid)
#' ## sequence of 100 km between distances, from 20 km to the length of the 
#' ##half diagonal of the bounding box.
#' ext <-bgRadio(xy = Oak_phylo2, bounding.coords = oak.extension, 
#' bg.absence = unsuitable.bg$absence, start = 0.166, by = 0.083, unit = "decimal degrees")
#' ## pseudo-absence generation at random
#' pa_random <-PseudoAbsences(xy = Oak_phylo2, bg.grids = ext, 
#' exclusion.buffer = 0.083, prevalence = 0.5, kmeans = FALSE)
#' presaus <-bindPresAbs(presences = Oak_phylo2, absences = pa_random)
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutíerrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 

#' 
#' 



bindPresAbs <- function (presences, absences){
  presaus<-list()
  
  if (class(absences[[1]])=="matrix"){
    absences <- list(absences)
  } else {absences <- absences}
  
  if (class(presences)=="data.frame"){
    presences <- list(presences)
  } else {presences <- presences}
  
  
  for (i in 1:length(presences)){
    pres<-presences[[i]]
    prau<-list()
    pr<-cbind(pres, rep(1, nrow(pres)))
    names(pr)<-c("x", "y", "v")
    au<-absences[[i]]
    for (j in 1:length(au)){
      aj<-cbind(as.data.frame(au[[j]]), rep(0,nrow(au[[j]])))
      names(aj)<-names(pr)
      prau[[j]]<-rbind(pr, aj)
    }
    names(prau)<-names(au)
    presaus[[i]]<-prau
    rm(aj, pr, au, prau)
  }
  if (length(presaus) == 1){
    presaus <-presaus[[1]]
  }else{
    names(presaus)<-names(presences)
  }
  return(presaus)
}
