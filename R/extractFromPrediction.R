#' @title Extrac objects from lists returned by function modelPredict
#' @description Extract values from objects or list of objects 
#' returned by function \code{\link[mopa]{mopaPredict}}
#' 
#' @param predictions Object or list of objects returned by \code{\link[mopa]{mopaPredict}}, 
#' i.e. list/s of rasters.
#' @param value Character pointing to the name of component/s in predictions
#' 

#' 
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' \dontrun{
#' data(Oak_phylo2)
#' data(biostack)
#' projection(biostack$baseline) <- CRS("+proj=longlat +init=epsg:4326")
#' r <- biostack$baseline[[1]]
#' 
#' ## Background of the whole study area
#' bg <- backgroundGrid(r)
#' 
#' 
#' ## considering an unique background extent
#' RS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg$xy, 
#' realizations = 10,
#' exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' fittingRS <- mopaTrain(y = RS_random, x = biostack$baseline, k = 10, 
#' algorithm = "glm", weighting = TRUE)
#' 
#' modsRS <- extractFromModel(models = fittingRS, value = "model")
#' 
#' #MODEL PREDICTION
#' prdTS.fut <- mopaPredict(models = modsTS, varstack = biostack$future)
#' prdTS.fut.sub <- extractFromPrediction(prdTS.fut, "realization01")
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 
extractFromPrediction <- function(predictions, value){
  if(class(predictions) == "list")  predictions <- stack(unlist(predictions))
  rasters <- predictions
  nm <- names(rasters)
  ind <- grep(value, nm)
  r <- rasters[[ind]]
  if(length(ind) > 1) names(r) <- nm[ind]
  
  return(r)
}

#end