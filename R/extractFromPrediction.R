#' @title Extract values from objects or list of objects
#' @description  Extrac objects from lists returned by function \code{\link[mopa]{mopaPredict}}.
#' 
#' 
#' @param predictions Listed lists of objects (e.g. as returned by \code{\link[mopa]{mopaPredict}})
#' @param value Character pointing to the name of component/s in the list
#' 

#' 
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' # SHORT EXAMPLE
#' destfile <- tempfile()
#' data.url <- "https://raw.githubusercontent.com/SantanderMetGroup/mopa/master/data/biostack.rda"
#' download.file(data.url, destfile)
#' load(destfile, verbose = TRUE)
#' 
#' ## Fitted models
#' data(mods)
#' ?mods
#' 
#' ## Model prediction
#' newClim <- lapply(1:4, function(x){
#' crop(biostack$future[[x]], extent(-10, 10, 35, 65))
#' })
#' names(newClim) <- names(biostack$future)[1:4]
#' prdRS.fut <- mopaPredict(models = mods, newClim = newClim)
#' names(newClim)
#' ## Extract predictions for climate projection CNRM
#' predsMPI <- extractFromPrediction(predictions = prdRS.fut, value = "CNRM")
#' spplot(predsMPI)
#' 
#' \donttest{
#' # FULL WORKED EXAMPLE
#' ## Load presence data
#' data(Oak_phylo2)
#' 
#' ## Load climate data
#' destfile <- tempfile()
#' data.url <- "https://raw.githubusercontent.com/SantanderMetGroup/mopa/master/data/biostack.rda"
#' download.file(data.url, destfile)
#' load(destfile, verbose = TRUE)
#' 
#' ## Spatial reference
#' r <- biostack$baseline[[1]]
#' ## Create background grid
#' bg <- backgroundGrid(r)
#' 
#' ## Generate pseudo-absences
#' RS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg$xy,
#'                            exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' ## Model training
#' fittedRS <- mopaTrain(y = RS_random, x = biostack$baseline,
#'                       k = 10, algorithm = "glm", weighting = TRUE)
#' ## Extract fitted models
#' mods <- extractFromModel(models = fittedRS, value = "model")
#' 
#' ## Model prediction
#' preds <- mopaPredict(models = mods, newClim = biostack$future)
#' ## Extract predictions for species climate projection MPI
#' predsMPI <- extractFromPrediction(predictions = preds, value = "MPI")
#' spplot(predsMPI)
#' }
#' 
#' @export
#' 
extractFromPrediction <- function(predictions, value){
  if(is.list(predictions))  predictions <- stack(unlist(predictions))
  rasters <- predictions
  nm <- names(rasters)
  ind <- grep(value, nm)
  r <- rasters[[ind]]
  if(length(ind) > 1) names(r) <- nm[ind]
  
  return(r)
}

#end