#' @title Model prediction 
#' @description Model projection into a RasterStack
#' 
#' @param models model class object (e.g. "glm") or list of model class objects, e.g. as returned by function \code{\link[mopa]{extractFromModel}}. 
#' @param varstack RasterStack or list of RasterStack objects with variables for projecting
#'  
#' 
#' @return RasterStack of the projected probabilities
#' @seealso \code{\link[mopa]{mopaTrain}}
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
#' ## inside different background extents
#' bg.extents <- backgroundRadios(xy = Oak_phylo2, background = bg$xy, 
#' start = 0.166, by = 0.083*10, unit = "decimal degrees")
#' TS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg.extents, 
#' exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' fittingTS <- mopaTrain(y = TS_random, x = biostack$baseline, k = 10, 
#' algorithm = "glm", weighting = TRUE, diagrams = T)
#' 
#' ## considering an unique background extent
#' RS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg$xy, 
#' exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' fittingRS <- mopaTrain(y = RS_random, x = biostack$baseline, k = 10, 
#' algorithm = "glm", weighting = TRUE)
#' 
#' modsTS <- extractFromModel(models = fittingTS, value = "model")
#' 
#' #MODEL PREDICTION
#' prdTS <- mopaPredict(models = modsTS, varstack = biostack$baseline)
#' spplot(stack(prdTS))
#' prdTS.fut <- mopaPredict(models = modsTS, varstack = biostack$future)
#' spplot(stack(prdTS.fut))
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 

mopaPredict <- function(models, varstack){
  if(class(varstack) != "list"){
    varstack <- list(varstack)
  }
  prd <- list()
  for(i in 1:length(models)){
    prd.var <- list()
    for(n in 1:length(varstack)){
      prd.var[[n]] <- mopaPredict0(models[[i]], varstack = varstack[[n]])
    }
    names(prd.var) <- names(varstack)
    if(length(prd.var)==1)  prd.var <- prd.var[[1]]
    prd[[i]] <- prd.var
  }
  names(prd) <- names(models)
  return(prd)
}

#end


#' @title Model prediction 
#' @description Model projection into a RasterStack
#' 
#' @param models model class object (e.g. "glm") or list of model class objects, e.g. as returned by function \code{\link[mopa]{extractFromModel}}. 
#' @param varstack RasterStack or list of RasterStacks of variables for projecting. If list, named lists are
#' recommended
#'  
#' 
#' @return RasterStack of the projected probabilities
#' @seealso \code{\link[mopa]{mopaTrain}}
#' 
#' @author M. Iturbide 
#' 
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' 

mopaPredict0 <- function(models, varstack){
  suppressWarnings(if(class(models) != "list") models <- list(models))
  b1 <- cbind(coordinates(varstack), rep(1, nrow(coordinates(varstack))))
  projenviro <- biomat(data = b1, varstack = varstack)[,-1]
  pro <- rep(NA, nrow(projenviro))
  projenviro2 <- projenviro[which(!is.na(projenviro[,1])),]
  projectionland <- cbind(coordinates(varstack), projenviro)
  ras <- list()
  for (i in 1:length(models)){
    alg <- models[[i]]
    algorithm <- class(alg)[1]
    if(algorithm != "MaxEnt"){
      if (algorithm == "rpart") {
            pro <- predict(alg, projenviro)
      }else if (algorithm=="cart.tree"){
            pro <- predict(alg, projenviro)
      }else if(algorithm == "ranger"){
            pro[which(!is.na(projenviro[,1]))] <- predict(alg, projenviro2)$predictions
      }else{
            pro <- predict(alg, projenviro, type="response")
      }
      pro[which(pro > 1)] <- 1
      pro[which(pro < 0)] <- 0
      ras[[i]] <- raster(SpatialPixelsDataFrame(coordinates(varstack), as.data.frame(pro)))
    }else{
      ras[[i]] <- predict(alg, varstack)
    }
  }
  names(ras) <- names(models)
  if(length(ras) == 1) ras <- ras[[1]]
  return(ras)
}

#end