#' @title Model prediction 
#' @description Model projection into a RasterStack
#' 
#' @param models model class object (e.g. "glm") or list of model class objects, e.g. as returned by function \code{\link[mopa]{extractFromModel}}. 
#' @param newClim RasterStack or list of RasterStack objects with variables for projecting
#'  
#' 
#' @return RasterStack of the projected probabilities
#' @seealso \code{\link[mopa]{mopaTrain}}
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' data(Oak_phylo2)
#' data(biostack)
#' projection(biostack$baseline) <- CRS("+proj=longlat +init=epsg:4326")
#' r <- biostack$baseline[[1]]
#' 
#' ## Background of the whole study area
#' bg <- backgroundGrid(r)
#' 
#' ## inside different background extents
#' bg.extents <- backgroundRadius(xy = Oak_phylo2, background = bg$xy, 
#' start = 0.166, by = 0.083*10, unit = "decimal degrees")
#' TS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg.extents, 
#' exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' fittingTS <- mopaTrain(y = TS_random, x = biostack$baseline, k = 10, 
#' algorithm = "glm", weighting = TRUE, diagrams = TRUE)
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
#' prdTS <- mopaPredict(models = modsTS, newClim = biostack$baseline)
#' spplot(stack(prdTS$H11$PA01$glm$baseClim01))
#' prdTS.fut <- mopaPredict(models = modsTS, newClim = biostack$future)
#' spplot(stack(prdTS.fut$H11$PA01$glm$baseClim01$MPI))
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export

mopaPredict <- function(models, newClim){
  if(class(newClim) != "list"){
    newClim <- list(newClim)
  }
  # ml <- depth(models)-1
  # repk <- 5-ml
  # while(repk != 0){
  #   repk <- repk - 1
  #   models <- list(models)
  # }
  pred <- list()
  for(l in 1:length(models)){
    prd <- list()
    for(i in 1:length(models[[l]])){
      prd0 <- list()
      for(k in 1:length(models[[l]][[i]])){
        prd1 <- list()
        for(h in 1:length(models[[l]][[i]][[k]])){
          prd.var <- list()
          nm <- character()
          for(n in 1:length(newClim)){
            prd.var[[n]] <- mopaPredict0(models[[l]][[i]][[k]][[h]], newClim = newClim[[n]])
            if(n < 10){
              nm[n] <- paste0("0", n)
            }else{
              nm[n] <- as.character(n)
            }
          }
          if(is.null(names(newClim))) names(newClim) <- paste0("newClim", nm)
          names(prd.var) <- names(newClim)
          # if(length(prd.var)==1)  prd.var <- prd.var[[1]]
          prd1[[h]] <- prd.var
        }
        names(prd1) <- names(models[[l]][[i]][[k]])
        prd0[[k]] <- prd1
      }
      names(prd0) <- names(models[[l]][[i]])
      # if(length(prd0)==1)  prd0 <- prd0[[1]]
      prd[[i]] <- prd0
    }
    names(prd) <- names(models[[l]])
    # if(length(prd)==1)  prd <- prd[[1]]
    pred[[l]] <- prd
  }
  names(pred) <- names(models)
  # if(length(pred) == 1) pred <- pred[[1]]
  return(pred)
}



#end


#' @title Internal function for model prediction 
#' @description Internal function for model projection into a RasterStack
#' 
#' @param models model class object (e.g. "glm") or list of model class objects, e.g. as returned by function \code{\link[mopa]{extractFromModel}}. 
#' @param newClim RasterStack or list of RasterStacks of variables for projecting. If list, named lists are
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

mopaPredict0 <- function(models, newClim){
  suppressWarnings(if(class(models) != "list") models <- list(models))
  b1 <- cbind(coordinates(newClim), rep(1, nrow(coordinates(newClim))))
  projenviro <- biomat(data = b1, varstack = newClim)[,-1]
  pro <- rep(NA, nrow(projenviro))
  projenviro2 <- projenviro[which(!is.na(projenviro[,1])),]
  projectionland <- cbind(coordinates(newClim), projenviro)
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
      ras[[i]] <- raster(SpatialPixelsDataFrame(coordinates(newClim), as.data.frame(pro)))
    }else{
      ras[[i]] <- predict(alg, newClim)
    }
  }
  names(ras) <- names(models)
  if(length(ras) == 1) ras <- ras[[1]]
  return(ras)
}

#end


#' @title Level depth in a list
#' @description Level depth in a list 
#' 
#' @param this list
#'  
#' @return number of nesting lists
#' @author M. Iturbide 

depth <- function(this){
  that <- this
  i <- 0
  while(is.list(that)){
    i <- i + 1
    that <- that[[1]]
  }
  return(i+1)
}
#end