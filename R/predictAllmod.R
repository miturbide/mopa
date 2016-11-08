#' @title Model prediction
#' @description Model projection into a RasterStack
#' 
#' @param model model class object/s, e.g. as returned by function `loadDefinitiveModel`. 
#' @param varstack RasterStack of variables for projecting
#'  
#' 
#' @return RasterStack of the projected probabilities
#'   
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' \dontrun{
#' data(presaus)
#' data(biostackENSEMBLES)
#' ##modeling
#' modirs <-allModeling(data = presaus, varstack = biostackENSEMBLES$baseline, k = 10, "mars") 
#' ##loading
#' auc_mars <-loadTestValues(models = modirs, test = "auc") 
#' ind <- indextent(testmat = auc_mars, diagrams = F)
#' def <- loadDefinitiveModel(defExtents = ind, slot = "allmod")
#' ##predicting
#' p <- predictAllmod(def, biostackENSEMBLES$future_71_100_SMHI_ECHAM5)
#' spplot(p$H11)
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 

predictAllmod <- function(model, varstack){
  b1 <- cbind(coordinates(environment), rep(1, nrow(coordinates(environment))))
  projenviro <- biomat(data = b1, varstack = environment)[,-1]
  projectionland <- cbind(coordinates(environment), projenviro)
  ras <- list()
  for (i in 1:length(model)){
    alg <- model[[i]]
    algorithm <- class(alg)[1]
    if (algorithm == "cart.rpart") {
      pro <- predict(alg, projenviro)
    }else if (algorithm=="cart.tree"){
      pro <- predict(alg, projenviro)
    }else{
      pro <- predict(alg, projenviro, type="response")
    }
    pro[which(pro > 1)] <- 1
    pro[which(pro < 0)] <- 0
    ras[[i]] <- raster(SpatialPixelsDataFrame(coordinates(environment), as.data.frame(pro)))
  }
  names(ras) <- names(model)
  return(ras)
}
