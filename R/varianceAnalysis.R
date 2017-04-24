#' @title Variance analysis of RasterStack objects
#' @description Extract componets of lists of objects 
#' (as returned by function \code{\link[mopa]{mopaPredict}}) 
#' and perform variance analysis to obtain raster objects of the 
#' contribution of each component to the observed varaibility.
#' 
#' 
#' @param predictions listed lists of raster objects as returned by \code{\link[mopa]{mopaPredict}}
#' @param component1 Character of the names in the list (e.g.\code{predictions}) that correspond to the first component in 
#' the variance analysis.
#' @param component2 Character of the names in the list (e.g.\code{predictions}) that correspond to the second component in 
#' the variance analysis.
#' @param stick Character of the component names corresponding to the components that are not being 
#' analyzed (component 3,..). One name for each component must be provided, 
#' e.g. if there is only one component left (i.e. component 3) 
#' stick must be a single character and if there are two (i.e. component 3 and 4) stick 
#' must be a character string of length 2. 
#' 
#' @details Rasters are extracted using function \code{\link[base]{grep}}, by matching
#' names in the lists and characters in \code{componen1} and \code{componen2}. 
#' The contribution of componen1 in front component2 to the spread (uncertainty) of the projected probabilities 
#' in \code{predictions} is here assessed using a simple analysis of variance approach, where the total 
#' variance (V) can be decomposed as the summation of the variance explained by component1 (Vcomp1), 
#' component2 (Vcomp2) and the combination of the previous two (Vcomp12):
#' 
#' \eqn{V = Vcomp1 + Vcomp2 + Vcomp2}.
#' 
#' @return A list of two RasterStack objects, the first containing the global mean and standard deviation and the 
#' second containing the percentage of variance correponding to each component in the analysis 
#' (component1, component2 and components 1 and 2). 
#' 
#' 
#' @author M. Iturbide 
#' 
#' @references San Martin, D., Manzanas, R., Brands, S., Herrera, S., & Gutierrez, J.M. (2016) Reassessing 
#' Model Uncertainty for Regional Projections of Precipitation with an Ensemble of Statistical Downscaling Methods. 
#' Journal of Climate 30, 203-223.
#' 
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
#' 
#' fittingRS <- mopaTrain(y = RS_random, x = biostack$baseline, k = 10, 
#' algorithm = "glm", weighting = TRUE)
#' 
#' modsRS <- extractFromModel(models = fittingRS, value = "model")
#' 
#' #MODEL PREDICTION AND ANALYSIS OF THE VARIABILITY IN PROJECTIONS
#' prdRS.fut <- mopaPredict(models = modsRS, varstack = biostack$future)
#' component2 <- names(prdRS.fut)
#' component1 <- names(prdRS.fut$realization01$H11)
#' result <- varianceAnalysis(prdRS.fut, component1, component2, stick = "H11")
#' spplot(result$variance, col.regions = rev(get_col_regions()))
#' }
#' 
#' 
#' @export
#' @importFrom stats sd

varianceAnalysis <- function(predictions, component1, component2, stick = NULL){
  d <- depth(predictions)
  comp2 <- list()
  for(i in 1:length(component2)){
    comp2[[i]] <- extractFromPrediction(predictions, component2[i])
  }
  names(comp2) <- component2
  
  comp1 <- list()
  for(i in 1:length(component1)){
    comp1[[i]] <- extractFromPrediction(comp2, component1[i])
  }
  
  if(d > 2){
    a <- d -3
    i <- 0
    if(is.null(stick) | length(stick) != a){
      stop("Available components in prediction is more than 2 or 3, 
           set argument stick to select the component member/s that is/are kept 
           constant in the analysis")
    }
      while(a!= 0){
        a <- a-1
        i <- i+1
        comp10 <- stack(unlist(comp1))
        comp1 <- extractFromPrediction(comp10, stick[i])
      }
    }
  
  
  bothcomp <- stack(unlist(comp1))
  names(bothcomp)
  if(length(component1)*length(component2)!= nlayers(bothcomp)) stop("speciefied component names do not completely 
                                                                     match with names in predictions")
  datos <- array(NA, dim = c(length(component1)*length(component2), ncell(bothcomp)), dimnames = list(names(bothcomp)))
  for(i in 1:nlayers(bothcomp)){
    datos[i, ] <- bothcomp[[i]]@data@values
  }
  
  mediaGlobal <- apply(datos, FUN = "mean", MARGIN=2, na.rm = TRUE)
  varGlobal <- apply(datos, FUN = "sd", MARGIN=2, na.rm = TRUE)*sqrt((nrow(datos)-1)/nrow(datos))
  
  
  mediacomp1 <- matrix(data = NA, nrow = length(component1), ncol = ncol(datos))
  for(i in 1:length(component1)){
    indcomp1 <- ((i-1)*length(component2)+1):(i*length(component2))
    mediacomp1[i,] <- apply(datos[indcomp1,], FUN = "mean", MARGIN=2, na.rm = TRUE)
  }
  rownames(mediacomp1) <-component1
  
  mediacomp2 <- matrix(data = NA, nrow = length(component2), ncol = ncol(datos))
  for(i in 1:length(component2)){
    mediacomp2[i,] <- apply(datos[seq(i,nrow(datos),length(component2)),], FUN = mean, MARGIN=2, na.rm = TRUE)
  }
  dos <- apply(mediacomp2, FUN = "sd", MARGIN=2, na.rm = TRUE)*sqrt((length(component2)-1)/length(component2))
  uno <- apply(mediacomp1, FUN = "sd", MARGIN=2, na.rm = TRUE)*sqrt((length(component1)-1)/length(component1))
  
  uno.dos <- matrix(data = NA, nrow = nrow(datos), ncol = ncol(datos))
  for(i in 1:length(component1)){
    for(j in 1:length(component2)){
      uno.dos[(i-1)*length(component2)+j,] <- (datos[(i-1)*length(component2)+j,] - mediacomp2[j,] - mediacomp1[i,] + mediaGlobal)^2
    }
  }
  
  uno.dos <- apply(uno.dos, FUN = "mean", MARGIN=2, na.rm = TRUE)
  
  plot(dos^2+uno^2+uno.dos, typ = "l")
  lines(varGlobal^2, col = "red")
  sd(dos^2+uno^2+uno.dos - varGlobal, na.rm = T)
  mean(dos^2+uno^2+uno.dos - varGlobal, na.rm = T)
  
  uno100 <- uno^2 *100 / (uno^2+dos^2+uno.dos)
  dos100 <- dos^2 *100 / (uno^2+dos^2+uno.dos)
  uno.dos100 <- uno.dos *100 / (uno^2+dos^2+uno.dos)
  nan.ind <- which(uno100 == "NaN")
  dos100[nan.ind] <- 0
  uno100[nan.ind] <- 0
  uno.dos100[nan.ind] <- 0
  
  
  bothcomp[[1]]@data@values <- uno100
  bothcomp[[2]]@data@values <- dos100
  bothcomp[[3]]@data@values <- uno.dos100
  bothcomp[[4]]@data@values <- varGlobal
  bothcomp[[5]]@data@values <- mediaGlobal
  l1 <- stack(bothcomp[[5]], bothcomp[[4]])
  l2 <- stack(bothcomp[[1]], bothcomp[[2]], bothcomp[[3]])
  
  names(l1) <- c("mean", "sd")
  onename <- deparse(substitute(component1))
  twoname <- deparse(substitute(component2))
  names(l2) <- c(onename, twoname, paste(onename, "and", twoname))
  
  return(list("mean" = l1, "variance" = l2))
}


#end