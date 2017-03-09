#' @title Easy species distribution modeling and cross validation
#' @description Species distribution modeling and k-fold cross validation 
#' for a set of presence/absence data per species, also considering different background 
#' extents (optional). Algorithms supported are "glm", "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' and "cart.tree" 
#' 
#' @param x Object returned by function \code{link[mopa]{pseudoAbsences}} or list/s of data frames with coordinates
#'  in the first two columns and presence/absence (1=presence, 0=absence) in the third column. 
#' @param y RasterStack of variables for modelling
#' 
#' @param k Integer. Number of folds for cross validation. Default is 10
#' @param algorithm Any character of the following: "glm", "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' or "cart.tree"
#' @param weighting Logical for "glm", "mars" and "randomForest" fitting with weighted presence/absences-s. Default is FALSE.
#' @param threshold Cut value between 0 and 1 to calculate the confusion matrix. Default is NULL (see Details).
#' @param diagrams logical. Only applied if \code{x} contains data for different background extents 
#' (see \code{\link[mopa]{backgroundRadios}} and \code{\link[mopa]{pseudoAbsences}}). Should diagrams of 
#' AUC extent fitting be printed? default is FALSE. 
#' 
#'  
#' 
#' @return Named Rdata objects are stored in the 
#' specified path. Each Object is given a name indicating the algorithm, background 
#' extent, and species in this order (if a single species is provided no name is given 
#' for the species). The object returned by the function is a list of characters with listed files, 
#' used algorithm, species names and background extents.
#' Each Rdata consists of a list with six components:
#' 
#'  \item{model }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{fold.models }{fitted model with partitioned data}
#'  \item{ObsPred }{cross model prediction}
#'  
#' 
#' @details This function calculates the AUC with the function \code{\link[PresenceAbsence]{auc}} from package 
#' \pkg{PresenceAbsence}. \strong{Note:} Package \pkg{SDMTools} must be detached. If \code{threshold} is not specified the value
#' that maximisez the TSS (true skill statistic) is used to calculate the accuracy.
#' 
#' @seealso \code{\link[mopa]{mopaPredict}}, \code{\link[mopa]{pseudoAbsences}}, \code{\link[mopa]{backgroundGrid}}, 
#' \code{\link[mopa]{OCSVMprofiling}}, \code{\link[mopa]{backgroundRadios}}
#' @author M. Iturbide 
#' 
#' @examples
#' 
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
#' fittingTS <- mopaFitting(y = TS_random, x = biostack$baseline, k = 10, 
#' algorithm = "glm", weighting = TRUE, diagrams = T)
#' 
#' ## considering an unique background extent
#' RS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg$xy,
#'  exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' fittingRS <- mopaFitting(y = RS_random, x = biostack$baseline, 
#' k = 10, algorithm = "glm", weighting = TRUE)
#' 
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 
#' @import raster
#' @import sp
#' @import PresenceAbsence
#' @importFrom sampling strata  
#' @importFrom e1071 best.svm 
#' @importFrom dismo maxent
#' @importFrom earth earth
#' @importFrom rpart rpart
#' @importFrom tree tree
#' @importFrom randomForest randomForest
#' 


mopaFitting <- function(y, 
                        x, 
                        k = 10, 
                        algorithm = c("glm", "svm", "maxent", "mars", "randomForest", "cart.rpart", "cart.tree"), 
                        weighting = FALSE,
                        threshold = NULL,
                        diagrams = FALSE){
  algorithm <- match.arg(algorithm, choices = c("glm", "svm", "maxent", "mars", "randomForest", "cart.rpart", "cart.tree"))
  data <- y
  biostack <- x
  if (class(data[[1]]) != "list"){
    data <- list(data)
    names(data) <- "species"
  }
  extents <- list()
  aucmat <- array(dim = c(length(data), max(unlist(lapply(data, length)))), dimnames = list(c(names(data))))
  for (i in 1:length(data)){
    sp_01 <- data[[i]]
    if(is.null(names(sp_01))) names(sp_01) <- NA 
    extents[[i]] <- names(sp_01)
    for(j in 1:length(sp_01)){
      #print(paste("running model for species", i, "considering pseudo-absences inside the extent of", names (sp_01)[j]))
      sp.bio <- biomat(sp_01[[j]], biostack)
      x <- kfold(k, df = sp.bio)
      xx <- leaveOneOut(x)
      # mod <- tryCatch({modelo(kdata = xx, data=sp.bio, algorithm = algorithm, weighting = weighting, threshold = threshold)},
      #                 error = function(err){xxx = list(rep(NA, k), NA, NA)})
      aucmat[i, j] <- modelo(kdata = xx, data=sp.bio, algorithm = algorithm, weighting = weighting, threshold = threshold)$auc
      rm(xx, x, sp.bio)
    }
  }
  extents2 <- as.integer(unique(sub(unlist(extents), pattern = "km", replacement = "")))
  ind <- rep(1, length(data))
  if(ncol(aucmat) > 1)  ind <- AUCextentFit(aucmat, extents = extents2, diagrams = diagrams)
  mod <- list()
  for (i in 1:length(data)){
    sp_01 <- data[[i]]
    sp.bio <- biomat(sp_01[[ind[i]]], biostack)
    x <- kfold(k, df = sp.bio)
    xx <- leaveOneOut(x)
    # mod <- tryCatch({modelo(kdata = xx, data=sp.bio, algorithm = algorithm, weighting = weighting, threshold = threshold)},
    #                 error = function(err){xxx = list(rep(NA, k), NA, NA)})
    mod[[i]] <- modelo(kdata = xx, data=sp.bio, algorithm = algorithm, weighting = weighting, threshold = threshold)
    mod[[i]]$extent <- extents[[i]][ind[i]]
    rm(xx, x, sp.bio)
  }
  names(mod) <- names(data)
  return(mod)
}
