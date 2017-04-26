#' @title Extrac objects from lists returned by function \code{\link[mopa]{mopaTrain}}
#' @description Extract values returned by function \code{\link[mopa]{mopaTrain}}
#' 
#' @param models Object returned by \code{\link[mopa]{mopaTrain}}.
#' @param value Any character of the following: "model", "auc", "kappa", "tss", "fold.models", "ObsPred"
#' 
#'  
#' 
#' @return Depending on the specified value:
#'  \item{model }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{fold.models }{fitted model with partitioned data}
#'  \item{ObsPred }{Observed and prediced (cross model prediction) values}
#'  
#' 
#' @details \code{ObsPred} allows to calculate further accuracy measures. 
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
#' ## inside different background extents
#' bg.extents <- backgroundRadios(xy = Oak_phylo2, background = bg$xy, 
#' start = 0.166, by = 0.083*10, unit = "decimal degrees")
#' TS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg.extents, 
#' exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' fittingTS <- mopaTrain(y = TS_random, x = biostack$baseline, 
#' k = 10, algorithm = "glm", weighting = TRUE, diagrams = T)
#' 
#' ## considering an unique background extent
#' RS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg$xy, 
#' exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' fittingRS <- mopaTrain(y = RS_random, x = biostack$baseline, 
#' k = 10, algorithm = "glm", weighting = TRUE)
#' 
#' modsTS <- extractFromModel(models = fittingTS, value = "model")
#' ObsPredTS <- extractFromModel(models = fittingTS, value = "ObsPred")
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export


extractFromModel <- function(models, value = c("model", "auc", "kappa", "tss", "fold.models", "ObsPred")){
  extmod <- list()  
  for(i in 1:length(models)){
    extmod1 <- list()
    for(l in 1:length(models[[i]])){
      extmod2 <- list()
      for(n in 1:length(models[[i]][[l]])){
        extmod3 <- list()
        for(h in 1:length(models[[i]][[l]][[n]])){
          extmod3[[h]] <- extractFromModel0(models = models[[i]][[l]][[n]][[h]], value = value)
        }
        names(extmod3) <- names(models[[i]][[l]][[n]])
        extmod2[[n]] <- extmod3
      }
      names(extmod2) <- names(models[[i]][[l]])
      extmod1[[l]] <- extmod2
    }
    names(extmod1) <- names(models[[i]])
    extmod[[i]] <- extmod1
  }
  names(extmod) <- names(models)
  return(extmod)
}


#end


#' @title Internal function to extrac objects from lists returned by function \code{\link[mopa]{mopaTrain}}
#' @description Internal function to extract values returned by function \code{\link[mopa]{mopaTrain}}
#' 
#' @param models Object returned by \code{\link[mopa]{mopaTrain}}.
#' @param value Any character of the following: "model", "auc", "kappa", "tss", "fold.models", "ObsPred"
#' 
#'  
#' 
#' @return Depending on the specified value:
#'  \item{model }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{fold.models }{fitted model with partitioned data}
#'  \item{ObsPred }{Observed and prediced (cross model prediction) values}
#'  
#' 
#' @details \code{ObsPred} allows to calculate further accuracy measures. 
#' 
#' 
#' 
#' @author M. Iturbide 
#' 
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 


extractFromModel0 <- function(models, value = c("model", "auc", "kappa", "tss", "fold.models", "ObsPred")){
  value <- match.arg(value, choices = c("model", "auc", "kappa", "tss", "fold.models", "ObsPred"))
  if (value == "model"){modelvalue <- models$model
    } else if (value == "auc"){modelvalue <- models$auc
    } else if (value == "kappa"){modelvalue <-models$kappa
    } else if (value == "tss"){modelvalue<-models$tss
    } else if (value == "fold.models"){modelvalue <-models$fold.models
    } else if (value == "ObsPred"){modelvalue <-models$ObsPred}
  return(modelvalue)
}

#end

