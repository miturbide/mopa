#' @title Easy species distribution modelling and cross validation 
#' into backgrounds of different extent
#' @description Species distribution modelling and k-fold cross validation 
#' for a set of presence/absence data per species corresponding to a different background 
#' extent. Algorithms supported are "glm", "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' and "cart.tree" 
#' 
#' @param data Object returned by function bindPresAbs or list/s of data frames with coordinates in the first two columns and presence/absence 
#' (1=presence, 0=absence) in the third column. 
#' @param varstack RasterStack of variables for modelling
#' 
#' @param k Integer. Number of folds for cross validation. Default is 10
#' @param algorithm Any character of the following: "glm", "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' or "cart.tree"
#' @param weighting Logical for "glm" and "mars" fitting with weighted presence/absences-s. Default is FALSE.
#' @param threshold Cut value between 0 and 1 to calculate the confusion matrix. Default is 0.5.
#' @param destdir Character of the output path
#' @param projection Object of class CRS with the coordinate reference system. Default is 
#' CRS("+proj=longlat +init=epsg:4326") 
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
#'  \item{allmod }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{mod }{fitted model with partitioned data}
#'  \item{p }{cross model prediction}
#'  
#' 
#' @details This function calculates the AUC with the function \code{\link[PresenceAbsence]{auc}} from package 
#' \pkg{PresenceAbsence}. \strong{Note:} Package \pkg{SDMTools} must be detached.
#' 
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' 
#' \dontrun{
#' data(presaus)
#' data(biostackENSEMBLES)
#' ##modeling
#' modirs <- allModeling(data = presaus, varstack = biostackENSEMBLES$baseline, k = 10, "mars") 
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
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


allModeling <- function(data, 
                        varstack, 
                        k = 10, 
                        algorithm = c("glm", "svm", "maxent", "mars", "randomForest", "cart.rpart", "cart.tree"), 
                        weighting = FALSE,
                        threshold = 0.5,
                        destdir =getwd(), 
                        projection = CRS("+proj=longlat +init=epsg:4326")){
  algorithm <- match.arg(algorithm, choices = c("glm", "svm", "maxent", "mars", "randomForest", "cart.rpart", "cart.tree"))
  biostack <- varstack
  if (class(data[[1]]) != "list"){
    data <- list(data)
    names(data) <- "species"
  } 
  extents <- list()
  dirsmain <- list()
  for (i in 1:length(data)){
    sp_01 <- data[[i]]
    dirs <- list()
    if(is.null(names(sp_01))) names(sp_01) <- NA 
    extents[[i]] <- names(sp_01)
    for(j in 1:length(sp_01)){
      #print(paste("running model for species", i, "considering pseudo-absences inside the extent of", names (sp_01)[j]))
      destfile <- names(sp_01)[j]
      sp.bio <- biomat(sp_01[[j]], biostack)
      x <- kfold(k, df = sp.bio)
      xx <- leaveOneOut(x)
      # mod <- tryCatch({modelo(kdata = xx, data=sp.bio, algorithm = algorithm, weighting = weighting, threshold = threshold)},
      #                 error = function(err){xxx = list(rep(NA, k), NA, NA)})
      mod <- modelo(kdata = xx, data=sp.bio, algorithm = algorithm, weighting = weighting, threshold = threshold)
      if (length(data)==1){
        dirs[[j]] <- paste(destdir, "/", algorithm,"_", destfile, ".rda",sep="")
        save(list=c("mod"), file = dirs[[j]])
      }else{
        dirs[[j]] <- paste(destdir, "/", algorithm,"_", destfile, "_hg",names(data)[i], ".rda",sep="")
        save(list=c("mod"), file = dirs[[j]])
      }
      rm(mod, xx, x, sp.bio)
    }
    dirsmain[[i]] <- unlist(dirs)
  }
  names(dirsmain) <- names(data)
  names(extents) <- names(data)
  #collect information
  modirs <- list("dirs" = dirsmain,
                 "algorithm" = algorithm,
                 "species" = names(data),
                 "extents" = extents)
  return(modirs)
}

#end



#' @title Extract spatial background from raster
#' @description Obtain a spatial points object of the background given by a raster object.
#' 
#' @param raster Raster object (typically corresponding to a variable used for modelling species distributions).
#' @param projection Object of class CRS with the coordinate reference system. Default is 
#' CRS("+proj=longlat +init=epsg:4326") 
#' 
#'  
#' 
#' @return SpatialPoints object.
#' 
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' 
#' \dontrun{
#' data(biostack)
#' grid <- background(biostack$baseline$bio2) 
#' plot(grid)
#' }

background <- function(raster, projection = CRS("+proj=longlat +init=epsg:4326")){
  ac <- xyFromCell(raster, 1:ncell(raster))
  ex <- extract(raster, ac)
  
  # Convert to a Spatial object and define projection
  sp_grid <- SpatialPoints(ac[-which(is.na(ex)), ])
  projection(sp_grid) <- projection
  return(sp_grid)
}

#end




#' @title Bind presences and absences  
#' @description Binds presence and absence data for each background extension 
#' 
#' @param presences Data frame or list of data frames with coordinates for presence data 
#' (each row is a point)
#' @param absences Object returned by function \code{\link[mopa]{pseudoAbsences}}. 
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
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
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

#end







#' @title Load Rdata files storing definitive fitted model 
#' @description Load information from Rdata generated with function allModelling based 
#' in an index to select the definitive fitted model/s, ideally, the index returned by 
#' function indextent should be used.
#' 
#' @param models Character of the rdata files to be loaded. Including, at least,
#'  the attribute "species" (a character of names corresponding to the modelled species/populations).
#' @param slot Any character of the following: "allmod", auc", "kappa", "tss", "mod", "p"
#' 
#'  
#' 
#' @return Depending on the specified slot:
#'  \item{allmod }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{mod }{fitted model with partitioned data}
#'  \item{p }{cross model prediction}
#'  
#' 
#' @details detail.
#' 
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
#' def <- loadDefinitiveModel(models = ind, slot = "allmod")
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 


loadDefinitiveModel<-function(models, slot = c("allmod", "auc", "kappa", "tss", "mod", "p")){
  slot <- match.arg(slot, choices = c("allmod", "auc", "kappa", "tss", "mod", "p"))
  species <- attr(models, "species") 
  modelslot<-list()
  for (i in 1:length(species)){
    g <- species[i]
    mod <- get(load(models[i]))
    if (slot == "allmod"){modelslot[[i]]<-mod$allmod
    } else if (slot == "auc"){modelslot[[i]]<-mod$auc
    } else if (slot == "kappa"){modelslot[[i]]<-mod$kappa
    } else if (slot == "tss"){modelslot[[i]]<-mod$tss
    } else if (slot == "mod"){modelslot[[i]]<-mod$mod
    } else if (slot == "p"){modelslot[[i]]<-mod$p}
  } 
  names(modelslot) <- species
  suppressWarnings(
    if(length(species) == 1) modelslot <- modelslot[[1]]
  )
  return (modelslot)
}



#end







#' @title Load model validation statistics from Rdata generated with function allModelling
#' @description Load model validation statistics from Rdata generated with function allModelling
#' 
#' @param models Object with the same structure as the object returned by function allModelling.
#' @param test Any character of the following: "auc", "kappa", "tss".
#' "cart.rpart" or "cart.tree"
#' 
#'  
#' 
#' @return Matrix with values for each species and background extents.
#' 
#' #@details detail.
#' 
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
#' auc_mars <- loadTestValues(models = modirs, test = "auc") 
#'  
#' library(lattice)
#' levelplot(auc_mars, aspect = 5, 
#'           scales = list(y = list(cex = 0.8, 
#'                                  at = c(1, 49, ncol(auc_mars)), 
#'                                  labels = c(colnames(auc_mars)[1],
#'                                             colnames(auc_mars)[49],
#'                                             colnames(auc_mars)[ncol(auc_mars)]))), 
#'           at = seq(0.6, 1, 0.01), col.regions = bpy.colors,
#'           xlab = "Haplogroups", ylab = "Background extent (km)", main = "AUC")
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 



loadTestValues <-function(models, test = c("auc", "kappa", "tss")){
  test <- match.arg(test, choices = c("auc", "kappa", "tss"))
  extents <- models$extents
  extent.lengths <- lengths(extents)
  max.extents <- extents[[which(extent.lengths == max(extent.lengths))[1]]]
  algorithm <- models$algorithm
  species <- models$species
  val_mat <- matrix(NA, nrow= length(species), ncol = length(max.extents), dimnames=list(species, max.extents))
  for (i in 1:length(species)){
    message(paste("loading values for species", species[i]))
    for (k in 1:length(extents[[i]])){
      mod <- get(load(models$dirs[[i]][k]))
      if (test == "auc"){
        val_mat[i,k] <- mod$auc
      }else if (test == "kappa"){
        val_mat[i,k]<-mod$kappa
      }else if (test == "tss"){
        val_mat[i,k]<-mod$tss
      }
    } 
  }
  message("Done.")
  attr(val_mat, "dirs") <- models$dirs
  attr(val_mat, "algorithm") <- algorithm
  attr(val_mat, "species") <- species
  attr(val_mat, "extents") <- extents
  return (val_mat)
}

#end