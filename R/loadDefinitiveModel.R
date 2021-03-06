


#' @title Load Rdata files storing definitive fitted model 
#' @description Load information from Rdata generated with function allModelling based 
#' in an index to select the definitive fitted model/s, ideally, the index returned by 
#' function indextent should be used.
#' 
#' @param data Object with the same structure as the object returned by function bindPresAbs. 
#' @param extents Named integer returned by function indextent (a named index 
#' corresponding to the definitive extents to be considered)
#' @param slot Any character of the following: "allmod", auc", "kappa", "tss", "mod", "p"
#' @param algorithm Any character of the following: "glm", "svm", "maxent", "mars", "randomForest", 
#' "cart.rpart" or "cart.tree"
#' @param sourcedir Character of the path where Rdata objects are 
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
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' data(presaus)
#' data(biostack)
#' ##modeling
#' modirs <-allModeling(data = presaus, varstack = biostack, k = 10, "mars") 
#' ##loading#'  
#' auc_mars <-loadTestValues(data = presaus, "auc", "mars") 
#' ind <- indextent(testmat = auc_mars, diagrams = TRUE)
#' 
#' def <-loadDefinitiveModel(data = presaus, extents = ind, slot = "allmod", algorithm = "mars")
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export


loadDefinitiveModel<-function(data, extents, slot=c("allmod", "auc", "kappa", "tss", "mod", "p"),
                    algorithm = c("glm","svm","maxent","mars","randomForest","cart.rpart",
                                  "cart.tree"), sourcedir = getwd()){
  
  if (class(data[[1]]) != "list"){
    dat<-list(data)
  }else{dat<-data}
  
  modelslot<-list()
  for (i in 1:length(data)){
    g<-names(data)[i]
    if (class(data[[1]]) != "list"){
    load(paste(sourcedir,"/", algorithm, "_bg", names(extents)[i],".Rdata",sep=""))
    } else {
    load(paste(sourcedir,"/", algorithm, "_bg", names(extents)[i],"_hg",g,".Rdata",sep=""))
    }
      if (slot == "allmod"){modelslot[[i]]<-mod$allmod
      } else if (slot == "auc"){modelslot[[i]]<-mod$auc
      } else if (slot == "kappa"){modelslot[[i]]<-mod$kappa
      } else if (slot == "tss"){modelslot[[i]]<-mod$tss
      } else if (slot == "mod"){modelslot[[i]]<-mod$mod
      } else if (slot == "p"){modelslot[[i]]<-mod$p}
  } 
  names(modelslot)<-names(data)
  return (modelslot)
}
