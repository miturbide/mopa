
#' @title Load model validation statistics from Rdata generated with function allModelling
#' @description Load model validation statistics from Rdata generated with function allModelling
#' 
#' @param data Object with the same structure as the object returned by function bindPresAbs. 
#' @param test Any character of the following: "auc", "kappa", "tss".
#' @param algorithm Any character of the following: "glm", "svm", "maxent", "mars", "randomForest", 
#' "cart.rpart" or "cart.tree"
#' @param sourcedir Character of the path where Rdata objects are 
#' 
#'  
#' 
#' @return Matrix with values for each species and background extents.
#' 
#' @details detail.
#' 
#' 
#' 
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples
#' \donttest{
#' ##delimit study area
#' data(Oak_phylo2)
#' data(sp_grid)
#' data(presaus)
#' data(biostack)
#' oak.extension<-boundingCoords(Oak_phylo2)
#' box.grid<-delimit(oak.extension, sp_grid, names(Oak_phylo2))
#' ##modelling
#' modirs <-allModeling(data = presaus, varstack = biostack, k = 10, "mars") 
#' ##loading
#' auc_mars <-loadTestValues(data = presaus, test = "auc", algorithm = "mars")
#'  
#' library(lattice)
#' levelplot(auc_mars ,aspect=5 ,at =seq(0.5,1,0.01), col.regions=bpy.colors,
#' xlab="Haplogroups" ,ylab="Background extent (km)", main = "AUC")
#' }
#' 
#' @export



loadTestValues <-function(data, test = c("auc", "kappa", "tss"), algorithm = c("glm","svm","maxent","mars","randomForest","cart.rpart","cart.tree"), sourcedir=getwd()){
  
  if (class(data[[1]]) != "list"){
    dat<-list(data)
  }else{dat<-data}
  
  l<-numeric()
  for( j in 1:length(dat)){
    l[j]<-length(dat[[j]])
  }
  
  cols<-max(l)
  maxind<-which(l==max(l))[1]
  extents<-names(dat[[maxind]])
  val_mat<-matrix(NA, nrow=length(dat), ncol=cols, dimnames=list(names(dat),extents))
  
  for (i in 1:length(dat)){
    g<-names(dat)[i]
    bgs<-length(dat[[i]])
    r<-names(dat[[i]])
    print(paste("loading values for species", i))
    for (k in 1:bgs){
      b<-k
      if (class(data[[1]]) != "list"){
      load(paste(sourcedir,"/", algorithm, "_bg", r[k],".Rdata",sep=""))
      }else {
        load(paste(sourcedir,"/", algorithm, "_bg", r[k],"_hg",g,".Rdata",sep=""))
      }
      
      if (test == "auc"){val_mat[i,k]<-mod$auc
      } else if (test == "kappa"){val_mat[i,k]<-mod$kappa
      } else if (test == "tss"){val_mat[i,k]<-mod$tss}
    } 
    
  }
  return (val_mat)
}