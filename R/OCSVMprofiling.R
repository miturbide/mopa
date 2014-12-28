
#' @title Environmental profiling with one-classification support vector machine
#' @description Presence-only modelling and classification of coordinates predicted as 
#' presence and absence
#' 
#' @param xy Data frame or list of data frames with coordinates (each row is a point)
#' @param varstack RasterStack of variables for modelling
#' @param bbs.grid Object derived from function \code{\link[mopa]{delimit}}. 
#' Matrix or list of matrixes of the background xy coordinates in columns. 
#' @param nu parameter needed for one-classification \code{\link[e1071]{svm}}. 
#' Default is 0.5
#'
#' @seealso \code{\link[e1071]{svm}} 
#' 
#' @return  A list with two components:
#'   \item{absence }{Matrix or list of matrixes with xy coordinates predicted 
#' as absence, each matrix correspond to a different background extent}
#' \item{presence }{Matrix or list of matrixes with xy coordinates predicted as presence, 
#' each matrix correspond to a different background extent}
#' 
#' @details This function constitutes the 
#' first step from a three-step proccess to generate pseudo-absences, and is aimed at 
#' separating from background (see \code{\link[mopa]{delimit}}) those areas (xy records) 
#' suitable for the species. 
#' 
#' 
#' 
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples
#' ##delimit study area
#' data(Oak_phylo2)
#' data(sp_grid)
#' oak.extension<-boundingCoords(Oak_phylo2)
#' box.grid<-delimit(oak.extension, sp_grid, names(Oak_phylo2))
#' ## environmental profiling
#' data(biostack)
#' unsuitable.bg <-OCSVMprofiling(xy = Oak_phylo2, varstack = biostack, 
#' bbs.grid = box.grid$bbs.grid)

#' ##plot
#' plot(unsuitable.bg$absence$H11, pch="*")
#' points(unsuitable.bg$presence$H11, pch="*", col= "pink4")
#' 
#' @export
#' 
#' @import raster
#' @import sp
#' @importFrom e1071 svm
#'  
#' 
#' 

OCSVMprofiling<-function(xy, varstack, bbs.grid, nu=.5){
  
  if (class(bbs.grid) != "list"){
    coords<-list(bbs.grid)
  }else{coords <-bbs.grid}
  
  if (class(xy) != "list"){
    pres<- list(xy)
  } else {pres <- xy}
  
  bioclim <-varstack
  absence<-list()    
  presence<-list()
  
    for(i in 1:length(pres)){
      length(absence)<-i
      coo<-coords[[i]]
      mat<- cbind(pres[[i]], rep(1, nrow(pres[[i]])))
      mat<-biomat(mat, bioclim)
      mod<-svm(mat[,-1], y=NULL, type='one-classification', nu=nu)
                  
          proj<-biomat(cbind(coo,rep(1,nrow(coo))), bioclim)
          pre<-predict(mod, proj[,-1])
          absence[[i]]<-coo[(which(pre==0)),]
          presence[[i]]<-coo[(which(pre!=0)),]

    }
   
  if (length(absence) == 1){
     absence <-absence[[1]]
     presence <-presence[[1]]
   }else{
     names(absence)<-names(pres)
     names(presence)<-names(pres)
   }
  return(list("absence"=absence, "presence"=presence))
}
