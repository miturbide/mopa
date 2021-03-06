

#' @title Matrix with variables for modelling
#' @description Prepares matrix with variables for modelling
#' 
#' @param data Data frame with coordinates in the first two columns and presence/absence 
#' (1=presence, 0=absence) in the third column.
#' @param varstack RasterStack of variables for modelling
#' @param projection Object of class CRS with the coordinate reference system. Default is 
#' CRS("+proj=longlat +init=epsg:4326") 
#' 
#' @return  2D matrix with the dependent variable (presence/absence) in the first column
#' and the independent variables (extracted from varstack) in the rest.
#' 
#' 
#' @details Details.
#' 
#' 
#' 
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' data(biostack)
#' data(Oak_phylo2)
#' dfp <-cbind(Oak_phylo2[[1]], "pa"= rep(1,nrow(Oak_phylo2[[1]])))
#' dfa <-cbind(Oak_phylo2[[2]], "pa"= rep(0,nrow(Oak_phylo2[[2]])))
#' df3 <-rbind(dfp, dfa)
#' str(df3)
#' mat <-biomat(df3, biostack)
#' str(mat)
#' }
#' 
#' @export
#' @keywords internal
#' @import sp
#' @import raster

biomat<-function(data, varstack, projection=CRS("+proj=longlat +init=epsg:4326")) {
  bio<-varstack
  coord.esp<-data[,1:2]
  sp.coord.esp<-SpatialPoints(coord.esp)
  proj4string(sp.coord.esp)<-projection
  proj4string(bio)<-projection
  z<-extract(bio,sp.coord.esp)
  bio.mat<-cbind(data[,3],z)
  bio.df<-as.data.frame(bio.mat)
  return(bio.df)
}
