
#' @title Create background coordinates from raster object
#' @description Creates the background coordinates used for defining the study area and for generating pseudo-absences.
#' @param raster raster object with projection (\code{?raster::projection} and \code{?crs}) 
#' from which to extract the the point grid (used as background coordinates for generating pseudo-absences)
#' @param spatial.subset object of class extent (see \code{\link[raster]{extent}}) or a two
#' column data.frame (or matrix) of coordinates (xy, each row is a point).
#' 
#' @details If a data.frame, matrix or a list of the previous is passed to \code{spatial.subset} the bounding coordinates
#' are extracted to delimit the background. For example, to bound the study area to the spatial distribution of a species
#' 
#' @return A list with a SpatialPolygons object and a matrix of the background coordinates 
#' 
#' @author M. Iturbide 
#' 
#' @importFrom abind abind
#' 
#' @examples 
#' \dontrun{
#' data(Oak_phylo2)
#' data(biostack)
#' projection(biostack$baseline) <- CRS("+proj=longlat +init=epsg:4326")
#' r <- biostack$baseline[[1]]
#' # Background around a set of coordinates
#' bg.species <- backgroundGrid(r, Oak_phylo2)
#' # Background of a subdomain of the study area
#' bg.subdomain <- backgroundGrid(r, extent(c(-10, 30, 35, 65)))
#' # Background of the whole study area
#' bg <- backgroundGrid(r)
#' }

#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' @export


backgroundGrid <- function(raster, spatial.subset = NULL){
  if(is.na(projection(raster))) stop("raster has no projection arguments. Type ?raster::projection")
  if(is.null(spatial.subset)) spatial.subset <- extent(raster)
  if(class(spatial.subset) == "Extent"){
      raster <- crop(raster, spatial.subset)
  }else{
    a <- abind(boundingCoords(spatial.subset), along = 2)
    b <- c(min(a[1,]), max(a[2,]), min(a[3,]), max(a[4,]))
    raster <- crop(raster, extent(b))
  }
  ac <- xyFromCell(raster, 1:ncell(raster))
  ex <- extract(raster, ac)
  # Convert to a Spatial object and define projection
  sp_grid <- SpatialPoints(ac[-which(is.na(ex)), ])
  projection(sp_grid) <- projection(raster)
  if(class(spatial.subset) != "Extent"){
      bcoord <- boundingCoords(spatial.subset)
  }else{
      bcoord <- boundingCoords(coordinates(sp_grid))
  }
  del <- delimit(bounding.coords = bcoord, grid = sp_grid)
  return(del)
}
  
    
 