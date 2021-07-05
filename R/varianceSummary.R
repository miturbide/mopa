#' @title Summary of the variance analysis results.
#' @description Illustrates the results of the variance analysis performed with 
#' function \code{\link[mopa]{varianceAnalysis}} based on a spatial subsetting of the 
#' study area.
#' 
#' 
#' @param ... A single or multiple objects returned by function \code{\link[mopa]{varianceAnalysis}}.
#' @param component Integer indicating the number of the component analyzed in function 
#' \code{\link[mopa]{varianceAnalysis}}. Choices are 1 (Component1), 2 (Component2) or 3 (Component12).
#' @param drawBoxplot Logical. If FALSE, points of the spatial mean are plotted.
#' @param regions Either a \code{SpatialPolygons} object (see \code{\link[sp]{SpatialPolygons}}) or integer
#' of length 2 indicating the number of chunks in which to divide longitudes and latitudes. 
#' @param addWorld Logical that indicates if the lines of the world should be drawn on the map.
#' @param outline Same as in \code{\link[graphics]{boxplot}}.
#' @param parMap  List of graphical parameters affecting the output map.
#' @param parGraph  List of graphical parameters affecting the output boxplot/pointplot.
#' 
#' @details When parameter 'regions' is an integer (not a \code{SpatialPolygons} object) a grid 
#' of polygons is created using function \code{\link[raster]{union}}
#' from package \pkg{raster}. 
#' 
#' Boxes/points of the total standard deviation in each region is are also plotted (over 100)
#' with colored background.
#' 
#' The output boxplot shows the spatial spread of the results in each region. On the contrary, if paramter
#' \code{drawBoxplot = FALSE}, the output graph only shows points of the spatial mean for each region.
#' However, might be useful when multiple results are compared, as the resulting graph has a cleaner 
#' appearance.
#' 
#' @return This function returns a map of the spatial subset performed and a graph of the variance 
#'  results for each spatial subset.  Additionally, returns a matrix (or list of matrixes) with the 
#'  summary of the results.
#' 
#' @importFrom raster stack extract
#' @importFrom sp Polygon Polygons SpatialPolygons bbox
#' @importFrom graphics par segments axis points
#' @importFrom grDevices colors rgb col2rgb
#' @author M. Iturbide
#' @seealso \code{\link[mopa]{varianceAnalysis}}
#' @export
#' @examples
#' ## Load climate data
#' destfile <- tempfile()
#' data.url <- "https://raw.githubusercontent.com/SantanderMetGroup/mopa/master/data/biostack.rda"
#' download.file(data.url, destfile)
#' load(destfile, verbose = TRUE)
#' 
#' ## Fitted models
#' data(mods)
#' ?mods
#' 
#' ## Model prediction and analysis of the variability in projections
#' newClim <- lapply(1:4, function(x){
#' crop(biostack$future[[x]], extent(-10, 5, 35, 60))
#' })
#' 
#' prdRS.fut <- mopaPredict(models = mods, newClim = newClim)
#' result <- varianceAnalysis(prdRS.fut, "PA", "newClim")
#' \donttest{
#' ## Summary of the results
#' varianceSummary(result, component = 2, regions = c(4, 6))
#' 
#' ## Use a SpatialPolygons object for spatial subsetting
#' 
#' destfile <- tempfile()
#' data.url <- "https://github.com/SantanderMetGroup/visualizeR/raw/devel/data/PRUDENCEregions.rda"
#' download.file(data.url, destfile)
#' load(destfile, verbose = TRUE)
#' varianceSummary("mars" = result, component = 2, regions = PRUDENCEregions)
#' }


varianceSummary <- function(..., 
                            component = 1L, 
                            drawBoxplot = TRUE,
                            regions = c(6L, 6L),
                            addWorld = TRUE, 
                            outline = FALSE, 
                            parMap = list(mar = c(0,0,0,0)),
                            parGraph = list(mar = c(2.5,2.5,3,1))){
  obj.list <- list(...)
  if (is.null(names(obj.list))) {
    nmes <- as.list(substitute(list(...)))[-1L]
    names(obj.list) <- as.character(nmes)
  }
  for (i in 1:length(obj.list)){
    if (!any(names(obj.list[[i]]) %in% c("mean", "variance"))) stop("Input object must have the structure of the object returned by function varianceAnalysis")
    if (!any(class(obj.list[[i]]$mean) != c("RasterStack", "RasterBrick"))) stop("Input object must have the structure of the object returned by function varianceAnalysis")
    if (!any(class(obj.list[[i]]$variance) != c("RasterStack", "RasterBrick"))) stop("Input object must have the structure of the object returned by function varianceAnalysis")
    if (nlayers(obj.list[[i]]$variance) != 3) stop("Input object must have the structure of the object returned by function varianceAnalysis")
    if (nlayers(obj.list[[i]]$mean) != 2) stop("Input object must have the structure of the object returned by function varianceAnalysis")
    if (!any(c(1,2,3) == component)) stop("Choices for argument 'component' are 1, 2 or 3")
  }
  if (is.integer(regions) | is.numeric(regions)) {
    ras <- obj.list[[1]]$mean$sd
    y <- seq(bbox(ras)[2,1], bbox(ras)[2,2], length.out = regions[2] + 1)
    x <- seq(bbox(ras)[1,1], bbox(ras)[1,2], length.out = regions[1] + 1)
    pols <- list()
    for(i in 1:(length(x)-1)){
      polsi <- list()
      x1 <- x[i]; x2 <- x[i + 1]
      for(l in 1:(length(y) - 1)) {
        polsi[[l]] <- Polygons(list(Polygon(rbind(c(x1, y[l]), c(x1, y[l+1]), c(x2, y[l+1]), c(x2, y[l])))),
                               ID = paste0(i,l))
      }
      pols[[i]] <- polsi
    }
    polist <- unlist(pols, recursive = T)
    regions <- SpatialPolygons(polist)
  }
  opar <- par()
  par(mfrow = c(2,1))
  do.call("par", parMap)
  plot(regions)
  if (addWorld) {
    wd <- get(load(paste0(find.package("mopa"), "/data/wrld.rda")))
    plot(wd, add = T, border = "red")
  }
  invisible(text(coordinates(regions), labels = names(regions), cex = 0.5))
  cols <- character()
  do.call("par", parGraph)
  summaries <- list()
  for (i in 1:length(obj.list)) {
    cols[i] <- colors()[i*23]
    allmaps <- stack(unlist(obj.list[[i]]))
    uncert <- extract(allmaps$sd*100, regions)
    vares <- extract(allmaps[[component + 2]], regions)
    names(vares) <- names(regions)
    if (i == 1) {
      if (drawBoxplot) {
        suppressWarnings(boxplot(NA, axes = FALSE, ylim = c(0,100), outline = outline, xlim = c(0.8, length(regions)+.2),
                                 cex.lab = 0.65))
      } else {
        plot(NA, axes = FALSE, ylim = c(0,100), xlim = c(0.8, length(regions) + .2),
                ylab = "",
                cex.lab = 0.65, xlab = "")
      }
      segments(x0 = 1:length(regions), y0 = rep(0,length(regions)), x1 = 1:length(regions), y1 = rep(100,length(regions)), 
               lty = 1, lwd = 2, col = "gray90")
      segments(x0 = 1, y0 = seq(0,100,10), x1 = length(regions), y1 = seq(0,100,10), 
               lty = 2, lwd = 1, col = "gray90")
      axis(side = 1, at = 1:length(regions), labels = names(regions), las = 2, cex.axis = 0.7)
      axis(side = 2, at = seq(0,100,10), labels = as.character(seq(0,100,10)), cex.axis = 0.7, las = 2, cex.lab = 0.5)
    }
    if (drawBoxplot) {
      b <- boxplot(uncert, add = TRUE, col = "#00000000", outline = outline, axes = FALSE,
                               notch = T, border = cols[i], plot = F)
      segments(1:length(regions), b$stats[2,], 1:length(regions), b$stats[4,], lwd = 7, col = do.call("rgb", as.list(c(col2rgb(cols[i])/255, 0.5))))
      segments(1:length(regions), b$stats[1,], 1:length(regions), b$stats[5,], lwd = 1, col = do.call("rgb", as.list(c(col2rgb(cols[i])/255, 0.5))))
      suppressWarnings(boxplot(vares, add = TRUE, border = cols[i], col = "#00000000", outline = outline, axes = FALSE, boxwex = 0.5))
      points(1:length(regions), b$stats[3,], col = "white", pch = ".", cex = 2)
      legend.pch <- c(0 , 15)
    } else {
      points(1:length(regions), lapply(vares, FUN = "mean", na.rm = T), col = cols[i], lwd = 1, pch = 1, cex = 0.7)
      points(1:length(regions), lapply(uncert, FUN = "mean", na.rm = T), col = cols[i], lwd = 0.5, pch = 16, cex = 0.7)
      legend.pch <- c(1 , 16)
    }
    m <- suppressWarnings(rbind(lapply(vares, FUN = "mean", na.rm = T)))
    s <- suppressWarnings(rbind(lapply(vares, FUN = "sd", na.rm = T)))
    sdd <- suppressWarnings(rbind(lapply(uncert, FUN = "mean", na.rm = T)))
    summ <- rbind(m, s, sdd)
    rownames(summ) <- c("var%.spatial.mean", "var%.spatial.sd", "uncertainty.spatial.mean")
    summaries[[i]] <- summ
  }
  add_legend(5,45, legend = names(obj.list), fill = cols, horiz = TRUE, bty ='n', cex = 0.65, border = "white")
  add_legend(5,48, legend = c(paste0("Variance proportion (%) explained by comp. ", names(allmaps[[component + 2]])), 
                              "Total standard deviation (%)"), pch = legend.pch, horiz = TRUE, bty = 'n', cex = 0.65)
  suppressWarnings(par(opar))
  names(summaries) <- names(obj.list)
  return(summaries)
}


#end

#' @title Add legend to the figure returned by \code{\link[mopa]{varianceSummary}}.
#' @description Add legend to the figure returned by \code{\link[mopa]{varianceSummary}}.
#' 
#' @param ... Named arguments to be passed to function \code{\link[graphics]{legend}}.
#' @importFrom graphics legend par
#' @author M. Iturbide

add_legend <- function(...) {
  opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
              mar = c(0, 0, 0, 0), new = TRUE)
  on.exit(par(opar))
  plot(0:100, 0:100, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend(...)
}

#end
