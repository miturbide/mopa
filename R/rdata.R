#' @title World map
#' @description A dataset of class SpatialPolygonsDataFrame of the World map
#' @name wrld
#' @docType data
#' @author M. Iturbide 
#' @keywords map
NULL


#' @title Oak distribution
#' @description A dataset consisting of a list with two
#' data frames with xy coordinates corresponding to the distribution of two 
#' phylogenetic groups of oaks (H11 and H1)
#' @name Oak_phylo2
#' @docType data
#' @keywords xy records
#' @source  Oak_phylo2 is a modified subset of the \strong{Quercus sp Europe Petit 2002} 
#' database (Petit et al., 2002b), which is available in the \strong{Georeferenced Database 
#' of Genetic Diversity} or \strong{(GD)^2}. 
#' @references Petit, R. J. \emph{et al} 2002. Chloroplast DNA variation in european 
#' white oaks: Phylogeography and patterns of diversity based on data from over 2600 populations. 
#' Forest Ecology and Management 156 (1-3), 5-26.
#' 
#' Evolution of Trees and Forest Communities \cr 
#' Ten years of the EVOLTREE network \cr
#' Evoltree E-Lab - An information system for forest genetics \cr
#' ISBN: 978-2-9519296-3-9
NULL

#' @title Quercus pubsencens distribution
#' @description A data frame with xy coordinates of Quercus pubesnces distribution
#' @name Q_pubescens
#' @docType data
#' @keywords xy records
#' @source  Q_pubescens is a modified subset of occurrences obtained from \strong{GBIF.org} 
#' @references GBIF.org (14th March 2017) GBIF Occurrence Download http://doi.org/10.15468/dl.4ss6vr
NULL


#' @title Fitted models 
#' @description List of fitted models as returned by functions \code{\link{mopaTrain}} and 
#' \code{\link{extractFromModel}}.
#' @name mods
#' @docType data
#' @keywords xy records
#' @source  # RS_random is the result of running the following code:
#' 
#' data(Oak_phylo2)
#' 
#' presences <- Oak_phylo2$H11
#' 
#' destfile <- tempfile()
#' 
#' data.url <- "https://raw.githubusercontent.com/SantanderMetGroup/mopa/master/data/biostack.rda"
#' 
#' download.file(data.url, destfile)
#' 
#' load(destfile, verbose = TRUE)
#' 
#' r <- biostack$baseline[[1]]
#' 
#' ## Background of the whole study area
#' bg <- backgroundGrid(r)
#' 
#' ## Considering an unique background extent
#' 
#' RS_random <-pseudoAbsences(xy = presences, background = bg$xy, 
#'                            realizations = 5,
#'                            exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' 
#' fittedModels <- mopaTrain(y = RS_random, x = biostack$baseline, k = 10, 
#'                     algorithm = "mars")
#'                     
#' mods <- extractFromModel(models = fittedModels, value = "model")
NULL
