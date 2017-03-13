#' @title Easy species distribution modeling and cross validation
#' @description Species distribution modeling and k-fold cross validation 
#' for a set of presence/absence data per species, also considering different background 
#' extents (optional). Algorithms supported are "glm", "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' and "cart.tree" 
#' 
#' @param y Object returned by function \code{\link[mopa]{pseudoAbsences}} or list/s of data frames with coordinates
#'  in the first two columns and presence/absence (1=presence, 0=absence) in the third column. 
#' @param x RasterStack of variables for modelling
#' 
#' @param k Integer. Number of folds for cross validation. Default is 10
#' @param algorithm Any character of the following: "glm", "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' or "cart.tree"
#' @param weighting Logical for "glm", "mars" and "randomForest" fitting with weighted presence/absences-s. 
#' Default is FALSE.
#' The processing time is considerably increased if weighting option is selected when applying 
#' "mars" algorithm (see \code{\link[earth]{earth}}.
#' @param threshold Cut value between 0 and 1 to calculate the confusion matrix. Default is NULL (see Details).
#' @param diagrams Logical. Only applied if \code{x} contains data for different background extents 
#' (see \code{\link[mopa]{backgroundRadios}} and \code{\link[mopa]{pseudoAbsences}}). Should diagrams of 
#' AUC extent fitting be printed? default is FALSE. 
#' 
#'  
#' 
#' @return A list of six components is returned for each species in \code{x}:
#' 
#'  \item{model }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{fold.models }{fitted models of each data partition for cross validation}
#'  \item{ObsPred }{cross model prediction (e.g. for further assessment of model accuracy)}
#'  
#' 
#' @details This function calculates the AUC with the function \code{\link[PresenceAbsence]{auc}} from package 
#' \pkg{PresenceAbsence}. \strong{Note:} Package \pkg{SDMTools} must be detached. If \code{threshold} is not specified the value
#' that maximisez the TSS (true skill statistic) is used to calculate the accuracy.
#' 
#' 
#' \code{mopaTrain} uses the algorithm implementations of the following R packages: 
#' \pkg{earth} for "mars"; \pkg{ranger} for "randomForest"; 
#' \pkg{dismo} for "maxent"; \pkg{rpart} for "cart.rpart";
#' \pkg{e1071} for "svm"; \pkg{tree} for "cart.tree";
#' \pkg{stats} for "glm".
#' 
#' If \code{y} contains data for different background extents (see \code{\link[mopa]{backgroundRadios}} and
#' \code{\link[mopa]{pseudoAbsences}}), \code{mopaFitting} performs the species distribution modelling for 
#' each different background extent, and fits obtained AUCs (corresponding to different background extents) 
#' to three non linear models (Michaelis-Menten, exponential2 and exponential3). 
#' The model that scores the lowest error is automatically selected to extract the Vm coefficient (equation 1 in 
#' Iturbide et al., 2015). Then, the minimum extent at which the AUC surpasses the Vm value is selected 
#' as the threshold extent (see Figure 3 in Iturbide et al., 2015), being the corresponding fitted SDM the 
#' one returned by \code{mopaFitting}. If argument \code{diagrams} is set to TRUE, A fitted model plot 
#' (as in Fig. 3 in Iturbide et al., 2015) is printed in the plotting environment.
#' 
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
#' fittingTS <- mopaTrain(y = TS_random, x = biostack$baseline, k = 10, 
#' algorithm = "glm", weighting = TRUE, diagrams = T)
#' 
#' ## considering an unique background extent
#' RS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg$xy,
#'  exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' glmRS <- mopaTrain(y = RS_random, x = biostack$baseline, 
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
#' @importFrom sampling strata  


mopaTrain <- function(y, 
                        x, 
                        k = 10, 
                        algorithm = c("glm", "svm", "maxent", "mars", "randomForest", "cart.rpart", "cart.tree"), 
                        weighting = FALSE,
                        threshold = NULL,
                        diagrams = FALSE){
  mfit <- list()
  for(i in 1:length(y)){
    message("[", Sys.time(), "] realization ", i)
    mfit[[i]] <- mopaTrain0(y = y[[i]], 
                              x = x, 
                              k = k, 
                              algorithm = algorithm, 
                              weighting = weighting,
                              threshold = threshold,
                              diagrams = diagrams)
  }
  names(mfit) <- names(y)
  return(mfit)
}


#end

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
#' @return A list of six components is returned for each species in \code{x}:
#' 
#'  \item{model }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{fold.models }{fitted model with partitioned data}
#'  \item{ObsPred }{cross model prediction}
#'  
#' 
#' 
#' 
#' 
#' @seealso \code{\link[mopa]{mopaPredict}}, \code{\link[mopa]{pseudoAbsences}}, \code{\link[mopa]{backgroundGrid}}, 
#' \code{\link[mopa]{OCSVMprofiling}}, \code{\link[mopa]{backgroundRadios}}
#' 
#' @author M. Iturbide 
#' 
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' 
#' @import raster
#' @import sp
#' @importFrom sampling strata  



mopaTrain0 <- function(y, 
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
    names(data) <- "unnamed"
  }
  extents <- array(dim = c(length(data), max(unlist(lapply(data, length)))), dimnames = list(c(names(data))))
  aucmat <- array(dim = c(length(data), max(unlist(lapply(data, length)))), dimnames = list(c(names(data))))
  for (i in 1:length(data)){
    sp_01 <- data[[i]]
    if(is.null(names(sp_01))) names(sp_01) <- "maximum extent"
    extnames <- names(sp_01)
    extents[i, 1:length(extnames)] <- as.integer(sub(unlist(extnames), pattern = "km", replacement = ""))
    for(j in 1:length(sp_01)){
      #print(paste("running model for species", i, "considering pseudo-absences inside the extent of", names (sp_01)[j]))
      sp.bio <- mopa:::biomat(sp_01[[j]], biostack)
      x <- kfold(k, df = sp.bio)
      xx <- leaveOneOut(x)
      # mod <- tryCatch({modelo(kdata = xx, data=sp.bio, algorithm = algorithm, weighting = weighting, threshold = threshold)},
      #                 error = function(err){xxx = list(rep(NA, k), NA, NA)})
      aucmat[i, j] <- modelo(kdata = xx, data=sp.bio, algorithm = algorithm, weighting = weighting, threshold = threshold)$auc
      rm(xx, x, sp.bio)
    }
  }
  ind <- rep(1, length(data))
  if(ncol(aucmat) > 1)  ind <- AUCextentFit(aucmat, extents = extents, diagrams = diagrams)
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

#end




#' @title Species distribution modelling and cross validation
#' @description Species distribution modelling with k-fold cross validation. 
#' Algorithms supported are "glm", "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' and "cart.tree" 
#' @param kdata Object returned by function leaveOneOut
#' @param data Object returned by function biomat. 2D matrix with the dependent variable 
#' (presence/absence) in the first column and the independent variables in the rest 
#' (extracted from varstack) 
#' @param algorithm Any character of the following: \code{"glm"}, "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' or "cart.tree"
#' @param weighting Logical for model fitting with weighted presence/absences-s. Applicable for algorithms "glm", "mars", 
#' "randomForest" and "cart.rpart". Default is FALSE.
#' @param threshold Cut value between 0 and 1 to calculate the confusion matrix. Default is 0.5.
#'
#'  
#' 
#' @return  A list with six components:
#'  \item{model }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{fold.models }{fitted model with partitioned data}
#'  \item{ObsPred }{cross model prediction}
#' 
#' @details This function calculates the AUC with the function "auc" from package 
#' "PresenceAbsence". Package SDMTools must be detached.
#' 
#' 
#' 
#' @author M. Iturbide 
#' 
#'
#' @keywords internal
#' @import PresenceAbsence
#' @importFrom e1071 svm
#' @importFrom sampling strata  
#' @importFrom e1071 best.svm 
#' @importFrom dismo maxent
#' @importFrom earth earth
#' @importFrom rpart rpart
#' @importFrom tree tree
#' @importFrom randomForest randomForest tuneRF
#' @importFrom ranger ranger
#' @importFrom stats glm binomial na.omit 



modelo <- function(kdata, data, algorithm = c("glm","svm","maxent","mars","randomForest","cart.rpart","cart.tree"), weighting = FALSE, threshold = 0.5){
  mod <- list()
  pmod <- list()
  algorithm <- as.character(algorithm)
  orig<-list()
  alltrn<-na.omit(data)
  colnames(alltrn)<-c("V1", colnames(alltrn)[2:ncol(alltrn)])
  if(weighting){
    Force.weights <- TRUE
    w.factor <- round(length(which(alltrn$V1 == 0))/length(which(alltrn$V1 == 1))*10)
    all.weights <- alltrn$V1*w.factor
    all.weights[(which(alltrn$V1 == 0))] <- 10
    all.weights.rf <- rep(length(which(alltrn$V1 == 1)), 2)
  }else{
    Force.weights <- FALSE
    all.weights <- rep(10, length(alltrn$V1))
    all.weights.rf <- nrow(alltrn)
  }
  for(i in 1:length(kdata)) {
    length(mod) <- i
    length(pmod) <- i
    dftst <- na.omit(kdata[[i]]$test)
    colnames(dftst)<-c("V1", colnames(dftst)[2:ncol(dftst)])   
    dftrn <- na.omit(kdata[[i]]$train)
    colnames(dftrn)<-c("V1", colnames(dftrn)[2:ncol(dftrn)])  
    orig[[i]]<-na.omit(dftst[,1])
    if(weighting){
      k.w.factor <- round(length(which(dftrn$V1 == 0))/length(which(dftrn$V1 == 1))*10)
      weights <- dftrn$V1*k.w.factor
      weights[(which(dftrn$V1 == 0))] <- 10
      pres <- length(which(dftrn$V1 == 1))
      aus <- length(which(dftrn$V1 == 0))
      if(pres>aus){
        weights.rf <- rep(aus, 2)
      }else{
        weights.rf <- rep(pres, 2) 
      }
    }else{
      weights <- rep(10, length(dftrn$V1))
      weights.rf <- nrow(dftrn)
    }
    # Training
    if(algorithm == "glm") {
      mod[[i]] <- glm(V1 ~., data=dftrn, family=binomial(link="logit"), weights = weights)
    } else if(algorithm == "svm"){
      mod[[i]] <- best.svm(V1 ~., data=dftrn)
    } else if(algorithm == "maxent"){
      mod[[i]] <- maxent(x = dftrn[ ,-1], p=dftrn[,1])
    } else if(algorithm == "mars"){
      mod[[i]] <- earth(x = dftrn[ ,-1], y = dftrn[,1], weights = weights, Force.weights = Force.weights)
    } else if(algorithm=="cart.rpart"){
      mod[[i]] <- rpart(V1 ~., data=dftrn, weights = weights)
    } else if(algorithm == "cart.tree"){
      mod[[i]]<- tree(V1 ~., data=dftrn, weights = weights)
    } else if(algorithm == "randomForest"){
      # dftrn$V1 <- as.factor(dftrn$V1)
      # mod[[i]]<- randomForest(V1 ~., data=dftrn, strata = dftrn$V1, sampsize = weights.rf)
      mtry1 <- tuneRF(dftrn[,-1], dftrn[,1], plot = F)
      mtry <- mtry1[which(mtry1[,2] == min(mtry1[,2])),1]
      mod[[i]]<- ranger(V1 ~., data=dftrn, case.weights = weights, num.trees = 300, mtry = mtry)
      
    }
    # Predictions
    if (algorithm == "cart.rpart") {
      pmod[[i]] <- predict(mod[[i]], dftst[,-1])
    }else if (algorithm=="cart.tree"){
      pmod[[i]] <- predict(mod[[i]], dftst[,-1])
    }else if(algorithm == "randomForest"){
      pmod[[i]] <- predict(mod[[i]], dftst[ ,-1])$predictions
    }else{
      pmod[[i]] <- predict(mod[[i]], dftst[ ,-1], type="response")
    }
    
  }
  
  
  # allTraining
  if(algorithm == "glm") {
    allmod <- glm(V1 ~., data=alltrn, family=binomial(link="logit"), weights = all.weights)
  } else if(algorithm == "svm"){
    allmod <- best.svm(V1 ~., data=alltrn)
  } else if(algorithm == "maxent"){
    allmod <- maxent(x = alltrn[ ,-1], p=alltrn[,1])
  } else if(algorithm == "mars"){
    allmod <- earth(x = alltrn[ ,-1], y = alltrn[,1], weights = all.weights)
  } else if(algorithm=="cart.rpart"){
    allmod <- rpart(V1 ~., data=alltrn, weights = all.weights)
  } else if(algorithm == "cart.tree"){
    allmod <- tree(V1 ~., data=alltrn, weights = all.weights)
  } else if(algorithm == "randomForest"){
    # alltrn$V1 <- as.factor(alltrn$V1)
    # allmod <- randomForest(V1 ~., data=alltrn, strata = alltrn$V1, sampsize = all.weights.rf)
    mtry1 <- tuneRF(alltrn[,-1], alltrn[,1], plot = F)
    mtry <- mtry1[which(mtry1[,2] == min(mtry1[,2])),1]
    allmod <- ranger(V1 ~., data=alltrn, case.weights = all.weights, num.trees = 300)
    
  }
  
  
  ori <- unlist(orig)
  p <- unlist(pmod)
  dat <- as.data.frame(cbind("id" = as.integer(c(1:length(ori))), "obs" = ori, "pred" = p))
  if(is.null(threshold)) threshold <- kappaRepet(Obs = dat$obs, Fit = dat$pred, TSS = T)$CutOff 
  p.auc <- p
  p.auc[which(p>=threshold)] <- 1
  p.auc[which(p.auc!=1)] <- 0
  dat.auc <- as.data.frame(cbind(as.integer(c(1:length(ori))), ori, p.auc))
  mod.auc <- auc(dat.auc, st.dev = FALSE)
  mod.kappa <- Kappa(cmx(dat, threshold = threshold), st.dev = FALSE)
  mod.tss <- sensitivity(cmx(dat, threshold = threshold),st.dev = FALSE) + 
    specificity (cmx(dat, threshold = threshold),st.dev = FALSE) - 1
  rm(ori, p)
  names(mod) <- paste0("fold", 1:length(kdata))
  return(list("model"=allmod, "auc" = mod.auc, "kappa"=mod.kappa, "tss"= mod.tss ,"fold.models" = mod, "ObsPred" = dat)) # behar dira ere modeloak
}

#end




#' @title Cut value of the max TSS
#' @param Obs Observed values
#' @param Fit fitted values
#' @param TSS Default is TRUE
#'   
#' @return threshold that of the max TSS


kappaRepet <-  function(Obs, Fit, TSS=TRUE){
  if(sum(Obs)==0) stop("\n The observed data only contains 0")
  tab <- as.data.frame(matrix(0, nrow=101, ncol=2))  ### 
  
  if(length(unique(Fit))==1){
    Misc<-table(as.vector(Fit) >= as.numeric(unique(Fit)), Obs) 
    if(TSS!=TRUE) a <- KappaStat(Misc)
    else a <- TSS.Stat(Misc)
    TP <- Misc[4]
    TN <- Misc[1]
    ca0 <- (TN * 100)/sum(Misc[,1])
    ca1 <- (TP * 100)/sum(Misc[,2])
    if(is.na(ca0)) ca0<-0
    if(is.na(ca1)) ca1<-0
    if(TSS!=TRUE) invisible(list(Kappa=0, CutOff=unique(Fit), TP=TP, se=ca1, TN=TN, sp=ca0))
    else invisible(list(TSS=0, CutOff=unique(Fit), TP=TP, se=ca1, TN=TN, sp=ca0))
  }
  else{
    Quant <- quantile(Fit)
    for(j in 0:100){
      Seuil <- Quant[1] + (j*((Quant[5] - Quant[1])/100))
      Misc<-table(Fit >= Seuil, Obs)
      if(TSS!=TRUE) a <- KappaStat(Misc) else a <- TSS.Stat(Misc)
      if(!is.na(a)) if(a > 0) {tab[j+1, 1] <- Seuil; tab[j+1, 2] <- a}
      rm(Misc, Seuil)
    }
    
    t <- max(tab[,2],na.rm=TRUE)
    seuil <- tab[tab[,2]==t,1]   
    if(t > 0) {
      Misc<-table(Fit >= seuil[1], Obs)
      TP <- Misc[4]
      TN <- Misc[1]
      ca0 <- (TN * 100)/sum(Misc[,1])
      ca1 <- (TP * 100)/sum(Misc[,2])
      if(TSS!=TRUE) invisible(list(Kappa = t, CutOff = seuil[1], TP = TP, se = ca1, TN = TN, sp = ca0))
      else invisible(list(TSS = t, CutOff = seuil[1], TP = TP, se = ca1, TN = TN, sp = ca0))
    }
    else {
      if(TSS!=TRUE) invisible(list(Kappa = 0, CutOff = 0, TP = 0, se = 0, TN = 0, sp = 0))
      else invisible(list(TSS = 0, CutOff = 0, TP = 0, se = 0, TN = 0, sp = 0))
    }
  }
}

#end



#' @title Internarl function for KappaRepet 
#' @param Misc data


TSS.Stat <-
  function(Misc)
  {
    if(dim(Misc)[1]==1){
      if(row.names(Misc)[1]=="FALSE") Misc<-rbind(Misc, c(0,0))
      else {
        a<-Misc
        Misc<-c(0,0)
        Misc<-rbind(Misc, a)
        
      }
    }
    n <- sum(Misc)
    a <- Misc[1,1]
    b <- Misc[1,2]
    c <- Misc[2,1]
    d <- Misc[2,2]
    sens<-a/(a+c)
    spec<-d/(b+d)
    K <- (sens + spec) - 1        #TSS
    return(K)
  }

#end
