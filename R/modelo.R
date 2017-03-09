
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
#'  \item{allmod }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{mod }{fitted model with partitioned data}
#'  \item{p }{cross model prediction}
#' 
#' @details This function calculates the AUC with the function "auc" from package 
#' "PresenceAbsence". Package SDMTools must be detached.
#' 
#' 
#' 
#' @author M. Iturbide 
#' 
#'
#' @export
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



modelo<-function(kdata, data, algorithm = c("glm","svm","maxent","mars","randomForest","cart.rpart","cart.tree"), weighting = FALSE, threshold = 0.5){
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
    print(sum(all.weights))
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
