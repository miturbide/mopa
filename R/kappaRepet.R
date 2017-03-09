


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
