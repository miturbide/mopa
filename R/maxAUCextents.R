

#' @title Michaelis-Mentel model fitting for extracting an index of threshold background extents 
#' @description Given a matrix as returned by function loadTextValues
#' extracts an index of the extents at which the AUC (or other test statistic) overcomes the Vm coefficient 
#' obtained from the Michaelis-Menten model fitting.
#' 
#' @param testmat Object returned by function loadTestValues (matrix with extents 
#' in columns and species in rows)
#' @param diagrams Logical. If TRUE diagrams of the fitted Michaelis-Menten model are returned
#'   
#' @return Named integer. Index of the extents at which the AUC, kappa or tss 
#' statistics obtain the maximum score
#' 
#' 
#' 
#' 
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples
#' data(presaus)
#' data(biostack)
#' ##modeling
#' modirs <-allModeling(data = presaus, varstack = biostack, k = 10, "mars") 
#' ##loading#'  
#' auc_mars <-loadTestValues(data = presaus, "auc", "mars") 
#' 
#' ind<-indextent(auc_mars)
#' 
#' @export
#' @import lattice

indextent<-function (testmat, diagrams = FALSE) {
  e <- list()
  dat<-list()
  
  for (i in 1:nrow(testmat)) {
    y<-na.omit(testmat[i,])
    x<-as.integer(sub("km", x=names(y), replacement= ""))[1:length(y)]
    group<-rep(rownames(testmat)[i], length(y))
    
    fm0 <- lm(y ~ x) # simpler model to get better starting values
    ken<-(fm0$fitted.values-coef(fm0)[1])
    micmen <- nls(y~Vm*x/(K+x), start=list(Vm=max(na.omit(y)), K=x[which(ken==min(ken))])) 
    dat[[i]]<- data.frame(x,y, group)
    
    #--------------------------
    a<-coef(micmen)[1] 
    qs <- which(na.omit(testmat[i, ]) > a)
    e[[i]] <- qs[which(qs == min(qs))]
    
  }
  
  if (diagrams == TRUE){
    
    dfa<-do.call("rbind", dat)
    
    pl<-xyplot( y ~ x|group, data = dfa, ylab="AUC", xlab="Background extent",
                panel=function(x, y){
                  ## add lines to the graph which denote means of x and y
                  panel.xyplot(x,y)
                  panel.lines(x, predict(nls(y~Vm*x/(K+x), 
                                             start=list(Vm=max(na.omit(y)), 
                                                        K=x[which(lm(y ~ x)$fitted.values-coef(lm(y ~ x))[1]==
                                                                    min(lm(y ~ x)$fitted.values-coef(lm(y ~ x))[1]))])), x),
                              col="black", lwd = 1.5)
                  panel.abline(coef(nls(y~Vm*x/(K+x), 
                                        start=list(Vm=max(na.omit(y)), 
                                                   K=x[which(lm(y ~ x)$fitted.values-coef(lm(y ~ x))[1]==
                                                               min(lm(y ~ x)$fitted.values-coef(lm(y ~ x))[1]))])))[1],
                               lty = 5, col = "deeppink3", lwd = 1.5)
                  
                  panel.abline(v = min(x[which(y > coef(nls(y~Vm*x/(K+x), 
                                                            start=list(Vm=max(na.omit(y)), 
                                                                       K=x[which(lm(y ~ x)$fitted.values-coef(lm(y ~ x))[1]==
                                                                                   min(lm(y ~ x)$fitted.values-coef(lm(y ~ x))[1]))])))[1])]),
                               lty = 5, col = "darkorange1", lwd = 1.5)
                  
                  
                  
                })
    
    print(pl)
  }
  return(unlist(e))
}