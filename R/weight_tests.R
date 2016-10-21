
###################### domain_est ############################################
#
#-----------------------  Domain Estimate with Variance -------------------
#
#'domain_est
#'
#' @param t1 rpms object
#' @param domdat data frame rows inside desired domain estimate
#' @param weights vector containing the sample weights for each observation
#' @param theta string of either "mean" or "total"
#' 
#' @return list containing F-statistic and corresponding p-value
#' 
#' @description  Produces estimates of domain mean or total along with the
#'               variance estimate based on estimated regression tree
#' 
#' @aliases domain_est
#'
# #' @export
domain_est<-function(t1, domdat, weights=rep(1,nrow(domdat)), theta="mean"){
  
  if(theta=="mean"){
    if(length(t1$ln_split)<2) 
      X<-as.matrix(sum(weights*covariates(splits=t1$ln_split, domdat)))/sum(weights)
    X<-as.matrix(colSums(weights*covariates(splits=t1$ln_split, domdat)))/sum(weights)
  }
  else if(theta=="total"){
    if(length(t1$ln_split)<2) 
      X<-as.matrix(sum(weights*covariates(splits=t1$ln_split, domdat)))
    X<-as.matrix(colSums(weights*covariates(splits=t1$ln_split, domdat)))
  }
  else return("Error: theta must be 'mean' or 'total'")
  
  y <- domdat[,all.vars(t1$e_equ)[1]]
  mse <- stats::weighted.mean((y-stats::predict(t1, domdat))^2, w=weights)
  
  return(list(theta=t(X)%*%t1$ln_coef, 
              var=t(X)%*%t1$ln_coef_cov%*%X, 
              tot_var=t(X)%*%t1$ln_coef_cov%*%X + mse/nrow(domdat)))
}
#--------------- End var_hat ----------------------------


########################### t-test ###########################################
#
#----------------- t-test for Weights on Domain Estimate -----------------
#
#'ttest
#'
#' @param t1 rpms object
#' @param data data frame
#' @param dom vector of indices of the data frame that define the domain
#' @param weights vector of sample weights for each observation 
#' 
#' @return list containing the overall effect of the weights on the domain 
#'         estimate, the t-statistic, and the corresponding p-value
#' 
#' @description  Conducts a t-test for the signifcance of the weights on 
#'               on the domain estimate
#' 
#' @aliases ttest
#'
# #' @export
#' 
ttest<-function(t1, data, dom, weights){
  # require("MASS")
  y<-data[,all.vars(as.formula(t1$e_equ))[1]]
  X<-covariates(t1$ln_split, data)
  if(length(t1$ln_split)<2) 
    Ad<-as.matrix(sum((weights*X)[dom,])/sum(weights[dom]))
  else Ad<-as.matrix(colSums((weights*X)[dom,])/sum(weights[dom]))
  n=nrow(X)
  k=ncol(X)
  Z<-cbind(X, weights*X)
  
  mod<-survLm(y, Z)
  V22<-mod$covM[(k+1):(2*k), (k+1):(2*k)]
  
  delt2=as.matrix(mod$coef[(k+1):(2*k)])
  
  #iV22<-solve(V22)
  iV22<-chol2inv(V22)
  
  mu<-t(Ad)%*%delt2
  se_mu<-sqrt(t(Ad)%*%V22%*% Ad)
  tstat=mu/se_mu
  result<-list(Effect=mu, se=se_mu, tstat=mu/se_mu, Prob=1-pt(abs(tstat),(n-2*k-1)))
  
  #result<-list(mean=Ad%*%delta2, se=sqrt(Ad%*%iV22%*%Ad), tstat=tstat, Prob=(1-pf(Fstat, df1=k, df2=(n-2*k-1)))))
  #names(result)<-c("Fstat", "Prob > ")
  #  print(cbind(Coef=mod$coefficients[1:k],se=mod$stderr[1:k]))
  
  #  cat("\n\n Weighted Beta\n")
  
  return(result)
}

########################### end t-test #######################################



############################################## Wtest #################################################

#' Wtest (Wald Test for sample weights)
#' 
#' @param t1 rpms object
#' @param data data.frame containing the data
#' @param weights vector of sample weights for each observation 
#' @param strata vector of strata labels
#' @param clusters vector of cluster labels
#' @return list containing F-statistic and corresponding p-value
#' 
#' @description  Conducts a Wald's test for the significance of the weights on the splits
#' 
#' @import stats
#' 
#' @aliases Wald test
#' 
Wtest<-function(t1, data, weights=rep(1, nrow(data)),
                strata=rep(1L, nrow(data)), clusters=(1L:nrow(data))){
  # require("MASS")
  y<-data[,all.vars(stats::as.formula(t1$e_equ))[1]]
  X<-covariates(t1$ln_split, data)
  n=nrow(X)
  k=ncol(X)
  Z<-cbind(X, weights*X)
  
  mod<-survLm_model(y, Z, rep(1, n), strata, clusters)
  
  delt2=mod$coef[(k+1):(2*k)]
  
  #iV22<-solve(mod$covM[(k+1):(2*k), (k+1):(2*k)])
  iV22<-chol2inv(mod$covM[(k+1):(2*k), (k+1):(2*k)])
  
  Fstat=(delt2%*%iV22%*%delt2)/k
  
  # cat("\n------------Wald Test -------------------\n\n")
  
  #   cat("  Beta \n")
  #   print(cbind(Coef=mod$coefficients[1:k],se=mod$stderr[1:k]))
  #   
  #   cat("\n\n Weighted Beta\n")
  #   
  #   print(cbind(Coef=mod$coefficients[(k+1):(2*k)],se=mod$stderr[(k+1):(2*k)]))
  #   
  #  cat("\n\n")
  
  #  cat("\n Test Null Hypothesis: Weighted Beta = 0 \n")
  #cat("\n Test Statistic and Probability \n\n")
  result<-as.data.frame(cbind(Fstat=Fstat, 
                              Prob=(1-pf(Fstat, df1=k, df2=(n-1-2*k)))))
  names(result)<-c("Fstat", "Prob > ")
  
  return(result)
}

########################### end Wtest #####################################################






