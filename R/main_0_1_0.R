#Write a Recursive Partition Algorithm with Modeling
#Version 0.6.0
#Daniell Toth 09/04/2015
#Handles one-stage stratified and clusters samples with unequal weights
#Does not have plot function
#Uses C++ functions through Rcpp
#Uses variable select; M=1000 permutations 

###############################################################################
#
#                                            Programs for RPMS
#
###############################################################################



#################  Front End to survLm  ###################
#' survLm
#' 
#' @param e_equ formula representing the equation to fit
#' @param data data.frame 
#' @param weights formula or vector of sample weights for each observation 
#' @param strata formula or vector of strata labels
#' @param clusters formula or vector of cluster labels
#' 
#' @return list containing coefficients, covariance matrix and the residuals
#' 
#' @description  wrapper function for the C++ function survLm_model 
#' 
#' @aliases survLm
#'
survLm<-function(e_equ, data, weights=rep(1, nrow(data)), strata=rep(1L, nrow(data)), clusters=(1L:nrow(data))){
  
  if(length(all.vars(e_equ))> nrow(data)) stop("Number of parameters p > n")
  
  y<-data[,all.vars(e_equ)[1]]
  mX<-model.matrix(e_equ, data)

  survLm_model(y=as.matrix(y), X=mX, weights=weights, strata=strata, clusters=clusters)
  
}
####################################################################################################################


################################ Variable Selection Function ##############################################
#' var_select
#' 
#' @param e_equ formula object of the estimating equation
#' @param data dataframe containing partitioning and estimating variables
#'        as well as any used sample design variables
#' @param e_fn error function
#' @param weights the sample weights for each observation 
#' @param strata lable of the strata containing the observation
#' @param clusters lable of the clusters containing the observation
#' @param X, string vector containing the names of pontential variables 
#' @param perm_reps integer specifying the number of permuations
#' @param pval numeric p-value used to reject null hypothesis in permutation 
#'        test  
#' 
#' @return vector of variable names with lowest estimated p-value
#' 
#' @description  selects variables likely to result in the best split
#' 
#' 
#' @aliases var_select

var_select<-function(e_equ, data, e_fn=survLm_fit, weights=rep(1, nrow(data)), 
                     strata=rep(1, nrow(data)), clusters=(1:nrow(data)), 
                     X, perm_reps, pval){
  
  if(perm_reps < 1) return(NULL)
  
  
  y<-data[,all.vars(e_equ)[1]]
  mX<-model.matrix(e_equ, data)
  
  res <- survLm_fit(y, mX, weights)$residuals
  #res <- eval(call(e_fn, y, mX, weights))$residuals
  
  pres <- clus_perm(y=res, weights=weights, clus=clusters, M=perm_reps)
  
  if(length(X)==1) cat_vec <- is.factor(data[,X])
  else
    cat_vec <- sapply(data[, X], FUN= function(x) is.factor(x))
  
  #make categories integers for Cpp
  if(length(which(cat_vec))==1) 
    data[,X[which(cat_vec)]] <- as.integer(data[,X[which(cat_vec)]])
  else
   data[,X[which(cat_vec)]] <- sapply(data[,X[which(cat_vec)]], as.integer)

  pvec <- get_pvec(p_scores=pres, mX=mX,  vars=as.matrix(data[, X]), 
                   cat_vec=cat_vec)
  #print(X)
  #print(pvec)  
 
  ix <- which.min(pvec)

  if(pvec[ix] <= pval)   return(X[ix]) else return(NULL)
  
  #if(min(pvec)<Inf) return(X[which(pvec==min(pvec))]) else return(NULL)

}
##################################  End var.select #############################

#---example ----------
# var_select(e_equ=Women04~Women03, data=belg, X=c("Men04", "Men03"))
# var_select(e_equ=Women04~1, data=belg, e_fn="lm", 
#           X=c("Women03", "Men04", "Province", "Var1", "Var2"))
# 
# var_select(e_equ=Women04~1, data=belg, e_fn="survLm_fit", 
#            X=c("Women03", "Men04", "Province", "Var1", "Var2"))


################################################################################################################################################


################################################## The Split Algorithm ##########################################
#' split
#' 
#' @param node is an integer label of the node
#' @param data this is the data
#' @param e_equ is the estimating equation used to model the data in each node
#' @param e_fn is the function that should be used to fit the model in each
#'        node.  The function should  
#' @param l_fn function taking a numeric vector and returning the estimated loss
#' @param X vector of string names for variables to consider for splitting
#' @param weights the sample weights for each observation 
#' @param strata lable of the strata containing the observation
#' @param clusters lable of the clusters containing the observation
#' @param bin_size integer specifying the minimum number of observations each
#'        node must contain 
#' @param perm_reps integer specifying the number of permuations
#' @param pval numeric p-value used to reject null hypothesis in permutation 
#'        test 
#' 
#' @return list containing split information
#' 
#' @description  internal split function for RPMS
#' 
split<-function(node, data, e_equ, e_fn=survLm, l_fn, X, 
                weights, strata, clusters, bin_size, perm_reps, pval){

  n<-nrow(data)
  
  if(n<=2*bin_size || length(all.vars(e_equ))> nrow(data)) return(NULL)
  else {
    
    x <- var_select(e_equ, data, e_fn, weights, 
                    strata, clusters, X, perm_reps=perm_reps, pval=pval) 
    
    if(is.null(x)) return(NULL)
    
    y<-data[,all.vars(e_equ)[1]]
    mX<-model.matrix(e_equ, data) 
    
    #------------------------ numeric variable  -------------------------------
      if(is.numeric(data[,x])) {
        
        x.val<-sort(unique(data[, x])) 
        
        if(e_fn=="survLm"){
          l.vec<-get_loss(data[,x], y, as.matrix(mX), weights, M=bin_size)
        } #end if e_fn=survLm
        
        else {return(NULL)
        }#end else not survLm
        
        #------ Find the miminimum split and add left and right nodes to frame 
        if(min(l.vec)<Inf) {
          
          min.i <- which.min(l.vec)

          if(is.integer(data[,x])) mxval <- as.integer(x.val[min.i])
             else mxval <- round((x.val[min.i] + x.val[min.i + 1])/2, 2)
          
          # ------- the Left node ---------------------------
          
          s <- which(data[,x]<=mxval)
          
          
          
          sp1 <- data.frame(node=2*node, cat=0, var=x, xval=mxval, n=length(s))

          modfit <- survLm_model(y[s], as.matrix(mX[s,]), weights=weights[s], 
                                strata=strata[s], clusters=clusters[s])
         
          sp1$value <- list(modfit$coefficients)
          sp1$loss <- (length(s)/n)*l_fn(modfit$residuals)
          sp1$cvar <- list(diag(modfit$covM))
          sp1$mean <- mean(y[s]-modfit$residuals)
          sp.L <- sp1
          
          # ------- the Right node ---------
          sp1 <- data.frame(node=2*node+1, cat=0, var=x,  xval=mxval, n=(n-length(s)))

          modfit <- survLm_model(y[-s], as.matrix(mX[-s,]), weights=weights[-s], 
                                strata=strata[-s], clusters=clusters[-s])

          sp1$value <- list(modfit$coefficients)
          sp1$loss <- ((n-length(s))/n)*l_fn(modfit$residuals)
          sp1$cvar <- list(diag(modfit$covM))
          sp1$mean <- mean(y[-s]-modfit$residuals)
          sp.R <- sp1
          
        } #end if < Inf
        else return(NULL)
      } # -----end  numeric--       
      
    #-------------------Categorical variables ---------------------------
      else { 
        x.cats<-unique(data[,x])
        p<-length(x.cats)  
        
        if(e_fn=="survLm"){
          power.set<-as.matrix(expand.grid(as.data.frame(matrix(rep(c(0, 1), p), 2, p))))
          l.vec<-t(get_loss_cat(power.set, x.cats, data[,x], y, mX, weights, M=bin_size))
        } #end if survLm
        
        else{return(NULL)}# end else slow loop
        
        #------ Find the miminimum split and add left and right nodes to frame 
        if(min(l.vec)<Inf){
          
          min.i <- which.min(l.vec)

          #-------- Left node (categorical) -----------------------
          s <- which(data[,x] %in% x.cats[which(power.set[min.i,]==1)])
          sp1<-data.frame(node=2*node, cat=1, var=x)
          sp1$xval<-list(x.cats[which(power.set[min.i,]==1)])
          sp1$n<-length(s)
        
          modfit <- survLm_model(y[s], as.matrix(mX[s,]), weights=weights[s], 
                                strata=strata[s], clusters=clusters[s])
      
          sp1$value <- list(modfit$coefficients)
          sp1$loss <- (length(s)/n)*l_fn(modfit$residuals)
          sp1$cvar <- list(diag(modfit$covM))
          sp1$mean <- mean(y[s]-modfit$residuals)
          sp.L<-sp1
          
          #-------- the Right node (categorical) -------------------------
          sp1<-data.frame(node=2*node+1, cat=1, var=x)
          sp1$xval<-list(x.cats[which(power.set[min.i,]==0)])
          sp1$n<-(n-length(s))

          modfit <- survLm_model(y[-s], as.matrix(mX[-s,]), weights=weights[-s], 
                                strata=strata[-s], clusters=clusters[-s])

          sp1$value <- list(modfit$coefficients)
          sp1$loss <- ((n-length(s))/n)*l_fn(modfit$residuals)
          sp1$cvar <- list(diag(modfit$covM))
          sp1$mean <- mean(y[-s]-modfit$residuals)
          sp.R<- sp1
        } #end if < Inf
        
        else return(NULL)
        
      } #end else categorical
      
   
      return(rbind(sp.L, split(sp.L$node, data[s,], e_equ, e_fn, l_fn, X, 
                               weights[s], strata[s], clusters[s], 
                               bin_size, perm_reps, pval), 
                   sp.R, split(sp.R$node, data[-s,], e_equ, e_fn, l_fn, X, 
                               weights[-s], strata[-s], clusters[-s], 
                               bin_size, perm_reps, pval))) 
      
    
  }# end else n>2*M
} ##################################################  End Split Algorithm #######################################



########################################## Linearize ###########################################
#takes frame (f1) from tree and outputs linear logical splits
#recursive function starts with node number 1 
linearize<-function(f1, n_num=1, sp=NULL){
  
  if(!(n_num %in% f1$node)) return(NULL) #end
  
  i <-which(f1$node==n_num)
  and<-ifelse(is.null(sp),"", "&") #No & for first split
  
  if(n_num==1) #first split
    return(rbind(linearize(f1, 2), linearize(f1, 3))) 
  
  if((n_num %% 2)==0) #left side 
    if(f1$cat[i]==0){
      #  print(paste("sp is", sp, "in cat"))
      sp<-paste(sp, and, paste(f1$var[i], "<=", f1$xval[i], sep=" "))
      return(rbind(sp, linearize(f1, 2*n_num, sp), linearize(f1, 2*n_num+1, sp)))
    }  
  else {
    # print(paste("sp is", sp, "in cat"))
    sp<-paste(sp, and, paste(f1$var[i], " %in% c('", paste(f1$xval[i][[1]], collapse="','"), "')", sep=""))
    return(rbind(sp, linearize(f1, 2*n_num, sp), linearize(f1, 2*n_num+1, sp)))
  }
  else #right side
    if((2*n_num) %in% f1$node)
      if(f1$cat[i]==0){
        sp<-paste(sp, and,paste(f1$var[i], ">", f1$xval[i], sep=" "))
        return(rbind(linearize(f1, 2*n_num, sp), linearize(f1, 2*n_num+1, sp)))
      }        
  else {
    sp<-paste(sp, and, paste(f1$var[i], " %in%  c('", paste(f1$xval[i][[1]], collapse="','"), "')", sep=""))
    return(rbind(linearize(f1, 2*n_num, sp), linearize(f1, 2*n_num+1, sp)))
  }
  else return(NULL)
}
##################################### End Linearize #########################################



####################### Mark the end nodes as ends in Frame ######################
#returns the frame with "E" at the end of each frame
mark_ends<-function(frame){
  frame$end <- "" #sets new colum to blank string vector
  
  get_ends<-function(frame, n_num=1){
    if((2*n_num) %in% frame$node) 
      return(c(get_ends(frame, 2*n_num), get_ends(frame, 2*n_num+1))) 
    else return(n_num)
  }
  
  frame$end[which(frame$node %in% get_ends(frame))]="E"

  return(frame)
}


######################### end mark.ends ###########################################

#################################  covariates ##################################
#---gets regressors for simple function

covariates <- function(splits, data){
  x<-rep(0,length=nrow(data)) #first column = x satifies no splits

  if(length(splits)<=0) return(x) else
    for(i in 1:length(splits))
      x<-cbind(x, #if x satisfies split 1 in that column else 0 in that colum
               eval(parse(text=paste("ifelse(", splits[i], ", 1, 0)")), data),
               deparse.level=0)
  
  #colnames(x)<-paste("E", seq(0:(ncol(x)-1)), sep="")
  return(as.matrix(x[,-1])) 
}
##################################### End covariates ##########################


################################## get_endnodes ###############################
#
# function returns a list containing: the node_index; coef; and coef_se
# the node index is a matrix with two rows
#     1.) the binary representaion of each end node using the linear splits
#     2.) and the corresponding node number
# coef is a list containing the e_equ model coefficients for each node
# coef_se a list containing the standard error coefficients for each node
#

get_endnodes<- function(ln_sp, e_fn, e_equ, data, frame){
  two_vec <- 2^(0:(nrow(ln_sp)-1))
  end_vals <- covariates(splits=ln_sp, data)%*%two_vec
  uends <- sort(unique(end_vals))
  
  enodes <- frame$node[which(frame$end=="E")]

   #function returns order of linear splits
   get_tdir <- 
    function(cn, frame=frame){
     if((2*cn) %in% frame$node)  
         return(c(2*cn, get_tdir(2*cn, frame), get_tdir(2*cn+1, frame)))
     else return(NULL)  
     
    }

  tdir <- get_tdir(1, frame)
  
  #function returns linear splits used to get to a given node
  get_turns<-function(node, tdir=tdir){
    if(node==1) return(NULL) 
    
    if(node %% 2 == 0) { #even
      if(node %in% tdir) return(c(get_turns(node/2, tdir), node)) 
      else return(c(get_turns(node/2,tdir), NULL))
    } #end even
    else #odd
      return(c(get_turns((node-1)/2, tdir), NULL))
  }#end functio get_turns

  
  get_endval<-function(node, tdir=tdir){
    bvec <- rep(0, length(tdir))
    index<-match(get_turns(node, tdir), tdir)
    bvec[index]<-1
    bvec<-c(1,bvec)
    return(as.numeric(bvec %*% two_vec))
  } 
  
  coef<-list()
  node_index<-rbind(uends, 0)
  row.names(node_index) <- c("binary index", "node")

  for(i in 1:length(enodes)){
    index <- which(uends==get_endval(enodes[i], tdir))
    node_index[2,index] <- enodes[i]
    coef[[index]] <- unlist(frame[which(frame$node==enodes[i]), "value"])  
  } #end loop
 
   
  
 return(list(node_index=node_index, coef=coef))

}

################################### End get.endnodes ###########################


################################################################################
#                                                                              #
#                               MAIN Function rpm                              #
#                                                                              #
################################################################################
#
#' rpms  
#' 
#' @description main function producing a regression tree using variables
#'  from rp_equ to partition the data and fit the model e_equ on each node.
#'  Currently only uses data with complete cases.
#'
#' @param rp_equ formula containing all variables for partitioning
#' @param data data.frame that includes variables used in rp_equ, e_equ, 
#'        and design information
#' @param weights formula or vector of sample weights for each observation 
#' @param strata formula or vector of strata labels
#' @param clusters formula or vector of cluster labels
#' @param e_equ formula for modeling data in each node
#' @param e_fn string name of function to use for modeling 
#'        (only "survLm" is operational)
#' @param l_fn loss function (does nothing yet)
#' @param bin_size numeric minimum number of observations in each node
#' @param perm_reps integer specifying the number of permuations
#' @param pval numeric p-value used to reject null hypothesis in permutation 
#'        test 
#' 
#' @details The dependent variable must be the same for rp_equ and e_equ, 
#'  but recommended that the independent variables be different.  Categorical 
#'  variables with many categories in the rp_equ will cause the algorithm to 
#'  take a long time to run.
#' 
#' @return object of class "rpms"
#' 
#' @examples
#' {
#' # model mean of retirement contributions with a binary tree while accounting 
#' # for clusterd data
#' 
#' rpms(FINDRETX~EDUC_REF+AGE_REF+BLS_URBN+REGION, data = CE,  clusters=~CID)
#'      
#'      
#' # model linear fit between retirement contributions and amount of income
#' # with a regression tree while accounting for clusterd data
#' 
#' rpms(FINDRETX~EDUC_REF+AGE_REF+BLS_URBN+REGION, data=CE,
#'      e_equ=FINDRETX~FINCBTAX, clusters=~CID)     
#' }
#' 
#' @export
#' @aliases rpms

rpms<-function(rp_equ, data, weights=~1, strata=~1, clusters=~1, 
               e_equ=~1, e_fn="survLm", l_fn=NULL, 
               bin_size=NULL, perm_reps=500L, pval=.05){
  
  varlist <- unique(c(all.vars(rp_equ), all.vars(e_equ), all.vars(weights), 
                        all.vars(strata), all.vars(clusters)))
  
  
  if(all(varlist %in% names(data)))
     data <- na.omit(data[,varlist]) #this can be removed once we handle missing
     else {
       e1ind <- which(!(varlist %in% names(data)))[1]
       stop(paste(varlist[e1ind], " not in dataset"))
     } 
  
    #capture design information if equations for use in graphing
  des<-list(weights=~1, strata=~1, clusters=~1)
  
  if(is.null(l_fn)) l_fn <- function(x){sum(x^2)}
    else if(!(class(l_fn)=="function")) stop("l_fn is not of type function")
  
  if(is.null(bin_size)) bin_size <- ceiling(nrow(data)^(5/8))
  
  if(length(all.vars(rp_equ))==0) 
    stop("no variable for recursive partition given in formula") else
      if(!all(all.vars(rp_equ) %in% names(data)))
        stop("e_equ contains variables that are not in the data set")
  
  if(length(all.vars(e_equ))==0)
    e_equ<-formula(paste(all.vars(rp_equ)[1], 1, sep="~")) else
      if(all.vars(e_equ)[1]!=all.vars(rp_equ)[1]) 
        stop("Dependent variable must be same for rp_equ as e_equ") else
          if(!all(all.vars(rp_equ) %in% names(data)))
            stop("e_equ contains variables that are not in the data set")
  

  if(is.numeric(weights)) stopifnot(length(weights)==nrow(data)) else
    if(length(all.vars(weights))==0) {weights <- rep(1, nrow(data))} else
      if(all.vars(weights)[1] %in% names(data)){
        des$weights=weights
        weigths <- as.numeric(data[,all.vars(weights)[1]])} else 
          stop(paste("Problem with weights:",
                     all.vars(weights)[1], "not in data set"))
  
  if(is.numeric(strata) | is.factor(strata)){
    strata<-as.integer(strata)
    if(length(strata)!=nrow(data)) 
      stop("Number of strata labels != rows of data")
  } else
    if(length(all.vars(strata))==0) {strata <- rep(1L, nrow(data))} else
      if(all.vars(strata)[1] %in% names(data)) {
        des$strata=strata
        strata <- as.integer(data[,all.vars(strata)[1]])} else 
          stop(paste("Problem with strata:",
                     all.vars(strata)[1], "not in data set"))
  
  if(is.numeric(clusters) | is.factor(clusters)){
    clusters<-as.integer(clusters)
    if(length(clusters)!=nrow(data)) 
      stop("Number of cluster labels != rows of data")
  } else
    if(length(all.vars(clusters))==0) {clusters <- seq(1L:nrow(data))} else
      if(all.vars(clusters)[1] %in% names(data)) {
        des$clusters<-clusters
        clusters <- as.integer(data[,all.vars(clusters)[1]])} else
          stop(paste("Problem with clusters:",
                     all.vars(clusters)[1], "not in data set")) 
  
  
  y <- data[,all.vars(e_equ)[1]]
  nas <- which(is.na(y))
  if(length(nas) > 0){
    ndat <- data[-nas, varlist]
    warning(paste0("Response variable had ", length(nas), " missing values, 
                   which have been removed."))
    
  }
  else  ndat <- data[ ,varlist]
  

  if(length(which(sapply(ndat[,varlist], 
                         function(x) !(is.factor(x)|is.numeric(x)))))>0)
    stop("RPMS works only on numeric or factor data types.")
 
  vX=all.vars(rp_equ)[-1]
#  p<-length(all.vars(rp_equ)[-1])
  
  mX<-model.matrix(e_equ, ndat)
  if(det(t(mX)%*%mX)==0) stop("Model matrix is singular")
  
  
  #--- get root node -------#  
  fn<-survLm_model(y=as.matrix(y), X=mX, weights=weights, strata=strata, 
                   clusters=clusters)
  
  frame<-data.frame(node=1, cat=NA, var="Root", xval=NA, n=nrow(data)) 
  frame$value <- list(fn$coefficients)
  frame$loss <- l_fn(fn$residuals)
  frame$cvar <- list(diag(fn$covM))
  frame$mean=mean(y)
 
  #-----  starts recursive spliting ----
  frame <- rbind(frame, 
                 split(1, ndat, e_equ, e_fn, l_fn, vX, 
                       weights, strata, clusters, bin_size, perm_reps, pval))
  
  frame <- mark_ends(frame) # puts an 'E' after each end node on the frame
  
  lt<-rbind(" 1", linearize(frame))
  row.names(lt)<-NULL
  
  endnodes<-get_endnodes(lt, e_fn, e_equ, ndat, frame)
  
  
  t1<-list(rp_equ=rp_equ, e_equ=e_equ, frame=frame, ln_split=lt, 
           end_nodes=endnodes, coef=endnodes$coef,
           survey_design=des)
  
  #fit<-survLm_fit(predict.rpm(t1, data), covariates(lt, data), s.weights)
  fit<-survLm_model(y, covariates(lt, ndat), weights, strata, clusters)
  
  t1$ln_coef<-as.numeric(fit$coefficients)
  t1$ln_coef_cov<-fit$covM
  
  class(t1)<-c("rpms")
  
  return(t1)
  
}
#############################################  End rpms ####################################################
