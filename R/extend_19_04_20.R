#######################################################################################################
#
# Code to produce random forest from rpms function 
# 
#  functions: forest, predict.forest
#  
#  Last updated: 4/13/2019
#
########################################################################################################


##############################################################################################################
#                       internal function f_rpms Fast rpms for use with rpms_forest 
################################################################################################################

f_rpms<-function(rp_equ, data, weights=~1, strata=~1, clusters=~1, 
               e_equ=~1, e_fn="survLm", l_fn=NULL,
               bin_size=NULL, perm_reps=100L, pval=.25,
               cat_vec, n_cats, cat_table, vX, des_ind){

  #-----------------------------------------------------------------------
  # get model matrix and variable data in matrix form
  #-----------------------------------------------------------------------
  mX<-model.matrix(e_equ, data)
  if(det(t(mX)%*%mX)==0) stop("Model matrix is singular")
 
  #------ recurisive partitioning variables ------
  X<- data.matrix(as.data.frame(data[,vX]))
  y <- data[,all.vars(e_equ)[1]]
  #-------------------------------------------------------------------
  
  
  ##################################  Recursive Partitioning  ####################################################
  #          calls C++ funtions get_node and split_rpms 
  ################################################################################################################
  frame <-
    rbind_splits(get_node(node=1, cat=as.integer(NA), vname="Root", y=as.matrix(y), mxval=as.matrix(NA), s=as.matrix(0), 
                          modfit=survLm_model(y=as.matrix(y), X=mX, weights=weights, strata=strata, clusters=clusters)), 
                 split_rpms(node=1, y, mX, X, vX, cat_vec, weights=weights, strata=strata, clusters=clusters, des_ind, 
                            bin_size=bin_size, perm_reps=perm_reps, pval=pval))
  ################################################################################################################
  
  
  #################  formating resulting tree frame for R  #######################################################
  frame <-  make_nice(frame)
  frame <- mark_ends(frame) # puts an 'E' after each end node on the frame
  
  # return to original the factor levels in the data 
  if(n_cats>0){
    for(i in 1:n_cats){
      #---- return to original levels -------
      levels(data[,vX[which(cat_vec)[i]]]) <- unlist(cat_table[[i]])
    }
  }
  
  # put original level names in the frame
  cat_splits<-which(frame$cat==1)
  for(i in cat_splits){
    vis <- which(vX[cat_vec]==frame$var[[i]]) #variable location in cat_table
    frame$xval[i] <- list(cat_table[[vis]][unlist(frame$xval[i])]) #replace numbers with factor names
  }
  
  
  lt<-rbind(" 1", linearize(frame))
  row.names(lt)<-NULL
  
  endnodes<-get_endnodes(lt, e_fn, e_equ, data, frame)
  
  
  t1 <- list(rp_equ=rp_equ, e_equ=e_equ, frame=frame, ln_split=lt, 
             end_nodes=endnodes, coef=endnodes$coef)
  
  
  fit<-survLm_fit(y, covariates(lt, data), weights)
  
  t1$ln_coef<-as.numeric(fit$coefficients)
  
  class(t1)<-c("rpms")
  
  return(t1)
  
}
#############################################  End f_rpms ####################################################


###############################################################################################
#                                                                                             #
#                                            forest                                           #
#                                                                                             #
###############################################################################################
#' rpms_forest
#' 
#' @description produces a random forest using rpms to create the individual
#'              trees.
#'
#' @param rp_equ formula containing all variables for partitioning
#' @param data data.frame that includes variables used in rp_equ, e_equ, 
#'        and design information
#' @param weights formula or vector of sample weights for each observation 
#' @param strata formula or vector of strata labels
#' @param clusters formula or vector of cluster labels
#' @param e_equ formula for modeling data in each node
#' @param bin_size numeric minimum number of observations in each node 
#' @param perm_reps integer specifying the number of permuations
#' @param pval numeric p-value used to reject null hypothesis in permutation 
#'        test 
#' @param f_size integer specifying the number of trees in the forest        
#' @param cores integer number of cores to use in parallel if > 1 (not implemented)
#' 
#' @return object of class "rpms"
#' 
#' @aliases forest
#' @export
rpms_forest <- function(rp_equ, data, weights=~1, strata=~1, clusters=~1, 
                  e_equ=~1, bin_size = NULL, perm_reps=100, pval=.25, 
                  f_size=200, cores=1){

  
  if(is.null(bin_size)) bin_size <- ceiling(nrow(data)^(1/2))
  else 
    if(bin_size<20) {
      warning("bin_size set to 20")  
      bin_size <- 20}
  
  #============= format data set ================================================
  varlist <- unique(c(all.vars(rp_equ), all.vars(e_equ), all.vars(weights), 
                      all.vars(strata), all.vars(clusters)))  
  
  #------- check variables against data set ------
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
  
  #---------- check all variables are in data set --------------
  if(!all(varlist %in% names(data))){
    e1ind <- which(!(varlist %in% names(data)))[1]
    stop(paste(varlist[e1ind], " not in dataset"))
  }
  
  #-----check all variables are numeric or categorical
  if(length(which(sapply(data[,varlist], 
                         function(x) !(is.factor(x)|is.numeric(x)))))>0){
    stop("RPMS works only on numeric or factor data types.")
  }
  
  #------ recurisive partitioning variables ------
  vX=all.vars(rp_equ)[-1]
  
  
  #--------- reduce data to only y with observations and required variables 
  data <- na.omit(data[,varlist])

  #-----------------------------------------------------------------------
  
  #===================== Design Information =================================

  des_ind <- c(FALSE, FALSE, FALSE)
  des<-list(weights=~1, strata=~1, clusters=~1)
  
  #------ Sample Weights ----------------------------
  if(is.numeric(weights)) { 
    if(length(weights)!=nrow(data)) stop("Number of strata labels != rows of data")
  } else
    if(length(all.vars(weights))==0) {
      weights <- rep(1, nrow(data))
    } else 
      if(all.vars(weights)[1] %in% names(data))
      {
        des$weights=weights
        weights <- as.numeric(data[,all.vars(weights)[1]])
        if(var(weights)>0) {
          des_ind[1] <- TRUE 
        } 
      } else {stop(paste("Problem with weights:",
                         all.vars(weights)[1], "not in data set"))}
  
  #------ Strata Labels ----------------------------
  if(is.numeric(strata) | is.factor(strata)){
    strata<-as.integer(strata)
    if(length(strata)!=nrow(data)) 
      stop("Number of strata labels != rows of data")
  } else
    if(length(all.vars(strata))==0) {strata <- rep(1L, nrow(data))} else
      if(all.vars(strata)[1] %in% names(data)) {
        des$strata=strata
        des_ind[2] <- TRUE
        strata <- as.integer(data[,all.vars(strata)[1]])} else 
          stop(paste("Problem with strata:",
                     all.vars(strata)[1], "not in data set"))
  
  #------ Cluster Labels ---------------------------- 
  if(is.numeric(clusters) | is.factor(clusters)){
    clusters<-as.integer(clusters)
    if(length(clusters)!=nrow(data)) 
      stop("Number of cluster labels != rows of data")
  } else
    if(length(all.vars(clusters))==0) {clusters <- seq(1L:nrow(data))} else
      if(all.vars(clusters)[1] %in% names(data)) {
        des$clusters <- clusters
        des_ind[3] <- TRUE
        clusters <- as.integer(data[,all.vars(clusters)[1]])} else
          stop(paste("Problem with clusters:",
                     all.vars(clusters)[1], "not in data set")) 
  #===================================================================================  
  
  des_string <- if(sum(des_ind)==0) "Simple Random Sample"
  else paste(c("unequal probability of selection", "stratified", "cluster")[which(des_ind)], "sample design")
  

  y <- all.vars(rp_equ)[1] 
  N<-nrow(data)
  xp <- length(vX) # number of variables
  
  
  tree <- vector(mode="list", length=f_size) # forest
  
  get_trees<-function(){
    
    # ----- randomize tree growth ------------------------------
    
    #--randomize data set ----- 
    s <- sample(N, size=N, replace = TRUE) #bootstap sample
    #-------------------------------------------------------
    
    #---randomize variables used  
    fvars <- paste(vX[sample(xp, size=sample(xp, size=1))], collapse="+")
    f_equ <- as.formula(paste0(y, "~", fvars))
    #---------------------------------------------------------
    fvX=all.vars(f_equ)[-1] 
    ####################   handle categorical variables for C++  ########################
    
    # ------------------identify categorical variable ------
    # need to handle length 1 separately
    if(length(fvX)==1) cat_vec <- is.factor(data[,fvX]) 
    else
      cat_vec <- sapply(data[, fvX], FUN=is.factor)
    
    n_cats<-length(which(cat_vec))
    
    #---------- There are categorical variables ------------
    if(n_cats>0){
      
      # ----- function to turn NA into ? category
      fn_q <- function(x){
        #x<-as.integer(x)
        #x[which(is.na(x))]<- (min(x, na.rm=TRUE)-1)
        nas <- which(is.na(x))
        if(length(nas>0)){
          x <- factor(x, levels=c(levels(x), "?"))
          x[nas]<- factor("?")
        }
        return(as.factor(x))
      } # end internal function fn
      # 
      
      # ---------- turn each NA into ? category
      if(length(which(cat_vec))==1) {
        data[,fvX[which(cat_vec)]] <- fn_q(data[,fvX[which(cat_vec)]])
      }
      else{
        data[,fvX[which(cat_vec)]] <- lapply(data[,fvX[which(cat_vec)]], fn_q)
      }
      
      # ----- store original levels for each categorical variable ------------
      cat_table <- list()
      
      # ----- function to turn categories into integers for C++ 
      for(i in 1:n_cats){
        cat_table[[i]] <- levels(data[,fvX[which(cat_vec)[i]]]) # --store original levels in cat_table
        #---- replace original levels with integer -------
        levels(data[,fvX[which(cat_vec)[i]]]) <- (1:length(levels(data[,fvX[which(cat_vec)[i]]])))
      }
      
    } # done with if categorical variables
    #------------------------------------------------------------------------------------
    
    #-------end randomize ------------------------------------
    
    return(f_rpms(rp_equ=f_equ, data[s,], weights[s], strata[s], clusters[s],
                  e_equ=e_equ, perm_reps=perm_reps, pval=pval, bin_size = bin_size,
                  cat_vec=cat_vec, n_cats = n_cats, cat_table = cat_table, vX=fvX,
                  des_ind=des_ind))
    
    }#---end get_trees ---------------------------
  
  tree<-replicate(f_size, get_trees(), simplify = FALSE)
  
 f1<-list(tree=tree)
  class(f1)<-c("rpms_forest")
  
  return(f1)
  
} #end forest

################ end forest ######################################

### ------------------------------- END FOREST -----------------------------------------------------------------------------------



###################################### predict.rpms_forest ###########################################################
#' predict.rpms_forest
#' 
#' @param object  Object inheriting from  \code{rpms_forest}
#' @param newdata data frame with variables to use for predicting new values. 
#' @param ...	further arguments passed to or from other methods.
#' 
#' @description  Gets predicted values given new data based on \code{rpms_forest} model.
#' 
#' @return vector of predicticed values for each row of newdata
#'
#' @export
predict.rpms_forest <- function(object, newdata, ...){
  
  pars<-list(...)
  
  if(!("etype" %in% names(pars))) etype <- "mean" else
    etype <- pars$etype
  
  estM <-sapply(object$tree, predict.rpms, newdata=newdata)
  if(etype=="median")
    return(apply(estM, 1, median))
  else  return(rowMeans(estM)) 
  
} 

########################### end predict.rpms_forest ############################

