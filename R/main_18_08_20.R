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


#################  Front End to survLm_model  ###################
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
#' @keywords internal
survLm<-function(e_equ, data, weights=rep(1, nrow(data)), strata=rep(1L, nrow(data)), clusters=(1L:nrow(data))){
  
  if(length(all.vars(e_equ))> nrow(data)) stop("Number of parameters p > n")
  
  y<-data[,all.vars(e_equ)[1]]
  mX<-model.matrix(e_equ, data)

  survLm_model(y=as.matrix(y), X=mX, weights=weights, strata=strata, clusters=clusters)
  
}
####################################################################################################################

##############################################################################################
#
#   functions to make output of rpms more readable
#
##############################################################################################

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
#for internal use
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


#################################################################################################
#
#  make nice puts the List returned from Cpp into a data.frame and puts nice names on the labels
#
#################################################################################################
make_nice<-function(frame){
  
  # test is a List of 9
  
  nice.frame<-as.data.frame(frame[-c(4,7,8)])
  names(nice.frame)<-c("node", "cat", "var", "n", "loss", "mean")
  
  nice.frame$xval=frame$xval
  
  nice.frame$value=frame$value
  names(nice.frame$value)<-"value"
  
  nice.frame$cvar=frame$cvar
  names(nice.frame$cvar)<-"cvar"
  
  return(nice.frame[c("node", "cat", "var", "xval", "n", "loss", "value", "cvar", "mean")])
  
}#end function 
######################  end make nice #############################################################




#####################################################################################################
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
#' @param l_fn loss function (ignored)
#' @param bin_size numeric minimum number of observations in each node
#' @param perm_reps integer specifying the number of thousands of permuation 
#'        replications to use to estimate p-value
#' @param pval numeric p-value used to reject null hypothesis in permutation 
#'        test 
#' 
#' @return object of class "rpms"
#' 
#' 
#' @examples
#' {
#' # model mean of retirement account value for households with reported 
#' # retirment account values > 0 using a binary tree while accounting for 
#' # clusterd data and sample weights.
#' 
#' s1<- which(CE$IRAX > 0)
#' rpms(IRAX~EDUCA+AGE+BLS_URBN, data=CE[s1,], weights=~FINLWT21, clusters=~CID)
#'
#'                  
#' # model linear fit between retirement account value and amount of income
#' # conditioning on education and accounting for clusterd data for households 
#' # with reported retirment account values > 0
#' 
#' rpms(IRAX~EDUCA, e_equ=IRAX~FINCBTAX, data=CE[s1,], weights=~FINLWT21, clusters=~CID)    
#' 
#' }
#' @export
#' @aliases rpms
#' 
rpms<-function(rp_equ, data, weights=~1, strata=~1, clusters=~1, 
                  e_equ=~1, e_fn="survLm", l_fn=NULL,
                  bin_size=NULL, perm_reps=1000L, pval=.05){
  


  # Unnecessary l_fn is ignored ##########
  # if(is.null(l_fn)) l_fn <- function(x){sum(x^2)}
  # else if(!(class(l_fn)=="function")) stop("l_fn is not of type function")
  
  if(is.null(bin_size)) bin_size <- ceiling(nrow(data)^(11/20))
  else 
    if(bin_size<2) {
      warning("bin_size set to 2")  
      bin_size <- 2}

  
  #============= format data set ================================================
  varlist <- unique(c(all.vars(rp_equ), all.vars(e_equ), all.vars(weights), 
                      all.vars(strata), all.vars(clusters)))  
  
  #------- check variables against data set ------
  if(length(all.vars(rp_equ))==0) 
    stop("no variable for recursive partition given in formula") else
      if(!all(all.vars(rp_equ) %in% names(data)))
        stop("rp_equ contains variables that are not in the data set")
  
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
    
    
    ##########################   handle categorical variables for C++  ########################
    
    # ------------------identify categorical variable ------
    # need to handle length 1 separately
    if(length(vX)==1) cat_vec <- is.factor(data[,vX]) 
    else
      cat_vec <- sapply(data[, vX], FUN=is.factor)
    
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
      if(n_cats==1) {
        data[,vX[which(cat_vec)]] <- fn_q(data[,vX[which(cat_vec)]])
      }
      else{
        data[,vX[which(cat_vec)]] <- lapply(data[,vX[which(cat_vec)]], fn_q)
      }
      
      # ----- store original levels for each categorical variable ------------
      cat_table <- list()
      
      # ----- function to turn categories into integers for C++ 
      for(i in 1:n_cats){
        
        cat_table[[i]] <- levels(data[,vX[which(cat_vec)[i]]]) # --store original levels in cat_table
        
        #---- replace original levels with integer -------
        levels(data[,vX[which(cat_vec)[i]]]) <- (1:length(levels(data[,vX[which(cat_vec)[i]]])))

      } #end for loop
      
    } # done with if categorical variables
    ######################################################################################################
    
    
    #--------- reduce data to only y with observations and required variables 
    # n_o<-nrow(data)
    data <- na.omit(data[,varlist]) #remove any other missing
    #  nas<-abs(n_o -nrow(data))
    #  if(nas > 0) 
    #    warning(paste0("Response variable had ", nas, " missing values, which have been removed."))
 
  #-----------------------------------------------------------------------
 
#===================== Design Information =================================
  #capture design information if equations for use in graphing
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

  #-----------------------------------------------------------------------
  # get model matrix and variable data in matrix form
  #-----------------------------------------------------------------------
  mX<-model.matrix(e_equ, data)
  if(det(t(mX)%*%mX)==0) stop("Model matrix is singular")
  
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
  
  #return(frame)
  endnodes<-get_endnodes(lt, e_fn, e_equ, data, frame)
  
  
  des_string <- if(sum(des_ind)==0) "Simple Random Sample"
  else paste(c("unequal probability of selection", "stratified", "cluster")[which(des_ind)], "sample design")
  
  callvals <-list(des_string=des_string, perm_reps=perm_reps, pval=pval)
  
  
  
  t1 <- list(rp_equ=rp_equ, e_equ=e_equ, frame=frame, ln_split=lt, 
             end_nodes=endnodes, coef=endnodes$coef, callvals=callvals)
  
  
  fit<-survLm_model(y, covariates(lt, data), weights, strata, clusters)
  
  t1$ln_coef<-as.numeric(fit$coefficients)
  t1$ln_coef_cov<-fit$covM  
  
  
  class(t1)<-c("rpms")
  
  return(t1)
  
}
#############################################  End rpms ####################################################
