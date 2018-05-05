#File Containing All support functions for RPM that are intended to be made public
#
# Version 0.1 
#
# Contains: end_nodes, in_node, node_plot, predict, print, (plot)
# 
# 

###################################### end_nodes ################################################################
#' end_nodes
#' 
#' @param t1 \code{rpms} object
#' 
#' @description  Get vector of end-node labels
#' 
#' @return vector of lables for each end-node.
#' 
#' @aliases rpms::end_nodes
#' 
#' @examples 
#' {
#' # model mean of retirement account value for households with reported 
#' # retirment account values > 0 using a binary tree while accounting for 
#' # clusterd data and sample weights.
#' 
#' s1<- which(CE$IRAX > 0)
#' r1 <-rpms(IRAX~EDUCA+AGE+BLS_URBN, data = CE[s1,],  weights=~FINLWT21, clusters=~CID) 
#'  
#'  end_nodes(r1)
#' }
#' 
#' @export
end_nodes<-function(t1){
 return(t1$end_nodes$node_index["node",]) 
}
# ################################### End end_nodes #################################################################


###################################### in_node ################################################################
#' in_node
#' 
#' @param node integer label of the desired end-node.
#' @param t1 \code{rpms} object
#' @param data dataframe containing the variables used for the recursive 
#'       partitioning. 
#' 
#' @description  Get index of elements in dataframe that are in the specified 
#'               end-node of an \code{rpms} object.  A "which" function for end-nodes.
#' 
#' @return vector of indexes for observations in the end-node. 
#' 
#' @aliases rpms::in_node
#' 
#' @examples
#' {
#' # model mean of retirement account value for households with reported 
#' # retirment account values > 0 using a binary tree while accounting for 
#' # clusterd data and sample weights.
#' 
#' s1<- which(CE$IRAX > 0)
#' r1 <-rpms(IRAX~EDUCA+AGE+BLS_URBN, data = CE[s1,],  weights=~FINLWT21, clusters=~CID) 
#' 
#' # Get summary statistics of CUTENURE for households in end-nodes 7 and 8 of the tree
#'
#' if(7 %in% end_nodes(r1)) 
#'   summary(CE$CUTENURE[in_node(node=7, r1, data=CE[s1,])])
#' if(8 %in% end_nodes(r1)) 
#'   summary(CE$CUTENURE[in_node(node=8, r1, data=CE[s1,])])
#' }
#' 
#' @export 
in_node<-function(node, t1, data){
  
  bin<-2^(0:(nrow(t1$ln_split)-1)) #vector 1, 2, 4, 8 ... 2^(p-1)
  ends<-as.matrix(rowSums(t(t(covariates(splits=t1$ln_split, data))*bin)))
   #gets the binary index value for each observation
  
  
  endnodes <- t1$end_nodes$node_index["node",] #all end-nodes
  
  index <- match(node, endnodes) #which node was requestd
  
  if(is.na(index)){
    stop(paste(node, " is not an end-node of the tree"))
    return(NULL)}
  else return(which(ends==t1$end_nodes$node_index["binary index",index]))
  #returns rows number that have binary index matching the requested index
}

###################################### End in_node ############################

################################## node_plot ##################################
#' node_plot
#' 
#' @param t1 \code{rpms} object
#' @param node integer label of the desired end-node. 
#' @param data data.frame that includes variables used in rp_equ, e_equ, 
#'        and design information
#' @param variable string name of variable in data to use as x-axis in plot
#' @param ...	further arguments passed to plot function.      
#' 
#' 
#' @description  plots end-node of object of class \code{rpms}
#' 
#' @aliases rpms::node_plot
#' 
#' @import graphics
#' 
#' @examples{
#' 
#' # model mean of retirement account value for households with reported 
#' # retirment account values > 0 using a binary tree while accounting for 
#' # clusterd data and sample weights.
#' 
#' s1<- which(CE$IRAX > 0)
#' r1 <-rpms(IRAX~EDUCA+AGE+BLS_URBN, data = CE[s1,],  weights=~FINLWT21, clusters=~CID) 
#' 
#' # plot node 6 if it is an end-node of the tree
#' if(6 %in% end_nodes(r1))
#'   node_plot(t1=r1, node=6, data=CE[s1,])
#' 
#' # plot node 6 if it is an end-node of the tree
#' if(8 %in% end_nodes(r1))
#'   node_plot(t1=r1, node=8, data=CE[s1,])
#' 
#' }
#'
#' @export
#'
node_plot<-function(t1, node, data, variable=NA, ...){
  
  #if no variable provided use first variable in estimating equation
  #as x axis to plot on
  if(is.na(variable)){
    if(length(all.vars(t1$e_equ))>1)
      variable <- all.vars(t1$e_equ)[2]
    else variable <- all.vars(t1$rp_equ)[2]
    
  } #if variable provided check to make sure it is in the dataset
  else  if(!(variable %in% names(data) && is.character(variable)))
    stop(paste("variable", paste(variable, collapse=""), 
               "is not a name in the data set"))

  
  # y variable is the y from the estimating equation 
  yvariable<-all.vars(t1$e_equ)[1]
  
  #which element has coefficients for that node
  nind <- which(t1$end_nodes$node_index["node",]==node)
  
  # if the variable chosen is numeric produce the best fit line in the plot
  if(is.numeric(data[nind, variable])) coline <-TRUE else coline <-FALSE
  
  #get index of observations contained in end node
  eindx<-in_node(node=node, t1=t1, data=data)
  plot(data[eindx, c(variable, yvariable)])
  title(main=paste0("Node ", node), ylab=yvariable, xlab=variable)
  
  #find variable in estimating equation being used to graph against
  vindx <- match(variable, all.vars(t1$e_equ))
  
  if(coline && !is.na(vindx) && vindx>1)
      abline(coef=t1$end_nodes$coef[[nind]][c(1,vindx)],
             col=2)
  
  
} #end node_plot
################################### End plot.rpms #################################################################


###################################### predict.rpms ################################################################
#' predict.rpms
#' 
#' @param object  Object inheriting from  \code{rpms}
#' @param newdata data frame with variables to use for predicting new values. 
#' @param ...	further arguments passed to or from other methods.
#' 
#' @description  Predicted values based on \code{rpms} object
#' 
#' @return vector of predicticed values for each row of newdata
#' 
#' @aliases predict
#'
#' @examples{
#' 
#' # get rpms model of mean Soc Security income for families headed by a 
#' # retired person by several factors
#' r1 <-rpms(SOCRRX~EDUCA+AGE+BLS_URBN+REGION, 
#'           data=CE[which(CE$INCNONWK==1),], clusters=~CID) 
#' 
#' r1
#' 
#' # first 10 predicted means
#' predict(r1, CE[10:20, ])
#' 
#' }
#'
#' @export
predict.rpms<-function(object, newdata, ...){
  
 # pars<-list(...)  #currently does not take other parameters
  
  t1<-object
  
  new_equ <- t1$e_equ[-2]
  vX=all.vars(t1$rp_equ)[-1]
  varlist <- unique(c(all.vars(t1$rp_equ[-1]), all.vars(t1$e_equ)[-1]))  

  newdata<-newdata[,varlist, drop=FALSE]
  
  # ------------------identify categorical variable ------
  # need to handle length 1 separately
  if(length(vX)==1) {
    cat_vec <- is.factor(newdata[,vX])
    } 
  else
    cat_vec <- sapply(newdata[,vX], FUN=is.factor)
  
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
    
    # ---------- now turn each NA into ? category
    if(length(which(cat_vec))==1) {
      newdata[,vX[which(cat_vec)]] <- fn_q(newdata[,vX[which(cat_vec)]])
    }
    else{
      newdata[,vX[which(cat_vec)]] <- lapply(newdata[,vX[which(cat_vec)]], fn_q)
    }
  } #end if n_cats>0
  
  #newdata <- na.omit(newdata) #remove any other missing
  if(length(which(is.na(newdata)))>0) stop("Remove any missing in numeric variables")
  
  y<-numeric(nrow(newdata)) #vector for new y values
  
  bin <- 2^(0:(nrow(t1$ln_split)-1)) #vector 1, 2, 4, 8 ... 2^(p-1)
  ends <- as.matrix(rowSums(t(t(covariates(splits=t1$ln_split, newdata))*bin)))

  get_y<-function(end){  
    
    beta <- unlist(t1$coef[which(t1$end_nodes$node_index[1,]==end)])
    X <- model.matrix(new_equ, newdata[which(ends==end), ,drop=FALSE]) 
    y[ends==end]<-X %*% beta
    
    return(y)
  }

  
  return(rowSums(apply(unique(ends), 1, get_y)))                
  
}

################################### End predict #################################################################


################################## print.rpms ###################################################################
#' print.rpms
#' 
#' @param x \code{rpms} object
#' @param ...	further arguments passed to or from other methods.
#' 
#' @description  print method for class \code{rpms}
#' 
#' @aliases print
#'
#' @export
print.rpms<-function(x, ...){

  if(class(x)!="rpms") stop("argument must be of class rpms") else t1 <- x
  
  cat("\n")
  cat("RPMS Recursive Partitioning Equation \n")
  print(t1$rp_equ, showEnv=FALSE)
  cat("\n")
  cat("Estimating Equation \n")
  print(t1$e_equ, showEnv=FALSE)
  cat("\n")

  if(length(all.vars(t1$e_equ))==1){ #no e_equ just show node means
    
    if("ln_coef_cov" %in% names(t1)){ #if it has a cov matrix
      sptab<-cbind(t1$ln_split, t1$ln_coef, sqrt(diag(t1$ln_coef_cov)))
      colnames(sptab)<-c("Splits", "Coefficients", "SE")
      print(sptab, quote=F, zero.print=".")
      cat("\n \n")
    }
    else { #does not have cov mat
      sptab<-cbind(t1$ln_split, t1$ln_coef)
      colnames(sptab)<-c("Splits", "Coefficients")
      print(sptab, quote=F, zero.print=".")
      cat("\n \n")
    }

  
    } else { #e_equ is provided
      #make matrix of estimated coefficients (colums) for each end node (row)
      ends<-which(t1$frame$end=="E")
      coef_names<-list(node=t1$frame$node[ends], 
                       coefficients=c("1", paste(all.vars(t1$e_equ)[-1])))
      coef_mat<-matrix(data=unlist(t1$frame$value[ends]),
                       nrow=length(ends), 
                       ncol=length(coef_names$coefficients), 
                       byrow=TRUE, dimnames=coef_names)
      
      sptab<-t1$ln_split
      colnames(sptab)<-c("Splits")
      print(sptab, quote=F, zero.print=".")
      cat("\n")
      print(coef_mat, quote=F, zero.print=".")
      cat("\n \n")
  
  } #end if e_equ provided
  
}#end print.rpm

################################### End print.rpm #################################################################

################################## plot.rpms ###################################################################
# #' plot.rpms
# #'
# #' @param x \code{rpms} object
# #' @param ...	further arguments passed to or from other methods.
# #'
# #' @description  plott method for class \code{rpms}
# #'
# #' @aliases plot
# #' 
# plot.rpms<-function(x, ...){
# 
# 
# } #end plot.rpms
# 
# ################################### End plot.rpms #################################################################


# ########################### group_rpms ######################################################
# 
# group_rpms<-function(tree, newdata){
#   bin<-2^(0:(nrow(tree$ln_split)-1))
#   grp<-as.numeric(rowSums(t(t(covariates(splits=tree$ln_split, newdata))*bin)))
#   glab<-rep(0, length(grp))
#   for(i in unique(grp))
#   
#     glab[which(grp==i)]<-which(sort(unique(grp))==i)
#   
#   return(glab)
#   
# }
# 
# ########################### End group.rpm ######################################################





########################################### function ################################################################################

##################################################################################################################################

########################################### function ################################################################################

##############################################################################################################################

########################################### function ################################################################################

##############################################################################################################################


