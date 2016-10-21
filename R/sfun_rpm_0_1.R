#File Containing All support functions for RPM that are intended to be made public
#
# Version 0.1 
#
# Contains: node_plot, predict, print, plot, in_node,
# 
# 
#

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
#   params <- list(...)
#   
#   if("variable" %in% names(params)) variable <-params$variable
#   else variable <- NULL
#   
#   if(is.null(variable)){
#     if(length(all.vars(x$e_equ))>1)
#       variable <- all.vars(x$e_equ)[2]
#     else variable <- all.vars(x$rp_equ)[2]
# 
#   }
#   else  if(!(variable %in% names(data) && is.character(variable)))
#     stop(paste("variable", paste(variable, collapse=""),
#                "is not a name in the data set"))
#   
#   if(is.numeric(variable)) coline <-TRUE else coline <-FALSE
# 
#   yvariable<-all.vars(x$e_equ)[1]
# 
#   #for each row get its associated end node
#   ends<-rowSums(t(t(covariates(x$ln_split, data))*(2^(0:(nrow(x$ln_split)-1)))))
#   uends<-sort(unique(ends))
# 
#   #  par(mfrow=c(length(uends), 1), mar=c(2, 2, 2, 2))
# 
#   for(i in 1:length(uends)){
# 
#     if(i%%4==1){
#      # if(i>1) dev.next()
#       par(mfrow=c(min(4, length(uends)-i+1), 1), mar=c(3, 3, 3, 3))
#     }
#     nind <- which(x$end_nodes$node_index["binary index",]==uends[i])
#     plot(data[ends==uends[i],c(variable, yvariable)])
#     title(main=paste0("Node ", x$end_nodes$node_index[2, nind]),
#           ylab=yvariable, xlab=variable)
# 
#     if(coline)
#       abline(coef=x$end_nodes$coef[[nind]][c(1,2)],
#              col=2)
# 
#   }
# 
# 
# } #end plot.rpms
# 
# ################################### End plot.rpms #################################################################


###################################### in_node ################################################################
#' in_node
#' 
#' @param node integer label of the desired end-node.
#' @param t1 \code{rpms} object
#' @param data dataframe containing the variables used for the recursive 
#'       partitioning. 
#' 
#' @description  Get index of elements in dataframe that are in the specified 
#'               end-node of an \code{rpms} object
#' 
#' @return vector of indexes for observations in the end-node.
#' 
#' @aliases rpms::in_node
#' 
#' @examples
#' # model linear fit between retirement contributions and amount of income 
#' r1 <-rpms(FINDRETX~EDUC_REF+AGE_REF+BLS_URBN+REGION, data=CE,
#'      e_equ=FINDRETX~FINCBTAX, clusters=~CID) 
#' 
#' summary(CE$FSALARYX[in_node(node=2, r1, data=CE)])
#'
#' summary(CE$FSALARYX[in_node(node=6, r1, data=CE)])
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

################################## plot_node ##################################
#' node_plot
#' 
#' 
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
#' # model linear fit between retirement contributions and amount of income 
#' r1 <-rpms(FINDRETX~EDUC_REF+AGE_REF+BLS_URBN+REGION, data=CE,
#'      e_equ=FINDRETX~FINCBTAX, clusters=~CID) 
#' 
#' # plot node 2
#' node_plot(r1, node=2, data=CE)
#' 
#' #' # plot node 7
#' node_plot(r1, node=7, data=CE)
#' 
#' 
#' \dontrun{node_plot(r1, node=11, data=CE)}
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
#' @examples
#' {
#' # get rpms model of mean retirement contribution by several factors
#' r1 <- rpms(FINDRETX~EDUC_REF+AGE_REF+BLS_URBN+REGION, data = CE)
#' 
#' # first 10 predicted means
#' predict(r1, CE[1:10, ])
#' 
#' }
#'
#' @export
predict.rpms<-function(object, newdata, ...){
  
 # pars<-list(...)  #currently does not take other parameters
  
#  if("newdata" %in%)
  t1<-object
  
  y<-numeric(nrow(newdata)) #vector for new y values
  bin<-2^(0:(nrow(t1$ln_split)-1)) #vector 1, 2, 4, 8 ... 2^(p-1)
  ends<-as.matrix(rowSums(t(t(covariates(splits=t1$ln_split, newdata))*bin)))
  
  get_y<-function(end){  
    beta<-unlist(t1$coef[which(t1$end_nodes$node_index[1,]==end)])
    X<-model.matrix(t1$e_equ, newdata[which(ends==end), ]) 
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
    sptab<-cbind(t1$ln_split, t1$ln_coef, sqrt(diag(t1$ln_coef_cov)))
    colnames(sptab)<-c("Splits", "Coefficients", "SE")
    print(sptab, quote=F, zero.print=".")
    cat("\n \n")
  
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


########################### group_rpms ######################################################

group_rpms<-function(tree, newdata){
  bin<-2^(0:(nrow(tree$ln_split)-1))
  grp<-as.numeric(rowSums(t(t(covariates(splits=tree$ln_split, newdata))*bin)))
  glab<-rep(0, length(grp))
  for(i in unique(grp))
  
    glab[which(grp==i)]<-which(sort(unique(grp))==i)
  
  return(glab)
  
}

########################### End group.rpm ######################################################





########################################### function ################################################################################

##################################################################################################################################

########################################### function ################################################################################

##############################################################################################################################

########################################### function ################################################################################

##############################################################################################################################


