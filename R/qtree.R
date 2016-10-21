

get.qtree<-function(f1, n.num=1, title="name", digits=2){
  # takes a frame f1, split number n.num, root name, and rounding digits
  i <-which(f1$node==n.num)
  
  if(n.num==1) #first split
    if(f1$end[i]=="E"){
      val=round(f1$mean[i], digits=digits)
      return(cat("\\", "Tree [.{", title, "} ", "\\", "\\ ", "\\", 
                 "fbox{", val, "} ]", sep="")) 
    }
  else
    return(cat("\\", "Tree [.{", title, "} ", 
               get.qtree(f1, 2), get.qtree(f1, 3), "]", sep="")) 
  
  
  if((n.num %% 2)==0) # left hand split 
    if(f1$end[i]=="E"){ # end
      val=round(f1$mean[i], digits=digits)
      if(f1$cat[i]==0) # numeric split
        return(paste("[.{$", f1$var[i], " \\", "leq ", 
                     round(unlist(f1$xval[i]),digits), " $} ",
                     "{", "\\", "fbox{node ", n.num, "} ", "\\", "\\ ",
                     "value ", val, "}} ]", sep="")) 
      else  #categorical split
        return(paste("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                     paste(unlist(f1$xval[i]),collapse=","), "\\", "}$} ",
                     "{", "\\", "fbox{node ", n.num, "}", " \\", "\\ ",
                     "\\", "fbox{value ", val, "}} ]", sep="")) 
    }  
  else #not end node
    if(f1$cat[i]==0)
      return(paste("[.{$", f1$var[i], " \\", "leq ", 
                   round(unlist(f1$xval[i]),digits), " $} ",
                   get.qtree(f1, 2*n.num), 
                   get.qtree(f1, 2*n.num+1), "]", sep="")) 
  else  #in cat
    return(paste("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                 paste(unlist(f1$xval[i]),collapse=","), "\\", "}$} ",
                 get.qtree(f1, 2*n.num), 
                 get.qtree(f1, 2*n.num+1), "]", sep="")) 
  
  else #right hand split
    if(f1$end[i]=="E"){ #end
      val=round(f1$mean[i], digits=digits)
      if(f1$cat[i]==0) #numeric split
        return(paste("[.{$", f1$var[i], " > ", 
                     round(unlist(f1$xval[i]),digits), " $} ",
                     "{", "\\", "fbox{node ", n.num, "}", " \\", "\\ ",
                     "\\", "fbox{value ", val, "}} ]", sep="")) 
      else  #categorical split
        return(paste("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                     paste(unlist(f1$xval[i]),collapse=","), "\\", "}$} ",
                     "{", "\\", "fbox{node ", n.num, "}", " \\", "\\ ",
                     "\\", "fbox{value ", val, "}} ]", sep=""))
    }  
  else #not end node
    if(f1$cat[i]==0)
      return(paste("[.{$", f1$var[i], " > ", round(unlist(f1$xval[i]),digits), " $} ",
                   get.qtree(f1, 2*n.num), get.qtree(f1, 2*n.num+1), "]", sep=""))  
  else  #in cat
    return(paste("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                 paste(unlist(f1$xval[i]), collapse=","), "\\","}$} ",
                 get.qtree(f1, 2*n.num), get.qtree(f1, 2*n.num+1), "]", sep="")) 
  
  return("Error in qtree")
  
}

################################################################################
#' qtree
#' 
#' Code to write a latex qtree plot
#' takes a rpm frame and returns latex code to produce qtree
#' uses linearize as a guide
#' Produces text code to produce tree structure in tex document
#' Requires using LaTex packages and the following commands in preamble of 
#' LaTex doc: 
#' usepackage{lscape}
#` usepackage{qtree}
#'
#' @param t1 rpms object created by rpms function
#' @param title text for the top node of the tree
#' @param digits integer number of displayed digits 
#' @examples
#' {
#' # get rpms model of mean retirement contribution by several factors
#' r1 <-rpms(FINDRETX~EDUC_REF+AGE_REF+BLS_URBN+REGION, data=CE,
#'      e_equ=FINDRETX~FINCBTAX, clusters=~CID) 
#' 
#' # get Latex code
#' qtree(r1)
#' 
#' }
#' @export
#' @aliases rpms::qtree

qtree<-function(t1, title="rpms", digits=2){
  ls<-ifelse(nrow(t1$frame)>8, 1, 0)
  if(ls==1) cat("\\", "begin{landscape} \n", sep="")
  cat("\\", "begin{figure}[ht] \n", sep="")
  cat("\\", "centering \n", sep="")
  cat("\\", "footnotesize \n", sep="")
  
  get.qtree(t1$frame, n.num=1, title=title, digits=digits)
  
  cat("\n")  
  cat("\\", "end{figure} \n", sep="")
  if(ls==1) cat("\\", "end{landscape} \n", sep="")
}

###################################### end qtree  ##############################################################################