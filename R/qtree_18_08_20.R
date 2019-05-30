

get.qtree<-function(f1, n.num, title=NA, digits, s_size=TRUE){
  # takes a frame f1, split number n.num, root name, and rounding digits
  i <-which(f1$node==n.num)
  if(n.num==1) #first split
    if(f1$end[i]=="E"){
      val=round(f1$mean[i], digits=digits)
      return(cat("\\", "Tree [.{", title, "} ", "\\", "\\ ", 
                 ifelse(s_size, paste0("\\", "fbox{n = ", f1$n[i], "} ", "\\", "\\ "),""),
                 "\\", "fbox{", paste(round(val, digits)), "} ]", sep="")) 
    }
    else #first split not end
      return(cat("\\", "Tree [.{", title, "} ", 
                 get.qtree(f1, n.num=2, digits=digits, s_size=s_size), 
                 get.qtree(f1, n.num=3, digits=digits, s_size=s_size), "]", sep="")) 
  
  
  if((n.num %% 2)==0) # left hand split (LHS)
    if(f1$end[i]=="E"){ # end
      val=round(f1$mean[i], digits=digits)
      if(f1$cat[i]==0) # numeric split
        return(paste("[.{$", f1$var[i], " \\", "leq ", 
                     round(unlist(f1$xval[i]),digits), " $} ",
                     "{", "\\", "fbox{node ", n.num, "} ", "\\", "\\ ",
                     ifelse(s_size, paste0("\\", "fbox{n = ", f1$n[i], "} ", "\\", "\\ "),""),
                     "\\", "fbox{value ", paste(round(val, digits)), "}} ]", sep="")) 
      else  #categorical split
        return(paste("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                     paste(unlist(f1$xval[i]),collapse=","), "\\", "}$} ",
                     "{", "\\", "fbox{node ", n.num, "}", " \\", "\\ ",
                     ifelse(s_size, paste0("\\", "fbox{n = ", f1$n[i], "} ", "\\", "\\ "),""),
                     "\\", "fbox{value ", paste(round(val, digits)), "}} ]", sep="")) 
    }  
    else #LHS not end node
      if(f1$cat[i]==0)
        return(paste("[.{$", f1$var[i], " \\", "leq ", 
                     round(unlist(f1$xval[i]),digits), " $} ",
                     get.qtree(f1, n.num=2*n.num, digits=digits, s_size=s_size), 
                     get.qtree(f1, n.num=2*n.num+1, digits=digits, s_size=s_size), "]", sep="")) 
      else  #categorical
        return(paste("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                     paste(unlist(f1$xval[i]),collapse=","), "\\", "}$} ",
                     get.qtree(f1, n.num=2*n.num, digits=digits, s_size=s_size), 
                     get.qtree(f1, n.num=2*n.num+1, digits=digits, s_size=s_size), "]", sep="")) 
  
  else #right hand split (RHS)
    if(f1$end[i]=="E"){ #end
      val=round(f1$mean[i], digits=digits)
      if(f1$cat[i]==0) #numeric split
        return(paste("[.{$", f1$var[i], " > ", 
                     round(unlist(f1$xval[i]),digits), " $} ",
                     "{", "\\", "fbox{node ", n.num, "} ", "\\", "\\ ",
                     ifelse(s_size, paste0("\\", "fbox{n = ", f1$n[i], "} ", "\\", "\\ "),""),
                     "\\", "fbox{value ", paste(round(val, digits)), "}} ]", sep="")) 
      else  #categorical split
        return(paste("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                     paste(unlist(f1$xval[i]),collapse=","), "\\", "}$} ",
                     "{", "\\", "fbox{node ", n.num, "}", " \\", "\\ ",
                     ifelse(s_size, paste0("\\", "fbox{n = ", f1$n[i], "} ", "\\", "\\ "),""),
                     "\\", "fbox{value ", paste(round(val, digits)), "}} ]", sep=""))
    }
    else #RHS not end node
      if(f1$cat[i]==0)
        return(paste("[.{$", f1$var[i], " > ", round(unlist(f1$xval[i]),digits), " $} ",
                     get.qtree(f1, n.num=2*n.num, digits=digits, s_size=s_size), 
                     get.qtree(f1, n.num=2*n.num+1, digits=digits, s_size=s_size), "]", sep=""))  
      else #RHS is categorical
        return(paste("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                     paste(unlist(f1$xval[i]), collapse=","), "\\","}$} ",
                     get.qtree(f1, n.num=2*n.num, digits=digits, s_size), 
                     get.qtree(f1, n.num=2*n.num+1, digits=digits, s_size), "]", sep="")) 
  
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
#' usepackage\{lscape\}
#' usepackage\{tikz-qtree\}
#'
#' @param t1 rpms object created by rpms function
#' @param title string for the top node of the tree
#' @param label string used for labeling the tree figure
#' @param caption string used for caption
#' @param digits integer number of displayed digits
#' @param s_size boolean indicating whether or not to include sample size
#' @param scale numeric factor for scaling size of tree
#' @param lscape boolean to display tree in landscape mode
#' @examples
#' {
#' # model mean of retirement account value for households with reported 
#' # retirment account values > 0 using a binary tree while accounting for 
#' # clusterd data and sample weights.
#' 
#' s1<- which(CE$IRAX > 0)
#' r1 <-rpms(IRAX~EDUCA+AGE+BLS_URBN, data = CE[s1,],  weights=~FINLWT21, clusters=~CID) 
#' 
#' # get Latex code
#' qtree(r1)
#' 
#' }
#' @export
#' @aliases rpms::qtree

qtree<-function(t1, title="rpms", label=NA, caption="", digits=2, s_size=TRUE, scale=1, lscape=FALSE){
  if(lscape) cat("\\", "begin{landscape} \n", sep="")

  cat("\\", "begin{figure}[htb] \n", sep="")
  cat("\\", "centering \n", sep="")
    
   cat("\\", "begin{tikzpicture}[scale=", scale, ", ] \n", sep="")
   cat("\\", "tikzset{every tree node/.style={align=center,anchor=north}} \n", sep="")
     get.qtree(t1$frame, n.num=1, title=title, digits=digits, s_size=s_size)
     
   cat("\\", "end{tikzpicture} \n", sep="")
    
  cat("\n")
  cat("\\", "caption{","\\", "small ", caption, "} \n", sep="")
  if(!is.na(label))
    cat("\\", "label{", label, "} \n", sep="")

  cat("\\", "end{figure} \n", sep="")
  if(lscape) cat("\\", "end{landscape} \n", sep="")
}

###################################### end qtree  ##############################################################################