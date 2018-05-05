// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include<Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//##########################################################################################################
//                    External Function Declarations: 
//##########################################################################################################

// functions defined in model_functions 
List survLm_fit(arma::colvec y, arma::mat X, arma::colvec weights); 
List survLm_model(arma::colvec y, arma::mat X, arma::colvec weights, arma::uvec strata, arma::uvec clusters);

//functions defined in variable_tests
int select_var(arma::colvec y, arma::mat mX, arma::mat vars, arma::uvec cat_vec, 
               arma::colvec weights, arma::uvec clusters, arma::uvec des_ind,
               Rcpp::StringVector vnames, arma::uword perm_reps, float pval);

//########################################################################################################



//############################################################################################################
//
//            Definitions of loss functions one for numeric and one for categorical
//
//###################### get_loss (numeric) ################################################################
// [[Rcpp::export]] 
arma::vec get_loss(arma::vec x_val,  arma::vec uq_xs, arma::colvec y, arma::mat mX, 
                   arma::colvec weights, double M) {
  
  double n = y.n_elem;
  //  arma::vec uq_xs = unique(x_val);  // sorted unique x values done in split now
  uword k = uq_xs.n_elem;       // number of unique values to check
  arma::vec l_vec(k);               
  l_vec.fill(datum::inf);     // loss vector same as unique values datum::inf
  
  for(uword i=0; i<k; ++i){
    uvec s=find(x_val<= uq_xs[i]);
    uvec not_s=find(x_val> uq_xs[i]);
    double sn=s.n_elem;
    
    if(sn>=M && (n-sn)>=M){
      
      arma::colvec res_s = survLm_fit(y.elem(s), mX.rows(s), weights.elem(s))["residuals"];
      arma::colvec res_ns = survLm_fit(y.elem(not_s), mX.rows(not_s), weights.elem(not_s))["residuals"];
      l_vec[i]=(sn/n)*dot(res_s,res_s) + ((n-sn)/n)*dot(res_ns,res_ns);
    } //end if
    
    
  } //end for i loop
  
  return l_vec;
}
//******** end get_loss numeric *****

//********************************* get_loss categorical  *****************************************************

// [[Rcpp::export]] 
arma::vec get_loss_cat(arma::imat sets, arma::vec cats, arma::vec x_val, arma::colvec y, arma::mat mX, 
                       arma::colvec weights, double M) {
  
  int n = y.n_elem;
  int nsets= sets.n_rows;
  arma::vec l_vec(nsets);               
  l_vec.fill(datum::inf);     // loss vector same as unique values datum::inf
  
  for(double j=0; j<nsets; ++j){
    
    uvec sx(n);
    
    for(int i=0; i<n; ++i){ //loop to get s
      if(any(cats.elem(find(sets.row(j)==1))==x_val(i)))
        sx(i)=1; 
      else sx(i)=0;    
    }      
    
    uvec s= find(sx);
    uvec not_s= find(ones(n)-sx);  //not_s is complement of s
    
    double sn=s.n_elem;
    
    if(sn>=M && (n-sn)>=M){
      arma::colvec res_s = survLm_fit(y.elem(s), mX.rows(s), weights.elem(s))["residuals"];
      arma::colvec res_ns = survLm_fit(y.elem(not_s), mX.rows(not_s), weights.elem(not_s))["residuals"];
      l_vec[j]=(sn/n)*dot(res_s,res_s) + ((n-sn)/n)*dot(res_ns,res_ns);
    } //end if
    
    
  } //end for j loop thru power set
  
  return l_vec;
}

//################################### end loss functions  #######################################################################



//####################################### get_set_grid  ########################################################################
// makes a grid of {0, 1} for each possible combination splitting the category values into two groups
imat get_set_grid(int C){
  
  int p=(int)pow(2, C-1)-1;
  
  imat set_grid(p, C);
  
  irowvec sets(C); //reps a row of matrix
  sets.zeros();   //change as scroll through rows
  
  for(int j=0; j<p; ++j){
    
    int i=0; //start at first element
    
    // loop through sets until find first 0
    while(i<C){
      
      if(sets(i)==0){
        sets(i)=1;
        i = C;
      }
      else{
        sets(i)=0;
        ++i;
      }
    } //end while-loop
    
    set_grid.row(j)=sets;
  } // end for-loop
  
  return set_grid; 
}
//####################################### End get_set_grid  ########################################################################


//################################################ null_split ##########################################
List null_split(){
  
  return(
    List::create(
      Named("node") =  NA_INTEGER, 
      Named("cat") =   NA_INTEGER, 
      Named("var") =   NA_STRING,
      Named("xval") =  NA_REAL,
      Named("n") =     NA_INTEGER,
      Named("loss")=   NA_REAL,
      Named("value") = NA_REAL,
      Named("cvar")=   NA_REAL,
      Named("mean")=   NA_REAL
    )
  );
  
} // end function

//############################### end null_splits ######################################################################

//############################### rbind_splits ######################################################################

// [[Rcpp::export]]
List rbind_splits(List split1, List split2) {
  
  //------------- handling if one of the splits is the null-split ---------------------------
  if(any(is_na(as<NumericVector>(split1["node"])))){
    if(any(is_na(as<NumericVector>(split2["node"])))){
      List split0=null_split();
      return(split0);
    } // split1 is empty but split2 is not
    else{return(split2);}
  } // if split1 not empty
  else{if(any(is_na(as<NumericVector>(split2["node"])))){return(split1);}}
  //------------------- done handling null ------------------------------------------------- 
  
  arma::vec node = join_vert(as<arma::vec>(split1["node"]), as<arma::vec>(split2["node"]));
  arma::vec cat = join_vert(as<arma::vec>(split1["cat"]), as<arma::vec>(split2["cat"]));
  
  //------------- get var -------------------------
  CharacterVector var1 =  split1["var"];
  CharacterVector var2 =  split2["var"];
  var1 = as< vector< string > >(var1);
  var2 = as< vector< string > >(var2);
  for (int i=0;i<var2.length();i++) {
    var1.push_back(var2(i));
  }
  
  //-------------- get xval ------------------
  vector< vec > xval1=split1["xval"];
  vector< vec > xval2=split2["xval"];
  for (uword i=0;i<xval2.size();i++) {
    xval1.push_back(xval2.at(i));
  }
  
  arma::vec n = join_vert(as<arma::vec>(split1["n"]), as<arma::vec>(split2["n"]));
  
  arma::vec loss = join_vert(as<arma::vec>(split1["loss"]), as<arma::vec>(split2["loss"]));
  
  //---------- get value -----------------------
  vector< vec > val1=split1["value"];
  vector< vec > val2=split2["value"];
  for (uword i=0;i<val2.size();i++) {
    val1.push_back(val2.at(i));
  }
  
  //-------------- get cvar -----------------------------
  vector< vec > cvar1=split1["cvar"];
  vector< vec > cvar2=split2["cvar"];
  for (uword i=0;i<cvar2.size();i++) {
    cvar1.push_back(cvar2.at(i));
  }
  
  arma::vec mean = join_vert(as<arma::vec>(split1["mean"]), as<arma::vec>(split2["mean"]));
  
  return List::create(
    Named("node") =  node, 
    Named("cat") =   cat, 
    Named("var") =   var1,
    Named("xval") =  xval1,
    Named("n") =     n,
    Named("loss")=   loss,
    Named("value") = val1,
    Named("cvar")=   cvar1,
    Named("mean")=   mean
  );
} // end function

//############################### end rbind_splits ######################################################################




//####################################### get_node  ########################################################################
// [[Rcpp::export]]
List get_node(arma::uword node, int cat, std::string vname, arma::colvec y, arma::vec mxval, arma::uvec s, List modfit){
  
  arma::uword n_y=y.n_elem;  //y is all observations in parent node
  
  if(node==1) { s=find(y<datum::inf); }
  
  arma::colvec resid=modfit["residuals"];
  double loss = ((double)(s.n_elem)/(double)(n_y))*sum(resid.t()*resid); //weighted squared residuals
  
  vector <rowvec> xval(1);
  xval.at(0)=mxval.t();
  //cout<<"Got xval"<<endl;
  
  arma::rowvec value1 = modfit["coefficients"];
  vector <rowvec> value;
  value.push_back(value1);
  //cout<<"Got value1"<<endl;
  
  arma::mat spL_covM=modfit["covM"];
  arma::rowvec  cvar1 = spL_covM.diag().t();
  vector <rowvec> cvar;
  cvar.push_back(cvar1);
  //cout<<"Got cvar"<<endl;
  
  double spL_mean = mean(y(s)-resid);
  
  return List::create(
    Named("node") =  node, 
    Named("cat") =   cat, 
    Named("var") =   vname,
    Named("xval") =  xval,
    Named("n") =     s.n_elem,
    Named("loss")=   loss,
    Named("value") = value,
    Named("cvar")=   cvar,
    Named("mean")=   spL_mean
  ); //end create
  
}//end function

//####################################### End get_node  ########################################################################


//################################################# split_rpms #######################################################################
//# main split function () uses mean squared error)
//######################################################################################################################################
// [[Rcpp::export]]
List split_rpms(arma::uword node, arma::colvec y, arma::mat mX, arma::mat X, Rcpp::StringVector vnames, arma::uvec cat_vec, arma::colvec weights, 
                arma::uvec strata, arma::uvec clusters, arma::uvec des_ind,  arma::uword bin_size, arma::uword perm_reps, float pval){
  
  arma::uword n = y.n_elem;
  
  if(n<=2*bin_size || n < mX.n_cols) return(null_split());
  
  //cout<<"going to select a variable"<<endl;
  int x = select_var(y, mX, X, cat_vec, weights, clusters, des_ind, vnames, perm_reps, pval);
  if(x ==-1) return(null_split());
  
  arma::vec uvar= unique(X.col(x)); // sorted unique x values
  
  arma::uword cat=cat_vec.at(x);
  std::string var=as<string>(vnames.at(x));
  
  //#------------------------ numeric variable  -------------------------------
  if(cat==0) {
    // cout<<"variable is numeric"<<endl;
    
    arma::vec l_vec=get_loss(X.col(x), uvar, y, mX, weights, bin_size);
    
    double loss= min(l_vec);
    
    if(loss == datum::inf){ return(null_split()); }  //end if no split
    else{
      arma::uword  min_idx = min(find(l_vec==loss));
      rowvec mxval(1); //value of x rep as vector of size one
      mxval.fill(uvar.at(min_idx)); 
      
      // can't do this here, might have to pass an integer indicator vector from R
      //if(is.integer(X.col(x))==TRUE) int mxval = (int)(x_val[min_i]);
      //else mxval <- round((x.val[min.i] + x.val[min.i + 1])/2, 2)
      
      //------- the Left node ---------------------------
      arma::uvec sL = find(X.col(x)<=mxval.at(0));
      
      List modfit = survLm_model(y.elem(sL), mX.rows(sL), weights.elem(sL), strata.elem(sL), clusters.elem(sL));
      
      List SpL=get_node(2*node, cat=0, var, y, mxval, sL, modfit);
      
      List SpL2=split_rpms(SpL["node"], y(sL), mX.rows(sL), X.rows(sL), vnames, cat_vec,
                           weights(sL), strata(sL), clusters(sL), des_ind, bin_size, perm_reps, pval);
      //-------------------------------------------------------------------------------------------------------
      
      //------- the Right node -----------------------------------------------------------------------------
      uvec sR = find(X.col(x)>mxval.at(0));
      
      modfit = survLm_model(y.elem(sR), mX.rows(sR), weights.elem(sR), strata.elem(sR), clusters.elem(sR));
      
      List SpR=get_node(2*node+1, cat=0, var, y, mxval, sR, modfit);
      
      List SpR2=split_rpms(SpR["node"], y(sR), mX.rows(sR), X.rows(sR), vnames, cat_vec,
                           weights(sR), strata(sR), clusters(sR), des_ind, bin_size, perm_reps, pval);
      
      //---------------------------------------------------------------------------------------------
      
      return rbind_splits(rbind_splits(SpL, SpL2), rbind_splits(SpR, SpR2)); 
      
    } //end else
    
  } // end numeric
  
  //#------------------------ categorical variable  -------------------------------  
  if(cat==1){
    
    //uvar is "x.cats" unique set of categories
    uword p = uvar.n_elem;  
    
    if(p<2) return(null_split());  //end if no split
    
    // get vector of uvecs for each category
    vector<uvec> cat_ind(p); //vector to hold indexs for each category
    
    for(uword c=0; c<p; ++c){
      uvec ci = find(X.col(x)==uvar(c)); //find index of each category label 
      cat_ind.at(c)= ci;
    }
    // got the vector
    
    imat power_set;
    //check if doing means
    if(mX.n_cols==1){
      
      power_set=imat(p-1, p);
      power_set.zeros();
      
      vec cat_means(p);
      cat_means.fill(datum::inf);
      
      for(uword i=0; i<p; ++i){
        //  uvec cs = find(X.col(x)==uvar(i));
        uvec cs = cat_ind.at(i);
        cat_means(i) = mean(y(cs));
      } 
      
      for(uword i=0; i<(p-1); ++i){
        uword cmin = cat_means.index_min();
        power_set.rows(i, p-2).col(cmin).fill(1); //all rows below 
        cat_means(cmin) = datum::inf;  //set min to inf then go find next smallest mean
      }
      
    } //end if doing means
    else{
     // if(p>12) warning("Categorical variable ", x, " has > than 12 categories.  This may take a long time when fitting models on each node.");
      power_set = get_set_grid(p);
    } 
    
    // NEED TO pass cind matrix to get_loss make it faster
    arma::vec l_vec=get_loss_cat(power_set, uvar, X.col(x), y, mX, weights, bin_size);
    
    uword  min_idx = l_vec.index_min();
    double loss= l_vec.at(min_idx);
    
    if(loss == datum::inf) return(null_split());  //end if no split
    else{
      
      //------- the Left node ---------------------------
      
      uvec left_cats = find(power_set.row(min_idx)==1);
      uword ncats_L = left_cats.n_elem;
      vec mxval_L = uvar(left_cats);
      
      //get index vector of all data with category left side
      uvec sL = cat_ind.at(left_cats(0));
      for(uword i=1; i<ncats_L; ++i) sL= join_vert(sL,cat_ind.at(left_cats(i)));
      
      List modfit = survLm_model(y.elem(sL), mX.rows(sL), weights.elem(sL), strata.elem(sL), clusters.elem(sL));
      
      List SpL=get_node(2*node, cat=1, var, y, mxval_L, sL, modfit);
      
      List SpL2=split_rpms(SpL["node"], y(sL), mX.rows(sL), X.rows(sL), vnames, cat_vec,
                           weights(sL), strata(sL), clusters(sL), des_ind, bin_size, perm_reps, pval);
      //-------------------------------------------------------------------------------------------------------
      
      //------- the Right node -----------------------------------------------------------------------------
      
      uvec right_cats = find(power_set.row(min_idx)==0);
      uword ncats_R = right_cats.n_elem;
      
      vec mxval_R = uvar(right_cats);
      
      //get index vector of all data with category right side
      uvec sR = cat_ind.at(right_cats(0));
      for(uword i=1; i<ncats_R; ++i) sR= join_vert(sR, cat_ind.at(right_cats(i)));
      
      modfit = survLm_model(y.elem(sR), mX.rows(sR), weights.elem(sR), strata.elem(sR), clusters.elem(sR));
      
      List SpR=get_node(2*node+1, cat=1, var, y, mxval_R, sR, modfit);
      
      List SpR2=split_rpms(SpR["node"], y(sR), mX.rows(sR), X.rows(sR), vnames, cat_vec,
                           weights(sR), strata(sR), clusters(sR), des_ind, bin_size, perm_reps, pval);
      
      //---------------------------------------------------------------------------------------------
      
      return rbind_splits(rbind_splits(SpL, SpL2), rbind_splits(SpR, SpR2)); 
      
    } //end else finite loss
    
    
  } //end categorical
  
  return(null_split());
  
} //end rpms_split

//################################# End rpms_split ##############################################################


//                             End of File
//###############################################################################################################
