#include <RcppArmadillo.h>
#include<Rcpp.h>

using namespace Rcpp;
using namespace arma;

List survLm_fit(arma::colvec y, arma::mat X, arma::colvec weights); //declaration: defined in model_functions

// ********************** For Splitting algorithm *************************** 

//Code to do numeric split loop

// [[Rcpp::export]] 
arma::vec get_loss(arma::vec x_val, arma::colvec y, arma::mat mX, 
                   arma::colvec weights, double M) {
   
   double n = y.n_elem;
   arma::vec uq_xs = unique(x_val);  // unique x values
   int k = uq_xs.n_elem;       // number of unique values to check
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
//***************************** end numeric split ****************************************************************


//***********************************************************************************************************************


//********************************* categorical split ******************************************************
//Code to do categorical split loop

// [[Rcpp::export]] 
arma::vec get_loss_cat(arma::mat sets, arma::vec cats, arma::vec x_val, arma::colvec y, arma::mat mX, 
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



