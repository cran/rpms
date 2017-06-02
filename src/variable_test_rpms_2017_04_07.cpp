
// Functions defined:
//
//    double peak_num
//    double peak_cat
//
//    List var_test
//    List get_pvec
//  


//can't include if using armadillo #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//#include <boost/random.hpp>

using namespace Rcpp;
using namespace arma;
 


//********************* peak functions for num and cat *************************
//Returns the highest value of the cumulative score function ordered
//by the categories to find model variable instability

double peak_num(arma::vec score, arma::vec var){

  uvec sindx = sort_index(var);
  
  double peak=max(abs(cumsum(score(sindx))));
  
  return peak;
}


//Returns the highest value of the cumulative score function ordered
//by the numeric variable to find model variable instability

//[[Rcpp::export]]
double peak_cat(arma::vec score, arma::vec var){
  
  arma::vec cat_id= unique(var); // get the unique categories
  uword C=cat_id.n_elem;
  
  arma::vec cat_sums(C); // make vector to hold category sums
  
  for(uword i=0; i<C; ++i){
    arma::uvec s = find(var == cat_id(i));
    cat_sums(i)=sum(score(s));
  }
  
  
  uvec pos = find(cat_sums >= 0);
  uvec neg = find(cat_sums < 0);
  
  double peak, a, b;
  
  a=fabs(sum(cat_sums(neg)));
  b=sum(cat_sums(pos));
  
  if(a > b) peak=a; else peak=b;
  
  return peak;
}
//*********************** end peak functions ********************************

// [[Rcpp::export]] 
List var_test(arma::mat p_scores, arma::mat mX, arma::vec var, int cat){
  
  uword p = mX.n_cols;
  uword m = p_scores.n_cols-1;
  arma::vec pvals(p);
  double max_peak=0;
  
  // do for each variable in the model matrix
  for(uword j=0; j <p; ++j){
    
    arma::mat chi = p_scores;
    chi.each_col() %= mX.col(j); 
    
    arma::vec pks(m);
    double a;

    // ****** categorical variable *******
    if(cat==1) {
      
      a=peak_cat(chi.col(0), var);
      
      // if only one category return p-value = 1
      if(arma::vec(unique(var)).n_elem <= 1) {
        return List::create(
          Named("pval") = 1,
          Named("peak") = a
        );
      } // end if only 1 category

      for(uword i=0; i <m; ++i){
        pks(i)=peak_cat(chi.col(i+1), var);
      } // end for-loop
      
    } // end if categorical
    
    
    // ******* numeric variable ********* 
    else {  
      a=peak_num(chi.col(0), var);
     
      for(uword i=0; i <m; ++i)
        pks(i)=peak_num(chi.col(i+1), var);
      
    } //end else numeric
    
    // compute p-value
    uvec s = find(pks>=a);
    pvals(j) = s.n_elem/(float)m;
    
    if(max_peak<a) max_peak=a;
    
  }//end model matrix loop
  
  //return(min(pvals));
  return List::create(
    Named("pval") = min(pvals),
    Named("peak") = max_peak
  );
}



// [[Rcpp::export]] 
List get_pvec(arma::mat p_scores, arma::mat mX, 
                       arma::mat vars, arma::ivec cat_vec){
                         
  uword p = vars.n_cols;
  
  arma::vec pvec(p), peaks(p); 
  pvec.fill(datum::inf);
  peaks.fill(0);
  
  // loop over each variable
  for(uword i=0; i<p; ++i){
   // arma::colvec var vars.col(i);
    List test_results = var_test(p_scores, mX, vars.col(i), cat_vec[i]);
    pvec(i) = test_results["pval"];
    peaks(i) = test_results["peak"];
  }
  
  return List::create(
    Named("pvec") = pvec,
    Named("peaks") = peaks
  );
  
}
