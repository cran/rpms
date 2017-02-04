// Functions defined:
//    double peak_num
//    double peak_cat
//    list get_clus_effect
//    mat clus_perm
//    var_test
//  
// version 2 handles single psus seperately... 
// because we can handle them faster.

//can't include if using armadillo #include <Rcpp.h>
#include <RcppArmadillo.h>
//#include <boost/random.hpp>

using namespace Rcpp;
using namespace arma;
 

 //#################### clus_effect #####################################
 //This is a support function for sample_perm
 //The function clus_effect returns a matrix row 1 is a vector of the unique 
 //cluster labels and row 2 is the vector of estimated cluster effect for each
 //cluster.
 
 //DO NOT EXPORT FUNCTION AFTER TESTING ( # // [[Rcpp::export]] )

 List get_clus_effect(arma::vec resid, arma::vec weights, arma::uvec clus){
   
   double mu= dot(resid, weights)/sum(weights); //est of pop mean
   
   resid -= mu;  //center resid so has mean 0
   
   arma::uvec clus_id= unique(clus); //ids of clusters
   uword C=clus_id.n_elem; //C is number of clusters
   
   arma::ivec clus_size(C); //vector to hold cluster sizes
   arma::vec clus_eff(C); // vector to hold cluster effects
   
   //find cluster effect for each cluster 
   for(uword c=0; c<C; ++c){
     //find elements in cluster c 
     arma::uvec ci = find(clus==clus_id(c));
     clus_eff(c)= dot(resid(ci), weights(ci))/sum(weights(ci));
     clus_size(c)=ci.n_elem;
   }
   
   return List::create(
     Named("clusters") = clus_id,
     Named("effects") = clus_eff,
     Named("sizes") = clus_size
   );
 }

//#################### resid_clus_perm ######################################
//Main Function returns a n x (M+1) matrix with  
//containing the original weighted residuals and the next M values are the
// M permutations of the weighted residuals 
//the group lables and permuted cluster effects for clustered data.
//The second component is the p value computed from the M+1 vector

// [[Rcpp::export]] 
arma::mat clus_perm(arma::vec y, arma::vec weights,
                       arma::uvec clus, arma::uword M) {

  int  n = y.n_elem;
  
  arma::mat r_vals(n, M+1); //matrix to hold M reordered values
  r_vals.col(0)=y; // first column is original order

  List clus_eff=get_clus_effect(y, weights, clus); // get clus effects

  
  arma::vec effs=clus_eff["effects"];
  arma::uvec clusters=clus_eff["clusters"];
  arma::vec n_c=clus_eff["sizes"];
  
  arma::uvec singles = find(n_c == 1);
  arma::uvec mults = find(n_c > 1);
  

  for(uword m=1; m<=M; ++m){ // m starts at 1 because first element is original
    
    arma::vec ny=y;         // vector for the new y values 
    arma::vec reffs=shuffle(effs);   // shuffled cluster effects
    

    //handle single clusters
    
    ny(singles) -= effs(singles); //subtract old effects
    ny(singles) += reffs(singles); //add reshuffled effect
    
    
    //loop over each multi-cluster

    uword C = mults.n_elem;
    
    for(uword c=0; c<C; ++c){
      
      uvec j=find(clus==clusters(mults(c))); //index of each observation in cluster c
      ny(j)=shuffle(y(j)); //shuffle residual values within cluster
      
      ny(j) -= effs(mults(c)); //subtract average old effect
      ny(j) += reffs(mults(c)); //add average reshuffled effect

    } //end mult-cluster loop
    
    r_vals.col(m)=ny;  //mth column gets permuted value

  } //end M perms loop
  
//  cout<<effs<<endl;
//  cout<<singles<<endl;
//  cout<<mults<<endl;
  return(r_vals);

}

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

double peak_cat(arma::vec score, arma::vec var){
  
  
  arma::vec cat_id= unique(var);
  uword C=cat_id.n_elem;
  
  arma::vec cat_sums(C); // category sums
  
  for(uword i=0; i<C; ++i){
    arma::uvec s = find(var == cat_id(i));
    cat_sums(i)=sum(score(s));
  }
  
  uvec pos = find(cat_sums > 0);
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
  
  for(uword j=0; j <p; ++j){
    
    arma::mat chi = p_scores;
    chi.each_col() %= mX.col(j); 
    
    arma::vec pks(m);
    double a;

    if(cat==1) {
      a=peak_cat(chi.col(0), var);
      
      for(uword i=0; i <m; ++i)
        pks(i)=peak_cat(chi.col(i+1), var);
    }
    
    else {
      a=peak_num(chi.col(0), var);
      
      for(uword i=0; i <m; ++i)
        pks(i)=peak_num(chi.col(i+1), var);
    }
    
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
  
  //loop over each variable
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


