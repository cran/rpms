// Functions defined:
//
//    list get_clus_effect
//    mat clus_perm
//    var_test
//  
// version 2 handles single psus seperately... 
// because we can handle them faster.

//can't include if using armadillo #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//#include <boost/random.hpp>

using namespace Rcpp;
using namespace arma;


//#################### get_clus_effect #####################################
//This is a support function for sample_perm
//The function clus_effect returns a matrix row 1 is a vector of the unique 
//cluster labels and row 2 is the vector of estimated cluster effect for each
//cluster.

//DO NOT EXPORT FUNCTION AFTER TESTING ( # // [[Rcpp::export]]

List get_clus_effect(arma::vec rij, arma::uvec clus){
  
  // double mu= dot(resid, weights)/sum(weights); //est of pop mean
  
  // resid -= mu;  //center resid so has mean 0  should be zero
  
  arma::uvec clus_id= unique(clus); //ids of clusters
  uword C=clus_id.n_elem; //C is number of clusters
  
  arma::ivec clus_size(C); //vector to hold cluster sizes
  arma::vec clus_eff(C); // vector to hold cluster effects
  
  //find cluster effect for each cluster 
  for(uword c=0; c<C; ++c){
    //find elements in cluster c 
    arma::uvec ci = find(clus==clus_id(c));
    clus_eff(c)= mean(rij(ci));
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
arma::mat clus_perm(arma::vec res, arma::vec weights,
                    arma::uvec clus, arma::uword M) {
  
  int  n = res.n_elem;
  
  res%=weights; //considering weighted residuals rij
  
  arma::mat r_vals(n, M+1); //matrix to hold M reordered weighted residuals
 
  r_vals.col(0)=res; // first column is original order of weighted residuals
  
  
  List clus_eff=get_clus_effect(res, clus); // get clus effects
  
  arma::uvec clusters=clus_eff["clusters"];
  arma::vec effs=clus_eff["effects"];
  arma::vec n_c=clus_eff["sizes"];
  
  // begin permutaions
  
  //perm loop
  for(uword m=1; m<=M; ++m){ // m starts at 1 because first element is original
    
    arma::vec nr=res;         // vector for the new weighted residual values 
    
    arma::vec reffs=shuffle(effs);   // shuffled cluster effects
    
    
    uword C = clusters.n_elem;
    
    for(uword i=0; i<C; ++i){
      
      uvec j=find(clus==clusters(i)); //index of each observation in cluster i
     
      
      nr(j) -= effs(i); //subtract average old effect
      nr(j) += reffs(i); //add average reshuffled effect
      
      uvec sj = shuffle(j); //shuffle labels j within cluster
      
      nr(j)=nr(sj); // new value is shuffled values
      
      
    } //end mult-cluster loop
    
    r_vals.col(m)=nr;  //mth column gets permuted value of weighted residual
    
  } //end M perms loop
  
  //  cout<<effs<<endl;
  //  cout<<singles<<endl;
  //  cout<<mults<<endl;
  return(r_vals);
  
}

//********************************* unclustered perm ***************************
// if design does not have clusters 
//use faster unclustered permuation

// [[Rcpp::export]] 
arma::mat perm(arma::vec res, arma::vec weights, arma::uword M) {
  
  int  n = res.n_elem;
  
  res%=weights; //considering weighted residuals rij
  
  arma::mat r_vals(n, M+1); //matrix to hold M reordered weighted residuals
  
  r_vals.col(0)=res; // first column is original order of weighted residuals
  
  
  // begin permutaions
  
  //perm loop
  for(uword m=1; m<=M; ++m){ // m starts at 1 because first element is original
    
    r_vals.col(m)=shuffle(res);  //mth column gets permuted value of weighted residual
    
  } //end M perms loop
  
  return(r_vals);
  
}
