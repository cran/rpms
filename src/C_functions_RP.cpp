#include <RcppArmadillo.h>
#include<Rcpp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]] 
List survLm_fit(arma::colvec y, arma::mat X, arma::colvec weights){
   
  //weights = sqrt(weights);
  arma::mat wX = X;
  wX.each_col() %= sqrt(weights);  // multiply each column by sqrt(weights)
 
  arma::mat iA = arma::trans(wX)*wX;

    // check for colinearitiy in X before inverting
  if(arma::det(iA)==0){ 
  
    arma::colvec coef(X.n_cols, fill::zeros);
    arma::colvec resid(y.n_elem); resid.fill(datum::inf);
    
    return List::create(Rcpp::Named("coefficients") = coef,
                        Rcpp::Named("residuals")    = resid); 
    
  } // end if singular
  else{
    
    iA = inv_sympd(iA);
    
    vec wy=y;
    wy %=sqrt(weights);
    
    arma::colvec coef = iA*(arma::trans(wX)*wy);
    
    // residuals-- these are unweighted
    arma::colvec resid = y - X*coef;   
    
    
    return List::create(Rcpp::Named("coefficients") = coef,
                        Rcpp::Named("residuals")    = resid);
  }
    

}


//**************** getting design based estimates of std.error for coefficient estimates ***


arma::mat survLM_covM(arma::colvec resid, arma::mat X, 
                      arma::colvec weights, arma::ivec strata, arma::ivec clusters) {
  
  size_t n = X.n_rows, p = X.n_cols; 
  
  arma::colvec sweights = sqrt(weights);
  X.each_col() %= sweights;  // multiply each column by sqrt(weights)
  
  arma::mat iA = arma::trans(X)*X;
  if(arma::det(iA)==0){
    arma::mat covM(p, p);
    covM.fill(datum::inf);
    return covM;
  }
  else{
    iA = inv(iA);

    arma::mat D = X;
    D.each_col() %= resid;  //has weight Wt since X and resid have sqrt(Wt)
    
    arma::mat G(p, p, fill::zeros);
    
    arma::ivec h_labs=unique(strata);
    //cout<<"h_labs is "<< endl << h_labs << endl;
    uword H = h_labs.n_elem;
    
    for(uword h=0; h<H; ++h){
      arma::mat G_h(p, p, fill::zeros);
      uvec uh = find(strata==h_labs[h]);
      //cout<<"h is "<<h<<endl;      
      //cout<<"current strat label is "<<h_labs[h]<<endl;      
      //cout<<"strata index is : \n" << uh << endl;
      
      arma::ivec ids = unique(clusters.elem(find(strata==h_labs[h])));
      uword nh = ids.n_elem;
      
      for(uword i=0; i<nh; ++i ){
        uvec s=find(clusters.elem(uh) == ids[i]);  //in strata h and cluster i 

        arma::rowvec dhi = sum(D.rows(s));
        
        G_h += (trans(dhi)*dhi);
        
      }
  
      if(nh <= 1) {G+=datum::inf;}  // single psu
      else G += (nh/(nh-1))*G_h; // 
      
    }
    
    G=((n-1)/(n-p))*G;
    
    arma::mat covM = iA*G*iA;
    
    return covM;
  }

} //End get_CovM



//' Fit a linear model using data collected from a complex sample
//' 
//' @param y A vector of values
//' @param X The design matrix of the linear model  
//' @param weights A vector of sample weights for each observation 
//' @param strata A vector of strata labels
//' @param clusters A vector of cluster labels
//' 
//' @return list containing coefficients, covariance matrix and the residuals
//' 
// [[Rcpp::export]] 
List survLm_model(arma::colvec y, arma::mat X, arma::colvec weights, 
                  arma::ivec strata, arma::ivec clusters) {
  
    List Lnmod = survLm_fit(y, X, weights);

    arma:mat covM = survLM_covM(Lnmod["residuals"], X,
                                weights, strata, clusters);
    
    return List::create(Rcpp::Named("coefficients") = Lnmod["coefficients"],
                        Rcpp::Named("covM")         = covM,
                        Rcpp::Named("residuals")    = Lnmod["residuals"]);
  
}


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


List get_wald(arma::colvec y, arma::mat X) {
  
  arma::mat iXX = inv(trans(X)*X);
  arma::colvec coef = iXX*(trans(X)*y);
  arma::colvec resid = y - X*coef; 
  int n=y.n_elem;
  double A=0;
  
  for(int i=0; i<n; ++i){
    A += as_scalar(X.row(i)*trans(X.row(i)));
  }
  
  return(A);
  
}