#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List rescaled( arma::mat beta_temp, arma::mat xi_temp, arma::vec Ustar_dems){
  int N = beta_temp.n_rows;
  int P = Ustar_dems.size() - 1;
  
  for( int i = 0; i < N; ++i ){
    
    for( int p = 0; p < P; ++p ){
      int Rp = Ustar_dems[ p + 1 ] - Ustar_dems[ p ]; 
      double summed = 0;
      
      for( int r = Ustar_dems[ p ]; r < Ustar_dems[ p + 1 ]; ++r){
        summed += std::abs( xi_temp( i, r ) );
      }
      
      if( summed > 0 ){
        for( int r = Ustar_dems[ p ]; r < Ustar_dems[ p + 1 ]; ++r){
          xi_temp( i, r ) = ( Rp/summed )*xi_temp( i, r );
        }
        
        beta_temp( i, 3*p ) = ( summed/Rp )*beta_temp( i, 3*p );
      }

    }
    
  }
  
  List rescaled( 2 );
  rescaled[ 0 ] = beta_temp;
  rescaled[ 1 ] = xi_temp;
  return rescaled;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# ustar_dems = c(0,2,4)
# beta_temp = matrix( data = c(rep(0.1,6), rep(0.2,6)), nrow = 2 )
# xi_temp = matrix( data = c(rep(c(3,4),4)), nrow = 2 ) 
# rescaled(beta_temp,xi_temp,ustar_dems)
*/

// 

// beta_temp
//   
//      [n,1] [,2] [,3] [n,4] [,5] [,6]
// [1,]  0.1  0.1  0.1  0.2  0.2  0.2
// [2,]  0.1  0.1  0.1  0.2  0.2  0.2
//
// xi_temp
//   
//   [,1] [,2] [,3] [,4]
// [1,]    3    3    3    3
// [2,]    4    4    4    4
//
// Ustar
//   
//      [,1] [,2] [,3] [,4]
// [1,]    1    1    2    2
// [2,]    1    1    2    2
// [3,]    1    1    2    2
// 