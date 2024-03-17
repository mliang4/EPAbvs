#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Mathematical constants
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

namespace help{

////////////////////////////////////
// Helper functions from PGBVS /////
////////////////////////////////////

// Generate exponential distribution random variates
double exprnd(double mu)
{
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

// Function a_n(x) defined in equations (12) and (13) of
// Bayesian inference for logistic models using Polya-Gamma latent variables
// Nicholas G. Polson, James G. Scott, Jesse Windle
// arXiv:1205.0310
//
// Also found in the PhD thesis of Windle (2013) in equations
// (2.14) and (2.15), page 24
double aterm(int n, double x, double t)
{
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }    
  return (double)exp(f);
}

// Generate inverse gaussian random variates
double randinvg(double mu)
{
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
  
  if(R::runif(0.0,1.0) > mu /(mu+out)) {    
    out = mu*mu / out; 
  }    
  return out;
}

// Sample truncated gamma random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables 
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
double truncgamma()
{
  double c = MATH_PI_2;
  double X, gX;
  
  bool done = false;
  while(!done)
  {
    X = help::exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
    
    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }
  
  return X;  
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
double tinvgauss(double z, double t)
{
  double X, u;
  double mu = 1.0/z;
  
  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma 
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / help::truncgamma();
      
      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }  
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = help::randinvg(mu);
    }
  }    
  return X;
}


// Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
double samplepg(double z)
{
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;
  
  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;
  
  // Compute p, q and the ratio q / (q + p)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
  double logK = (double)std::log(K);
  double Kt = K * t;
  double w = (double)std::sqrt(MATH_PI_2);
  
  double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q); 
  
  double u, X;
  
  // Main sampling loop; page 130 of the Windle PhD thesis
  while(1) 
  {
    // Step 1: Sample X ? g(x|z)
    u = R::runif(0.0,1.0);
    if(u < ratio) {
      // truncated exponential
      X = t + help::exprnd(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = help::tinvgauss(z, t);
    }
    
    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = help::aterm(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;
    
    while(1) 
    {
      Sn = Sn + asgn * help::aterm(i, X, t);
      
      // Accept if n is odd
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }
      
      // Return to step 1 if n is even
      if(even && (U > Sn)) {
        break;
      }
      
      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}

// Simulate MVT normal data
arma::mat mvrnormArma( int n, arma::vec mu, arma::mat sigma ) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn( n, ncols );
  return arma::repmat( mu, 1, n ).t()  + Y * arma::chol( sigma );
}

////////////////////////////////////////////
// Helper functions created for MCMC   /////
////////////////////////////////////////////

// Function :: remove the gaps within the integer vector (e.g. [0, 1, 3, 5] --> [0, 1, 2, 3])
// Also let the vector start from 0 
arma::vec tighten(arma::vec somevec) {
  arma::vec uni = unique(somevec); // unique elements in the vector
  int unisize = uni.size();
  int start = 0;
  
  for(int i = 0; i < unisize; ++i){
    uvec pos = find( somevec == uni(i) ); // position in the vector that equals to value uni[i]
    arma::vec rpl( pos.size() );
    rpl.fill( start + i );
    somevec( pos ) = rpl;
  }
  
  return somevec;
}

// Function :: make sX matrix of [sum_i sum_j ij] rows and [3P] columns.
arma::mat make_sX( arma::vec subject, arma::mat Xbar, arma::mat Ustar, arma::vec Ustar_dems, arma::mat xi ){
  int P = Ustar_dems.size() - 1;
  int obs = Ustar.n_rows;
  
  arma::mat sX( obs, 3*P, fill::zeros );              
  
  for( int p = 0; p < P; ++p ){
    arma::vec sum_p( obs, fill::zeros );
    
    // Get the range of columns corresponding to covariate p 
    arma::mat UstarXi = Ustar.cols( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1) ;
    
    for( int ij = 0; ij < obs ; ++ij  ){
      arma::mat UstarXi_ind = UstarXi.row( ij ) * xi.submat( subject[ ij ], Ustar_dems[ p ], subject[ ij ], Ustar_dems[ p + 1 ] - 1 ).t(); 
      sum_p[ ij ] +=  UstarXi_ind[ 0 ];
    }
    
    // Order the covariates 
    sX.col( 3*p ) = sum_p % Xbar.col( 2*p + 1 );
    sX.col( 3*p + 1 ) = Xbar.col( 2*p );
    sX.col( 3*p + 2 ) = Xbar.col( 2*p + 1 ); 
    
  }
  return sX; 
}


arma::mat make_dsX( arma::vec subject, arma::vec subject_dems, arma::mat Xbar, 
                    arma::mat Ustar, arma::vec Ustar_dems, arma::mat beta_temp  ){
  
  int sum_rp = Ustar.n_cols; // sum_rp
  int P = Ustar_dems.size() - 1; // number of smooth functions
  int obs = Ustar.n_rows; // number of total observations
  int N = subject_dems.size() - 1; // number of subjects
  
  arma::mat X( obs, P, fill::zeros ); // temporary space for all x_ijp
  arma::mat beta_n( N, P, fill::zeros ); // temporary space for all beta_n
  for( int p = 0; p < P; ++p ){
    X.col( p ) = Xbar.col( 2*p + 1 );
    beta_n.col( p ) = beta_temp.col( 3*p );
  }
  
  
  arma::mat BnX( obs, P, fill::zeros ); // make beta_n * X, where each row (ij) is beta_n^(i) * x_ij
  int count = 0;
  for( int i = 0; i < N; ++i ){
    for( int j = subject_dems[i]; j < subject_dems[i+1]; ++j){
      BnX.row( count ) = X.row( count ) % beta_n.row( i );
      count++; // for each row of X
    }
  }
  
  arma::mat dsX( obs, sum_rp, fill::zeros ); // make dot{sX}, where each column (p,r_p) is beta_np * x_ijp * Ustar_pr_p
  count = 0;
  for( int p = 0; p < P; ++p ){
    for( int r = Ustar_dems[p]; r < Ustar_dems[p+1]; ++r ){
      dsX.col( count ) = Ustar.col( count ) % BnX.col( p );
      count++; // for each column of Ustar
    }
  }
  
  return dsX; 
}


// Make a distance function parameterized by eta
arma::mat make_lambda_dist_p_logis( arma::mat dist, double eta_p ){
  arma::mat numer = exp( -dist + eta_p );
  arma::mat denom = 1 + exp( -dist + eta_p );
  arma::mat out = numer/denom;
  return out;
}

arma::mat make_lambda_dist_p_window( arma::mat dist, double eta_p ){
  int n = dist.n_cols;
  arma::mat out( n, n, fill::zeros );
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if( dist(i,j) < eta_p ){
        out(i,j) = 1;
      }else{
        out(i,j) = 0.01;
      }
    }
  }
  return out;
}

arma::mat make_lambda_dist_p_exp( arma::mat dist, double eta_p ){
  return exp( -eta_p * dist );
}

// Sample from an integrer vector 
int sample_cpp( arma::vec x ){
  // Calling sample()
  Function f( "sample" );
  IntegerVector sampled = f( x, Named( "size" ) = 1 );
  return sampled[ 0 ];
}

// Sample according to a vector of probability
int sample_prob_cpp( arma::vec x, arma::vec prob ){
  // Calling sample()
  Function f( "sample" );
  IntegerVector sampled = f( x, Named( "size" ) = 1, Named( "prob" ) = prob );
  return sampled[ 0 ];
}

// Function :: sampler for inverse-gamma distribution with shape (a) and scale (b) parameters
double rinvgamma_cpp( double a, double b ){
  double invg_return = 1 / Rcpp::rgamma( 1, a, 1/b )[0];
  return invg_return;
}

// Function :: Calculate normal log-density ( univariate )
double log_normal_cpp( double value, double mu, double sigma2 ){
  double log_normal = -0.50*log( 2*M_PI*sigma2 ) - 1/( 2*sigma2 )*pow( value - mu ,2 );
  
  // Return output
  return log_normal; 
}

// Function :: Calculate Beta log-density ( univariate )
double log_beta_cpp( double value, double a, double b ){
  double log_beta = (a-1) * log(value) + (b-1) * log(1-value) - lgamma(a) - lgamma(b) + lgamma(a+b);
  
  // Return output
  return log_beta; 
}

// Function :: Calculate beta-binomial log-density ( univariate )
double log_beta_binomial_cpp( double indicate, double a, double b ){
  double post_a = indicate + a;
  double post_b = 1 - indicate + b;
  double log_indicator = lgamma( post_a ) + lgamma( post_b ) - lgamma( post_a + post_b ) - ( lgamma( a ) + lgamma( b ) - lgamma( a + b ) );
  
  // Return output
  return log_indicator ;
}

// Function :: Calculate multivariate normal log-density only with diagonal covariance matrix
double log_mvnorm_cpp( arma::vec value, arma::vec mu, arma::vec diag_sig2 ){
  int n = value.size();
  double d = 0;
  
  for (int i = 0; i < n; ++i){
    d += help::log_normal_cpp( value[i], mu[i], diag_sig2[i] );
  }
  
  return d;
}

// Function :: Log-likelihood for a single subject: h_ij = k_ij/w_ij and h ~ N( psi, Omega), k_ij = y_ij - 1/2
double log_like_logit_single( 
    int subj_id, // The subject that we are looking at, 
    arma::vec subject,
    arma::vec Y, 
    arma::mat Xbar, 
    arma::mat Ustar, 
    arma::vec Ustar_dems, 
    arma::mat beta_temp, 
    arma::mat xi_temp
){
  
  int obs = subject.size();
  int S = beta_temp.n_cols; // 3P
  
  // Make sX from Ustar, xi, Xbar
  arma::mat sX( obs, S, fill::zeros ); 
  sX = help::make_sX( subject, Xbar, Ustar, Ustar_dems, xi_temp );
  
  double log_like = 0;
  
  // Find rows that belongs to subject i
  uvec pos = find( subject == subj_id );
  int nrow_i = pos.size(); // The number of rows for i
  
  // Make psi and calculate the log-likelihood contribution 
  for(int i = 0; i < nrow_i; ++i){
    arma::mat psi_val = sX.row( pos[i] )*beta_temp.row( subj_id ).t();
    log_like += Y[ pos[i] ]*psi_val[0] - log( 1 + exp( psi_val[0] ) );
    // log_like += -0.50*log( 2*M_PI*Winv ) - 1/( 2*Winv )*pow( H[ pos[i] ] - psi_val[0], 2 );
  }
  
  // Return output
  return log_like;
} 


// Function :: log probability of subject being in a given cluster.
double log_prob_clust_single(
    int i,                    // chosen subject
    int p,                    // current smooth function
    arma::mat s_temp,         // a given cluster membership indicator matrix
    arma::vec alpha,          // Vector of current alphas. Indexed by p
    arma::vec delta,
    arma::uvec sigma,
    arma::mat lambda_dist_p
){
  
  // Return the probability of a subject belongs to the cluster given the cluster membership indicator matrix
  // i.e. p( i \in s_temp[i,p] | ... )
  
  arma::vec sp = s_temp.col( p ); // Current cluster membership indicator for p
  arma::vec lambda_dist_pi = lambda_dist_p.col( i ); // distance between i and other subjects
  double a_p = alpha( p );
  double d_p = delta( p );
  
  // Find the position of i in sigma
  uvec t = find( sigma == i ); 
  int sigma_i = t[ 0 ]; // sigma_i (cpp index) = i - 1 (number of items before i)
  
  double p_return = 1; // Initialize the log probability for return
  
  if ( sigma_i > 0 ){ 
    // Calculate probability if sigma_i is not the first one in the permutation. 
    // If sigma_i is the first item in permutation then return log(1) = 0
    
    // Get sigma_1 to sigma_i-1, if t is not the first element in sigma
    uvec sigma_vec = sigma( span( 0, sigma_i - 1) );
    arma::vec unique_clust = unique( sp( sigma_vec ) );
    int q = unique_clust.size(); // q_p,i-1
    
    // Find elements with the same cluster indicator as i
    uvec j = find( sp( sigma_vec ) == sp( i ) );
    
    if ( j.size() > 0 ){ // If there is a match in the cluster membership
      p_return = ( sigma_i - d_p*q ) / ( a_p + sigma_i ) * sum( lambda_dist_pi( sigma_vec(j) ) ) / sum( lambda_dist_pi( sigma_vec ) );
    }else{
      p_return = ( a_p + d_p*q ) / ( a_p + sigma_i ); 
    }
    
    p_return += 1e-30;
    
  } // End of if( sigma_i > 0 )
  
  return log(p_return);
  
}

// Function :: Rescale beta_temp and xi_temp
List rescaled( arma::mat beta_temp, arma::mat xi_temp, arma::vec Ustar_dems ){
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



// FUNCTION :: Update W for each ij

arma::vec update_w( 
    arma::mat Xbar, 
    arma::mat Ustar, 
    arma::vec Ustar_dems, 
    arma::mat beta_temp, // current slice of beta
    arma::mat xi_temp, // current slice of xi
    arma::vec subject 
){
  int obs = subject.size();
  int S = beta_temp.n_cols; // 3P
  
  // Make a home for w updates
  arma::vec updated_w( obs, fill::zeros ); 
  
  // Make sX from Ustar, xi, Xbar
  arma::mat sX( obs, S, fill::zeros ); 
  sX = help::make_sX( subject, Xbar, Ustar, Ustar_dems, xi_temp );
  
  // Update each w individually
  for( int j = 0; j < obs; ++j ){
    int sub = subject[ j ];
    arma::mat psi_j(1,1);
    psi_j = sX.row( j ) * beta_temp.row( sub ).t();
    updated_w[ j ] = help::samplepg( psi_j[ 0 ] );
  }
  
  return updated_w;
  
}

// Function :: Update lambda2_pk
arma::vec update_lambda2pk(
    arma::vec beta_pk,
    arma::vec nu_lambda_pk,
    double tau2_pk
){
  
  arma::vec lambda2pk_return( 3, fill::zeros );
  
  if( sum(beta_pk) == 0 ){
    
    for( int nlm = 0; nlm < 3; nlm++ ){
      lambda2pk_return[ nlm ] = rinvgamma_cpp( 0.5, 1/nu_lambda_pk[ nlm ] );
    }
    
  }else{
    
    for( int nlm = 0; nlm < 3; nlm++ ){
      lambda2pk_return[ nlm ] = rinvgamma_cpp( 1, 1 / nu_lambda_pk[ nlm ] + 0.5 * pow( beta_pk[ nlm ], 2 ) / tau2_pk );
    }
    
  } 
  
  return lambda2pk_return;
}

int getgamma(arma::mat s_temp, arma::mat gamma_temp, int p, int k) {
  
  uvec temp_pos = find( s_temp.col(p) == k );
  int temp_i = temp_pos[ 0 ]; // one of the subjects in cluster uk
  int temp_g = gamma_temp(temp_i, p);
  
  return temp_g;
}

// Function :: Put all zero-coefficient subjects in one cluster
arma::vec relabel_zero(arma::vec sp, arma::vec gp) {
  
  uvec tmp1 = find( gp == 0 );
  sp( tmp1 ).fill( max(sp) + 1 );
  
  return sp;
}





/// PY helper functions ///


// Function :: calculate the log of the bracket function in proposition 9 of Pitman 1995
double log_bracket(double x, int m, double a){
  double res = 0;
  if(m > 0){
    for(int mm = 0; mm < m; mm++){
      res += log( x + mm*a );
    }
  }
  return res;
}

// Function :: calculate the log of the exchangeable probability function (EPF) (for 1 smooth function) by proposition 9 in Pitman 1995
double log_epf(arma::vec sp, double alpha_p, double delta_p){
  
  int N = sp.size();
  arma::vec usp = unique( sp );
  int K = usp.size(); // number of unique clusters
  
  double res = 0;
  for(int i = 0; i < K; i++){
    uvec fvec = find(sp == usp[i]);
    int ni = fvec.size();
    res += log_bracket( 1-delta_p, ni-1, 1 );
  }
  
  res += log_bracket( alpha_p+delta_p, K-1, delta_p);
  res -= log_bracket( alpha_p+1, N-1, 1);
  
  return res;
  
}

arma::vec get_spmi(arma::vec sp, int i){
  
  arma::vec spmi(sp.size()-1, fill::zeros);
  int qp = sp.size();
  int j = 0;
  
  for(int k = 0; k < qp; k++){
    if(k != i){
      spmi[j++] = sp[k];
    }
  }
  
  return spmi;
  
}

/// PY helper functions ///



////////////////////////////////////////////////
//              Main functions             /////
////////////////////////////////////////////////

List between_step(
    int i, // chosen subject
    int p, // current smooth function
    arma::vec subject,
    arma::vec Y,
    arma::mat Xbar,
    arma::mat Ustar,
    arma::vec Ustar_dems,
    arma::mat beta_temp, 
    arma::mat xi_temp,
    arma::mat s_temp, // current cluster membership indicator
    arma::mat gamma_temp,
    arma::vec mu,
    arma::mat lambda2_temp,
    arma::vec tau2_temp,
    double a_gamma,
    double b_gamma
){
  
  
  // perform a between step for everybody that are in the same cluster as subject i
  arma::vec sp = s_temp.col( p ); // Current cluster membership indicator
  uvec member = find( sp == sp( i ) );
  int spk_size = member.size();

  arma::mat beta_proposal = beta_temp; // Create a temporary space for beta_proposal (N-by-3P matrix)
  
  arma::rowvec beta_temp_ip = beta_temp( i, span( 3*p, 3*p+2 ) ); // Current beta_p (beta_np, beta_lp, beta_mp) for subject i
  arma::rowvec beta_proposal_ip( 3, fill::zeros ); 
  // Create a temporary space for beta_proposal (beta_np, beta_lp, beta_mp) proposal for subject i
  
  double tau2_p = tau2_temp[ p ];
  arma::vec zero_vec3( 3, fill::zeros ); // Make a 3-d vector of 0s
  arma::vec lam2_prime(3, fill::zeros ); // Make a temporary space for the diagonal elements of the lambda2_prime in the add step
  arma::vec add_var_diag(3, fill::zeros ); // Make a temporary space for the diagonal elements of the prior variance in the add step
  
  arma::vec lambda2_pk = lambda2_temp( i, span( 3*p, 3*p+2 ) ).t();
  arma::vec var_diag = tau2_p * lambda2_pk; // diagonal elements of the prior variance var(beta) in the add/delete step
  arma::mat one_mat3 = eye( 3, 3 ); // Make a 3-d identity matrix for proposal variance
  
  
  if( gamma_temp( i, p ) == 0 ){
    
    // Add
    
    beta_proposal_ip = help::mvrnormArma( 1, zero_vec3, 5*one_mat3 );
    for(int j = 0; j < spk_size; j++){
      beta_proposal( member(j), span( 3*p, 3*p+2 ) ) = beta_proposal_ip;
    }
    
    double r = help::log_mvnorm_cpp(beta_proposal_ip.t(), zero_vec3, var_diag)
      + help::log_beta_binomial_cpp(1, a_gamma, b_gamma)
      - help::log_beta_binomial_cpp(0, a_gamma, b_gamma);
      
    for(int j = 0; j < spk_size; j++){
      r += help::log_like_logit_single(member(j), subject, Y, Xbar, Ustar, Ustar_dems, beta_proposal, xi_temp) -
        help::log_like_logit_single(member(j), subject, Y, Xbar, Ustar, Ustar_dems, beta_temp, xi_temp); 
    }
      
    
      
    double a  = log( runif( 1 )[ 0 ] );
      
    // Determine acceptance
    if( a < r ){
      beta_temp = beta_proposal;
      for(int j = 0; j < spk_size; j++){
        s_temp( member(j), p ) = max(sp) + 1; // assign this subject to a new cluster
        gamma_temp( member(j), p ) = 1;
      }
    }
      
  } else {
    
    // Delete
    
    beta_proposal_ip = zero_vec3.t();
    for(int j = 0; j < spk_size; j++){
      beta_proposal( member(j), span( 3*p, 3*p+2 ) ) = beta_proposal_ip;
    }
    
    double r = help::log_beta_binomial_cpp(0, a_gamma, b_gamma)
      - help::log_beta_binomial_cpp(1, a_gamma, b_gamma)
      - help::log_mvnorm_cpp(beta_temp_ip.t(), zero_vec3, var_diag);
      
    for(int j = 0; j < spk_size; j++){
      r += help::log_like_logit_single(member(j), subject, Y, Xbar, Ustar, Ustar_dems, beta_proposal, xi_temp) -
        help::log_like_logit_single(member(j), subject, Y, Xbar, Ustar, Ustar_dems, beta_temp, xi_temp); 
    }
      
    double a  = log( runif( 1 )[ 0 ] );
    
    // Determine acceptance
    if( a < r ){
      beta_temp = beta_proposal;
      for(int j = 0; j < spk_size; j++){
        s_temp( member(j), p ) = max(sp) + 1; // assign this subject to a new cluster
        gamma_temp( member(j), p ) = 0;
      }
    }
  }
  
  List between_clust( 3 );
  between_clust[ 0 ] = beta_temp;
  between_clust[ 1 ] = s_temp;
  between_clust[ 2 ] = gamma_temp;
  
  return between_clust;
  
} 

arma::mat within_beta_pk(
    int p,
    int k, // the cluster indicator, must be a non-trivial cluster (i.e. greater than 0)
    arma::vec subject,
    arma::vec subject_dems,
    arma::vec Y,
    arma::mat Xbar, 
    arma::mat Ustar, 
    arma::vec Ustar_dems, 
    arma::mat beta_temp,
    arma::mat xi_temp,
    arma::mat s_temp,
    arma::vec w,
    arma::mat lambda2_temp,
    arma::vec tau2_temp
){
  
  // Find subjects and the number of subjects belong to kth cluster
  arma::vec sp = s_temp.col( p );
  uvec pos = find( sp == k ); // The vector of all subjects in cluster k
  // Rcout << "The vector of all subjects in cluster k" << pos << std::endl;
  int spk_size = pos.size();
  int i = pos[ 0 ]; // one of the subjects in cluster k
  
  int obs = subject.size();
  int S = beta_temp.n_cols; // 3P
  double tau2_p = tau2_temp[ p ];
  arma::vec lambda2_pi = lambda2_temp( i, span( 3*p, 3*p+2 ) ).t();
  
  // Make sX from Ustar, xi, Xbar
  arma::mat sX( obs, S, fill::zeros );
  sX = help::make_sX( subject, Xbar, Ustar, Ustar_dems, xi_temp ); // sum_ij-by-3P
  
  // Make sX_p, sX_mp (i.e. sX_{-p})
  arma::mat sX_p( obs, 3, fill::zeros );
  arma::mat sX_mp( obs, S-3, fill::zeros );
  
  sX_p = sX.cols( 3*p, 3*p+2 );
  
  if(p > 0){
    sX_mp.cols( 0, 3*p-1 ) = sX.cols( 0, 3*p-1 );
  }
  if(3*p+3 < S){
    sX_mp.cols( 3*p, S-4 ) = sX.cols( 3*p+3, S-1 );
  }
  
  // Make h_ij = k_ij/w_ij, k_ij = y_ij - 1/2
  arma::vec H( obs, fill::zeros );
  H = (Y - 0.5)/w;
  
  // Make the prior variance matrix
  arma::vec tau2_three (3);
  tau2_three.fill(tau2_p);
  arma::vec s2_diag = tau2_three % lambda2_pi;
  // s2_diag += 0.000000000001;
  
  arma::mat Sigma_star_inv = diagmat( 1 / s2_diag );
  
  arma::mat Sigma_bkp_star = Sigma_star_inv; // Starting point of calculating Sigma_{beta^k_p}*
  
  arma::vec mu_bkp_star( 3, fill::zeros );
  
  for(int idx = 0; idx < spk_size; ++idx){
    
    i = pos[ idx ]; // subject id
    
    // Make beta_ip, beta_imp (i.e. beta^i_{-p})
    arma::vec beta_i = beta_temp.row( i ).t();
    arma::vec beta_ip = beta_i( span( 3*p, 3*p+2 ) ); // beta_i^p is the same across p, but beta_i^{-p} is different
    arma::vec beta_imp( S-3, fill::zeros );
    if(p > 0){
      beta_imp( span( 0, 3*p-1 ) ) = beta_i( span( 0, 3*p-1 ) );
    }
    if(3*p+3 < S){
      beta_imp( span( 3*p, S-4 ) ) = beta_i( span( 3*p+3, S-1 ) );
    }
    
    arma::mat sX_ip = sX_p.rows( subject_dems[i], subject_dems[i + 1] - 1 ); // sum_j-by-3
    arma::mat sX_imp = sX_mp.rows( subject_dems[i], subject_dems[i + 1] - 1 ); // sum_j-by-(3P-3)
    arma::vec sZ_ip = H( span( subject_dems[i], subject_dems[i + 1] - 1 ) ); // sum_j
    sZ_ip -= sX_imp * beta_imp; //
    
    arma::vec omega_i = w( span( subject_dems[i], subject_dems[i + 1] - 1 ) );
    arma::mat Omega_i = diagmat( omega_i );
    
    Sigma_bkp_star += sX_ip.t() * Omega_i * sX_ip;
    mu_bkp_star += sX_ip.t() * Omega_i * sZ_ip;
  }
  
  Sigma_bkp_star = inv( Sigma_bkp_star ); // inverse of a 3-by-3 matrix, shouldn't take too long
  mu_bkp_star = Sigma_bkp_star * mu_bkp_star;
  
  
  arma::vec beta_new = help::mvrnormArma( 1, mu_bkp_star, Sigma_bkp_star ).t();
  
  for(int idx = 0; idx < spk_size; ++idx){
    
    int i = pos[ idx ]; // subject id
    beta_temp( i, span(3*p, 3*p+2) ) = beta_new.t();
    
  }
  
  return beta_temp;
}

arma::mat within_xi_pk(
    int p,
    int k, // the cluster indicator
    arma::vec subject,
    arma::vec subject_dems,
    arma::vec Y,
    arma::mat Xbar, 
    arma::mat Ustar, 
    arma::vec Ustar_dems, 
    arma::mat beta_temp,
    arma::mat xi_temp,
    arma::mat s_temp,
    arma::vec w,
    arma::vec mu
){
  
  int obs = subject.size(); // sum_ij
  int S = Ustar.n_cols; // sum_rp
  int rp = Ustar_dems[p+1] - Ustar_dems[p];
  int N = subject_dems.size() - 1;
  int P = Ustar_dems.size() - 1;
  
  // Make dot_sX
  arma::mat dsX( obs, S, fill::zeros );
  dsX = help::make_dsX( subject, subject_dems, Xbar, Ustar, Ustar_dems, beta_temp );
  
  // make dsX_(m)p, xi_(m)p
  arma::mat dsX_mp( obs, S - rp, fill::zeros );
  arma::mat xi_mp( N, S - rp, fill::zeros );
  
  if( Ustar_dems[p] > 0 ){
    dsX_mp.cols( 0, Ustar_dems[p] - 1 ) = dsX.cols(     0, Ustar_dems[p] - 1 );
    xi_mp.cols(  0, Ustar_dems[p] - 1 ) = xi_temp.cols( 0, Ustar_dems[p] - 1 );
  }
  
  if( Ustar_dems[p] < S - rp ){
    dsX_mp.cols( Ustar_dems[p], S - rp - 1 ) = dsX.cols(    Ustar_dems[ p + 1 ], S - 1 );
    xi_mp.cols(  Ustar_dems[p], S - rp - 1 ) = xi_temp.cols( Ustar_dems[ p + 1 ], S - 1 );
  }
  
  arma::mat dsX_p = dsX.cols(    Ustar_dems[p], Ustar_dems[ p + 1 ] - 1 );
  arma::mat xi_p = xi_temp.cols( Ustar_dems[p], Ustar_dems[ p + 1 ] - 1 );
  
  // Make beta_{-n} (with beta_l and beta_m only)
  arma::mat beta_mn( N, 2*P, fill::zeros );
  for( int p = 0; p < P; ++p ){
    beta_mn.col( 2 * p ) = beta_temp.col( 3 * p + 1 );
    beta_mn.col( 2 * p + 1 ) = beta_temp.col( 3 * p + 2 );
  }
  
  // Make h_ij = k_ij/w_ij, k_ij = y_ij - 1/2
  arma::vec H( obs, fill::zeros );
  H = (Y - 0.5)/w;
  
  // Find subjects and the number of subjects belong to kth cluster
  arma::vec sp = s_temp.col( p );
  uvec pos = find( sp == k ); // The vector of all subjects in cluster k
  int spk_size = pos.size();
  
  // Make the prior mean
  
  arma::vec mu_p = mu( span( Ustar_dems[p], Ustar_dems[ p + 1 ] - 1 ) );
  
  arma::mat Sigma_xkp_star = eye( rp, rp ); // Use the identity matrix as the Starting point of calculating Sigma_xi^k_p_star
  
  arma::vec mu_xkp_star( rp, fill::zeros );
  
  for(int idx = 0; idx < spk_size; ++idx){
    
    int i = pos[ idx ]; // subject id
    
    arma::vec xi_imp = xi_mp.row( i ).t();
    arma::vec beta_imn = beta_mn.row( i ).t();
    
    arma::mat dsX_ip = dsX_p.rows( subject_dems[i], subject_dems[i + 1] - 1 );
    arma::mat dsX_imp = dsX_mp.rows( subject_dems[i], subject_dems[i + 1] - 1 );
    arma::mat ddsX_ip = Xbar.rows( subject_dems[i], subject_dems[i + 1] - 1 );
    
    arma::vec dsZ_ip = H( span( subject_dems[i], subject_dems[i + 1] - 1 ) );
    
    dsZ_ip = dsZ_ip - dsX_imp * xi_imp - ddsX_ip * beta_imn;
    
    arma::vec omega_i =  w( span( subject_dems[i], subject_dems[i + 1] - 1 ) );
    arma::mat Omega_i = diagmat( omega_i );
    Sigma_xkp_star += dsX_ip.t() * Omega_i * dsX_ip;
    mu_xkp_star += dsX_ip.t() * Omega_i * dsZ_ip;
  }
  
  Sigma_xkp_star = inv( Sigma_xkp_star );
  mu_xkp_star = Sigma_xkp_star * mu_xkp_star;
  
  
  arma::vec xi_new = help::mvrnormArma( 1, mu_xkp_star, Sigma_xkp_star ).t();
  
  for(int idx = 0; idx < spk_size; ++idx){
    
    int i = pos[ idx ]; // subject id
    xi_temp( i, span( Ustar_dems[p], Ustar_dems[ p + 1 ] - 1 ) ) = xi_new.t();
    
  }
  
  return xi_temp;
}

// Update xi_pk for subjects in the excluded clusters, xi_p0|. ~ N(mu_p, I)
arma::mat within_xi_p0(
    int p,
    int k, // the cluster indicator
    arma::vec Ustar_dems, 
    arma::mat xi_temp,
    arma::mat s_temp,
    arma::vec mu
){
  
  int rp = Ustar_dems[p+1] - Ustar_dems[p];
  
  // Find subjects and the number of subjects belong to trivial cluster
  arma::vec sp = s_temp.col( p );
  uvec pos = find( sp == k ); // The vector of all subjects in cluster
  int sp0_size = pos.size();
  
  arma::vec mu_p = mu( span( Ustar_dems[p], Ustar_dems[ p + 1 ] - 1 ) ); // prior mean
  arma::mat I_rp = eye( rp, rp ); // identity matrix
  
  arma::vec xi_new = help::mvrnormArma( 1, mu_p, I_rp ).t();
  
  for(int idx = 0; idx < sp0_size; ++idx){
    
    int i = pos[ idx ]; // subject id
    xi_temp( i, span( Ustar_dems[p], Ustar_dems[ p + 1 ] - 1 ) ) = xi_new.t();
    
  }
  
  return xi_temp;
}

// Update the cluster membership according to EPA distribution
List update_spi_EPA(
    int i, // chosen subject 
    int p, // current smooth function
    arma::vec between_covas,
    arma::vec subject,
    arma::vec subject_dems,
    arma::vec Y,
    arma::mat Xbar,
    arma::mat Ustar,
    arma::vec Ustar_dems,
    arma::mat beta_temp, 
    arma::mat xi_temp,
    arma::mat s_temp,             // current cluster membership indicator
    arma::mat gamma_temp,
    arma::vec w,
    arma::vec alpha,              // Vector of current alphas. Indexed by p
    arma::vec delta,
    arma::uvec sigma,
    arma::vec mu,
    arma::mat lambda2_temp,
    arma::vec tau2_temp,
    arma::mat lambda_dist_p
){
  
  int N = sigma.size(); // total number of subjects
  int rp = Ustar_dems[ p + 1 ] - Ustar_dems[ p ]; // Dimension of xi_temp_ip
  arma::vec mu_beta( 3, fill::zeros ); // Make a 3-d vector of 0s
  arma::vec mu_xi = mu( span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) ); // get mu_rp
  arma::mat one_beta = eye( 3, 3 ); // Make a 3-d identity matrix
  arma::mat one_xi = eye( rp, rp ); // Make a rp-d identity matrix
  arma::rowvec ltmp( 3, fill::zeros ); // Make a temporary space to put the lambdas
  
  // Get the number and the indexes of all existing clusters
  arma::vec sp = s_temp.col( p );              // the cluster indicator s_p
  arma::uvec uclust_pos = find_unique( sp );   // positions for the unique clusters, one for each cluster
  int q_i = uclust_pos.size();                 // number of existing unique clusters
  
  // Make a copy of the parameters (beta_temp/xi_temp should be unchanged until the cluster to move to is decided)
  arma::mat beta_temp2 = beta_temp;
  arma::mat xi_temp2 = xi_temp;
  arma::mat s_temp2 = s_temp;
  arma::mat gamma_temp2 = gamma_temp;
  arma::mat lambda2_temp2 = lambda2_temp;
  
  // make space to save the change for moving to each cluster
  arma::mat beta_candidate( q_i + 2, 3, fill::zeros );
  arma::mat xi_candidate( q_i + 2, rp, fill::zeros );
  arma::vec s_candidate( q_i + 2, fill::zeros );
  arma::vec gamma_candidate( q_i + 2, fill::zeros );
  arma::mat lambda2_candidate( q_i + 2, 3, fill::zeros );
  arma::vec v_pi( q_i + 2, fill::zeros ); // Make space to store each probability of moving i to a updated cluster

  // Fill the candidates with parameters from existing clusters and Propose 1 non-zero new cluster from N(0,5) and 1 zero new cluster
  for(int k = 0; k < q_i + 2; ++k){ // propose a move to cluster k
    
    double v_pik = 0;
    
    if(k == q_i + 1){ // moving to a new beta!=0 cluster
      
      s_temp2( i, p ) = max( sp ) + 2; // Assume the indicator of the new cluster to be greater than all indicators in s_p 
      beta_temp2( i, span( 3*p, 3*p+2 ) ) = help::mvrnormArma( 1, mu_beta, 5*one_beta );
      xi_temp2( i, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) ) = help::mvrnormArma( 1, mu_xi, one_xi );
      gamma_temp2( i, p ) = 1;
      ltmp = Rcpp::rgamma(3,5,1);
      lambda2_temp2( i, span( 3*p, 3*p+2 ) ) = ltmp; // sample something greater than 1
      
    }else if( k == q_i ){ // moving to a new beta=0 cluster
      
      s_temp2( i, p ) = max( sp ) + 1; // Assume the indicator of the new cluster to be greater than all indicators in s_p 
      beta_temp2( i, span( 3*p, 3*p+2 ) ).zeros();
      xi_temp2( i, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) ) = help::mvrnormArma( 1, mu_xi, one_xi );
      gamma_temp2( i, p ) = 0;
      ltmp = Rcpp::rgamma(3,1,0.2);
      lambda2_temp2( i, span( 3*p, 3*p+2 ) ) = ltmp; // sample something close to 0
      
    }else{ // if moving to one of the previous non-trivial clusters
      
      int j = uclust_pos( k ); // find a subject j in cluster k and mark down its parameters
      beta_temp2( i, span( 3*p, 3*p+2 ) ) = beta_temp( j, span( 3*p, 3*p+2 ) );
      xi_temp2( i, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) ) = xi_temp( j, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) );
      s_temp2( i, p ) = s_temp( j, p ); 
      gamma_temp2( i, p ) = gamma_temp( j, p );
      lambda2_temp2( i, span( 3*p, 3*p+2 ) ) = lambda2_temp( j, span( 3*p, 3*p+2 ) );
      
    } // end of if-elseif-else (k = 0,1,...,q,q+1) loop
    
    // Rcout << "move to cluster: " << s_temp2( i, p ) << " with beta: " << beta_temp2( i, span( 3*p, 3*p+2 ) ) << endl;
    
    // Calculate full conditionals by summing the log partition probability and the log likelihood with new proposed parameters
    
    // Calculate the sum of the log partition probability for all j' >= i, j' < i can be cancelled out 
    // Find the position of i in sigma
    uvec tmp1 = find( sigma == i );
    int sigma_i = tmp1[0];
    uvec sigma_vec2 = sigma( span( sigma_i, N - 1 ) );     // Get sigma_i to sigma_N, if t is not the first element in sigma
    int n_vec2 = sigma_vec2.size(); // the number of elements in sigma_vec2
    
    for(int idx = 0; idx < n_vec2; ++idx){
      int j_prime = sigma_vec2[idx]; // j_prime is the original index
      v_pik += help::log_prob_clust_single( j_prime, p, s_temp2, alpha, delta, sigma, lambda_dist_p );
    }
    
    // Rcout << "log partition probability: " << v_pik << endl;
    
    // add likelihood
    v_pik += help::log_like_logit_single( i, subject, Y, Xbar, Ustar, Ustar_dems, beta_temp2, xi_temp2 );
    v_pi[ k ] = v_pik;
    
    // Rcout << "log loglik: " << help::log_like_logit_single( i, subject, Y, Xbar, Ustar, Ustar_dems, beta_temp2, xi_temp2 ) << endl;
    
    // save parameters candidates
    beta_candidate.row( k ) = beta_temp2( i, span( 3*p, 3*p+2 ) );
    xi_candidate.row( k ) = xi_temp2( i, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) );
    s_candidate[ k ] = s_temp2( i, p );
    gamma_candidate[ k ] = gamma_temp2( i, p );
    lambda2_candidate.row( k ) = lambda2_temp2( i, span( 3*p, 3*p+2 ) );
    
  } // End for k = 0 to q+1
  
  // Sample the new cluster from all the full conditionals
  v_pi = exp(v_pi - max(v_pi));
  
  bool zero_cluster = any(between_covas == p);
  if (!zero_cluster){
    v_pi[q_i] = 0;
  }
  
  int new_clust = help::sample_prob_cpp( s_candidate, v_pi );
  
  int index = 0;
  uvec indextemp = find( s_candidate == new_clust ); // which s_candidate is selected?
  index = indextemp[0];
  
  beta_temp( i, span( 3*p, 3*p+2 ) ) = beta_candidate.row( index );
  xi_temp( i, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) ) = xi_candidate.row( index );
  s_temp( i, p ) = new_clust;
  gamma_temp( i, p ) = gamma_candidate( index );
  lambda2_temp( i, span( 3*p, 3*p+2 ) ) = lambda2_candidate.row( index );
  
  // Rcout << "moving to cluster: " << new_clust << endl;
  // Rcout << "new parameter: " << beta_temp( i, span( 3*p, 3*p+2 ) ) << endl;
  
  List vpi_list( 5 );
  vpi_list[ 0 ] = beta_temp;
  vpi_list[ 1 ] = xi_temp;
  vpi_list[ 2 ] = s_temp;
  vpi_list[ 3 ] = gamma_temp;
  vpi_list[ 4 ] = lambda2_temp;
  
  return vpi_list;
  
}

arma::vec update_mu_rp( int p, arma::vec Ustar_dems, arma::mat xi_temp, arma::mat s_temp, arma::vec mu ){
  
  uvec uclust_ind = find_unique( s_temp.col(p) ); // Find one subject from each of the clusters
  double exp_temp = 0; // the amount within the exponential
  double prob = 0; // p(mu_pr = 1 | cdot)
  for( int r = Ustar_dems[ p ]; r < Ustar_dems[ p+1 ]; ++r ){
    
    arma::vec xi_i = xi_temp.col( r ); // get xi_pr for all subject i
    
    exp_temp = -2 * sum( xi_i( uclust_ind ) );
    
    prob = 1 / ( 1 + exp( exp_temp ) );
    
    int mu_ind = rbinom( 1, 1, prob )[0];
    
    if( mu_ind == 1){
      mu[ r ] = 1;
    }else{
      mu[ r ] = -1;
    }
    
  } // end for
  
  return mu;
  
}

arma::mat update_lambda( arma::mat beta_temp, arma::mat s_temp, arma::mat gamma_temp, arma::vec tau2_temp, arma::mat nu_lambda  ){
  
  int P = tau2_temp.size();
  int N = beta_temp.n_rows;
  arma::mat lambda2_return( N, 3*P, fill::zeros );
  // double ltemp = 0;
  double aterm = 0;
  double bterm = 0;
  arma::vec lambda2pk_temp( 3 ); // a temporary space to put lambda2_npk, lambda2_lpk, lambda2_mpk
  
  for( int p = 0; p < P; ++p ){
    
    arma::vec uclust = unique( s_temp.col(p) ); // unique cluters
    
    int n_unique = uclust.size();
    
    for( int k = 0; k < n_unique; ++k ){
      
      int uk = uclust[k];
      int included = help::getgamma( s_temp, gamma_temp, p, uk);
      
      if( included == 1 ){ // If in included cluster
        
        uvec pos = find( s_temp.col(p) == uk ); // pos is the set of subjects in cluster k
        int spk_size = pos.size();
        int idx = pos[ 0 ]; // get one of the subjects in cluster k
        
        for( int nlm = 0; nlm < 3; ++nlm ){ // 0: _n term, 1: _l term, 2: _m term
          aterm = 1;
          bterm = 1/nu_lambda( idx, 3*p + nlm ) +  0.5*pow( beta_temp( idx, 3*p + nlm ), 2 ) / tau2_temp[ p ];
          lambda2pk_temp[ nlm ] = rinvgamma_cpp(aterm, bterm);
        } // end for nlm
        
        for(int idx2 = 0; idx2 < spk_size; ++idx2){ // put updated lambda2 back to each subject within cluster k
          int i = pos[ idx2 ]; // subject id
          lambda2_return( i, span(3*p, 3*p+2) ) = lambda2pk_temp.t();
        } // end for idx
        
      }else{ // If in excluded cluster
        
        uvec pos = find( s_temp.col(p) == uk ); // pos is the set of subjects in cluster k
        int spk_size = pos.size();
        int idx = pos[ 0 ]; // get one of the subjects in cluster k
        
        for( int nlm = 0; nlm < 3; ++nlm ){ // 0: _n term, 1: _l term, 2: _m term
          aterm = 0.5;
          bterm = 1/nu_lambda( idx, 3*p + nlm );
          lambda2pk_temp[ nlm ] = rinvgamma_cpp(aterm, bterm);
        } // end for nlm
        
        for(int idx2 = 0; idx2 < spk_size; ++idx2){ // put updated lambda2 back to each subject within cluster k
          int i = pos[ idx2 ]; // subject id
          lambda2_return( i, span(3*p, 3*p+2) ) = lambda2pk_temp.t();
        } // end for idx
        
      } // end if-else in trivial cluster
      
    } // end for k
    
  } // end for p
  
  return lambda2_return;
  
}


arma::mat update_nu_lambda( arma::mat s_temp, arma::mat lambda2_temp  ){
  
  int P = s_temp.n_cols;
  int N = lambda2_temp.n_rows;
  arma::mat nu_lambda_return( N, 3*P, fill::zeros );
  // double ltemp = 0;
  double bterm = 0;
  arma::vec nu_pk_temp( 3 ); // a temporary space to put nu_lambda2_npk, nu_lambda2_lpk, nu_lambda2_mpk
  
  for( int p = 0; p < P; ++p ){
    
    arma::vec uclust = unique( s_temp.col( p ) ); // unique clusters
    int n_unique = uclust.size();
    
    for( int k = 0; k < n_unique; ++k ){
      
      uvec pos = find( s_temp.col( p ) == uclust[ k ] ); // pos is the set of subjects in cluster k
      int spk_size = pos.size();
      int idx = pos[ 0 ]; // get one of the subjects in cluster k
      
      for( int nlm = 0; nlm < 3; ++nlm ){
        bterm = 1 + 1 / lambda2_temp( idx, 3*p + nlm );
        nu_pk_temp[ nlm ] = rinvgamma_cpp(1, bterm);
      } // end for nlm
      
      for( int idx2 = 0; idx2 < spk_size; ++idx2 ){
        
        int i = pos[ idx2 ]; // subject id
        nu_lambda_return( i, span(3*p, 3*p+2) ) = nu_pk_temp.t(); // put the sampled value in all subjects within the same cluster
        
      } // end for idx
      
    } // end for k
    
  } // end for p
  
  return nu_lambda_return;
  
}

arma::vec update_tau( arma::mat beta_temp,
                      arma::mat s_temp,
                      arma::mat gamma_temp,
                      arma::mat lambda2_temp,
                      arma::vec nu_tau){
  
  int P = nu_tau.size();
  
  arma::vec tau2_return( P, fill::zeros );
  
  double aterm = 0;
  double bterm = 0;
  // double ltemp = 0;
  
  for( int p = 0; p < P; ++p ){
    
    arma::vec uclust = unique( s_temp.col(p) ); // unique cluters
    int n_unique = uclust.size();
    aterm = 0.5;
    bterm = 1/nu_tau[ p ]; // first start with 1 over nu_tau_p for the b term.
    
    for( int k = 0; k < n_unique; ++k ){
      
      int uk = uclust[k];
      int included = help::getgamma( s_temp, gamma_temp, p, uk);
      
      if( included == 1 ){ // If excluded cluster add nothing to aterm and bterm
        
        uvec pos = find( s_temp.col(p) == uclust[ k ] ); // pos is the set of subjects in cluster k
        
        int idx = pos[ 0 ]; // get one of the subjects in cluster k
        
        aterm += 1.5; // for each of the non-trivial cluster in smooth function p, add one to size
        
        for( int nlm = 0; nlm < 3; ++nlm ){
          bterm += 0.5 * pow( beta_temp( idx, 3*p + nlm ), 2) / lambda2_temp( idx, 3*p + nlm ); // add 0.5 * three beta^2/lambda2 terms
        } // end for nlm
        
      } // end for uclust
      
    } // end for k
    
    // ltemp = rgamma( 1, aterm, 1/bterm )[ 0 ];
    // tau2_return[ p ] = 1/ltemp;
    tau2_return[ p ] = rinvgamma_cpp(aterm, bterm);
    
  } // end for p
  
  return tau2_return;
  
}

arma::vec update_nu_tau( arma::vec tau2_temp ){
  
  int P = tau2_temp.size(); 
  arma::vec nu_tau_return( P, fill::ones );
  // double ltemp = 0;
  double bterm = 0;
  
  for( int p = 0; p < P; ++p ){
    bterm = 1 + 1 / tau2_temp( p );
    // ltemp = rgamma( 1, 1, 1/bterm )[0];
    // nu_tau_return( p ) = 1 / ltemp;
    nu_tau_return( p ) = rinvgamma_cpp(1, bterm);
  }
  
  return nu_tau_return;
  
}

arma::vec update_alpha_gibbs(
    arma::mat s_temp,
    arma::vec alpha,
    double a_alpha,
    double b_alpha){
  
  int N = s_temp.n_rows;
  int P = s_temp.n_cols;
  
  for(int p = 0; p < P; p++){
    // Get the number and the indexes of all existing clusters
    arma::vec sp = s_temp.col( p );              // the cluster indicator s_p
    arma::uvec uclust_pos = find_unique( sp );   // positions for the unique clusters, one for each cluster
    int qpn = uclust_pos.size();                 // number of existing unique clusters
    
    double ap = alpha( p );
    double tp = Rcpp::rbeta( 1, ap+1, N )[0];
    
    double p_shape = a_alpha + qpn;
    double p_rate = b_alpha - log(tp);
    double p_scale = 1 / p_rate;
    
    double A = Rcpp::rgamma( 1, p_shape, p_scale )[0];
    double B = Rcpp::rgamma( 1, p_shape - 1, p_scale )[0];
    
    double M = ( p_shape - 1 ) / p_rate / N;
    double pi_tp = M / ( 1 + M );
    
    double p_alpha = Rcpp::rbinom( 1, 1, pi_tp )[0];
    double ap_new = p_alpha*A + (1-p_alpha)*B;
    
    alpha(p) = ap_new;
    
  } // end for p
  
  return alpha;
}

arma::vec update_eta_MH(
    arma::mat dist, // the original distance matrix, where dist(i,i) = 0
    arma::mat s_temp, // current cluster membership indicator
    arma::vec alpha,  // Vector of current alphas. Indexed by p
    arma::vec delta,
    arma::vec eta,
    arma::uvec sigma,
    double a_eta,
    double b_eta, // rate parameter for gamma distribution
    int ftype // distance function type: 1-window, 2-exp, 3-logistic
){
  
  int N = s_temp.n_rows;
  int P = s_temp.n_cols;
  arma::mat lambda_dist_p_prop(N, N, fill::zeros);
  arma::mat lambda_dist_p(N, N, fill::zeros);
  double eta_prop = 0;
  
  double r = 0;
  double a = 0;
  double propose_lik = 0;
  double origin_lik = 0;
  
  for(int p = 0; p < P; ++p){
    
    eta_prop = Rcpp::rgamma( 1, a_eta, 1 / b_eta )[ 0 ];
    
    if(ftype == 1){
      lambda_dist_p_prop = help::make_lambda_dist_p_window( dist, eta_prop );
      lambda_dist_p = help::make_lambda_dist_p_window( dist, eta[p] );
    }else if(ftype == 2){
      lambda_dist_p_prop = help::make_lambda_dist_p_exp( dist, eta_prop );
      lambda_dist_p = help::make_lambda_dist_p_exp( dist, eta[p] );
    }else{
      lambda_dist_p_prop = help::make_lambda_dist_p_logis( dist, eta_prop );
      lambda_dist_p = help::make_lambda_dist_p_logis( dist, eta[p] );
    }
    
    propose_lik = 0;
    origin_lik = 0;
    
    for(int i = 0; i < N; ++i){
      propose_lik += help::log_prob_clust_single( i, p, s_temp, alpha, delta, sigma, lambda_dist_p_prop ); 
      origin_lik += help::log_prob_clust_single( i, p, s_temp, alpha, delta, sigma, lambda_dist_p );
    } // end for i
    
    r = propose_lik - origin_lik;
    a  = log( runif( 1 )[ 0 ] );
    if( a < r ){ // accept proposal if a < r
      eta[p] = eta_prop;
    }
    
  } // end for p
  
  return eta;
  
}


arma::uvec update_sigma_MH(
    int n_shuffle,
    arma::mat s_temp, // current cluster membership indicator
    arma::vec alpha,  // Vector of current alphas. Indexed by p
    arma::vec delta,
    arma::vec eta,
    arma::uvec sigma,
    arma::mat lambda_dist_p
){
  
  int N = s_temp.n_rows;
  int P = s_temp.n_cols;
  
  arma::uvec sigma_prop = sigma;
  arma::uvec shuf = randperm( N, n_shuffle); // randomly pick n_sh positions
  arma::uvec sorted_shuf = sort( shuf ); // sorted order of these positions
  sigma_prop( sorted_shuf ) = sigma_prop( shuf );
  
  double r = 0;
  double a = 0;
  double propose_lik = 0;
  double origin_lik = 0;
  
  for(int p = 0; p < P; ++p){
    for(int i = 0; i < N; ++i){
      propose_lik += help::log_prob_clust_single( i, p, s_temp, alpha, delta, sigma_prop, lambda_dist_p );
      origin_lik += help::log_prob_clust_single( i, p, s_temp, alpha, delta, sigma, lambda_dist_p );
    } // end for i
  } // end for p
  
  r = propose_lik - origin_lik;
  a  = log( runif( 1 )[ 0 ] );
  if( a < r ){ // accept proposal if a < r
    sigma = sigma_prop;
  }
  
  return sigma;
  
}



/// PY main functions ///


// Update the cluster membership according to EPA distribution
List update_cluster(
    int i, // chosen subject (will be greater than 0 here)
    int p, // current smooth function
    arma::vec between_covas,
    arma::vec subject,
    arma::vec subject_dems,
    arma::vec Y,
    arma::mat Xbar,
    arma::mat Ustar,
    arma::vec Ustar_dems,
    arma::mat beta_temp, 
    arma::mat xi_temp,
    arma::mat s_temp,             // current cluster membership indicator
    arma::mat gamma_temp,
    arma::vec w,
    arma::vec alpha,              // Vector of current alphas. Indexed by p
    arma::vec delta,
    arma::vec mu,
    arma::mat lambda2_temp,
    arma::vec tau2_temp
){
  
  int N = beta_temp.n_rows; // total number of subjects
  int rp = Ustar_dems[ p + 1 ] - Ustar_dems[ p ]; // Dimension of xi_temp_ip
  arma::vec mu_beta( 3, fill::zeros ); // Make a 3-d vector of 0s
  arma::vec mu_xi = mu( span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) ); // get mu_rp
  arma::mat one_beta = eye( 3, 3 ); // Make a 3-d identity matrix
  arma::mat one_xi = eye( rp, rp ); // Make a rp-d identity matrix
  arma::rowvec ltmp( 3, fill::zeros ); // Make a temporary space to put the lambdas
  
  // Get the number and the indexes of all existing clusters
  arma::vec sp = s_temp.col( p );              // the cluster indicator s_p
  arma::vec spmi = get_spmi(sp, i);            // cluster indicators s_p without i 
  arma::uvec uclust_pos = find_unique( spmi );   // positions for the unique clusters, one for each cluster
  int q_i = uclust_pos.size();                 // number of existing unique clusters
  
  // Make a copy of the parameters (beta_temp/xi_temp should be unchanged until the cluster to move to is decided)
  arma::mat beta_temp2 = beta_temp;
  arma::mat xi_temp2 = xi_temp;
  arma::mat s_temp2 = s_temp;
  arma::mat gamma_temp2 = gamma_temp;
  arma::mat lambda2_temp2 = lambda2_temp;
  
  // make space to save the change for moving to each cluster
  arma::mat beta_candidate( q_i + 2, 3, fill::zeros );
  arma::mat xi_candidate( q_i + 2, rp, fill::zeros );
  arma::vec s_candidate( q_i + 2, fill::zeros );
  arma::vec gamma_candidate( q_i + 2, fill::zeros );
  arma::mat lambda2_candidate( q_i + 2, 3, fill::zeros );
  arma::vec v_pi( q_i + 2, fill::zeros ); // Make space to store each probability of moving i to a updated cluster
  
  // Rcout << "--------------------------------------------"  << endl;
  // Rcout << "subject: " << i << " at smooth function:" << p << endl;
  // Rcout << "current clustering: " << sp << endl;
  // Rcout << "current selection: " << gamma_temp.col( p ) << endl;
  
  // Fill the candidates with parameters from existing clusters and Propose 1 non-zero new cluster from N(0,5) and 1 zero new cluster
  for(int k = 0; k < q_i + 2; ++k){ // propose a move to cluster k
    
    double v_pik = 0;
    
    if(k == q_i + 1){ // moving to a new beta!=0 cluster
      
      s_temp2( i, p ) = max( sp ) + 2; // Assume the indicator of the new cluster to be greater than all indicators in s_p 
      beta_temp2( i, span( 3*p, 3*p+2 ) ) = help::mvrnormArma( 1, mu_beta, 5*one_beta );
      xi_temp2( i, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) ) = help::mvrnormArma( 1, mu_xi, one_xi );
      gamma_temp2( i, p ) = 1;
      ltmp = Rcpp::rgamma(3,5,1);
      lambda2_temp2( i, span( 3*p, 3*p+2 ) ) = ltmp; // sample something greater than 1
      
    }else if( k == q_i ){ // moving to a new beta=0 cluster
      
      s_temp2( i, p ) = max( sp ) + 1; // Assume the indicator of the new cluster to be greater than all indicators in s_p 
      beta_temp2( i, span( 3*p, 3*p+2 ) ).zeros();
      xi_temp2( i, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) ) = help::mvrnormArma( 1, mu_xi, one_xi );
      gamma_temp2( i, p ) = 0;
      ltmp = Rcpp::rgamma(3,1,0.2);
      lambda2_temp2( i, span( 3*p, 3*p+2 ) ) = ltmp; // sample something close to 0
      
    }else{ // if moving to one of the previous non-trivial clusters
      
      int j = uclust_pos( k ); // find a subject j in cluster k and mark down its parameters
      
      if(j >= i){
        j++;
      } // add one back to the index so that j is the correct index in the original sp
      
      beta_temp2( i, span( 3*p, 3*p+2 ) ) = beta_temp( j, span( 3*p, 3*p+2 ) );
      xi_temp2( i, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) ) = xi_temp( j, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) );
      s_temp2( i, p ) = s_temp( j, p ); 
      gamma_temp2( i, p ) = gamma_temp( j, p );
      lambda2_temp2( i, span( 3*p, 3*p+2 ) ) = lambda2_temp( j, span( 3*p, 3*p+2 ) );
      
    } // end of if-elseif-else (k = 0,1,...,q,q+1) loop
    
    // Rcout << "move to cluster: " << s_temp2( i, p ) << " with beta: " << beta_temp2( i, span( 3*p, 3*p+2 ) ) << endl;
    
    // Update clustering using Neal algorithm 8
    uvec tmp1 = find( spmi == s_temp2( i, p ) );
    int n_kmi = tmp1.size(); // n_k minus i: the cardinality of cluster k after taking i off of s_p
    if( k < q_i ){
      v_pik = log( ( n_kmi - delta[p] ) / ( alpha[p] + N - 1 ) ); // moving to an existing cluster
    }else{
      v_pik = log( 0.5 * ( alpha[p] + delta[p] * q_i ) / ( alpha[p] + N - 1 ) ); // move to one of the 2 new clusters
    }
    
    // Rcout << "log partition probability: " << v_pik << endl;
    
    // add likelihood
    v_pik += help::log_like_logit_single( i, subject, Y, Xbar, Ustar, Ustar_dems, beta_temp2, xi_temp2 );
    v_pi[ k ] = v_pik;
    
    // Rcout << "log loglik: " << help::log_like_logit_single( i, subject, Y, Xbar, Ustar, Ustar_dems, beta_temp2, xi_temp2 ) << endl;
    
    // save parameters candidates
    beta_candidate.row( k ) = beta_temp2( i, span( 3*p, 3*p+2 ) );
    xi_candidate.row( k ) = xi_temp2( i, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) );
    s_candidate[ k ] = s_temp2( i, p );
    gamma_candidate[ k ] = gamma_temp2( i, p );
    lambda2_candidate.row( k ) = lambda2_temp2( i, span( 3*p, 3*p+2 ) );
    
  } // End for k = 0 to q+1
  
  // Sample the new cluster from all the full conditionals
  v_pi = exp(v_pi - max(v_pi));
  
  bool zero_cluster = any(between_covas == p);
  if (!zero_cluster){
    v_pi[q_i] = 0;
  }
  
  int new_clust = help::sample_prob_cpp( s_candidate, v_pi );
  
  int index = 0;
  uvec indextemp = find( s_candidate == new_clust ); // which s_candidate is selected?
  index = indextemp[0];
  
  beta_temp( i, span( 3*p, 3*p+2 ) ) = beta_candidate.row( index );
  xi_temp( i, span( Ustar_dems[ p ], Ustar_dems[ p + 1 ] - 1 ) ) = xi_candidate.row( index );
  s_temp( i, p ) = new_clust;
  gamma_temp( i, p ) = gamma_candidate( index );
  lambda2_temp( i, span( 3*p, 3*p+2 ) ) = lambda2_candidate.row( index );
  
  // Rcout << "moving to cluster: " << new_clust << endl;
  // Rcout << "new parameter: " << beta_temp( i, span( 3*p, 3*p+2 ) ) << endl;
  
  List vpi_list( 5 );
  vpi_list[ 0 ] = beta_temp;
  vpi_list[ 1 ] = xi_temp;
  vpi_list[ 2 ] = s_temp;
  vpi_list[ 3 ] = gamma_temp;
  vpi_list[ 4 ] = lambda2_temp;
  
  return vpi_list;
  
}

// An MH step to update alpha_p, proposal from gamma distribution 
double update_ap_MH(
    arma::vec sp, // current cluster membership indicator
    double ap,  // Vector of current alphas. Indexed by p
    double dp,
    double a_alpha,
    double b_alpha
){
  
  double a  = log( runif( 1 )[ 0 ] );
  double ap2 = R::rgamma(a_alpha,b_alpha);
  
  double orig_ll = help::log_epf( sp, ap, dp);
  double prop_ll = help::log_epf( sp, ap2, dp);
  // Rcout << "propose loglike: " << prop_ll << endl;
  // Rcout << "origin loglike: " << orig_ll << endl;
  double r = prop_ll - orig_ll;
  
  if( a < r ){ // accept proposal if a < r
    ap = ap2;
  }
  
  return ap;  
}

// An MH within step to update delta_p, proposal from beta distribution  
double update_dp_MH(
    arma::vec sp, // current cluster membership indicator
    double ap,  // Vector of current alphas. Indexed by p
    double dp,
    double a_delta,
    double b_delta
){
  
  double a = log( runif( 1 )[ 0 ] );
  double dp2 = R::rbeta( a_delta, b_delta);
  
  double orig_ll = help::log_epf( sp, ap, dp);
  double prop_ll = help::log_epf( sp, ap, dp2);
  // Rcout << "propose loglike: " << prop_ll << endl;
  // Rcout << "origin loglike: " << orig_ll << endl;
  double r = prop_ll - orig_ll;
  
  if( a < r ){ // accept proposal if a < r
    dp = dp2;
  }
  
  return dp;
}

/// PY main functions ///

}  // For namespace 'help'

// Function :: MCMC algorithm
// [[Rcpp::export]]
List bvsEPAcpp(
    int iterations, // number of iterations
    int n_clust_betn, // number of subjects to perform a cluster between move 
    int n_shuffle, // number of subjects to shuffle in the permutation order
    arma::vec between_covas, // columns index for covariates allowed to go through between steps
    arma::vec cluster_covas, // column indices for covariates to go through updating cluster steps
    arma::vec subject,        
    arma::vec subject_dems,   
    arma::vec Y, 
    arma::mat Xbar, 
    arma::mat Ustar,
    arma::vec Ustar_dems,  
    arma::mat dist,
    arma::vec subj_clust_free, // A vector of subjects that are open to between steps
    arma::mat beta_init, 
    arma::mat xi_init, 
    arma::mat s_init,
    arma::mat gamma_init,
    arma::vec w, 
    arma::vec alpha,
    arma::vec delta,
    arma::vec gamma_delta,
    arma::vec eta,
    arma::uvec sigma,
    arma::vec mu,
    arma::mat lambda2_init,
    arma::vec tau2_init,
    arma::mat nu_lambda_init,
    arma::vec nu_tau_init,
    double a_alpha,
    double b_alpha,
    double a_delta,
    double b_delta,
    double a_gamma,
    double b_gamma,
    double a_eta,
    double b_eta,
    int ftype // distance function type: 1-window, 2-exp, 3-logistic
){
  
  int P = tau2_init.size(); // Get the total number of smooth function
  int N = dist.n_cols;
  int nclust = 0;
  int sum_rp = mu.size();
  
  // Initialization
  arma::cube beta_cube(N, 3*P, iterations, fill::zeros);        // beta_temp - Matrix of coefficients for fixed effects. x axis indexed by 3P (beta_n1, beta_l1, beta_m1,...,beta_nP, beta_lP, beta_mP), y axis indexed by i
  arma::cube xi_cube(N, sum_rp, iterations, fill::zeros);       // xi_temp - Matrix of parameter expansion for beta. x axis indexed by sum r_p over p, y axis indexed by i
  arma::cube s_cube(N, P, iterations, fill::zeros);             // s_temp - Matrix of cluster membership indicators. x axis indexed by p, y axis indexed by i
  arma::cube gamma_cube(N, P, iterations, fill::zeros);         // gamma_temp - Matrix of inclusion indicators for fixed effects. x axis indexed by p, y axis indexed by i
  arma::cube lambda2_cube(N, 3*P, iterations, fill::zeros);     // lambda_temp - the local shrinkage parameter. x axis indexed by 3P, y axis indexed by i
  arma::mat tau2(iterations, P, fill::zeros);                   // tau_temp - the global shrinkage parameter, indexed by p
  arma::cube nu_lambda_cube(N, 3*P, iterations, fill::zeros);
  arma::mat nu_tau_mat(iterations, P, fill::zeros); 
  
  arma::mat beta_temp = beta_init;
  arma::mat xi_temp = xi_init;
  arma::mat s_temp = s_init;
  arma::mat gamma_temp = gamma_init;
  arma::mat lambda2_temp = lambda2_init;
  arma::vec tau2_temp = tau2_init;     
  arma::mat nu_lambda = nu_lambda_init;
  arma::vec nu_tau = nu_tau_init;
  
  arma::mat lambda_dist_p = dist;
  arma::mat alpha_mat(iterations, P, fill::zeros);
  arma::mat delta_mat(iterations, P, fill::zeros);
  arma::mat eta_mat(iterations, P, fill::zeros);
  arma::umat sigma_mat(iterations, N, fill::zeros);
  
  arma::mat dist_mat_window = help::make_lambda_dist_p_window( dist, eta[0] ); 
  // make a fixed distance matrix for the window decay (since we don't need to update eta)
  
  // Looping over the number of iterations specified by user
  
  for( int iter = 0; iter < iterations; ++iter ){
    
    for( int p = 0; p < P; ++p ){
      
      // FUNCTION :: Update w, for each ij
      w = help::update_w( Xbar, Ustar, Ustar_dems, beta_temp, xi_temp, subject );
      
      // FUNCTION :: between step
      bool perform_between = any(between_covas == p);
      
      if ( perform_between ){
        for ( int t = 0; t < n_clust_betn; ++t ){
          
        int subj_clust_between = help::sample_cpp( subj_clust_free ); // Find a subject to perform a between step
        List return_clust_between( 3 ); // assign space to save updated cluster membership indicator (s) and parameters (beta, xi)
        return_clust_between = help::between_step( subj_clust_between, p, subject, Y, Xbar, Ustar, Ustar_dems, beta_temp, xi_temp, s_temp, gamma_temp, mu, lambda2_temp, tau2_temp, a_gamma, b_gamma );
        beta_temp = as<arma::mat>( return_clust_between[ 0 ] );
        s_temp = as<arma::mat>( return_clust_between[ 1 ] );
        gamma_temp = as<arma::mat>( return_clust_between[ 2 ] );
          
        } // end for n_clust_betn
      }
      
      // Function :: Make distance matrix (given eta_p)
      if(ftype == 1){
        lambda_dist_p = dist_mat_window;
      }else if(ftype == 2){
        lambda_dist_p = help::make_lambda_dist_p_exp( dist, eta[p] );
      }else{
        lambda_dist_p = help::make_lambda_dist_p_logis( dist, eta[p] );
      }
      
      // FUNCTION :: update cluster membership for all subjects
      bool perform_cluster = any(cluster_covas == p);
      if (perform_cluster){
        for ( int i = 0; i < N; ++i ){
          
          List return_update_spi( 5 );
          return_update_spi = help::update_spi_EPA(sigma(i), p, between_covas, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, beta_temp, xi_temp, s_temp, gamma_temp, w, alpha, delta, sigma, mu, lambda2_temp, tau2_temp, lambda_dist_p);
          beta_temp = as<arma::mat>( return_update_spi[ 0 ] );
          xi_temp = as<arma::mat>( return_update_spi[ 1 ] );
          s_temp = as<arma::mat>( return_update_spi[ 2 ] );
          gamma_temp = as<arma::mat>( return_update_spi[ 3 ] );
          lambda2_temp = as<arma::mat>( return_update_spi[ 4 ] );
        }
      }
      
      // Put zero-coefficient subjects in one cluster and tight S_PK to prevent crashing
      s_temp.col( p ) = help::relabel_zero( s_temp.col(p), gamma_temp.col(p) );
      s_temp.col( p ) = help::tighten( s_temp.col(p) );
      
      arma::vec uclust_p = unique( s_temp.col(p) );
      nclust = uclust_p.size(); // size of all clusters: 1+q_p
      
      // FUNCTION :: within step
      for(int k = 0; k < nclust; ++k){
        
        int uk = uclust_p[k];
        int temp_g = help::getgamma( s_temp, gamma_temp, p, uk);
        
        if(temp_g == 1){
          beta_temp = help::within_beta_pk(p, uk, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, beta_temp, xi_temp, s_temp, w, lambda2_temp, tau2_temp);
          xi_temp = help::within_xi_pk(p, uk, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, beta_temp, xi_temp, s_temp, w, mu);
        }else if(temp_g == 0){
          xi_temp = help::within_xi_p0(p, uk, Ustar_dems, xi_temp, s_temp, mu);
        }
        
      } // end for k
    } // end for P
    
    // Rescale
    List return_rescaled = help::rescaled( beta_temp, xi_temp, Ustar_dems);
    beta_temp = as<arma::mat>( return_rescaled[ 0 ] );
    xi_temp = as<arma::mat>( return_rescaled[ 1 ] );
    
    // FUNCTION :: update_mu, for all r
    for( int p = 0; p < P; ++p ){
      mu = help::update_mu_rp( p, Ustar_dems, xi_temp, s_temp, mu );
    }
    
    // FUNCTION :: update tau, lambda, nu_lambda, nu_tau
    lambda2_temp = help::update_lambda( beta_temp, s_temp, gamma_temp, tau2_temp, nu_lambda );
    nu_lambda = help::update_nu_lambda( s_temp, lambda2_temp );
    tau2_temp = help::update_tau( beta_temp, s_temp, gamma_temp, lambda2_temp, nu_tau );
    nu_tau = help::update_nu_tau( tau2_temp );
    
    alpha = help::update_alpha_gibbs(s_temp, alpha, a_alpha, b_alpha);
    if(ftype != 1){
      eta = help::update_eta_MH(dist, s_temp, alpha, delta, eta, sigma, a_eta, b_eta, ftype);
    }
    sigma = help::update_sigma_MH(n_shuffle, s_temp, alpha, delta, eta, sigma, lambda_dist_p);
    
    // Rcout << sigma << endl;
    
    beta_cube.slice( iter ) = beta_temp;           // beta_temp - Matrix of coefficients for fixed effects. x axis indexed by 3P (beta_n1, beta_l1, beta_m1,...,beta_nP, beta_lP, beta_mP), y axis indexed by i
    xi_cube.slice( iter ) = xi_temp;               // xi_temp - Matrix of parameter expansion for beta. x axis indexed by sum r_p over p, y axis indexed by i
    s_cube.slice( iter ) = s_temp;                 // s_temp - Matrix of cluster membership indicators. x axis indexed by p, y axis indexed by i
    gamma_cube.slice( iter ) = gamma_temp;         // gamma_temp - Matrix of inclusion indicators for fixed effects. x axis indexed by p, y axis indexed by i
    lambda2_cube.slice( iter ) = lambda2_temp;     // lambda_temp - the local shrinkage parameter. x axis indexed by 3P, y axis indexed by i
    tau2.row( iter ) = tau2_temp.t();              // tau_temp - the global shrinkage parameter, indexed by p
    nu_lambda_cube.slice( iter ) = nu_lambda;
    nu_tau_mat.row( iter ) = nu_tau.t();
    
    alpha_mat.row( iter ) = alpha.t();
    delta_mat.row( iter ) = delta.t();
    eta_mat.row( iter ) = eta.t();
    sigma_mat.row( iter ) = sigma.t();
    
    // Print out progress
    double printer = iter % 50;
    if( printer == 0 ){
      Rcpp::Rcout << "Iteration = " << iter << std::endl;
    }
    
  } // end for iterations
  
  // Return output
  List output( 12 );
  output[ 0 ] = beta_cube;
  output[ 1 ] = xi_cube;
  output[ 2 ] = s_cube;
  output[ 3 ] = gamma_cube;
  output[ 4 ] = lambda2_cube;
  output[ 5 ] = tau2;
  output[ 6 ] = nu_lambda_cube;
  output[ 7 ] = nu_tau_mat;
  output[ 8 ] = alpha_mat;
  output[ 9 ] = delta_mat;
  output[ 10 ] = eta_mat;
  output[ 11 ] = sigma_mat;
  
  return output;
  
}

// [[Rcpp::export]]
List bvsPYcpp(
    int iterations, // number of iterations
    int n_clust_betn, // number of subjects to perform a cluster between move 
    int n_shuffle, // number of subjects to shuffle in the permutation order
    arma::vec between_covas, // columns index for covariates allowed to go through between steps
    arma::vec cluster_covas,
    arma::vec subject,        
    arma::vec subject_dems,   
    arma::vec Y, 
    arma::mat Xbar, 
    arma::mat Ustar,
    arma::vec Ustar_dems,  
    arma::vec subj_clust_free, // A vector of subjects that are open to between steps
    arma::mat beta_init, 
    arma::mat xi_init, 
    arma::mat s_init,
    arma::mat gamma_init,
    arma::vec w, 
    arma::vec alpha,
    arma::vec delta,
    arma::vec mu,
    arma::mat lambda2_init,
    arma::vec tau2_init,
    arma::mat nu_lambda_init,
    arma::vec nu_tau_init,
    double a_gamma, // indicator hyperparameter
    double b_gamma,
    double a_alpha, // concentration hyperparameter
    double b_alpha,
    double a_delta, // discount hyperparameter
    double b_delta,
    bool DP
){
  
  int P = s_init.n_cols; // Get the total number of smooth function
  int N = s_init.n_rows;
  int nclust = 0;
  int sum_rp = mu.size();
  
  // Initialization
  arma::cube beta_cube(N, 3*P, iterations, fill::zeros);        // beta_temp - Matrix of coefficients for fixed effects. x axis indexed by 3P (beta_n1, beta_l1, beta_m1,...,beta_nP, beta_lP, beta_mP), y axis indexed by i
  arma::cube xi_cube(N, sum_rp, iterations, fill::zeros);       // xi_temp - Matrix of parameter expansion for beta. x axis indexed by sum r_p over p, y axis indexed by i
  arma::cube s_cube(N, P, iterations, fill::zeros);             // s_temp - Matrix of cluster membership indicators. x axis indexed by p, y axis indexed by i
  arma::cube gamma_cube(N, P, iterations, fill::zeros);         // gamma_temp - Matrix of inclusion indicators for fixed effects. x axis indexed by p, y axis indexed by i
  arma::cube lambda2_cube(N, 3*P, iterations, fill::zeros);     // lambda_temp - the local shrinkage parameter. x axis indexed by 3P, y axis indexed by i
  arma::mat tau2(iterations, P, fill::zeros);                   // tau_temp - the global shrinkage parameter, indexed by p
  arma::cube nu_lambda_cube(N, 3*P, iterations, fill::zeros);
  arma::mat nu_tau_mat(iterations, P, fill::zeros); 
  
  arma::mat beta_temp = beta_init;
  arma::mat xi_temp = xi_init;
  arma::mat s_temp = s_init;
  arma::mat gamma_temp = gamma_init;
  arma::mat lambda2_temp = lambda2_init;
  arma::vec tau2_temp = tau2_init;     
  arma::mat nu_lambda = nu_lambda_init;
  arma::vec nu_tau = nu_tau_init;
  
  arma::mat alpha_mat(iterations, P, fill::zeros);
  arma::mat delta_mat(iterations, P, fill::zeros);
  
  // Looping over the number of iterations specified by user
  
  for( int iter = 0; iter < iterations; ++iter ){
    
    // Rcpp::Rcout << "Iteration = " << iter << std::endl;
    
    for( int p = 0; p < P; ++p ){
      
      // FUNCTION :: Update w, for each ij
      w = help::update_w( Xbar, Ustar, Ustar_dems, beta_temp, xi_temp, subject );
      
      // FUNCTION :: between step
      // Rcpp::Rcout << "between step, p=" << p << std::endl;
      
      bool perform_between = any(between_covas == p);
      
      if ( perform_between ){
        for ( int t = 0; t < n_clust_betn; ++t ){
          
          int subj_clust_between = help::sample_cpp( subj_clust_free ); // Find a subject to perform a between step
          
          List return_clust_between( 3 ); // assign space to save updated cluster membership indicator (s) and parameters (beta, xi)
          
          return_clust_between = help::between_step( subj_clust_between, p, subject, Y, Xbar, Ustar, Ustar_dems, beta_temp, xi_temp, s_temp, gamma_temp, mu, lambda2_temp, tau2_temp, a_gamma, b_gamma );
          beta_temp = as<arma::mat>( return_clust_between[ 0 ] );
          s_temp = as<arma::mat>( return_clust_between[ 1 ] );
          gamma_temp = as<arma::mat>( return_clust_between[ 2 ] );
          
        } // end fot n_clust_betn
      }
      
      
      // FUNCTION :: update cluster membership for all subjects
      bool perform_cluster = any(cluster_covas == p);
      if (perform_cluster){
        for ( int i = 0; i < N; ++i ){
          List return_update_spi( 5 );
          
          return_update_spi = help::update_cluster(i, p, between_covas, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, beta_temp, xi_temp, s_temp, gamma_temp, w, alpha, delta, mu, lambda2_temp, tau2_temp);
          
          beta_temp = as<arma::mat>( return_update_spi[ 0 ] );
          xi_temp = as<arma::mat>( return_update_spi[ 1 ] );
          s_temp = as<arma::mat>( return_update_spi[ 2 ] );
          gamma_temp = as<arma::mat>( return_update_spi[ 3 ] );
          lambda2_temp = as<arma::mat>( return_update_spi[ 4 ] );
        }
      }
      
      // Put zero-coefficient subjects in one cluster and tight S_PK to prevent crashing
      s_temp.col( p ) = help::relabel_zero( s_temp.col(p), gamma_temp.col(p) );
      s_temp.col( p ) = help::tighten( s_temp.col(p) );
      
      arma::vec uclust_p = unique( s_temp.col(p) );
      nclust = uclust_p.size(); // size of all clusters: 1+q_p
      
      // Rcpp::Rcout << "p = " << p << std::endl;
      // Rcpp::Rcout << "before within step, beta = " << beta_temp << std::endl;
      
      // FUNCTION :: within step
      for(int k = 0; k < nclust; ++k){
        
        int uk = uclust_p[k];
        int temp_g = help::getgamma( s_temp, gamma_temp, p, uk);
        
        if(temp_g == 1){
          beta_temp = help::within_beta_pk(p, uk, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, beta_temp, xi_temp, s_temp, w, lambda2_temp, tau2_temp);
          xi_temp = help::within_xi_pk(p, uk, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, beta_temp, xi_temp, s_temp, w, mu);
        }
        
        if(temp_g == 0){
          xi_temp = help::within_xi_p0(p, uk, Ustar_dems, xi_temp, s_temp, mu);
        }
        
      } // end for k
      
      // Rcpp::Rcout << "p = " << p << std::endl;
      // Rcpp::Rcout << "after within step, beta = " << beta_temp << std::endl;
      
    } // end for P
    
    // Rcpp::Rcout << "before rescale, beta = " << beta_temp << std::endl;
    // Rcpp::Rcout << "before rescale, xi = " << xi_temp << std::endl;
    // Rescale
    List return_rescaled = help::rescaled( beta_temp, xi_temp, Ustar_dems);
    beta_temp = as<arma::mat>( return_rescaled[ 0 ] );
    xi_temp = as<arma::mat>( return_rescaled[ 1 ] );
    
    // Rcpp::Rcout << "after rescale, beta = " << beta_temp << std::endl;
    // Rcpp::Rcout << "after rescale, xi = " << xi_temp << std::endl;
    
    // FUNCTION :: update_mu, for all r
    for( int p = 0; p < P; ++p ){
      mu = help::update_mu_rp( p, Ustar_dems, xi_temp, s_temp, mu );
    }
    
    // Rcpp::Rcout << "before update, lambda2 = " << lambda2_temp << std::endl;
    // Rcpp::Rcout << "before update, tau2 = " << tau2_temp << std::endl;
    
    // FUNCTION :: update tau, lambda, nu_lambda, nu_tau
    lambda2_temp = help::update_lambda( beta_temp, s_temp, gamma_temp, tau2_temp, nu_lambda );
    nu_lambda = help::update_nu_lambda( s_temp, lambda2_temp );
    tau2_temp = help::update_tau( beta_temp, s_temp, gamma_temp, lambda2_temp, nu_tau );
    nu_tau = help::update_nu_tau( tau2_temp );
    
    if(DP){
      alpha = help::update_alpha_gibbs(s_temp, alpha, a_alpha, b_alpha);
    }else{
      for( int p = 0; p < P; ++p ){
        alpha[p] = help::update_ap_MH(s_temp.col(p), alpha[p], delta[p], a_alpha, b_alpha);
        delta[p] = help::update_dp_MH(s_temp.col(p), alpha[p], delta[p], a_delta, b_delta);
      }
    }
    
    // Rcpp::Rcout << "after update, lambda2 = " << lambda2_temp << std::endl;
    // Rcpp::Rcout << "after update, tau2 = " << tau2_temp << std::endl;
    
    beta_cube.slice( iter ) = beta_temp;           // beta_temp - Matrix of coefficients for fixed effects. x axis indexed by 3P (beta_n1, beta_l1, beta_m1,...,beta_nP, beta_lP, beta_mP), y axis indexed by i
    xi_cube.slice( iter ) = xi_temp;               // xi_temp - Matrix of parameter expansion for beta. x axis indexed by sum r_p over p, y axis indexed by i
    s_cube.slice( iter ) = s_temp;                 // s_temp - Matrix of cluster membership indicators. x axis indexed by p, y axis indexed by i
    gamma_cube.slice( iter ) = gamma_temp;         // gamma_temp - Matrix of inclusion indicators for fixed effects. x axis indexed by p, y axis indexed by i
    lambda2_cube.slice( iter ) = lambda2_temp;     // lambda_temp - the local shrinkage parameter. x axis indexed by 3P, y axis indexed by i
    tau2.row( iter ) = tau2_temp.t();              // tau_temp - the global shrinkage parameter, indexed by p
    nu_lambda_cube.slice( iter ) = nu_lambda;
    nu_tau_mat.row( iter ) = nu_tau.t();
    
    alpha_mat.row( iter ) = alpha.t();
    delta_mat.row( iter ) = delta.t();
    
    // Print out progress
    double printer = iter % 50;
    if( printer == 0 ){
      Rcpp::Rcout << "Iteration = " << iter << std::endl;
    }
    
  } // end for iterations
  
  // Return output
  List output( 10 );
  output[ 0 ] = beta_cube;
  output[ 1 ] = xi_cube;
  output[ 2 ] = s_cube;
  output[ 3 ] = gamma_cube;
  output[ 4 ] = lambda2_cube;
  output[ 5 ] = tau2;
  output[ 6 ] = nu_lambda_cube;
  output[ 7 ] = nu_tau_mat;
  output[ 8 ] = alpha_mat;
  output[ 9 ] = delta_mat;
  
  return output;
  
}

