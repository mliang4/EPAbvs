genData = function(
  N = 60, # number of subjects
  P = 15, # number of smooth functions (must be greater than 3)
  n_i = 60, # a number or a vector represent observations of each subject
  cor = 0.3, # correlation between continuous covariates
  dist_x = NULL,
  seed = 121113
){
  
  library( mvtnorm )
  library( spikeSlabGAM )
  
  # Adjust the number of observations per person if necessary
  
  if( length( n_i ) == 1 ){ n_i <- rep( n_i, N ) }
  set.seed( seed )
  
  # smooth functions
  
  f1 = function(x){
    pi*x*cos(5*pi*x) - 1.2*x
  }
  f2 = function(x){
    pi*sin(3*pi*x) + 1.4*x - 1.6
  }
  f3 = function(x){
    pi*cos(2*pi*x) + 1.6
  }
  f4 = function(x){
    -pi*cos(2*pi*x) + 1.6*x
  }
  f5 = function(x){
    pi*sin(5*x) - cos(pi*x)
  }
  f6 = function(x){
    rep(0,length(x))
  }
  
  # Make X - cor is the correlation between covariates
  
  sigma2 <- diag( P )
  for( i in 1:P ){
    for( j in 1:P ){
      if( i != j ){
        sigma2[ i , j ] = cor^abs(i - j)
      }
    }
  }
  
  # Simulate a covariate for each subject
  X <- rmvnorm( N, rep( 0, nrow( sigma2 ) ), sigma2 )
  
  # Replicate some and make the others wander
  X_ij <- numeric()
  cols <- sample( seq(1,P), floor( P/2 ) )
  
  # cols <- c(2,3,4,5)
  for( i in 1:N ){
    # Replicate by number of observations
    x_ij <- matrix( rep( X[ i, ], n_i[ i ] ), ncol = P, byrow = TRUE )
    # Give 1/2 of the columns random noise
    x_ij[, cols ] <- x_ij[, cols ] + matrix( rnorm( length( cols )*n_i[ i ] ), ncol = length( cols ), byrow = TRUE )
    # Append observations
    X_ij <- rbind( X_ij, x_ij )
  }
  
  # Add intercept term to X_ij
  X_ij <- scale( X_ij )
  X_ij[,1] <- 1
  
  # Create 4 covariates for the distance matrix
  if( is.null( dist_x ) ){
    binary_temp = matrix(data=0, nrow=N, ncol=2)
    contin_temp = matrix(data=0, nrow=N, ncol=2)
    for(i in 1:(N/3)){
      binary_temp[i,] = rbinom(2,1,0.9)
      contin_temp[i,] = rnorm(2,-3,1)
    }
    for(i in (N/3+1):(2*N/3)){
      binary_temp[i,] = rbinom(2,1,0.1)
      contin_temp[i,] = rnorm(2,0,1)
    }
    for(i in (2*N/3+1):N){
      binary_temp[i,] = rbinom(2,1,0.5)
      contin_temp[i,] = rnorm(2,3,1)
    }
    dist_x = cbind(binary_temp,contin_temp)
  }
  
  dist = as.matrix(dist(dist_x))
  
  # Make U - Assumes that the input into spline function is the same
  T <- numeric()
  for( j in 1:length( n_i ) ){
    T <- c( T, sort( runif( n_i[ j ], 0, 1 )))
  }
  U <- matrix( rep( T, ncol( X_ij) ), ncol = ncol( X_ij) )
  
  # Ustar - Matrix of spline functions. Rows indexed by ij. Columns indexed by sum r_p over p. ( Should be a list if r_p != r_p' )
  # Ustar_dems - Input vector of dimensions for each spline function. Element is equal to starting indicies for corresponding p. Last term is number of columns. Ustar_dems[ 0 ] = 0 and length = P + 1
  # Allows for only one input for U
  
  # if( ncol( as.matrix( U ) ) == 1 ){
  #   Ustar <- numeric()
  #   Ustar_dems <- c( 0 )
  #   tmp <- sm( U, rankZ = 0.995 )
  #   for( p in 1:P ){
  #     Ustar <- cbind( Ustar, tmp )
  #     Ustar_dems <- c( Ustar_dems, ( Ustar_dems[ p ] + ncol( tmp ) ) )
  #   }
  #   U <- matrix( rep( U, P ), ncol = P , nrow = length( Y ) )
  # }else{
  #   Ustar <- numeric()
  #   Ustar_dems <- c( 0 )
  #   for( p in 1:ncol( U ) ){
  #     tmp <- sm( U[ , p ], rankZ = 0.995 )
  #     Ustar <- cbind( Ustar, tmp )
  #     Ustar_dems <- c( Ustar_dems, ( Ustar_dems[ p ] + ncol( tmp ) ) )
  #   }
  # }
  
  # Matrix of barX. Rows indexed by ij. Columns indexed by 2P ( x1u, x1, x2u, x2,...xPu, xP )
  Xbar <- numeric()
  for( p in 1:ncol( U ) ){
    tmp <- X_ij[ , p]*U[ , p ]
    Xbar <- cbind( Xbar, tmp, X_ij[ , p ] )
  }
  
  # Initialize the log-odds
  psi <- rep( 0 , sum( n_i ) )
  
  # main effects, linear interactions and non-linear interactions
  
  ## non-linear interactions
  
  o2sind = NULL
  for(i in 1:N){
    o2sind = c(o2sind,rep(i,n_i[i])) # a set of index mapping obs to subject
  }
  
  # smooth function/clustering assignment: sm1: 1-20 f2, 21-60 f5;
  # sm2: 1-20 f3, 21-40 f4, 41-60 f3
  # sm3: 1-40 f1, 41-60 f6
  # sm4-smP (if exist): f6
  
  #---------------------sm1-----------------------------
  
  smOut = numeric(length(o2sind))
  sub1 = 1:(N/3)
  sub2 = (N/3+1):N
  
  obstemp = which(o2sind %in% sub1)
  t1 = U[obstemp,1]
  smOut[obstemp] = f2(t1)
  
  obstemp = which(o2sind %in% sub2)
  t2 = U[obstemp,1]
  smOut[obstemp] = f5(t2)
  
  # par(mfrow=c(2,3))
  # plot(smOut[(1-1)*ni+ (1:ni)],type="l",col="red")
  # plot(smOut[(21-1)*ni+ (1:ni)],type="l",col="red")
  
  psi = psi + smOut * Xbar[,2]
  
  #--------------------sm2------------------------------
  
  smOut = numeric(length(o2sind))
  sub1 = c(1:(N/3),(2*N/3+1):N)
  sub2 = (N/3+1):(2*N/3)
  
  obstemp = which(o2sind %in% sub1)
  t3 = U[obstemp,2]
  smOut[obstemp] = f3(t3)
  
  obstemp = which(o2sind %in% sub2)
  t4 = U[obstemp,2]
  smOut[obstemp] = f4(t4)
  
  psi = psi + smOut * Xbar[,4]
  
  # par(mfrow=c(2,3))
  # plot(smOut[(1-1)*ni+ (1:ni)],type="l",col="red")
  # plot(smOut[(16-1)*ni+ (1:ni)],type="l",col="red")
  # plot(smOut[(31-1)*ni+ (1:ni)],type="l",col="red")
  
  #--------------------sm3------------------------------
  
  smOut = numeric(length(o2sind))
  sub1 = 1:(2*N/3)
  sub2 = (2*N/3+1):N
  
  obstemp = which(o2sind %in% sub1)
  t5 = U[obstemp,3]
  smOut[obstemp] = f1(t5)
  
  obstemp = which(o2sind %in% sub2)
  t6 = U[obstemp,3]
  smOut[obstemp] = f6(t6)
  
  psi = psi + smOut * Xbar[,6]
  
  # par(mfrow=c(2,3))
  # plot(smOut[(1-1)*ni+ (1:ni)],type="l",col="red")
  # plot(smOut[(11-1)*ni+ (1:ni)],type="l",col="red")
  
  prob <- exp( psi )/( 1 + exp( psi ) )
  Y <- rbinom( length(prob) , 1 , prob )
  
  data_sim <- list( 
    "f1" = f1,
    "f2" = f2,
    "f3" = f3,
    "f4" = f4,
    "f5" = f5,
    "f6" = f6,
    "Y" = Y,
    "X" = X_ij,
    "dist" = dist,
    "n_i" = n_i,
    "T" = T
     )
  return( data_sim ) 
}