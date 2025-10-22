library ( lpSolveAPI )
library ( gtools )
library( codalm )
library( CUB )
library( compositions )

solve_simplex_lp <- function ( P_list , P_prime_list , a_weights ) {
  N <- length ( P_list )
  n_dim <- length ( P_list[[1]])
  m_dim <- length ( P_prime_list[[1]])
  if ( length ( a_weights ) != m_dim - 1) stop ( "Weight vector length must be m-1." )
  
  num_vars_A <- m_dim * n_dim
  num_vars_t <- N * ( m_dim - 1)
  lp_model <- make.lp (0 , num_vars_A + num_vars_t )
  lp.control ( lp_model , sense = "min" )
  
  obj_coeffs <- c ( rep (0 , num_vars_A ) , rep ( a_weights , N ) )
  set.objfn ( lp_model , obj_coeffs )
  
  get_A_idx <- function (r , c ) ( c - 1) * m_dim + r
  get_t_idx <- function (i , k ) num_vars_A + ( i - 1) * ( m_dim - 1) + k
  
  set.bounds ( lp_model , lower = rep (0 , num_vars_A ) , columns = 1: num_vars_A )
  for ( j in 1: n_dim ) {
    indices <- sapply (1: m_dim , get_A_idx , c = j )
    add.constraint ( lp_model , rep (1 , m_dim ) , "=" , 1 , indices = indices )
  }
  
  CDF_targets <- lapply ( P_prime_list , cumsum )
  for ( i in 1: N ) {
    P_i <- P_list[[ i ]]
    for ( k in 1:( m_dim - 1) ) {
      t_idx <- get_t_idx (i , k )
      A_coeffs <- numeric ( num_vars_A ) 
      for ( j in 1: n_dim ) {
        if ( P_i[ j ] > 0) {
          for ( s in 1: k ) {
            A_coeffs [ get_A_idx (s , j ) ] <- A_coeffs [ get_A_idx (s , j ) ] + P_i[ j ] 
          }
        }
      }
      
      nz_idx <- which( A_coeffs != 0)
      rhs <- CDF_targets[[ i ]][ k ]
      add.constraint ( lp_model , c ( A_coeffs[ nz_idx ] , -1) , "<=" , rhs , c ( nz_idx , t_idx ) )
      add.constraint ( lp_model , c ( - A_coeffs[ nz_idx ] , -1) , "<=" , -rhs , c ( nz_idx ,t_idx ) )
    }
  }
  
  solve( lp_model )
  A_opt <- matrix ( get.variables ( lp_model ) [1: num_vars_A ] , nrow = m_dim , ncol = n_dim )
  return ( list ( A = A_opt , loss = get.objective ( lp_model ) ) )
}


####################distance functions##############################
wd<-function(weights, P, Q){
  n<-length(P)
  p<-vector()
  p[1]=P[1]
  for(i in 2:n){
    p[i]<-p[i-1]+P[i]
  }
  q<-vector()
  q[1]=Q[1]
  for(i in 2:n){
    q[i]<-q[i-1]+Q[i]
  }
  dist<-sum(weights*abs(p-q))
  return(dist)
}

kld<-function(P,Q){
  n<-length(P)
  dist=0
  for(i in 1:n){
    if(P[i]==0){
      P[i]=0.0001
    }
    dist=dist+P[i]*(log(Q[i])-log(P[i]))
  }
  dist = -dist
  return(dist)
}

###########################################################
inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
