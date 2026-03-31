source(file="Functions_simplex.R")
#Case I
run_simulation_caseI <- function(
    iter = 1000,
    N = 100,
    Cx = 5,
    Cy = 8,
    dgp_type = "laplace",
    lambda = 2,
    sigma = 1.2,
    alpha_dir = 80,
    equidistance = TRUE,
    gamma = 1.6,
    latent_thresholds = NULL,
    tol_equal = 1e-8
){
  
  result <- matrix(NA, nrow = iter, ncol = 3)
  
  colnames(result) <- c(
    "W1",
    "OPI",
    "R2"
  )
  
  for(it in 1:iter){
    
    data <- generate_ordinal_data(
      N = N, Cx = Cx, Cy = Cy,
      dgp_type = dgp_type,
      lambda = lambda,
      sigma = sigma,
      alpha_dir = alpha_dir,
      equidistance = equidistance,
      gamma = gamma,
      latent_thresholds = latent_thresholds
    )
    
    X <- data$X
    Y <- data$Y
    
    weights_unit <- rep(1, Cy-1)
    
    wsum_u <- sum(weights_unit)
    
    
    X_list <- lapply(1:N, function(i) X[i,])
    Y_list <- lapply(1:N, function(i) Y[i,])
    
    A_unit <- solve_simplex_lp(X_list,Y_list,weights_unit)$A
    
    # ==========================
    # W1 (SSE locale)
    # ==========================
    
    SSE_u <- 0
    pred_u<-matrix(0,nrow=N,ncol=Cy)
    for(i in 1:N){
      pred_u[i,] <- as.vector(A_unit %*% X[i,])
      }
    
    for(i in 1:N){
      SSE_u <- SSE_u + wd(weights_unit,Y[i,],pred_u[i,])
      }
    
    SSE_u <- SSE_u
    
    result[it,"W1"] <- SSE_u/(wsum_u*N)
     
    
    result[it,"R2"] <- compute_R2(Y,X,A_unit,weights_unit)$R2
    
    
    ind_u  <- 0
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        if(opiwd(weights_unit,Y[i,],Y[j,])==
           opiwd(weights_unit,pred_u[i,],pred_u[j,]))
          ind_u <- ind_u+1
      }
    }
    
    result[it,"OPI"] <- ind_u/choose(N,2)
  }
  
  
  list(
    mean_res = colMeans(result,na.rm=TRUE),
    sd_res   = apply(result,2,sd,na.rm=TRUE),
    raw_results = result
  )
}

s1 <- run_simulation_caseI(
  iter = 500,
  N = 50,
  Cx = 3,
  Cy = 3,
  dgp_type = "laplace",
  equidistance = TRUE
)

#Case III
run_simulation_caseIII <- function(
    iter = 500,
    N = 100,
    Cx1 = 5,
    Cx2 = 5,
    Cy = 8,
    dgp_type = "laplace",
    lambda = 2,
    sigma = 1.2,
    alpha_dir = 70,
    equidistance = TRUE,
    gamma = 1.5,
    w1 = 0.6,
    w2 = 0.4
){
  
  weights <- rep(1, Cy-1)
  weight_sum <- sum(weights)
  
  result <- matrix(0, nrow=iter, ncol=3)
  colnames(result) <- c(
    "W1","OPI", "R2"
  )
  
  it <- 1
  
  while(it <= iter){
    
    success <- tryCatch({
      
      data <- generate_ordinal_data_2x(
        N=N,Cx1=Cx1,Cx2=Cx2,Cy=Cy,
        dgp_type=dgp_type,
        lambda=lambda,sigma=sigma,
        alpha_dir=alpha_dir,
        equidistance=equidistance,
        gamma=gamma,w1=w1,w2=w2
      )
      
      X1 <- data$X1; X2 <- data$X2
       Y2 <- data$Y2
      
      x1list <- lapply(1:N,function(i) X1[i,])
      x2list <- lapply(1:N,function(i) X2[i,])
     
      y2list <- lapply(1:N,function(i) Y2[i,])
     
      
      Z_full <- lapply(1:N,function(i)
        tensor_product(comps=list(X1[i,], X2[i,]),list(rep(1,Cx1-1),rep(1,Cx2-1)),return_indices = FALSE)$product)
      
     
      A_full_Y2 <- solve_simplex_lp(Z_full,y2list,weights)$A
      
      # ======================
      # W1 (LASCIATO INVARIATO)
      # ======================
      
      W1_fun <- function(A,Y,x){
        mean(sapply(1:N,function(i)
          wd(weights,Y[i,],A%*%x[[i]]))) /sum(weights)
      }
      
      
      W1_full_Y2 <- W1_fun(A_full_Y2,Y2,Z_full)
      
      # ======================
      # OPI (INVARIATO)
      # ======================
      
      OPI_fun <- function(A,Y,x){
        count <- 0
        for(a in 1:(N-1)){
          for(b in (a+1):N){
            if(opiwd(weights,Y[a,],Y[b,]) ==
               opiwd(weights,A%*%x[[a]],A%*%x[[b]])){
              count <- count + 1
            }
          }
        }
        count/choose(N,2)
      }
      
      
      OPI_full_Y2 <- OPI_fun(A_full_Y2,Y2,Z_full)
      Z_mat<-do.call(rbind, Z_full)
      v2 <- compute_R2(Y2,Z_mat,A_full_Y2,weights)$R2
      
      # SALVA
      result[it,"W1"] <- W1_full_Y2
      result[it,"OPI"] <- OPI_full_Y2
        result[it,"R2"] <- v2
      
      TRUE
      
    }, error=function(e){FALSE})
    
    it <- it + 1
  }
  
  list(
    mean_res = colMeans(result,na.rm=TRUE),
    sd_res   = apply(result,2,sd,na.rm=TRUE),
    raw_results = result
  )
}

s2 <- run_simulation_caseIII(
  iter = 500,
  N = 50,
  Cx1 = 5,
  Cx2 = 5,
  Cy  = 5,
  dgp_type = "laplace",
  w1 = 0.6,
  w2 = 0.4
)
