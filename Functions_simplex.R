if(!require(lpSolveAPI)) install.packages("lpSolveAPI")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(codalm)) install.packages("codalm") 
if(!require(gtools)) install.packages("gtools") 
if(!require(CUB)) install.packages("CUB") 
if(!require(compositions)) install.packages("compositions") 

library(lpSolveAPI)
library(gtools)
library(codalm)
library(CUB)
library(compositions)
library(ggplot2)

solve_simplex_lp <- function(P_list, P_prime_list, a_weights, lambda = 0, dims_grid = NULL) {
  N <- length(P_list); n_cols_A <- length(P_list[[1]]); m_rows_A <- length(P_prime_list[[1]])
  num_vars_A <- m_rows_A * n_cols_A; num_vars_t <- N * (m_rows_A - 1)
  
  use_reg <- (lambda > 0 && !is.null(dims_grid))
  edges <- list(); num_vars_reg <- 0
  
  if (use_reg) {
    
    dims <- dims_grid
    D <- length(dims)
    
    coord_list <- lapply(dims, seq_len)
    coords <- as.matrix(expand.grid(coord_list))
    
    perm_order <- do.call(order, c(list(rowSums(coords)), as.data.frame(coords)))
    
    map_vec <- integer(nrow(coords))
    for (k in seq_along(perm_order)) {
      map_vec[perm_order[k]] <- k
    }
    
    coord_to_index <- function(coord) {
      row_idx <- 1
      mult <- 1
      for (d in seq_len(D)) {
        row_idx <- row_idx + (coord[d] - 1) * mult
        mult <- mult * dims[d]
      }
      map_vec[row_idx]
    }
    
    edges <- list()
    
    for (cell in seq_len(nrow(coords))) {
      base <- coords[cell, ]
      for (d in seq_len(D)) {
        if (base[d] < dims[d]) {
          neigh <- base
          neigh[d] <- neigh[d] + 1
          edges[[length(edges) + 1]] <- c(
            coord_to_index(base),
            coord_to_index(neigh)
          )
        }
      }
    }
    
    num_vars_reg <- length(edges) * (m_rows_A - 1)
  }
  
  total_vars <- num_vars_A + num_vars_t + num_vars_reg
  lp_model <- make.lp(0, total_vars); lp.control(lp_model, sense = "min", verbose="neutral") 
  
  obj_t <- rep(a_weights, N); obj_reg <- if(use_reg) rep(lambda, num_vars_reg) else numeric(0)
  set.objfn(lp_model, c(rep(0, num_vars_A), obj_t, obj_reg))
  set.bounds(lp_model, lower = rep(0, total_vars), columns = 1:total_vars)
  
  get_A_idx <- function(r, c) (c - 1) * m_rows_A + r
  get_t_idx <- function(obs, k) num_vars_A + (obs - 1) * (m_rows_A - 1) + k
  get_u_idx <- function(edge, k) num_vars_A + num_vars_t + (edge - 1) * (m_rows_A - 1) + k
  
  for (j in 1:n_cols_A) {
    indices <- sapply(1:m_rows_A, get_A_idx, c = j)
    add.constraint(lp_model, rep(1, m_rows_A), "=", 1, indices = indices)
  }
  
  CDF_targets <- lapply(P_prime_list, cumsum)
  for (i in 1:N) {
    P_in <- P_list[[i]]; nz_idx <- which(P_in > 0)
    for (k in 1:(m_rows_A - 1)) {
      t_idx <- get_t_idx(i, k); col_indices <- numeric(); col_values <- numeric()
      for(j in nz_idx) for(r in 1:k) {
        col_indices <- c(col_indices, get_A_idx(r, j)); col_values <- c(col_values, P_in[j])
      }
      rhs <- CDF_targets[[i]][k]
      add.constraint(lp_model, c(col_values, -1), "<=", rhs, indices = c(col_indices, t_idx))
      add.constraint(lp_model, c(col_values, 1), ">=", rhs, indices = c(col_indices, t_idx))
    }
  }
  
  if (use_reg) {
    for (e in 1:length(edges)) {
      c1 <- edges[[e]][1]; c2 <- edges[[e]][2]
      for (k in 1:(m_rows_A - 1)) {
        u_idx <- get_u_idx(e, k)
        idx_c1 <- sapply(1:k, function(r) get_A_idx(r, c1)); idx_c2 <- sapply(1:k, function(r) get_A_idx(r, c2))
        vals <- c(rep(1, length(idx_c1)), rep(-1, length(idx_c2)), -1)
        add.constraint(lp_model, vals, "<=", 0, indices = c(idx_c1, idx_c2, u_idx))
        vals2 <- c(rep(1, length(idx_c1)), rep(-1, length(idx_c2)), 1)
        add.constraint(lp_model, vals2, ">=", 0, indices = c(idx_c1, idx_c2, u_idx))
      }
    }
  }
  
  res_code <- solve(lp_model)
  if(res_code != 0) return(NULL) 
  vars <- get.variables(lp_model)
  A_opt <- matrix(vars[1:num_vars_A], nrow = m_rows_A, ncol = n_cols_A)
  return(list(A = A_opt))
}


#LEAVE-ONE-OUT CROSS-VALIDATION (LOOCV)
do_loocv <- function(Z, Y, w, dims, lambda_seq) {
  
  results <- data.frame(lambda = lambda_seq, error = NA, se = NA)
  N_obs <- length(Z)
  
  #cat(paste("Avvio LOOCV su", N_obs, "osservazioni per", length(lambda_seq), "valori di Lambda...\n"))
  #pb <- txtProgressBar(min = 0, max = length(lambda_seq), style = 3)
  
  for(i in seq_along(lambda_seq)) {
    lam <- lambda_seq[i]
    errors <- numeric(N_obs)
    
    for(k in 1:N_obs) {
      
      train_Z <- Z[-k]
      train_Y <- Y[-k]
      test_Z <- Z[[k]]
      test_Y <- Y[[k]]
      
      model <- solve_simplex_lp(train_Z, train_Y, w, lambda=lam, dims_grid=dims)
      
      if(!is.null(model)) {
        pred <- as.vector(model$A %*% test_Z)
        errors[k] <- wd(w, pred, test_Y)
      } else {
        errors[k] <- NA 
      }
    }
    
    results$error[i] <- mean(errors, na.rm=TRUE)
    # Standard Error 
    results$se[i] <- sd(errors, na.rm=TRUE) / sqrt(N_obs)
    
    #setTxtProgressBar(pb, i)
  }
  #close(pb)
  return(results)
}

#Grid
lambda_grid <- c(0, 0.002, 0.005, 0.008, 0.01, 0.012, 0.015, 0.018, 0.02, 0.025, 0.03, 0.05, 0.1)
                                                                              

####################distance functions##############################
wd <- function(weights, P, Q){
  n <- length(P)
  sum(weights * abs(cumsum(P)[1:(n-1)] - cumsum(Q)[1:(n-1)]))
}

kld<-function(P,Q){
  n<-length(P)
  dist=0
  for(i in 1:n){
    if(P[i]==0){
      P[i]=0.0001
    }
    if(Q[i]==0){
      Q[i]=0.0001
    }
    dist=dist+P[i]*(log(Q[i])-log(P[i]))
  }
  dist = -dist
  return(dist)
}

# --- Global Wasserstein distance between Matrixes ---
matrix_wasserstein_dist <- function(A, B, weights) {
  m <- nrow(A); n <- ncol(A)
  total_dist <- 0
  for(j in 1:n) {
    col_A <- A[, j]; col_B <- B[, j]
    cdf_A <- cumsum(col_A); cdf_B <- cumsum(col_B)
    diffs <- abs(cdf_A[1:(m-1)] - cdf_B[1:(m-1)])
    total_dist <- total_dist + sum(diffs * weights)
  }
  return(total_dist)
}

###########################################################
inv_logit <- function(logit) exp(logit) / (1 + exp(logit))

opiwd <- function(weights, P, Q, tol = 1e-8) {
  n <- length(P)
  
  for (k in 1:n) {
    ek <- rep(0, n)
    ek[k] <- 1
    
    dP <- wd(weights, P, ek)
    dQ <- wd(weights, Q, ek)
    
    if (dP > dQ + tol) return(-1)  # P < Q
    if (dP < dQ - tol) return(1)   # P > Q
  }
  
  return(0)  # P == Q
}

opi <- function(P, Q, tol = 1e-8) {
  n <- length(P)
  
  cP <- cumsum(P)[1:(n-1)]
  cQ <- cumsum(Q)[1:(n-1)]
  
  leq <- cP <= cQ + tol
  geq <- cP >= cQ - tol
  
  if (all(leq) && any(cP < cQ - tol)) {
    return(-1)   # P < Q
  }
  
  if (all(geq) && any(cP > cQ + tol)) {
    return(1)    # P > Q
  }
  
  if (all(abs(cP - cQ) < tol)) {
    return(2)    # equivalenti
  }
  
  return(0)     # non confrontabili
}

####################################à
tensor_product_ordered <- function(...) {
  comps <- list(...)
  lengths <- sapply(comps, length)
  
  # all combinations of indexes
  index_grid <- do.call(expand.grid, lapply(lengths, function(n) 1:n))
  
  # product 
  values <- apply(index_grid, 1, function(idx) {
    prod(mapply(`[`, comps, idx))
  })
  
  index_sum <- rowSums(index_grid)
  
  # first order: sum of indexes
  # second order: natural order
  ord <- do.call(order, c(list(index_sum), as.data.frame(index_grid)))
  
  list(
    #indices    = index_grid[ord, , drop = FALSE],
    #index_sum  = index_sum[ord],
    product    = values[ord]
  )
}

#############plot functions################
# Ternary transformation (for plot) ---
# (p1, p2, p3) -> (x, y) in the triangol
ternary_to_xy <- function(p_mat) {
  p2 <- p_mat[,2]
  p3 <- p_mat[,3]
  data.frame(
    x = 0.5 * p2 + 1.0 * p3,
    y = (sqrt(3) / 2) * p2
  )
}


get_cloud_points <- function(matrix_list, col_idx) {
  mat_col <- do.call(rbind, lapply(matrix_list, function(M) M[, col_idx]))
  df_xy <- ternary_to_xy(mat_col)
  df_xy$Column <- paste0("C", col_idx) 
  return(df_xy)
}
