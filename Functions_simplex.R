if(!require(lpSolveAPI)) install.packages("lpSolveAPI")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(codalm)) install.packages("codalm") 
if(!require(gtools)) install.packages("gtools")  
if(!require(compositions)) install.packages("compositions") 
if(!require(dplyr)) install.packages("dplyr") 
if(!require(tidyr)) install.packages("tidyr") 
if(!require(reshape2)) install.packages("reshape2") 

library(lpSolveAPI)
library(gtools)
library(codalm)
library(compositions)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

safe_dirichlet <- function(alpha, scale = 10, eps = 1e-6){
  alpha <- scale * alpha
  alpha[alpha < eps] <- eps
  alpha[!is.finite(alpha)] <- eps
  return(rdirichlet(1, alpha))
}


solve_simplex_lp<-function(
    P_list, P_prime_list, a_weights,
    lambda1 = 0,
    lambda2 = 0
){

  N <- length(P_list)
  n_cols_A <- length(P_list[[1]])
  m_rows_A <- length(P_prime_list[[1]])

  num_vars_A <- m_rows_A * n_cols_A
  num_vars_t <- N * (m_rows_A - 1)

  use_reg1 <- (lambda1 > 0)
  use_reg2 <- (lambda2 > 0)

  # ---- EDGES primo ordine ----
  edges1 <- list()
  num_vars_reg1 <- 0

  if(use_reg1 && n_cols_A > 1){
    edges1 <- lapply(1:(n_cols_A - 1), function(j) c(j, j + 1))
    num_vars_reg1 <- length(edges1) * (m_rows_A - 1)
  }

  # ---- EDGES secondo ordine ----
  edges2 <- list()
  num_vars_reg2 <- 0

  if(use_reg2 && n_cols_A > 2){
    edges2 <- lapply(2:(n_cols_A - 1), function(j) c(j-1, j, j+1))
    num_vars_reg2 <- length(edges2) * (m_rows_A - 1)
  }

  total_vars <- num_vars_A + num_vars_t + num_vars_reg1 + num_vars_reg2

  lp_model <- make.lp(0, total_vars)
  lp.control(lp_model, sense = "min", verbose = "neutral")

  # ---- OBJECTIVE ----
  obj <- c(
    rep(0, num_vars_A),
    rep(a_weights, N),
    rep(lambda1 * a_weights, length(edges1)),
    rep(lambda2 * a_weights, length(edges2))
  )

  set.objfn(lp_model, obj)
  set.bounds(lp_model, lower = rep(0, total_vars), columns = 1:total_vars)

  # ---- INDEX HELPERS ----
  get_A_idx <- function(r, c) (c - 1) * m_rows_A + r
  get_t_idx <- function(obs, k) num_vars_A + (obs - 1) * (m_rows_A - 1) + k

  get_u1_idx <- function(e, k)
    num_vars_A + num_vars_t + (e - 1) * (m_rows_A - 1) + k

  get_u2_idx <- function(e, k)
    num_vars_A + num_vars_t + num_vars_reg1 + (e - 1) * (m_rows_A - 1) + k

  # ---- SIMPLEX ----
  for(j in 1:n_cols_A){
    idx <- sapply(1:m_rows_A, get_A_idx, c = j)
    add.constraint(lp_model, rep(1, m_rows_A), "=", 1, indices = idx)
  }

  # ---- WASSERSTEIN FIT ----
  CDF_targets <- lapply(P_prime_list, cumsum)

  for(i in 1:N){
    P_in <- P_list[[i]]
    nz_idx <- which(P_in > 0)

    for(k in 1:(m_rows_A - 1)){
      t_idx <- get_t_idx(i, k)

      col_indices <- c()
      col_values  <- c()

      for(j in nz_idx){
        for(r in 1:k){
          col_indices <- c(col_indices, get_A_idx(r, j))
          col_values  <- c(col_values, P_in[j])
        }
      }

      rhs <- CDF_targets[[i]][k]

      add.constraint(lp_model,
                     c(col_values, -1),
                     "<=",
                     rhs,
                     indices = c(col_indices, t_idx))

      add.constraint(lp_model,
                     c(col_values, 1),
                     ">=",
                     rhs,
                     indices = c(col_indices, t_idx))
    }
  }

  # ---- REG 1 (NOW WEIGHTED) ----
  if(use_reg1){
    for(e in seq_along(edges1)){
      c1 <- edges1[[e]][1]
      c2 <- edges1[[e]][2]

      for(k in 1:(m_rows_A - 1)){
        u_idx <- get_u1_idx(e, k)

        idx1 <- sapply(1:k, function(r) get_A_idx(r, c1))
        idx2 <- sapply(1:k, function(r) get_A_idx(r, c2))

        w <- a_weights[k]

        add.constraint(lp_model,
                       c(rep(w, k), rep(-w, k), -1),
                       "<=",
                       0,
                       indices = c(idx1, idx2, u_idx))

        add.constraint(lp_model,
                       c(rep(w, k), rep(-w, k), 1),
                       ">=",
                       0,
                       indices = c(idx1, idx2, u_idx))
      }
    }
  }

  # ---- REG 2 (NOW WEIGHTED) ----
  if(use_reg2){
    for(e in seq_along(edges2)){
      c1 <- edges2[[e]][1]
      c2 <- edges2[[e]][2]
      c3 <- edges2[[e]][3]

      for(k in 1:(m_rows_A - 1)){
        u_idx <- get_u2_idx(e, k)

        idx1 <- sapply(1:k, function(r) get_A_idx(r, c1))
        idx2 <- sapply(1:k, function(r) get_A_idx(r, c2))
        idx3 <- sapply(1:k, function(r) get_A_idx(r, c3))

        w <- a_weights[k]

        add.constraint(lp_model,
                       c(rep(w, k), rep(-2*w, k), rep(w, k), -1),
                       "<=",
                       0,
                       indices = c(idx1, idx2, idx3, u_idx))

        add.constraint(lp_model,
                       c(rep(w, k), rep(-2*w, k), rep(w, k), 1),
                       ">=",
                       0,
                       indices = c(idx1, idx2, idx3, u_idx))
      }
    }
  }

  res_code <- solve(lp_model)
  if(res_code != 0) return(NULL)

  vars <- get.variables(lp_model)

  A_opt <- matrix(vars[1:num_vars_A],
                  nrow = m_rows_A,
                  ncol = n_cols_A)

  return(list(A = A_opt))
}

                       
select_lambda <- function(
    P_list,
    P_prime_list,
    a_weights,
    lambda1_grid,
    lambda2_grid,
    method = c("cv", "gcv"),
    K = NULL,
    seed = 123,
    verbose = FALSE
){

  method <- match.arg(method)

  N <- length(P_list)
  m <- length(P_prime_list[[1]])

  lambda1_grid <- sort(unique(lambda1_grid))
  lambda2_grid <- sort(unique(lambda2_grid))

  
  compute_fit_single <- function(A, P, Pp){

    F_mix <- rep(0, m)

    for(j in 1:length(P)){
      F_mix <- F_mix + P[j] * cumsum(A[,j])
    }

    fit <- sum(
      a_weights * abs(
        F_mix[-m] - cumsum(Pp)[-m]
      )
    )

    if(is.na(fit) || is.infinite(fit)){
      return(Inf)
    }

    fit
  }

  compute_fit <- function(A, P_list, P_prime_list){

    fit <- 0

    for(i in 1:length(P_list)){

      F_mix <- rep(0, m)

      for(j in 1:length(P_list[[i]])){
        F_mix <- F_mix +
          P_list[[i]][j] * cumsum(A[,j])
      }

      fit <- fit +
        sum(
          a_weights[1:(m-1)] *
            abs(
              F_mix[-m] -
                cumsum(P_prime_list[[i]])[-m]
            )
        )
    }

    fit
  }

  # =====================================================
  # CROSS-VALIDATION
  # =====================================================
  if(method == "cv"){

    # vera CV = LOOCV
    if(is.null(K)){
      K <- N
    }

    if(K > N){
      stop("K must be lower then the number of observations.")
    }

    set.seed(seed)

    if(K == N){

      folds <- seq_len(N)

    } else {

      folds <- sample(rep(1:K, length.out = N))

    }

    cv_values <- matrix(
      Inf,
      nrow = length(lambda1_grid),
      ncol = length(lambda2_grid)
    )

    rownames(cv_values) <- paste0("lam1_", lambda1_grid)
    colnames(cv_values) <- paste0("lam2_", lambda2_grid)

    for(i1 in seq_along(lambda1_grid)){

      lam1 <- lambda1_grid[i1]

      for(i2 in seq_along(lambda2_grid)){

        lam2 <- lambda2_grid[i2]

        if(verbose){

          cat("\n============================\n")
          cat(
            "lambda1 =", lam1,
            "| lambda2 =", lam2, "\n"
          )
          cat("============================\n")
        }

        cv_error <- 0
        valid <- TRUE

        for(k in unique(folds)){

          train_idx <- which(folds != k)
          test_idx  <- which(folds == k)

          if(length(train_idx) < 1){
            valid <- FALSE
            break
          }

          res <- tryCatch({

            solve_simplex_lp(
              P_list[train_idx],
              P_prime_list[train_idx],
              a_weights,
              lambda1 = lam1,
              lambda2 = lam2
            )

          }, error = function(e) NULL)

          if(is.null(res)){
            valid <- FALSE
            break
          }

          A_hat <- res$A

          for(t in test_idx){

            err <- compute_fit_single(
              A_hat,
              P_list[[t]],
              P_prime_list[[t]]
            )

            if(is.infinite(err)){
              valid <- FALSE
              break
            }

            cv_error <- cv_error + err
          }

          if(!valid) break
        }

        if(valid){

          cv_values[i1, i2] <- cv_error

          if(verbose){
            cat("CV error =", cv_error, "\n")
          }

        } else {

          cv_values[i1, i2] <- Inf

          if(verbose){
            cat("Configuration discarded\n")
          }
        }
      }
    }

    if(all(is.infinite(cv_values))){
      stop("All configurations were discarded.")
    }

    idx <- which(
      cv_values == min(cv_values, na.rm = TRUE),
      arr.ind = TRUE
    )[1,]

    return(list(
      method = "cv",
      best_lambda1 = lambda1_grid[idx[1]],
      best_lambda2 = lambda2_grid[idx[2]],
      cv_values = cv_values,
      lambda1_grid = lambda1_grid,
      lambda2_grid = lambda2_grid,
      K = K,
      seed = seed
    ))
  }

  # =====================================================
  # GCV
  # =====================================================
  if(method == "gcv"){

    n1 <- length(lambda1_grid)
    n2 <- length(lambda2_grid)

    n_total <- n1 * n2

    fit_values <- rep(Inf, n_total)
    score_values <- rep(Inf, n_total)

    cc <- 0

    for(i1 in seq_along(lambda1_grid)){

      lam1 <- lambda1_grid[i1]

      for(i2 in seq_along(lambda2_grid)){

        lam2 <- lambda2_grid[i2]

        cc <- cc + 1

        res <- tryCatch(

          solve_simplex_lp(
            P_list,
            P_prime_list,
            a_weights,
            lambda1 = lam1,
            lambda2 = lam2
          ),

          error = function(e) NULL
        )

        if(is.null(res)){
          next
        }

        A_hat <- res$A

        fit <- compute_fit(
          A_hat,
          P_list,
          P_prime_list
        )

        fit_values[cc] <- fit

        svd_res <- svd(A_hat)

        sigma <- svd_res$d
        sigma <- sigma[sigma > 1e-6]

        df <- sum(
          sigma^2 /
            (sigma^2 + lam1 + lam2)
        )

        n <- length(P_list) *
          (nrow(A_hat) - 1)

        score_values[cc] <-
          fit / (max(n - df, 1)^2)
      }
    }

    best_idx <- which.min(score_values)

    best_i1 <- (best_idx - 1) %/% n2 + 1
    best_i2 <- (best_idx - 1) %% n2 + 1

    return(list(
      method = "gcv",
      best_lambda1 = lambda1_grid[best_i1],
      best_lambda2 = lambda2_grid[best_i2],
      score_values = score_values,
      fit_values = fit_values,
      lambda1_grid = lambda1_grid,
      lambda2_grid = lambda2_grid
    ))
  }
}
                          
OCC <- function(P, Pprime) {
  
  # check dimensioni
  if(nrow(P) != nrow(Pprime)) {
    stop("The 2 matrices must have the same numbers of observations")
  }
  
  N <- nrow(P)
  
  # center of mass
  cmass_P <- as.vector(P %*% (1:ncol(P)))
  cmass_Pprime <- as.vector(Pprime %*% (1:ncol(Pprime)))
  
  # Spearman correlation
  occ <- cor(cmass_P, cmass_Pprime, method = "spearman")
  
  return(occ)
}


tensor_product <- function(comps, weights = NULL, return_indices = FALSE) {
  K <- length(comps)
  if (is.null(weights)) {
    weights <- lapply(comps, function(x) rep(1, length(x)-1))
  }
  
  # calcola costi cumulativi
  cumcosts <- lapply(1:K, function(k) {
    w <- weights[[k]]
    if(length(w) != length(comps[[k]])-1) stop("Weights must have length = length(composition)-1")
    cumsum(c(0, w))
  })
  
  # tutte le combinazioni di indici
  index_grid <- do.call(expand.grid, lapply(comps, function(x) seq_along(x)))
  
  # calcola prodotto e costo
  costs <- apply(index_grid, 1, function(idx) sum(sapply(1:K, function(k) cumcosts[[k]][idx[k]])))
  values <- apply(index_grid, 1, function(idx) prod(sapply(1:K, function(k) comps[[k]][idx[k]])))
  
  # ordine: prima per costi, poi lessicografico
  ord <- do.call(order, c(list(costs), as.list(index_grid)))
  
  out <- list(
    product = values[ord]
  )
  if(return_indices) out$indices <- index_grid[ord, , drop=FALSE]
  return(out)
}


calc_wrps <- function(p_hat, y_true){
  Fhat <- cumsum(p_hat)
  Fobs <- cumsum(y_true)
  sum((Fhat - Fobs)^2)
}

# Brier multicategoria
calc_brier <- function(p_hat, y_true){
  sum((p_hat - y_true)^2)
}

# Wasserstein W1 per ordinale
calc_wass <- function(p_hat, y_true){
  Fhat <- cumsum(p_hat)
  Fobs <- cumsum(y_true)
  sum(abs(Fhat - Fobs))
}                                                                          

#Grid
lambda_grid <- c(0, 0.005, 0.01, 0.02, 0.05, 0.08, 0.1)
                                                                              

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

compute_R2 <- function(Y, X, A, weights){

  N <- nrow(Y)
  Cy <- ncol(Y)
  
  # --- funzione mediana Wasserstein ---
  compute_wmedian <- function(mat, weights){
    distot <- rep(0, nrow(mat))
    
    for(i in 1:nrow(mat)){
      for(j in 1:nrow(mat)){
        if(i != j){
          if(opiwd(weights, mat[i,], mat[j,]) == 1){
            distot[i] <- distot[i] + 1
          }
        }
      }
    }
    
    ord <- order(distot)
    
    if(nrow(mat) %% 2 == 1){
      return(mat[ord[round(nrow(mat)/2)], ])
    } else {
      return((mat[ord[nrow(mat)/2], ] +
              mat[ord[nrow(mat)/2 + 1], ]) / 2)
    }
  }
  
  ymed <- compute_wmedian(Y, weights)
  
  pred <- matrix(0, nrow = N, ncol = Cy)
  for(i in 1:N){
    pred[i, ] <- as.vector(A %*% X[i, ])
  }
  
  SSE <- 0
  SST <- 0
  
  for(i in 1:N){
    SSE <- SSE + wd(weights, Y[i, ], pred[i, ])
    SST <- SST + wd(weights, Y[i, ], ymed)
  }
  
  # --- R2 ---
  if(SST == 0){
    return(NA)  # evita divisione per zero
  }
  
  R2 <- 1 - SSE / SST
  
  return(list(
    R2 = R2,
    SSE = SSE,
    SST = SST,
    y_median = ymed,
    pred = pred
  ))
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

  # mat_col ha colonne: High, Medium, Low
  High   <- mat_col[, 1]
  Medium <- mat_col[, 2]
  Low    <- mat_col[, 3]

  # mapping coerente con:
  # High -> top
  # Medium -> left
  # Low -> right
  df_xy <- data.frame(
    x = Low + 0.5 * High,
    y = (sqrt(3)/2) * High
  )

  df_xy$Column <- paste0("C", col_idx)

  return(df_xy)
}

################sim functions###############################
generate_ordinal_data <- function(
    N = 200,
    Cx = 5,
    Cy = 8,
    dgp_type = "laplace",      # "laplace", "shift", "gaussian"
    lambda = 1.5,              # usato in laplace
    sigma = 1.2,               # usato in gaussian
    alpha_dir = 60,            # concentrazione Dirichlet
    equidistance = TRUE,       # default supporto equispaziato
    gamma = 1.5,               # usato se equidistance = FALSE
    latent_thresholds = NULL   # nuovo: vettore per salti non equispaziati
){
  
  # ----- coordinate ordinali -----
  zx <- 1:Cx
  
  if(!is.null(latent_thresholds)){
    if(length(latent_thresholds) != Cy){
      stop("latent_thresholds must have length Cy")
    }
    zy <- latent_thresholds
  } else if(equidistance){
    zy <- 1:Cy
  } else {
    zy <- (1:Cy)^gamma
  }
  
  # ----- inizializzo matrici -----
  x0 <- matrix(0, N, Cx)
  y  <- matrix(0, N, Cy)
  
  for(n in 1:N){
    # X generato da Dirichlet uniforme
    x0[n,] <- rdirichlet(1, rep(1, Cx))
    
    # score continuo latente
    s <- sum(zx * x0[n,])
    
    if(dgp_type == "laplace"){
      prob <- exp(-lambda * abs(zy - s))
      prob <- prob / sum(prob)
      y[n,] <- rdirichlet(1, alpha_dir * prob)
      
    } else if(dgp_type == "shift"){
      s_round <- round(sum((1:Cx) * x0[n,]))
      shift <- sample(c(-2,-1,0,1,2),
                      1,
                      prob = c(0.05,0.2,0.5,0.2,0.05))
      target <- min(max(s_round + shift, 1), Cy)
      #prob <- exp(-abs((1:Cy) - target))#kernel esponenziale
      sigma_loc <- 0.8
      prob <- exp(-((1:Cy - target)^2)/(2*sigma_loc^2))
      prob <- prob / sum(prob)
      y[n,] <- rdirichlet(1, alpha_dir * prob + 1)
      
    } else if(dgp_type == "gaussian"){
      prob <- exp(-((zy - s)^2)/(2*sigma^2))
      prob <- prob / sum(prob)
      y[n,] <- rdirichlet(1, alpha_dir * prob)
      
    } else {
      stop("dgp_type must be 'laplace', 'shift', or 'gaussian'")
    }
  }
  
  return(list(X = x0, Y = y))
}


generate_ordinal_data_2x <- function(
    N = 200,
    Cx1 = 4,
    Cx2 = 4,
    Cy = 6,
    dgp_type = "gaussian",
    lambda = 2,
    sigma = 1.2,
    alpha_dir = 80,
    equidistance = TRUE,
    gamma = 1.6,
    w1 = 0.6,        
    w2 = 0.4,        
    w2_y1 = 0.05     
){
  
  if(equidistance){
    zx1 <- 1:Cx1
    zx2 <- 1:Cx2
    zy  <- 1:Cy
  } else {
    zx1 <- (1:Cx1)^gamma
    zx2 <- (1:Cx2)^gamma
    zy  <- (1:Cy)^gamma
  }
  
  X1 <- matrix(0, N, Cx1)
  X2 <- matrix(0, N, Cx2)
  Y1 <- matrix(0, N, Cy)
  Y2 <- matrix(0, N, Cy)
  Y3 <- matrix(0, N, Cy)
  
  build_prob <- function(s){
    if(dgp_type == "gaussian"){
      p <- exp(-((zy - s)^2)/(2*sigma^2))
    } else if(dgp_type == "laplace"){
      p <- exp(-lambda * abs(zy - s))
    } else if(dgp_type == "shift"){
      s_round <- round(s)
      #s_round <- round(sum((1:Cx) * x0[n,]))
      
      shift <- sample(c(-2,-1,0,1,2),
                      1,
                      prob=c(0.05,0.2,0.5,0.2,0.05))
      
      target <- min(max(s_round + shift, 1), Cy)
      
      #prob <- rep(0, Cy)
      #prob[target] <- 1
      p <- exp(-abs((1:Cy) - target))#kernel esponenziale
      #sigma_loc=0.8
      #prob <- exp(-((1:Cy - target)^2)/(2*sigma_loc^2)) #kernel gaussiano
      #prob <- prob / sum(prob)
    }
    p/sum(p)
  }
  
  for(n in 1:N){
    
    X1[n,] <- rdirichlet(1, rep(1, Cx1))
    X2[n,] <- rdirichlet(1, rep(1, Cx2))
    
    s1 <- sum(zx1 * X1[n,])
    s2 <- sum(zx2 * X2[n,])
    
    # -------------------------
    # Y1 → X1 + X2 (X2 almonst irrilevant)
    # -------------------------
    s_y1 <- (1 - w2_y1) * s1 + w2_y1 * s2
    
    # -------------------------
    # Y2 → X1 + X2 both relevant
    # -------------------------
    s_y2 <- w1 * s1 + w2 * s2
    
    # -------------------------
    # Y3 → SOLO X1
    # -------------------------
    s_y3 <- s1
    
    prob1 <- build_prob(s_y1)
    prob2 <- build_prob(s_y2)
    prob3 <- build_prob(s_y3)
    
    Y1[n,] <- rdirichlet(1, alpha_dir * prob1 + 1e-6)
    Y2[n,] <- rdirichlet(1, alpha_dir * prob2 + 1e-6)
    Y3[n,] <- rdirichlet(1, alpha_dir * prob3 + 1e-6)
  }
  
  return(list(
    X1 = X1,
    X2 = X2,
    Y1 = Y1,
    Y2 = Y2,
    Y3 = Y3
  ))
}


