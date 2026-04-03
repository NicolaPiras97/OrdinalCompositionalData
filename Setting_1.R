source(file="Functions_simplex.R")

run_simulation <- function(
    iter = 1000,
    N = 100,
    Cx = 5,
    Cy = 8,
    dgp_type = "laplace",
    lambda = 2,
    sigma = 1.2,
    alpha_dir = 80,
    train_ratio = 0.7,
    equidistance = TRUE,    
    gamma = 1.6,
    prop_zero = 0
){
  
  result <- matrix(0, nrow = iter, ncol = 14)
  colnames(result) <- c(
    "WRPS_OT_OT","WRPS_COD_OT",
    "Brier_OT_OT","Brier_COD_OT",
    "WRPS_OT_COD","WRPS_COD_COD",
    "Brier_OT_COD","Brier_COD_COD",
    "Win_WRPS_OT_DGP","Win_Brier_OT_DGP",
    "Win_WRPS_COD_DGP","Win_Brier_COD_DGP",
    "Win_WRPS_TRUE_DGP","Win_Brier_TRUE_DGP"
  )
  
  for(it in 1:iter){
    
    # --- GENERAZIONE DGP ---
    data <- generate_ordinal_data(
      N = N,
      Cx = Cx,
      Cy = Cy,
      dgp_type = dgp_type,
      lambda = lambda,
      sigma = sigma,
      alpha_dir = alpha_dir,
      equidistance = equidistance,
      gamma = gamma
    )
    
    x0 <- data$X
    y  <- data$Y
    
    # COD matrix from DGP
    Ac <- codalm(y, x0)
    
    weights <- rep(1, Cy-1)
    
    # OT matrix from DGP
    ydata <- lapply(1:N, function(i) y[i,])
    xdata <- lapply(1:N, function(i) x0[i,])
    sol   <- solve_simplex_lp(xdata, ydata, weights)
    mat   <- sol$A
    
    # --- GENERO NUOVO CAMPIONE ---
    x  <- matrix(0,N,Cx)
    y1 <- matrix(0,N,Cy)
    y2 <- matrix(0,N,Cy)
    
    for(n in 1:N){
      x[n,]  <- rdirichlet(1, rep(1,Cx))
      y1[n,] <- safe_dirichlet(mat %*% x[n,])
      y2[n,] <- safe_dirichlet(t(Ac) %*% x[n,])
    }
    
    # --- TRAIN/VALID SPLIT ---
    tr <- trunc(train_ratio * N)
    va <- N - tr
    
    idx  <- sample(1:N, N)
    intr <- idx[1:tr]
    inva <- idx[(tr+1):N]

    # --- FORZO PERCENTUALE DI ZERI NEL VALIDATION SET ---
    if(prop_zero > 0){
      total_cells <- length(inva) * Cy
      n_zero <- floor(prop_zero * total_cells)
      if(n_zero > 0){
        val_matrix <- expand.grid(row = inva, col = 1:Cy)
        
        # y1
        values_y1 <- mapply(function(r,c) y1[r,c], val_matrix$row, val_matrix$col)
        ord_idx <- order(values_y1)
        selected <- val_matrix[ord_idx[1:n_zero], ]
        for(k in 1:nrow(selected)){
          y1[selected$row[k], selected$col[k]] <- 0
        }
        # y2
        values_y2 <- mapply(function(r,c) y2[r,c], val_matrix$row, val_matrix$col)
        ord_idx <- order(values_y2)
        selected <- val_matrix[ord_idx[1:n_zero], ]
        for(k in 1:nrow(selected)){
          y2[selected$row[k], selected$col[k]] <- 0
        }
        # rinormalizzazione
        for(r in inva){
          if(sum(y1[r,]) > 0) y1[r,] <- y1[r,] / sum(y1[r,])
          if(sum(y2[r,]) > 0) y2[r,] <- y2[r,] / sum(y2[r,])
        }
      }
    }
                    
    
    Plist  <- lapply(1:N, function(i) x[i,])
    P1list <- lapply(1:N, function(i) y1[i,])
    P2list <- lapply(1:N, function(i) y2[i,])
    
    Plist_tr  <- lapply(intr, function(i) x[i,])
    P1list_tr <- lapply(intr, function(i) y1[i,])
    P2list_tr <- lapply(intr, function(i) y2[i,])
    
    # --- STIMA SU TRAIN ---
    A1 <- solve_simplex_lp(Plist_tr, P1list_tr, weights)$A
    A2 <- solve_simplex_lp(Plist_tr, P2list_tr, weights)$A
    
    # COD stimato
    eps <- 1e-8
    Ymat1 <- do.call(rbind, P1list_tr)
    Ymat2 <- do.call(rbind, P2list_tr)
    Xmat <- do.call(rbind, Plist_tr)
    Ymat1[Ymat1 < eps] <- eps
    Ymat2[Ymat2 < eps] <- eps
    Xmat[Xmat < eps] <- eps
    Ymat1 <- Ymat1 / rowSums(Ymat1)
    Ymat2 <- Ymat2 / rowSums(Ymat2)
    Xmat <- Xmat / rowSums(Xmat)
    A1_cod <- codalm(Ymat1, Xmat)
    A2_cod <- codalm(Ymat2, Xmat)
    
    # --- ERRORI ---
    calc_errors <- function(A_ot, A_cod, y_true, x_mat){
      err_wrps_ot <- err_wrps_cod <- err_b_ot <- err_b_cod <- err_w1_ot <- err_w1_cod <- 0
      for(i in inva){
        pred_ot  <- as.vector(A_ot  %*% x_mat[i,])
        pred_cod <- as.vector(t(A_cod) %*% x_mat[i,])
        err_wrps_ot  <- err_wrps_ot  + calc_wrps(pred_ot,  y_true[i,])
        err_wrps_cod <- err_wrps_cod + calc_wrps(pred_cod, y_true[i,])
        err_b_ot     <- err_b_ot  + calc_brier(pred_ot,  y_true[i,])
        err_b_cod    <- err_b_cod + calc_brier(pred_cod, y_true[i,])
      }
      return(c(
        err_wrps_ot/va, err_wrps_cod/va,
        err_b_ot/va,    err_b_cod/va
      ))
    }
    
    errs1 <- calc_errors(A1, A2_cod, y1, x)
    errs2 <- calc_errors(A2, A1_cod, y2, x)
    result[it,1:4]  <- errs1
    result[it,5:8] <- errs2
    
    # --- WINS ---
    result[it,9:10] <- as.numeric(errs1[c(1,3)] < errs1[c(2,4)])
    result[it,11:12] <- as.numeric(errs2[c(1,3)] < errs2[c(2,4)])
    result[it,13:14] <- as.numeric(errs1[c(1,3)] < errs2[c(2,4)])
  }
  
  list(
    mean_res  = colMeans(result),
    sd_res    = apply(result,2,sd),
    win_rates = colMeans(result[,9:14]),
    raw_results = result
  )
}

s1 <- run_simulation(
  iter = 500,
  N = 50,
  Cx = 5,
  Cy = 3,
  dgp_type = "laplace"
)
s1$win_rates         

#Plot Setting 4.2.1
                        library(dplyr)
library(tidyr)
library(ggplot2)

# =========================
# DATA COMPLETI (N = 100)
# =========================

data <- data.frame(
  Cx = rep(c(3,7,5,3,5,5,7), each = 5),
  Cy = rep(c(3,7,5,5,3,7,5), each = 5),
  prop_zero = rep(c(0,0.1,0.2,0.3,0.4), 7),
  
  OT = c(
    # (3,3)
    0.530,0.492,0.550,0.704,0.802,
    # (7,7)
    0.742,0.698,0.750,0.744,0.732,
    # (5,5)
    0.579,0.602,0.574,0.636,0.656,
    # (3,5)
    0.558,0.538,0.536,0.550,0.590,
    # (5,3)
    0.548,0.544,0.542,0.668,0.804,
    # (5,7)
    0.596,0.538,0.606,0.592,0.578,
    # (7,5)
    0.730,0.738,0.718,0.754,0.754
  ),
  
  COD = c(
    # (3,3)
    0.430,0.454,0.538,0.568,0.658,
    # (7,7)
    0.654,0.592,0.654,0.638,0.650,
    # (5,5)
    0.480,0.430,0.490,0.528,0.580,
    # (3,5)
    0.382,0.358,0.414,0.350,0.462,
    # (5,3)
    0.526,0.534,0.546,0.538,0.556,
    # (5,7)
    0.446,0.460,0.474,0.448,0.486,
    # (7,5)
    0.616,0.662,0.620,0.644,0.632
  )
)

# =========================
# LONG FORMAT
# =========================

data_long <- data %>%
  pivot_longer(cols = c(OT, COD),
               names_to = "DGP",
               values_to = "WinRate")

# =========================
# PLOT
# =========================

p1<-ggplot(data_long, aes(x = prop_zero, y = WinRate,
                      color = DGP, group = DGP)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  facet_grid(Cx ~ Cy, labeller = label_both) +
  labs(
    x = "Proportion of zeros",
    y = "Win rate (RPS)",
    title = "Predictive performance under increasing zero inflation",
    color = "DGP"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")
p1

