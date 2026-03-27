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
    gamma = 1.6
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
