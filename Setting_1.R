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
    "RPS_OT_OT","RPS_COD_OT",
    "Brier_OT_OT","Brier_COD_OT",
    "RPS_OT_COD","RPS_COD_COD",
    "Brier_OT_COD","Brier_COD_COD",
    "Win_RPS_OT_DGP","Win_Brier_OT_DGP",
    "Win_RPS_COD_DGP","Win_Brier_COD_DGP",
    "Win_RPS_TRUE_DGP","Win_Brier_TRUE_DGP"
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

    # --- TRAIN/VALID SPLIT ---
    tr <- trunc(train_ratio * N)
    va <- N - tr
    
    idx  <- sample(1:N, N)
    intr <- idx[1:tr]
    inva <- idx[(tr+1):N]

    # --- FORZO PERCENTUALE DI ZERI NEL VALIDATION SET ---
    if(prop_zero > 0){

  # numero di righe da modificare
  n_rows <- length(inva)
  n_zero_rows <- floor(prop_zero * n_rows)

  if(n_zero_rows > 0){

    # selezione casuale delle righe nel validation set
    selected_rows <- sample(inva, n_zero_rows, replace = FALSE)

    for(r in selected_rows){

      # indice della componente più grande
      j_max <- which.max(y[r,])

      y[r, j_max] <- 0

      if(sum(y[r,]) > 0){
        y[r,] <- y[r,] / sum(y[r,])
      }
    }
  }
}
    
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
  N = 100,
  Cx = 5,
  Cy = 5
)
s1$win_rates         


# =========================
# DATA  (N = 100)
# =========================

data <- data.frame(
  Cx = rep(c(3,7,5,3,5,5,7), each = 4),
  Cy = rep(c(3,7,5,5,3,7,5), each = 4),
  prop_zero = rep(c(0,0.2,0.4,0.6), 7),
  
  OT = c(
    # (3,3)
    0.530,0.535,0.586,0.680,
    # (7,7)
    0.742,0.774,0.758,0.814,
    # (5,5)
    0.579,0.596,0.658,0.700,
    # (3,5)
    0.538,0.541,0.569,0.622,
    # (5,3)
    0.558,0.728,0.852,0.922,
    # (5,7)
    0.596,0.608,0.672,0.728,
    # (7,5)
    0.730,0.760,0.778,0.792
  ),
  
  COD = c(
    # (3,3)
    0.430,0.506,0.560,0.638,
    # (7,7)
    0.654,0.658,0.712,0.726,
    # (5,5)
    0.480,0.518,0.556,0.640,
    # (3,5)
    0.526,0.528,0.543,0.596,
    # (5,3)
    0.382,0.652,0.848,0.919,
    # (5,7)
    0.446,0.530,0.588,0.650,
    # (7,5)
    0.616,0.644,0.690,0.745
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


data_rps <- data.frame(
  config = rep(c("(3,3)","(7,7)","(5,5)","(3,5)","(5,3)","(5,7)","(7,5)"), each = 4),
  prop_zero = rep(c(0,0.2,0.4,0.6), 7),
  
  OT = c(
    # (3,3)
    0.0284,0.0291,0.0287,0.0293,
    # (7,7)
    0.0420,0.0427,0.0436,0.0443,
    # (5,5)
    0.0390,0.0389,0.0391,0.0403,
    # (3,5)
    0.0321,0.0329,0.0336,0.0331,
    # (5,3)
    0.0215,0.0221,0.0218,0.0224,
    # (5,7)
    0.0396,0.0399,0.0395,0.0410,
    # (7,5)
    0.0385,0.0381,0.0385,0.0388
  ),
  
  COD = c(
    # (3,3)
    0.0286,0.0305,0.0320,0.0342,
    # (7,7)
    0.0512,0.0537,0.0555,0.0566,
    # (5,5)
    0.0410,0.0429,0.0444,0.0469,
    # (3,5)
    0.0328,0.0341,0.0354,0.0389,
    # (5,3)
    0.0239,0.0271,0.0325,0.0387,
    # (5,7)
    0.0424,0.0449,0.0466,0.0481,
    # (7,5)
    0.0447,0.0462,0.0484,0.0492
  )
)

# long format
data_long <- melt(data_rps, id.vars = c("config","prop_zero"),
                  variable.name = "Method", value.name = "RPS")

# plot
p2<-ggplot(data_long, aes(x = prop_zero, y = RPS, color = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~config, scales = "free_y") +
  labs(
    title = "Mean RPS with Zero Inflation (N = 100)",
    x = "Proportion of zeros",
    y = "Mean RPS"
  ) +
  theme_minimal()+
scale_color_manual(values = c("OT" = "red", "COD" = "blue"))
p2



