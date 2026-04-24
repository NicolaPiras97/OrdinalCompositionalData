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
    prop_zero_max = 0,
    prop_zero_max2 = 0,
    prop_zero =0,
    prop_inv=0
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
    if(prop_zero_max > 0){

  # numero di righe da modificare
  n_rows <- length(inva)
  n_zero_rows <- floor(prop_zero_max* n_rows)

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

if(prop_zero_max2 > 0){

  # numero di righe da modificare
  n_rows <- length(inva)
  n_zero_rows <- floor(prop_zero_max2 * n_rows)

  if(n_zero_rows > 0){

    # selezione casuale delle righe nel validation set
    selected_rows <- sample(inva, n_zero_rows, replace = FALSE)

    for(r in selected_rows){

      # controllo numero colonne
      if(ncol(y) >= 5){
        # indici delle due componenti più grandi
        j_max <- order(y[r, ], decreasing = TRUE)[1:2]
      } else {
        # indice della componente più grande
        j_max <- which.max(y[r, ])
      }

      # azzera le componenti selezionate
      y[r, j_max] <- 0

      # normalizzazione se la somma è positiva
      if(sum(y[r,]) > 0){
        y[r,] <- y[r,] / sum(y[r,])
      }
    }
  }
}

if(prop_zero > 0){

  n_rows <- length(inva)
  n_zero_rows <- floor(prop_zero * n_rows)

  if(n_zero_rows > 0){

    selected_rows <- sample(inva, n_zero_rows, replace = FALSE)

    for(r in selected_rows){

      n_comp <- ncol(y)

      # quante componenti annullare
      k <- if(n_comp >= 5) 2 else 1

      # indici candidati (solo quelli > 0 per evitare inutilità)
      non_zero_idx <- which(y[r, ] > 0)

      # se ci sono abbastanza componenti non nulle
      if(length(non_zero_idx) >= k){

        j_sel <- sample(non_zero_idx, k, replace = FALSE)

      } else {
        # fallback: prendi da tutte le colonne
        j_sel <- sample(seq_len(n_comp), k, replace = FALSE)
      }

      # annulla le componenti selezionate
      y[r, j_sel] <- 0

      # rinormalizza
      s <- sum(y[r,])
      if(s > 0){
        y[r,] <- y[r,] / s
      } else {
        # fallback opzionale (evita riga tutta zero)
        y[r,] <- rep(1/ncol(y), ncol(y))
      }
    }
  }
}

 if(prop_inv > 0){

  n_rows <- length(inva)
  n_inv_rows <- floor(prop_inv * n_rows)

  if(n_inv_rows > 0){

    selected_rows <- sample(inva, n_inv_rows, replace = FALSE)

    for(r in selected_rows){

      # inversione delle componenti
      y[r, ] <- rev(y[r, ])
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


data_rps <- data.frame(
     config = rep(c("(3,3)","(7,7)","(5,5)","(3,5)","(5,3)","(5,7)","(7,5)"), each = 5),
     prop_zero = rep(c(0,0.05,0.10,0.15,0.20), 7),
     
     WD = c(
         # (3,3)
         0.0284,0.0287,0.0284,0.0286,0.0282,
         # (7,7)
         0.0420,0.0420,0.0414,0.0415,0.0421,
         # (5,5)
         0.0390,0.0388,0.0380,0.0387,0.0382,
         # (3,5)
         0.0321,0.0317,0.0320,0.0320,0.0323,
         # (5,3)
         0.0215,0.0217,0.0217,0.0213,0.0219,
         # (5,7)
         0.0396,0.0392,0.0392,0.0384,0.0385,
         # (7,5)
         0.0385,0.0379,0.0378,0.0373,0.0375
     ),
     
     COD = c(
         # (3,3)
         0.0286,0.0288,0.0288,0.0293,0.0299,
         # (7,7)
         0.0512,0.0518,0.0519,0.0532,0.0539,
         # (5,5)
         0.0410,0.0419,0.0428,0.0426,0.0432,
         # (3,5)
         0.0328,0.0329,0.0330,0.0336,0.0351,
         # (5,3)
         0.0239,0.0247,0.0257,0.0259,0.0270,
         # (5,7)
         0.0424,0.0422,0.0433,0.0433,0.0442,
         # (7,5)
         0.0447,0.0455,0.0462,0.0468,0.0478
     )
 )
 
 # long format
 data_long <- melt(data_rps, id.vars = c("config","prop_zero"),
                   variable.name = "Method", value.name = "RPS")
 
 # plot
 p1<-ggplot(data_long, aes(x = prop_zero, y = RPS, color = Method)) +
     geom_line(size = 1.2) +
     geom_point(size = 2) +
     facet_wrap(~config, scales = "free_y") +
     labs(
         title = "Mean RPS with Zeri inflation (N = 100)",
         x = "Proportion of zeros",
         y = "Mean RPS"
     ) +
     theme_minimal()+
     scale_color_manual(values = c("WD" = "red", "COD" = "blue"))
 p1


data_rps <- data.frame(
     config = rep(c("(3,3)","(7,7)","(5,5)","(3,5)","(5,3)","(5,7)","(7,5)"), each = 5),
     prop_zero = rep(c(0,0.05,0.10,0.15,0.20), 7),
     
     WD = c(
         # (3,3)
         0.0284,0.0289,0.0292,0.0296,0.0305,
         # (7,7)
         0.0420,0.0429,0.0445,0.0449,0.0459,
         # (5,5)
         0.0390,0.0392,0.0397,0.0403,0.0399,
         # (3,5)
         0.0321,0.0325,0.0317,0.0329,0.0326,
         # (5,3)
         0.0215,0.0219,0.0221,0.0226,0.0223,
         # (5,7)
         0.0396,0.0397,0.0410,0.0412,0.0419,
         # (7,5)
         0.0385,0.0374,0.0385,0.0389,0.0376
     ),
     
     COD = c(
         # (3,3)
         0.0286,0.0298,0.0311,0.0324,0.0343,
         # (7,7)
         0.0512,0.0554,0.0608,0.0644,0.0693,
         # (5,5)
         0.0410,0.0438,0.0470,0.0494,0.0541,
         # (3,5)
         0.0328,0.0355,0.0394,0.0420,0.0480,
         # (5,3)
         0.0239,0.0259,0.0296,0.0345,0.0391,
         # (5,7)
         0.0424,0.0448,0.0515,0.0566,0.0614,
         # (7,5)
         0.0447,0.0481,0.0522,0.0557,0.0612
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
         title = "Mean RPS with Zeri inflation (N = 100)",
         x = "Proportion of zeros",
         y = "Mean RPS"
     ) +
     theme_minimal()+
     scale_color_manual(values = c("WD" = "red", "COD" = "blue"))
 p2


data_rps <- data.frame(
     config = rep(c("(3,3)","(7,7)","(5,5)","(3,5)","(5,3)","(5,7)","(7,5)"), each = 5),
     prop_inv = rep(c(0,0.1,0.2,0.3,0.4), 7),
     
     WD = c(
         # (3,3)
         0.0284,0.0282,0.0279,0.0282,0.0275,
         # (7,7)
         0.0420,0.0420,0.0420,0.0411,0.0397,
         # (5,5)
         0.0390,0.0384,0.0380,0.0369,0.0367,
         # (3,5)
         0.0321,0.0320,0.0323,0.0322,0.0322,
         # (5,3)
         0.0215,0.0217,0.0217,0.0222,0.0219,
         # (5,7)
         0.0396,0.0395,0.0394,0.0392,0.0391,
         # (7,5)
         0.0385,0.0381,0.0382,0.0379,0.0378
     ),
     
     COD = c(
         # (3,3)
         0.0286,0.0286,0.0290,0.0296,0.0299,
         # (7,7)
         0.0512,0.0526,0.0538,0.0542,0.0577,
         # (5,5)
         0.0410,0.0413,0.0425,0.0433,0.0443,
         # (3,5)
         0.0328,0.0365,0.0430,0.0506,0.0600,
         # (5,3)
         0.0239,0.0268,0.0324,0.0384,0.0451,
         # (5,7)
         0.0424,0.0470,0.0528,0.0597,0.0691,
         # (7,5)
         0.0447,0.0506,0.0563,0.0632,0.0706
     )
 )
 # long format
 data_long <- melt(data_rps, id.vars = c("config","prop_inv"),
                   variable.name = "Method", value.name = "RPS")
 
 # plot
 p3<-ggplot(data_long, aes(x = prop_inv, y = RPS, color = Method)) +
     geom_line(size = 1.2) +
     geom_point(size = 2) +
     facet_wrap(~config, scales = "free_y") +
     labs(
         title = "Mean RPS with Zeri inflation (N = 100)",
         x = "Proportion of inversions",
         y = "Mean RPS"
     ) +
     theme_minimal()+
     scale_color_manual(values = c("WD" = "red", "COD" = "blue"))
