source(file="Functions_simplex.R")
quality<-read.table("Data_life_quality.txt",sep=",",header=T)
#Variable names
#4 modalities VU:Very unsatisfied; RU:Rather unsatisfied; RS=Rather satisfied; VS=Very satisfied
#City
#DH:Health care services, doctors and hospitals.
#GS:Green spaces such as parks and gardens.
#QA:The quality of the air
#SL:I'm satisfied to live in my city.
#GJ:It is easy to find a good job in my city.
#SN:I feel safe walking alone at night in my city.
#GH:It is easy to find good housing in my city at a reasonable price.
#A:Age 15-24;25-39;40-54;55+

satisf <- as.matrix(quality[,14:17]); ya <- satisf/rowSums(satisf)
safety <- as.matrix(quality[,22:25]); x1a <- safety/rowSums(safety)
air <- as.matrix(quality[,10:13]); x2a <- air/rowSums(air)
green <- as.matrix(quality[,6:9]); x3a <- green/rowSums(green)

N <- dim(ya)[1]; Cy <- dim(ya)[2]; Cx1 <- dim(x1a)[2]; Cx2 <- dim(x2a)[2]
y<-matrix(0,nrow=N,ncol=Cy)
x1<-matrix(0,nrow=N,ncol=Cx1)
x2<-matrix(0,nrow=N,ncol=Cx2)
for(i in 1:N){
  y[i,1]<-ya[i,4];y[i,2]<-ya[i,3];y[i,3]<-ya[i,2];y[i,4]<-ya[i,1]
  x1[i,1]<-x1a[i,4];x1[i,2]<-x1a[i,3];x1[i,3]<-x1a[i,2];x1[i,4]<-x1a[i,1]
  x2[i,1]<-x2a[i,4];x2[i,2]<-x2a[i,3];x2[i,3]<-x2a[i,2];x2[i,4]<-x2a[i,1]
}

xdata1 <- list()
xdata2 <- list()
ydata <- list(); Z_data <- list()

weights<-select_weights(x1, x2, ydata, a, b, lambda_grid)$best$weights_raw

for(i in 1:N){ 
       xdata1[[i]]<-x1[i,]
       xdata2[[i]]<-x2[i,]
       ydata[[i]] <- y[i,]
       Z_data[[i]] <- tensor_product(comps=list(x1[i,], x2[i,]),list(a,b),return_indices = TRUE)$product
   }
res_svd <- select_lambda_loocv(Z_data, ydata, weights, lambda_grid)
#res_svd$best_lambda  # λ1 selezionato con GCV basato su SVD
final_model <- solve_simplex_lp(Z_data, ydata, weights, lambda = res_svd$best_lambda)
A_hat<-final_model$A
indices_map<-tensor_product(comps=list(x1[1,], x2[1,]),list(a,b),return_indices=TRUE)$indices
colnames(A_hat) <- apply(indices_map, 1, function(v) paste0("(", v[1], ",", v[2], ")"))
print(round(A_hat, 4))


#########################################################
# Indexes
#########################################################

cat("\n--- OCC(x1,y) ---\n")
OCC(x1,y)
cat("\n--- OCC(x2,y) ---\n")
OCC(x2,y)
z <- do.call(rbind, Z_data)
cat("\n--- OCC(x1*x2,y) ---\n")                         
OCC(z,y)

indporwd=0
 for(i in 1:(N-1)){
     for(j in (1+i):N){
         if(opiwd(weights,y[i,],y[j,])==opiwd(weights,A_hat%*%Z_data[[i]],A_hat%*%Z_data[[j]])){
             indporwd=indporwd+1
         }
     }
 }
OPIwd=indporwd/choose(N,2)
cat("\n--- ORDINAL PRESERVATION INDEX (OPI) ---\n")
print(OPIwd)

R2W<-compute_R2(y,z,A_hat,weights)$R2
cat("\n--- WASSERSTEIN R-SQUARED ---\n")
print(R2W)

# ==============================================================================
# 3. BOOTSTRAP 
# ==============================================================================

B_iter <- 1000
distances_W <- numeric(B_iter)
A_boot_list <- list()

cat(paste(" Bootstrap (", B_iter, " iteration)...\n", sep=""))
pb <- txtProgressBar(min = 0, max = B_iter, style = 3)
set.seed(123) 

for(b in 1:B_iter) {
  indices <- sample(1:N, N, replace=TRUE)
  xdataB1 <- list();xdataB2 <- list(); ydataB <- list()
  
  for(k in 1:N) {
    idx <- indices[k]
    xdataB1[[k]] <- xdata1[[idx]]
    xdataB2[[k]] <- xdata2[[idx]]
    ydataB[[k]] <- ydata[[idx]]
  }
  
  y_mat_B <- do.call(rbind, ydataB)
 
  weightsB <- weights
  
  # Solver
  Z_dataB <-lapply(1:N, function(i) tensor_product(comps=list(xdataB1[[i]],xdataB2[[i]]),list(a,b) )$product)
  res_svdB <- select_lambda_loocv(Z_dataB, ydataB, weightsB, lambda_grid)
  solB <- solve_simplex_lp(Z_dataB, ydataB, weightsB, lambda = res_svdB$best_lambda)
  A_boot <- solB$A
  A_boot_list[[b]] <- A_boot
  
  distances_W[b] <- matrix_wasserstein_dist(A_hat, A_boot, weights_orig)
  setTxtProgressBar(pb, b)
}
close(pb)


radius_95 <- quantile(distances_W, 0.95)
cat("\n\nRadius R95:", radius_95/sum(weights), "\n")

valid_indices <- which(distances_W <= radius_95)
valid_boots <- A_boot_list[valid_indices]

A_Freshet_W <- matrix(0, nrow=Cy, ncol=Cx1*Cx2)
dimnames(A_Freshet_W) <- dimnames(A_hat)

for(j in 1:(Cx1*Cx2)) {
  cols_j_list <- lapply(valid_boots, function(M) M[, j])
  mat_cols_j  <- do.call(cbind, cols_j_list) 
  A_Freshet_W[, j] <- rowMeans(mat_cols_j)
}
colnames(A_Freshet_W) <- apply(indices_map, 1, function(v) paste0("(", v[1], ",", v[2], ")"))
cat("\n--- FRÉCHET MEAN ---\n")
print(round(A_Freshet_W, 4))
