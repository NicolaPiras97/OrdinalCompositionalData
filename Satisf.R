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

satisf <- as.matrix(quality[,14:17]) 
y <- satisf/rowSums(satisf)
safety <- as.matrix(quality[,22:25]) 
x1 <- safety/rowSums(safety)
air<- as.matrix(quality[,10:13]) 
x2 <- air/rowSums(air)

quality <- read.table("Data_life_quality.txt", sep=",", header=T)
satisf <- as.matrix(quality[,14:17]); y <- satisf/rowSums(satisf)
safety <- as.matrix(quality[,22:25]); x1 <- safety/rowSums(safety)
air <- as.matrix(quality[,10:13]); x2 <- air/rowSums(air)

N <- dim(y)[1]; Cy <- dim(y)[2]; Cx1 <- dim(x1)[2]; Cx2 <- dim(x2)[2]
xdata1 <- list()
xdata2 <- list()
ydata <- list(); Z_data <- list()
for(i in 1:N){ 
  xdata1[[i]]<-x1[i,]
  xdata2[[i]]<-x2[i,]
  ydata[[i]] <- y[i,]
  Z_data[[i]] <- tensor_product_ordered(x1[i,], x2[i,])$product
}

ycdf_mat <- t(apply(y, 1, cumsum))
weights_orig <- numeric(Cy - 1)
for(k in 1:(Cy-1)) weights_orig[k] <- median(ycdf_mat[,k+1]) - median(ycdf_mat[,k])
if(any(weights_orig <= 0)) weights_orig <- rep(1, Cy-1)


# LEAVE-ONE-OUT CROSS-VALIDATION (LOOCV)
cv_results <- do_loocv(
  Z = Z_data, 
  Y = ydata, 
  w = weights_orig, 
  dims = c(Cx1, Cx2), 
  lambda_seq = lambda_grid
)


best_idx <- which.min(cv_results$error)
min_lambda <- cv_results$lambda[best_idx]
min_err <- cv_results$error[best_idx]
target_err <- min_err + cv_results$se[best_idx]

# Plot
plot(cv_results$lambda, cv_results$error, type="b", pch=19, main="LOOCV Error Curve",
     xlab="Lambda", ylab="Mean Wasserstein Error", ylim=c(min(cv_results$error)*0.9, max(cv_results$error)*1.1))
abline(v=min_lambda, col="red", lty=2)
abline(h=target_err, col="gray", lty=3)
legend("topright", legend=c("Min"), col=c("red"), lty=2)


FINAL_LAMBDA <- min_lambda 

cat(paste("\n\nLambda =", FINAL_LAMBDA, "...\n"))
final_model <- solve_simplex_lp(Z_data, ydata, weights_orig, lambda = FINAL_LAMBDA, dims_grid = c(Cx1, Cx2))
print(round(final_model$A, 4))
A_hat<-final_model$A

#########################################################
# Indexes
#########################################################

ymedianw<-vector()
for(j in 1:Cy){
  ymedianw[j]=0
}

distot<-rep(0,N)
for(i in 1:N){
  for(j in 1:N){
    if(i!=j){
      if(porwd(weights_orig,y[i,],y[j,])==2){
        distot[i]=distot[i]+1
      }
    }
  }
}
i=1
ordtot<-vector()
for(j in 1:N){
  if(length(which(distot==sort(distot,decreasing = F)[i]))==1){
    ordtot[i]<-which(distot==sort(distot,decreasing = F)[i])
  }
  if(length(which(distot==sort(distot,decreasing = F)[i]))>1){
    ordtot[i:(i+length(which(distot==sort(distot,decreasing = F)[i]))-1)]<-which(distot==sort(distot,decreasing = F)[i])
  }
  i=i+length(which(distot==sort(distot,decreasing = F)[i]))
  if(length(ordtot==N)){
    stop
  }
}
if(N%%2==1){
  ymedianw<-y[ordtot[round(N/2)],]
}
if(N%%2==0){
  ymedianw<-apply(cbind(y[ordtot[N/2],],y[ordtot[N/2+1],]),1,mean)
}


denw<-0
for(i in 1:N){
  denw=denw+wd(weights_orig,ymedianw,ydata[[i]])
}

indporwd=0
for(i in 1:(N-1)){
  for(j in (1+i):N){
    if(porwd(weights_orig,y[i,],y[j,])==porwd(weights_orig,A_hat%*%tensor_product_ordered(xdata1[[i]],xdata2[[i]])$product,A_hat%*%tensor_product_ordered(xdata1[[j]],xdata2[[j]])$product)){
      indporwd=indporwd+1
    }
  }
}
OPIwd=indporwd/choose(N,2)
cat("\n--- ORDINAL PRESERVATION INDEX (OPI) ---\n")
print(OPIwd*((N-2*(Cy-1))/N))

SSRw<-0
for(i in 1:N){
  SSRw=SSRw+wd(weights_orig,ymedianw,A_hat%*%tensor_product_ordered(xdata1[[i]],xdata2[[i]])$product)
}
R2W<-SSRw/denw
cat("\n--- WASSERSTEIN R-SQUARED ---\n")
print(R2W*((N-2*(Cy-1))/N))

# ==============================================================================
# 3. BOOTSTRAP 
# ==============================================================================

B_iter <- 20000
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
  ycdfB <- t(apply(y_mat_B, 1, cumsum))
  weightsB <- c(median(ycdfB[,2])-median(ycdfB[,1]), median(ycdfB[,3])-median(ycdfB[,2]), median(ycdfB[,4])-median(ycdfB[,3]))
  
  # Solver
  Z_dataB <-lapply(1:N, function(i) tensor_product_ordered(xdataB1[[i]],xdataB2[[i]] )$product)
  solB <- solve_simplex_lp(Z_dataB, ydataB, weightsB)
  A_boot <- solB$A
  A_boot_list[[b]] <- A_boot
  
  distances_W[b] <- matrix_wasserstein_dist(A_hat, A_boot, weights_orig)
  setTxtProgressBar(pb, b)
}
close(pb)


radius_95 <- quantile(distances_W, 0.95)
cat("\n\nRadius R95:", radius_95, "\n")

valid_indices <- which(distances_W <= radius_95)
valid_boots <- A_boot_list[valid_indices]

A_Freshet_W <- matrix(0, nrow=Cy, ncol=Cx1*Cx2)
dimnames(A_Freshet_W) <- dimnames(A_hat)

for(j in 1:(Cx1*Cx2)) {
  cols_j_list <- lapply(valid_boots, function(M) M[, j])
  mat_cols_j  <- do.call(cbind, cols_j_list) 
  A_Freshet_W[, j] <- rowMeans(mat_cols_j)
}
cat("\n--- FRÉCHET MEAN ---\n")
print(round(A_Freshet_W, 4))
