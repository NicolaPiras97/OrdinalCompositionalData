source(file="Functions_simplex.R")
data("educFM")
father <- as.matrix(educFM[,2:4])
y <- father/rowSums(father)
mother <- as.matrix(educFM[,5:7])
x <- mother/rowSums(mother)

N <- dim(y)[1]
Cy <- dim(y)[2]
Cx <- dim(x)[2]
ydata <- list(); xdata <- list()
for(i in 1:N){ ydata[[i]]<-y[i,]; xdata[[i]]<-x[i,] }

weights_orig <- c(1,1)

res <- select_lambda(xdata, ydata, weights_orig, lambda_grid)
final_model <- solve_simplex_lp(xdata, ydata, weights_orig, lambda = res$best_lambda)
A_hat<-final_model$A

cat("Original estimate A_hat\n")
print(round(A_hat,4))

AC<-(codalm(y,x))

# ==============================================================================
# BOOTSTRAP (Confidence regions)
# ==============================================================================

B <- 1000 
distances_W <- numeric(B)
A_boot_list <- list()
cat(paste(" Bootstrap (", B, " iteration)...\n", sep=""))
pb <- txtProgressBar(min = 0, max = B, style = 3)
set.seed(123) 

for(b in 1:B) {
  indices <- sample(1:N, N, replace=TRUE)
  xdataB <- list(); ydataB <- list()
  
  for(k in 1:N) {
    idx <- indices[k]
    xdataB[[k]] <- xdata[[idx]]
    ydataB[[k]] <- ydata[[idx]]
  }
  
  y_mat_B <- do.call(rbind, ydataB)
  ycdfB <- t(apply(y_mat_B, 1, cumsum))
  weightsB <- c(1,1)
  
  # Solver
  resB <- select_lambda(xdataB, ydataB, weightsB, lambda_grid)
  solB <- solve_simplex_lp(xdataB, ydataB, weightsB, lambda = resB$best_lambda)
  A_boot <- solB$A
  A_boot_list[[b]] <- A_boot
  
  # Distance
  distances_W[b] <- matrix_wasserstein_dist(A_hat, A_boot, weights_orig)
  setTxtProgressBar(pb, b)
}
close(pb)

radius_95 <- quantile(distances_W, 0.95)
cat("\n\nRadius R95:", radius_95/sum(weights_orig), "\n")

valid_indices <- which(distances_W <= radius_95)
valid_boots <- A_boot_list[valid_indices]

# ==================================================
# PLOT
# ==================================================

# clouds
cloud_c1 <- get_cloud_points(valid_boots, 1)
cloud_c2 <- get_cloud_points(valid_boots, 2)
cloud_c3 <- get_cloud_points(valid_boots, 3)
all_clouds <- rbind(cloud_c1, cloud_c2, cloud_c3)

# Triangol definition (Low left, Med top, High right)
triangle_border <- data.frame(
  x = c(0, 0.5, 1, 0),
  y = c(0, sqrt(3)/2, 0, 0)
)

p_clean <- ggplot() +
  geom_path(data=triangle_border, aes(x,y), color="black", linewidth=0.8) +
  
  annotate("text", x=-0.05, y=0, label="Low", fontface="bold") +
  annotate("text", x=1.05, y=0, label="High", fontface="bold") +
  annotate("text", x=0.5, y=sqrt(3)/2 + 0.05, label="Medium", fontface="bold") +
  
  geom_point(data=all_clouds, aes(x=x, y=y, color=Column), 
             size=0.8, alpha=0.5) +
  
  scale_color_manual(values=c("red", "green", "blue"), 
                     name="Confidence regions (95%):",
                     labels=c("Column 1 (Low)", "Column 2 (Medium)", "Column 3 (High)")) +
  coord_fixed() +
  theme_void() +
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
  labs(title = "Bootstrap Confidence Clouds") 

print(p_clean)

A_Freshet_W <- matrix(0, nrow=Cy, ncol=Cx)
dimnames(A_Freshet_W) <- dimnames(A_hat)

for(j in 1:Cx) {
  cols_j_list <- lapply(valid_boots, function(M) M[, j])
  mat_cols_j  <- do.call(cbind, cols_j_list) 
  A_Freshet_W[, j] <- rowMeans(mat_cols_j)
}
cat("\n--- MATRIX FRÉCHET MEDIAN (Wasserstein center) ---\n")
print(round(A_Freshet_W, 4))

#########################################################
# Indexes
#########################################################

#OCC
cat("\n--- ORDINAL COMPOSITIONAL CORRELATION (OCC) ---\n")
print(OCC(x,y))


#R2
cat("\n--- WASSERSTEIN R-SQUARED ---\n")
print(compute_R2(y,x,A_hat,weights_orig)$R2)
                        
#OPI
indopiwd=0
for(i in 1:(N-1)){
  for(j in (1+i):N){
    if(opiwd(weights_orig,y[i,],y[j,])==opiwd(weights_orig,A_hat%*%x[i,],A_hat%*%x[j,])){
      indopiwd=indopiwd+1
    }
  }
}
OPIwd=indopiwd/choose(N,2)
cat("\n--- ORDINAL PRESERVATION INDEX (OPI) ---\n")
print(OPIwd)

##############BOOTSTRAP USING WASSERSTEIN ORDER#################
weightsB<-list()
Btot<-list()
vartotc<-rep(0,Cy)
for(b in 1:B){
xB<-matrix(0,nrow=N,ncol=Cx)
yB<-matrix(0,nrow=N,ncol=Cy)
rid<-matrix(0,nrow=N,ncol=Cy)
indexB<-sample(1:N,N,replace=T)
for(i in 1:N){
xB[i,]<-x[indexB[i],]
}
for(i in 1:N){
  yB[i,]<-y[indexB[i],]
}

weightsB[[b]]<-c(1,1)

ydataB<-list()
xdataB<-list()
for(i in 1:N){
  ydataB[[i]]<-yB[i,]
}
for(i in 1:N){
  xdataB[[i]]<-xB[i,]
}
resB <- select_lambda(xdataB, ydataB, weightsB[[b]], lambda_grid)
solB<-solve_simplex_lp( xdataB , ydataB , weightsB[[b]], lambda = resB$best_lambda )
Btot[[b]]<-solB$A

for(j in 1:Cy){
  vartotc[j]=vartotc[j]+wd(weightsB[[b]],A_hat[,j],Btot[[b]][,j])
}

}

#95% confidence region
distotB<-matrix(0,nrow=B,ncol=Cy)
for(i in 1:B){
  for(j in 1:B){
    if(i!=j){
    for(k in 1:Cy){
    if(opiwd(weightsB[[b]],Btot[[i]][,k],Btot[[j]][,k])==2){
      distotB[i,k]=distotB[i,k]+1
    }
    }
    }
  }
}
ordtotB<-matrix(0,nrow=B,ncol=Cy)
for(k in 1:Cy){
  i=1
for(j in 1:B){
  if(length(which(distotB[,k]==sort(distotB[,k],decreasing = F)[i]))==1){
    ordtotB[i,k]<-which(distotB[,k]==sort(distotB[,k],decreasing = F)[i])
  }
  if(length(which(distotB[,k]==sort(distotB[,k],decreasing = F)[i]))>1){
    ordtotB[i:(i+length(which(distotB[,k]==sort(distotB[,k],decreasing = F)[i]))-1),k]<-which(distotB[,k]==sort(distotB[,k],decreasing = F)[i])
  }
  i=i+length(which(distotB[,k]==sort(distotB[,k],decreasing = F)[i]))
  if(length(ordtotB[,k]==B)){
    stop
  }
}
}
bmedianw<-matrix(0,nrow=Cx,ncol=Cy)
if(B%%2==1){
  for(k in 1:Cy){
     bmedianw[,k]<-Btot[[ordtotB[round(B/2),k]]][,k]
  }
}
if(B%%2==0){
  for(k in 1:Cy){
  bmedianw[,k]<-apply(cbind(Btot[[ordtotB[B/2,k]]][,k],Btot[[ordtotB[B/2+1,k]]][,k]),1,mean)
  }
}

binfw<-matrix(0,nrow=Cx,ncol=Cy)
if(B%%2==1){
  for(k in 1:Cy){
    binfw[,k]<-Btot[[ordtotB[round(B/40),k]]][,k]
  }
}
if(B%%2==0){
  for(k in 1:Cy){
    binfw[,k]<-apply(cbind(Btot[[ordtotB[B/40,k]]][,k],Btot[[ordtotB[B/40+1,k]]][,k]),1,mean)
  }
}
bsupw<-matrix(0,nrow=Cx,ncol=Cy)
if(B%%2==1){
  for(k in 1:Cy){
    bsupw[,k]<-Btot[[ordtotB[round(B/1.02555),k]]][,k]
  }
}
if(B%%2==0){
  for(k in 1:Cy){
    bsupw[,k]<-apply(cbind(Btot[[ordtotB[B/1.02555,k]]][,k],Btot[[ordtotB[B/1.02555+1,k]]][,k]),1,mean)
  }
}

bq1w<-matrix(0,nrow=Cx,ncol=Cy)
if(B%%2==1){
  for(k in 1:Cy){
    bq1w[,k]<-Btot[[ordtotB[round(B/4),k]]][,k]
  }
}
if(B%%2==0){
  for(k in 1:Cy){
    bq1w[,k]<-apply(cbind(Btot[[ordtotB[B/4,k]]][,k],Btot[[ordtotB[B/4+1,k]]][,k]),1,mean)
  }
}

bq3w<-matrix(0,nrow=Cx,ncol=Cy)
if(B%%2==1){
  for(k in 1:Cy){
    bq3w[,k]<-Btot[[ordtotB[round(B/1.33333),k]]][,k]
  }
}
if(B%%2==0){
  for(k in 1:Cy){
    bq3w[,k]<-apply(cbind(Btot[[ordtotB[B/1.33333,k]]][,k],Btot[[ordtotB[B/1.33333+1,k]]][,k]),1,mean)
  }
}

cat("\n--- Bootstrap median ---\n")
print(round(bmedianw, 4))
cat("\n--- Bootstrap inf ---\n")
print(round(binfw, 4))
cat("\n--- Bootstrap sup ---\n")
print(round(bsupw, 4))
cat("\n--- Bootstrap q1 ---\n")
print(round(bq1w, 4))
cat("\n--- Bootstrap q3 ---\n")
print(round(bq3w, 4))
                        
