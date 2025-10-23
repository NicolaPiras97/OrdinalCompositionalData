source(file="Functions_simplex.R")

m<-100
Ni<-50

b01 <- 1 
b02 <- 0.05
b03 <- -0.05
b04 <- -1
Cx1=3
Cx2=4
Cy=5
b1b <- 2
b1c <- 6
b2b <- 0.5 #0 for Setting 1
b2c <- 3 #0 for Setting 1
b2d<- 5 #0 for Setting 1
x1<-matrix(0,nrow=Ni,ncol=Cx1)
x2<-matrix(0,nrow=Ni,ncol=Cx2)
y<-matrix(0,nrow=Ni,ncol=Cy)
ridit<-matrix(0,nrow=N,ncol=Cy)

for(n in 1:Ni){
    
      xa<-c()
      xa <- sample(c("A","B","C"), m, replace=TRUE, prob=rdirichlet(1,rep(1,Cx1)))
      x1[n,]<-table(factor(xa,levels=c("A","B","C")))/m
      xb<-c()
      xb <- sample(c("A","B","C","D"), m, replace=TRUE, prob=rdirichlet(1,rep(1,Cx2)))
      x2[n,]<-table(factor(xb,levels=c("A","B","C","D")))/m
      
      logodds1 <- b01 + b1b * (xa=="B") + b1c * (xa=="C") + b2b * (xb=="B") + b2c * (xb=="C") + b2d * (xb=="D")  
      logodds2 <- b02 + b1b * (xa=="B") + b1c * (xa=="C") + b2b * (xb=="B") + b2c * (xb=="C") + b2d * (xb=="D")
      logodds3 <- b03 + b1b * (xa=="B") + b1c * (xa=="C") + b2b * (xb=="B") + b2c * (xb=="C") + b2d * (xb=="D")
      logodds4 <- b04 + b1b * (xa=="B") + b1c * (xa=="C") + b2b * (xb=="B") + b2c * (xb=="C") + b2d * (xb=="D")
      
      prob_2to5 <- inv_logit(logodds1)
      prob_3to5 <- inv_logit(logodds2)
      prob_4to5 <- inv_logit(logodds3)
      prob_5 <- inv_logit(logodds4)
      
      prob_1 <- 1 - prob_2to5
      prob_2 <- prob_2to5 - prob_3to5
      prob_3 <- prob_3to5 - prob_4to5
      prob_4 <- prob_4to5 - prob_5
      
      y1 <- c()
      for (i in 1:m) {
        y1[i] <- sample(1:Cy, 1, prob = c(prob_1[i], prob_2[i], prob_3[i], prob_4[i], prob_5[i]) )
      }
      y[n,]<-table(factor(y1,levels=seq(1:Cy)))/m
}

rid<-matrix(0,nrow=dim(y)[1],ncol=dim(y)[2])
for(i in 1:dim(y)[1]){
  rid[i,1]<-0.5*y[i,1]
  for(c in 2:dim(y)[2]){
    rid[i,c]<-0.5*y[i,c]
    for(j in 1:(c-1)){
      rid[i,c]<-rid[i,c]+y[i,c-j]
    }
  }
}
weights<-c()
for(i in 1:(Cy-1)){
  weights[i]<-mean(rid[,i+1]-rid[,i])
}

ydata<-list()
x1data<-list()
x2data<-list()
for(i in 1:dim(y)[1]){
  ydata[[i]]<-y[i,]
}
for(i in 1:dim(y)[1]){
  x1data[[i]]<-x1[i,]
}
for(i in 1:dim(y)[1]){
  x2data[[i]]<-x2[i,]
}
Z_data <-lapply(1:Ni, function(i) as.vector(x2data[[i]] %o%x1data[[i]] ))
sol<-solve_simplex_lp( Z_data , ydata , weights )
mat<-sol$A

sol2<-solve_simplex_lp(x1data , ydata , weights )
mat2<-sol2$A

#Simulation
iter=100
N=50 #chainge for different sizes
x1<-matrix(0,nrow=N,ncol=Cx1)
x2<-matrix(0,nrow=N,ncol=Cx2)
y1<-matrix(0,nrow=N,ncol=Cy)
y2<-matrix(0,nrow=N,ncol=Cy)
ridit1<-matrix(0,nrow=N,ncol=Cy)
ridit2<-matrix(0,nrow=N,ncol=Cy)
result<-matrix(0,nrow=iter,ncol=7)
colnames(result)<-c("R2adj","Error validation-set","Error training-set","R2 1 var","Error validation-set 1var","Error training-set 1 var","R2ds")
weightstot1<-matrix(0,nrow=iter,ncol=Cy-1)
weightstot2<-matrix(0,nrow=iter,ncol=Cy-1)

for(it in 1:iter){
  
  for(n in 1:N){
    
    x1[n,]<-rdirichlet(1,rep(1,Cx1))
    x2[n,]<-rdirichlet(1,rep(1,Cx2))
    y1[n,]<-mat%*%(lapply(n, function(i) as.vector(x2[n,] %o%x1[n,] )))[[1]]
    y1[n,]<-rdirichlet(1,10*y1[n,])
    
    ridit1[n,1]<-0.5*y1[n,1]
    for(c in 2:Cy){
      ridit1[n,c]<-0.5*y1[n,c]
      for(i in 1:(c-1)){
        ridit1[n,c]<-ridit1[n,c]+y1[n,c-i]
      }
    }
    
    y2[n,]<-mat2%*%x1[n,]
    y2[n,]<-rdirichlet(1,10*y2[n,])
    
    ridit2[n,1]<-0.5*y2[n,1]
    for(c in 2:Cy){
      ridit2[n,c]<-0.5*y2[n,c]
      for(i in 1:(c-1)){
        ridit2[n,c]<-ridit2[n,c]+y2[n,c-i]
      }
    }
    
  }
  
  weights1<-c()
  for(i in 1:(Cy-1)){
    weights1[i]<-mean(ridit1[,i+1]-ridit1[,i])
  }
  weights2<-c()
  for(i in 1:(Cy-1)){
    weights2[i]<-mean(ridit2[,i+1]-ridit2[,i])
  }
  
  tr<-trunc(0.7*N)
  va<-N-tr
  intot<-sample(1:N,N,replace=F)
  intr<-intot[1:tr]
  inva<-intot[tr+1:va]
  
  Plist1 <- list()
  for(i in 1:N){
    Plist1[[i]]<-x1[i,]
  }
  Plist2 <- list()
  for(i in 1:N){
    Plist2[[i]]<-x2[i,]
  }
  Pprime1<-list()
  for(i in 1:N){
    Pprime1[[i]]<-y1[i,]
  }
  
  Plistx1 <- list()
  for(i in 1:tr){
    Plistx1[[i]]<-Plist1[[intr[i]]]
  }
  Plistx2 <- list()
  for(i in 1:tr){
    Plistx2[[i]]<-Plist2[[intr[i]]]
  }
  
  Pprimey1<-list()
  for(i in 1:tr){
    Pprimey1[[i]]<-Pprime1[[intr[i]]]
  }
  
  Z_data <-lapply(1:tr, function(i) as.vector(Plistx2[[i]] %o%Plistx1[[i]] ))
  resestim <- solve_simplex_lp( Z_data , Pprimey1 , weights1 )
  matrixcoeff1<-resestim$A
  
  ymean<-vector()
  for(j in 1:Cy){
    ymean[j]=0
  }
  for(j in 1:Cy){
    for(i in 1:tr){
      ymean[j]<-ymean[j]+Pprimey1[[i]][j]
    }
  }
  ymean=ymean/tr
  den<-0
  for(i in 1:tr){
    den=den+wd(c(weights1,1),ymean,Pprimey1[[i]])
  }
  R2<-1-((tr-1)/(tr-2-1))*(resestim$loss/den)
  R2ds<-1-(((tr+1)/(tr))*((tr-2)/(tr-2-2))*((tr-1)/(tr-2-1)))*(resestim$loss/den) #ds
  result[it,1]<-R2
  result[it,7]<-R2ds
  
  yestim1<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    yestim1[i,]<-matrixcoeff1%*%(lapply(va, function(i) as.vector(Plist2[[inva[i]]] %o%Plist1[[inva[i]]] )))[[1]]
  }
  yerr1<-rep(0,va)
  for(i in 1:va){
    yerr1[i]<-wd(c(weights1,1),yestim1[i,],Pprime1[[inva[i]]])
  }
  yestim2<-matrix(0,nrow=tr,ncol=Cy)
  for(i in 1:tr){
    yestim2[i,]<-matrixcoeff1%*%(lapply(tr, function(i) as.vector(Plist2[[intr[i]]] %o%Plist1[[intr[i]]] )))[[1]]
  }
  yerr2<-rep(0,tr)
  for(i in 1:tr){
    yerr2[i]<-wd(c(weights1,1),yestim2[i,],Pprime1[[intr[i]]])
  }
  
  result[it,2]<-sum(yerr1)/va
  result[it,3]<-sum(yerr2)/tr
  
  weightstot1[it,]<-weights1
  
  Pprime2<-list()
  for(i in 1:N){
    Pprime2[[i]]<-y2[i,]
  }
  
  Pprimey2<-list()
  for(i in 1:tr){
    Pprimey2[[i]]<-Pprime2[[intr[i]]]
  }
  
  resestim2 <- solve_simplex_lp( Plistx1 , Pprimey2 , weights2 )
  matrixcoeff2<-resestim2$A
  
  ymean<-vector()
  for(j in 1:Cy){
    ymean[j]=0
  }
  for(j in 1:Cy){
    for(i in 1:tr){
      ymean[j]<-ymean[j]+Pprimey2[[i]][j]
    }
  }
  ymean=ymean/tr
  den<-0
  for(i in 1:tr){
    den=den+wd(c(weights2,1),ymean,Pprimey2[[i]])
  }
  R2<-1-resestim2$loss/den
  result[it,4]<-R2
  
  yestim3<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    yestim3[i,]<-matrixcoeff2%*%as.matrix(Plist1[[inva[i]]])
  }
  yerr3<-rep(0,va)
  for(i in 1:va){
    yerr3[i]<-wd(c(weights2,1),yestim3[i,],Pprime2[[inva[i]]])
  }
  yestim4<-matrix(0,nrow=tr,ncol=Cy)
  for(i in 1:tr){
    yestim4[i,]<-matrixcoeff2%*%as.matrix(Plist1[[intr[i]]])
  }
  yerr4<-rep(0,tr)
  for(i in 1:tr){
    yerr4[i]<-wd(c(weights2,1),yestim4[i,],Pprime2[[intr[i]]])
  }
  
  result[it,5]<-sum(yerr3)/va
  result[it,6]<-sum(yerr4)/tr
  
  weightstot2[it,]<-weights2
  
}

apply(result,2,mean)
apply(result,2,sd)
apply(weightstot1,2,mean)
apply(weightstot1,2,sd)
apply(weightstot2,2,mean)
apply(weightstot2,2,sd)
