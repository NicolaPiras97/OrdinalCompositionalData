source(file="Functions_simplex.R")

#setting 1
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
b2b <- 0
b2c <- 0
b2d<- 0
x1<-matrix(0,nrow=Ni,ncol=Cx1)
x2<-matrix(0,nrow=Ni,ncol=Cx2)
y<-matrix(0,nrow=Ni,ncol=Cy)

#setting 2
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
b2b <- 0.5
b2c <- 3
b2d<- 5
x1<-matrix(0,nrow=Ni,ncol=Cx1)
x2<-matrix(0,nrow=Ni,ncol=Cx2)
y<-matrix(0,nrow=Ni,ncol=Cy)

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

ycdf<-matrix(0,nrow=dim(y)[1],ncol=dim(y)[2])
for(i in 1:dim(y)[1]){
  ycdf[i,]<-cumsum(y[i,])
}
weights<-c()
for(i in 1:(Cy-1)){
  weights[i]<-median(ycdf[,i+1])-median(ycdf[,i])
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
Z_data <-lapply(1:Ni, function(i) tensor_product_ordered(x1data[[i]],x2data[[i]])$product)
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
result<-matrix(0,nrow=iter,ncol=6)
colnames(result)<-c("R2W_X1X2","OPI_X1X2","Error validation-set_X1X2","R2W_X1","OPI_X1","Error validation-set_X1")
weightstot1<-matrix(0,nrow=iter,ncol=Cy-1)
weightstot2<-matrix(0,nrow=iter,ncol=Cy-1)

for(it in 1:iter){
  
  for(n in 1:N){
    
    x1[n,]<-rdirichlet(1,rep(1,Cx1))
    x2[n,]<-rdirichlet(1,rep(1,Cx2))
    y1[n,]<-mat%*%ordprod(x1[n,],x2[n,])
    y1[n,]<-rdirichlet(1,10*y1[n,])
    
    y2[n,]<-mat2%*%x1[n,]
    y2[n,]<-rdirichlet(1,10*y2[n,])
    
  }
  
  ycdf<-matrix(0,nrow=dim(y1)[1],ncol=dim(y1)[2])
  for(i in 1:dim(y1)[1]){
    ycdf[i,]<-cumsum(y1[i,])
  }
  weightsa<-c()
  for(i in 1:(Cy-1)){
    weightsa[i]<-median(ycdf[,i+1])-median(ycdf[,i])
  }
  ycdf<-matrix(0,nrow=dim(y2)[1],ncol=dim(y2)[2])
  for(i in 1:dim(y2)[1]){
    ycdf[i,]<-cumsum(y2[i,])
  }
  weightsb<-c()
  for(i in 1:(Cy-1)){
    weightsb[i]<-median(ycdf[,i+1])-median(ycdf[,i])
  }
  
  
  tr<-trunc(0.7*N)
  va<-N-tr
  intot<-sample(1:N,N,replace=F)
  intr<-intot[1:tr]
  inva<-intot[tr+1:va]
  
  ycdf<-matrix(0,nrow=tr,ncol=dim(y1)[2])
  for(i in 1:tr){
    ycdf[i,]<-cumsum(y1[intr[i],])
  }
  weights1<-c()
  for(i in 1:(Cy-1)){
    weights1[i]<-median(ycdf[,i+1])-median(ycdf[,i])
  }
  ycdf<-matrix(0,nrow=tr,ncol=dim(y2)[2])
  for(i in 1:tr){
    ycdf[i,]<-cumsum(y2[intr[i],])
  }
  weights2<-c()
  for(i in 1:(Cy-1)){
    weights2[i]<-median(ycdf[,i+1])-median(ycdf[,i])
  }
  
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
  
  Z_data <-lapply(1:tr, function(i) tensor_product_ordered(Plistx1[[i]],Plistx2[[i]])$product)
  resestim <- solve_simplex_lp( Z_data , Pprimey1 , weights1 )
  matrixcoeff1<-resestim$A
  
  ymedian<-vector()
  for(j in 1:Cy){
    ymedian[j]=0
  }

  distot<-rep(0,N)
  for(i in 1:N){
    for(j in 1:N){
      if(i!=j){
        if(porwd(c(weightsa,1),Pprime1[[i]],Pprime1[[j]])==2){
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
    ymedian<-Pprime1[[ordtot[round(N/2)]]]
  }
  if(N%%2==0){
    ymedian<-apply(cbind(Pprime1[[ordtot[N/2]]],Pprime1[[ordtot[N/2+1]]]),1,mean)
  }

  denw<-0
  for(i in 1:N){
    denw=denw+wd(c(weightsa,1),ymedian,Pprime1[[i]])
  }

  SSRw<-0
  for(i in 1:N){
    SSRw=SSRw+wd(c(weightsa,1),ymedian,matrixcoeff1%*%tensor_product_ordered(Plist1[[i]],Plist2[[i]])$product)
    }
  R2W<-SSRw/denw

  indporwd=0
  for(i in 1:(N-1)){
    for(j in (1+i):N){
      if(porwd(c(weightsa,1),Pprime1[[i]],Pprime1[[j]])==porwd(c(weightsa,1),matrixcoeff1%*%tensor_product_ordered(Plist1[[i]],Plist2[[i]])$product,matrixcoeff1%*%tensor_product_ordered(Plist1[[j]],Plist2[[j]])$product)){
          indporwd=indporwd+1
      }
    }
  }
  PORwd=indporwd/choose(N,2)

  result[it,1]<-R2W*((N-2*(Cy-1))/N)
  result[it,2]<-PORwd*((N-2*(Cy-1))/N)

  yestim1<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    yestim1[i,]<-matrixcoeff1%*%tensor_product_ordered(Plist1[[inva[i]]],Plist2[[inva[i]]])$product
  }
  yerr1<-rep(0,va)
  for(i in 1:va){
    yerr1[i]<-wd(c(weights1,1),yestim1[i,],Pprime1[[inva[i]]])
  }
  
  result[it,3]<-sum(yerr1)/va
  
  
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
  
  ymedian<-vector()
  for(j in 1:Cy){
    ymedian[j]=0
  }

  distot<-rep(0,N)
  for(i in 1:N){
    for(j in 1:N){
      if(i!=j){
        if(porwd(c(weightsb,1),y2[i,],y2[j,])==2){
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
    ymedian<-y2[ordtot[round(N/2)],]
  }
  if(N%%2==0){
    ymedian<-apply(cbind(y2[ordtot[N/2],],y2[ordtot[N/2+1],]),1,mean)
  }
  denw<-0
  for(i in 1:N){
    denw=denw+wd(c(weightsb,1),ymedian,Pprime2[[i]])
  }

  SSRw<-0
  for(i in 1:N){
    SSRw=SSRw+wd(c(weightsb,1),ymedian,matrixcoeff2%*%Plist1[[i]])
  }
  R2W<-SSRw/denw

  indporwd=0
  for(i in 1:(N-1)){
    for(j in (1+i):N){
      if(porwd(c(weightsb,1),y2[i,],y2[j,])==porwd(c(weightsb,1),matrixcoeff2%*%Plist1[[i]],matrixcoeff2%*%Plist1[[j]])){
        indporwd=indporwd+1
      }
    }
  }
  PORwd=indporwd/choose(N,2)

  result[it,4]<-R2W*((N-(Cy-1))/N)
  result[it,5]<-PORwd*((N-(Cy-1))/N)

  yestim3<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    yestim3[i,]<-matrixcoeff2%*%as.matrix(Plist1[[inva[i]]])
  }
  yerr3<-rep(0,va)
  for(i in 1:va){
    yerr3[i]<-wd(c(weights2,1),yestim3[i,],Pprime2[[inva[i]]])
  }
  
  result[it,6]<-sum(yerr3)/va
  
  weightstot2[it,]<-weights2
  
}

apply(result,2,mean)
apply(result,2,sd)
apply(weightstot1,2,mean)
apply(weightstot1,2,sd)
apply(weightstot2,2,mean)
apply(weightstot2,2,sd)


