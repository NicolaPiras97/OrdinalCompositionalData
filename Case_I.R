source(file="Functions_simplex.R")

#Scenario 1
Cx=3
Cy=3
data("educFM")
father <- as.matrix(educFM[,2:4])
y <- father/rowSums(father)
mother <- as.matrix(educFM[,5:7])
x0 <- mother/rowSums(mother)
Ac<-codalm(y,x0)
ycdf<-matrix(0,nrow=dim(y)[1],ncol=dim(y)[2])
for(i in 1:dim(y)[1]){
  ycdf[i,]<-cumsum(y[i,])
}
weights<-c(median(ycdf[,2])-median(ycdf[,1]),median(ycdf[,3])-median(ycdf[,2]))

ydata<-list()
xdata<-list()
for(i in 1:dim(y)[1]){
  ydata[[i]]<-y[i,]
}
for(i in 1:dim(y)[1]){
  xdata[[i]]<-x0[i,]
}
sol<-solve_simplex_lp( xdata , ydata , weights )
mat<-sol$A

#Scenario 2
m<-100
Ni<-50
generation<-"CUB"
pi<-0.8
Cx=5
Cy=5
x0<-matrix(0,nrow=Ni,ncol=Cx)
y<-matrix(0,nrow=Ni,ncol=Cy)
xi2<-c(0.8,0.7,0.5,0.3,0.2)

#Scenario 3
m<-100
Ni<-50
generation<-"CUB"
pi<-0.8
Cx=3
Cy=3
xi2<-c(0.9,0.5,0.1)
x0<-matrix(0,nrow=Ni,ncol=Cx)
y<-matrix(0,nrow=Ni,ncol=Cy)

#Scenario 4
m<-100
Ni<-50
generation<-"Proportional_odds"
Cx=5
Cy=5
b01 <- 1 
b02 <- 0.05
b03 <- -0.05
b04 <- -1
b1b <- 0.5
b1c <- 2
b1d <- 4
b1e <- 6
x0<-matrix(0,nrow=Ni,ncol=Cx)
y<-matrix(0,nrow=Ni,ncol=Cy)

#Scenario 5
m<-100
Ni<-50
generation<-"Proportional_odds"
Cx=3
Cy=5
b01 <- 1 
b02 <- 0.05
b03 <- -0.05
b04 <- -1
b1b <- 2
b1c <- 6
x0<-matrix(0,nrow=Ni,ncol=Cx)
y<-matrix(0,nrow=Ni,ncol=Cy)

for(n in 1:Ni){
  
  if(generation=="CUB"){
    
    x0[n,]<-rdirichlet(1,rep(1,Cx))
    
    for(i in 1:Cx){
      if(which.max(x0[n,])==i){
          xi<-xi2[i]
      }
    }  
    
    classiy<-simcub(m,Cy,pi,xi)
    y[n,]<-table(factor(classiy,levels=seq(1:Cy)))/m
    
  }
  
  if(generation=="Proportional_odds"){
    
    if(Cx==3){
      x2<-c()
      x2 <- sample(c("A","B","C"), m, replace=TRUE, prob=rdirichlet(1,rep(1,Cx)))
      x0[n,]<-table(factor(x2,levels=c("A","B","C")))/m
      logodds1 <- b01 + b1b * (x2=="B") + b1c * (x2=="C")  
      logodds2 <- b02 + b1b * (x2=="B") + b1c * (x2=="C") 
      logodds3 <- b03 + b1b * (x2=="B") + b1c * (x2=="C") 
      logodds4 <- b04 + b1b * (x2=="B") + b1c * (x2=="C") 
      
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
    
    if(Cx==5){
      x2<-c()
      x2 <- sample(c("A","B","C","D","E"), m, replace=TRUE, prob=rdirichlet(1,rep(1,Cx)))
      x0[n,]<-table(factor(x2,levels=c("A","B","C","D","E")))/m
      logodds1 <- b01 + b1b * (x2=="B") + b1c * (x2=="C") + b1d * (x2=="D") + b1e * (x2=="E") 
      logodds2 <- b02 + b1b * (x2=="B") + b1c * (x2=="C") + b1d * (x2=="D") + b1e * (x2=="E")
      logodds3 <- b03 + b1b * (x2=="B") + b1c * (x2=="C") + b1d * (x2=="D") + b1e * (x2=="E")
      logodds4 <- b04 + b1b * (x2=="B") + b1c * (x2=="C") + b1d * (x2=="D") + b1e * (x2=="E")
      
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
    
  }
  
}

Ac<-codalm(y,x0)
# weights<-rep(1,Cy-1) #for case A-1
ycdf<-matrix(0,nrow=dim(y)[1],ncol=dim(y)[2])
for(i in 1:dim(y)[1]){
  ycdf[i,]<-cumsum(y[i,])
}
weights<-c()
for(i in 1:(Cy-1)){
  weights[i]<-median(ycdf[,i+1])-median(ycdf[,i])
}

ydata<-list()
xdata<-list()
for(i in 1:dim(y)[1]){
  ydata[[i]]<-y[i,]
}
for(i in 1:dim(y)[1]){
  xdata[[i]]<-x0[i,]
}
sol<-solve_simplex_lp( xdata , ydata , weights )
mat<-sol$A

#Simulation
iter=1000
N=50 #chainge for different sample sizes
x<-matrix(0,nrow=N,ncol=Cx)
y1<-matrix(0,nrow=N,ncol=Cy)
y2<-matrix(0,nrow=N,ncol=Cy)
result<-matrix(0,nrow=iter,ncol=8)
colnames(result)<-c("OCC","R2W","R2Codalm","OPI","Error (A A) validation-set","Error (A B) validation-set","Error (B A) validation-set","Error (B B) validation-set")
weightstot1<-matrix(0,nrow=iter,ncol=Cy-1)
weightstot2<-matrix(0,nrow=iter,ncol=Cy-1)

for(it in 1:iter){
  
  for(n in 1:N){
      
      x[n,]<-rdirichlet(1,rep(1,Cx))
      y1[n,]<-mat%*%x[n,]
      y2[n,]<-t(Ac)%*%x[n,]
      
      y2[n,]<-rdirichlet(1,10*y2[n,])
      y1[n,]<-rdirichlet(1,10*y1[n,])
  }
  
  #weights1<-rep(1,Cy-1) #for case A-1
  #weights2<-rep(1,Cy-1) #for case A-1
  ycdf<-matrix(0,nrow=N,ncol=dim(y1)[2])
  for(i in 1:N){
    ycdf[i,]<-cumsum(y1[i,])
  }
  weights<-c()
  for(i in 1:(Cy-1)){
    weights[i]<-median(ycdf[,i+1])-median(ycdf[,i])
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
  
  Plist <- list()
  for(i in 1:N){
    Plist[[i]]<-x[i,]
  }
  Pprime1<-list()
  for(i in 1:N){
    Pprime1[[i]]<-y1[i,]
  }
  Pprime2<-list()
  for(i in 1:N){
    Pprime2[[i]]<-y2[i,]
  }
  
  Plistx <- list()
  for(i in 1:tr){
    Plistx[[i]]<-Plist[[intr[i]]]
  }
  Pprimey1<-list()
  for(i in 1:tr){
    Pprimey1[[i]]<-Pprime1[[intr[i]]]
  }
  Pprimey2<-list()
  for(i in 1:tr){
    Pprimey2[[i]]<-Pprime2[[intr[i]]]
  }
  
  resestim <- solve_simplex_lp( Plistx , Pprimey1 , weights1 )
  matrixcoeff1<-resestim$A
  resestim2 <- solve_simplex_lp( Plistx , Pprimey2 , weights2 )
  matrixcoeff2<-resestim2$A
  
  Px<-matrix(0,nrow=N-va,ncol=Cx)
  for(i in 1:(N-va)){
    Px[i,]<-Plistx[[i]]
  }
  Py1<-matrix(0,nrow=N-va,ncol=Cy)
  for(i in 1:(N-va)){
    Py1[i,]<-Pprimey1[[i]]
  }
  Py2<-matrix(0,nrow=N-va,ncol=Cy)
  for(i in 1:(N-va)){
    Py2[i,]<-Pprimey2[[i]]
  }
  
  Ex<-rep(0,N)
  for(i in 1:N){
    for(k in 1:Cx){
      Ex[i]=Ex[i]+k*x[i,k]
    }
  }
  Ey<-rep(0,N)
  for(i in 1:N){
    for(k in 1:Cy){
      Ey[i]=Ey[i]+k*y1[i,k]
    }
  }
  
  A1<-(codalm(Py1,Px))
  A2<-(codalm(Py2,Px))
  
  ymean<-vector()
  for(j in 1:Cy){
    ymean[j]=0
  }
  for(j in 1:Cy){
    for(i in 1:N){
      ymean[j]<-ymean[j]+Pprime2[[i]][j]
    }
  }
  ymean=ymean/N
  
  ymedian<-vector()
  for(j in 1:Cy){
    ymedian[j]=0
  }
  
  #median with Wasserstein order
  distot<-rep(0,N)
  for(i in 1:N){
    for(j in 1:N){
      if(i!=j){
        if(porwd(weights,y1[i,],y1[j,])==2){
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
    ymedian<-y1[ordtot[round(N/2)],]
  }
  if(N%%2==0){
    ymedian<-apply(cbind(y1[ordtot[N/2],],y1[ordtot[N/2+1],]),1,mean)
  }
  
  denw<-0
  for(i in 1:N){
    denw=denw+wd(weights,ymedian,Pprime1[[i]])
  }
  denC<-0
  for(i in 1:N){
    denC=denC+kld(y2[i,],ymean)
  }
  SSRw<-0
  for(i in 1:N){
    SSRw=SSRw+wd(weights,ymedian,matrixcoeff1%*%Plist[[i]])
  }
  R2W<-SSRw/denw
  SSRc<-0
  for(i in 1:N){
    SSRc=SSRc+kld(t(A2)%*%Plist[[i]],ymean)
  }
  R2C<-SSRc/denC
  
  indporwd=0
  for(i in 1:(N-1)){
    for(j in (1+i):N){
      if(porwd(weights,y1[i,],y1[j,])==porwd(weights,matrixcoeff1%*%x[i,],matrixcoeff1%*%x[j,])){
        indporwd=indporwd+1
      }
    }
  }
  OPI=indporwd/choose(N,2)
  
  result[it,1]<-cov(Ex,Ey)/(sd(Ex)*sd(Ey))*((N-1)/N)
  result[it,2]<-R2W*((N-(Cy-1))/N)
  result[it,3]<-R2C*((N-(Cy-1))/N)
  result[it,4]<-OPI*((N-(Cy-1))/N)
  
  yestim1<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    yestim1[i,]<-matrixcoeff1%*%as.matrix(Plist[[inva[i]]])
  }
  yerr1<-rep(0,va)
  for(i in 1:va){
    yerr1[i]<-wd(weights1,yestim1[i,],Pprime1[[inva[i]]])
  }
  yestim2<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    yestim2[i,]<-matrixcoeff2%*%as.matrix(Plist[[inva[i]]])
  }
  yerr2<-rep(0,va)
  for(i in 1:va){
    yerr2[i]<-wd(weights2,yestim2[i,],Pprime2[[inva[i]]])
  }
  
  result[it,5]<-sum(yerr1)/va
  result[it,7]<-sum(yerr2)/va
  
  yestim5<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    yestim5[i,]<-t(A1)%*%as.matrix(Plist[[inva[i]]])
  }
  yerr5<-rep(0,va)
  for(i in 1:va){
    yerr5[i]<-kld(Pprime1[[inva[i]]],yestim5[i,])
  }
  yestim6<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    yestim6[i,]<-t(A2)%*%as.matrix(Plist[[inva[i]]])
  }
  yerr6<-rep(0,va)
  for(i in 1:va){
    yerr6[i]<-kld(Pprime2[[inva[i]]],yestim6[i,])
  }
  
  result[it,6]<-sum(yerr5)/va
  result[it,8]<-sum(yerr6)/va
  
  weightstot1[it,]<-weights1
  weightstot2[it,]<-weights2
  
}

apply(result,2,mean)
apply(result,2,sd)

apply(weightstot1,2,mean)
apply(weightstot1,2,sd)
apply(weightstot2,2,mean)
apply(weightstot2,2,sd)
