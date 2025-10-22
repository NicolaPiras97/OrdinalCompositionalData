source(file="Functions_simplex.R")

#Scenario 1
Cx=3
Cy=3
data("educFM")
father <- as.matrix(educFM[,2:4])
y <- father/rowSums(father)
mother <- as.matrix(educFM[,5:7])
x <- mother/rowSums(mother)
Ac<-codalm(y,x)
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
weights<-c(mean(rid[,2]-rid[,1]),mean(rid[,3]-rid[,2]))

ydati<-list()
xdati<-list()
for(i in 1:dim(y)[1]){
  ydati[[i]]<-y[i,]
}
for(i in 1:dim(y)[1]){
  xdati[[i]]<-x[i,]
}
sol<-solve_simplex_lp( xdati , ydati , weights )
mat<-sol$A

#Scenario 2
m<-100
Ni<-50
generation<-"CUB"
pi<-0.8
Cx=5
Cy=5
x<-matrix(0,nrow=Ni,ncol=Cx)
y<-matrix(0,nrow=Ni,ncol=Cy)
ridit<-matrix(0,nrow=N,ncol=Cy)
xi2<-c(0.8,0.7,0.5,0.3,0.2)

#Scenario 3
m<-100
Ni<-50
generation<-"CUB"
pi<-0.8
Cx=3
Cy=3
xi2<-c(0.9,0.5,0.1)
x<-matrix(0,nrow=Ni,ncol=Cx)
y<-matrix(0,nrow=Ni,ncol=Cy)
ridit<-matrix(0,nrow=N,ncol=Cy)

#Scenario 4
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
x<-matrix(0,nrow=Ni,ncol=Cx)
y<-matrix(0,nrow=Ni,ncol=Cy)
ridit<-matrix(0,nrow=N,ncol=Cy)

#Scenario 5
Cx=3
Cy=5
b01 <- 1 
b02 <- 0.05
b03 <- -0.05
b04 <- -1
b1b <- 2
b1c <- 6
x<-matrix(0,nrow=Ni,ncol=Cx)
y<-matrix(0,nrow=Ni,ncol=Cy)
ridit<-matrix(0,nrow=N,ncol=Cy)

for(n in 1:Ni){
  
  if(generation=="CUB"){
    
    x[n,]<-rdirichlet(1,rep(1,Cx))
    
    rb<-rbinom(1,1,prob=c(0.5,0.5))
    rb<-1
    for(i in 1:Cx){
      if(which.max(x[n,])==i){
        if(rb==0){
          xi<-xi2[Cx+1-i]
        }
        if(rb==1){
          xi<-xi2[i]
        }
      }
    }  
    
    classiy<-simcub(m,Cy,pi,xi)
    y[n,]<-table(factor(classiy,levels=seq(1:Cy)))/m
    
  }
  
  if(generation=="Proportional_odds"){
    
    if(Cx==3){
      x2<-c()
      x2 <- sample(c("A","B","C"), m, replace=TRUE, prob=rdirichlet(1,rep(1,Cx)))
      x[n,]<-table(factor(x2,levels=c("A","B","C")))/m
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
      x[n,]<-table(factor(x2,levels=c("A","B","C","D","E")))/m
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
  
  ridit[n,1]<-0.5*y[n,1]
  for(c in 2:Cy){
    ridit[n,c]<-0.5*y[n,c]
    for(i in 1:(c-1)){
      ridit[n,c]<-ridit[n,c]+y[n,c-i]
    }
  }
}

Ac<-codalm(y,x)
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
#weights<-rep(1,Cy-1) #for case A-1

ydati<-list()
xdati<-list()
for(i in 1:dim(y)[1]){
  ydati[[i]]<-y[i,]
}
for(i in 1:dim(y)[1]){
  xdati[[i]]<-x[i,]
}
sol<-solve_simplex_lp( xdati , ydati , weights )
mat<-sol$A

#Simulation
iter=100
N=50 #chainge for different sample sizes
x<-matrix(0,nrow=N,ncol=Cx)
y1<-matrix(0,nrow=N,ncol=Cy)
y2<-matrix(0,nrow=N,ncol=Cy)
ridit1<-matrix(0,nrow=N,ncol=Cy)
ridit2<-matrix(0,nrow=N,ncol=Cy)
result<-matrix(0,nrow=iter,ncol=11)
colnames(result)<-c("OCC","R2 (A)","R2 (B)","Error (A A) validation-set","Error (B A) validation-set","Error (A A) training-set","Error (B A) training-set","Error (A B) validation-set","Error (B B) validation-set","Error (A B) training-set","Error (B B) training-set")
weightstot1<-matrix(0,nrow=iter,ncol=Cy-1)
weightstot2<-matrix(0,nrow=iter,ncol=Cy-1)

for(it in 1:iter){
  
  for(n in 1:N){
      
      x[n,]<-rdirichlet(1,rep(1,Cx))
      y1[n,]<-mat%*%x[n,]
      y2[n,]<-t(Ac)%*%x[n,]
      
      y2[n,]<-rdirichlet(1,10*y2[n,])
      y1[n,]<-rdirichlet(1,10*y1[n,])
      
    ridit1[n,1]<-0.5*y1[n,1]
    for(c in 2:Cy){
      ridit1[n,c]<-0.5*y1[n,c]
      for(i in 1:(c-1)){
        ridit1[n,c]<-ridit1[n,c]+y1[n,c-i]
      }
    }
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
  #weights1<-rep(1,Cy-1) #for case A-1
  #weights2<-rep(1,Cy-1) #for case A-1
  
  tr<-trunc(0.7*N)
  va<-N-tr
  intot<-sample(1:N,N,replace=F)
  intr<-intot[1:tr]
  inva<-intot[tr+1:va]
  
  Plista <- list()
  for(i in 1:N){
    Plista[[i]]<-x[i,]
  }
  Pprimo1<-list()
  for(i in 1:N){
    Pprimo1[[i]]<-y1[i,]
  }
  Pprimo2<-list()
  for(i in 1:N){
    Pprimo2[[i]]<-y2[i,]
  }
  
  Plistax <- list()
  for(i in 1:tr){
    Plistax[[i]]<-Plista[[intr[i]]]
  }
  Pprimoy1<-list()
  for(i in 1:tr){
    Pprimoy1[[i]]<-Pprimo1[[intr[i]]]
  }
  Pprimoy2<-list()
  for(i in 1:tr){
    Pprimoy2[[i]]<-Pprimo2[[intr[i]]]
  }
  
  resstima <- solve_simplex_lp( Plistax , Pprimoy1 , weights1 )
  matricecoeff1<-resstima$A
  resstima <- solve_simplex_lp( Plistax , Pprimoy2 , weights2 )
  matricecoeff2<-resstima$A
  
  Px<-matrix(0,nrow=N-va,ncol=Cx)
  for(i in 1:(N-va)){
    Px[i,]<-Plistax[[i]]
  }
  Py1<-matrix(0,nrow=N-va,ncol=Cy)
  for(i in 1:(N-va)){
    Py1[i,]<-Pprimoy1[[i]]
  }
  Py2<-matrix(0,nrow=N-va,ncol=Cy)
  for(i in 1:(N-va)){
    Py2[i,]<-Pprimoy2[[i]]
  }
  
  Ex<-rep(0,tr)
  for(i in 1:tr){
    for(k in 1:Cx){
      Ex[i]=Ex[i]+k*x[i,k]
    }
  }
  Ey<-rep(0,tr)
  for(i in 1:tr){
    for(k in 1:Cy){
      Ey[i]=Ey[i]+k*y1[i,k]
    }
  }
  
  ymedio<-vector()
  for(j in 1:Cy){
    ymedio[j]=0
  }
  for(j in 1:Cy){
    for(i in 1:tr){
      ymedio[j]<-ymedio[j]+Pprimoy1[[i]][j]
    }
  }
  ymedio=ymedio/tr
  den<-0
  for(i in 1:tr){
    den=den+wd(c(weights1,1),ymedio,Pprimoy1[[i]])
  }
  R2<-1-resstima$loss/den
  result[it,2]<-R2
  
  A1<-(codalm(Py1,Px))
  A2<-(codalm(Py2,Px))
  b<-matrix(0,nrow=N,ncol=dim(A2)[2])
  lo<-0
  for(i in 1:tr){
    for(k in 1:dim(A2)[2]){
      for(j in 1:dim(A2)[1]){
        b[i,k]=b[i,k]+A2[j,k]*Plistax[[i]][j]
      }
    }
    lo=lo+kld(Pprimoy2[[i]],b[i,])
  }
  
  ymedio<-vector()
  for(j in 1:Cy){
    ymedio[j]=0
  }
  for(j in 1:Cy){
    for(i in 1:tr){
      ymedio[j]<-ymedio[j]+Pprimoy2[[i]][j]
    }
  }
  ymedio=ymedio/tr
  den<-0
  for(i in 1:tr){
    den<-den+kld(Pprimoy2[[i]],ymedio)
  }
  
  R2<-1-lo/den
  result[it,3]<-R2
  
  result[it,1]<-cov(Ex,Ey)/(sd(Ex)*sd(Ey))*((N-va-1)/(N-va))
  
  ystime1<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    ystime1[i,]<-matricecoeff1%*%as.matrix(Plista[[inva[i]]])
  }
  yerr1<-rep(0,va)
  for(i in 1:va){
    yerr1[i]<-wd(c(weights1,1),ystime1[i,],Pprimo1[[inva[i]]])
  }
  ystime2<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    ystime2[i,]<-matricecoeff2%*%as.matrix(Plista[[inva[i]]])
  }
  yerr2<-rep(0,va)
  for(i in 1:va){
    yerr2[i]<-wd(c(weights2,1),ystime2[i,],Pprimo2[[inva[i]]])
  }
  
  ystime3<-matrix(0,nrow=tr,ncol=Cy)
  for(i in 1:tr){
    ystime3[i,]<-matricecoeff1%*%as.matrix(Plista[[intr[i]]])
  }
  yerr3<-rep(0,tr)
  for(i in 1:tr){
    yerr3[i]<-wd(c(weights1,1),ystime3[i,],Pprimo1[[intr[i]]])
  }
  ystime4<-matrix(0,nrow=tr,ncol=Cy)
  for(i in 1:tr){
    ystime4[i,]<-matricecoeff2%*%as.matrix(Plista[[intr[i]]])
  }
  yerr4<-rep(0,tr)
  for(i in 1:tr){
    yerr4[i]<-wd(c(weights2,1),ystime4[i,],Pprimo2[[intr[i]]])
  }
  
  result[it,4]<-sum(yerr1)/va
  result[it,5]<-sum(yerr2)/va
  result[it,6]<-sum(yerr3)/tr
  result[it,7]<-sum(yerr4)/tr
  
  
  ystime5<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    ystime5[i,]<-t(A1)%*%as.matrix(Plista[[inva[i]]])
  }
  yerr5<-rep(0,va)
  for(i in 1:va){
    yerr5[i]<-kld(Pprimo1[[inva[i]]],ystime5[i,])
  }
  ystime6<-matrix(0,nrow=va,ncol=Cy)
  for(i in 1:va){
    ystime6[i,]<-t(A2)%*%as.matrix(Plista[[inva[i]]])
  }
  yerr6<-rep(0,va)
  for(i in 1:va){
    yerr6[i]<-kld(Pprimo2[[inva[i]]],ystime6[i,])
  }
  
  ystime7<-matrix(0,nrow=tr,ncol=Cy)
  for(i in 1:tr){
    ystime7[i,]<-t(A1)%*%as.matrix(Plista[[intr[i]]])
  }
  yerr7<-rep(0,tr)
  for(i in 1:tr){
    yerr7[i]<-kld(Pprimo1[[intr[i]]],ystime7[i,])
  }
  ystime8<-matrix(0,nrow=tr,ncol=Cy)
  for(i in 1:tr){
    ystime8[i,]<-t(A2)%*%as.matrix(Plista[[intr[i]]])
  }
  yerr8<-rep(0,tr)
  for(i in 1:tr){
    yerr8[i]<-kld(Pprimo2[[intr[i]]],ystime8[i,])
  }
  
  result[it,8]<-sum(yerr5)/va
  result[it,9]<-sum(yerr6)/va
  result[it,10]<-sum(yerr7)/tr
  result[it,11]<-sum(yerr8)/tr
  
  weightstot1[it,]<-weights1
  weightstot2[it,]<-weights2
  
}

apply(result,2,mean)
apply(result,2,sd)

apply(weightstot1,2,mean)
apply(weightstot1,2,sd)
apply(weightstot2,2,mean)
apply(weightstot2,2,sd)

