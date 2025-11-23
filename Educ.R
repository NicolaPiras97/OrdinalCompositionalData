data("educFM")
father <- as.matrix(educFM[,2:4])
y <- father/rowSums(father)
mother <- as.matrix(educFM[,5:7])
x <- mother/rowSums(mother)
N=dim(y)[1]
Cy=dim(y)[2]
Cx=dim(x)[2]

ydata<-list()
xdata<-list()
for(i in 1:N){
  ydata[[i]]<-y[i,]
}
for(i in 1:N){
  xdata[[i]]<-x[i,]
}
ycdf<-matrix(0,nrow=N,ncol=Cy)
for(i in 1:N){
  ycdf[i,]<-cumsum(ydata[[i]])
}

weights<-c(median(ycdf[,2])-median(ycdf[,1]),median(ycdf[,3])-median(ycdf[,2]))

sol<-solve_simplex_lp(xdata,ydata,weights)
A<-sol$A
AC<-(codalm(y,x))

Ex<-rep(0,N)
for(i in 1:N){
  for(k in 1:Cx){
     Ex[i]=Ex[i]+k*x[i,k]
  }
}
Ey<-rep(0,N)
for(i in 1:N){
  for(k in 1:Cy){
    Ey[i]=Ey[i]+k*y[i,k]
  }
}
#OCC
cov(Ex,Ey)/(sd(Ex)*sd(Ey))*((N-1)/N)


ymean<-vector()
for(j in 1:Cy){
  ymean[j]=0
}

ymedianw<-vector()
for(j in 1:Cy){
  ymedianw[j]=0
}

 for(j in 1:Cy){
   for(i in 1:N){
     ymean[j]<-ymean[j]+ydata[[i]][j]
   }
 }
 ymean=ymean/N


distot<-rep(0,N)
for(i in 1:N){
  for(j in 1:N){
    if(i!=j){
    if(porwd(c(weights,1),y[i,],y[j,])==2){
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
  denw=denw+wd(c(weights,1),ymedianw,ydata[[i]])
}


denC<-0
for(i in 1:N){
    denC=denC+kld(y[i,],ymean)
}

indporwd=0
for(i in 1:(N-1)){
  for(j in (1+i):N){
    if(porwd(c(weights,1),y[i,],y[j,])==porwd(c(weights,1),A%*%x[i,],A%*%x[j,])){
      indporwd=indporwd+1
    }
  }
}
OPIwd=indporwd/choose(N,2)
OPIwd*((N-2)/N)

SSRw<-0
for(i in 1:N){
  SSRw=SSRw+wd(c(weights,1),ymedianw,A%*%xdata[[i]])
}
R2W<-SSRw/denw
R2W*((N-2)/N)

SSRc<-0
for(i in 1:N){
  SSRc=SSRc+kld(t(AC)%*%xdata[[i]],ymean)
}
R2C<-SSRc/denC
R2C*((N-2)/N)


#bootstrap
B=1000
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

ycdf<-matrix(0,nrow=N,ncol=Cy)
for(i in 1:N){
  ycdf[i,]<-cumsum(yB[i,])
}
weightsB[[b]]<-c(median(ycdf[,2])-median(ycdf[,1]),median(ycdf[,3])-median(ycdf[,2]))

ydataB<-list()
xdataB<-list()
for(i in 1:N){
  ydataB[[i]]<-yB[i,]
}
for(i in 1:N){
  xdataB[[i]]<-xB[i,]
}
solB<-solve_simplex_lp( xdataB , ydataB , weightsB[[b]] )
Btot[[b]]<-solB$A

for(j in 1:Cy){
  vartotc[j]=vartotc[j]+wd(c(weightsB[[b]],1),A[,j],Btot[[b]][,j])
}

}

#sd with wd
sqrt(vartotc/B)

#95% confidence region
distotB<-matrix(0,nrow=B,ncol=Cy)
for(i in 1:B){
  for(j in 1:B){
    if(i!=j){
    for(k in 1:Cy){
    if(porwd(c(weightsB[[b]],1),Btot[[i]][,k],Btot[[j]][,k])==2){
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

A
bmedianw
binfw
bsupw
bq1w
bq3w
