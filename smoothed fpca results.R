##################################################################
##################################################################
#Regression analysis on the the smoothed FPCA
##################################################################
##################################################################
load("/Users/Selene/Desktop/menu+rfh/complete consecutive profiles.rdata")
load("/Users/Selene/Desktop/menu+rfh/J.rdata")
load("/Users/Selene/Desktop/menu+rfh/smooth fpca/Gb.rdata")
load("/Users/Selene/Desktop/menu+rfh/smooth fpca/Gw.rdata")
setwd("/Users/Selene/Desktop/menu+rfh/smooth fpca")

n=nrow(gooddf)/720
x=gooddf$activity
y=matrix(x,nrow=n,byrow=TRUE)
t=1:720
mu=apply(y, 2,mean)
resd=matrix(0, nrow=n, ncol=720) 
resd=t( t(y) - mu ) 
gooddf$identifier=as.character(gooddf$identifier)
M=length(unique(gooddf$identifier))
Jsum=J*(J-1)
SUM=sum(J*(J-1))
N=720

#eigen-decomposition of the covariance matrices
e1 <- eigen(Gb)
e2 <- eigen(Gw)
fpca1.value <- e1$values 
fpca2.value <- e2$values 
fpca1.value <- ifelse(fpca1.value>=0, fpca1.value, 0)
fpca2.value <- ifelse(fpca2.value>=0, fpca2.value, 0)
percent1 <- (fpca1.value)/sum(fpca1.value)
percent2 <- (fpca2.value)/sum(fpca2.value)
K1 <- max( which(cumsum(percent1) < 0.9 | percent1 > 1/N ) + 1)
K2 <- max( which(cumsum(percent2) < 0.9 | percent2 > 1/N ) + 1)
rho=sum(fpca1.value)/(sum(fpca1.value)+sum(fpca2.value))
#K1=9 K2=12 rho=0.23906

fpca1.vectors <- e1$vectors[, 1:K1]
fpca2.vectors <- e2$vectors[, 1:K2]

for(i in 1:K1) {
  v2 <- fpca1.vectors[,i]
  tempsign <- sum(v2)
  fpca1.vectors[,i] <- ifelse(tempsign<0, -1,1) * v2
}
for(i in 1:K2) {
  v2 <- fpca2.vectors[,i]
  tempsign <- sum(v2)
  fpca2.vectors[,i] <- ifelse(tempsign<0, -1,1) * v2
}
save(fpca1.vectors,file='fpca1.vectors.rdata')
save(fpca2.vectors,file='fpca2.vectors.rdata')

#calculate principal component scores using the projection method introduced in Di's paper
cross.integral <- matrix(0, K1, K2)
for(i in 1:K1)
  for(j in 1:K2) 
    cross.integral[i,j] <- sum(fpca1.vectors[,i]* fpca2.vectors[,j]) 

n=nrow(gooddf)/720
int1 <- matrix(0, n, K1)
int2 <- matrix(0, n, K2)
for(i in 1:n)   {
  for(j in 1:K1)  {
    int1[ i ,j] <- sum( resd[i,] * fpca1.vectors[,j] ) 
  }
  for(j in 1:K2) {
    int2[ i ,j] <- sum( resd[i,] * fpca2.vectors[,j] )    
  }
}

s1 <- matrix(0, n, K1)
s2 <- matrix(0, n, K2)
library(MASS)
design.xi <- ginv( diag(rep(1,K1)) - cross.integral %*% t(cross.integral) )
for(m in 1:M) {
  resid <- rep(0, K1)
  for(j in 1:J[m]) {
    index <-  ifelse(m==1,0,sum(J[1:(m-1)])) + j
    resid <- resid + ( int1[index,] - drop(cross.integral %*% int2[index,]) )/J[m]
  }
  index.m <- ( ifelse(m==1,0,sum(J[1:(m-1)])) + 1 ) : (sum(J[1:m]))
  xi.temp <- design.xi %*% resid
  s1[index.m,] <- matrix(rep(xi.temp, each=J[m]), nrow=J[m])
  for(j in 1:J[m]) {
    index <- ifelse(m==1,0,sum(J[1:(m-1)])) + j
    s2[index,] <- int2[index,] - drop( t(cross.integral) %*% xi.temp )
  }
}
coef.1=s1
save(coef.1,file='coef.1.rdata')
coef.2=s2
save(coef.2,file='coef.2.rdata')


#menu vs rfh (no cancer vs cancer) analysis
n=nrow(gooddf2)/720
c.m=coef.1[1:n,]
nc.m=coef.1[(n+1):nrow(coef.1),]
c.m=unique(c.m)
nc.m=unique(nc.m)
t.test(c.m[,1],nc.m[,1]) #p=0.2018
t.test(c.m[,2],nc.m[,2]) #p=0.7728
t.test(c.m[,3],nc.m[,3]) #p=0.3929

c.l2=rep(0,nrow(c.m))
nc.l2=rep(0,nrow(nc.m))
for(i in 1:nrow(c.m)){
	v=c.m[i,]
	c.l2[i]=sqrt(sum(v^2))
}
for(i in 1:nrow(nc.m)){
	v=nc.m[i,]
	nc.l2[i]=sqrt(sum(v^2))
}
t.test(c.l2,nc.l2) #p=0.1667


#biomarkers and qol
load("/Users/Selene/Desktop/menu+rfh/df.rdata")
com=matrix(nrow=M,ncol=ncol(coef.1))
for(i in 1:M){
	row=ifelse(i==1,0,sum(J[1:(i-1)]))+1
	com[i,]=coef.1[row,]
}
com.l2=rep(0,M)
for(i in 1:M){
	v=com[i,]
	com.l2[i]=sqrt(sum(v^2))
}
com=com/1000
df$PC1=com[,1]
df$PC2=com[,2]
df$PC3=com[,3]
df$PC4=com[,4]
df$PC5=com[,5]
df$PC6=com[,6]
df$PC7=com[,7]
df$PC8=com[,8]
df$PC9=com[,9]
df$l2=com.l2/1000

mod1.2=lm(insulin~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9, data=df)
mod2.2=lm(CRP~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9, data=df)
mod3.2=lm(QOLp~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9, data=df)
mod4.2=lm(QOLm~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9, data=df)
mod5.2=lm(insulin~age+education+bmi+smoke+cancer+l2, data=df)
mod6.2=lm(CRP~age+education+bmi+smoke+cancer+l2, data=df)
mod7.2=lm(QOLp~age+education+bmi+smoke+cancer+l2, data=df)
mod8.2=lm(QOLm~age+education+bmi+smoke+cancer+l2, data=df)
#nothing is significant


#total variation of the derivative of the intensity curve
pred=fpca1.vectors %*% t(coef.1)
com=matrix(nrow=720,ncol=M)
for(i in 1:M){
	col=ifelse(i==1,0,sum(J[1:(i-1)]))+1
	com[,i]=pred[,col]
}
totv=rep(NA,M)
for (i in 1: M){
	z=com[,i]
	v1=z[1:719]
	v2=z[2:720]
	totv[i]=sum((v2-v1)^2)
}
df$totv=totv

mod5.2=lm(insulin~age+education+bmi+smoke+cancer+totv, data=df)
mod6.2=lm(CRP~age+education+bmi+smoke+cancer+totv, data=df)
mod7.2=lm(QOLp~age+education+bmi+smoke+cancer+totv, data=df)
mod8.2=lm(QOLm~age+education+bmi+smoke+cancer+totv, data=df)
#nothing is significant



























