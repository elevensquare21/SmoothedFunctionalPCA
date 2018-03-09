##################################################################
##################################################################
#Add smoothness when estimating intensity curve by:
#1. smooth original accelerometer data
#2. smooth covariance function matrix
##################################################################
##################################################################
load("/Users/Selene/Desktop/menu+rfh/complete consecutive profiles.rdata")

##################################################################
#smooth original accelerometer data
##################################################################

#plot different smoothing methods
x=gooddf[1:720,3]
plot(x,type='l')
plot(lowess(x,f=0.01)$y,type='l')
plot(ksmooth(1:720,x,kernel='normal',bandwidth=100)$y, type='l')

#start with lowess with f=0.01
n=nrow(gooddf)/720
for(i in 1:n){
	x=gooddf[(i-1)*720+1:720,3]
	gooddf[(i-1)*720+1:720,3]=lowess(x,f=0.01)$y
}
save(gooddf, file='/Users/Selene/Desktop/menu+rfh/lowess gooddf.rdata')
#perform multilevel FPCA as usual
save(Gb,file='lowess Gb.rdata')
save(G,file='lowess G.rdata')
save(Gw,file='lowess Gw.rdata')
#K1=33 K2=91 rho=0.11536

#lowess with f=0.05
n=nrow(gooddf)/720
for(i in 1:n){
	x=gooddf[(i-1)*720+1:720,3]
	gooddf[(i-1)*720+1:720,3]=lowess(x,f=0.05)$y
}
save(gooddf, file='/Users/Selene/Desktop/menu+rfh/lowess2 gooddf.rdata')
#perform multilevel FPCA as usual
save(Gb,file='lowess2 Gb.rdata')
save(G,file='lowess2 G.rdata')
save(Gw,file='lowess2 Gw.rdata')
#K1=16 K2=34 rho=0.15425

#kernel smooth with bandwidth=10
n=nrow(gooddf)/720
for(i in 1:n){
	x=gooddf[(i-1)*720+1:720,3]
	gooddf[(i-1)*720 + 1:720, 3] = ksmooth(1:720, x, kernel = 'normal', bandwidth =10)$y
}
save(gooddf, file='/Users/Selene/Desktop/menu+rfh/ksmooth gooddf.rdata')
#perform multilevel FPCA as usual
save(Gb,file='ksmooth Gb.rdata')
save(G,file='ksmooth G.rdata')
save(Gw,file='ksmooth Gw.rdata')
#K1=27 K2=66 rho=0.13369

#kernel smooth with bandwidth=20
n=nrow(gooddf)/720
for(i in 1:n){
	x=gooddf[(i-1)*720+1:720,3]
	gooddf[(i-1)*720 + 1:720, 3] = ksmooth(1:720, x, kernel = 'normal', bandwidth =20)$y
}
save(gooddf, file='/Users/Selene/Desktop/menu+rfh/ksmooth2 gooddf.rdata')
#perform multilevel FPCA as usual
save(Gb,file='ksmooth2 Gb.rdata')
save(G,file='ksmooth2 G.rdata')
save(Gw,file='ksmooth2 Gw.rdata')
#K1=19 K2=42 rho=0.15179

#kernel smooth with bandwidth=100
n=nrow(gooddf)/720
for(i in 1:n){
	x=gooddf[(i-1)*720+1:720,3]
	gooddf[(i-1)*720 + 1:720, 3] = ksmooth(1:720, x, kernel = 'normal', bandwidth =100)$y
}
save(gooddf, file='/Users/Selene/Desktop/menu+rfh/ksmooth3 gooddf.rdata')
#perform multilevel FPCA as usual
save(Gb,file='ksmooth3 Gb.rdata')
save(G,file='ksmooth3 G.rdata')
save(Gw,file='ksmooth3 Gw.rdata')
#K1=9 K2=14 rho=0.24597


##################################################################
#smooth covariance function matrix
##################################################################
smoothing=TRUE
#try Chong-zhi's method: smooth covariance function
if(smoothing==TRUE) {

    ### drop the diagonal elements in G2 in smoothing step (the reason of doing so is
    ### explained in the paper)
        gw.temp <- Gw
        diag(gw.temp) <- rep(NA, N)
       
        data.gb <- data.frame(gb=as.vector(Gb),gw=as.vector(gw.temp),
                    x1=rep(seq(0,1,length=N), N), x2=rep(seq(0,1,length=N), each=N))
        attach(data.gb)
        
    ### load the "SemiPar" package which contains functions for semiparametric smoothing 
    ### using linear mixed models
        library(SemiPar)
        
        ####??      Define my own knots, but there seems to be a bug
             #myknots <- list( myknots=data.frame(x1=rep(seq(0,1,length=15),each=15), x2=rep(seq(0,1,length=15),15)) )
             #attach(myknots)
             #myknots <- data.frame(x1=rep(seq(0,1,length=15),each=15), x2=rep(seq(0,1,length=15),15)) 
             
             #takes too long with self defined knots, let the package do it on its own
             fit1<- spm(gb ~ f(x1, x2))
             fit<- spm(gw ~ f(x1, x2),omit.missing=T)        
        
        pred1 <- predict(fit1,newdata)
        newdata <- data.frame(x1=x1,x2=x2)
        pred <- predict(fit,newdata)

    
    ### obtain the smoothed covariance functions Gw(s,t) and Gb(s,t)
        var.noise <- mean( diag(Gw) - diag(matrix(pred,720)) )
        s.gw <-matrix(pred,720) 
        Gw <- (s.gw + t(s.gw) )/2 
        s.gb <- matrix(pred1, 720)
        Gb <- (s.gb + t(s.gb))/2
}

#K1=9 K2=12 rho=0.23906
setwd("/Users/Selene/Desktop/menu+rfh/smooth fpca")
save(Gb,file='Gb.rdata')
save(Gw,file='Gw.rdata')

#plot smoothed covariance function matrices
library('plotly')
p <- plot_ly(z=Gw, type="surface",showscale=TRUE)
p
p <- plot_ly(z=Gb, type="surface",showscale=TRUE)
p

#plot the first four level 1 eigen functions
par(mfrow=c(2,2))

mu.sm=ksmooth(1:720,mu,kernel='normal',bandwidth=100)
pc1=fpca1.vectors[,1]*1000
plus1=mu.sm$y+pc1$y
minus1=mu.sm$y-pc1$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 1 (30.5%)',ylim=c(200,450))
lines(plus1,lty=2,col='red')
lines(minus1,lty=2,col='blue')

pc2=fpca1.vectors[,2]*1000
plus2=mu.sm$y+pc2$y
minus2=mu.sm$y-pc2$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 2 (10.2%)',ylim=c(200,450))
lines(plus2,lty=2,col='red')
lines(minus2,lty=2,col='blue')

pc3=fpca1.vectors[,3]*1000
plus3=mu.sm$y+pc3$y
minus3=mu.sm$y-pc3$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 3 (7.6%)',ylim=c(200,450))
lines(plus3,lty=2,col='red')
lines(minus3,lty=2,col='blue')

pc4=fpca1.vectors[,4]*1000
plus4=mu.sm$y+pc4$y
minus4=mu.sm$y-pc4$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 4 (4.7%)',ylim=c(200,450))
lines(plus4,lty=2,col='red')
lines(minus4,lty=2,col='blue')














