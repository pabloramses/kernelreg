library(kernelreg)
library(locfit)
data<-data(ethanol)

#estimaciones
predict<- RegKer(sort(ethanol$E), ethanol$E, ethanol$NOx, Gauss, 0.05)
plot(ethanol$E, ethanol$NOx)
par(new=TRUE)
plot(sort(ethanol$E), predict,type = 'l', col = 'red')

predict2<- RegKer(sort(ethanol$E), ethanol$E, ethanol$NOx, Naive, 0.05)
plot(ethanol$E, ethanol$NOx)
par(new=TRUE)
plot(sort(ethanol$E), predict2,type = 'l', col = 'red')

#COMPARACIONES
###############################MDOELO 1########################################
k<-5; sigma<-0.5
estecm<-rep(0,k)
for (i in 1:k){n<- 10^i
               samplecm<-replicate(1000,{x1<-rnorm(n,0,1)
                                         y1<- exp(x1) + rnorm(n,0,sigma)
                                         pred<- RegKer(0, x1, y1, Gauss, n^(-1/5))
                                         (pred-exp(0))^2})
               estecm[i]<-mean(samplecm)}
##valores reales
f<-function(x){Gauss(x)}
f1<-function(x){-Gauss(x)*x}
mu<-function(x){exp(x)}
mu1<-function(x){exp(x)}
mu2<-function(x){exp(x)}


ecm_asintotico<-function(x,ker,n,h,sigma){((h^4)/4)*((mu2(x)+2*mu1(x)*f1(x)/f(x))^2)*(dk(ker)^2)+((sigma^2)*(ck(ker)^2))/(n*h*f(x))}
erreal<-rep(0,k)
for (i in 1:k){n<- 10^i
               h<-n^(-1/5)
               erreal[i]<-ecm_asintotico(0,'Gauss', n,h,sigma)}

plot(estecm, col= 'red', type='l',xlab='n=10^x',ylab='',main='Modelo 1 x=0')
lines(erreal,type='l')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)
###############################MDOELO 1 x=0.7########################################
k<-5; sigma<-0.5
estecm2<-rep(0,k)
for (i in 1:k){n<- 10^i
samplecm<-replicate(1000,{x1<-rnorm(n,0,1)
y1<- exp(x1) + rnorm(n,0,sigma)
pred<- RegKer(0.7, x1, y1, Gauss, n^(-1/5))
(pred-exp(0.7))^2})
estecm2[i]<-mean(samplecm)}
##valores reales
f<-function(x){Gauss(x)}
f1<-function(x){-Gauss(x)*x}
mu<-function(x){exp(x)}
mu1<-function(x){exp(x)}
mu2<-function(x){exp(x)}

erreal2<-rep(0,k)
for (i in 1:k){n<- 10^i
h<-n^(-1/5)
erreal2[i]<-ecm_asintotico(0.7,'Gauss', n,h,sigma)}
######PLOT
par(mfrow=c(1,2))
plot(estecm, col= 'red', type='l',xlab='n=10^x',ylab='',main='Modelo 1 x=0')
lines(erreal,type='l')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)

plot(estecm2, col= 'red', type='l',xlab='n=10^x',ylab='',main='Modelo 1 x=0.7')
lines(erreal2,type='l')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)
###############################MDOELO 2########################################
k<-5; sigma<-0.2;
estecm<-rep(0,k)
for (i in 1:k){n<- 10^i
               samplecm<-replicate(1000,{x1<-runif(n,-1,1)
                                         y1<- x1^2 + rnorm(n,0,sigma)
                                         pred<- RegKer(0, x1, y1, Gauss, n^(-1/5))
                                         (pred-0)^2})
               estecm[i]<-mean(samplecm)}
##valores reales
f<-function(x){1/2}
f1<-function(x){0}
mu<-function(x){x^2}
mu1<-function(x){2*x}
mu2<-function(x){2}


ecm_asintotico<-function(x,ker,n,h,sigma){((h^4)/4)*(mu2(x)+2*mu1(x)*f1(x)/f(x))*(dk(ker)^2)+((sigma^2)*(ck(ker)^2))/(n*h*f(x))}
erreal<-rep(0,k)
for (i in 1:k){n<- 10^i
h<-n^(-1/5)
erreal[i]<-ecm_asintotico(0,'Gauss', n,h)}



###############################MDOELO 2 x=0.7########################################
k<-5; sigma<-0.2;
estecm2<-rep(0,k)
for (i in 1:k){n<- 10^i
               samplecm<-replicate(1000,{x1<-runif(n,-1,1)
                                         y1<- x1^2 + rnorm(n,0,sigma)
                                         pred<- RegKer(0.7, x1, y1, Gauss, n^(-1/5))
                                         (pred-0.7^2)^2})
               estecm2[i]<-mean(samplecm)}
##valores reales
f<-function(x){1/2}
f1<-function(x){0}
mu<-function(x){x^2}
mu1<-function(x){2*x}
mu2<-function(x){2}



erreal2<-rep(0,k)
for (i in 1:k){n<- 10^i
               h<-n^(-1/5)
               erreal2[i]<-ecm_asintotico(0.7,'Gauss', n,h,sigma)}
##########PLOT
par(mfrow=c(1,2))
plot(estecm, col= 'red', type='l',xlab='n=10^x',ylab='',main='Modelo 2 x=0')
lines(erreal,type='l')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)

plot(estecm2, col= 'red', type='l',xlab='n=10^x',ylab='',main='Modelo 2 x=0.7')
lines(erreal2,type='l')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)
###############################MDOELO 3 x=1########################################
k<-5; sigma<-0.1;
estecm<-rep(0,k)
for (i in 1:k){n<- 10^i
               samplecm<-replicate(1000,{x1<-rexp(n,1)
                                         y1<- sin(x1) + rnorm(n,0,sigma)
                                         pred<- RegKer(1, x1, y1, Gauss, n^(-1/5))
                                         (pred-sin(1))^2})
               estecm[i]<-mean(samplecm)}
##valores reales
f<-function(x){exp(-x)}
f1<-function(x){-exp(-x)}
mu<-function(x){sin(x)}
mu1<-function(x){cos(x)}
mu2<-function(x){-sin(x)}



erreal<-rep(0,k)
for (i in 1:k){n<- 10^i
h<-n^(-1/5)
erreal[i]<-ecm_asintotico(0,'Gauss', n,h,sigma)}



###############################MDOELO 3 x=1.7########################################
k<-5; sigma<-0.1;
estecm2<-rep(0,k)
for (i in 1:k){n<- 10^i
               samplecm<-replicate(1000,{x1<-rexp(n,1)
                                         y1<- sin(x1) + rnorm(n,0,sigma)
                                         pred<- RegKer(1.7, x1, y1, Gauss, n^(-1/5))
                                         (pred-sin(1.7))^2})
               estecm2[i]<-mean(samplecm)}
##valores reales
erreal2<-rep(0,k)
for (i in 1:k){n<- 10^i
               h<-n^(-1/5)
               erreal2[i]<-ecm_asintotico(1.7,'Gauss', n,h,sigma)}
##########PLOT
par(mfrow=c(1,2))
plot(estecm, col= 'red', type='l',xlab='n=10^x',ylab='',main='Modelo 3 x=1')
lines(erreal,type='l')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)

plot(estecm2, col= 'red', type='l',xlab='n=10^x',ylab='',main='Modelo 3 x=1.7')
lines(erreal2,type='l')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)
###############################MDOELO 2 x=1########################################
k<-5; sigma<-0.2;
estecm2<-rep(0,k)
for (i in 1:k){n<- 10^i
               samplecm<-replicate(1000,{x1<-runif(n,-1,1)
                                         y1<- (x1^2)*sin(1/x1) + rnorm(n,0,sigma)
                                         pred<- RegKer(0.02, x1, y1, Gauss, n^(-1/5))
                                        (pred-(0.02^2)*sin(1/0.02))^2})
               estecm2[i]<-mean(samplecm)}


##valores reales
f<-function(x){1/2}
f1<-function(x){0}
mu<-function(x){(x^2)*sin(1/x)*(x!=0)}
mu1<-function(x){(2*x*sin(1/x)-cos(1/x))*(x!=0)}
mu2<-function(x){(2*sin(1/x)+(2/x)*cos(1/x)-sin(1/x)*(1/x^2))*(x!=0)}



erreal2<-rep(0,k)
for (i in 1:k){n<- 10^i
h<-n^(-1/5)
erreal2[i]<-ecm_asintotico(0.02,'Gauss', n,h,sigma)}
##########PLOT
par(mfrow=c(1,2))
plot(estecm, col= 'red', type='l',xlab='n=10^x',ylab='',main='Modelo 2 x=0')
lines(erreal,type='l')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)

plot(erreal2, col= 'black', type='l',xlab='n=10^x',ylab='',main='Modelo 2 x=0.7')
lines(estecm2,type='l',col='red')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)
erreal2
estecm2
erreal2

