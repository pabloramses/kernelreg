library(kernelreg)
k<-5
###############################MODELO 1########################################
####x=0####
sigma2<-0.5
estecm<-rep(0,k)
for (i in 1:k){n<- 10^i
               samplecm<-replicate(1000,{x1<-rnorm(n,0,1)
                                         y1<- exp(x1) + rnorm(n,0,sqrt(sigma))
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


estecm2<-rep(0,k)
for (i in 1:k){n<- 10^i
               samplecm<-replicate(1000,{x1<-rnorm(n,0,1)
                                         y1<- exp(x1) + rnorm(n,0,sigma)
                                         pred<- RegKer(0.7, x1, y1, Gauss, n^(-1/5))
                                         (pred-exp(0.7))^2})
               estecm2[i]<-mean(samplecm)}
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

###############################MDOELO 2 ########################################
####x=1####
sigma2<-0.1;
estecm<-rep(0,k)
for (i in 1:k){n<- 10^i
               samplecm<-replicate(1000,{x1<-rexp(n,1)
                                         y1<- sin(x1) + rnorm(n,0,sqrt(sigma))
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

####x=1.7####
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
plot(erreal, type='l',xlab='n=10^x',ylab='',main='Modelo 2 x=1')
lines(estecm,type='l',col= 'red')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)

plot(estecm2, col= 'red', type='l',xlab='n=10^x',ylab='',main='Modelo 2 x=1.7')
lines(erreal2,type='l')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)
###############################MDOELO 3 x=0.02########################################
sigma2<-0.2;
estecm<-rep(0,k)
for (i in 1:k){n<- 10^i
               samplecm<-replicate(1000,{x1<-runif(n,-1,1)
                                         y1<- (x1^2)*sin(1/x1) + rnorm(n,0,sqrt(sigma))
                                         pred<- RegKer(0.02, x1, y1, Gauss, n^(-1/5))
                                        (pred-(0.02^2)*sin(1/0.02))^2})
               estecm[i]<-mean(samplecm)}

##valores reales
f<-function(x){1/2}
f1<-function(x){0}
mu<-function(x){(x^2)*sin(1/x)*(x!=0)}
mu1<-function(x){(2*x*sin(1/x)-cos(1/x))*(x!=0)}
mu2<-function(x){(2*sin(1/x)+(2/x)*cos(1/x)-sin(1/x)*(1/x^2))*(x!=0)}

erreal<-rep(0,k)
for (i in 1:k){n<- 10^i
               h<-n^(-1/5)
               erreal[i]<-ecm_asintotico(0.02,'Gauss', n,h,sigma)}
##########PLOT
par(mfrow = c(1,1))
plot(erreal, col= 'black', type='l',xlab='n=10^x',ylab='',main='Modelo 3 x=0.02')
lines(estecm,type='l',col='red')
legend('topright',legend=c('ECM simulado', 'ECM asintótico'),col=c('red','black'),lty=1)

