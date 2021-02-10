library(kernelreg)
library(locfit)
data(ethanol)

x<-ethanol$E
y<-ethanol$NOx


n<-100
x<-runif(n,0,1)
y<-1-x+exp(-200*(x-1/2)^2) + rnorm(n,0,1)
varxi<-function(x){1+2*x}
wl<-function(x){(2<x)*(x<10)}


################PENALIZING FUNCTIONS############################################
G<-function(x,y,Ker,h){result<-0
                         for (i in 1:length(x)){
                         result<-result + ((y[i]-RegKer(x[i],x,y,Ker,h))^2)*varxi((Ker((x[i]-x[i])/h))/(sum(Ker((x[i]-x)/h))))*wl(x[i])}
                         return((1/n)*result)}

Gvec<-Vectorize(G, vectorize.args = 'h')



xord<-sort(x)
l<-10
hmin<-function(x){xord<-sort(x)
                  m<-sum(wl(xord))
                  h<-rep(0,m)
                  j<-1
                  for(i in which(wl(xord)==1)){h[j]<-min(xord[i]-xord[i-1],xord[i+1]-xord[i])
                                            j<-j+1}
                  return(h)}
hmean<-max(hmin(x))
h<-seq(hmean,8,0.01)
gh<-Gvec(x,y,Gauss,h)
hminG<-h[which(min(gh)==gh)]
plot(h,gh,type='l',xaxt='n',ylab='G(h)', main='Función penalizadora G(h)')
points(hmean,G(x,y,Gauss,hmean), col='red',pch=4)
points(hminG,G(x,y,Gauss,hminG), col='blue',pch=4)
axis(1,at=c(hmean,1,2,3,4,5,6,7,8),labels=c(round(hmean,1),1,2,3,4,5,6,7,8))
axis(1,at=hminG, labels=round(hminG,2))
legend('bottomright',legend=c('hmin', 'hG'), col=c('red','blue'), pch=4)
#################ONE-LEFT METHOD################################################
n<-100
x<-runif(n,0,1)
y<-1-x+exp(-200*(x-1/2)^2) + rnorm(n,0,1)

x<-0.5*(rlnorm(n,0,1)+rnorm(n,10,4))
y<-sin(0.4*x)+2*rnorm(n,0,1)
n<-88
plot(x,y)
CV<-function(x,y,Ker,h){result<-0
                        for (i in 1:length(x)){result<-result + ((y[i]-RegKer(x[i],x[-i],y[-i],Ker,h))^2)*wl(x[i])}
                        return((1/n)*result)}
CVvec<-Vectorize(CV, vectorize.args = 'h')
h<-seq(hmean,8,0.01)
cv<-CVvec(x,y,Gauss,h)
hminCV<-h[which(min(cv)==cv)]
plot(h,cv,type='l',xaxt='n',ylab='CV(h)', main='Función CV(h)')
points(c(hmean,hminCV),c(CV(x,y,Gauss,hmean),CV(x,y,Gauss,hminCV)), col=c('red','blue'),pch=4)
axis(1,at=c(hmean,2,3,4,5,6,7,8),labels=c(round(hmean,1),2,3,4,5,6,7,8))
axis(1,at=hminCV, labels=round(hminCV,2))
legend('topright',legend=c('hmin', 'hCV'), col=c('red','blue'), pch=4)


#################OSCV#################################
#ONE SIDED KERNELS
L1<-function(x){(15/8)*((1-x^2)^2)*(x>=0)*(x<=1)}
cL<-10/7
dL<-1/7
l<-10

OSCV<-function(x,y,Ker,l,b){result<-0
                            xw<-x[(l+1):(n-l)]; yw<-y[(l+1):(n-l)]
                            for (i in 1:length(xw)){result<-result + ((yw[i]-RegKer(xw[i],xw[-i],yw[-i],Ker,b))^2)}
                            return((1/(n-2*l))*result)}

OSCV<-Vectorize(OSCV, vectorize.args = 'b')
b<-seq(0,8,0.01)
os<-OSCV(x,y,Gauss,10,b)

bmin<-function(x){index<-seq((l+1),(n-l),1)
                  m<-n-2*l
                  h<-rep(0,m)
                  j<-1
                  xord<-sort(x)
                  for(i in index){h[j]<-min(xord[i]-xord[i-1],xord[i+1]-xord[i])
                                               j<-j+1}
                  return(h)}
bmean<-max(bmin(x))
bminOSCV<-b[which(min(os[6:length(os)])==os)]
plot(b,os,type='l',xaxt='n',ylab='OSCV(b)', main='Función OSCV(b)')
points(c(bmean,bminOSCV),c(OSCV(x,y,Gauss,10,bmean),OSCV(x,y,Gauss,10,bminOSCV)), col=c('red','blue'),pch=4)
axis(1,at=c(bmean,2,3,4,5,6,7,8),labels=c(round(bmean,1),2,3,4,5,6,7,8))
axis(1,at=bminOSCV, labels=round(bminOSCV,2))
legend('topright',legend=c('bmin', 'bOSCV'), col=c('red','blue'), pch=4)



#######ETHANOL##############
x<-ethanol$E
y<-ethanol$NOx
library(kernelreg)
plot(x,y)
G<-function(x,y,Ker,h){result<-0
for (i in 1:length(x)){
  result<-result + ((y[i]-RegKer(x[i],x,y,Ker,h))^2)*varxi((Ker((x[i]-x[i])/h))/(sum(Ker((x[i]-x)/h))))*wl(x[i])}
return((1/n)*result)}

Gvec<-Vectorize(G, vectorize.args = 'h')
wl<-function(x){(0.7<x)*(x<1.1)}
hmin<-function(x){xord<-sort(x)
                  m<-sum(wl(xord))
                  h<-rep(0,m)
                  j<-1
                  for(i in which(wl(xord)==1)){h[j]<-min(xord[i]-xord[i-1],xord[i+1]-xord[i])
                  j<-j+1}
                  return(h)}
hmean<-max(hmin(x))
h<-seq(0,0.8,0.01)
gh<-Gvec(x,y,Gauss,h)
hminG<-h[which(min(gh[2:length(gh)])==gh)]
plot(h,gh,type='l',xaxt='n',ylab='G(h)', main='Función penalizadora G(h) sobre los datos de Ethanol')
points(hmean,G(x,y,Gauss,hmean), col='red',pch=4)
points(hminG,G(x,y,Gauss,hminG), col='blue',pch=4)
axis(1,at=c(hmean,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c(round(hmean,2),0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))
legend('bottomright',legend=c('hmin', 'hG'), col=c('red','blue'), pch=4)

CV<-function(x,y,Ker,h){result<-0
for (i in 1:length(x)){result<-result + ((y[i]-RegKer(x[i],x[-i],y[-i],Ker,h))^2)*wl(x[i])}
return((1/n)*result)}
CVvec<-Vectorize(CV, vectorize.args = 'h')
h<-seq(0,1,0.01)
cv<-CVvec(x,y,Gauss,h)
hminCV<-h[which(min(cv[2:length(cv)])==cv)]
plot(h,cv,type='l',xaxt='n',ylab='CV(h)', main='Función CV(h) sobre los datos de Ethanol')
points(c(hmean,hminCV),c(CV(x,y,Gauss,hmean),CV(x,y,Gauss,hminCV)), col=c('red','blue'),pch=4)
axis(1,at=c(hmean,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),labels=c(round(hmean,2),0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
legend('bottomright',legend=c('hmin', 'hCV'), col=c('red','blue'), pch=4)

prediction1<-RegKer(sort(ethanol$E), ethanol$E, ethanol$NOx, Gauss, 0.022)
prediction2<-RegKer(sort(ethanol$E), ethanol$E, ethanol$NOx, Gauss, 0.05)
plot(ethanol$E,ethanol$NOx,  ylab = 'NOx', xlab = 'E', main='Ethanol')
lines(sort(ethanol$E),prediction1, col='red')
lines(sort(ethanol$E),prediction2, col='blue')
legend('topright',legend=c('Predicción h=hmin','Predicción h=0.05'),col=c('red','blue'),lty=1)


