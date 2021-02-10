library(locfit)
data(ethanol)
library(kernelreg)
RegKerKn<-function(x,X,Y,Ker,alpha){n<-length(x); k<-floor(alpha*n)
                                    rk<-rep(0,n)
                                    for (i in 1:n){dist<-abs(x[i]-X)
                                                   h<-sort(dist)[k]
                                                   rk[i]<-(sum(Ker((x[i]-X)/h)*Y))/(sum(Ker((x[i]-X)/h)))}
                                    return(rk)}
x<-ethanol$E
y<-ethanol$NOx
pred1<-RegKerKn(sort(x),x,y,Epane,0.2)
pred2<-RegKerKn(sort(x),x,y,Epane,0.4)
pred3<-RegKerKn(sort(x),x,y,Epane,0.6)
pred4<-RegKerKn(sort(x),x,y,Epane,0.9)
plot(x,y, xlab='E',ylab='NOx', main='Elección del ancho de banda a través del método de k-vecinos')
lines(sort(x),pred4,col='yellow')
lines(sort(x),pred3,col='green')
lines(sort(x),pred2,col='blue')
lines(sort(x),pred1,col='red')
legend('topright', legend=c('alpha=0.2','alpha=0.4','alpha=0.6','alpha=0.8'), col=c('red','blue','green','yellow'),lty=1)
