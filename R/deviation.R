library(kernelreg)
library(locfit)
data(ethanol)

varest<-function(x,X,Y,Ker,h){rk<-rep(0,length(x))
                              for (i in 1:length(x)){rk[i]<-(sum(Ker((x[i]-X)/h)*(Y-RegKer(x,X,Y,Ker,h)[i])^2)/(sum(Ker((x[i]-X)/h))))}
                              return(rk)}
densest<-function(x,X,Ker,h){rk<-rep(0,length(x))
                               n<-length(x)
                               for (i in 1:length(x)){rk[i]<-(1/(n*h))*sum(Ker((x[i]-X)/h))}
                               return(rk)}


x<-ethanol$E
y<-ethanol$NOx
n<-length(x)
delta<-8
alpha<-0.05





maxdev<-function(x,X,Y,Ker,delta,alpha){n<-length(X)
                                  h<-0.05
                                  cl<-rep(0,n)
                                  cup<-rep(0,n)
                                  ck<-3/5
                                  dn<-((2*log(1/h))^(1/2))+((2*log(1/h))^(-1/2))*log(5/(4*pi))
                                  calpha<--log((1/2)*log((1/(1-alpha))))
                                  muh<-RegKer(x,X,Y,Ker,h)
                                  vare<-varest(x,X,Y,Ker,h)
                                  dens<-densest(x,X,Ker,h)
                                  for (i in 1:length(x)){cl[i]<-muh[i]-((calpha/((2*log(1/h))^(1/2)+dn))*(((vare[i]^2)/(dens[i]*n*h))^(1/2))/ck)
                                                         cup[i]<-muh[i]+((calpha/((2*log(1/h))^(1/2)+dn))*(((vare[i]^2)/(dens[i]*n*h))^(1/2))/ck)}
                                  return(list('cl'=cl,'cup'=cup))}

int<-maxdev(sort(x),x,y,Epane,delta,alpha)

plot(x,y, xlab='E',ylab='NOx',main='Banda de confianza asintótica')
lines(sort(x),int$cl,col='blue')
lines(sort(x),int$cup,col='blue')
lines(sort(x),RegKer(sort(x),x,y,Epane,0.05),col='red')
legend('topright', legend=c('Estimación','Banda de confianza'),col=c('red','blue'), lty=1)
legend(0.8, 0.65, "Amplitud esperada:0.126", box.col = "lightblue", bg = "lightblue", adj = 0.1)

h<-n^(-delta)
delta<-1/2
dist<-mean(int$cup-int$cl)

