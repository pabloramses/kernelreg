library(kernelreg)
library(locfit)
data(ethanol)

x<-ethanol$E
y<-ethanol$NOx
n<-length(x)
h<-0.05
pred<-RegKer(x,x,y,Epane,0.05)

###INTERVALO DE CONFIANZA ASINTOTICO
alpha<-0.05
z<-qnorm(c(0.975), mean=0, sd=1, lower.tail=TRUE)
intervalosasintoticos<-function(x,X,Y,Ker,ker,h){cl<-rep(0,length(x))
                              cup<-rep(0,length(x))
                              for (i in 1:length(x)){cl[i]<-RegKer(x,X,Y,Ker,h)[i]-z*sqrt(ck(ker)/(n*h*sum(Ker((x[i]-X)/h))))}
                              for (i in 1:length(x)){cup[i]<-RegKer(x,X,Y,Ker,h)[i]+z*sqrt(ck(ker)/(n*h*sum(Ker((x[i]-X)/h))))}
                              return(list('cl'=cl, 'cup'=cup))}
int<-intervalosasintoticos(x,x,y,Epane,'Epa',h)
up<-int$cup
low<-int$cl
dist<-mean(up-low)
dist
plot(x,y,ylab = 'NOx', xlab = 'E', main='Estimación por intervalo asintótico de la función de regresión')
points(sort(x),RegKer(sort(x),x,y,Epane,0.05),col='red',pch=4)
points(x,int$cl,col='blue',pch=4)
points(x,int$cup,col='blue',pch=4)
legend('topright', legend=c('Estimacion', 'Intervalo de confianza'), col=c('red','blue'), pch=4)
legend(0.8, 0.66, "Amplitud media: 0.578", box.col = "lightblue", bg = "lightblue", adj = 0.1)
###INTERVALO DE CONFIANZA BOOTSTRAP

intervalosbootstrap<-function(x,y){clboot<-rep(0,length(x))
                                   cupboot<-rep(0,length(x))
                                   for(i in 1:length(x)){bootmu<-replicate(1000,{samp<-sample(1:88,88,replace=TRUE)
                                                                                 bootx<-x[samp]; booty<-y[samp]
                                                                                 RegKer(x[i],bootx,booty, Epane, 0.05)})
                                                         cupboot[i]<-unname(quantile(bootmu,0.975,na.rm=TRUE))
                                                         clboot[i]<-unname(quantile(bootmu,0.025,na.rm=TRUE))}
                                   return(list('cl'=clboot, 'cup'=cupboot))}

int2<-intervalosbootstrap(x,y)
dist2<-mean(int2$cup-int2$cl)
plot(x,y,ylab = 'NOx', xlab = 'E', main='Estimación por intervalo bootstrap de la función de regresión')
points(sort(x),RegKer(sort(x),x,y,Epane,0.05),col='red',pch=4)
points(x,int2$cl,col='blue',pch=4)
points(x,int2$cup,col='blue',pch=4)
legend('topright', legend=c('Estimacion', 'Intervalo de confianza'), col=c('red','blue'), pch=4)
legend(0.8, 0.66, "Amplitud media: 0.409", box.col = "lightblue", bg = "lightblue", adj = 0.1)
dist2
#WILD BOOTSTRAPPING
wildboots<-function(x,y,h,g,dist){prediction<-RegKer(x,x,y, Epane, h)
                                  residuals<-y-prediction
                                  clwild<-rep(0,length(x))
                                  cupwild<-rep(0,length(x))
                                  if (dist=='Mammen'){p<-(sqrt(5)+1)/(2*sqrt(5))
                                                      for (i in 1:length(residuals)){booteps<-rep(0,length(residuals))
                                                                                     muboot<-replicate(1000,{for (j in 1:length(residuals)){booteps[j]<-sample(c(-(sqrt(5)-1)*residuals[j]/2,(sqrt(5)+1)*residuals[j]/2),1,p=c(p,1-p))}
                                                      yboot<-RegKer(x,x,y,Epane,g)+booteps
                                                      RegKer(x[i],x,yboot,Epane,h)})
                                                      clwild[i]<-unname(quantile(muboot,0.025))
                                                      cupwild[i]<-unname(quantile(muboot,0.975))}}

                                  else if (dist=='Rademacher'){for (i in 1:length(residuals)){booteps<-rep(0,length(residuals))
                                                                                              muboot<-replicate(1000,{for (j in 1:length(residuals)){booteps[j]<-sample(c(-residuals[j],residuals[j]),1,p=c(0.5,0.5))}
                                                               yboot<-RegKer(x,x,y,Epane,g)+booteps
                                                               RegKer(x[i],x,yboot,Epane,h)})
                                                               clwild[i]<-unname(quantile(muboot,0.025))
                                                               cupwild[i]<-unname(quantile(muboot,0.975))}}
                                  return(list('cl'=clwild,'cup'=cupwild))}
intwild1<-wildboots(x,y,0.05,0.1,'Mammen')
intwild2<-wildboots(x,y,0.05,0.1,'Rademacher')
par(mfrow=c(1,1))
dist3<-mean(intwild1$cup-intwild1$cl)
plot(x,y,ylab = 'NOx', xlab = 'E',main = 'Intervalo Wildbootstrap con distribución de Mammen')
points(sort(x),RegKer(sort(x),x,y,Epane,0.1),col='red',pch=4)
points(x,intwild1$cl,col='blue',pch=4)
points(x,intwild1$cup,col='blue',pch=4)
legend('topright', legend=c('Estimacion', 'Int.confianza'), col=c('red','blue'), pch=4)
legend(0.8, 0.66, "Amplitud media: 0.308", box.col = "lightblue", bg = "lightblue", adj = 0.1)
dist4<-mean(intwild2$cup-intwild2$cl)
plot(x,y,ylab = 'NOx', xlab = 'E',main = 'Intervalo Wildbootstrap con distribución de Rademacher')
points(sort(x),RegKer(sort(x),x,y,Epane,0.1),col='red',pch=4)
points(x,intwild2$cup,col='purple',pch=4)
points(x,intwild2$cl,col='purple',pch=4)
legend('topright', legend=c('Estimacion', 'Int.confianza'), col=c('red','blue'), pch=4)
legend(0.8, 0.66, "Amplitud media: 0.301", box.col = "lightblue", bg = "lightblue", adj = 0.1)
#######BANDA POR BONFERRONI#####################
alpha<-0.1

bandabonf<-function(x,X,Y,Ker,ker,h,alpha){alphan<-alpha/n
                                           n<-length(x)
                                           z<-qnorm(c(1-(alphan/2)), mean=0, sd=1, lower.tail=TRUE)
                                           cl<-rep(0,length(x)); cup<-rep(0,length(x))
                                           for (i in 1:length(x)){cl[i]<-RegKer(x,X,Y,Ker,h)[i]-z*sqrt(ck(ker)/(n*h*sum(Ker((x[i]-X)/h))))}
                                           for (i in 1:length(x)){cup[i]<-RegKer(x,X,Y,Ker,h)[i]+z*sqrt(ck(ker)/(n*h*sum(Ker((x[i]-X)/h))))}
                                           return(list('cl'=cl,'cup'=cup))}
bandabonferroni<-bandabonf(sort(x),x,y,Epane,'Epa',h,alpha)

dist5<-mean(bandabonferroni$cup-bandabonferroni$cl)
plot(x,y,ylab = 'NOx', xlab = 'E', main='Banda de confianza de Bonferroni de la función de regresión')
lines(sort(x),RegKer(sort(x),x,y,Epane,0.05),col='red',pch=4)
lines(sort(x),bandabonferroni$cup,col='blue',pch=4)
lines(sort(x),bandabonferroni$cl,col='blue',pch=4)
legend('topright', legend=c('Estimacion', 'Banda de confianza'), col=c('red','blue'), lty=1)
legend(0.8, 0.66, "Amplitud esperada: 0.961", box.col = "lightblue", bg = "lightblue", adj = 0.1)
########BANDA BOOTSTRAP########################

bandasbootstrap<-function(x,y){n<-length(x)
                               clboot<-rep(0,n)
                               cupboot<-rep(0,n)
                               prediction<-RegKer(x,x,y,Epane,0.05)
                               bootmu<-replicate(1000,{samp<-sample(1:n,n,replace=TRUE)
                                                       bootx<-x[samp]; booty<-y[samp]
                                                       max(abs(RegKer(x,bootx,booty, Epane, 0.05)-prediction))})
                               distboot<-unname(quantile(bootmu,0.95,na.rm=TRUE))
                               return(list('cup'=RegKer(sort(x),x,y,Epane,0.05)+distboot,'cl'=RegKer(sort(x),x,y,Epane,0.05)-distboot))}

int<-bandasbootstrap(x,y)

plot(x,y,ylab = 'NOx', xlab = 'E', main='Banda de confianza bootstrap de la función de regresión', ylim=c(min(int$cl),max(int$cup)))
lines(sort(x),RegKer(sort(x),x,y,Epane,0.05),col='red',lty=1)
lines(sort(x),int$cup,col='blue',lty=1)
lines(sort(x),int$cl,col='blue',lty=1)
legend('topright', legend=c('Estimacion', 'Banda conf.'), col=c('red','blue'), lty=1)
dist=mean(int$cup-int$cl)
legend(0.8, 0.6, "Amplitudm esperada:1.070", box.col = "lightblue", bg = "lightblue", adj = 0.1)


