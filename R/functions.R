
##REGRESION KERNEL####
#DEFINICION DE LAS FUNCIONES NÚCLEO
Gauss<-function(x){(1/sqrt(2*pi))*exp((-x^2)/2)}
Naive<-function(x){rep(0.5,length(x))*(abs(x)<=1)}
Triang<-function(x) {(1-abs(x))*(abs(x)<=1)}
Epane<-function(x) {(3/4)*(1-x^2)*(abs(x)<=1)}
Biw<-function(x) {((15/16)*(1-x^2)^2)*(abs(x)<=1)}
triw<-function(x) {((35/32)*(1-x^2)^3)*(abs(x)<=1)}
cosen<-function(x) {(pi/4)*cos((pi/2)*x)*(abs(x)<=1)}
#FUNCION QUE PROPORCIONA ESTIMACION DE LA FREGRESION A TRAVES DEL ESTIMADOR DE NADARAYA-WATSON
RegKer<-function(x,X,Y,Ker,h){rk<-rep(0,length(x))
                              for (i in 1:length(x)){rk[i]<-(sum(Ker((x[i]-X)/h)*Y))/(sum(Ker((x[i]-X)/h)))}
                              return(rk)}


#VALORES DE LAS CONSTANTES DK Y CK PARA LOS DISTINTOS NUCLEOS
dk<-function(x){if(x=='Naive'){1/3}
                else if(x=='Gauss'){1}
                else if(x=='Triangular'){1/6}
                else if(x=='Epane'){1/5}
                else if(x=='biw'){1/7}
                else if(x=='tri'){1/9}
                else if(x=='cos'){1-8/(pi^2)}}
ck<-function(x){if(x=='Naive'){1/2}
                else if(x=='Gauss'){1/(2*sqrt(pi))}
                else if(x=='Triangular'){2/3}
                else if(x=='Epa'){3/5}
                else if(x=='biw'){5/7}
                else if(x=='tri'){350/429}
                else if(x=='cos'){(pi^2)/16}}

#VARIANCE ESTIMATOR
varest<-function(x,X,Y,Ker,h){rk<-rep(0,length(x))
                             for (i in 1:length(x)){rk[i]<-(sum(Ker((x[i]-X)/h)*(Y-RegKer(x,X,Y,Ker,h)[i])^2)/(sum(Ker((x[i]-X)/h))))}
                             return(rk)}
#INTERVALOS DE CONFIANZA
intervalosasintoticos<-function(x,X,Y,Ker,ker,h){cl<-rep(0,length(x))
                                                 cup<-rep(0,length(x))
                                                 for (i in 1:length(x)){cl[i]<-RegKer(x,X,Y,Ker,h)[i]-z*sqrt(ck(ker)/(n*h*sum(Ker((x[i]-X)/h))))}
                                                 for (i in 1:length(x)){cup[i]<-RegKer(x,X,Y,Ker,h)[i]+z*sqrt(ck(ker)/(n*h*sum(Ker((x[i]-X)/h))))}
                                                 return(list('cl'=cl, 'cup'=cup))}

intervalosbootstrap<-function(x,y){clboot<-rep(0,length(x))
                                   cupboot<-rep(0,length(x))
                                   for(i in 1:length(x)){bootmu<-replicate(1000,{samp<-sample(1:88,88,replace=TRUE)
                                                                                 bootx<-x[samp]; booty<-y[samp]
                                                                                 RegKer(x[i],bootx,booty, Epane, 0.1)})
                                                         cupboot[i]<-unname(quantile(bootmu,0.975))
                                                         clboot[i]<-unname(quantile(bootmu,0.025))}
                                                         return(list('cl'=clboot, 'cup'=cupboot))}


wildboots<-function(x,y,h,g,dist){prediction<-RegKer(x,x,y, Epane, h)
residuals<-y-prediction
clwild<-rep(0,length(x))
cupwild<-rep(0,length(x))

if (dist=='Mammen'){p<-(sqrt(5)+1)/(2*sqrt(5))
for (i in 1:length(residuals)){muboot<-replicate(1000,{for (j in 1:length(residuals)){booteps[j]<-sample(c(-(sqrt(5)-1)*residuals[j]/2,(sqrt(5)+1)*residuals[j]/2),1,p=c(p,1-p))}
  yboot<-RegKer(x,x,y,Epane,g)+booteps
  RegKer(x[i],x,yboot,Epane,h)})
clwild[i]<-unname(quantile(muboot,0.025))
cupwild[i]<-unname(quantile(muboot,0.975))}}

else if (dist=='Rademacher'){for (i in 1:length(residuals)){muboot<-replicate(1000,{for (j in 1:length(residuals)){booteps[j]<-sample(c(-residuals[j],residuals[j]),1,p=c(0.5,0.5))}
  yboot<-RegKer(x,x,y,Epane,g)+booteps
  RegKer(x[i],x,yboot,Epane,h)})
clwild[i]<-unname(quantile(muboot,0.025))
cupwild[i]<-unname(quantile(muboot,0.975))}}
return(list('cl'=clwild,'cup'=cupwild))}

bandabonf<-function(x,X,Y,Ker,ker,h,alpha){alphan<-alpha/n
n<-length(x)
z<-qnorm(c(1-(alphan/2)), mean=0, sd=1, lower.tail=TRUE)
cl<-rep(0,length(x)); cup<-rep(0,length(x))
for (i in 1:length(x)){cl[i]<-RegKer(x,X,Y,Ker,h)[i]-z*sqrt(ck(ker)/(n*h*sum(Ker((x[i]-X)/h))))}
for (i in 1:length(x)){cup[i]<-RegKer(x,X,Y,Ker,h)[i]+z*sqrt(ck(ker)/(n*h*sum(Ker((x[i]-X)/h))))}
return(list('cl'=cl,'cup'=cup))}

bandasbootstrap<-function(x,y){n<-length(x)
clboot<-rep(0,n)
cupboot<-rep(0,n)
prediction<-RegKer(x,x,y,Epane,0.15)
bootmu<-replicate(1000,{samp<-sample(1:n,n,replace=TRUE)
bootx<-x[samp]; booty<-y[samp]
max(abs(RegKer(x,bootx,booty, Epane, 0.15)-prediction))})
distboot<-unname(quantile(bootmu,0.95))
return(list('cup'=RegKer(sort(x),x,y,Epane,0.15)+distboot,'cl'=RegKer(sort(x),x,y,Epane,0.15)-distboot))}
