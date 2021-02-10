library(kernelreg)

#parametros del modelo:
f<-function(x){0.5*(abs(x)<=1)}
f1<-function(x){0}
mu<-function(x){x^{3}}
mu1<-function(x){3*x^{2}}
mu2<-function(x){6*x}

k<-4
sigma<-0.5
colores<-c('red','blue','green')
t<-seq(-7,7,0.01)
plot(t, pnorm(t, mean=0, sd=1), type = 'l', xlab='x', ylab='Fn(x)', main= 'Convergencia en Ley de la distribución del estimador NW')
legend('topleft',legend=c('n=10', 'n=100', 'n=1000', 'F.distr N(0,1)'),col=c('red','blue','green', 'black' ),lty=1)
for (i in 1:k){n<-10^i
               h<-n^(-1/5)
               x<-runif(n,-1,1)
               y<-mu(x) + rnorm(n,0,sqrt(sigma))
               sample<-replicate(1000, {x<-runif(n,-1,1)
                                        y<-mu(x) + rnorm(n,0,sqrt(sigma))
                                        RegKer(0,x,y,Epane,h)*sqrt((n*h)/(sigma*ck('Epa')/0.5))})
              cdf<-ecdf(sample)
              lines(t,cdf(t), col=colores[i])}


pred<-RegKer(0,x,y,Epane,h)
sigma(0,x,y,Epane,h)
pred

Y<-y
X<-x
x<-0
(3/4)*(1-x^2)*(abs(x)<1)
(sum(Epane((x[1]-X)/h)*(Y-c[1])^2)/(sum(Epane((x[1]-X)/h))))
Epane((x-X)/h)
Y-c[1]
