library(kernelreg)
#########################cO2#########################
data(co2)
#tendencia
dataco2<-monco2$CO2[385:744]

months<-seq_len(length(dataco2))
plot(months,dataco2,type='l',xaxt='n',ylab='emisions',main='CO2 mensual en Hawaii')
pred<-RegKer(months,months,dataco2,Gauss,7)
lines(months,pred,col='red')
axis(1,at=seq(0,350,50),labels=c('1-1991','2-1995','4-1999','6-2003','8-2007','10-2011','12-2015','2-2020'))
legend('bottomright',legend=c('CO2 emisions','Trend'),col=c('black','red'),lty=1)


detend<-dataco2-pred
predses<-RegKer(months,months,detend,Gauss,1)
plot(months,detend,type='l',ylab='Detrended Serie',xlab='Months',xaxt='n',main='Serie sin tendencia')
lines(predses,col='red')
axis(1,at=seq(0,350,50),labels=c('1-1991','2-1995','4-1999','6-2003','8-2007','10-2011','12-2015','2-2020'))
legend('topleft',legend='Suavizado', col='red',lty=1)

#destendenciao papi
mesfactor<-levels(as.factor(monco2$Mn[385:744]))
media<-c()
for (mes in mesfactor){media<-append(media,mean(predses[monco2$Mn[385:744]==mes]))}
plot(monco2$Mn[385:744],predses,main='Datos agrupados por meses',ylab='Emisiones',xlab='Meses',xaxt='n')
axis(1,at=seq_len(12),labels=c('J','F','M','A','M','J','J','A','S','O','N','D'))
lines(media)
legend('topright',legend='Media por mes',col='black', lty=1)
#estacionalidad todo wapa niño
estacionalidad<-rep(media,length(levels(as.factor(monco2$Yr[385:744]))))
plot(estacionalidad,type='l',xaxt='n',main='Componente estacional de la serie',xlab='Months')
axis(1,at=seq(0,350,50),labels=c('1-1991','2-1995','4-1999','6-2003','8-2007','10-2011','12-2015','2-2020'))

plot(dataco2,type='l',xaxt='n',xlab='Months',ylab='Emisiones co2',main='Serie temporal y estimación STL')
lines(estacionalidad+pred,col='red')
axis(1,at=seq(0,350,50),labels=c('1-1991','2-1995','4-1999','6-2003','8-2007','10-2011','12-2015','2-2020'))
legend('bottomright',legend=c('Datos','Estimación'),col=c('black','red'),lty=1)

plot(dataco2-estacionalidad-pred,type='l',main='Estimación del ruido',xlab='months',ylab='residuals')
axis(1,at=seq(0,350,50),labels=c('1-1991','2-1995','4-1999','6-2003','8-2007','10-2011','12-2015','2-2020'))


