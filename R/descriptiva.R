library(locfit)
data(ethanol)
plot(ethanol$E,ethanol$NOx,xlab='E', ylab='NOx', main='E frente a NOx')
boxplot(ethanol$E, main='E')
boxplot(ethanol$NOx, main='NOx')


plot(anos,dataco2, ylab='months', xlab='CO2', main='CO2 frente a Yr')
plot(meses,dataco2, ylab='months', xlab='CO2', main='CO2 frente a Mn')
plot(dataco2, main='Serie temporal CO2')
boxplot(dataco2,main='CO2')

plot(AABA$High,ylab='Highest daily value', xlab='day',main='Datos AABA')
boxplot(AABA$High,main='High')
