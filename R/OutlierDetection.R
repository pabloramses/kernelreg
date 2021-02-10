library(kernelreg)
library(locfit)
#ONE SIDED KERNELS
L1<-function(x){(15/8)*((1-x^2)^2)*(x>=0)*(x<=1)}
#simulacion outliers etanol
data(ethanol)
x<-ethanol$E
y<-ethanol$NOx
n<-length(x)
h<-0.05
#DATA SPLITTING
#split into normal data and data with anomalies
train_size<-floor(n*0.6)
val_size<-ceiling(n*0.2)
test_size<-n-train_size-val_size
train_index<-sample(seq_len(n),train_size)
xtrain<-x[train_index]; ytrain<-y[train_index]
xnotrain<-x[-train_index]; ynotrain<-y[-train_index]
val_index<-sample(seq_len(n-train_size),val_size)
xval<-xnotrain[val_index]; yval<-ynotrain[val_index]
xtest<-xnotrain[-val_index]; ytest<-ynotrain[-val_index]
#simulated outliers
w<-sample(1:length(ytest),floor(length(ytest)*0.4))
ytest[w]<-ytest[w]+(-1)^(w)*rnorm(floor(length(ytest)*0.4),1.5,0.1)

#DETECTION
prediction<-RegKer(xval,xtrain,ytrain,Gauss,0.05)
val_loss<-prediction-yval
threshold<-max(abs(val_loss))
####OUTLIER DETECTION:
testpred<-RegKer(xtest,xtrain,ytrain,Gauss,0.05)
test_loss<-abs(testpred-ytest)
outlier<-which(test_loss>threshold)
#representacion
plot(c(xtest,xtrain,xval),c(ytest,ytrain,yval),xlab='E',ylab='NOx',main='Representación datos Ethanol con outliers simulados')
points(xtest[outlier],ytest[outlier],col='red')
legend('topright',legend='Outlier', col='red', pch=1)

#MATRIZ DE CONFUSION
confmat<-function(outliersdetected,realoutliers,n){falsopositivo<-0
truepositivo<-0
for (i in 1:length(outliersdetected)){contador<-0
for (j in 1:length(realoutliers)){if (realoutliers[j]==outliersdetected[i]){contador<-contador+1}}
if(contador==1){truepositivo<-truepositivo+1}
else if (contador == 0) {falsopositivo<-falsopositivo+1}}
return(matrix(c(truepositivo,falsopositivo,length(realoutliers)-length(outliersdetected)+falsopositivo, n-length(realoutliers)-falsopositivo),nrow=2,byrow=FALSE))}
#ESPEC Y SENS
espsens<-function(mat){c(mat[1,1]/apply(mat,1,sum)[1],(mat[2,2]-mat[1,2])/apply(mat,1,sum)[2])}
#####TRAINING DATA
n<-500
ntrain<-floor(n*0.6)
nval<-ceiling(n*0.2)
ntest<-n-ntrain-ntest
x<-runif(n,-1,1)
y<-x^2+rnorm(n,0,0.2)
resultado3<-replicate(1000,{train_index<-sample(seq_len(n),train_size)
                            xtrain<-x[train_index]; ytrain<-y[train_index]
                            xnotrain<-x[-train_index]; ynotrain<-y[-train_index]
                            val_index<-sample(seq_len(n-train_size),val_size)
                            xval<-xnotrain[val_index]; yval<-ynotrain[val_index]
                            xtest<-xnotrain[-val_index]; ytest<-ynotrain[-val_index]
                            w<-sample(1:ntest,floor(ntest*0.2))
                            ytest[w]<-ytest[w]+(-1)^(w)*rnorm(floor(ntest*0.2),0.75,0.1)
                            ###CALCULO EL THRESHOLD:
                            prediction<-RegKer(xval,xtrain,ytrain,Gauss,0.1)
                            val_loss<-prediction-yval
                            threshold<-max(abs(val_loss[5:length(val_loss)]))
                            ####OUTLIER DETECTION:
                            testpred<-RegKer(xtest,xtrain,ytrain,Gauss,0.1)
                            test_loss<-abs(testpred-ytest)
                            outlier<-which(test_loss>threshold)
                            if (length(outlier)>0){mat<-confmat(outlier,w,n)
                                                   espsens(mat)}
                            else {c(0,1)}})
apply(resultado3,1,mean)
plot(xtest,ytest,main='Detección de Outliers')
points(xtest[outlier],ytest[outlier],col='red')
legend('topright',legend=c('outlier'),col='red',pch=1)

outlier
w
#SERIES TEMPORALES
#####TRAINING DATA
n<-201
x<-seq(-1,1,0.01)
resultados<-replicate(1000,{ytrain<-x^2+rnorm(n,0,0.2)
                            #####VALIDATION DATAxtrain<-runif(n,-1,1)
                            yval<-x^2+rnorm(n,0,0.2)
                            #####TEST DATA
                            ytest<-x^2+rnorm(n,0,0.2)
                            w<-sample(1:n,floor(n*0.1))
                            ytest[w]<-ytest[w]+(-1)^(w)*rnorm(floor(n*0.1),0.75,0.1)
                            ###CALCULO EL THRESHOLD:
                            prediction<-RegKer(x,x,ytrain,L1,0.5)
                            val_loss<-prediction-yval
                            threshold<-max(abs(val_loss))
                            ####OUTLIER DETECTION:
                            testpred<-RegKer(x,x,ytrain,L1,0.5)
                            test_loss<-abs(testpred-ytest)
                            outlier<-which(test_loss>threshold)
                            if (length(outlier)>0){mat<-confmat(outlier,w,n)
                                                   espsens(mat)}
                            else {c(0,1)}})

apply(resultados,1,mean)
plot(x,ytest, main='Detección de Outlier en Serie Temporal')
points(x[outlier],ytest[outlier],col='red')
legend('bottomright',legend='outlier',col='red',pch=1)
