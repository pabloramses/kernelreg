#DEFINICION DE LAS FUNCIONES NÚCLEO
Gauss<-function(x){(1/sqrt(2*pi))*exp((-x^2)/2)}
Naive<-function(x){rep(0.5,length(x))*(abs(x)<=1)}
Triang<-function(x) {(1-abs(x))*(abs(x)<=1)}
Epane<-function(x) {(3/4)*(1-x^2)*(abs(x)<=1)}
Biw<-function(x) {((15/16)*(1-x^2)^2)*(abs(x)<=1)}
triw<-function(x) {((35/32)*(1-x^2)^3)*(abs(x)<=1)}
cosen<-function(x) {(pi/4)*cos((pi/2)*x)*(abs(x)<=1)}


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
