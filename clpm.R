library(lavaan)

lsay<-read.table("D:/LSAY.dat",sep="\t",header=F)
View(lsay)

CLPM<-
'V11~~V5
 V12~~V6
 V7~~V13
 V8~~V14
 V9~~V15
 V10~~V16
 
 V16~X*V15+Y*V9
 V15~X*V14+Y*V8
 V14~X*V7+Y*V13
 V7~X*V12+Y*V6
 V12~X*V11+Y*V5
 V10~A*V9+B*V15
 V9~A*V8+B*V14
 V8~A*V13+B*V7
 V13~A*V6+B*V12
 V6~A*V5+B*V11
'

fit<-sem(CLPM,data=lsay,estimator="ML")
summary(fit)
install.packages("lme")
m1<-lme(V12~V6+(1|V4),data=lsay)
library(mediation)
View(lsay)