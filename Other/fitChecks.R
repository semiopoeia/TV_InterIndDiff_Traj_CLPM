library(readxl)
ls<-read_excel("D:/LSAY.xlsx",sheet="LSAY")
ls<-filter(df,x1>0&x2>0&x3>0&x4>0&x5>0&x6>0&y1>0&y2>0&y3>0&y4>0&y5>0&y6>0)

df<-seeDatasets$datasets[[195]]

#fit the model
tic()
fit_gclm<-sem(LGCM_SR_SlpVar,data=SR_df,estimator="ML")
toc()

fitMeasures(fit_gclm)

View(summary(fit_gclm))


fact_mod<-
'
Sx=~x6+x5+x4+x3+x2+x1
Sy=~y6+y5+y4+y3+y2+y1

Sx~~Sy
'

lcm_mod<-
'
Ix=~1*x1+1*x2+1*x3+1*x4+1*x5+1*x6
Sx=~0*x1+1*x2+lx3*x3+lx4*x4+lx5*x5+lx6*x6
Iy=~1*y1+1*y2+1*y3+1*y4+1*y5+1*y6
Sy=~0*y1+1*y2+ly3*y3+ly4*y4+ly5*y5+ly6*y6
'
