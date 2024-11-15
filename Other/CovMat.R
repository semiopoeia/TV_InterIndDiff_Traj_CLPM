library(MASS)
library(lavaan)
library(tidyverse)
set.seed(1507)

##########Structured Residual Generating Model#############
#set mean structure
mu<-c(50.36,3.34,50.35,3.29,rep(0,12))

#set correlation matrix
R<-matrix(c(
  1,	0.524,	0.912,	0.561,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0.524,	1,	0.424,	0.811,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0.912	,0.424,	1,	0.617,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0.561	,0.811,	0.617,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0,	0	,0,	0,	1,	0,	0,	0,	0,	0,	0.25,	0,	0,	0,	0,	0,
  0,	0	,0,	0,	0,	1,	0,	0,	0,	0,	0,	0.25,	0,	0,	0,	0,
  0,	0,0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0.25,	0,	0,	0,
  0,	0	,0,	0	,0,	0,	0,	1,	0,	0,	0,	0,	0,	0.25,	0,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0.25,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0.25,
  0,	0	,0,	0,	0.25,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,
  0,	0	,0,	0,	0,	0.25,	0,	0	,0	,0	,0	,1	,0	,0	,0	,0,
  0,	0	,0,	0,	0,	0,	0.25,	0,	0,	0,	0,	0,	1,	0,	0,	0,
  0,	0	,0,	0,	0,	0,	0,	0.25,	0,	0,	0,	0,	0,	1,	0,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	0.25,	0,	0,	0,	0,	0,	1,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	0,	0.25,	0,	0,	0,	0,	0,	1
  
),16,16)

#create covariance from correlation
#Variances
V<-diag(c(
	73.06^.5,1.84^.5,79.34^.5,2.15^.5,
	30.96^.5,28.74^.5,22.53^.5,22.30^.5,15.68^.5,22.6^.5,
	24.52^.5,25.15^.5,23.47^.5,22.39^.5,23.01^.5,23.99^.5))

Sigma<-V%*%R%*%V

rownames(Sigma)<-c("Ix","Sx","Iy","Sy",
                   "Vx1","Vx2","Vx3","Vx4","Vx5","Vx6",
                   "Vy1","Vy2","Vy3","Vy4","Vy5","Vy6")
colnames(Sigma)<-c("Ix","Sx","Iy","Sy",
                   "Vx1","Vx2","Vx3","Vx4","Vx5","Vx6",
                   "Vy1","Vy2","Vy3","Vy4","Vy5","Vy6")

#raw dataframe from covariance matrix & mean structure
df<-mvrnorm(3109, mu, Sigma, tol = 1e-6, empirical = FALSE)

#create structured residuals dataframe
library(tidyverse)
SR_df<-
df%>%as_tibble()%>%
	mutate(x1=Ix+Vx1) %>%
	mutate(y1=Iy+Vy1) %>%
	
	mutate(x2=(Ix+Sx)+(0.35*Vx1+0.05*Vy1)+Vx2) %>%
	mutate(y2=(Iy+Sy)+(0.32*Vy1+0.05*Vx1)+Vy2) %>%
	
	mutate(x3=(Ix+2.2*Sx)+(0.35*Vx2+0.05*Vy2)+Vx3) %>%
	mutate(y3=(Iy+2.3*Sy)+(0.32*Vy2+0.05*Vx2)+Vy3) %>%
	
	mutate(x4=(Ix+2.84*Sx)+(0.35*Vx3+0.05*Vy3)+Vx4) %>%
	mutate(y4=(Iy+3.55*Sy)+(0.32*Vy3+0.05*Vx3)+Vy4) %>%
	
	mutate(x5=(Ix+3.46*Sx)+(0.35*Vx4+0.05*Vy4)+Vx5) %>%
	mutate(y5=(Iy+4.31*Sy)+(0.32*Vy4+0.05*Vx4)+Vy5) %>%
	
	mutate(x6=(Ix+3.78*Sx)+(0.35*Vx5+0.05*Vy5)+Vx6) %>%
	mutate(y6=(Iy+4.59*Sy)+(0.32*Vy5+0.05*Vx5)+Vy6) %>%

	subset(select=c(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6))

write.table(SR_df, "C:/Users/pws5/Documents/DPM_simul_Dat/SR_df_i.dat",sep="\t",
		col.names=FALSE,row.names=FALSE)


###############Observed Variable Generating MOdel############
omu<-c(29.17,28.62,50.35,50.35,rep(0,10))

#set correlation matrix
oR<-matrix(c(
  1,0.821,0.789,0.701,0,0,0,0,0,0,0,0,0,0,
  0.821,1,0.692,0.819,0,0,0,0,0,0,0,0,0,0,
  0.789,0.692,1,0.716,0,0,0,0,0,0,0,0,0,0,
  0.701,0.819,0.716,1,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,0.28,0,0,0,0,
  0,0,0,0,0,1,0,0,0,0,0.252,0,0,0,
  0,0,0,0,0,0,1,0,0,0,0,0.194,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,0.195,0,
  0,0,0,0,0,0,0,0,1,0,0,0,0,0.306,
  0,0,0,0,0.28,0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,0.252,0,0,0,0,1,0,0,0,
  0,0,0,0,0,0,0.194,0,0,0,0,1,0,0,
  0,0,0,0,0,0,0,0.195,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,0.306,0,0,0,0,1
),14,14)

#create covariance from correlation
#Variances
oV<-diag(c(
  29.78^.5,36.21^.5,104.283^.5,103.81^.5,
  30.16^.5,23.36^.5,22.63^.5,16.6^.5,24.7^.5,
  27.22^.5,25^.5,24.61^.5,24.34^.5,27.13^.5))

oSigma<-oV%*%oR%*%oV

rownames(oSigma)<-c("Bx","By","x1","y1",
                   "Vx2","Vx3","Vx4","Vx5","Vx6",
                   "Vy2","Vy3","Vy4","Vy5","Vy6")
colnames(oSigma)<-c("Bx","By","x1","y1",
                "Vx2","Vx3","Vx4","Vx5","Vx6",
                "Vy2","Vy3","Vy4","Vy5","Vy6")

#raw dataframe from covariance matrix & mean structure
odf<-mvrnorm(3109, omu, oSigma, tol = 1e-6, empirical = FALSE)

#create structured residuals dataframe
library(tidyverse)
OV_df<-
  odf%>%as_tibble()%>%
  mutate(x1=x1) %>%
  mutate(y1=y1) %>%
  
  mutate(x2=(Bx)+(0.394*x1+0.09*y1)+Vx2) %>%
  mutate(y2=(By)+(0.405*y1+0.09*x1)+Vy2) %>%
  
  mutate(x3=(1.086*Bx)+(0.394*x2+0.09*y2)+Vx3) %>%
  mutate(y3=(1.094*By)+(0.405*y2+0.09*x2)+Vy3) %>%
  
  mutate(x4=(1.093*Bx)+(0.394*x3+0.09*y3)+Vx4) %>%
  mutate(y4=(1.170*By)+(0.405*y3+0.09*x3)+Vy4) %>%
  
  mutate(x5=(1.121*Bx)+(0.394*x4+0.09*y4)+Vx5) %>%
  mutate(y5=(1.197*By)+(0.405*y4+0.09*x4)+Vy5) %>%
  
  mutate(x6=(1.119*Bx)+(0.394*x5+0.09*y5)+Vx6) %>%
  mutate(y6=(1.182*By)+(0.405*y5+0.09*x5)+Vy6) %>%
  
  subset(select=c(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6))

write.table(SR_df, "C:/Users/scottpw/Desktop/RR_Simul&Files/OV_df.dat",sep="\t",
            col.names=FALSE,row.names=FALSE)
