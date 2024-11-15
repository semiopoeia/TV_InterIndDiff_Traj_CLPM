#load in libraries
library(MASS)
library(tidyverse)
library(lavaan)
library(MplusAutomation)
library(texreg)
library(furrr)
library(tictoc)

###########################################################
##########LCM (deterministic pure) Generating Model########
###########################################################
tic()
SxSyR<-c(0.1,0.5)
VSx<-c(3.44,13.76)
VSy<-c(4.22,16.88)
CLxy<-0
ARx<-0
ARy<-0
Niter<-1:1000

cond<-crossing(SxSyR,VSx,VSy,CLxy,ARx,ARy,Niter)	

LCM_gen<-function(SxSyR,VSx,VSy,CLxy,ARx,ARy,Niter){

		#set mean structure
mu<-c(50.36,3.32,50.35,3.27,rep(0,12))

#set correlation matrix
R<-matrix(c(
  1,	0.3,	0.5,	0.1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0.3,	1,	0.1,	SxSyR,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0.5	, 0.1,	1,	0.3,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0.1	,SxSyR,	0.3,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0,	0	,0,	0,	1,	0,	0,	0,	0,	0,	0.2,	0,	0,	0,	0,	0,
  0,	0	,0,	0,	0,	1,	0,	0,	0,	0,	0,	0.2,	0,	0,	0,	0,
  0,	0,0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0.2,	0,	0,	0,
  0,	0	,0,	0	,0,	0,	0,	1,	0,	0,	0,	0,	0,	0.2,	0,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0.2,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0.2,
  0,	0	,0,	0,	0.2,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,
  0,	0	,0,	0,	0,	0.2,	0,	0	,0	,0	,0	,1	,0	,0	,0	,0,
  0,	0	,0,	0,	0,	0,	0.2,	0,	0,	0,	0,	0,	1,	0,	0,	0,
  0,	0	,0,	0,	0,	0,	0,	0.2,	0,	0,	0,	0,	0,	1,	0,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	0.2,	0,	0,	0,	0,	0,	1,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	0,	0.2,	0,	0,	0,	0,	0,	1
  
),16,16)

#create covariance from correlation
#Variances
V<-diag(c(
	72.08^.5,VSx^.5,78.90^.5,VSy^.5,
	28^.5,28^.5,28^.5,28^.5,28^.5,28^.5,
	25^.5,25^.5,25^.5,25^.5,25^.5,25^.5))

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

LCM_df<-
df%>%as_tibble()%>%
	mutate(x1=Ix+Vx1) %>%
	mutate(y1=Iy+Vy1) %>%
	
	mutate(x2=(Ix+Sx)+(ARx*Vx1+CLxy*Vy1)+Vx2) %>%
	mutate(y2=(Iy+Sy)+(ARy*Vy1+0*Vx1)+Vy2) %>%
	
	mutate(x3=(Ix+2.2*Sx)+(ARx*Vx2+CLxy*Vy2)+Vx3) %>%
	mutate(y3=(Iy+2.3*Sy)+(ARy*Vy2+0*Vx2)+Vy3) %>%
	
	mutate(x4=(Ix+2.86*Sx)+(ARx*Vx3+CLxy*Vy3)+Vx4) %>%
	mutate(y4=(Iy+3.57*Sy)+(ARy*Vy3+0*Vx3)+Vy4) %>%
	
	mutate(x5=(Ix+3.49*Sx)+(ARx*Vx4+CLxy*Vy4)+Vx5) %>%
	mutate(y5=(Iy+4.35*Sy)+(ARy*Vy4+0*Vx4)+Vy5) %>%
	
	mutate(x6=(Ix+3.8*Sx)+(ARx*Vx5+CLxy*Vy5)+Vx6) %>%
	mutate(y6=(Iy+4.62*Sy)+(ARy*Vy5+0*Vx5)+Vy6) %>%

	subset(select=c(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6))
return(LCM_df)
}
toc()

tic()
LCM_simDats<-
  cond %>%
  group_by(SxSyR,VSx,VSy,CLxy,ARx,ARy,Niter) %>%
  mutate(
    datasets=pmap(list(SxSyR,VSx,VSy,CLxy,ARx,ARy,Niter),
                  possibly(LCM_gen,NA))
  )
toc()

tic()
tbl<-
LCM_simDats%>%
  mutate(GenMod="LCM")%>%
	mutate(LSlpCov=factor(SxSyR,levels=c(0.1,0.5),labels=c("LoCov","HiCov")))%>%
	mutate(LVSx=factor(VSx,levels=c(3.44,13.76),labels=c("LoVSx","HiVSx")))%>%
	mutate(LVSy=factor(VSy,levels=c(4.22,16.88),labels=c("LoVSy","HiVSy")))%>%
	mutate(Cell=paste(LSlpCov,LVSx,LVSy))%>%
	group_by(Cell)

DataList<-group_split(tbl)
toc()

#######################################################################
##########Write out data and input files for MPlus to use##############
#######################################################################
tic()
for(k in 1:8){
  dir.create(paste0("C:/LCM/ConditionSet",k))}
for (j in 1:8){
      for (i in 1:1000){
        write.table(DataList[[j]]$datasets[i],
        paste0("C:/LCM/ConditionSet",j,"/","Rep",i,".dat"),
        sep="\t",row.names=FALSE,col.names=FALSE)
        write.table((paste0("Rep",i,".dat")),
                    paste0("C:/LCM/ConditionSet",j,"/Rep.dat"),
                    row.names=FALSE,col.names=FALSE,quote=FALSE)
      }
  write.table((paste0("Rep",1:1000,".dat")),
              paste0("C:/LCM/ConditionSet",j,"/Rep.dat"),
              row.names=FALSE,col.names=FALSE,quote=FALSE)

write.table((paste0(
    "TITLE: 
LGCM_SR fit to ",DataList[[j]]$GenMod[1],"  ",DataList[[j]]$Cell[1],
    
    "

DATA:
    FILE = C:/LCM/ConditionSet",j,"/Rep.dat;
  TYPE = MONTECARLO;
  
  VARIABLE:
    names=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;
  usevar=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;
  analysis: estimator=ml;
  model:
    !Set Intercept and Slope
  Xint Xslp | x1@0 x2@1 x3 x4 x5 x6; 
  Yint Yslp | y1@0 y2@1 y3 y4 y5 y6; 
  
  !set measurement error in observed vars to 0
  x1-x6@0; 
  y1-y6@0;
  
  !create residuals to structure
  Xres1 by x1@1; 
  Xres2 by x2@1; 
  Xres3 by x3@1; 
  Xres4 by x4@1;
  Xres5 by x5@1;
  Xres6 by x6@1;
  
  Yres1 by y1@1;
  Yres2 by y2@1;
  Yres3 by y3@1;
  Yres4 by y4@1;
  Yres5 by y5@1;
  Yres6 by y6@1;
  
  !AR structured resids
  Xres2-Xres6 pon Xres1-Xres5 (1);
  Yres2-Yres6 pon Yres1-Yres5 (2);
  
  !CL structured resids
  Xres2-Xres6 pon Yres1-Yres5 (3); 
  Yres2-Yres6 pon Xres1-Xres5 (4);
  
  !Innovation structured resids
  Xres1-Xres6 pwith Yres1-Yres6 (5);
  
  
  !Zero out Between & Within-level Exogenous Covariances
  Xint WITH Xres1@0 Yres1@0;
  Xslp WITH Xres1@0 Yres1@0;
  Yint WITH Xres1@0 Yres1@0;
  Yslp WITH Xres1@0 Yres1@0;
  
  output: 
    sampstat stdyx;
  "
)),paste0("C:/LCM/ConditionSet",j,"/LGCM_SR.inp"),
    row.names=FALSE,col.names=FALSE,quote=FALSE)

write.table((paste0(
  "TITLE: 
GCLM fit to ",DataList[[j]]$GenMod[1],"  ",DataList[[j]]$Cell[1],
  
  "

DATA:
    FILE = C:/LCM/ConditionSet",j,"/Rep.dat;
  TYPE = MONTECARLO;
  
 VARIABLE:
   names=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6 ;
   usevar= x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;
analysis: estimator=ml;
model:
!Set unit effect
 EtaX BY x6 x5 x4 x3 x2 x1;
 EtaY BY y6 y5 y4 y3 y2 y1; 
 
!set measurement error in observed vars to 0
 x1-x6@0; 
 y1-y6@0;
 
 !create innovations
 Xres1 by x1; 
 Xres2 by x2; 
 Xres3 by x3; 
 Xres4 by x4;
 Xres5 by x5;
 Xres6 by x6;
 
 Yres1 by y1;
 Yres2 by y2;
 Yres3 by y3;
 Yres4 by y4;
 Yres5 by y5;
 Yres6 by y6;

 !AR 
x2-x6 pon x1-x5 (1);
y2-y6 pon y1-y5 (2);

!CL 
x2-x6 pon y1-y5 (3);
y2-y6 pon x1-x5 (4);


!Zero out Between & Within-level Exogenous Covariances
Xres1-Xres6 WITH EtaX@0 EtaY@0;
Yres1-Yres6 WITH EtaX@0 EtaY@0;
Xres1-Xres6 WITH Xres1-Xres6@0;
Yres1-Yres6 WITH Yres1-Yres6@0;
Xres1-Xres6 WITH Yres1-Yres6@0;

!comovements
Xres1-Xres6 PWITH Yres1-Yres6;

output: 
sampstat stdyx;
 "
)),paste0("C:/LCM/ConditionSet",j,"/GCLM.inp"),
  row.names=FALSE,col.names=FALSE,quote=FALSE)

write.table((paste0(
  "TITLE: 
ALT fit to ",DataList[[j]]$GenMod[1],"  ",DataList[[j]]$Cell[1],
  
  "

DATA:
    FILE = C:/LCM/ConditionSet",j,"/Rep.dat;
  TYPE = MONTECARLO;
  
VARIABLE:
   names=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;

   usevar= x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;

analysis: estimator=ml;

model:
!Set Intercept and Slope
Xb BY x2@1 x3 x4 x5 x6;
Yb BY y2@1 y3 y4 y5 y6;
 
 !Means
 [Xb Yb];
 [x2-x6@0 y2-y6@0];
 
   !AR structured resids
  x2-x6 pon x1-x5 (1);
  y2-y6 pon y1-y5 (2);

  !CL structured resids
   x2-x6 pon y1-y5 (3);
   y2-y6 pon x1-x5 (4);

  !covariances
  x1 with y1;
  Xb with Yb;
  Xb with x1 y1;
  Yb with x1 y1;
  x2-x6 pwith y2-y6 ;

output: 
sampstat stdyx;
"
)),paste0("C:/LCM/ConditionSet",j,"/ALT.inp"),
row.names=FALSE,col.names=FALSE,quote=FALSE)

write.table((paste0(
  "TITLE: 
RI_CLPM fit to ",DataList[[j]]$GenMod[1],"  ",DataList[[j]]$Cell[1],
  
  "

DATA:
    FILE = C:/LCM/ConditionSet",j,"/Rep.dat;
  TYPE = MONTECARLO;
  
  VARIABLE:
    names=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;
  usevar=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;
  analysis: estimator=ml;
  model:
    !Set Intercept and Slope
  Xint Xslp | x1@0 x2@1 x3 x4 x5 x6; 
  Yint Yslp | y1@0 y2@1 y3 y4 y5 y6; 
  Xslp@0; 
  Yslp@0;
  
  !set measurement error in observed vars to 0
  x1-x6@0; 
  y1-y6@0;
  
  !create residuals to structure
  Xres1 by x1@1; 
  Xres2 by x2@1; 
  Xres3 by x3@1; 
  Xres4 by x4@1;
  Xres5 by x5@1;
  Xres6 by x6@1;
  
  Yres1 by y1@1;
  Yres2 by y2@1;
  Yres3 by y3@1;
  Yres4 by y4@1;
  Yres5 by y5@1;
  Yres6 by y6@1;
  
  !AR structured resids
  Xres2-Xres6 pon Xres1-Xres5 (1);
  Yres2-Yres6 pon Yres1-Yres5 (2);
  
  !CL structured resids
  Xres2-Xres6 pon Yres1-Yres5 (3); 
  Yres2-Yres6 pon Xres1-Xres5 (4);
  
  !Innovation structured resids
  Xres1-Xres6 pwith Yres1-Yres6 (5);
  
  
  !Zero out Between & Within-level Exogenous Covariances
  Xint WITH Xres1@0 Yres1@0;
  Xslp WITH Xres1@0 Yres1@0 Xint@0 Yint@0 Yslp@0;
  Yint WITH Xres1@0 Yres1@0;
  Yslp WITH Xres1@0 Yres1@0 Yint@0 Xint@0;
  
  output: 
    sampstat stdyx;
  ")),paste0("C:/LCM/ConditionSet",j,"/RI_CLPM.inp"),
  row.names=FALSE,col.names=FALSE,quote=FALSE)


write.table((paste0(
  "TITLE: 
GCLM Mean Stationary fit to ",DataList[[j]]$GenMod[1],"  ",DataList[[j]]$Cell[1],
  
  "

DATA:
    FILE = C:/LCM/ConditionSet",j,"/Rep.dat;
  TYPE = MONTECARLO;
  
 VARIABLE:
   names=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6 ;
   usevar= x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;
analysis: estimator=ml;
model:
!Set unit effect
 EtaX BY x1@1 x2@1 x3@1 x4@1 x5@1 x6;
 EtaY BY y1@1 y2@1 y3@1 y4@1 y5@1 y6; 
 
!set measurement error in observed vars to 0
 x1-x6@0; 
 y1-y6@0;
 
 !create innovations
 Xres1 by x1; 
 Xres2 by x2; 
 Xres3 by x3; 
 Xres4 by x4;
 Xres5 by x5;
 Xres6 by x6;
 
 Yres1 by y1;
 Yres2 by y2;
 Yres3 by y3;
 Yres4 by y4;
 Yres5 by y5;
 Yres6 by y6;

 !AR 
x2-x6 pon x1-x5 (1);
y2-y6 pon y1-y5 (2);

!CL 
x2-x6 pon y1-y5 (3);
y2-y6 pon x1-x5 (4);


!Zero out Between & Within-level Exogenous Covariances
Xres1-Xres6 WITH EtaX@0 EtaY@0;
Yres1-Yres6 WITH EtaX@0 EtaY@0;
Xres1-Xres6 WITH Xres1-Xres6@0;
Yres1-Yres6 WITH Yres1-Yres6@0;
Xres1-Xres6 WITH Yres1-Yres6@0;
EtaX WITH EtaY;

!comovements
Xres1-Xres6 PWITH Yres1-Yres6;

output: 
sampstat stdyx;
 "
)),paste0("C:/LCM/ConditionSet",j,"/GCLM_MeanStationary.inp"),
row.names=FALSE,col.names=FALSE,quote=FALSE)

write.table((paste0(
  "TITLE: 
ALT_noSlpVar fit to ",DataList[[j]]$GenMod[1],"  ",DataList[[j]]$Cell[1],
  
  "

DATA:
    FILE = C:/LCM/ConditionSet",j,"/Rep.dat;
  TYPE = MONTECARLO;
  
VARIABLE:
   names=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;

   usevar= x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;

analysis: estimator=ml;

model:
!Set Intercept and Slope
Xb BY x2@1 x3 x4 x5 x6;
Yb BY y2@1 y3 y4 y5 y6;
Xb@0; Yb@0; 
 !Means
 [Xb Yb];
 [x2-x6@0 y2-y6@0];
 
   !AR structured resids
  x2-x6 pon x1-x5 (1);
  y2-y6 pon y1-y5 (2);

  !CL structured resids
   x2-x6 pon y1-y5 (3);
   y2-y6 pon x1-x5 (4);

  !covariances
  x1 with y1;
  Xb with Yb@0;
  Xb with x1@0 y1@0;
  Yb with x1@0 y1@0;
  x2-x6 pwith y2-y6 ;

output: 
sampstat stdyx;
"
)),paste0("C:/LCM/ConditionSet",j,"/ALT_noSlpVar.inp"),
row.names=FALSE,col.names=FALSE,quote=FALSE)
  
write.table((paste0(
  "TITLE: 
CLPM fit to ",DataList[[j]]$GenMod[1],"  ",DataList[[j]]$Cell[1],
  
  "

DATA:
    FILE = C:/LCM/ConditionSet",j,"/Rep.dat;
  TYPE = MONTECARLO;
  
VARIABLE:
   names=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;

   usevar= x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;

analysis: estimator=ml;

model:
   !AR structured resids
  x2-x6 pon x1-x5 (1);
  y2-y6 pon y1-y5 (2);

  !CL structured resids
   x2-x6 pon y1-y5 (3);
   y2-y6 pon x1-x5 (4);

  !covariances
x1-x6 pwith y1-y6 ;

output: 
sampstat stdyx;
"
)),paste0("C:/LCM/ConditionSet",j,"/CLPM.inp"),
row.names=FALSE,col.names=FALSE,quote=FALSE)

write.table((paste0(
  "TITLE: 
LGCM fit to ",DataList[[j]]$GenMod[1],"  ",DataList[[j]]$Cell[1],
  
  "

DATA:
    FILE = C:/LCM/ConditionSet",j,"/Rep.dat;
  TYPE = MONTECARLO;
  
  VARIABLE:
    names=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;
  usevar=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;
  analysis: estimator=ml;
  model:
    !Set Intercept and Slope
  Xint Xslp | x1@0 x2@1 x3 x4 x5 x6; 
  Yint Yslp | y1@0 y2@1 y3 y4 y5 y6; 
  
  output: 
    sampstat stdyx;
  "
  )),paste0("C:/LCM/ConditionSet",j,"/LGCM.inp"),
  row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
toc()
##################################################################
############Automate reading of batches from MPlus files##########
##################################################################
tic()
runModels("C:/LCM/", recursive=TRUE)
toc()
#################################################################
###########Bring in MPlus out file estimates into R##############
#################and prepare for analysis########################
#################################################################

tbl<-
cond%>%	
mutate(GenMod="LCM")%>%
	mutate(LSlpCov=factor(SxSyR,levels=c(0.1,0.5),labels=c("LoCov","HiCov")))%>%
	mutate(LVSx=factor(VSx,levels=c(3.44,13.76),labels=c("LoVSx","HiVSx")))%>%
	mutate(LVSy=factor(VSy,levels=c(4.22,16.88),labels=c("LoVSy","HiVSy")))%>%
	mutate(Cell=paste(LSlpCov,LVSx,LVSy))%>%
  group_by(Cell)
DataList<-group_split(tbl)

tic()
for (m in 1:8){

###ALT###
read_alt<-readModels(paste0("C:/LCM/ConditionSet",m,"/alt.out"))
alt_params<-as.tibble(read_alt$parameters$unstandardized)
alt_ARCLest<-alt_params[c(11:12,21:22),c(1:2,4:6)]
alt_ARCLest<-alt_ARCLest%>%
            unite(Param,"paramHeader","param")%>%
            mutate(Param=factor(Param,
                                levels=c("X2.ON_X1","X2.ON_Y1","Y2.ON_Y1","Y2.ON_X1"),
                                labels=c("ARx","CLxy","ARy","CLyx")))%>%
            unite(pooledparams,"average","population_sd","average_se")%>%
          spread(Param,pooledparams)%>%
          separate(ARx,sep="_",into=c("ARx_est","ARx_sd","ARx_se"),convert=T)%>%
          separate(CLxy,sep="_",into=c("CLxy_est","CLxy_sd","CLxy_se"),convert=T)%>%
          separate(ARy,sep="_",into=c("ARy_est","ARy_sd","ARy_se"),convert=T)%>%
          separate(CLyx,sep="_",into=c("CLyx_est","CLyx_sd","CLyx_se"),convert=T)%>%
          mutate(FitMod="alt")%>%
          mutate(GenMod="LCM")
##do similar for fit stats (RMSEA,CFI,SRMR,BIC)
ifelse(read_alt$summaries[1,10]>0,
alt_FitInd<-as.tibble(read_alt$summaries[1,c(17,21,33,36,39)]),
alt_FitInd<-tibble(LL_NumComputations=0,
			    CFI_Mean="NA", aBIC_Mean="NA",
			    RMSEA_Mean="NA",SRMR_Mean="NA")
)

#cbind fit and ARCL
alt_FitInd_ARCL<-cbind(alt_FitInd,alt_ARCLest)

#cbind (fit&ARCL) with conditons & parameter values
alt_ModSum<-cbind(alt_FitInd_ARCL,
              DataList[[m]][1,c(4:6,9:12)])

#calculate Relative Bias for est & se to produce
#final 1 X 35 object per condition/loop
alt_ModSum<-
  alt_ModSum %>%
  mutate(bias_est_ARx=(ARx_est-ARx)/ARx,
        bias_est_ARy=(ARy_est-ARy)/ARy,
        bias_se_ARx=(ARx_se-ARx_sd)/ARx_sd,
        bias_se_ARy=(ARy_se-ARy_sd)/ARy_sd,
        bias_est_CLxy=(CLxy_est-CLxy)/CLxy,
        bias_est_CLyx=(CLyx_est-0.05)/0.05,
        bias_se_CLxy=(CLxy_se-CLxy_sd)/CLxy_sd,
        bias_se_CLyx=(CLyx_se-CLyx_sd)/CLyx_sd,
	  ConvRate=LL_NumComputations/10)

###LGCM-SR###
read_lgcm_sr<-readModels(paste0("C:/LCM/ConditionSet",m,"/lgcm_sr.out"))
lgcm_sr_params<-as.tibble(read_lgcm_sr$parameters$unstandardized)
lgcm_sr_ARCLest<-lgcm_sr_params[c(37:38,47:48),c(1:2,4:6)]
lgcm_sr_ARCLest<-lgcm_sr_ARCLest%>%
  unite(Param,"paramHeader","param")%>%
  mutate(Param=factor(Param,
                      levels=c("XRES2.ON_XRES1","XRES2.ON_YRES1","YRES2.ON_YRES1","YRES2.ON_XRES1"),
                      labels=c("ARx","CLxy","ARy","CLyx")))%>%
  unite(pooledparams,"average","population_sd","average_se")%>%
  spread(Param,pooledparams)%>%
  separate(ARx,sep="_",into=c("ARx_est","ARx_sd","ARx_se"),convert=T)%>%
  separate(CLxy,sep="_",into=c("CLxy_est","CLxy_sd","CLxy_se"),convert=T)%>%
  separate(ARy,sep="_",into=c("ARy_est","ARy_sd","ARy_se"),convert=T)%>%
  separate(CLyx,sep="_",into=c("CLyx_est","CLyx_sd","CLyx_se"),convert=T)%>%
  mutate(FitMod="lgcm_sr")%>%
  mutate(GenMod="LCM")
##do similar for fit stats (RMSEA,CFI,SRMR,BIC)
ifelse(read_lgcm_sr$summaries[1,10]>0,
lgcm_sr_FitInd<-as.tibble(read_lgcm_sr$summaries[1,c(17,21,33,36,39)]),
lgcm_sr_FitInd<-tibble(LL_NumComputations=0,
			    CFI_Mean="NA", aBIC_Mean="NA",
			    RMSEA_Mean="NA",SRMR_Mean="NA")
)

#cbind fit and ARCL
lgcm_sr_FitInd_ARCL<-cbind(lgcm_sr_FitInd,lgcm_sr_ARCLest)

#cbind (fit&ARCL) with conditons & parameter values
lgcm_sr_ModSum<-cbind(lgcm_sr_FitInd_ARCL,
                      DataList[[m]][1,c(4:6,9:12)])

#calculate Relative Bias for est & se to produce
#final 1 X 35 object per condition/loop
lgcm_sr_ModSum<-
  lgcm_sr_ModSum %>%
  mutate(bias_est_ARx=(ARx_est-ARx)/ARx,
         bias_est_ARy=(ARy_est-ARy)/ARy,
         bias_se_ARx=(ARx_se-ARx_sd)/ARx_sd,
         bias_se_ARy=(ARy_se-ARy_sd)/ARy_sd,
         bias_est_CLxy=(CLxy_est-CLxy)/CLxy,
         bias_est_CLyx=(CLyx_est-0.05)/0.05,
         bias_se_CLxy=(CLxy_se-CLxy_sd)/CLxy_sd,
         bias_se_CLyx=(CLyx_se-CLyx_sd)/CLyx_sd,
	   ConvRate=LL_NumComputations/10)

###GCLM###
read_gclm<-readModels(paste0("C:/LCM/ConditionSet",m,"/gclm.out"))
gclm_params<-as.tibble(read_gclm$parameters$unstandardized)
gclm_ARCLest<-gclm_params[c(25:26,35:36),c(1:2,4:6)]
gclm_ARCLest<-gclm_ARCLest%>%
  unite(Param,"paramHeader","param")%>%
  mutate(Param=factor(Param,
                      levels=c("X2.ON_X1","X2.ON_Y1","Y2.ON_Y1","Y2.ON_X1"),
                      labels=c("ARx","CLxy","ARy","CLyx")))%>%
  unite(pooledparams,"average","population_sd","average_se")%>%
  spread(Param,pooledparams)%>%
  separate(ARx,sep="_",into=c("ARx_est","ARx_sd","ARx_se"),convert=T)%>%
  separate(CLxy,sep="_",into=c("CLxy_est","CLxy_sd","CLxy_se"),convert=T)%>%
  separate(ARy,sep="_",into=c("ARy_est","ARy_sd","ARy_se"),convert=T)%>%
  separate(CLyx,sep="_",into=c("CLyx_est","CLyx_sd","CLyx_se"),convert=T)%>%
  mutate(FitMod="gclm")%>%
  mutate(GenMod="LCM")
##do similar for fit stats (RMSEA,CFI,SRMR,BIC)
ifelse(read_gclm$summaries[1,10]>0,
gclm_FitInd<-as.tibble(read_gclm$summaries[1,c(17,21,33,36,39)]),
gclm_FitInd<-tibble(LL_NumComputations=0,
			    CFI_Mean="NA", aBIC_Mean="NA",
			    RMSEA_Mean="NA",SRMR_Mean="NA")
)

#cbind fit and ARCL
gclm_FitInd_ARCL<-cbind(gclm_FitInd,gclm_ARCLest)

#cbind (fit&ARCL) with conditons & parameter values
gclm_ModSum<-cbind(gclm_FitInd_ARCL,
                   DataList[[m]][1,c(4:6,9:12)])

#calculate Relative Bias for est & se to produce
#final 1 X 35 object per condition/loop
gclm_ModSum<-
  gclm_ModSum %>%
  mutate(bias_est_ARx=(ARx_est-ARx)/ARx,
         bias_est_ARy=(ARy_est-ARy)/ARy,
         bias_se_ARx=(ARx_se-ARx_sd)/ARx_sd,
         bias_se_ARy=(ARy_se-ARy_sd)/ARy_sd,
         bias_est_CLxy=(CLxy_est-CLxy)/CLxy,
         bias_est_CLyx=(CLyx_est-0.05)/0.05,
         bias_se_CLxy=(CLxy_se-CLxy_sd)/CLxy_sd,
         bias_se_CLyx=(CLyx_se-CLyx_sd)/CLyx_sd,
	   ConvRate=LL_NumComputations/10)

read_alt_noslpvar<-readModels(paste0("C:/LCM/ConditionSet",m,"/alt_noslpvar.out"))
alt_noslpvar_params<-as.tibble(read_alt_noslpvar$parameters$unstandardized)
alt_noslpvar_ARCLest<-alt_noslpvar_params[c(11:12,21:22),c(1:2,4:6)]
alt_noslpvar_ARCLest<-alt_noslpvar_ARCLest%>%
            unite(Param,"paramHeader","param")%>%
            mutate(Param=factor(Param,
                                levels=c("X2.ON_X1","X2.ON_Y1","Y2.ON_Y1","Y2.ON_X1"),
                                labels=c("ARx","CLxy","ARy","CLyx")))%>%
            unite(pooledparams,"average","population_sd","average_se")%>%
          spread(Param,pooledparams)%>%
          separate(ARx,sep="_",into=c("ARx_est","ARx_sd","ARx_se"),convert=T)%>%
          separate(CLxy,sep="_",into=c("CLxy_est","CLxy_sd","CLxy_se"),convert=T)%>%
          separate(ARy,sep="_",into=c("ARy_est","ARy_sd","ARy_se"),convert=T)%>%
          separate(CLyx,sep="_",into=c("CLyx_est","CLyx_sd","CLyx_se"),convert=T)%>%
          mutate(FitMod="alt_noslpvar")%>%
          mutate(GenMod="LCM")
##do similar for fit stats (RMSEA,CFI,SRMR,BIC)
ifelse(read_alt_noslpvar$summaries[1,10]>0,
alt_noslpvar_FitInd<-as.tibble(read_alt_noslpvar$summaries[1,c(17,21,33,36,39)]),
alt_noslpvar_FitInd<-tibble(LL_NumComputations=0,
			    CFI_Mean="NA", aBIC_Mean="NA",
			    RMSEA_Mean="NA",SRMR_Mean="NA")
)

#cbind fit and ARCL
alt_noslpvar_FitInd_ARCL<-cbind(alt_noslpvar_FitInd,alt_noslpvar_ARCLest)

#cbind (fit&ARCL) with conditons & parameter values
alt_noslpvar_ModSum<-cbind(alt_noslpvar_FitInd_ARCL,
              DataList[[m]][1,c(4:6,9:12)])

#calculate Relative Bias for est & se to produce
#final 1 X 35 object per condition/loop
alt_noslpvar_ModSum<-
  alt_noslpvar_ModSum %>%
  mutate(bias_est_ARx=(ARx_est-ARx)/ARx,
        bias_est_ARy=(ARy_est-ARy)/ARy,
        bias_se_ARx=(ARx_se-ARx_sd)/ARx_sd,
        bias_se_ARy=(ARy_se-ARy_sd)/ARy_sd,
        bias_est_CLxy=(CLxy_est-CLxy)/CLxy,
        bias_est_CLyx=(CLyx_est-0.05)/0.05,
        bias_se_CLxy=(CLxy_se-CLxy_sd)/CLxy_sd,
        bias_se_CLyx=(CLyx_se-CLyx_sd)/CLyx_sd,
	  ConvRate=LL_NumComputations/10)

###RI-CLPM###
read_ri_clpm<-readModels(paste0("C:/LCM/ConditionSet",m,"/ri_clpm.out"))
ri_clpm_params<-as.tibble(read_ri_clpm$parameters$unstandardized)
ri_clpm_ARCLest<-ri_clpm_params[c(37:38,47:48),c(1:2,4:6)]
ri_clpm_ARCLest<-ri_clpm_ARCLest%>%
  unite(Param,"paramHeader","param")%>%
  mutate(Param=factor(Param,
                      levels=c("XRES2.ON_XRES1","XRES2.ON_YRES1","YRES2.ON_YRES1","YRES2.ON_XRES1"),
                      labels=c("ARx","CLxy","ARy","CLyx")))%>%
  unite(pooledparams,"average","population_sd","average_se")%>%
  spread(Param,pooledparams)%>%
  separate(ARx,sep="_",into=c("ARx_est","ARx_sd","ARx_se"),convert=T)%>%
  separate(CLxy,sep="_",into=c("CLxy_est","CLxy_sd","CLxy_se"),convert=T)%>%
  separate(ARy,sep="_",into=c("ARy_est","ARy_sd","ARy_se"),convert=T)%>%
  separate(CLyx,sep="_",into=c("CLyx_est","CLyx_sd","CLyx_se"),convert=T)%>%
  mutate(FitMod="ri_clpm")%>%
  mutate(GenMod="LCM")
##do similar for fit stats (RMSEA,CFI,SRMR,BIC)
ifelse(read_ri_clpm$summaries[1,10]>0,
ri_clpm_FitInd<-as.tibble(read_ri_clpm$summaries[1,c(17,21,33,36,39)]),
ri_clpm_FitInd<-tibble(LL_NumComputations=0,
			    CFI_Mean="NA", aBIC_Mean="NA",
			    RMSEA_Mean="NA",SRMR_Mean="NA")
)

#cbind fit and ARCL
ri_clpm_FitInd_ARCL<-cbind(ri_clpm_FitInd,ri_clpm_ARCLest)

#cbind (fit&ARCL) with conditons & parameter values
ri_clpm_ModSum<-cbind(ri_clpm_FitInd_ARCL,
                   DataList[[m]][1,c(4:6,9:12)])

#calculate Relative Bias for est & se to produce
#final 1 X 35 object per condition/loop
ri_clpm_ModSum<-
  ri_clpm_ModSum %>%
  mutate(bias_est_ARx=(ARx_est-ARx)/ARx,
         bias_est_ARy=(ARy_est-ARy)/ARy,
         bias_se_ARx=(ARx_se-ARx_sd)/ARx_sd,
         bias_se_ARy=(ARy_se-ARy_sd)/ARy_sd,
         bias_est_CLxy=(CLxy_est-CLxy)/CLxy,
         bias_est_CLyx=(CLyx_est-0.05)/0.05,
         bias_se_CLxy=(CLxy_se-CLxy_sd)/CLxy_sd,
         bias_se_CLyx=(CLyx_se-CLyx_sd)/CLyx_sd,
	   ConvRate=LL_NumComputations/10)

###GCLM Mean Stationary###
read_gclm_meanstationary<-readModels(paste0("C:/LCM/ConditionSet",m,"/gclm_meanstationary.out"))
gclm_meanstationary_params<-as.tibble(read_gclm_meanstationary$parameters$unstandardized)
gclm_meanstationary_ARCLest<-gclm_meanstationary_params[c(25:26,35:36),c(1:2,4:6)]
gclm_meanstationary_ARCLest<-gclm_meanstationary_ARCLest%>%
  unite(Param,"paramHeader","param")%>%
  mutate(Param=factor(Param,
                      levels=c("X2.ON_X1","X2.ON_Y1","Y2.ON_Y1","Y2.ON_X1"),
                      labels=c("ARx","CLxy","ARy","CLyx")))%>%
  unite(pooledparams,"average","population_sd","average_se")%>%
  spread(Param,pooledparams)%>%
  separate(ARx,sep="_",into=c("ARx_est","ARx_sd","ARx_se"),convert=T)%>%
  separate(CLxy,sep="_",into=c("CLxy_est","CLxy_sd","CLxy_se"),convert=T)%>%
  separate(ARy,sep="_",into=c("ARy_est","ARy_sd","ARy_se"),convert=T)%>%
  separate(CLyx,sep="_",into=c("CLyx_est","CLyx_sd","CLyx_se"),convert=T)%>%
  mutate(FitMod="gclm_meanstationary")%>%
  mutate(GenMod="LCM")
##do similar for fit stats (RMSEA,CFI,SRMR,BIC)
ifelse(read_gclm_meanstationary$summaries[1,10]>0,
gclm_meanstationary_FitInd<-as.tibble(read_gclm_meanstationary$summaries[1,c(17,21,33,36,39)]),
gclm_meanstationary_FitInd<-tibble(LL_NumComputations=0,
			    CFI_Mean="NA", aBIC_Mean="NA",
			    RMSEA_Mean="NA",SRMR_Mean="NA")
)

#cbind fit and ARCL
gclm_meanstationary_FitInd_ARCL<-cbind(gclm_meanstationary_FitInd,gclm_meanstationary_ARCLest)

#cbind (fit&ARCL) with conditons & parameter values
gclm_meanstationary_ModSum<-cbind(gclm_meanstationary_FitInd_ARCL,
                   DataList[[m]][1,c(4:6,9:12)])

#calculate Relative Bias for est & se to produce
#final 1 X 35 object per condition/loop
gclm_meanstationary_ModSum<-
  gclm_meanstationary_ModSum %>%
  mutate(bias_est_ARx=(ARx_est-ARx)/ARx,
         bias_est_ARy=(ARy_est-ARy)/ARy,
         bias_se_ARx=(ARx_se-ARx_sd)/ARx_sd,
         bias_se_ARy=(ARy_se-ARy_sd)/ARy_sd,
         bias_est_CLxy=(CLxy_est-CLxy)/CLxy,
         bias_est_CLyx=(CLyx_est-0.05)/0.05,
         bias_se_CLxy=(CLxy_se-CLxy_sd)/CLxy_sd,
         bias_se_CLyx=(CLyx_se-CLyx_sd)/CLyx_sd,
	   ConvRate=LL_NumComputations/10)

###CLPM###
read_clpm<-readModels(paste0("C:/LCM/ConditionSet",m,"/clpm.out"))
clpm_params<-as.tibble(read_clpm$parameters$unstandardized)
clpm_ARCLest<-clpm_params[c(1:2,11:12),c(1:2,4:6)]
clpm_ARCLest<-clpm_ARCLest%>%
  unite(Param,"paramHeader","param")%>%
  mutate(Param=factor(Param,
                      levels=c("X2.ON_X1","X2.ON_Y1","Y2.ON_Y1","Y2.ON_X1"),
                      labels=c("ARx","CLxy","ARy","CLyx")))%>%
  unite(pooledparams,"average","population_sd","average_se")%>%
  spread(Param,pooledparams)%>%
  separate(ARx,sep="_",into=c("ARx_est","ARx_sd","ARx_se"),convert=T)%>%
  separate(CLxy,sep="_",into=c("CLxy_est","CLxy_sd","CLxy_se"),convert=T)%>%
  separate(ARy,sep="_",into=c("ARy_est","ARy_sd","ARy_se"),convert=T)%>%
  separate(CLyx,sep="_",into=c("CLyx_est","CLyx_sd","CLyx_se"),convert=T)%>%
  mutate(FitMod="clpm")%>%
  mutate(GenMod="LCM")
##do similar for fit stats (RMSEA,CFI,SRMR,BIC)
ifelse(read_clpm$summaries[1,10]>0,
clpm_FitInd<-as.tibble(read_clpm$summaries[1,c(17,21,33,36,39)]),
clpm_FitInd<-tibble(LL_NumComputations=0,
			    CFI_Mean="NA", aBIC_Mean="NA",
			    RMSEA_Mean="NA",SRMR_Mean="NA")
)


#cbind fit and ARCL
clpm_FitInd_ARCL<-cbind(clpm_FitInd,clpm_ARCLest)

#cbind (fit&ARCL) with conditons & parameter values
clpm_ModSum<-cbind(clpm_FitInd_ARCL,
                   DataList[[m]][1,c(4:6,9:12)])

#calculate Relative Bias for est & se to produce
#final 1 X 35 object per condition/loop
clpm_ModSum<-
  clpm_ModSum %>%
  mutate(bias_est_ARx=(ARx_est-ARx)/ARx,
         bias_est_ARy=(ARy_est-ARy)/ARy,
         bias_se_ARx=(ARx_se-ARx_sd)/ARx_sd,
         bias_se_ARy=(ARy_se-ARy_sd)/ARy_sd,
         bias_est_CLxy=(CLxy_est-CLxy)/CLxy,
         bias_est_CLyx=(CLyx_est-0.05)/0.05,
         bias_se_CLxy=(CLxy_se-CLxy_sd)/CLxy_sd,
         bias_se_CLyx=(CLyx_se-CLyx_sd)/CLyx_sd,
	   ConvRate=LL_NumComputations/10)

#rbind all 7 model entries, resulting in  7 X 35 object
GeneralModSum<-rbind(
		lgcm_sr_ModSum,alt_ModSum,gclm_ModSum,
		ri_clpm_ModSum, alt_noslpvar_ModSum,gclm_meanstationary_ModSum,
		clpm_ModSum)		
write.table(GeneralModSum,
			paste0("C:/LCM/ConditionSet",m,"/Set4Analysis.dat"),
			sep="\t",row.names=FALSE)
}
toc()

#append tables from each loop for object with 7*8=56 rows X 37 columns
#use pmap to read external files into a tibble
file<-paste0("C:/LCM/ConditionSet",1:8,"/Set4Analysis.dat")
sep<-rep("\t",8)
header<-rep(TRUE,8)
args<-tibble(file,sep,header)
MastLCM<-
args%>%
	pmap(read.table)
#append all of the datasets together
MastDat<-do.call(rbind,MastLCM)

write.table(MastDat,"C:/LCM/SimResult_LCM.dat",
		sep="\t",row.names=FALSE)

#this object can then be rbinded to other generating model products

othermods<-read.table("C:/Users/pws5/Desktop/GCLM_LGCMSR_CLPM_res.dat",
				sep="\t",header=TRUE)

AnalyzeLCM<-read.table("C:/LCM/SimResult_LCM.dat",sep="\t",header=TRUE)
AnalyzeLCM<-AnalyzeLCM%>%
	mutate(LARx="NA")%>%
	mutate(LARy="NA")%>%
	mutate(LCL="NA")%>%
	rename(LCov=LSlpCov)%>%
	select(-Cell)

AllModRes<-rbind(othermods,AnalyzeLCM)
write.table(AllModRes,"C:/Users/pws5/desktop/all_mod_results.dat",sep="\t",row.names=FALSE)
bringin<-read.table("C:/Users/pws5/desktop/all_mod_results.dat",sep="\t",header=TRUE)

##final object should be 1008 by 37

###############################################################
########################Analyze MCMC results###################
###############################################################

