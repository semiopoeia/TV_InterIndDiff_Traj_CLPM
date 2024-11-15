			####GCLM#####
	###Generating###
GCLM_gen<-  
'
#impulse variances
Vx1~~33.75*Vx1
Vx2~~30.13*Vx2
Vx3~~25.11*Vx3
Vx4~~24.80*Vx4
Vx5~~18.10*Vx5
Vx6~~23.23*Vx6

Vy1~~29.31*Vy1
Vy2~~27.22*Vy2
Vy3~~24.59*Vy3
Vy4~~25.86*Vy4
Vy5~~26.59*Vy5
Vy6~~25.37*Vy6

#unit effect variances
EtaX~~51.65*EtaX
EtaY~~75.03*EtaY

#unit covariance
EtaX~~xy_eta_cov*EtaY

#co-movements
Vx1~~6.742*Vy1
Vx2~~6.742*Vy1
Vx3~~6.742*Vy1
Vx4~~6.742*Vy1
Vx5~~6.742*Vy1
Vx6~~6.742*Vy1

#observed vars
x1~50.352*1+Lx1*EtaX+1*Vx1
y1~50.352*1+Ly1*EtaY+1*Vy1

x2~30.61*1+Lx2*EtaX+ARx21*x1+CLxy21*y1+1*Vx2
y2~31.45*1+Ly2*EtaY+ARy21*y1+CLyx21*x1+1*Vy2

x3~33.27*1+Lx3*EtaX+ARx32*x2+CLxy32*y2+1*Vx3
y3~34.06*1+Ly3*EtaY+ARy32*y2+CLyx32*x2+1*Vy3

x4~33.51*1+Lx4*EtaX+ARx43*x3+CLxy43*y3+1*Vx4
y4~36.50*1+Ly4*EtaY+ARy43*y3+CLyx43*x3+1*Vy4

x5~34.33*1+Lx5*EtaX+ARx54*x4+CLxy54*y4+1*Vx5
y5~37.47*1+Ly5*EtaY+ARy54*y4+CLyx54*x4+1*Vy5

x6~34.20*1+Lx6*EtaX+ARx65*x5+CLxy65*y5+1*Vx6
y6~36.73*1+Ly6*EtaY+ARy65*y5+CLyx65*x5+1*Vy6
'

#####LGCM-SR#######
  ###Generating######
LGCM_SR_gen<-
'
	#zero resid on observed vars
#x1~~0*x1
#x2~~0*x2
#x3~~0*x3
#x4~~0*x4
#x5~~0*x5
#x6~~0*x6
#y1~~0*y1
#y2~~0*y2
#y3~~0*y3
#y4~~0*y4
#y5~~0*y5
#y6~~0*y6

#Vx1=~1*x1
#Vx2=~1*x2
#Vx3=~1*x3
#Vx4=~1*x4
#Vx5=~1*x5
#Vx6=~1*x6
#Vy1=~1*y1
#Vy2=~1*y2
#Vy3=~1*y3
#Vy4=~1*y4
#Vy5=~1*y5
#Vy6=~1*y6

#growth curve params
#Ix=~1*x1+1*x2+1*x3+1*x4+1*x5+1*x6
#Sx=~0*x1+1*x2+2.3*x3+3.6*x4+4.3*x5+4.6*x6
#Iy=~1*y1+1*y2+1*y3+1*y4+1*y5+1*y6
#Sy=~0*y1+1*y2+2.2*y3+2.8*y4+3.5*y5+3.8*y6

	#residual variance/impulse values
		Vx1~~24.5*Vx1
		Vy1~~31*Vy1
		Vx2~~25.1*Vx2
		Vy2~~28.7*Vy2
		Vx3~~23.5*Vx3
		Vy3~~22.5*Vy3
		Vx4~~22.4*Vx4
		Vy4~~22.3*Vy4
		Vx5~~23*Vx5
		Vy5~~15.7*Vy5
		Vx6~~24*Vx6
		Vy6~~22.6*Vy6
		
	#residual covariance/co-movements
		Vx1~~5.12*Vy1
		Vx2~~5.12*Vy2
		Vx3~~5.12*Vy3
		Vx4~~5.12*Vy4
		Vx5~~5.12*Vy5
		Vx6~~5.12*Vy6
		
	#ARCL on SR	
		Vx2~ARx21*Vx1+CLxy21*Vy1
		Vy2~ARy21*Vy1+CLyx21*Vx1
		Vx3~ARx32*Vx2+CLxy32*Vy2
		Vy3~ARy32*Vy2+CLyx32*Vx2
		Vx4~ARx43*Vx3+CLxy43*Vy3
		Vy4~ARy43*Vy3+CLyx43*Vx3
		Vx5~ARx54*Vx4+CLxy54*Vy4
		Vy5~ARy54*Vy4+CLyx54*Vx4
		Vx6~ARx65*Vx5+CLxy65*Vy5
		Vy6~ARy65*Vy5+CLyx65*Vx5

  #Growth Curve Parameter Means
    Ix~50.3*1
    Sx~3.3*1
    Iy~50.4*1
    Sy~3.3*1
    
	#Growth Curve Parameter Var/Cov
		Ix~~79.3*Ix
		Sx~~x_slp_var*Sx
		Iy~~73.1*Iy
		Sy~~y_slp_var*Sy
		
		Ix~~XiXs_cov*Sx
		Iy~~YiYs_cov*Sy
		Ix~~XiYi_cov*Iy
		Sx~~xy_slp_cov*Sy
		Ix~~XiYs_cov*Sy
		Iy~~YiXs_cov*Sx
		Ix~~0*Vx1+0*Vy1
		Sx~~0*Vx1+0*Vy1
		Sy~~0*Vx1+0*Vy1
		Iy~~0*Vx1+0*Vy1
		
  #observed vars
  x1~0*1+1*Ix+0*Sx+1*Vx1
  y1~0*1+1*Iy+0*Sy+1*Vy1
  x2~0*1+1*Ix+1*Sx+1*Vx2
  y2~0*1+1*Iy+1*Sy+1*Vy2
  x3~0*1+1*Ix+2.3*Sx+1*Vx3
  y3~0*1+1*Iy+2.2*Sy+1*Vy3
  x4~0*1+1*Ix+3.6*Sx+1*Vx4
  y4~0*1+1*Iy+2.8*Sy+1*Vy4
  x5~0*1+1*Ix+4.3*Sx+1*Vx5
  y5~0*1+1*Iy+3.5*Sy+1*Vy5
  x6~0*1+1*Ix+4.6*Sx+1*Vx6
  y6~0*1+1*Iy+3.8*Sy+1*Vy6
	 '


ALT_gen<-
'

'

#########Fitting LGCM-SR slope var (Model E)#########
LGCM_SR_SlpVar<-
  '
#growth curve params
Ix=~1*x1+1*x2+1*x3+1*x4+1*x5+1*x6
Sx=~0*x1+1*x2+lx3*x3+lx4*x4+lx5*x5+lx6*x6
Iy=~1*y1+1*y2+1*y3+1*y4+1*y5+1*y6
Sy=~0*y1+1*y2+ly3*y3+ly4*y4+ly5*y5+ly6*y6

#zero out measurement errors on ov
x1~~0*x1
x2~~0*x2
x3~~0*x3
x4~~0*x4
x5~~0*x5
x6~~0*x6

y1~~0*y1
y2~~0*y2
y3~~0*y3
y4~~0*y4
y5~~0*y5
y6~~0*y6

#set resids as measurement
vx1=~1*x1      
vx2=~1*x2      
vx3=~1*x3      
vx4=~1*x4      
vx5=~1*x5      
vx6=~1*x6

vy1=~1*y1
vy2=~1*y2
vy3=~1*y3
vy4=~1*y4
vy5=~1*y5
vy6=~1*y6

#ARCL structured resids
vx2~arx*vx1+clxy*vy1
vx3~arx*vx2+clxy*vy2
vx4~arx*vx3+clxy*vy3
vx5~arx*vx4+clxy*vy4
vx6~arx*vx5+clxy*vy5

vy2~ary*vy1+clyx*vx1
vy3~ary*vy2+clyx*vx2
vy4~ary*vy3+clyx*vx3
vy5~ary*vy4+clyx*vx4
vy6~ary*vy5+clyx*vx5

#set covaraiances
vx1~~c*vy1          
vx2~~c*vy2          
vx3~~c*vy3          
vx4~~c*vy4          
vx5~~c*vy5          
vx6~~c*vy6

Ix~~Iy
Ix~~Sx
Ix~~Sy
Iy~~Sy
Iy~~Sx
Sx~~Sy

#estimate latent means
Ix~1
Iy~1
Sx~1
Sy~1

#constraints
vx1~~0*Ix
vx1~~0*Sx
vx1~~0*Iy
vx1~~0*Sy

vy1~~0*Ix
vy1~~0*Sx
vy1~~0*Iy
vy1~~0*Sy
'

#########Fitting LGCM-SR slope var (Model E)#########
LGCM_SR_NoSlpVar<-
  '
#growth curve params
Ix=~1*x1+1*x2+1*x3+1*x4+1*x5+1*x6
Sx=~0*x1+1*x2+lx3*x3+lx4*x4+lx5*x5+lx6*x6
Iy=~1*y1+1*y2+1*y3+1*y4+1*y5+1*y6
Sy=~0*y1+1*y2+ly3*y3+ly4*y4+ly5*y5+ly6*y6

#zero out measurement errors on ov
x1~~0*x1
x2~~0*x2
x3~~0*x3
x4~~0*x4
x5~~0*x5
x6~~0*x6

y1~~0*y1
y2~~0*y2
y3~~0*y3
y4~~0*y4
y5~~0*y5
y6~~0*y6

#set resids as measurement
vx1=~1*x1      
vx2=~1*x2      
vx3=~1*x3      
vx4=~1*x4      
vx5=~1*x5      
vx6=~1*x6

vy1=~1*y1
vy2=~1*y2
vy3=~1*y3
vy4=~1*y4
vy5=~1*y5
vy6=~1*y6

#ARCL structured resids
vx2~arx*vx1+clxy*vy1
vx3~arx*vx2+clxy*vy2
vx4~arx*vx3+clxy*vy3
vx5~arx*vx4+clxy*vy4
vx6~arx*vx5+clxy*vy5

vy2~ary*vy1+clyx*vx1
vy3~ary*vy2+clyx*vx2
vy4~ary*vy3+clyx*vx3
vy5~ary*vy4+clyx*vx4
vy6~ary*vy5+clyx*vx5

#set covaraiances
vx1~~vy1          
vx2~~vy2          
vx3~~vy3          
vx4~~vy4          
vx5~~vy5          
vx6~~vy6

Sx~~0*Sx
Sy~~0*Sy

Ix~~Iy
Ix~~0*Sx
Ix~~0*Sy
Iy~~0*Sy
Iy~~0*Sx
Sx~~0*Sy

#estimate latent means
Ix~1
Iy~1
Sx~1
Sy~1

#constraints
vx1~~0*Ix
vx1~~0*Sx
vx1~~0*Iy
vx1~~0*Sy

vy1~~0*Ix
vy1~~0*Sx
vy1~~0*Iy
vy1~~0*Sy

x1~0*1
y1~0*1
x2~0*1
y2~0*1
x3~0*1
y3~0*1
x4~0*1
y4~0*1
x5~0*1
y5~0*1
x6~0*1
y6~0*1
'


#####Fitting GCLM time-varying unit effect (Model C)#######
GCLM_FreeLoad<-
  '
  #unit effects
Sx=~x6+x5+x4+x3+x2+x1
Sy=~y6+y5+y4+y3+y2+y1

Sx~0*1
Sy~0*1

#zero out measurement errors
x1~~0*x1
x2~~0*x2
x3~~0*x3
x4~~0*x4
x5~~0*x5
x6~~0*x6

y1~~0*y1
y2~~0*y2
y3~~0*y3
y4~~0*y4
y5~~0*y5
y6~~0*y6

#set impulses
vx1=~x1
vx2=~x2
vx3=~x3
vx4=~x4
vx5=~x5
vx6=~x6

vy1=~y1
vy2=~y2
vy3=~y3
vy4=~y4
vy5=~y5
vy6=~y6

#ARCL
x2~arx*x1+clxy*y1
x3~arx*x2+clxy*y2
x4~arx*x3+clxy*y3
x5~arx*x4+clxy*y4
x6~arx*x5+clxy*y5

y2~ary*y1+clyx*x1
y3~ary*y2+clyx*x2
y4~ary*y3+clyx*x3
y5~ary*y4+clyx*x4
y6~ary*y5+clyx*x5

#constraints
vx1~~0*Sx+0*Sy
vy1~~0*Sx+0*Sy
vx2~~0*Sx+0*Sy
vy2~~0*Sx+0*Sy
vx3~~0*Sx+0*Sy
vy3~~0*Sx+0*Sy
vx4~~0*Sx+0*Sy
vy4~~0*Sx+0*Sy
vx5~~0*Sx+0*Sy
vy5~~0*Sx+0*Sy
vx6~~0*Sx+0*Sy
vy6~~0*Sx+0*Sy

vx1~~0*vx2+0*vy2+0*vx3+0*vy3+0*vx4+0*vy4+0*vx5+0*vy5+0*vx6+0*vy6
vx2~~0*vx3+0*vy1+0*vy3+0*vx4+0*vy4+0*vx5+0*vy5+0*vx6+0*vy6
vx3~~0*vx4+0*vx5+0*vx6+0*vy1+0*vy2+0*vy4+0*vy5+0*vy6
vx4~~0*vx5+0*vx6+0*vy1+0*vy2+0*vy3+0*vy5+0*vy6
vx5~~0*vx6+0*vy1+0*vy2+0*vy3+0*vy4+0*vy6
vx6~~0*vy1+0*vy2+0*vy3+0*vy4+0*vy5
vy1~~0*vy2+0*vy3+0*vy4+0*vy5+0*vy6
vy2~~0*vy3+0*vy4+0*vy5+0*vy6
vy3~~0*vy4+0*vy5+0*vy6
vy4~~0*vy5+0*vy6
vy5~~0*vy6

#set comovements
vx1~~vy1          
vx2~~vy2          
vx3~~vy3          
vx4~~vy4          
vx5~~vy5          
vx6~~vy6

Sx~~Sy
'

#######Fitting GCLM time-invariant unit effect (Model D)######
GCLM_FixLoad<-
  '  
#unit effects
Sx=~1*x6+1*x5+1*x4+1*x3+1*x2+x1
Sy=~1*y6+1*y5+1*y4+1*y3+1*y2+y1

Sx~0*1
Sy~0*1

#zero out measurement errors
x1~~0*x1
x2~~0*x2
x3~~0*x3
x4~~0*x4
x5~~0*x5
x6~~0*x6

y1~~0*y1
y2~~0*y2
y3~~0*y3
y4~~0*y4
y5~~0*y5
y6~~0*y6

#set impulses
vx1=~x1
vx2=~x2
vx3=~x3
vx4=~x4
vx5=~x5
vx6=~x6

vy1=~y1
vy2=~y2
vy3=~y3
vy4=~y4
vy5=~y5
vy6=~y6

#ARCL
x2~arx*x1+clxy*y1
x3~arx*x2+clxy*y2
x4~arx*x3+clxy*y3
x5~arx*x4+clxy*y4
x6~arx*x5+clxy*y5

y2~ary*y1+clyx*x1
y3~ary*y2+clyx*x2
y4~ary*y3+clyx*x3
y5~ary*y4+clyx*x4
y6~ary*y5+clyx*x5

#set comovements
vx1~~vy1          
vx2~~vy2          
vx3~~vy3          
vx4~~vy4          
vx5~~vy5          
vx6~~vy6

Sx~~Sy

#constraints
vx1~~0*Sx+0*Sy
vy1~~0*Sx+0*Sy
vx2~~0*Sx+0*Sy
vy2~~0*Sx+0*Sy
vx3~~0*Sx+0*Sy
vy3~~0*Sx+0*Sy
vx4~~0*Sx+0*Sy
vy4~~0*Sx+0*Sy
vx5~~0*Sx+0*Sy
vy5~~0*Sx+0*Sy
vx6~~0*Sx+0*Sy
vy6~~0*Sx+0*Sy

vx1~~0*vx2+0*vy2+0*vx3+0*vy3+0*vx4+0*vy4+0*vx5+0*vy5+0*vx6+0*vy6
vx2~~0*vx3+0*vy1+0*vy3+0*vx4+0*vy4+0*vx5+0*vy5+0*vx6+0*vy6
vx3~~0*vx4+0*vx5+0*vx6+0*vy1+0*vy2+0*vy4+0*vy5+0*vy6
vx4~~0*vx5+0*vx6+0*vy1+0*vy2+0*vy3+0*vy5+0*vy6
vx5~~0*vx6+0*vy1+0*vy2+0*vy3+0*vy4+0*vy6
vx6~~0*vy1+0*vy2+0*vy3+0*vy4+0*vy5
vy1~~0*vy2+0*vy3+0*vy4+0*vy5+0*vy6
vy2~~0*vy3+0*vy4+0*vy5+0*vy6
vy3~~0*vy4+0*vy5+0*vy6
vy4~~0*vy5+0*vy6
vy5~~0*vy6
'
GCLMALT_FreeLoad<-
  '
  #unit effects
Ix=~1*x6+1*x5+1*x4+1*x3+1*x2+1*x1
Iy=~1*y6+1*y5+1*y4+1*y3+1*y2+1*y1
Sx=~0*x1+1*x2+x3+x4+x5+x6
Sy=~0*y1+1*y2+y3+y4+y5+y6
Ix~1
Iy~1
Sx~1
Sy~1
x1~0*1
x2~0*1
x3~0*1
x4~0*1
x5~0*1
x6~0*1

y1~0*1
y2~0*1
y3~0*1
y4~0*1
y5~0*1
y6~0*1

#zero out measurement errors
x1~~0*x1
x2~~0*x2
x3~~0*x3
x4~~0*x4
x5~~0*x5
x6~~0*x6

y1~~0*y1
y2~~0*y2
y3~~0*y3
y4~~0*y4
y5~~0*y5
y6~~0*y6

#set impulses
vx1=~x1
vx2=~x2
vx3=~x3
vx4=~x4
vx5=~x5
vx6=~x6

vy1=~y1
vy2=~y2
vy3=~y3
vy4=~y4
vy5=~y5
vy6=~y6

#ARCL
x2~arx*x1+clxy*y1
x3~arx*x2+clxy*y2
x4~arx*x3+clxy*y3
x5~arx*x4+clxy*y4
x6~arx*x5+clxy*y5

y2~ary*y1+clyx*x1
y3~ary*y2+clyx*x2
y4~ary*y3+clyx*x3
y5~ary*y4+clyx*x4
y6~ary*y5+clyx*x5

#constraints
vx1~~0*Sx+0*Sy
vy1~~0*Sx+0*Sy
vx2~~0*Sx+0*Sy
vy2~~0*Sx+0*Sy
vx3~~0*Sx+0*Sy
vy3~~0*Sx+0*Sy
vx4~~0*Sx+0*Sy
vy4~~0*Sx+0*Sy
vx5~~0*Sx+0*Sy
vy5~~0*Sx+0*Sy
vx6~~0*Sx+0*Sy
vy6~~0*Sx+0*Sy

vx1~~0*vx2+0*vy2+0*vx3+0*vy3+0*vx4+0*vy4+0*vx5+0*vy5+0*vx6+0*vy6
vx2~~0*vx3+0*vy1+0*vy3+0*vx4+0*vy4+0*vx5+0*vy5+0*vx6+0*vy6
vx3~~0*vx4+0*vx5+0*vx6+0*vy1+0*vy2+0*vy4+0*vy5+0*vy6
vx4~~0*vx5+0*vx6+0*vy1+0*vy2+0*vy3+0*vy5+0*vy6
vx5~~0*vx6+0*vy1+0*vy2+0*vy3+0*vy4+0*vy6
vx6~~0*vy1+0*vy2+0*vy3+0*vy4+0*vy5
vy1~~0*vy2+0*vy3+0*vy4+0*vy5+0*vy6
vy2~~0*vy3+0*vy4+0*vy5+0*vy6
vy3~~0*vy4+0*vy5+0*vy6
vy4~~0*vy5+0*vy6
vy5~~0*vy6

#set comovements
vx1~~vy1          
vx2~~vy2          
vx3~~vy3          
vx4~~vy4          
vx5~~vy5          
vx6~~vy6

Ix~~Iy
Ix~~Sx
Ix~~Sy
Iy~~Sx
Iy~~Sy
Sx~~Sy
'
