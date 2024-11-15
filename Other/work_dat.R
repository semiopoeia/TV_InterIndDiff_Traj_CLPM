View(MCMC)
#write.table(MCMC,"C:/Users/scottpw/Desktop/RR_Simul&Files/MCMC_314.dat",
 #           sep="\t",row.names = FALSE)

#library(haven)
#write_sav(MCMC,"C:/Users/scottpw/Desktop/RR_Simul&Files/MCMC.sav")

#setting to get Bayes factor
prep4BF<-
  MCMC %>%
  filter(FitMod!="alt_noslpvar")%>%
  select(aBIC_Mean,GenMod,FitMod,LCov,LVSx,LVSy,LCL,LARx,LARy)%>%
  arrange(LCov,LVSx,LVSy,LCL,LARx,LARy)%>%
  mutate(Gen2Fit=paste0(GenMod,FitMod))
set4BF<-split.data.frame(prep4BF,prep4BF$Gen2Fit)

BF<-function(j,k){
  BICdiff<-(set4BF[[j]]$aBIC_Mean-set4BF[[k]]$aBIC_Mean)
  BFdf<-round(exp(BICdiff),4)
  BFdf<-cbind(set4BF[[k]],BFdf)
  BFdf$Favors[BFdf$BFdf<1/150]<-"Very Strongly j"
  BFdf$Favors[BFdf$BFdf>=1/150&BFdf$BFdf<(1/20)]<-"Strongly j"
  BFdf$Favors[BFdf$BFdf>=(1/20)&BFdf$BFdf<(1/3)]<-"Moderately j"
  BFdf$Favors[BFdf$BFdf>=(1/3)&BFdf$BFdf<1]<-"Weakly j"
  BFdf$Favors[BFdf$BFdf>=1&BFdf$BFdf<3]<-"Weakly k"
  BFdf$Favors[BFdf$BFdf>=3&BFdf$BFdf<20]<-"Moderately k"
  BFdf$Favors[BFdf$BFdf>=20&BFdf$BFdf<150]<-"Strongly k"
  BFdf$Favors[BFdf$BFdf>=150]<-"Very Strongly k"
  k<-print(k)
  j<-print(j)
  BF<-BFdf$BFdf
  Favors<-BFdf$Favors
  LCov<-BFdf$LCov
  LVSx<-BFdf$LVSx
  LVSy<-BFdf$LVSy
  LCL<-BFdf$LCL
  LARx<-BFdf$LARx
  LARy<-BFdf$LARy
return(tibble(k,j,BICdiff,BF,Favors,LCov,LVSx,LVSy,LCL,LARx,LARy))
}

chk<-BF("LGCM_SRri_clpm","LGCM_SRgclm_meanstaionary")
table(chk$Favors)
View(chk)

riclpm<-165385.025
lgcmsr<-165038.609
alt<-165083.672
gclmMS<-165422.877
gclm<-165082.313

exp(alt-gclm)
