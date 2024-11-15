#10A
set.seed(315340)
wA10<-rnorm(10,0,1)
Random_Walk10A<-cumsum(wA10)
Random_Walk_with_Drift10A<-cumsum(wA10+.6)
Linear_Trend10A<-.6*(1:10)


#100A
set.seed(315340)
wA100<-rnorm(100,0,1)
Random_Walk100A<-cumsum(wA100)
Random_Walk_with_Drift100A<-cumsum(wA100+.6)
Linear_Trend100A<-.6*(1:100)

#10B
set.seed(3140)
wB10<-rnorm(10,0,1)
Random_Walk10B<-cumsum(wB10)
Random_Walk_with_Drift10B<-cumsum(wB10+.6)
Linear_Trend10B<-.6*(1:10)


#100B
set.seed(3140)
wB100<-rnorm(100,0,1)
Random_Walk100B<-cumsum(wB100)
Random_Walk_with_Drift100B<-cumsum(wB100+.6)
Linear_Trend100B<-.6*(1:100)

par(mfrow=c(2,2))
#10A
plot.ts(Random_Walk10A,ylim=c(-1,6),main="Stochastic & Deterministic Trends (t=10,i=A)",ylab="Scores"); 
lines(Random_Walk_with_Drift10A,lty="longdash");
lines(Linear_Trend10A,lty="dotted");
legend(1,6,legend=c("RW","RW+Drift","Linear Trend"),
	lty=c(1,5,3), box.lty=0,cex=.8)

#100A
plot.ts(Random_Walk100A,ylim=c(-1,65),main="Stochastic & Deterministic Trends (t=100,i=A)",ylab="Scores"); 
lines(Random_Walk_with_Drift100A,lty="longdash");
lines(Linear_Trend100A,lty="dotted");
legend(1,60,legend=c("RW","RW+Drift","Linear Trend"),
	lty=c(1,5,3), box.lty=0,cex=.8)

#10B
plot.ts(Random_Walk10B,ylim=c(-6,10),main="Stochastic & Deterministic Trends (t=10,i=B)",ylab="Scores"); 
lines(Random_Walk_with_Drift10B,lty="longdash");
lines(Linear_Trend10B,lty="dotted");
legend(1,10,legend=c("RW","RW+Drift","Linear Trend"),
	lty=c(1,5,3), box.lty=0,cex=.8)

#100B
plot.ts(Random_Walk100B,ylim=c(-25,65),main="Stochastic & Deterministic Trends (t=100,i=B)",ylab="Scores"); 
lines(Random_Walk_with_Drift100B,lty="longdash");
lines(Linear_Trend100B,lty="dotted");
legend(0,65,legend=c("RW","RW+Drift","Linear Trend"),
	lty=c(1,5,3), box.lty=0,cex=.8)
