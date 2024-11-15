set.seed(315340)
w<-rnorm(10,0,1)
Random_Walk<-cumsum(w)
Random_Walk_with_Drift<-cumsum(w+.6)
Linear_Trend<-.6*(1:10)
plot.ts(Random_Walk,ylim=c(-1,9),main="Stochastic & Deterministic Trends (t=10)",ylab="Scores"); 
lines(Random_Walk_with_Drift,lty="longdash");
lines(Linear_Trend,lty="dotted");
legend(1,9,legend=c("RW","RW+Drift","Linear Trend"),
	lty=c(1,5,3), box.lty=0,cex=.6)
