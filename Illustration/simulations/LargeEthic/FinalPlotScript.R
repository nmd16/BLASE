# Explore creating a graph with CIs. 

boxplot(cbind(GTHolder,GHolder,GMEHolder),at=c(1,4,7,10,13,16,19,2,5,8,11,14,17,20,3,6,9,12,15,18,21),col=rep(c("red","white","blue"),each=7))
points( rep(truthCoefs[-1],each=3), pch = 19, col = "green", lwd = 2)
abline( v = c(3.5,6.5,9.5,12.5,15.5,18.5,21.5))
legend( "topright", c("PB","GM","BLASE"),fill = c("white","red","blue"), cex = .5)

GTMean = apply(GTHolder,2,mean)
GMean = apply(GHolder,2,mean)
GMEMean = apply(GMEHolder,2,mean)
plot(c(GTMean,GMean,GMEMean),pch=16)
points(rep(truthCoefs[-1],3),pch=21,col="green")
abline(v = c(7.5,14.5))

# Plot with CIs. 
plotCI(1:7,GTMean,ui = GTMean + 1.96*apply(GTHolder,2,sd),li = GTMean- 1.96*apply(GTHolder,2,sd))
# Now we need to figure out how to do it all at once 

ExtMean = apply(ExtHolder,2,median)

ui1 =  GTMean + 1.96*apply(GTHolder,2,sd)
ui2 =  GMean + 1.96*apply(GHolder,2,sd)
ui3 =  ExtMean + 1.96*apply(ExtHolder,2,sd)
li1 =  GTMean - 1.96*apply(GTHolder,2,sd)
li2 =  GMean - 1.96*apply(GHolder,2,sd)
li3 =  ExtMean - 1.96*apply(ExtHolder,2,sd)
MeanHolder = c(GTMean,GMean,ExtMean)
names(MeanHolder) = rep(c("PB","GM","BL"),7)
MeanHolder = as.matrix(MeanHolder)
plotCI(c(2,6,10,14,18,22,26,3,7,11,15,19,23,27,4,8,12,16,20,24,28),MeanHolder,ui = c(ui1,ui2,ui3),li = c(li1,li2,li3),xlim = c(1,30),xlab ="",ylab = "",ylim = c(-2.2,2.2),axes=FALSE)
points(x = c(1,5,9,13,17,21,25), y = truthCoefs[-1],pch = 8)
axis(side=3,tck=0,at=c(0:30),labels=rep(c(""),31))
axis(side = 3, at=c(2,6,10,14,18,22,26), xlab = "Parameter",labels=c("Math","SexM","EthB","EthH","Math:Male","Male:EthB","Male:EthH"))
axis(side = 2, at=seq(-2.5,2.5,0.5))
axis(side = 1, at=c(0,2,6,10,14,18,22,26,3,7,11,15,19,23,27,4,8,12,16,20,24,28),labels=c("",names(MeanHolder)))


PlotHolder = cbind(GTHolder,GHolder,ExtHolder)
colnames(PlotHolder) = rep(c("PB","GM","BL"),each=7)
boxplot(PlotHolder,at=c(2,6,10,14,18,22,26,3,7,11,15,19,23,27,4,8,12,16,20,24,28), xlim = c(1,30))#,col=rep(c("red","white","blue"),each=7))
points( x= c(1,5,9,13,17,21,25), y = rep(truthCoefs[-1]), pch = 8, lwd = 2)
axis(side = 3, at=c(2,6,10,14,18,22,26), xlab = "Parameter",labels=c("Math","SexM","EthB","EthH","Math:Male","Male:EthB","Male:EthH"))
#abline( v = c(3.5,6.5,9.5,12.5,15.5,18.5,21.5))
#legend( "topright", c("PB","GM","BLASE"),fill = c("white","red","blue"), cex = .5)

