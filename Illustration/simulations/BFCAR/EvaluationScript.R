## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##    Posterior Summaries    ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ## 

C2011 = File2011T[,1] - mean(File2011T[,1])
truthModel = lm(C2011 ~ Cmathscal + ethnic+sex, data = File2011T)
truthCoefs = truthModel$coefficients

SeedModel = lm( C2011[type1seeds] ~ Cmathscal[type1seeds] + ethnic+sex, data = File2011T[type1seeds,])
SeedCoefs = SeedModel$coefficients

load("N:/Simulations/NewDataRun/PreBLASE.RData")

load("N:/Simulations/BFCAR/StartSpace.RData")

##%%%%%%%%%%%%%%%%%%%%
## Load the BLASE Run

FullDataRunThetaOut <- read.csv("N:/Simulations/BFCAR/MainErrorOrigEthnicRunThetaOut.txt", quote="\"",header=F)

p= 6

Ext <- FullDataRunThetaOut 

Ext <- as.matrix( Ext)
Ext <- strsplit( Ext, " ")
ExtNew2 <- unlist(Ext)
ExtNew2 <- matrix( ExtNew2, ncol = p, byrow = T )

ExtHold <- ExtNew2

ExtNew2 <-as.numeric(ExtHold)
ExtNew2 <- matrix( ExtNew2, ncol = p, nrow = nrow(ExtHold), byrow = F )

rm(FullDataRunThetaOut, ExtHold,Ext)
##%%%%%%%%%%%%%%%%%%%%
## Load the GM Run

GMThetaOut <- read.csv("N:/Simulations/BFCAR/GMainErrorOrigEthnicRunThetaOut.txt", quote="\"",header=F)

p= 6

Ext <- GMThetaOut 

Ext <- as.matrix( Ext)
Ext <- strsplit( Ext, " ")
GM <- unlist(Ext)
GM <- matrix( GM, ncol = p, byrow = T )

ExtHold <- GM

GM <-as.numeric(ExtHold)
GM<- matrix( GM, ncol = p, nrow = nrow(ExtHold), byrow = F )

rm(GMThetaOut, ExtHold,Ext)

##%%%%%%%%%%%%%%%%%%%%
## Load the PB Run

PBThetaOut <- read.csv("N:/Simulations/BFNAR/GGTestDelibRunThetaOut.txt", quote="\"",header=F)

p= 6

Ext <- PBThetaOut 

Ext <- as.matrix( Ext)
Ext <- strsplit( Ext, " ")
PB <- unlist(Ext)
PB <- matrix( PB, ncol = p, byrow = T )

ExtHold <- PB

PB <-as.numeric(ExtHold)
PB<- matrix( PB, ncol = p, nrow = nrow(ExtHold), byrow = F )
PB[,1]=PB[,1]-mean(File1[,1])
rm(PBThetaOut, ExtHold,Ext)

## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##   CI Plots with PB        ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

ExtNew2[,1] = ExtNew2[,1] - 360.75
GM[,1] = GM[,1] - 360.75

burnin = 200
meanvec = c( apply(PB[-c(1:burnin),1:5],2,mean),apply(GM[-c(1:burnin),1:5],2,mean), apply(ExtNew2[-c(1:burnin),1:5],2,mean))
meanvec = meanvec[c(1,6,11,2,7,12,3,8,13,4,9,14,5,10,15)]
li1     = c( apply(PB[-c(1:burnin),1:5],2,function(x) quantile(x,0.025)), apply(GM[-c(1:burnin),1:5],2,function(x) quantile(x,0.025)), apply(ExtNew2[-c(1:burnin),1:5],2,function(x) quantile(x,0.025)) ) 
li1     = li1[c(1,6,11,2,7,12,3,8,13,4,9,14,5,10,15)]
ui1     = c( apply(PB[-c(1:burnin),1:5],2,function(x) quantile(x,0.975)),apply(GM[-c(1:burnin),1:5],2,function(x) quantile(x,0.975)), apply(ExtNew2[-c(1:burnin),1:5],2,function(x) quantile(x,0.975)) ) 
ui1     = ui1[c(1,6,11,2,7,12,3,8,13,4,9,14,5,10,15)]

pdf("N:/Simulations/BFCAR/Figures/CIPB_NoLines.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( PB[c(100:341),1],GM[c(100:341),1],ExtNew2[-c(1:50),1] , PB[c(100:341),2], GM[c(100:341),2],ExtNew2[-c(1:50),2],PB[c(100:341),3], GM[c(100:341),3],ExtNew2[-c(1:50),3],PB[c(100:341),4],GM[c(100:341),4],ExtNew2[-c(1:50),4],PB[c(100:341),5],GM[c(100:341),5],ExtNew2[-c(1:50),5],border="white",xlim=c(0,15),names=rep(c("PB","GM","BL"),5))
plotCI(1:15,meanvec,ui=ui1,li=li1,xlab ="",ylab = "Point Estimate", add=T)
#points(x = c(3,7,11,15,19), y = truthCoefs,pch = 15)
axis(side = 3, at=c(2,5,8,11,14), xlab = "Parameter",labels=c("Int","Math","EthB","EthH","Male"),outer=F,tick = c(0,4,7,10,13))
segments(.5,truthCoefs[1],3.5,truthCoefs[1])
segments(4,truthCoefs[2],6.25,truthCoefs[2])
segments(6.5,truthCoefs[3],9.25,truthCoefs[3])
segments(9.5,truthCoefs[4],12.25,truthCoefs[4])
segments(12.5,truthCoefs[5],15,truthCoefs[5])
dev.off()

pdf("N:/Simulations/BFCAR/Figures/CIPB_Intercept.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( PB[c(100:341),1],GM[c(100:341),1],ExtNew2[-c(1:50),1],border="white",xlim=c(0,4),names=rep(c("PB","GM","BL"),1),ylab="Intercept")
plotCI(1:3,meanvec[1:3],ui=ui1[1:3],li=li1[1:3],xlab ="",ylab = "Point Estimate", add=T)
#points(x = c(3,7,11,15,19), y = truthCoefs,pch = 15)
segments(1,truthCoefs[1],3,truthCoefs[1])
dev.off()

pdf("N:/Simulations/BFCAR/Figures/CIPB_Math.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( PB[c(100:341),2],GM[c(100:341),2],ExtNew2[-c(1:50),2],border="white",xlim=c(0,4),names=rep(c("PB","GM","BL"),1),ylab="Math")
plotCI(1:3,meanvec[4:6],ui=ui1[4:6],li=li1[4:6],xlab ="",ylab = "Point Estimate", add=T)
#points(x = c(3,7,11,15,19), y = truthCoefs,pch = 15)
segments(1,truthCoefs[2],3,truthCoefs[2])
dev.off()

pdf("N:/Simulations/BFCAR/Figures/CIPB_EthB.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( PB[c(100:341),3],GM[c(100:341),3],ExtNew2[-c(1:50),3],border="white",names=rep(c("PB","GM","BL"),1),ylab="Ethnicity B")
plotCI(1:3,meanvec[7:9],ui=ui1[7:9],li=li1[7:9],xlab ="",ylab = "Point Estimate", add=T)
#points(x = c(3,7,11,15,19), y = truthCoefs,pch = 15)
segments(1,truthCoefs[3],3,truthCoefs[3])
dev.off()

pdf("N:/Simulations/BFCAR/Figures/CIPB_EthH.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( PB[c(100:341),4],GM[c(100:341),4],ExtNew2[-c(1:50),4],border="white",names=rep(c("PB","GM","BL"),1),ylab="Ethnicity H")
plotCI(1:3,meanvec[10:12],ui=ui1[10:12],li=li1[10:12],xlab ="",ylab = "Point Estimate", add=T)
#points(x = c(3,7,11,15,19), y = truthCoefs,pch = 15)
segments(1,truthCoefs[4],3,truthCoefs[4])
dev.off()

pdf("N:/Simulations/BFCAR/Figures/CIPB_Male.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( PB[c(100:341),5],GM[c(100:341),5],ExtNew2[-c(1:50),5],border="white",names=rep(c("PB","GM","BL"),1),ylab="Male")
plotCI(1:3,meanvec[13:15],ui=ui1[13:15],li=li1[13:15],xlab ="",ylab = "Point Estimate", add=T)
#points(x = c(3,7,11,15,19), y = truthCoefs,pch = 15)
segments(1,truthCoefs[5],3,truthCoefs[5])
dev.off()


## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##         CI Plots          ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

GM[,1] = GM[,1]-mean(File2011T[,1])
ExtNew2[,1] = ExtNew2[,1]-mean(File2011T[,1])

meanvec = c( apply(GM[-c(1:burnin),1:5],2,mean), apply(ExtNew2[-c(1:burnin),1:5],2,mean))
meanvec = meanvec[c(1,6,2,7,3,8,4,9,5,10)]
li1     = c( apply(GM[-c(1:burnin),1:5],2,function(x) quantile(x,0.025)), apply(ExtNew2[-c(1:burnin),1:5],2,function(x) quantile(x,0.025)) ) 
li1     = li1[c(1,6,2,7,3,8,4,9,5,10)]
ui1     = c( apply(GM[-c(1:burnin),1:5],2,function(x) quantile(x,0.975)), apply(ExtNew2[-c(1:burnin),1:5],2,function(x) quantile(x,0.975)) ) 
ui1     = ui1[c(1,6,2,7,3,8,4,9,5,10)]
  
pdf("N:/Simulations/BFCAR/Figures/CI_All.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( GM[c(100:341),1],ExtNew2[-c(1:50),1] ,GM[c(100:341),2],ExtNew2[-c(1:50),2],GM[c(100:341),3],ExtNew2[-c(1:50),3],GM[c(100:341),4],ExtNew2[-c(1:50),4],GM[c(100:341),5],ExtNew2[-c(1:50),5],border="white",at=c(1,3,4,6,7,9,10,12,13,15),xlim=c(0,15),names=rep(c("GM","B"),5),ylim=c(-2,1.1),xlab = "BFCAR") 
plotCI(c(1,3,4,6,7,9,10,12,13,15),meanvec,ui=ui1,li=li1,xlab ="",ylab = "Point Estimate", add=T)
points(x = c(2,5,8,11,14), y = truthCoefs,pch = 15)
axis(side = 3, at=c(2,5,8,11,14), xlab = "Parameter",labels=c("Int","Math","EthB","EthH","Male"),outer=F,tick = c(0,4,7,10,13))
#abline( v= c(3.5,6.5,9.5,12.5))
dev.off()

# Int ONLY 
pdf("N:/Simulations/BFCAR/Figures/CI_Intercept.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( GM[c(200:441,1)], ExtNew2[-c(1:50),1], at=c(2,3),border="white",names=c("GM","BLASE"),ylim=c(0,.4),ylab = "Intercept")
points(x =2.5, y = truthCoefs[1], pch = 15)
plotCI(c(2,3),meanvec[1:2],ui=ui1[1:2],li=li1[1:2],xlab ="",ylab = "Point Estimate", add=T)
axis(side = 2, at=c(0.1), ylab = "Intercept",outer=F)
dev.off()

# MATH ONLY 
pdf("N:/Simulations/BFCAR/Figures/CI_Math.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( GM[c(200:441,2)], ExtNew2[-c(1:50),2], at=c(2,3),border="white",ylim=c(.86,.9),names=c("GM","BLASE"),ylab="Math")
points(x =2.5, y = truthCoefs[2], pch = 15)
plotCI(c(2,3),meanvec[3:4],ui=ui1[3:4],li=li1[3:4],xlab ="",ylab = "Point Estimate", add=T)
dev.off()

# EthB ONLY 
pdf("N:/Simulations/BFCAR/Figures/CI_EthB.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( GM[c(200:441),3], ExtNew2[-c(1:50),3], at=c(2,3),border = "white", names=c("GM","BLASE"),ylab = "Ethnic B",ylim =c(-1.6,-1))
points(x =2.5, y = truthCoefs[3], pch = 15)
axes(2,"EthnicB")
plotCI(c(2,3),meanvec[5:6],ui=ui1[5:6],li=li1[5:6],xlab ="",ylab = "Point Estimate", add=T)
dev.off()

# EthH ONLY 
pdf("N:/Simulations/BFCAR/FiguresCI_EthH.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( GM[c(200:441),4], ExtNew2[-c(1:50),4], at=c(2,3),border = "white", names=c("GM","BLASE"),ylab = "Ethnic H")
points(x =2.5, y = truthCoefs[4], pch = 15)
plotCI(c(2,3),meanvec[7:8],ui=ui1[7:8],li=li1[7:8],xlab ="",ylab = "Point Estimate", add=T)
dev.off()

# SexM ONLY 
pdf("N:/Simulations/BFCAR/Figures/CI_SexM.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( GM[c(200:441),5], ExtNew2[-c(1:50),5], at=c(2,3),border = "white", names=c("GM","BLASE"),ylab = "Sex M")
points(x =2.5, y = truthCoefs[5], pch = 15)
plotCI(c(2,3),meanvec[9:10],ui=ui1[9:10],li=li1[9:10],xlab ="",ylab = "Point Estimate", add=T)
dev.off()

## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##       Match Rates         ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

# Now we want to make the same kind of plots, but for the match rate.

BMatch <- read.csv("N:/Simulations/MainOnly/MainOrigEthnicRun/MainErrorOrigEthnicRunLambdaOut.txt", quote="\"",header=F)
BMatch = as.numeric(unlist(BMatch))
GMatch <- read.csv("N:/Simulations/MainOnly/MainOrigEthnicRun/GMainErrorOrigEthnicRunLambdaOut.txt", quote="\"",header=F)
GMatch = as.numeric(unlist(GMatch))

meanvec = c( mean(GMatch[-c(1:100)]), mean(BMatch[-c(1:100)]))/(N/2)
li1     = c(quantile(GMatch[-c(1:100)],0.025),quantile(BMatch[-c(1:100)],0.025) )/(N/2) 
ui1     = c(quantile(GMatch[-c(1:100)],0.975),quantile(BMatch[-c(1:100)],0.975) )/(N/2)

pdf("N:/Simulations/BFCAR/Figures/MatchPercent.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( GMatch[c(200:441)]/(N/2), BMatch[-c(1:50)]/(N/2), at=c(2,3),border = "white", names=c("GM","BLASE"),ylab = "Matches")
#points(x =2.5, y = truthCoefs[5], pch = 15)
plotCI(c(2,3),meanvec,ui=ui1,li=li1,xlab ="",ylab = "Point Estimate", add=T)
dev.off()

meanvec = c( mean(GMatch[-c(1:100)]), mean(BMatch[-c(1:100)]))
li1     = c(quantile(GMatch[-c(1:100)],0.025),quantile(BMatch[-c(1:100)],0.025) )
ui1     = c(quantile(GMatch[-c(1:100)],0.975),quantile(BMatch[-c(1:100)],0.975) )

pdf("N:/Simulations/BFCAR/Figures/Match.pdf",width = 11, height = 7 )
#par(mar=c(5,3,2,2)+0.1)
boxplot( GMatch[c(200:441)], BMatch[-c(1:50)], at=c(2,3),border = "white", names=c("GM","BLASE"),ylab = "Matches")
#points(x =2.5, y = truthCoefs[5], pch = 15)
plotCI(c(2,3),meanvec,ui=ui1,li=li1,xlab ="",ylab = "Point Estimate", add=T)
dev.off()

plot(BMatch[-c(1:200)])

