## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##  Load the starting space  ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ## 

setwd("N:/Simulations/LargeEthnic")

load("N:/Simulations/LargeEthnic/TruthData.RData")

truthModel = lm( File2011[,1] ~ Cmathscal + ethnic+sex, data = File2011)
truthCoefs = truthModel$coefficients

SeedsModel = lm( File2011[type1seeds,1] ~ Cmathscal[type1seeds] + ethnic+sex, data = File2011[type1seeds,])
SeedsCoefs = SeedsModel$coefficients
AllModel = lm( File1_data[,1] ~ File2_data[order(order(File2_data[,"lambda"])),1] + File1_data[,2] + File1_data[,3] + File1_data[,4])

##%%%%%%%%%%%%%%%%%%%%
## Load the BLASE Run

FullDataRunThetaOut <- read.csv("N:/Simulations/LargeEthnic/MainShortRunThetaOut.txt", quote="\"",header=F)

FullDataRunThetaOut <- read.csv("N:/Simulations/LargeEthnic/MainOnlyRun2ThetaOut.txt", quote="\"",header=F)

p= 6

Ext <- FullDataRunThetaOut 

Ext <- as.matrix( Ext)
Ext <- strsplit( Ext, " ")
ExtNew2 <- unlist(Ext)
ExtNew2 <- matrix( ExtNew2, ncol = p, byrow = T )

ExtHold <- ExtNew2

ExtNew2 <-as.numeric(ExtHold)
ExtNew2 <- matrix( ExtNew2, ncol = p, nrow = nrow(ExtHold), byrow = F )

##%%%%%%%%%%%%%%%%%%%%
## Load the PB Run

GTest <- read.csv("N:/Simulations/LargeEthnic/LargeEthnicRunTestThetaOut.txt", quote="\"",header=F)

p= 9

Ext <- GTest

Ext <- as.matrix( Ext)
Ext <- strsplit( Ext, " ")
GTest <- unlist(Ext)
GTest <- matrix( GTest, ncol = p, byrow = T )

ExtHold <- GTest

GTest <-as.numeric(ExtHold)
GTest <- matrix( GTest, ncol = p, nrow = nrow(ExtHold), byrow = F )

##%%%%%%%%%%%%%%%%%%%%
## Load the GM Run
##%%%%%%%%%%%%%%%%%%%%
GM <- read.csv("N:/Simulations/LargeEthnic/LargeEthnicRunGMThetaOut.txt", quote="\"",header=F)

GME <- read.csv("N:/Simulations/LargeEthnic/ErrorRun/MoreError/MoreErrorGMThetaOut.txt", quote="\"",header=F)

load("~/GOutMainOnly.RData")
write(G.outlist$out.coefs,file = "GMainOnlyRun2ThetaOut.txt")

p= 9

Ext <- GME

Ext <- as.matrix( Ext)
Ext <- strsplit( Ext, " ")
GME <- unlist(Ext)
GME <- matrix( GME, ncol = p, byrow = T )

ExtHold <- GME

GME <-as.numeric(ExtHold)
GME <- matrix( GME, ncol = p, nrow = nrow(ExtHold), byrow = F )
GM = GME


N = 2*N

posterior_summaries_neat(G.outlist$out.coefs,c("Intercept", "CMath","Ethnic B", "Ethnic H", "SexM"),truthCoefs)

posterior_summaries_neat(GTest$out.coefs,c("Intercept", "Math","Year Post 1996", "Ethnic A", "Ethnic B", "Ethnic H", "Ethnic I", "Ethnic M","Sex Male"),truthCoefs)

posterior_summaries_neat(ExtNew2[-c(1:500),],c("Intercept", "CMath","Ethnic B", "Ethnic H", "SexM"),truthCoefs)

## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##         Boxplot           ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

boxplot(G.outlist$out.coefs[-c(1:1000),-c(1,p)])#,names=c("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"))
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

boxplot(ExtNew2[-c(1:100),-c(1,p)])
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

# Okay, so we can now see the ethnicity coefficients are messy. Good :) Let's hope that our method can do better! 
EHolderBox = (ExtNew2[-c(1:200),-c(1,9)] - truthCoefs[-1])/truthCoefs[-1]
boxplot(EHolderBox)
abline(h=0)
boxplot(ExtNew2[-c(1:200),-c(1,20)])#,names=c("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"))
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

range.E = nrow(ExtNew)- 50
range.E = length(G.outlist$out.coefs[,1])- 200


boxplot(cbind(G.outlist$out.coefs[-c(1:622),-c(1,p)],ExtNew2[-c(1:500),-c(1,p)]),at=c(1,3,5,7,2,4,6,8) ,col=rep(c("red","blue"),each=4)) #,names=rep(c ("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"),2)
points( rep(truthCoefs[-1],each=2), pch = 19, col = "green", lwd = 2)
abline( v = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5,30.5,32.5,34.5))
# Why is everything so crazy with BLASE? 

GHolder = G.outlist$out.coefs[-c(1:622),]
BHolder = ExtNew2[-c(1:500),]
hist(BHolder[,1])
hist(BHolder[,2])
hist(BHolder[,3])
hist(BHolder[,4])

plot(BHolder[,1])
plot(BHolder[,2])
plot(BHolder[,3])
plot(BHolder[,4])

Mode1 = which(BHolder[,2]>=.8)
mean(BHolder[Mode1,1]) # . 85
Mode2 = which(BHolder[,2]<.8 & BHolder[,1] > 0.6)
Mode3 = which(BHolder[,2]<=0.6)
mean(BHolder[Mode2,2]) # . 68
plot(x = Mode1, y = BHolder[Mode1,2],ylim = c(.2,1),xlim = c(1,5000))
points(x=Mode2, y = BHolder[Mode2,2],col="blue")
points(x= Mode3,y = BHolder[Mode3,2],col="red")
acf(BHolder[,4])


plot(x = Mode1, y = BHolder[Mode1,3],ylim = c(-20,0),xlim =c(0,5000))
points(x= Mode2, y = BHolder[Mode2,3],col="blue")
points(x= Mode3, y = BHolder[Mode3,3],col="red")

summary(BHolder[Mode1,3])
summary(BHolder[Mode2,3])

plot(x= Mode1, y = BHolder[Mode1,4],ylim = c(-20,0),xlim =c(0,5000))
points(x= Mode2, y = BHolder[Mode2,4],col="blue")
points(x = Mode3, y = BHolder[Mode3,4],col="red")

summary(BHolder[Mode1,4])
summary(BHolder[Mode2,4])

plot(x = Mode2, y = BHolder[Mode2,5],ylim = c(-1,1),xlim =c(0,5000),col="blue")
points(x= Mode1, y = BHolder[Mode1,5],col="black")
points(x= Mode3, y = BHolder[Mode3,5],col="red")

# Let's predict for all three modes. 
ModelMatrix = model.matrix( lm( File2011T[,1] ~ File2012T$Cmathscal + ethnic+sex, data = File2011T))
predictMode1 = rnorm( n1, ModelMatrix%*%apply(BHolder[Mode1,c(1:(p-1))],2,mean), sqrt(BHolder[Mode1,p]))
summary(predictMode1)
summary(File2011T[,1])
sqrt( sum(predictMode1 -File2011T[,1])^2/n1 ) 
predictTrue = rnorm( n1, ModelMatrix%*%truthCoefs, 5.041)
sqrt( sum(predictTrue -File2011T[,1])^2/n1 )
summary(predictTrue)

predictMode2 = rnorm( n1, ModelMatrix%*%apply(BHolder[Mode2,c(1:(p-1))],2,mean), sqrt(BHolder[Mode2,p]))
summary(predictMode2)  

predictMode3 = rnorm( n1, ModelMatrix%*%apply(BHolder[Mode3,c(1:(p-1))],2,mean), sqrt(BHolder[Mode3,p]))
summary(predictMode3) 

predictG = rnorm( n1, ModelMatrix%*%apply(GHolder[,c(1:(p-1))],2,mean), sqrt(GHolder[,p]))
summary(predictG) 
#Let's check the root MSE
sqrt( sum(predictG -File2011T[,1])^2/n1 ) 

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExtHolderS <- (ExtNew2[-c(1:500),-c(1,p)] - truthCoefs[-1])/truthCoefs[-1]
GHolderS <- (GM[-c(1:500,3879:5000),-c(1,p)] - truthCoefs[-1])/truthCoefs[-1]
GTHolderS <- (GTest[-c(1:500,3879:5000),-c(1,p)] - truthCoefs[-1])/truthCoefs[-1]
GMEHolderS <- (GME[-c(1:500,3879:5000),-c(1,p)] - truthCoefs[-1])/truthCoefs[-1]

ExtHolder <- (ExtNew2[-c(1:500),-c(1,p)])
GHolder <- (GM[-c(1:500,3879:5000),-c(1,p)] )
GTHolder <- (GTest[-c(1:500,3879:5000),-c(1,p)])
GMEHolder <- (GME[-c(1:500,3879:5000),-c(1,p)])

boxplot(cbind(GTHolderS,GHolderS,ExtHolderS),at=c(1,4,7,10,13,16,19,2,5,8,11,14,17,20,3,6,9,12,15,18,21),ylim = c(-50,50))#,col=rep(c("red","white","blue"),each=7))
points( rep(truthCoefs[-1],each=3), pch = 19, col = "green", lwd = 2)
abline( v = c(3.5,6.5,9.5,12.5,15.5,18.5,21.5))
legend( "topright", c("PB","GM","BLASE"),fill = c("white","red","blue"), cex = .5)

# PB vs. BLASE
boxplot(cbind(GTHolder[,1:7],ExtHolder[,1:7]),at=c(1,3,5,7,9,11,13,2,4,6,8,10,12,14),ylab = "% Distance from Truth") 
abline( v = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5))
mtext(c("CMath","SexM","EthB","EthH","Math:Sex", "SexM:B", "SexM:H"),1,line=3,at = c(1.5,3.5,5.5,7.5,9.5,11.5,13.5))
abline(h=0,lty=2)

# GM vs. BLASE
boxplot(cbind(GHolder[,1:7],ExtHolder[,1:7]),at=c(1,3,5,7,9,11,13,2,4,6,8,10,12,14),ylab = "% Distance from Truth") 
abline( v = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5))
mtext(c("CMath","SexM","EthB","EthH","Math:Sex", "SexM:B", "SexM:H"),1,line=3,at = c(1.5,3.5,5.5,7.5,9.5,11.5,13.5))
abline(h=0,lty=2)


## %%%%%%%%%%%%%%%%%%%%%%%% ## 
##   Trace Plots : Check    ##
## %%%%%%%%%%%%%%%%%%%%%%%% ## 

burnoff = 500
plot( G.outlist$out.coefs[-c(1:burnoff),1],type="l")
abline(h = truthCoefs[1],col="green")
plot( G.outlist$out.coefs[-c(1:burnoff),2],type="l")
abline(h = truthCoefs[2],col="green")
plot( G.outlist$out.coefs[-c(1:burnoff),3],type="l")
abline(h = truthCoefs[3],col="green")
# 3 = YES
plot( G.outlist$out.coefs[-c(1:burnoff),4],type="l")
abline(h = truthCoefs[4],col="green")
plot( G.outlist$out.coefs[-c(1:burnoff),5],type="l")
abline(h = truthCoefs[5],col="green")
# 5 = YES
plot( G.outlist$out.coefs[-c(1:burnoff),6],type="l")
abline(h = truthCoefs[6],col="green")
plot( G.outlist$out.coefs[-c(1:burnoff),7],type="l")
abline(h = truthCoefs[7],col="green")
# 7 = Barely 
plot( G.outlist$out.coefs[-c(1:burnoff),8],type="l")
abline(h = truthCoefs[8],col="green")
# 8 = YES 
plot( G.outlist$out.coefs[-c(1:burnoff),9],type="l")
abline(h = truthCoefs[9],col="green")

ExtNew = ExtNew2
burnoff = 100
plot( ExtNew[-c(1:burnoff),1])
abline(h=truthCoefs[1])
plot( ExtNew[-c(1:burnoff),2], type = "l")
abline(h=truthCoefs[2])
plot( ExtNew[-c(1:burnoff),3], type = "l")
abline(h=truthCoefs[3],col="green")
# 3 = YES
plot( ExtNew[-c(1:burnoff),4], type = "l")
abline(h=truthCoefs[4],col="green")
plot( ExtNew[-c(1:burnoff),5], type = "l")
abline(h=truthCoefs[5],col="green")
plot( ExtNew[-c(1:burnoff),6], type = "l")
abline(h=truthCoefs[6],col="green")
plot( ExtNew[-c(1:burnoff),7], type = "l")
abline(h=truthCoefs[7],col="green")
# 7 = Barely 
plot( ExtNew[-c(1:burnoff),8], type = "l")
abline(h=truthCoefs[8],col="green")
# 8 = YES
plot( ExtNew[-c(1:burnoff),9], type = "l")
abline(h=truthCoefs[9],col="green")


plot( GTest$out.coefs[-c(1:100),1],type="l")
abline(h = truthCoefs[1],col="green")
plot( GTest$out.coefs[-c(1:100),2],type="l")
abline(h = truthCoefs[2],col="green")
plot( GTest$out.coefs[-c(1:100),3],type="l")
abline(h = truthCoefs[3],col="green")
plot( GTest$out.coefs[-c(1:100),4],type="l")
abline(h = truthCoefs[4],col="green")
plot( GTest$out.coefs[,5],type="l")
abline(h = truthCoefs[5],col="green")
plot( GTest$out.coefs[,6],type="l")
abline(h = truthCoefs[6],col="green")
plot( GTest$out.coefs[,7],type="l")
abline(h = truthCoefs[7],col="green")
plot( GTest$out.coefs[,8],type="l")
abline(h = truthCoefs[8],col="green")
plot( GTest$out.coefs[,9],type="l")
abline(h = truthCoefs[9],col="green")
plot( GTest$out.coefs[,10],type="l")
abline(h = truthCoefs[10],col="green")
# 10 = YES 
plot( GTest$out.coefs[,11],type="l")
abline(h = truthCoefs[11],col="green")
# 11 = YES
plot( GTest$out.coefs[,12],type="l")
abline(h = truthCoefs[12],col="green")
# 12 = YES 
plot( GTest$out.coefs[,13],type="l")
abline(h = truthCoefs[13],col="green")
# 13 = YES
plot( GTest$out.coefs[,14],type="l")
abline(h = truthCoefs[14],col="green")
# 14 = YES 
plot( GTest$out.coefs[,15],type="l")
abline(h = truthCoefs[15],col="green")
# 15 = YES 
plot( GTest$out.coefs[,16],type="l")
abline(h = truthCoefs[16],col="green")
plot( GTest$out.coefs[,17],type="l")
abline(h = truthCoefs[17],col="green")
plot( GTest$out.coefs[,18],type="l")
abline(h = truthCoefs[18],col="green")

## %%%%%%%%%%%%%%%%%%%%%%%% ## 
##       T-values           ##
## %%%%%%%%%%%%%%%%%%%%%%%% ##

names(truthModel)

summary(truthModel)

tvalues = c(7.487,263.988,6.234,-13.392,-10.361,-4.534,-5.437,-2.221)
tranks = order(abs(tvalues[-1]),decreasing=T)

#We want positive values to mean that BLASe does better 
GMvsBLASE = abs(GHolder)-abs(ExtHolder)
dim(GMvsBLASE)
pdf("N:/Simulations/SimulationOutput/Figures/Boxplot_ByRank.pdf")
boxplot(GMvsBLASE, at = order(tranks),names=order(tranks),xlab = "Ranks:T-values",ylab = "Abs Diff", main =" Gain in Absolute Percent Difference: BLASE vs.GM")
abline(h = 0 ,lty = 2)
dev.off()

G.mean = apply(GHolder,2,mean)
Ext.mean = apply(ExtHolder,2,mean)

pdf("Gains_BLASEvsGM_byRank.pdf")
plot(y = (abs(G.mean) - abs(Ext.mean))[order(tranks)], x = 1:(p-2),xlab = "Ranks:T-values",ylab = "Abs Diff", main =" Gain in Absolute Percent Difference: BLASE vs.GM")
abline(h = 0 ,lty = 2)
dev.off()

plot(y = ((abs(G.mean) - abs(Ext.mean))/(abs(G.mean)))[order(tranks)], x = 1:(p-2),xlab = "Ranks:T-values",ylab = "Abs Diff", main =" Gain in Absolute Percent Difference: BLASE vs.GM")
abline(h = 0 ,lty = 2)

## %%%%%%%%%%%%%%%%%%%%%%%% ## 
##    Error Distribution    ##
## %%%%%%%%%%%%%%%%%%%%%%%% ##

setcheck = setdiff( type1seeds, which(File2012[,"ethnic"]!=File2011[,"ethnic"]) %in% in.error)

table(File2012[in.error, "ethnic"],File2011[in.error,"ethnic"])

table(File2012[type1seeds, "ethnic"],File2011[type1seeds,"ethnic"])/n.type1


## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##       Match Rate          ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

LAM <- G.outlist$out.lambda
save(LAM, file = "MainRun2GoutLambda.txt")
LAM <- G.outlist$out.coefs
save(LAM, file = "GoutCoefs.txt")
LAM <- G.outlist$out.coefs2
save(LAM, file = "GoutCoefs2.txt")

# BLASE - Seeds and Non Seeds
FullDataRunAllLambdaPerc <- read.table("N:/Simulations/LargeEthnic/MainOnlyRun2LambdaOut.txt", quote="\"")
(summary(FullDataRunAllLambdaPerc[-c(1:100),])-n.type1)/(N/2)*100

plot(unlist(FullDataRunAllLambdaPerc),type="l")
sum(File2_data[,"Seed"]==1)
intersect(type1seeds,MoveCand_F2)

# BLASE - original records only 
FullLamOrig <- read.table("N:/Simulations/NewDataRun/NewDataRunLambdaOut.txt", quote="\"")
summary(FullLamOrig[-c(1:20),])/N

# GM - Original Only (1666 posterior summaries)
holder = apply( G.outlist$out.lambda,1, function(x) sum(x[1:N]==1:N))
summary(holder[-c(1:100)]/N)

# PB 
holderT = apply( GTest$out.lambda,1, function(x) sum(x[1:N]==1:N))
summary(holderT[-c(1:100)]/N)

# Compare the match rates considering ALL original records
N = 9999

# Gutman 
summary(holder[-c(1:400)]/dim(G.outlist$out.lambda)[2])

# Extension 
Run2AllLambdaPerc <- read.table("C:/Users/nmd16/Desktop/Final Runs/Run2/Run2AllLambdaPerc.txt", quote="\"")
summary(Run2AllLambdaPerc[-c(1:400),])

plot(as.matrix(Run2AllLambdaPerc), type="l")
plot(as.matrix(Run2AllLambdaPerc[-c(1:400),]), type="l")

## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##   Posterior for Gamma     ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

GammaOut <- read.table("C:/Users/nmd16/Desktop/Final Runs/Run1/Run1GammaOut1.txt", quote="\"")

summary(GammaOut)

AcceptOut <- read.table("N:/Simulations/LargeEthnic/MainShortRunAcceptOut.txt", quote="\"")

plot( as.matrix(AcceptOut[-c(1:10),]), type= "l")

SizeOut <- read.table("N:/Simulations/LargeEthnic/MainShortRunSizeOut.txt", quote="\"")

plot( as.matrix(SizeOut), type= "l")

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
##   Evaluating the Block Sizes   ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

# On the error free data 

Bhat1          = File2012[,2:7]
d = sapply(Bhat1,function(x) length(unique(x)) )

blocks <-NULL 

for( i in 1:(N)){
  blocks[i] = category.block(Bhat1[i,],d,6)
}

# On the error prone data 

load("C:/Users/nmd16/Desktop/Final Runs/Run2/File1Formatted.RData")

table(table(File1_data[,'Block']))


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
##   Evaluating the Block Moves   ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

Ext <- as.matrix( ExtNew, ncol = 10 )
ExtNew <- matrix( Ext, ncol = 10, nrow = length(ExtNew)/10, byrow = T )

LOut <- read.csv("C:/Users/nmd16/Desktop/Final Runs/Run2/Run2LambdaAll.csv", header=F, sep = "")

LM <-LOut 
LM = LM[!is.na(LOut)]

LM <- c(LOut[-nrow(LOut),ncol(LOut)])
LM <- matrix( LM, ncol = 9999, nrow = 100, byrow = T)

HolderL = LM[1,]
for( i in 2:2000){
  HolderL <-cbind(HolderL,LM[i,])
}
HolderL <-t(HolderL)
Holder <- HolderL[-nrow(HolderL),]

# Explore creating a graph with CIs. 

boxplot(cbind(GTHolder,GHolder,GMEHolder),at=c(1,4,7,10,13,16,19,2,5,8,11,14,17,20,3,6,9,12,15,18,21),col=rep(c("red","white","blue"),each=7))
points( rep(truthCoefs[-1],each=3), pch = 19, col = "green", lwd = 2)
abline( v = c(3.5,6.5,9.5,12.5,15.5,18.5,21.5))
legend( "topright", c("PB","GM","BLASE"),fill = c("white","red","blue"), cex = .5)

GTMean = apply(GTHolder,2,mean)
GMean = apply(GHolder[,-c(1,6)],2,mean)
GMEMean = apply(GMEHolder,2,mean)
plot(c(GTMean,GMean,GMEMean),pch=16)
points(rep(truthCoefs[-1],3),pch=21,col="green")
abline(v = c(7.5,14.5))

# Plot with CIs. 
plotCI(1:4,GMean,ui = GMean + 1.96*apply(GHolder[,-c(1,6)],2,sd),li = GMean- 1.96*apply(GHolder[,-c(1,6)],2,sd),xlim =c(0,6))
points(y = truthCoefs[-1],x=c(0.5,1.5,2.5,3.5))
# Now we need to figure out how to do it all at once 

ExtMean = apply(BHolder[,-c(1,6)],2,mean)
# Plot with CIs. 
plotCI(1:4,ExtMean,ui = ExtMean + 1.96*apply(BHolder[,-c(1,6)],2,sd),li = ExtMean- 1.96*apply(BHolder[,-c(1,6)],2,sd),xlim =c(0,6))
points(y = truthCoefs[-1],x=c(0.5,1.5,2.5,3.5))
# Now we need to figure out how to do it all at once 

ui1 =  GTMean + 1.96*apply(GTHolder,2,sd)
ui2 =  GMean + 1.96*apply(GHolder,2,sd)
ui3 =  ExtMean + 1.96*apply(ExtHolder,2,sd)
li1 =  GTMean - 1.96*apply(GTHolder,2,sd)
li2 =  GMean - 1.96*apply(GHolder,2,sd)
li3 =  ExtMean - 1.96*apply(ExtHolder,2,sd)
MeanHolder = c(GTMean,GMean,ExtMean)
names(MeanHolder) = rep(c("PB","GM","GME"),7)
plotCI(c(2,6,10,14,18,22,26,3,7,11,15,19,23,27,4,8,12,16,20,24,28),MeanHolder,ui = c(ui1,ui2,ui3),li = c(li1,li2,li3),xlim = c(1,30),xlab ="",ylab = "Point Estimate")
points(x = c(1,5,9,13,17,21,25), y = truthCoefs[-1],pch = "-")
# These are good plots, but they are a little damning. Let's look at the boxplots again. 

PlotHolder = cbind(GTHolder,GHolder,GMEHolder)
colnames(PlotHolder) = rep(c("PB","GM","BL"),each=7)
boxplot(PlotHolder,at=c(2,6,10,14,18,22,26,3,7,11,15,19,23,27,4,8,12,16,20,24,28), xlim = c(1,30))#,col=rep(c("red","white","blue"),each=7))
points( x= c(1,5,9,13,17,21,25), y = rep(truthCoefs[-1]), pch = 8, lwd = 2)
axis(side = 3, at=c(2,6,10,14,18,22,26), xlab = "Parameter",labels=c("Math","SexM","EthB","EthH","Math:Male","Male:EthB","Male:EthH"))
#abline( v = c(3.5,6.5,9.5,12.5,15.5,18.5,21.5))
legend( "topright", c("PB","GM","BLASE"),fill = c("white","red","blue"), cex = .5)


