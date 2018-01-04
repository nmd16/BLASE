## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##  Load the starting space  ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ## 

setwd("N:/Simulations/LargeEthnic")

load("N:/Simulations/LargeEthnic/TruthData.RData")

truthModel = lm(mathscal2011 ~ Cmathscal * sex2011 + sex2011 * ethnic2011, data = Data)
truthCoefs = truthModel$coefficients


##%%%%%%%%%%%%%%%%%%%%
## Load the BLASE Run

FullDataRunThetaOut <- read.csv("N:/Simulations/LargeEthnic/LargeEthnicRunThetaOut.txt", quote="\"",header=F)

p= 9

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

#GME <- read.csv("N:/Simulations/LargeEthnic/ErrorRun/MoreError/MoreErrorGMThetaOut.txt", quote="\"",header=F)

p= 9

Ext <- GM

Ext <- as.matrix( Ext)
Ext <- strsplit( Ext, " ")
GM <- unlist(Ext)
GM <- matrix( GM, ncol = p, byrow = T )

ExtHold <- GM

GM <-as.numeric(ExtHold)
GM <- matrix( GM, ncol = p, nrow = nrow(ExtHold), byrow = F )
#GM = GME


N = 2*N

## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##         Boxplot           ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

GHolderBox = (G.outlist$out.coefs[-c(1:1000),-c(1,p)] - truthCoefs[-1])/truthCoefs[-1]
boxplot(GHolderBox)
boxplot(G.outlist$out.coefs[-c(1:1000),-c(1,p)])#,names=c("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"))
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

boxplot(GTest$out.coefs[-c(1:100),-c(1,10)],names=c("Math","YearPost1996","EthB","EthH","EthI","EthM","EthW","Male"))
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

# Okay, so we can now see the ethnicity coefficients are messy. Good :) Let's hope that our method can do better! 
EHolderBox = (ExtNew2[-c(1:200),-c(1,9)] - truthCoefs[-1])/truthCoefs[-1]
boxplot(EHolderBox)
abline(h=0)
boxplot(ExtNew2[-c(1:200),-c(1,20)])#,names=c("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"))
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

range.E = nrow(ExtNew)- 50
range.E = length(G.outlist$out.coefs[,1])- 200

boxplot(cbind(G.outlist$out.coefs[-c(1:1404),-c(1,p)],ExtNew2[-c(1:500),-c(1,p)]),at=c(1,3,5,7,9,11,13,2,4,6,8,10,12,14) ,col=rep(c("red","blue"),each=7)) #,names=rep(c ("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"),2)
points( rep(truthCoefs[-1],each=2), pch = 19, col = "green", lwd = 2)
abline( v = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5,30.5,32.5,34.5))

ExtHolder <- (ExtNew2[-c(1:500),-c(1,p)] - truthCoefs[-1])/truthCoefs[-1]
GHolder <- (GTest$out.coefs[-c(1:500,2822:4576),-c(1,p)] - truthCoefs[-1])/truthCoefs[-1]

boxplot(cbind(GHolder,ExtHolder),at=c(1,3,5,7,9,11,13, 2,4,6,8,10,12,14) ,col=rep(c("red","blue"),each=(p-2))) #,names=rep(c ("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"),2)
#points( rep(truthCoefs[-1],each=2), pch = 19, col = "green", lwd = 2)
abline( v = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5,30.5,32.5,34.5))
abline(h=0,lty=2)

which(ExtHolder[,3]>50)  
hist(sort(ExtHolder[,3]),breaks=50)
which(ExtHolder[,4]>50)
hist(sort(ExtHolder[,4]),breaks=50)

which(ExtHolder[,7]>50)

summary(GHolder[,1])
summary(ExtHolder[,1])

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Math11Error = File1[in.error,1]
summary(Math11Error)

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
plot( ExtNew[-c(1:burnoff),1], type = "l")
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
# Now, plot the absolute differences between GM and BLASE. 
# The x-axis represents the ranks. 
pdf("N:/Simulations/SimulationOutput/Figures/Boxplot_ByRank.pdf")
boxplot(GMvsBLASE, at = order(tranks),names=order(tranks),xlab = "Ranks:T-values",ylab = "Abs Diff", main =" Gain in Absolute Percent Difference: BLASE vs.GM")
abline(h = 0 ,lty = 2)
dev.off()

G.mean = apply(GHolder,2,mean)
Ext.mean = apply(ExtHolder,2,mean)

# Let's do this same thing with just the means 
pdf("Gains_BLASEvsGM_byRank.pdf")
plot(y = (abs(G.mean) - abs(Ext.mean))[order(tranks)], x = 1:(p-2),xlab = "Ranks:T-values",ylab = "Abs Diff", main =" Gain in Absolute Percent Difference: BLASE vs.GM")
abline(h = 0 ,lty = 2)
dev.off()

# The problem with this plot is that it does not take into account the fact that the differences are relatively tiny. We can scale them.

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
save(LAM, file = "GoutLambda.txt")
LAM <- G.outlist$out.coefs
save(LAM, file = "GoutCoefs.txt")
LAM <- G.outlist$out.coefs2
save(LAM, file = "GoutCoefs2.txt")

# GM - Imputed and Non Imputed 
holderall = apply( G.outlist$out.lambda,1, function(x) sum(x[1:nrow.F1]==1:nrow.F1))
summary(holderall[-c(1:100)]/nrow.F1)


# BLASE - Imputed and Non Imputed
FullDataRunAllLambdaPerc <- read.table("N:/Simulations/NewDataRun/NewDataRunAllLambdaPerc.txt", quote="\"")
summary(FullDataRunAllLambdaPerc[-c(1:100),])
plot(FullDataRunAllLambdaPerc[-c(1:1740),],type="l")

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

# GTest does not change 

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
# Yeah, we are just drawing from the prior here. 

AcceptOut <- read.table("N:/Simulations/DataSet/FullDataRun/FullDataRunAcceptOut.txt", quote="\"")

plot( as.matrix(AcceptOut[-c(1:10),]), type= "l")

SizeOut <- read.table("N:/Simulations/DataSet/FullDataRun/FullDataRunSizeOut.txt", quote="\"")

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

# Now, we can see how many records in the incorrect block 

## %%%%%%%%%%%%%%%%%%%%%%%%%%%
#### CI PLots #####
###################

meanVec = c( mean(GHolder[,1]),mean(BHolder[,1]), mean(GTest[-c(1:10),2]),mean(ExtNew2[-c(1:100),2]), mean(GTest[-c(1:10),3]),mean(ExtNew3[-c(1:100),3]), mean(GTest[-c(1:10),4]),mean(ExtNew2[-c(1:100),4]),mean(GTest[-c(1:10),5]),mean(ExtNew2[-c(1:100),5]),mean(GTest[-c(1:10),6]),mean(ExtNew2[-c(1:100),6]),mean(GTest[-c(1:10),7]),mean(ExtNew2[-c(1:100),7]),mean(GTest[-c(1:10),8]),mean(ExtNew2[-c(1:100),8]))

lil  = c( quantile( GHolder[,1],0.025) , quantile( BHolder[,1],0.025),quantile( GTest[-c(1:100,3879:5000),2],0.025) , quantile( ExtNew2[-c(1:100),2],0.025),quantile( GTest[-c(1:100,3879:5000),3],0.025) , quantile( ExtNew2[-c(1:100),3],0.025),quantile( GTest[-c(1:100,3879:5000),4],0.025) , quantile( ExtNew2[-c(1:100),4],0.025),quantile( GTest[-c(1:100,3879:5000),5],0.025) , quantile( ExtNew2[-c(1:100),5],0.025),quantile( GTest[-c(1:100,3879:5000),6],0.025) , quantile( ExtNew2[-c(1:100),6],0.025),quantile( GTest[-c(1:100,3879:5000),7],0.025) , quantile( ExtNew2[-c(1:100),7],0.025),quantile( GTest[-c(1:100,3879:5000),8],0.025) , quantile( ExtNew2[-c(1:100),8],0.025))

uil  = c( quantile( GHolder[,1],0.975) , quantile(BHolder[,1],0.975),quantile( GTest[-c(1:100,3879:5000),2],0.975) , quantile( ExtNew2[-c(1:100),2],0.975),quantile( GTest[-c(1:100,3879:5000),3],0.975) , quantile( ExtNew2[-c(1:100),3],0.975),quantile( GTest[-c(1:100,3879:5000),4],0.975) , quantile( ExtNew2[-c(1:100),4],0.975),quantile( GTest[-c(1:100,3879:5000),5],0.975) , quantile( ExtNew2[-c(1:100),5],0.975),quantile( GTest[-c(1:100,3879:5000),6],0.975) , quantile( ExtNew2[-c(1:100),6],0.975),quantile( GTest[-c(1:100,3879:5000),7],0.975) , quantile( ExtNew2[-c(1:100),7],0.975),quantile( GTest[-c(1:100,3879:5000),8],0.975) , quantile( ExtNew2[-c(1:100),8],0.975))

GHolder = GTest[-c(1:100,3879:5000),-c(p)]
GHolder[,1] = GHolder[,1] - mean(File2011[,1])
BHolder = ExtNew2[-c(1:100),-c(p)]
BHolder[,1] = BHolder[,1] - mean(File2011[,1])

# Plot with CIs.
pdf('~/AlEthnicWithCI.pdf')
boxplot(cbind(GHolder,BHolder),at=c(1,4,7,10,13,16,19,22, 2,5,8,11,14,17,20,23) ,names=rep(c("GM","B"),each=8),border="white")
plotCI(x=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23),meanVec,ui = uil,li = lil,add=T)
points(y = truthCoefs,x=c(0.5,3.5,6.5,9.5,12.5,15.5,18.5,21.5),pch="x")
mtext(c("Int","Math","Male", "EthB", "EthH", "Male:Math", "Male:EthB","Male:EthH"),1,line=3,at = c(1.5,4.5,7.5,10.5,13.5,16.5,19.5,22.5))
dev.off()

meanVec = c( mean(GHolder[,1]),mean(BHolder[,1]), mean(GM[-c(1:10),2]),mean(ExtNew2[-c(1:100),2]), mean(GM[-c(1:10),3]),mean(ExtNew3[-c(1:100),3]), mean(GM[-c(1:10),4]),mean(ExtNew2[-c(1:100),4]),mean(GM[-c(1:10),5]),mean(ExtNew2[-c(1:100),5]),mean(GM[-c(1:10),6]),mean(ExtNew2[-c(1:100),6]),mean(GM[-c(1:10),7]),mean(ExtNew2[-c(1:100),7]),mean(GM[-c(1:10),8]),mean(ExtNew2[-c(1:100),8]))

lil  = c( quantile( GHolder[,1],0.025) , quantile( BHolder[,1],0.025),quantile( GM[-c(1:100,3879:5000),2],0.025) , quantile( ExtNew2[-c(1:100),2],0.025),quantile( GM[-c(1:100,3879:5000),3],0.025) , quantile( ExtNew2[-c(1:100),3],0.025),quantile( GM[-c(1:100,3879:5000),4],0.025) , quantile( ExtNew2[-c(1:100),4],0.025),quantile( GM[-c(1:100,3879:5000),5],0.025) , quantile( ExtNew2[-c(1:100),5],0.025),quantile( GM[-c(1:100,3879:5000),6],0.025) , quantile( ExtNew2[-c(1:100),6],0.025),quantile( GM[-c(1:100,3879:5000),7],0.025) , quantile( ExtNew2[-c(1:100),7],0.025),quantile( GM[-c(1:100,3879:5000),8],0.025) , quantile( ExtNew2[-c(1:100),8],0.025))

uil  = c( quantile( GHolder[,1],0.975) , quantile(BHolder[,1],0.975),quantile( GM[-c(1:100,3879:5000),2],0.975) , quantile( ExtNew2[-c(1:100),2],0.975),quantile( GM[-c(1:100,3879:5000),3],0.975) , quantile( ExtNew2[-c(1:100),3],0.975),quantile( GM[-c(1:100,3879:5000),4],0.975) , quantile( ExtNew2[-c(1:100),4],0.975),quantile( GM[-c(1:100,3879:5000),5],0.975) , quantile( ExtNew2[-c(1:100),5],0.975),quantile( GM[-c(1:100,3879:5000),6],0.975) , quantile( ExtNew2[-c(1:100),6],0.975),quantile( GM[-c(1:100,3879:5000),7],0.975) , quantile( ExtNew2[-c(1:100),7],0.975),quantile( GM[-c(1:100,3879:5000),8],0.975) , quantile( ExtNew2[-c(1:100),8],0.975))

GHolder = GM[-c(1:100,3879:5000),-c(p)]
GHolder[,1] = GHolder[,1] - mean(File2011[,1])
BHolder = ExtNew2[-c(1:100),-c(p)]
BHolder[,1] = BHolder[,1] - mean(File2011[,1])

# Plot with CIs.
pdf('~/AlEthnicWithCI.pdf')
boxplot(cbind(GHolder,BHolder),at=c(1,4,7,10,13,16,19,22, 2,5,8,11,14,17,20,23) ,names=rep(c("GM","B"),each=8),border="white")
plotCI(x=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23),meanVec,ui = uil,li = lil,add=T)
points(y = truthCoefs,x=c(0.5,3.5,6.5,9.5,12.5,15.5,18.5,21.5),pch="x")
mtext(c("Int","Math","Male", "EthB", "EthH", "Male:Math", "Male:EthB","Male:EthH"),1,line=3,at = c(1.5,4.5,7.5,10.5,13.5,16.5,19.5,22.5))
dev.off()

