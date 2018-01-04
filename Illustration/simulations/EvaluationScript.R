## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##    Posterior Summaries    ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ## 

File2012True = File2012
File2012True[,"ethnic"] = File2011[,"ethnic"]
(lm(File2011[, 1] ~ Cmathscal + sex*ethnic + YearInd*ethnic, data = File2012))$coefficients
truthCoefs = lm(File2011[, 1] ~ File2012[,"Cmathscal"] + sex * ethnic + YearInd *ethnic, data = File2011)$coefficients

truthCoefs = lm(File2011[, 1] ~ File2012$Cmathscal + sex * ethnic + YearInd *ethnic, data = File2011)$coefficients

load("N:/Simulations/NewDataRun/PreBLASE.RData")

##%%%%%%%%%%%%%%%%%%%%
## Load the BLASE Run

FullDataRunThetaOut <- read.csv("N:/Simulations/NewDataRun/NewDataRunThetaOut.txt", quote="\"",header=F)

p= 20

Ext <- FullDataRunThetaOut 

Ext <- as.matrix( Ext)
Ext <- strsplit( Ext, " ")
ExtNew2 <- unlist(Ext)
ExtNew2 <- matrix( ExtNew2, ncol = p, byrow = T )

ExtHold <- ExtNew2

ExtNew2 <-as.numeric(ExtHold)


rm(FullDataRunThetaOut, ExtHold,Ext)
##%%%%%%%%%%%%%%%%%%%%
## Load the PB Run
load("N:/Simulations/DataSet/NewDataRun/GOutTestNewDataRun.RData")
load("N:/Simulations/DataSet/FullDataRun/GTest.RData")
GTest <-G.outlist

##%%%%%%%%%%%%%%%%%%%%
## Load the GM Run
load("N:/Simulations/NewDataRun/GOutNewDataRun.RData")

N = 2*N

posterior_summaries_neat(G.outlist$out.coefs,c("Intercept", "CMath","SexM","Ethnic A", "Ethnic B", "Ethnic H", "Ethnic I", "Ethnic M","Year Ind", "SexM:EthA","SexM:EthB","SexM:EthH","SexM:EthI","SexM:EthM","EthA:Year","EthB:Year","EthH:Year","EthI:Year","EthM:Year"),truthCoefs)

posterior_summaries_neat(GTest$out.coefs,c("Intercept", "Math","Year Post 1996", "Ethnic A", "Ethnic B", "Ethnic H", "Ethnic I", "Ethnic M","Sex Male"),truthCoefs)

posterior_summaries_neat(ExtNew2[-c(1:20),],c("Intercept", "Math","Year Post 1996", "Ethnic A", "Ethnic B", "Ethnic H", "Ethnic I", "Ethnic M","Sex Male"),truthCoefs)

## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##         Boxplot           ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

Gnew = G.outlist
GHolderBoxNew = GHolderBox

GHolderBox = (G.outlist$out.coefs[-c(1:716),-c(1,20)] - truthCoefs[-1])/truthCoefs[-1]
boxplot(GHolderBox)
boxplot(G.outlist$out.coefs[-c(1:1000),-c(1,20)])#,names=c("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"))
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

boxplot(GTest$out.coefs[-c(1:100),-c(1,10)],names=c("Math","YearPost1996","EthB","EthH","EthI","EthM","EthW","Male"))
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

EHolderBox = (ExtNew2[-c(1:200),-c(1,20)] - truthCoefs[-1])/truthCoefs[-1]
boxplot(EHolderBox)
abline(h=0)
boxplot(ExtNew2[-c(1:200),-c(1,20)])#,names=c("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"))
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

range.E = nrow(ExtNew)- 50
range.E = length(G.outlist$out.coefs[,1])- 200

boxplot(cbind(G.outlist$out.coefs[c(250:1132),-c(1,20)],ExtNew2[c(250:1132),-c(1,20)]),at=c(1,3,5,7,9,11,13,15,17,19,21, 23, 25, 27, 29,31,33, 35, 2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36) ,col=rep(c("red","blue"),each=18)) #,names=rep(c ("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"),2)
points( rep(truthCoefs[-1],each=2), pch = 19, col = "green", lwd = 2)
abline( v = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5,30.5,32.5,34.5))

ExtHolder <- (ExtNew2[-c(1:250),-c(1,20)] - truthCoefs[-1])/truthCoefs[-1]
GHolder <- (G.outlist$out.coefs[-c(1:766),-c(1,20)] - truthCoefs[-1])/truthCoefs[-1]

boxplot(cbind(GHolder,ExtHolder),at=c(1,3,5,7,9,11,13,15,17,19,21, 23, 25, 27, 29,31,33, 35, 2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36) ,col=rep(c("red","blue"),each=18)) #,names=rep(c ("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"),2)
#points( rep(truthCoefs[-1],each=2), pch = 19, col = "green", lwd = 2)
abline( v = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5,30.5,32.5,34.5))
abline(h=0,lty=2)

# Creating plots for the paper 

# Look just at Main Effects 
colnames(GHolder)= rep("PB",18)
colnames(ExtHolder)=rep("B",18)
pdf("Box_MainEffects.pdf")
boxplot(cbind(GHolder[,1:8],ExtHolder[,1:8]),at=c(1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16),ylab = "% Distance from Truth",main = "Main Effects by Method") 
abline( v = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5))
mtext(c("CMath","SexM","EthA","EthB","EthH", "EthI", "EthM", "YearPost"),1,line=3,at = c(1.5,3.5,5.5,7.5,9.5,11.5,13.5,15.5))
abline(h=0,lty=2)
dev.off()

# Look just at Interactions 
pdf("Box_Interactions.pdf")
boxplot(cbind(GHolder[,9:18],ExtHolder[,9:18]),,at=c(1,3,5,7,9,11,13,15,17,19,2,4,6,8,10,12,14,16,18,20),ylab = "% Distance from Truth",main = "Interactions") 
abline( v = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5))
mtext(c("SexM:ethA","SexM:ethB","SexM:EthH","SexM:EthI","SexM:EthM","EthA:Year","EthB:Year","EthH:Year","EthI:Year","EthM:Year"),1,line=3,at = c(1.5,3.5,5.5,7.5,9.5,11.5,13.5,15.5,17.5,19.5))
abline(h=0,lty=2)
dev.off()

AllHolder = cbind(GHolder[,1:8],ExtHolder[,1:8])

AllHolder = rbind(GHolder[,1:8],ExtHolder[,1:8]) 
AllHolder = cbind(AllHolder,c(rep("PB",900),rep("BLASE",900)))
AllHolder = data.frame(AllHolder)
AllHolder[,1] = as.numeric(levels(AllHolder[,1])[AllHolder[,1]])
AllHolder[,2] = as.numeric(levels(AllHolder[,2])[AllHolder[,2]])
AllHolder[,3] = as.numeric(levels(AllHolder[,3])[AllHolder[,3]])
AllHolder[,4] = as.numeric(levels(AllHolder[,4])[AllHolder[,4]])
AllHolder[,5] = as.numeric(levels(AllHolder[,5])[AllHolder[,5]])
AllHolder[,6] = as.numeric(levels(AllHolder[,6])[AllHolder[,6]])
AllHolder[,7] = as.numeric(levels(AllHolder[,7])[AllHolder[,7]])
AllHolder[,8] = as.numeric(levels(AllHolder[,8])[AllHolder[,8]])


bwplot(V9 ~ AllHolder[,1],data=AllHolder,xlab = "% Distance from Truth",main="Intercept")

bwplot(V9 ~ AllHolder[,2],data=AllHolder,xlab = "% Distance from Truth",main="Math")

GHolderNew = c(GHolder[,1:8])
EHolderNew = c(ExtHolder[,1:8])
AllHolderNew = c(GHolderNew,EHolderNew)
AllHolderNew = cbind(AllHolderNew, c(rep("PB",7200),rep("BLASE",7200)))
nameholder = names(truthCoefs)[1:8]
nameholder[2] = "CMath"
AllHolderNew = cbind(AllHolderNew,c((rep(nameholder, each = 900))))
AllHolderNew = data.frame(AllHolderNew)
AllHolderNew[,1] = as.numeric(levels(AllHolderNew[,1])[AllHolderNew[,1]])
colnames(AllHolderNew)=c("Value","Method","Parameter")

bwplot( ~Value|Parameter, data=AllHolderNew,xlab = "% Distance from Truth")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExtHolder <- (ExtNew2[-c(1:250),-c(1,20)] - truthCoefs[-1])/truthCoefs[-1]
GHolder <- (G.outlist$out.coefs[-c(1:766),-c(1,20)] - truthCoefs[-1])/truthCoefs[-1]
GTHolder <- (GTest$out.coefs[-c(1:250,1151:1961),-c(1,20)] - truthCoefs[-1])/truthCoefs[-1]
colnames(GTHolder)= NULL

colnames(GHolder) = c("CMath", "sexM", "ethA", "ethB", "ethH", "ethI", "ethM", "YearPost", "sexM:ethA", "sexM:ethB", "sexM:ethH", "sexM:ethI", "sexM:ethM",      "ethA:YearPost","ethB:YearPost", "ethH:YearPost","ethI:YearPost","ethM:YearPost")

par(las=3) # Always make vertical labels 
boxplot(cbind(GTHolder,GHolder,ExtHolder),at=c(1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54),col=rep(c("white","red","blue"),each=18))
#points( rep(truthCoefs[-1],each=3), pch = 19, col = "green", lwd = 2)
abline( v = c(3.5,6.5,9.5,12.5,15.5,18.5,21.5,24.5,27.5,30.5,33.5,36.5,39.5,42.5,45.5,48.5,51.5))
legend( "topright", c("PB","GM","BLASE"),fill = c("white","red","blue"), cex = .5)

# Look just at Math 
pdf("Boxplot_Math.pdf")
boxplot(cbind(PB = GTHolder[,1],GM =GHolder[,1],BLASE = ExtHolder[,1]),col=rep(c("red","white","blue"),each=1),main = "Perc Distance from Truth: Math")
abline(h = 0, lty = 2 )
dev.off()

pdf("Figures/Boxplot_Sex.pdf")
boxplot(cbind(PB = GTHolder[,2],GM =GHolder[,2],BLASE = ExtHolder[,2]),col=rep(c("red","white","blue"),each=1),main = "Perc Distance from Truth: Sex M")
abline(h = 0, lty = 2 )
dev.off()

pdf("Figures/Boxplot_EthABHI.pdf")
boxplot(cbind(GTHolder[,3:6],GHolder[,3:6],ExtHolder[,3:6]),at=c(1,4,7,10,2,5,8,11,3,6,9,12),col=rep(c("red","white","blue"),each=4),names= rep(c("PB","GM","BLASE"),each = 4),main = "Perc Distance from Truth: Eth A, B, H,I")
abline(v =c(3.5,6.5,9.5))
abline( h = 0, lty = 2)
dev.off()

pdf("Figures/Boxplot_7to10.pdf")
boxplot(cbind(GTHolder[,7:10],GHolder[,7:10],ExtHolder[,7:10]),at=c(1,4,7,10,2,5,8,11,3,6,9,12),col=rep(c("red","white","blue"),each=4),names= rep(c("PB","GM","BLASE"),each = 4),main = "Perc Distance from Truth: Eth M, Year 1996+,Sex M:Eth A, SexM:Eth B")
abline(v =c(3.5,6.5,9.5))
abline( h = 0, lty = 2)
dev.off()

pdf("Figures/Boxplot_11to14.pdf")
boxplot(cbind(GTHolder[,11:14],GHolder[,11:14],ExtHolder[,11:14]),at=c(1,4,7,10,2,5,8,11,3,6,9,12),col=rep(c("red","white","blue"),each=4),names= rep(c("PB","GM","BLASE"),each = 4),main = "Perc Distance from Truth:SexM: EthH,Sex M: EthI, SexM:EthM, EthA:1996+")
abline(v =c(3.5,6.5,9.5))
abline(h=0, lty = 2)
dev.off()

pdf("Figures/Boxplot_15to18.pdf")
boxplot(cbind(PB = GTHolder[,15:18],GM = GHolder[,15:18],BLASE = ExtHolder[,15:18]),at=c(1,4,7,10,2,5,8,11,3,6,9,12),col=rep(c("red","white","blue"),each=4),names= rep(c("PB","GM","BLASE"),each = 4),main = "Perc Distance from Truth: EthB:1996+, EthH:1996+, EthI:1996+, EthM:1996+")
abline(v =c(3.5,6.5,9.5))
abline( h = 0, lty = 2)
dev.off()

## %%%%%%%%%%%%%%%%%%%%%%%% ## 
##   Trace Plots : Check    ##
## %%%%%%%%%%%%%%%%%%%%%%%% ## 

burnoff = 10
plot( G.outlist$out.coefs[-c(1:burnoff),1],type="l")
abline(h = truthCoefs[1],col="green")
plot( G.outlist$out.coefs[-c(1:burnoff),2],type="l")
abline(h = truthCoefs[2],col="green")
# 2 = YES 
plot( G.outlist$out.coefs[-c(1:burnoff),3],type="l")
abline(h = truthCoefs[3],col="green")
plot( G.outlist$out.coefs[-c(1:burnoff),4],type="l")
abline(h = truthCoefs[4],col="green")
plot( G.outlist$out.coefs[,5],type="l")
abline(h = truthCoefs[5],col="green")
plot( G.outlist$out.coefs[,6],type="l")
abline(h = truthCoefs[6],col="green")
plot( G.outlist$out.coefs[,7],type="l")
abline(h = truthCoefs[7],col="green")
plot( G.outlist$out.coefs[,8],type="l")
abline(h = truthCoefs[8],col="green")
plot( G.outlist$out.coefs[,9],type="l")
abline(h = truthCoefs[9],col="green")
plot( G.outlist$out.coefs[,10],type="l")
abline(h = truthCoefs[10],col="green")
# 10 = YES 
plot( G.outlist$out.coefs[,11],type="l")
abline(h = truthCoefs[11],col="green")
plot( G.outlist$out.coefs[,12],type="l")
abline(h = truthCoefs[12],col="green")
# 12 = YES 
plot( G.outlist$out.coefs[,13],type="l")
abline(h = truthCoefs[13],col="green")
# 13 = YES
plot( G.outlist$out.coefs[,14],type="l")
abline(h = truthCoefs[14],col="green")
# 14 = YES 
plot( G.outlist$out.coefs[,15],type="l")
abline(h = truthCoefs[15],col="green")
# 15 = YES 
plot( G.outlist$out.coefs[,16],type="l")
abline(h = truthCoefs[16],col="green")
plot( G.outlist$out.coefs[,17],type="l")
abline(h = truthCoefs[17],col="green")
plot( G.outlist$out.coefs[,18],type="l")
abline(h = truthCoefs[18],col="green")


ExtNew = ExtNew2
plot( ExtNew[,1], type = "l")
abline(h=truthCoefs[1])
plot( ExtNew[,2], type = "l")
abline(h=truthCoefs[2])
plot( ExtNew[,3], type = "l")
abline(h=truthCoefs[3])
# 3 = YES
plot( ExtNew[,4], type = "l")
abline(h=truthCoefs[4])
plot( ExtNew[,5], type = "l")
abline(h=truthCoefs[5])
# 5 = YES
plot( ExtNew[,6], type = "l")
abline(h=truthCoefs[6])
plot( ExtNew[,7], type = "l")
abline(h=truthCoefs[7])
plot( ExtNew[,8], type = "l")
abline(h=truthCoefs[8])
plot( ExtNew[,9], type = "l")
abline(h=truthCoefs[9])
plot( ExtNew[,10], type = "l")
abline(h=truthCoefs[10])
# 10 = YES 
plot( ExtNew[,11], type = "l")
abline(h=truthCoefs[11])
# 11 = YES 
plot( ExtNew[,12], type = "l")
abline(h=truthCoefs[12])
# 12 = YES 
plot( ExtNew[,13], type = "l")
abline(h=truthCoefs[13])
# 13 = YES 
plot( ExtNew[,14], type = "l")
abline(h=truthCoefs[14])
# 14 = YES 
plot( ExtNew[,15], type = "l")
abline(h=truthCoefs[15])
# 15 = YES 
plot( ExtNew[,16], type = "l")
abline(h=truthCoefs[16])
plot( ExtNew[,17], type = "l")
abline(h=truthCoefs[17])
plot( ExtNew[,18], type = "l")
abline(h=truthCoefs[18])
# 18 = YES 


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
##       Scatter Plot       ##
## %%%%%%%%%%%%%%%%%%%%%%%% ##

GT.mean   = apply(GTest$out.coefs[-c(1:100),-c(1,20)],2, mean)
G.mean   = apply( G.outlist$out.coefs[-c(1:100),-c(1,20)],2,mean)
Ext.mean = apply( ExtNew2[-c(1:100),-c(1,20)],2,mean)

G.diff  = (G.mean - truthCoefs[-1])/abs(truthCoefs[-1])
Ext.diff = (Ext.mean - truthCoefs[-1])/abs(truthCoefs[-1])

pdf("Figures/Scatter_GMandBLASE.pdf")
plot( G.diff~ Ext.diff, xlab = "Perc Difference: BLASE" , pch= c("1","2","3","4","5","6","7","8","9","A","B","C","D","E","F","G","H","J","K","L"), ylab = "Perc Difference: GM",xlim = c(-10,10),ylim= c(-10,10),main = "GM and BLASE vs. Truth")
abline( a=0, b = 1)
abline( h = 0, lty = 2)
abline( v = 0, lty = 2)

# Shade the areas that Gutman out performs the extension
cord.x <- 0
cord.y <- 0

cord.x <-c(cord.x, -10)
cord.y <-c(cord.y, 0)

cord.x <-c(cord.x, -10.1)
cord.y <-c(cord.y, -10.1)

polygon(cord.x, cord.y, col=rgb(1,0,0,0.25))

cord.x <- 0
cord.y <- 0

cord.x <-c(cord.x, 10.1)
cord.y <-c(cord.y, 0)

cord.x <-c(cord.x, 10.1)
cord.y <-c(cord.y, 10.1)

polygon(cord.x, cord.y, col=rgb(1,0,0,0.25))
dev.off()

G.diff.VT   = abs(G.mean - GT.mean)/abs(GT.mean)
Ext.diff.VT = abs(Ext.mean - GT.mean)/abs(GT.mean)

plot( G.diff.VT~ Ext.diff.VT, xlab = "Perc Difference: BLASE" , pch= c("1","2","3","4","5","6","7","8","9","A","B","C","D","E","F","G","H","J","K","L"), ylab = "Perc Difference: GM",xlim = c(-10,10),ylim= c(-10,10),main="GM and BLASE vs. PB")
abline( a=0, b = 1)
abline( h = 0, lty = 2)
abline( v = 0, lty = 2)
 
# Shade the areas that Gutman out performs the extension
cord.x <- 0
cord.y <- 0

cord.x <-c(cord.x, -10.3)
cord.y <-c(cord.y, 0)

cord.x <-c(cord.x, -10.3)
cord.y <-c(cord.y, -10.3)

polygon(cord.x, cord.y, col=rgb(1,0,0,0.25))

cord.x <- 0
cord.y <- 0

cord.x <-c(cord.x, 10.3)
cord.y <-c(cord.y, 0)

cord.x <-c(cord.x, 10.3)
cord.y <-c(cord.y, 10.3)

polygon(cord.x, cord.y, col=rgb(1,0,0,0.25))

## %%%%%%%%%%%%%%%%%%%%%%%% ## 
##       T-values           ##
## %%%%%%%%%%%%%%%%%%%%%%%% ##

names(truthModel)

summary(truthModel)

tvalues = c(10356.569,376.999,5.627,1.950,-14.234,-10.7,-2.04,-4.586,-10.046,-1.570,-4.016,-1.566,-0.731,0.614,-0.136,3.709,4,1.20,0.384)
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

# Let's do this same thing with just the means 
pdf("Gains_BLASEvsGM_byRank.pdf")
plot(y = (abs(G.mean) - abs(Ext.mean))[order(tranks)], x = 1:18,xlab = "Ranks:T-values",ylab = "Abs Diff", main =" Gain in Absolute Percent Difference: BLASE vs.GM")
abline(h = 0 ,lty = 2)
dev.off()

plot(y = ((abs(G.mean) - abs(Ext.mean))/(abs(G.mean)))[order(tranks)], x = 1:18,xlab = "Ranks:T-values",ylab = "Abs Diff", main =" Gain in Absolute Percent Difference: BLASE vs.GM")
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

# Compare the match rates ONLY considering the original records

setwd("N:/Simulations/DataSet/FullDataRun")
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


