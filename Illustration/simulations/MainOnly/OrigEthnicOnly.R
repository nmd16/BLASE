## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##  Load the starting space  ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ## 

setwd("N:/Simulations/LargeEthnic")

load("N:/Simulations/LargeEthnic/TruthData.RData")

truthModel = lm( File2011[,1] ~ Cmathscal + ethnic+sex, data = File2011)
truthCoefs = truthModel$coefficients

##%%%%%%%%%%%%%%%%%%%%
## Load the BLASE Run

FullDataRunThetaOut <- read.csv("N:/Simulations/MainOnly/MainAllEthnicRun/MainErrorOrigEthnicRunThetaOut.txt", quote="\"",header=F)

p= 6

Ext <- FullDataRunThetaOut 

Ext <- as.matrix( Ext)
Ext <- strsplit( Ext, " ")
ExtNew3 <- unlist(Ext)
ExtNew3 <- matrix( ExtNew3, ncol = p, byrow = T )

ExtHold <- ExtNew3

ExtNew3 <-as.numeric(ExtHold)
ExtNew3 <- matrix( ExtNew3, ncol = p, nrow = nrow(ExtHold), byrow = F )

##%%%%%%%%%%%%%%%%%%%%
## Load the GM Run

GTest <- read.csv("N:/Simulations/MainOnly/MainAllEthnicRun/GMainErrorOrigEthnicRunThetaOut.txt", quote="\"",header=F)

Ext <- GTest

Ext <- as.matrix( Ext)
Ext <- strsplit( Ext, " ")
GTest <- unlist(Ext)
GTest <- matrix( GTest, ncol = p, byrow = T )

ExtHold <- GTest

GTest <-as.numeric(ExtHold)
GTest <- matrix( GTest, ncol = p, nrow = nrow(ExtHold), byrow = F )

## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##         Boxplot           ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

GTest[,1] = GTest[,1] -mean(File2011[,1])
# Look at the main effects with random errors 
boxplot(cbind(GTest[-c(1:10),-c(1,p)],ExtNew3[-c(1:2436),-c(1,p)]),at=c(1,3,5,7,2,4,6,8) ,col=rep(c("red","blue"),each=4),names=c(rep(c("GM","B"),each=4))) #,names=rep(c ("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"),2)
points( rep(truthCoefs[-1],each=2), pch = 19, col = "green", lwd = 2)
abline( v = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5,30.5,32.5,34.5))

## %%%%%%%%%%%%%%%%%%%%%%%% ## 
##    Error Distribution    ##
## %%%%%%%%%%%%%%%%%%%%%%%% ##

setcheck = setdiff( type1seeds, which(File2012[,"ethnic"]!=File2011[,"ethnic"]) %in% in.error)

table(File2012[in.error, "ethnic"],File2011[in.error,"ethnic"])


table(File2012[type1seeds, "ethnic"],File2011[type1seeds,"ethnic"])/n.type1


## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##   Posterior for Gamma     ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

GammaOut <- read.table("C:/Users/nmd16/Desktop/Final Runs/Run1/Run1GammaOut1.txt", quote="\"")

summary(GammaOut)

AcceptOut <- read.table("N:/Simulations/LargeEthnic/MainShortRunAcceptOut.txt", quote="\"")

plot( as.matrix(AcceptOut[-c(1:10),]), type= "l")

SizeOut <- read.table("N:/Simulations/LargeEthnic/MainShortRunSizeOut.txt", quote="\"")

plot( as.matrix(SizeOut), type= "l")


##%%%%%%%%%%%%%%%
###    CIS  #####
## %%%%%%%%%%%%%%

# Explore creating a graph with CIs. 

myrepo = getOption("repos")
myrepo["CRAN"] = "http://archive.linux.duke.edu/cran/"
options(repos=myrepo)
install.packages("plotrix")
library(plotrix)
rm(myrepo)

meanVec = c( mean(GTest[-c(1:10),1]),mean(ExtNew3[-c(1:100),1]), mean(GTest[-c(1:10),2]),mean(ExtNew3[-c(1:100),2]), mean(GTest[-c(1:10),3]),mean(ExtNew3[-c(1:100),3]), mean(GTest[-c(1:10),4]),mean(ExtNew3[-c(1:100),4]),mean(GTest[-c(1:10),5]),mean(ExtNew3[-c(1:100),5]))

lil  = c( quantile( GTest[-c(1:10),1],0.025) , quantile( ExtNew3[-c(1:100),1],0.025),quantile( GTest[-c(1:10),2],0.025) , quantile( ExtNew3[-c(1:100),2],0.025),quantile( GTest[-c(1:10),3],0.025) , quantile( ExtNew3[-c(1:100),3],0.025),quantile( GTest[-c(1:10),4],0.025) , quantile( ExtNew3[-c(1:100),4],0.025),quantile( GTest[-c(1:10),5],0.025) , quantile( ExtNew3[-c(1:100),5],0.025) )

uil  = c( quantile( GTest[-c(1:10),1],0.975) , quantile( ExtNew3[-c(1:100),1],0.975),quantile( GTest[-c(1:10),2],0.975) , quantile( ExtNew3[-c(1:100),2],0.975),quantile( GTest[-c(1:10),3],0.975) , quantile( ExtNew3[-c(1:100),3],0.975),quantile( GTest[-c(1:10),4],0.975) , quantile( ExtNew3[-c(1:100),4],0.975),quantile( GTest[-c(1:10),5],0.975) , quantile( ExtNew3[-c(1:100),5],0.975) )

# Plot with CIs.
pdf('~/MainOnlyWithCI.pdf')
boxplot(cbind(GTest[-c(1:10),-c(p)],ExtNew3[-c(1:2436),-c(p)]),at=c(1,4,7,10,13,2,5,8,11,14),names=c(rep(c("GM","B"),each=5)),border="white")
plotCI(x=c(1,2,4,5,7,8,10,11,13,14),meanVec,ui = uil,li = lil,add=T)
points(y = truthCoefs2,x=c(0.5,3.5,6.5,9.5,12.5),pch="x")
mtext(c("Intercept","Math","EthB","EthH","Male"),1,line=3,at = c(1.5,4.5,7.5,10.5,13.5))
dev.off()

meanVec = c( mean(GTest[-c(1:10),1]),mean(ExtNew3[-c(1:100),1]), mean(GTest[-c(1:10),2]),mean(ExtNew3[-c(1:100),2]), mean(GTest[-c(1:10),3]),mean(ExtNew3[-c(1:100),3]), mean(GTest[-c(1:10),4]),mean(ExtNew3[-c(1:100),4]),mean(GTest[-c(1:10),5]),mean(ExtNew3[-c(1:100),5]))

lil  = c( quantile( GTest[-c(1:10),1],0.025) , quantile( ExtNew3[-c(1:100),1],0.025),quantile( GTest[-c(1:10),2],0.025) , quantile( ExtNew3[-c(1:100),2],0.025),quantile( GTest[-c(1:10),3],0.025) , quantile( ExtNew3[-c(1:100),3],0.025),quantile( GTest[-c(1:10),4],0.025) , quantile( ExtNew3[-c(1:100),4],0.025),quantile( GTest[-c(1:10),5],0.025) , quantile( ExtNew3[-c(1:100),5],0.025) )

uil  = c( quantile( GTest[-c(1:10),1],0.975) , quantile( ExtNew3[-c(1:100),1],0.975),quantile( GTest[-c(1:10),2],0.975) , quantile( ExtNew3[-c(1:100),2],0.975),quantile( GTest[-c(1:10),3],0.975) , quantile( ExtNew3[-c(1:100),3],0.975),quantile( GTest[-c(1:10),4],0.975) , quantile( ExtNew3[-c(1:100),4],0.975),quantile( GTest[-c(1:10),5],0.975) , quantile( ExtNew3[-c(1:100),5],0.975) )

# Plot with CIs.
pdf('~/MainOnlyWithCI.pdf')
boxplot(cbind(GTest[-c(1:10),-c(p)],ExtNew3[-c(1:2436),-c(p)]),at=c(1,4,7,10,13,2,5,8,11,14),names=c(rep(c("GM","B"),each=5)),border="white")
plotCI(x=c(1,2,4,5,7,8,10,11,13,14),meanVec,ui = uil,li = lil,add=T)
points(y = truthCoefs2,x=c(0.5,3.5,6.5,9.5,12.5),pch="x")
mtext(c("Intercept","Math","EthB","EthH","Male"),1,line=3,at = c(1.5,4.5,7.5,10.5,13.5))
dev.off()



