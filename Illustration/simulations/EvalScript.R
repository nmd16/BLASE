#Eval Full Data 

load("C:/Users/nmd16/Desktop/Final Runs/RunAll/PostGutman.RData")

setwd("C:/Users/nmd16/Desktop/Final Runs/RunAll")

load("C:/Users/nmd16/Desktop/Final Runs/RunAll/GOutRunAllTest.RData")
GTest <- G.outlist
load("C:/Users/nmd16/Desktop/Final Runs/RunAll/GOutRunAll.RData")

ExtNew <- read.table("C:/Users/nmd16/Desktop/Final Runs/RunAll/RunAllLongThetaOut.txt", quote="\"")
Ext <- as.matrix( ExtNew, ncol = 10)
ExtNew <- matrix( Ext, ncol = 10, nrow = length(ExtNew)/10, byrow = T )

## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##    Posterior Summaries    ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ## 

posterior_summaries_neat(G.outlist$out.coefs,c("Intercept", "Math","Year Post 1996", "Ethnic A", "Ethnic B", "Ethnic H", "Ethnic M","Sex Male"),truthCoefs)

posterior_summaries_neat(GTest$out.coefs,c("Intercept", "Math","Year Post 1996", "Ethnic A", "Ethnic B", "Ethnic H", "Ethnic I", "Ethnic M","Sex Male"),truthCoefs)

posterior_summaries_neat(ExtNew,c("Intercept", "Math","Year Post 1996", "Ethnic A", "Ethnic B", "Ethnic H", "Ethnic M","Sex Male"),truthCoefs)

## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##         Boxplot           ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

boxplot(G.outlist$out.coefs[,-c(1,10)],names=c("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"))
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

boxplot(GTest$out.coefs[,-c(1,10)],names=c("Math","YearPost1996","EthB","EthH","EthI","EthM","EthW","Male"))
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

boxplot(ExtNew[,-c(1,10)],names=c("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"))
points( truthCoefs[-1], pch = 19, col = "green", lwd = 2)

lExt = nrow(ExtNew)-50

boxplot(cbind(G.outlist$out.coefs[c(501:1000),-c(1,10)],ExtNew[-c(1:295),-c(1,10)]),at=c(1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16),names=rep(c ("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"),2),col=rep(c("red","blue"),each=8))
points( rep(truthCoefs[-1],each=2), pch = 19, col = "green", lwd = 2)
abline( v = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5))

boxplot(cbind(G.outlist$out.coefs[c(501:(1000)),-c(1,10)],GTest$out.coefs[c(501:(400+600)),-c(1,10)],ExtNew[-c(1:295),-c(1,10)]),at=c(1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,23,3,6,9,12,15,18,21,24),names=rep(c ("Math","YearPost1996","EthA","EthB","EthH","EthI","EthM","Male"),3),col=rep(c("red","white","blue"),each=8))
points( rep(truthCoefs[-1],each=3), pch = 19, col = "green", lwd = 2)
abline( v = c(3.5,6.5,9.5,12.5,15.5,18.5,21.5))
legend( "topright", c("GOnly","Matched","Dalzell"),fill = c("red","white","blue"), cex = .5)

## %%%%%%%%%%%%%%%%%%%%%%%% ## 
##   Trace Plots : Check    ##
## %%%%%%%%%%%%%%%%%%%%%%%% ## 

plot( G.outlist$out.coefs[,1],type="l")
plot( G.outlist$out.coefs[,2],type="l")
plot( G.outlist$out.coefs[,3],type="l")
plot( G.outlist$out.coefs[,4],type="l")
plot( G.outlist$out.coefs[,5],type="l")
plot( G.outlist$out.coefs[,6],type="l")
plot( G.outlist$out.coefs[,7],type="l")
plot( G.outlist$out.coefs[,8],type="l")
plot( G.outlist$out.coefs[,9],type="l")

plot( ExtNew[,1], type = "l")
abline(h=truthCoefs[1])
plot( ExtNew[,2], type = "l")
abline(h=truthCoefs[2])
plot( ExtNew[,3], type = "l")
abline(h=truthCoefs[3])
plot( ExtNew[,4], type = "l")
abline(h=truthCoefs[4])
plot( ExtNew[,5], type = "l")
abline(h=truthCoefs[5])
plot( ExtNew[,6], type = "l")
abline(h=truthCoefs[6])
plot( ExtNew[,7], type = "l")
abline(h=truthCoefs[7])
plot( ExtNew[,8], type = "l")
abline(h=truthCoefs[8])
plot( ExtNew[,9], type = "l")
abline(h=truthCoefs[9])

plot( GTest$out.coefs[,1],type="l")
plot( GTest$out.coefs[,2],type="l")
plot( GTest$out.coefs[,3],type="l")
plot( GTest$out.coefs[,4],type="l")
plot( GTest$out.coefs[,5],type="l")
plot( GTest$out.coefs[,6],type="l")
plot( GTest$out.coefs[,7],type="l")
plot( GTest$out.coefs[,8],type="l")
plot( GTest$out.coefs[,9],type="l")

new.burn = 200

make_traceplotsALL( ExtNew[-c(1:new.burn),1], truthS[1], "Intercept",G.outlist$out.coefs[c(401:(400+lExt)),1],option= "compare",GTest$out.coefs[c(401:(400+lExt)),1])

make_traceplotsALL( ExtNew[-c(1:new.burn),2], truthS[2], "Math",G.outlist$out.coefs[c(401:476),2],option= "compare",GTest$out.coefs[c(401:476),2])

make_traceplotsALL( ExtNew[-c(1:new.burn),3], truthS[3], "YearPost1996",G.outlist$out.coefs[c(401:476),3],option= "compare",GTest$out.coefs[c(401:476),3])

make_traceplotsALL( ExtNew[-c(1:new.burn),4], truthS[4], "Eth:A",G.outlist$out.coefs[-c(1:200),4],option= "compare",GTest$out.coefs[-c(1:200),4])

make_traceplotsALL( ExtNew[-c(1:new.burn),5], truthS[5], "Eth:B",G.outlist$out.coefs[c(401:476),5],option= "compare",GTest$out.coefs[c(401:476),5])

make_traceplotsALL( ExtNew[-c(1:new.burn),6], truthS[6], "Eth:H",G.outlist$out.coefs[c(401:476),6],option= "compare",GTest$out.coefs[c(401:476),6])

make_traceplotsALL( ExtNew[-c(1:new.burn),7], truthS[7], "Eth:M",G.outlist$out.coefs[c(401:476),7],option= "compare",GTest$out.coefs[c(401:476),7])

make_traceplotsALL( ExtNew[-c(1:new.burn),8], truthS[8], "Male",G.outlist$out.coefs[c(401:476),8],option= "compare",GTest$out.coefs[c(401:476),8])

## %%%%%%%%%%%%%%%%%%%%%%%% ## 
##       Scatter Plot       ##
## %%%%%%%%%%%%%%%%%%%%%%%% ##

GT.mean   = apply(GTest$out.coefs[-c(1:new.burn),-c(1,10)],2, mean)
G.mean   = apply( G.outlist$out.coefs[-c(1:new.burn),-c(1,10)],2,mean)
Ext.mean = apply( ExtNew[-c(1:new.burn),-c(1,10)],2,mean)

G.diff  = (G.mean - truthS[-1])/abs(truthS[-1])
Ext.diff = (Ext.mean - truthS[-1])/abs(truthS[-1])

plot( G.diff~ Ext.diff, pch = c("M","Y","A", "B","H","I","E","S"), main="Scatter: Gutman/Extension Versus Truth", xlab = "(1/|Matched|)* |Posterior Mean Extension - Matched|" , ylab = "(1/|Matched|)* |Posterior Mean Gutman - Matched|",xlim = c(-1,1),ylim = c(-1,1))
abline( a=0, b = 1)
abline( h = 0, lty = 2)
abline( v = 0, lty = 2)

# Shade the areas that Gutman out performs the extension
cord.x <- 0
cord.y <- 0

cord.x <-c(cord.x, -1.1)
cord.y <-c(cord.y, 0)

cord.x <-c(cord.x, -1.1)
cord.y <-c(cord.y, -1.1)

polygon(cord.x, cord.y, col=rgb(1,0,0,0.25))

cord.x <- 0
cord.y <- 0

cord.x <-c(cord.x, 1.1)
cord.y <-c(cord.y, 0)

cord.x <-c(cord.x, 1.1)
cord.y <-c(cord.y, 1.1)

polygon(cord.x, cord.y, col=rgb(1,0,0,0.25))

# Now see how we do compared to the best Gutman can do, rather than the truth.

G.diff.VT   = (G.mean - GT.mean)/abs(GT.mean)
Ext.diff.VT = (Ext.mean - GT.mean)/abs(GT.mean)

plot( G.diff.VT~ Ext.diff.VT, pch = c("M","Y","A", "B","H","I","E","S"),main="Scatter: Gutman/Extension Vs Matched", xlab = "(1/|GTest|)* |Posterior Mean Extension - GTest|" , ylab = "(1/|GTest|)* |Posterior Mean Gutman - GTest|",xlim = c(-1.2,1.2),ylim = c(-1.2,1.2))
abline( a=0, b = 1)
abline( h = 0, lty = 2)
abline( v = 0, lty = 2)

# Shade the areas that Gutman out performs the extension
cord.x <- 0
cord.y <- 0

cord.x <-c(cord.x, -1.2)
cord.y <-c(cord.y, 0)

cord.x <-c(cord.x, -1.2)
cord.y <-c(cord.y, -1.2)

polygon(cord.x, cord.y, col=rgb(1,0,0,0.25))

cord.x <- 0
cord.y <- 0

cord.x <-c(cord.x, 1.3)
cord.y <-c(cord.y, 0)

cord.x <-c(cord.x, 1.3)
cord.y <-c(cord.y, 1.3)

polygon(cord.x, cord.y, col=rgb(1,0,0,0.25))

## %%%%%%%%%%%%%%%%%%%%%%%% ## 
##    Error Distribution    ##
## %%%%%%%%%%%%%%%%%%%%%%%% ##

table(File2012[in.error, "ethnic"],File2011[in.error,"ethnic"])

table(File2012[matches[,1],"ethnic"])/length(matches[,1])

table(File2012[,"ethnic"])/N


(abs(Ext.diff)-abs(G.diff))*100


## %%%%%%%%%%%%%%%%%%%%%%%%% ## 
##       Match Rate          ##
## %%%%%%%%%%%%%%%%%%%%%%%%% ##

# Gutman 
holder = apply( G.outlist$out.lambda,1, function(x) sum(x[1:N]==1:N))
summary(holder/N)

# Gutman Test 
holderT = apply( GTest$out.lambda,1, function(x) sum(x[1:N]==1:N))
summary(holderT/N)

# Extension
LongRunLambdaOut <- read.table("RunAllLongLambdaOut.txt", quote="\"")
summary(LongRunLambdaOut[-c(1:new.burn),]/N)

MatchRate <- read.table("RunAllLongAllLambdaPerc.txt", quote="\"")
summary(MatchRate[-c(1:new.burn),])
plot(as.matrix(MatchRate[-c(1:10),]),type="l")

Size <- read.table("RunAllLongSizeOut.txt", quote="\"")
plot(as.matrix(Size[-c(1:10),]),type="l",ylab = "Number of Record Pairs")
abline(h = N)

AcceptOut <- read.table("RunAllLongAcceptOut.txt", quote="\"")
plot(as.matrix(AcceptOut[-c(1:new.burn),]),type="p",ylab = "Accepted Moves per Iteration")
abline(h = N)

CheckSpace = cbind( AcceptOut, Size)

RunAllLongGammaOut1 <- read.table("C:/Users/nmd16/Desktop/Final Runs/RunAll/RunAllLongGammaOut1.txt", quote="\"")
plot(as.matrix(RunAllLongGammaOut1[-c(1:50),]),type="l",ylab = "Error Percentage")
abline(h = 1)

# %%%%%%%%%%%%%%%%%%%%%%%%% #
#  Overlapping Histograms   #
# %%%%%%%%%%%%%%%%%%%%%%%%% # 

hist( G.outlist$out.coefs[ -c(1:new.burn), 2 ] - truthS[2], main = "Overlapping Histograms: Math",xlim = c(-0.005, 0.014) )
hist( ExtNew[-c(1:new.burn),2]-truthS[2], col = rgb(0.8,0.8,0.8,0.5), add = T)
abline( v = 0 )
legend( "topleft", c("Gutman","Dalzell"), col = c("white",rgb(0.8,0.8,0.8,0.5)),fill=c("white",rgb(0.8,0.8,0.8,0.5)))

(G.mean[1] - truthS[2])/truthS[2]*100
(Ext.mean[1] - truthS[2])/truthS[2]*100

Gpart = G.outlist$out.coefs[-c(1:new.burn),3]-truthS[3]
Epart = ExtNew[-c(1:295),3] - truthS[3]
min.part = min(min(Gpart),min(Epart))
max.part = max(max(Gpart),max(Epart))
name.hist = "Overlapping Histograms: Year Post 1996"

hist( Gpart, main =  name.hist,xlim = c(min.part,max.part))
hist( Epart, col = rgb(0.8,0.8,0.8,0.5), add = T)
abline( v = 0)
legend( "topleft", c("Gutman","Dalzell"), col = c("white",rgb(0.8,0.8,0.8,0.5)),fill=c("white",rgb(0.8,0.8,0.8,0.5)))

Gpart = G.outlist$out.coefs[-c(1:new.burn),4]-truthS[4]
Epart = ExtNew[-c(1:new.burn),4] - truthS[4]
min.part = min(min(Gpart),min(Epart))
max.part = max(max(Gpart),max(Epart))
name.hist = "Overlapping Histograms: Eth A"

hist( Gpart, main =  name.hist,xlim = c(min.part,max.part))
hist( Epart, col = rgb(0.8,0.8,0.8,0.5), add = T)
abline( v = 0)
legend( "topleft", c("Gutman","Dalzell"), col = c("white",rgb(0.8,0.8,0.8,0.5)),fill=c("white",rgb(0.8,0.8,0.8,0.5)))

Gpart = G.outlist$out.coefs[-c(1:new.burn),5]-truthS[5]
Epart = ExtNew[-c(1:new.burn),5] - truthS[5]
min.part = min(min(Gpart),min(Epart))
max.part = max(max(Gpart),max(Epart))
name.hist = "Overlapping Histograms: Eth B"

hist( Gpart, main =  name.hist,xlim = c(min.part,max.part),breaks=20)
hist( Epart, col = rgb(0.8,0.8,0.8,0.5), add = T)
abline( v = 0)
legend( "topleft", c("Gutman","Dalzell"), col = c("white",rgb(0.8,0.8,0.8,0.5)),fill=c("white",rgb(0.8,0.8,0.8,0.5)))

Gpart = G.outlist$out.coefs[-c(1:new.burn),6]-truthS[6]
Epart = ExtNew[-c(1:new.burn),6] - truthS[6]
min.part = min(min(Gpart),min(Epart))
max.part = max(max(Gpart),max(Epart))
name.hist = "Overlapping Histograms: Eth H"

hist( Gpart, main =  name.hist,xlim = c(min.part,max.part),breaks=20)
hist( Epart, col = rgb(0.8,0.8,0.8,0.5), add = T)
abline( v = 0)
legend( "topleft", c("Gutman","Dalzell"), col = c("white",rgb(0.8,0.8,0.8,0.5)),fill=c("white",rgb(0.8,0.8,0.8,0.5)))

Gpart = G.outlist$out.coefs[-c(1:new.burn),7]-truthS[7]
Epart = ExtNew[-c(1:new.burn),7] - truthS[7]
min.part = min(min(Gpart),min(Epart))
max.part = max(max(Gpart),max(Epart))
name.hist = "Overlapping Histograms: Eth I"

hist( Gpart, main =  name.hist,xlim = c(min.part,max.part),ylim = c(1,200))
hist( Epart, col = rgb(0.8,0.8,0.8,0.5), add = T)
abline( v = 0)
legend( "topleft", c("Gutman","Dalzell"), col = c("white",rgb(0.8,0.8,0.8,0.5)),fill=c("white",rgb(0.8,0.8,0.8,0.5)))

Gpart = G.outlist$out.coefs[-c(1:new.burn),8]-truthS[8]
Epart = ExtNew[-c(1:new.burn),8] - truthS[8]
min.part = min(min(Gpart),min(Epart))
max.part = max(max(Gpart),max(Epart))
name.hist = "Overlapping Histograms: Eth M"

hist( Gpart, main =  name.hist,xlim = c(min.part,max.part), ylim = c(0, 200))
hist( Epart, col = rgb(0.8,0.8,0.8,0.5), add = T)
abline( v = 0)
legend( "topleft", c("Gutman","Dalzell"), col = c("white",rgb(0.8,0.8,0.8,0.5)),fill=c("white",rgb(0.8,0.8,0.8,0.5)))

Gpart = G.outlist$out.coefs[-c(1:new.burn),9]-truthS[9]
Epart = ExtNew[-c(1:295),9] - truthS[9]
min.part = min(min(Gpart),min(Epart))
max.part = max(max(Gpart),max(Epart))
name.hist = "Overlapping Histograms: Sex M"

hist( Gpart, main =  name.hist,xlim = c(min.part,max.part))
hist( Epart, col = rgb(0.8,0.8,0.8,0.5), add = T)
abline( v = 0)
legend( "topleft", c("Gutman","Dalzell"), col = c("white",rgb(0.8,0.8,0.8,0.5)),fill=c("white",rgb(0.8,0.8,0.8,0.5)))


