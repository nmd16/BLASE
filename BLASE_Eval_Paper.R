## %%%%%%%%%%%%%%%%%%% ##
##    Paper Metrics    ##
## %%%%%%%%%%%%%%%%%%% ## 

require(foreign)
ml <-read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
GenCoefs = lm(read~math+prog,data = ml)$coefficients

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### High Fault, High Seed, Diffuse  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

setwd("BLASERuns/HF_HS_D")

num.sim  = 100
newburn = 1000 
indices  = 1:num.sim
lastcoef = 5 

load("OutSpace.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### High Fault, High Seed, CP ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

setwd("BLASERuns/HF_HS_CP")

load("OutSpace.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### High Fault, High Seed, CA ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

setwd("BLASERuns/HF_HS_CA")

load("OutSpace.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%##
####     HSHF : PMR     ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Posterior Match Rates 
summary(EMatchOut_HS_D_P)
summary(GMatchOut_HS_D_P)
summary(GTMatchOut_HS_D_P)

PMRCompare = (EMatchOut_HS_D_P/GTMatchOut_HS_D_P) - (GMatchOut_HS_D_P/GTMatchOut_HS_D_P)
summary(PMRCompare)
hist(PMRCompare)

summary(EMatchOut_HS_D_P-GTMatchOut_HS_D_P)
summary(GMatchOut_HS_D_P-GTMatchOut_HS_D_P)
summary(EMatchOut_HS_D_P-GMatchOut_HS_D_P)

# Option 1 
mean(EMatchOut_HS_CA_P-GMatchOut_HS_CA_P)*100
mean(EMatchOut_HS_CP_P-GMatchOut_HS_CP_P)*100
mean(EMatchOut_HS_D_P-GMatchOut_HS_D_P)*100

sd(EMatchOut_HS_CA_P-GMatchOut_HS_CA_P)*100
sd(EMatchOut_HS_CP_P-GMatchOut_HS_CP_P)*100
sd(EMatchOut_HS_D_P-GMatchOut_HS_D_P)*100

# All are significant 
t.test(EMatchOut_HS_CA_P,GMatchOut_HS_CA_P,paired=T)
t.test(EMatchOut_HS_CP_P,GMatchOut_HS_CP_P,paired=T)
t.test(EMatchOut_HS_D_P,GMatchOut_HS_D_P,paired=T)

# Option 2 
mean( (EMatchOut_HS_CA_P/GTMatchOut_HS_CA_P) - (GMatchOut_HS_CA_P/GTMatchOut_HS_CA_P) )*100
mean( (EMatchOut_HS_CP_P/GTMatchOut_HS_CP_P) - (GMatchOut_HS_CP_P/GTMatchOut_HS_CP_P) )*100
mean( (EMatchOut_HS_D_P/GTMatchOut_HS_D_P) - (GMatchOut_HS_D_P/GTMatchOut_HS_D_P) )*100

sd( (EMatchOut_HS_CA_P/GTMatchOut_HS_CA_P) - (GMatchOut_HS_CA_P/GTMatchOut_HS_CA_P) )*100
sd( (EMatchOut_HS_CP_P/GTMatchOut_HS_CP_P) - (GMatchOut_HS_CP_P/GTMatchOut_HS_CP_P) )*100
sd( (EMatchOut_HS_D_P/GTMatchOut_HS_D_P) - (GMatchOut_HS_D_P/GTMatchOut_HS_D_P) )*100

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSHF : Y1 Coefficients  SC ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BLSC= EMeanOut_HS_CA_P
Coefs_GSC = GMeanOut_HS_CA_P
Coefs_PBSC= GTMeanOut_HS_CA_P
theta = c(17.1,.65,2.02,-1.20)

# Option 1 
CompareBL = (Coefs_BLSC - Coefs_PBSC)
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = (Coefs_GSC - Coefs_PBSC)
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))
summary(Compare*100)
apply(Compare*100,2,mean)

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

# Option 2
CompareBL = sweep(Coefs_BLSC,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GSC,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(abs(CompareBL[,1]),abs(CompareG[,1]),paired=T) # BL 
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # BL 
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # GM 
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # BL 

t.test(Coefs_BLSC[,1],Coefs_GSC[,1],paired=T) # BL 
t.test(Coefs_BLSC[,2],Coefs_GSC[,2],paired=T) # BL 
t.test(Coefs_BLSC[,3],Coefs_GSC[,3],paired=T) # GM 
t.test(Coefs_BLSC[,4],Coefs_GSC[,4],paired=T) # BL 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSHF : Y1 Coefficients  D ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BLW= EMeanOut_HS_D_P
Coefs_GW = GMeanOut_HS_D_P
Coefs_PBW= GTMeanOut_HS_D_P

# Option 2
CompareBL = sweep(Coefs_BLW,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GW,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(abs(CompareBL[,1]),abs(CompareG[,1]),paired=T) # BL 
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # BL 
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # GM 
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # BL 

t.test(Coefs_BLW[,1],Coefs_GW[,1],paired=T) # 
t.test(Coefs_BLW[,2],Coefs_GW[,2],paired=T) # 
t.test(Coefs_BLW[,3],Coefs_GW[,3],paired=T) # 
t.test(Coefs_BLW[,4],Coefs_GW[,4],paired=T) #  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSHF : Y1 Coefficients  CP ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BLSI= EMeanOut_HS_CP_P
Coefs_GSI = GMeanOut_HS_CP_P
Coefs_PBSI= GTMeanOut_HS_CP_P

# Option 2
CompareBL = sweep(Coefs_BLSI,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GSI,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(CompareBL[,1],CompareG[,1],paired=T) # BL 
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # BL 
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # NOT SIGNIFICANT
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # BL 

t.test(Coefs_BLSI[,1],Coefs_GSI[,1],paired=T) # 
t.test(Coefs_BLSI[,2],Coefs_GSI[,2],paired=T) # 
t.test(Coefs_BLSI[,3],Coefs_GSI[,3],paired=T) # 
t.test(Coefs_BLSI[,4],Coefs_GSI[,4],paired=T) # 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSHF : Y2 Coefficients  CA ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

lm2 = lm(math~ prog + ses + female, data= ml )
eta = lm2$coefficients 

eta = apply(CoefsOut2_HS_CA_P,2,mean)

Coefs_BHSC= EMeanOut2_HS_CA_P[,1:6]
Coefs_GSC = GMeanOut2_HS_CA_P[,1:6]
Coefs_PBSC= GTMeanOut2_HS_CA_P[,1:6]

CompareBL = sweep(Coefs_BHSC,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")
CompareBL = abs(CompareBL)

CompareG = sweep(Coefs_GSC,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")
CompareG = abs(CompareG)

Compare = ((CompareBL) -(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(CompareBL[,2],CompareG[,2],paired=T) # BL
t.test(CompareBL[,3],CompareG[,3],paired=T) # BL
t.test(CompareBL[,4],CompareG[,4],paired=T) # BL 
t.test(CompareBL[,5],CompareG[,5],paired=T) # GM 
t.test(CompareBL[,6],CompareG[,6],paired=T) # BL  

# let's try something different 
CompareBL = apply(Coefs_BHSC,2,mean)
CompareBL = (CompareBL-eta)/eta 

CompareG = apply(Coefs_GSC,2,mean)
CompareG = (CompareG-eta)/eta 

(abs(CompareBL)-abs(CompareG))*100

# The values are very similiar, we jsut get no SDs this way. 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSHF : Y2 Coefficients  D ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

eta = apply(CoefsOut2_HS_D_P,2,mean)

Coefs_BLW= EMeanOut2_HS_D_P[,1:6]
Coefs_GW = GMeanOut2_HS_D_P[,1:6]
Coefs_PBW= GTMeanOut2_HS_D_P[,1:6]

# Option 2
CompareBL = sweep(Coefs_BLW,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")
CompareBL = abs(CompareBL)

CompareG = sweep(Coefs_GW,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")
CompareG = abs(CompareG)

Compare = ((CompareBL) -(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(CompareBL[,2],CompareG[,2],paired=T) # BL
t.test(CompareBL[,3],CompareG[,3],paired=T) # BL
t.test(CompareBL[,4],CompareG[,4],paired=T) # NOT SIGNIFICANT
t.test(CompareBL[,5],CompareG[,5],paired=T) # GM 
t.test(CompareBL[,6],CompareG[,6],paired=T) # BL  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSHF : Y2 Coefficients  CP ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

eta = apply(CoefsOut2_HS_CP_P,2,mean)

Coefs_BHSI= EMeanOut2_HS_CP_P[,1:6]
Coefs_GSI = GMeanOut2_HS_CP_P[,1:6]
Coefs_PBSI= GTMeanOut2_HS_CP_P[,1:6]

# Option 2
CompareBL = sweep(Coefs_BHSI,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")
CompareBL = abs(CompareBL)

CompareG = sweep(Coefs_GSI,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")
CompareG = abs(CompareG)

Compare = ((CompareBL) -(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(CompareBL[,2],CompareG[,2],paired=T) # BL
t.test(CompareBL[,3],CompareG[,3],paired=T) # BL
t.test(CompareBL[,4],CompareG[,4],paired=T) # NOT SIGNIFICANT
t.test(CompareBL[,5],CompareG[,5],paired=T) # GM 
t.test(CompareBL[,6],CompareG[,6],paired=T) # BL   


#%%%%%%%%%%%%%%%%%%%%%%%%%##
####     HSHF : RMSE     ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%##

RMSE_BLSC= sqrt(E_MSE_HS_CA_P)
RMSE_GSC= sqrt(G_MSE_HS_CA_P)
RMSE_PBSC= sqrt(GT_MSE_HS_CA_P)

RMSE_BLW= sqrt(E_MSE_HS_D_P)
RMSE_GW= sqrt(G_MSE_HS_D_P)
RMSE_PBW= sqrt(GT_MSE_HS_D_P)

RMSE_BLSI= sqrt(E_MSE_HS_CP_P)
RMSE_GSI= sqrt(G_MSE_HS_CP_P)
RMSE_PBSI= sqrt(GT_MSE_HS_CP_P)

# Option 1 
mean(RMSE_BLSC-RMSE_GSC)
mean(RMSE_BLW-RMSE_GW)
mean(RMSE_BLSI-RMSE_GSI)

sd(EMatchOut_HS_CA_P-GMatchOut_HS_CA_P)
sd(EMatchOut_HS_CP_P-GMatchOut_HS_CP_P)
sd(EMatchOut_HS_D_P-GMatchOut_HS_D_P)

# All are significant 
t.test(RMSE_BLSC,RMSE_GSC,paired=T)
t.test(RMSE_BLW,RMSE_GW,paired=T)
t.test(RMSE_BLSI,RMSE_GSI,paired=T)

# Option 2 
mean( (RMSE_BLSC/RMSE_PBSC) - (RMSE_GSC/RMSE_PBSC) )*100
mean( (RMSE_BLW/RMSE_PBW) - (RMSE_GW/RMSE_PBW) )*100
mean( (RMSE_BLSI/RMSE_PBSI) - (RMSE_GSI/RMSE_PBSI) )*100

sd( (RMSE_BLSC/RMSE_PBSC) - (RMSE_GSC/RMSE_PBSC) )*100
sd( (RMSE_BLW/RMSE_PBW) - (RMSE_GW/RMSE_PBW) )*100
sd( (RMSE_BLSI/RMSE_PBSI) - (RMSE_GSI/RMSE_PBSI) )*100

##%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
###  Plot: PMR  #### 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

# Plot for paper
pdf("Paper_HEHS_Match.pdf",width = 11, height = 7)
par(mar=c(5,3,2,2)+0.1)
boxplot((GTMatchOut_HS_CA_P-.6)/.4*100,(GTMatchOut_HS_D_P-.6)/.4*100,(GTMatchOut_HS_CP_P-.6)/.4*100,(GMatchOut_HS_CA_P-.6)/.4*100,(GMatchOut_HS_D_P-.6)/.4*100,(GMatchOut_HS_CP_P-.6)/.4*100,(EMatchOut_HS_CA_P-.6)/.4*100, (EMatchOut_HS_D_P-.6)/.4*100,(EMatchOut_HS_CP_P-.6)/.4*100,names=c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"),ylab="NonSeed Posterior Match Rate")
dev.off()

##%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
###  Plot: MSE  #### 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

pdf("MSE.pdf",width = 11, height = 7)
par(mar=c(5,3,2,2)+0.1)
boxplot(sqrt(GT_MSE_HS_CA_P),sqrt(GT_MSE_HS_D_P[-72]),sqrt(GT_MSE_HS_CP_P),sqrt(G_MSE_HS_CA_P), sqrt(G_MSE_HS_D_P),sqrt(G_MSE_HS_CP_P),sqrt(E_MSE_HS_CA_P),sqrt(E_MSE_HS_D_P[-72]),sqrt(E_MSE_HS_CP_P), ylab= "Root MSE", names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"))
dev.off()

pdf("Paper_MSE.pdf",width = 11, height = 7)
par(mar=c(5,3,2,2)+0.1)
boxplot(sqrt(GT_MSE_HS_CA_P),sqrt(GT_MSE_HS_D_P),sqrt(GT_MSE_HS_CP_P),sqrt(G_MSE_HS_CA_P), sqrt(G_MSE_HS_D_P),sqrt(G_MSE_HS_CP_P),sqrt(E_MSE_HS_CA_P),sqrt(E_MSE_HS_D_P),sqrt(E_MSE_HS_CP_P), ylab= "Root MSE", names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"))
dev.off()

##### Parameter Estimates: Plots ########

install.packages("plotrix")
require(plotrix)

#Step 1: The Mean Vector
meanVec = c( mean(GTMeanOut_HS_CA_P[,1]),mean(GTMeanOut_HS_D_P[,1]),mean(GTMeanOut_HS_CP_P[,1]),mean(GMeanOut_HS_CA_P[,1]),mean(GMeanOut_HS_D_P[,1]),mean(GMeanOut_HS_CP_P[,1]),mean(EMeanOut_HS_CA_P[,1]),mean(EMeanOut_HS_D_P[,1]),mean(EMeanOut_HS_CP_P[,1]))

li1 = c( quantile(GTMeanOut_HS_CA_P[,1],c(0.025)),quantile(GTMeanOut_HS_D_P[,1],c(0.025)),quantile(GTMeanOut_HS_CP_P[,1],c(0.025)),quantile(GMeanOut_HS_CA_P[,1],c(0.025)),quantile(GMeanOut_HS_D_P[,1],c(0.025)),quantile(GMeanOut_HS_CP_P[,1],c(0.025)),quantile(EMeanOut_HS_CA_P[,1],c(0.025)),quantile(EMeanOut_HS_D_P[,1],c(0.025)),quantile(EMeanOut_HS_CP_P[,1],c(0.025)))

ui1 = c( quantile(GTMeanOut_HS_CA_P[,1],c(0.975)),quantile(GTMeanOut_HS_D_P[,1],c(0.975)),quantile(GTMeanOut_HS_CP_P[,1],c(0.975)),quantile(GMeanOut_HS_CA_P[,1],c(0.975)),quantile(GMeanOut_HS_D_P[,1],c(0.975)),quantile(GMeanOut_HS_CP_P[,1],c(0.975)),quantile(EMeanOut_HS_CA_P[,1],c(0.975)),quantile(EMeanOut_HS_D_P[,1],c(0.975)),quantile(EMeanOut_HS_CP_P[,1],c(0.975)))

# Step 3: The Plot 
pdf('Paper_HFHS_All_Intercept_CIs.pdf',width = 8, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeHS and axes 
 boxplot(  GTMeanOut_HS_CA_P[,1], GTMeanOut_HS_D_P[,1], GTMeanOut_HS_CP_P[,1],GMeanOut_HS_CA_P[,1], GMeanOut_HS_D_P[,1], GMeanOut_HS_CP_P[,1], EMeanOut_HS_CA_P[,1], EMeanOut_HS_D_P[,1], EMeanOut_HS_CP_P[,1], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x=c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[1],lty = 2)
dev.off()

meanVec = c( mean(GTMeanOut_HS_CA_P[,2]),mean(GTMeanOut_HS_D_P[,2]),mean(GTMeanOut_HS_CP_P[,2]),mean(GMeanOut_HS_CA_P[,2]),mean(GMeanOut_HS_D_P[,2]),mean(GMeanOut_HS_CP_P[,2]),mean(EMeanOut_HS_CA_P[,2]),mean(EMeanOut_HS_D_P[,2]),mean(EMeanOut_HS_CP_P[,2]))

li1 = c( quantile(GTMeanOut_HS_CA_P[,2],c(0.025)),quantile(GTMeanOut_HS_D_P[,2],c(0.025)),quantile(GTMeanOut_HS_CP_P[,2],c(0.025)),quantile(GMeanOut_HS_CA_P[,2],c(0.025)),quantile(GMeanOut_HS_D_P[,2],c(0.025)),quantile(GMeanOut_HS_CP_P[,2],c(0.025)),quantile(EMeanOut_HS_CA_P[,2],c(0.025)),quantile(EMeanOut_HS_D_P[,2],c(0.025)),quantile(EMeanOut_HS_CP_P[,2],c(0.025)))

ui1 = c( quantile(GTMeanOut_HS_CA_P[,2],c(0.975)),quantile(GTMeanOut_HS_D_P[,2],c(0.975)),quantile(GTMeanOut_HS_CP_P[,2],c(0.975)),quantile(GMeanOut_HS_CA_P[,2],c(0.975)),quantile(GMeanOut_HS_D_P[,2],c(0.975)),quantile(GMeanOut_HS_CP_P[,2],c(0.975)),quantile(EMeanOut_HS_CA_P[,2],c(0.975)),quantile(EMeanOut_HS_D_P[,2],c(0.975)),quantile(EMeanOut_HS_CP_P[,2],c(0.975)))

pdf('Paper_HFHS_All_Math_CIs.pdf',width = 8, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeHS and axes 
 boxplot(  GTMeanOut_HS_CA_P[,2], GTMeanOut_HS_D_P[,2], GTMeanOut_HS_CP_P[,2],GMeanOut_HS_CA_P[,2], GMeanOut_HS_D_P[,2], GMeanOut_HS_CP_P[,2], EMeanOut_HS_CA_P[,2], EMeanOut_HS_D_P[,2], EMeanOut_HS_CP_P[,2], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white",ylim = c(0.45,0.7), at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x= c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[2],lty = 2)
dev.off()

meanVec = c( mean(GTMeanOut_HS_CA_P[,3]),mean(GTMeanOut_HS_D_P[,3]),mean(GTMeanOut_HS_CP_P[,3]),mean(GMeanOut_HS_CA_P[,3]),mean(GMeanOut_HS_D_P[,3]),mean(GMeanOut_HS_CP_P[,3]),mean(EMeanOut_HS_CA_P[,3]),mean(EMeanOut_HS_D_P[,3]),mean(EMeanOut_HS_CP_P[,3]))

li1 = c( quantile(GTMeanOut_HS_CA_P[,3],c(0.025)),quantile(GTMeanOut_HS_D_P[,3],c(0.025)),quantile(GTMeanOut_HS_CP_P[,3],c(0.025)),quantile(GMeanOut_HS_CA_P[,3],c(0.025)),quantile(GMeanOut_HS_D_P[,3],c(0.025)),quantile(GMeanOut_HS_CP_P[,3],c(0.025)),quantile(EMeanOut_HS_CA_P[,3],c(0.025)),quantile(EMeanOut_HS_D_P[,3],c(0.025)),quantile(EMeanOut_HS_CP_P[,3],c(0.025)))

ui1 = c( quantile(GTMeanOut_HS_CA_P[,3],c(0.975)),quantile(GTMeanOut_HS_D_P[,3],c(0.975)),quantile(GTMeanOut_HS_CP_P[,3],c(0.975)),quantile(GMeanOut_HS_CA_P[,3],c(0.975)),quantile(GMeanOut_HS_D_P[,3],c(0.975)),quantile(GMeanOut_HS_CP_P[,3],c(0.975)),quantile(EMeanOut_HS_CA_P[,3],c(0.975)),quantile(EMeanOut_HS_D_P[,3],c(0.975)),quantile(EMeanOut_HS_CP_P[,3],c(0.975)))

pdf('Paper_HFHS_All_Academic_CIs.pdf',width = 8, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeHS and axes 
 boxplot(  GTMeanOut_HS_CA_P[,3], GTMeanOut_HS_D_P[,3], GTMeanOut_HS_CP_P[,3],GMeanOut_HS_CA_P[,3], GMeanOut_HS_D_P[,3], GMeanOut_HS_CP_P[,3], EMeanOut_HS_CA_P[,3], EMeanOut_HS_D_P[,3], EMeanOut_HS_CP_P[,3], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x=c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[3],lty = 2)
dev.off()


meanVec = c( mean(GTMeanOut_HS_CA_P[,4]),mean(GTMeanOut_HS_D_P[,4]),mean(GTMeanOut_HS_CP_P[,4]),mean(GMeanOut_HS_CA_P[,4]),mean(GMeanOut_HS_D_P[,4]),mean(GMeanOut_HS_CP_P[,4]),mean(EMeanOut_HS_CA_P[,4]),mean(EMeanOut_HS_D_P[,4]),mean(EMeanOut_HS_CP_P[,4]))

li1 = c( quantile(GTMeanOut_HS_CA_P[,4],c(0.025)),quantile(GTMeanOut_HS_D_P[,4],c(0.025)),quantile(GTMeanOut_HS_CP_P[,4],c(0.025)),quantile(GMeanOut_HS_CA_P[,4],c(0.025)),quantile(GMeanOut_HS_D_P[,4],c(0.025)),quantile(GMeanOut_HS_CP_P[,4],c(0.025)),quantile(EMeanOut_HS_CA_P[,4],c(0.025)),quantile(EMeanOut_HS_D_P[,4],c(0.025)),quantile(EMeanOut_HS_CP_P[,4],c(0.025)))

ui1 = c( quantile(GTMeanOut_HS_CA_P[,4],c(0.975)),quantile(GTMeanOut_HS_D_P[,4],c(0.975)),quantile(GTMeanOut_HS_CP_P[,4],c(0.975)),quantile(GMeanOut_HS_CA_P[,4],c(0.975)),quantile(GMeanOut_HS_D_P[,4],c(0.975)),quantile(GMeanOut_HS_CP_P[,4],c(0.975)),quantile(EMeanOut_HS_CA_P[,4],c(0.975)),quantile(EMeanOut_HS_D_P[,4],c(0.975)),quantile(EMeanOut_HS_CP_P[,4],c(0.975)))

pdf('Paper_HFHS_All_Vocational_CIs.pdf',width = 8, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeHS and axes 
 boxplot(  GTMeanOut_HS_CA_P[,4], GTMeanOut_HS_D_P[,4], GTMeanOut_HS_CP_P[,4],GMeanOut_HS_CA_P[,4], GMeanOut_HS_D_P[,4], GMeanOut_HS_CP_P[,4], EMeanOut_HS_CA_P[,4], EMeanOut_HS_D_P[,4], EMeanOut_HS_CP_P[,4], names = c("PB","PB","PB","GM","GM","GM","B.CA","B.D","B.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x= c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[4],lty = 2)
dev.off()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Clear the Memory 
rm(list = ls() )

require(foreign)
ml <-read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
GenCoefs = lm(read~math+prog,data = ml)$coefficients

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### High Fault, Low Seed, D  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

setwd("BLASERuns/HE_LS_D")

load("OutSpace.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### High Error, Low Seed, SI ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# EDiffOut_LS_CP_P - FORMATTED 

setwd("/home/grad/nmd16/BLASERuns/HE_LS_CP")

load("OutSpace.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### High Error, Low Seed, SC ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# EDiffOut_LS_CA_P

setwd("/home/grad/nmd16/BLASERuns/HE_LS_CA")

load("OutSpace.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%##
####     LSHF : PMR     ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Option 1 
mean(EMatchOut_LS_CA_P-GMatchOut_LS_CA_P)*100
mean(EMatchOut_LS_CP_P-GMatchOut_LS_CP_P)*100
mean(EMatchOut_LS_D_P-GMatchOut_LS_D_P)*100

sd(EMatchOut_LS_CA_P-GMatchOut_LS_CA_P)*100
sd(EMatchOut_LS_CP_P-GMatchOut_LS_CP_P)*100
sd(EMatchOut_LS_D_P-GMatchOut_LS_D_P)*100

# All are significant 
t.test(EMatchOut_LS_CA_P,GMatchOut_LS_CA_P,paired=T)
t.test(EMatchOut_LS_CP_P,GMatchOut_LS_CP_P,paired=T)
t.test(EMatchOut_LS_D_P,GMatchOut_LS_D_P,paired=T)

# Option 2 
mean( (EMatchOut_LS_CA_P/GTMatchOut_LS_CA_P) - (GMatchOut_LS_CA_P/GTMatchOut_LS_CA_P) )*100
mean( (EMatchOut_LS_CP_P/GTMatchOut_LS_CP_P) - (GMatchOut_LS_CP_P/GTMatchOut_LS_CP_P) )*100
mean( (EMatchOut_LS_D_P/GTMatchOut_LS_D_P) - (GMatchOut_LS_D_P/GTMatchOut_LS_D_P) )*100

sd( (EMatchOut_LS_CA_P/GTMatchOut_LS_CA_P) - (GMatchOut_LS_CA_P/GTMatchOut_LS_CA_P) )*100
sd( (EMatchOut_LS_CP_P/GTMatchOut_LS_CP_P) - (GMatchOut_LS_CP_P/GTMatchOut_LS_CP_P) )*100 
sd( (EMatchOut_LS_D_P/GTMatchOut_LS_D_P) - (GMatchOut_LS_D_P/GTMatchOut_LS_D_P) )*100

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSHF : Y1 Coefficients  SC ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BLSC= EMeanOut_LS_CA_P
Coefs_GSC = GMeanOut_LS_CA_P
Coefs_PBSC= GTMeanOut_LS_CA_P
theta = c(17.1,.65,2.02,-1.20)

# Option 2
CompareBL = sweep(Coefs_BLSC,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GSC,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(CompareBL[,1],CompareG[,1],paired=T) # NOT
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # NOT
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # GM 
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # NOT 

t.test(Coefs_BLSC[,1],Coefs_GSC[,1],paired=T) # NOT
t.test(Coefs_BLSC[,2],Coefs_GSC[,2],paired=T) # NOT 
t.test(Coefs_BLSC[,3],Coefs_GSC[,3],paired=T) # GM 
t.test(Coefs_BLSC[,4],Coefs_GSC[,4],paired=T) # NOT

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSHF : Y1 Coefficients  W ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BLW= EMeanOut_LS_D_P
Coefs_GW = GMeanOut_LS_D_P
Coefs_PBW= GTMeanOut_LS_D_P

# Option 2
CompareBL = sweep(Coefs_BLW,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GW,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(CompareBL[,1],CompareG[,1],paired=T) # BL 
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # BL 
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # GM 
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # BL 

t.test(Coefs_BLW[,1],Coefs_GW[,1],paired=T) #  
t.test(Coefs_BLW[,2],Coefs_GW[,2],paired=T) # 
t.test(Coefs_BLW[,3],Coefs_GW[,3],paired=T) # 
t.test(Coefs_BLW[,4],Coefs_GW[,4],paired=T) # 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSHF : Y1 Coefficients  SI ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BLSI= EMeanOut_LS_CP_P
Coefs_GSI = GMeanOut_LS_CP_P
Coefs_PBSI= GTMeanOut_LS_CP_P

# Option 2
CompareBL = sweep(Coefs_BLSI,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GSI,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(CompareBL[,1],CompareG[,1],paired=T) # BL 
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # BL 
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # GM
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # BL 

t.test(Coefs_BLSI[,1],Coefs_GSI[,1],paired=T) #  
t.test(Coefs_BLSI[,2],Coefs_GSI[,2],paired=T) # 
t.test(Coefs_BLSI[,3],Coefs_GSI[,3],paired=T) # 
t.test(Coefs_BLSI[,4],Coefs_GSI[,4],paired=T) # 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSHF : Y2 Coefficients  SC ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

lm2 = lm(math~ prog + ses + female, data= ml )
eta = lm2$coefficients 

eta = apply(Coefs_Out2_LS_CA_P,2,mean)

Coefs_BLSC= EMeanOut2_LS_CA_P[,1:6]
Coefs_GSC = GMeanOut2_LS_CA_P[,1:6]
Coefs_PBSC= GTMeanOut2_LS_CA_P[,1:6]

# Option 2
CompareBL = sweep(Coefs_BLSC,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")

CompareG = sweep(Coefs_GSC,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # GM
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # BL 
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # BL
t.test(abs(CompareBL[,5]),abs(CompareG[,5]),paired=T) # BL
t.test(abs(CompareBL[,6]),abs(CompareG[,6]),paired=T) # BL  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSHF : Y2 Coefficients  W ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

eta = apply(CoefsOut2_LS_D_P,2,mean)

Coefs_BLW= EMeanOut2_LS_D_P[,1:6]
Coefs_GW = GMeanOut2_LS_D_P[,1:6]
Coefs_PBW= GTMeanOut2_LS_D_P[,1:6]

# Option 2
CompareBL = sweep(Coefs_BLW,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")

CompareG = sweep(Coefs_GW,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))
apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # GM
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # BL 
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # BL
t.test(abs(CompareBL[,5]),abs(CompareG[,5]),paired=T) # BL
t.test(abs(CompareBL[,6]),abs(CompareG[,6]),paired=T) # BL  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSHF : Y2 Coefficients  SI ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

eta = apply(CoefsOut2_LS_CP_P,2,mean)

Coefs_BLSI= EMeanOut2_LS_CP_P[,1:6]
Coefs_GSI = GMeanOut2_LS_CP_P[,1:6]
Coefs_PBSI= GTMeanOut2_LS_CP_P[,1:6]

# Option 2
CompareBL = sweep(Coefs_BLSI,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")

CompareG = sweep(Coefs_GSI,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))
apply(Compare*100,2,mean)
apply(Compare*100,2,sd) 

t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # GM
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # BL 
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # BL
t.test(abs(CompareBL[,5]),abs(CompareG[,5]),paired=T) # BL
t.test(abs(CompareBL[,6]),abs(CompareG[,6]),paired=T) # BL  

#%%%%%%%%%%%%%%%%%%%%%%%%%##
####     LSHF : RMSE     ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%##

RMSE_BLSC= sqrt(E_MSE_LS_CA_P)
RMSE_GSC= sqrt(G_MSE_LS_CA_P)
RMSE_PBSC= sqrt(GT_MSE_LS_CA_P)

RMSE_BLW= sqrt(E_MSE_LS_D_P)
RMSE_GW= sqrt(G_MSE_LS_D_P)
RMSE_PBW= sqrt(GT_MSE_LS_D_P)

RMSE_BLSI= sqrt(E_MSE_LS_CP_P)
RMSE_GSI= sqrt(G_MSE_LS_CP_P)
RMSE_PBSI= sqrt(GT_MSE_LS_CP_P)

# All are significant 
t.test(RMSE_BLSC,RMSE_GSC,paired=T)
t.test(RMSE_BLW,RMSE_GW,paired=T)
t.test(RMSE_BLSI,RMSE_GSI,paired=T)

# Option 2 
summary( (RMSE_BLSC/RMSE_PBSC) - (RMSE_GSC/RMSE_PBSC) )
summary( (RMSE_BLW/RMSE_PBW) - (RMSE_GW/RMSE_PBW) )
summary( (RMSE_BLSI/RMSE_PBSI) - (RMSE_GSI/RMSE_PBSI) )

sd( (RMSE_BLSC/RMSE_PBSC) - (RMSE_GSC/RMSE_PBSC) )
sd( (RMSE_BLW/RMSE_PBW) - (RMSE_GW/RMSE_PBW) )
sd( (RMSE_BLSI/RMSE_PBSI) - (RMSE_GSI/RMSE_PBSI) )

##%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
###  Plot: PMR #### 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

# Plot for paper
pdf("Paper_HFLS_Match.pdf",width = 11, height = 7)
par(mar=c(5,3,2,2)+0.1)
boxplot((GTMatchOut_LS_CA_P-.2)/.8*100,(GTMatchOut_LS_D_P-.2)/.8*100,(GTMatchOut_LS_CP_P-.2)/.8*100,(GMatchOut_LS_CA_P-.2)/.8*100,(GMatchOut_LS_D_P-.2)/.8*100,(GMatchOut_LS_CP_P-.2)/.8*100,(EMatchOut_LS_CA_P-.2)/.8*100, (EMatchOut_LS_D_P-.2)/.8*100,(EMatchOut_LS_CP_P-.2)/.8*100,names=c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"),ylab="NonSeed Posterior Match Rate")
dev.off()

##%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
###  Plot: MSE #### 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

pdf("Paper_MSE.pdf",width = 11, height = 7)
boxplot(sqrt(GT_MSE_LS_CA_P),sqrt(GT_MSE_LS_D_P),sqrt(GT_MSE_LS_CA_P),sqrt(G_MSE_LS_CA_P), sqrt(G_MSE_LS_D_P),sqrt(G_MSE_LS_CP_P),sqrt(E_MSE_LS_CA_P),sqrt(E_MSE_LS_D_P),sqrt(E_MSE_LS_CP_P), ylab= "Root MSE", names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"))
dev.off()

###### Plot: Parameter Estimates ########

install.packages("plotrix")
require(plotrix)

meanVec = c( mean(GTMeanOut_LS_CA_P[,1]),mean(GTMeanOut_LS_D_P[,1]),mean(GTMeanOut_LS_CA_P[,1]),mean(GMeanOut_LS_CA_P[,1]),mean(GMeanOut_LS_D_P[,1]),mean(GMeanOut_LS_CP_P[,1]),mean(EMeanOut_LS_CA_P[,1]),mean(EMeanOut_LS_D_P[,1]),mean(EMeanOut_LS_CP_P[,1]))

li1 = c( quantile(GTMeanOut_LS_CA_P[,1],c(0.025)),quantile(GTMeanOut_LS_D_P[,1],c(0.025)),quantile(GTMeanOut_LS_CA_P[,1],c(0.025)),quantile(GMeanOut_LS_CA_P[,1],c(0.025)),quantile(GMeanOut_LS_D_P[,1],c(0.025)),quantile(GMeanOut_LS_CP_P[,1],c(0.025)),quantile(EMeanOut_LS_CA_P[,1],c(0.025)),quantile(EMeanOut_LS_D_P[,1],c(0.025)),quantile(EMeanOut_LS_CP_P[,1],c(0.025)))

ui1 = c( quantile(GTMeanOut_LS_CA_P[,1],c(0.975)),quantile(GTMeanOut_LS_D_P[,1],c(0.975)),quantile(GTMeanOut_LS_CA_P[,1],c(0.975)),quantile(GMeanOut_LS_CA_P[,1],c(0.975)),quantile(GMeanOut_LS_D_P[,1],c(0.975)),quantile(GMeanOut_LS_CP_P[,1],c(0.975)),quantile(EMeanOut_LS_CA_P[,1],c(0.975)),quantile(EMeanOut_LS_D_P[,1],c(0.975)),quantile(EMeanOut_LS_CP_P[,1],c(0.975)))

pdf('Paper_HFLS_All_Intercept_CIs.pdf',width=8,height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labels and axes 
 boxplot(  GTMeanOut_LS_CA_P[,1], GTMeanOut_LS_D_P[,1], GTMeanOut_LS_CP_P[,1],GMeanOut_LS_CA_P[,1], GMeanOut_LS_D_P[,1], GMeanOut_LS_CP_P[,1], EMeanOut_LS_CA_P[,1], EMeanOut_LS_D_P[,1], EMeanOut_LS_CP_P[,1], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x= c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[1],lty = 2)
dev.off()

meanVec = c( mean(GTMeanOut_LS_CA_P[,2]),mean(GTMeanOut_LS_D_P[,2]),mean(GTMeanOut_LS_CP_P[,2]),mean(GMeanOut_LS_CA_P[,2]),mean(GMeanOut_LS_D_P[,2]),mean(GMeanOut_LS_CP_P[,2]),mean(EMeanOut_LS_CA_P[,2]),mean(EMeanOut_LS_D_P[,2]),mean(EMeanOut_LS_CP_P[,2]))

li1 = c( quantile(GTMeanOut_LS_CA_P[,2],c(0.025)),quantile(GTMeanOut_LS_D_P[,2],c(0.025)),quantile(GTMeanOut_LS_CP_P[,2],c(0.025)),quantile(GMeanOut_LS_CA_P[,2],c(0.025)),quantile(GMeanOut_LS_D_P[,2],c(0.025)),quantile(GMeanOut_LS_CP_P[,2],c(0.025)),quantile(EMeanOut_LS_CA_P[,2],c(0.025)),quantile(EMeanOut_LS_D_P[,2],c(0.025)),quantile(EMeanOut_LS_CP_P[,2],c(0.025)))

ui1 = c( quantile(GTMeanOut_LS_CA_P[,2],c(0.975)),quantile(GTMeanOut_LS_D_P[,2],c(0.975)),quantile(GTMeanOut_LS_CP_P[,2],c(0.975)),quantile(GMeanOut_LS_CA_P[,2],c(0.975)),quantile(GMeanOut_LS_D_P[,2],c(0.975)),quantile(GMeanOut_LS_CP_P[,2],c(0.975)),quantile(EMeanOut_LS_CA_P[,2],c(0.975)),quantile(EMeanOut_LS_D_P[,2],c(0.975)),quantile(EMeanOut_LS_CP_P[,2],c(0.975)))

pdf('Paper_HFLS_All_Math_CIs.pdf',width=8,height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labels and axes 
 boxplot(  GTMeanOut_LS_CA_P[,2], GTMeanOut_LS_D_P[,2], GTMeanOut_LS_CP_P[,2],GMeanOut_LS_CA_P[,2], GMeanOut_LS_D_P[,2], GMeanOut_LS_CP_P[,2], EMeanOut_LS_CA_P[,2], EMeanOut_LS_D_P[,2], EMeanOut_LS_CP_P[,2], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white",ylim = c(0.2,0.7), at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x= c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[2],lty = 2)
dev.off()

meanVec = c( mean(GTMeanOut_LS_CA_P[,3]),mean(GTMeanOut_LS_D_P[,3]),mean(GTMeanOut_LS_CP_P[,3]),mean(GMeanOut_LS_CA_P[,3]),mean(GMeanOut_LS_D_P[,3]),mean(GMeanOut_LS_CP_P[,3]),mean(EMeanOut_LS_CA_P[,3]),mean(EMeanOut_LS_D_P[,3]),mean(EMeanOut_LS_CP_P[,3]))

li1 = c( quantile(GTMeanOut_LS_CA_P[,3],c(0.025)),quantile(GTMeanOut_LS_D_P[,3],c(0.025)),quantile(GTMeanOut_LS_CP_P[,3],c(0.025)),quantile(GMeanOut_LS_CA_P[,3],c(0.025)),quantile(GMeanOut_LS_D_P[,3],c(0.025)),quantile(GMeanOut_LS_CP_P[,3],c(0.025)),quantile(EMeanOut_LS_CA_P[,3],c(0.025)),quantile(EMeanOut_LS_D_P[,3],c(0.025)),quantile(EMeanOut_LS_CP_P[,3],c(0.025)))

ui1 = c( quantile(GTMeanOut_LS_CA_P[,3],c(0.975)),quantile(GTMeanOut_LS_D_P[,3],c(0.975)),quantile(GTMeanOut_LS_CP_P[,3],c(0.975)),quantile(GMeanOut_LS_CA_P[,3],c(0.975)),quantile(GMeanOut_LS_D_P[,3],c(0.975)),quantile(GMeanOut_LS_CP_P[,3],c(0.975)),quantile(EMeanOut_LS_CA_P[,3],c(0.975)),quantile(EMeanOut_LS_D_P[,3],c(0.975)),quantile(EMeanOut_LS_CP_P[,3],c(0.975)))

pdf('Paper_HFLS_All_Academic_CIs.pdf',width=8,height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labels and axes 
 boxplot(  GTMeanOut_LS_CA_P[,3], GTMeanOut_LS_D_P[,3], GTMeanOut_LS_CP_P[,3],GMeanOut_LS_CA_P[,3], GMeanOut_LS_D_P[,3], GMeanOut_LS_CP_P[,3], EMeanOut_LS_CA_P[,3], EMeanOut_LS_D_P[,3], EMeanOut_LS_CP_P[,3], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x= c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[3],lty = 2)
dev.off()


meanVec = c( mean(GTMeanOut_LS_CA_P[,4]),mean(GTMeanOut_LS_D_P[,4]),mean(GTMeanOut_LS_CP_P[,4]),mean(GMeanOut_LS_CA_P[,4]),mean(GMeanOut_LS_D_P[,4]),mean(GMeanOut_LS_CP_P[,4]),mean(EMeanOut_LS_CA_P[,4]),mean(EMeanOut_LS_D_P[,4]),mean(EMeanOut_LS_CP_P[,4]))

li1 = c( quantile(GTMeanOut_LS_CA_P[,4],c(0.025)),quantile(GTMeanOut_LS_D_P[,4],c(0.025)),quantile(GTMeanOut_LS_CP_P[,4],c(0.025)),quantile(GMeanOut_LS_CA_P[,4],c(0.025)),quantile(GMeanOut_LS_D_P[,4],c(0.025)),quantile(GMeanOut_LS_CP_P[,4],c(0.025)),quantile(EMeanOut_LS_CA_P[,4],c(0.025)),quantile(EMeanOut_LS_D_P[,4],c(0.025)),quantile(EMeanOut_LS_CP_P[,4],c(0.025)))

ui1 = c( quantile(GTMeanOut_LS_CA_P[,4],c(0.975)),quantile(GTMeanOut_LS_D_P[,4],c(0.975)),quantile(GTMeanOut_LS_CP_P[,4],c(0.975)),quantile(GMeanOut_LS_CA_P[,4],c(0.975)),quantile(GMeanOut_LS_D_P[,4],c(0.975)),quantile(GMeanOut_LS_CP_P[,4],c(0.975)),quantile(EMeanOut_LS_CA_P[,4],c(0.975)),quantile(EMeanOut_LS_D_P[,4],c(0.975)),quantile(EMeanOut_LS_CP_P[,4],c(0.975)))

pdf('HFLS/Paper_HFLS_All_Vocational_CIs.pdf',width=8,height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labels and axes 
 boxplot(  GTMeanOut_LS_CA_P[,4], GTMeanOut_LS_D_P[,4], GTMeanOut_LS_CP_P[,4],GMeanOut_LS_CA_P[,4], GMeanOut_LS_D_P[,4], GMeanOut_LS_CP_P[,4], EMeanOut_LS_CA_P[,4], EMeanOut_LS_D_P[,4], EMeanOut_LS_CP_P[,4], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x=c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[4],lty = 2)
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#####         Low Error            #######
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

rm(list=ls())

require(foreign)
ml <-read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
GenCoefs = lm(read~math+prog,data = ml)$coefficients

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### Low Fault, High Seed, CP  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

setwd("/home/grad/nmd16/BLASERuns/LF_HS_CP")
num.sim = 100 
newburn = 1000 
indices  = 1:100
lastcoef = 5 

load("OutSpace.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
####  Low Error, High Seed, SC ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

setwd("/home/grad/nmd16/BLASERuns/LF_HS_CA")
load("OutSpace2.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### Low Error, High Seed, W ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

setwd("/home/grad/nmd16/BLASERuns/LF_HS_D")

load("OutSpace.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%##
####     HSLF : PMR     ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Option 1 
mean(EMatchOut_LF_HS_CA_P-GMatchOut_LF_HS_CA_P)*100
mean(EMatchOut_LF_HS_D_P-GMatchOut_LF_HS_D_P)*100
mean(EMatchOut_LF_HS_CP_P-GMatchOut_LF_HS_CP_P)*100

sd(EMatchOut_LF_HS_CA_P-GMatchOut_LF_HS_CA_P)*100
sd(EMatchOut_LF_HS_D_P-GMatchOut_LF_HS_D_P)*100
sd(EMatchOut_LF_HS_CP_P-GMatchOut_LF_HS_CP_P)*100

# All are significant 
t.test(EMatchOut_LF_HS_CA_P,GMatchOut_LF_HS_CA_P,paired=T)
t.test(EMatchOut_LF_HS_CP_P,GMatchOut_LF_HS_CP_P,paired=T)
t.test(EMatchOut_LF_HS_D_P,GMatchOut_LF_HS_D_P,paired=T)

# Option 2 
mean( (EMatchOut_LF_HS_CA_P/GTMatchOut_LF_HS_CA_P) - (GMatchOut_LF_HS_CA_P/GTMatchOut_LF_HS_CA_P) )*100
mean( (EMatchOut_LF_HS_CP_P/GTMatchOut_LF_HS_CP_P) - (GMatchOut_LF_HS_CP_P/GTMatchOut_LF_HS_CP_P) )*100
mean( (EMatchOut_LF_HS_D_P/GTMatchOut_LF_HS_D_P) - (GMatchOut_LF_HS_D_P/GTMatchOut_LF_HS_D_P) )*100

sd( (EMatchOut_LF_HS_CA_P/GTMatchOut_LF_HS_CA_P) - (GMatchOut_LF_HS_CA_P/GTMatchOut_LF_HS_CA_P) )*100
sd( (EMatchOut_LF_HS_CP_P/GTMatchOut_LF_HS_CP_P) - (GMatchOut_LF_HS_CP_P/GTMatchOut_LF_HS_CP_P) )*100
sd( (EMatchOut_LF_HS_D_P/GTMatchOut_LF_HS_D_P) - (GMatchOut_LF_HS_D_P/GTMatchOut_LF_HS_D_P) )*100

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSLF : Y1 Coefficients  SC ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BLSC= EMeanOut_LF_HS_CA_P
Coefs_GSC = GMeanOut_LF_HS_CA_P
Coefs_PBSC= GTMeanOut_LF_HS_CA_P
theta = c(17.1,.65,2.02,-1.20)

# Option 1 
CompareBL = (Coefs_BLSC - Coefs_PBSC)
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = (Coefs_GSC - Coefs_PBSC)
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

# Option 2
CompareBL = sweep(Coefs_BLSC,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GSC,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)
 

t.test(CompareBL[,1],CompareG[,1],paired=T) # BL
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # BL
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # BL
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # BL

t.test(Coefs_BLSC[,1],Coefs_GSC[,1],paired=T) # 
t.test(Coefs_BLSC[,2],Coefs_GSC[,2],paired=T) #  
t.test(Coefs_BLSC[,3],Coefs_GSC[,3],paired=T) # 
t.test(Coefs_BLSC[,4],Coefs_GSC[,4],paired=T) # 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSLF : Y1 Coefficients  W ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BLW= EMeanOut_LF_HS_D_P
Coefs_GW = GMeanOut_LF_HS_D_P
Coefs_PBW= GTMeanOut_LF_HS_D_P

# Option 1 
CompareBL = (Coefs_BLW - Coefs_PBW)
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = (Coefs_GW - Coefs_PBW)
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))
apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

# Option 2
CompareBL = sweep(Coefs_BLW,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GW,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(CompareBL[,1],CompareG[,1],paired=T) # NOT 
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # BL (barely)
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # BL 
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # BL 


t.test(Coefs_BLW[,1],Coefs_GW[,1],paired=T) # NOT 
t.test(Coefs_BLW[,2],Coefs_GW[,2],paired=T) # barely 
t.test(Coefs_BLW[,3],Coefs_GW[,3],paired=T) #  
t.test(Coefs_BLW[,4],Coefs_GW[,4],paired=T) # 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSLF : Y1 Coefficients  SI ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BHSI= EMeanOut_LF_HS_CP_P
Coefs_GSI = GMeanOut_LF_HS_CP_P
Coefs_PBSI= GTMeanOut_LF_HS_CP_P

# Option 2
CompareBL = sweep(Coefs_BHSI,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GSI,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)
# This is a better metric. It doesn't compare relative to PB, but it DOES provide a better indication of accuracy. 

t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # GM 
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # GM
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # GM 

t.test(Coefs_BHSI[,1],Coefs_GSI[,1],paired=T) #  
t.test(Coefs_BHSI[,2],Coefs_GSI[,2],paired=T) # 
t.test(Coefs_BHSI[,3],Coefs_GSI[,3],paired=T) # 
t.test(Coefs_BHSI[,4],Coefs_GSI[,4],paired=T) # 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSLF : Y2 Coefficients  SC ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

lm2 = lm(math~ prog + ses + female, data= ml )
eta = lm2$coefficients 

eta = apply(CoefsOut2_LF_HS_CA_P,2,mean)

Coefs_BHSC= EMeanOut2_LF_HS_CA_P[,1:6]
Coefs_GSC = GMeanOut2_LF_HS_CA_P[,1:6]
Coefs_PBSC= GTMeanOut2_LF_HS_CA_P[,1:6]

CompareBL = sweep(Coefs_BHSC,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")
CompareBL = abs(CompareBL)

CompareG = sweep(Coefs_GSC,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")
CompareG = abs(CompareG)

Compare = ((CompareBL) -(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(CompareBL[,2],CompareG[,2],paired=T) # BL
t.test(CompareBL[,3],CompareG[,3],paired=T) # GM 
t.test(CompareBL[,4],CompareG[,4],paired=T) # NOT SIGNIFICANT 
t.test(CompareBL[,5],CompareG[,5],paired=T) # NOT SIGNIFICANT 
t.test(CompareBL[,6],CompareG[,6],paired=T) # NOT SIGNIFICANT   

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSLF : Y2 Coefficients  W ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

eta = apply(CoefsOut2_LF_HS_D_P,2,mean)

Coefs_BLW= EMeanOut2_LF_HS_D_P[,1:6]
Coefs_GW = GMeanOut2_LF_HS_D_P[,1:6]
Coefs_PBW= GTMeanOut2_LF_HS_D_P[,1:6]

# Option 2
CompareBL = sweep(Coefs_BLW,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")
CompareBL = abs(CompareBL)

CompareG = sweep(Coefs_GW,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")
CompareG = abs(CompareG)

Compare = ((CompareBL) -(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(CompareBL[,2],CompareG[,2],paired=T) # BL
t.test(CompareBL[,3],CompareG[,3],paired=T) # GM 
t.test(CompareBL[,4],CompareG[,4],paired=T) # NOT SIGNIFICANT 
t.test(CompareBL[,5],CompareG[,5],paired=T) # NOT SIGNIFICANT 
t.test(CompareBL[,6],CompareG[,6],paired=T) # NOT SIGNIFICANT   


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### HSLF : Y2 Coefficients  SI ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

eta = apply(CoefsOut2_LF_HS_CP_P,2,mean)

Coefs_BHSI= EMeanOut2_LF_HS_CP_P[,1:6]
Coefs_GSI = GMeanOut2_LF_HS_CP_P[,1:6]
Coefs_PBSI= GTMeanOut2_LF_HS_CP_P[,1:6]

# Option 2
CompareBL = sweep(Coefs_BHSI,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")
CompareBL = abs(CompareBL)

CompareG = sweep(Coefs_GSI,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")
CompareG = abs(CompareG)

Compare = ((CompareBL) -(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(CompareBL[,2],CompareG[,2],paired=T) # GM
t.test(CompareBL[,3],CompareG[,3],paired=T) # GM 
t.test(CompareBL[,4],CompareG[,4],paired=T) # GM
t.test(CompareBL[,5],CompareG[,5],paired=T) # BL 
t.test(CompareBL[,6],CompareG[,6],paired=T) # GM    


#%%%%%%%%%%%%%%%%%%%%%%%%%##
####     HSLF : RMSE     ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%##

RMSE_BHSC= sqrt(E_MSE_LF_HS_CA_P)
RMSE_GSC= sqrt(G_MSE_LF_HS_CA_P)
RMSE_PBSC= sqrt(GT_MSE_LF_HS_CA_P)

RMSE_BLW= sqrt(E_MSE_LF_HS_D_P)
RMSE_GW= sqrt(G_MSE_LF_HS_D_P)
RMSE_PBW= sqrt(GT_MSE_LF_HS_D_P)

RMSE_BHSI= sqrt(E_MSE_LF_HS_CP_P)
RMSE_GSI= sqrt(G_MSE_LF_HS_CP_P)
RMSE_PBSI= sqrt(GT_MSE_LF_HS_CP_P)

t.test(RMSE_BHSC,RMSE_GSC,paired=T) # not  
t.test(RMSE_BLW,RMSE_GW,paired=T)
t.test(RMSE_BHSI,RMSE_GSI,paired=T)

# Option 2 
mean( (RMSE_BHSC/RMSE_PBSC) - (RMSE_GSC/RMSE_PBSC) )*100
mean( (RMSE_BLW/RMSE_PBW) - (RMSE_GW/RMSE_PBW) )*100
mean( (RMSE_BHSI/RMSE_PBSI) - (RMSE_GSI/RMSE_PBSI) )*100

sd( (RMSE_BHSC/RMSE_PBSC) - (RMSE_GSC/RMSE_PBSC) )*100
sd( (RMSE_BLW/RMSE_PBW) - (RMSE_GW/RMSE_PBW) )*100
sd( (RMSE_BHSI/RMSE_PBSI) - (RMSE_GSI/RMSE_PBSI) )*100

##%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
###  Low Error High Seed: Plots #### 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

##### Plots: Match Rates #### 

# Plot for paper
pdf("/home/grad/nmd16/Dropbox/Images/LEHS/Paper_LEHS_Match.pdf",width = 11, height = 7)
par(mar=c(5,3,2,2)+0.1)
boxplot((GTMatchOut_LF_HS_CA_N-.6)/.4*100,(GTMatchOut_LF_HS_D_N-.6)/.4*100,(GTMatchOut_LF_HS_CP_N-.6)/.4*100,(GMatchOut_LF_HS_CA_P-.6)/.4*100,(GMatchOut_LF_HS_D_P-.6)/.4*100,(GMatchOut_LF_HS_CP_P-.6)/.4*100,(EMatchOut_LF_HS_CA_P-.6)/.4*100, (EMatchOut_LF_HS_D_P-.6)/.4*100,(EMatchOut_LF_HS_CP_P-.6)/.4*100,names=c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"),ylab="NonSeed Posterior Match Rate")
dev.off()

##### Plots: MSE #####

pdf("/home/grad/nmd16/Dropbox/Images/LEHS/Paper_MSE.pdf",width = 11, height = 7)
par(mar=c(5,3,2,2)+0.1)
boxplot(sqrt(GT_MSE_LF_HS_CA_N),sqrt(GT_MSE_LF_HS_D_N),sqrt(GT_MSE_LF_HS_CP_N),sqrt(G_MSE_LF_HS_CA_P), sqrt(G_MSE_LF_HS_D_P),sqrt(G_MSE_LF_HS_CP_P),sqrt(E_MSE_LF_HS_CA_P),sqrt(E_MSE_LF_HS_D_P),sqrt(E_MSE_LF_HS_CP_P), ylab= "Root MSE", names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"))
dev.off()

###### Plots: Parameter Estimates ########

install.packages("plotrix")
require(plotrix)

meanVec = c( mean(GTMeanOut_LF_HS_CA_N[,1]),mean(GTMeanOut_LF_HS_D_N[,1]),mean(GTMeanOut_LF_HS_CP_N[,1]),mean(GMeanOut_LF_HS_CA_P[,1]),mean(GMeanOut_LF_HS_D_P[,1]),mean(GMeanOut_LF_HS_CP_P[,1]),mean(EMeanOut_LF_HS_CA_P[,1]),mean(EMeanOut_LF_HS_D_P[,1]),mean(EMeanOut_LF_HS_CP_P[,1]))

li1 = c( quantile(GTMeanOut_LF_HS_CA_N[,1],c(0.025)),quantile(GTMeanOut_LF_HS_D_N[,1],c(0.025)),quantile(GTMeanOut_LF_HS_CP_N[,1],c(0.025)),quantile(GMeanOut_LF_HS_CA_P[,1],c(0.025)),quantile(GMeanOut_LF_HS_D_P[,1],c(0.025)),quantile(GMeanOut_LF_HS_CP_P[,1],c(0.025)),quantile(EMeanOut_LF_HS_CA_P[,1],c(0.025)),quantile(EMeanOut_LF_HS_D_P[,1],c(0.025)),quantile(EMeanOut_LF_HS_CP_P[,1],c(0.025)))

ui1 = c( quantile(GTMeanOut_LF_HS_CA_N[,1],c(0.975)),quantile(GTMeanOut_LF_HS_D_N[,1],c(0.975)),quantile(GTMeanOut_LF_HS_CP_N[,1],c(0.975)),quantile(GMeanOut_LF_HS_CA_P[,1],c(0.975)),quantile(GMeanOut_LF_HS_D_P[,1],c(0.975)),quantile(GMeanOut_LF_HS_CP_P[,1],c(0.975)),quantile(EMeanOut_LF_HS_CA_P[,1],c(0.975)),quantile(EMeanOut_LF_HS_D_P[,1],c(0.975)),quantile(EMeanOut_LF_HS_CP_P[,1],c(0.975)))

pdf('/home/grad/nmd16/Dropbox/Images/LEHS/Paper_LF_HS_All_Intercept_CIs.pdf',width = 8, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeHS and axes 
 boxplot(  GTMeanOut_LF_HS_CA_N[,1], GTMeanOut_LF_HS_D_N[,1], GTMeanOut_LF_HS_CP_N[,1],GMeanOut_LF_HS_CA_P[,1], GMeanOut_LF_HS_D_P[,1], GMeanOut_LF_HS_CP_P[,1], EMeanOut_LF_HS_CA_P[,1], EMeanOut_LF_HS_D_P[,1], EMeanOut_LF_HS_CP_P[,1], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x= c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[1],lty = 2)
dev.off()

meanVec = c( mean(GTMeanOut_LF_HS_CA_N[,2]),mean(GTMeanOut_LF_HS_D_N[,2]),mean(GTMeanOut_LF_HS_CP_N[,2]),mean(GMeanOut_LF_HS_CA_P[,2]),mean(GMeanOut_LF_HS_D_P[,2]),mean(GMeanOut_LF_HS_CP_P[,2]),mean(EMeanOut_LF_HS_CA_P[,2]),mean(EMeanOut_LF_HS_D_P[,2]),mean(EMeanOut_LF_HS_CP_P[,2]))

li1 = c( quantile(GTMeanOut_LF_HS_CA_N[,2],c(0.025)),quantile(GTMeanOut_LF_HS_D_N[,2],c(0.025)),quantile(GTMeanOut_LF_HS_CP_N[,2],c(0.025)),quantile(GMeanOut_LF_HS_CA_P[,2],c(0.025)),quantile(GMeanOut_LF_HS_D_P[,2],c(0.025)),quantile(GMeanOut_LF_HS_CP_P[,2],c(0.025)),quantile(EMeanOut_LF_HS_CA_P[,2],c(0.025)),quantile(EMeanOut_LF_HS_D_P[,2],c(0.025)),quantile(EMeanOut_LF_HS_CP_P[,2],c(0.025)))

ui1 = c( quantile(GTMeanOut_LF_HS_CA_N[,2],c(0.975)),quantile(GTMeanOut_LF_HS_D_N[,2],c(0.975)),quantile(GTMeanOut_LF_HS_CP_N[,2],c(0.975)),quantile(GMeanOut_LF_HS_CA_P[,2],c(0.975)),quantile(GMeanOut_LF_HS_D_P[,2],c(0.975)),quantile(GMeanOut_LF_HS_CP_P[,2],c(0.975)),quantile(EMeanOut_LF_HS_CA_P[,2],c(0.975)),quantile(EMeanOut_LF_HS_D_P[,2],c(0.975)),quantile(EMeanOut_LF_HS_CP_P[,2],c(0.975)))

pdf('/home/grad/nmd16/Dropbox/Images/LEHS/Paper_LF_HS_All_Math_CIs.pdf',width = 8, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeHS and axes 
 boxplot(  GTMeanOut_LF_HS_CA_N[,2], GTMeanOut_LF_HS_D_P[,2], GTMeanOut_LF_HS_CP_P[,2],GMeanOut_LF_HS_CA_P[,2], GMeanOut_LF_HS_D_P[,2], GMeanOut_LF_HS_CP_P[,2], EMeanOut_LF_HS_CA_P[,2], EMeanOut_LF_HS_D_P[,2], EMeanOut_LF_HS_CP_P[,2], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white",ylim = c(0.5,0.7), at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x=c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[2],lty = 2)
dev.off()

meanVec = c( mean(GTMeanOut_LF_HS_CA_N[,3]),mean(GTMeanOut_LF_HS_D_N[,3]),mean(GTMeanOut_LF_HS_CP_N[,3]),mean(GMeanOut_LF_HS_CA_P[,3]),mean(GMeanOut_LF_HS_D_P[,3]),mean(GMeanOut_LF_HS_CP_P[,3]),mean(EMeanOut_LF_HS_CA_P[,3]),mean(EMeanOut_LF_HS_D_P[,3]),mean(EMeanOut_LF_HS_CP_P[,3]))

li1 = c( quantile(GTMeanOut_LF_HS_CA_N[,3],c(0.025)),quantile(GTMeanOut_LF_HS_D_N[,3],c(0.025)),quantile(GTMeanOut_LF_HS_CP_N[,3],c(0.025)),quantile(GMeanOut_LF_HS_CA_P[,3],c(0.025)),quantile(GMeanOut_LF_HS_D_P[,3],c(0.025)),quantile(GMeanOut_LF_HS_CP_P[,3],c(0.025)),quantile(EMeanOut_LF_HS_CA_P[,3],c(0.025)),quantile(EMeanOut_LF_HS_D_P[,3],c(0.025)),quantile(EMeanOut_LF_HS_CP_P[,3],c(0.025)))

ui1 = c( quantile(GTMeanOut_LF_HS_CA_N[,3],c(0.975)),quantile(GTMeanOut_LF_HS_D_N[,3],c(0.975)),quantile(GTMeanOut_LF_HS_CP_N[,3],c(0.975)),quantile(GMeanOut_LF_HS_CA_P[,3],c(0.975)),quantile(GMeanOut_LF_HS_D_P[,3],c(0.975)),quantile(GMeanOut_LF_HS_CP_P[,3],c(0.975)),quantile(EMeanOut_LF_HS_CA_P[,3],c(0.975)),quantile(EMeanOut_LF_HS_D_P[,3],c(0.975)),quantile(EMeanOut_LF_HS_CP_P[,3],c(0.975)))

pdf('/home/grad/nmd16/Dropbox/Images/LEHS/Paper_LF_HS_All_Academic_CIs.pdf',width = 8, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeHS and axes 
 boxplot(  GTMeanOut_LF_HS_CA_N[,3], GTMeanOut_LF_HS_D_N[,3], GTMeanOut_LF_HS_CP_N[,3],GMeanOut_LF_HS_CA_P[,3], GMeanOut_LF_HS_D_P[,3], GMeanOut_LF_HS_CP_P[,3], EMeanOut_LF_HS_CA_P[,3], EMeanOut_LF_HS_D_P[,3], EMeanOut_LF_HS_CP_P[,3], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x= c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[3],lty = 2)
dev.off()


meanVec = c( mean(GTMeanOut_LF_HS_CA_N[,4]),mean(GTMeanOut_LF_HS_D_N[,4]),mean(GTMeanOut_LF_HS_CP_N[,4]),mean(GMeanOut_LF_HS_CA_P[,4]),mean(GMeanOut_LF_HS_D_P[,4]),mean(GMeanOut_LF_HS_CP_P[,4]),mean(EMeanOut_LF_HS_CA_P[,4]),mean(EMeanOut_LF_HS_D_P[,4]),mean(EMeanOut_LF_HS_CP_P[,4]))

li1 = c( quantile(GTMeanOut_LF_HS_CA_N[,4],c(0.025)),quantile(GTMeanOut_LF_HS_D_N[,4],c(0.025)),quantile(GTMeanOut_LF_HS_CP_N[,4],c(0.025)),quantile(GMeanOut_LF_HS_CA_P[,4],c(0.025)),quantile(GMeanOut_LF_HS_D_P[,4],c(0.025)),quantile(GMeanOut_LF_HS_CP_P[,4],c(0.025)),quantile(EMeanOut_LF_HS_CA_P[,4],c(0.025)),quantile(EMeanOut_LF_HS_D_P[,4],c(0.025)),quantile(EMeanOut_LF_HS_CP_P[,4],c(0.025)))

ui1 = c( quantile(GTMeanOut_LF_HS_CA_N[,4],c(0.975)),quantile(GTMeanOut_LF_HS_D_N[,4],c(0.975)),quantile(GTMeanOut_LF_HS_CP_N[,4],c(0.975)),quantile(GMeanOut_LF_HS_CA_P[,4],c(0.975)),quantile(GMeanOut_LF_HS_D_P[,4],c(0.975)),quantile(GMeanOut_LF_HS_CP_P[,4],c(0.975)),quantile(EMeanOut_LF_HS_CA_P[,4],c(0.975)),quantile(EMeanOut_LF_HS_D_P[,4],c(0.975)),quantile(EMeanOut_LF_HS_CP_P[,4],c(0.975)))

pdf('/home/grad/nmd16/Dropbox/Images/LEHS/Paper_LF_HS_All_Vocational_CIs.pdf',width = 8, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeHS and axes 
 boxplot(  GTMeanOut_LF_HS_CA_N[,4], GTMeanOut_LF_HS_D_N[,4], GTMeanOut_LF_HS_CP_N[,4],GMeanOut_LF_HS_CA_P[,4], GMeanOut_LF_HS_D_P[,4], GMeanOut_LF_HS_CP_P[,4], EMeanOut_LF_HS_CA_P[,4], EMeanOut_LF_HS_D_P[,4], EMeanOut_LF_HS_CP_P[,4], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x= c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[4],lty = 2)
dev.off()


### %%%%%%%%%%%%%%%%%%%%%%%%%%%

# Clear the memory 

rm(list=ls())

require(foreign)
ml <-read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
GenCoefs = lm(read~math+prog,data = ml)$coefficients

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### Low Error, Low Seed, SI ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

setwd("/home/grad/nmd16/BLASERuns/LE_LS_SI")
load("OutSpace2.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
####  Low Error, Low Seed, SC ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

setwd("/home/grad/nmd16/BLASERuns/LE_LS_SC")

load("OutSpace2.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### Low Error, Low Seed, W  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

setwd("/home/grad/nmd16/BLASERuns/LE_LS_W")

load("OutSpace2.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%##
####     LSLF : PMR     ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Option 1 
mean(EMatchOut_LF_LS_CA_P-GMatchOut_LF_LS_CA_P)*100
mean(EMatchOut_LF_LS_D_P-GMatchOut_LF_LS_D_P)*100
mean(EMatchOut_LF_LS_CP_P-GMatchOut_LF_LS_CP_P)*100

sd(EMatchOut_LF_LS_CA_P-GMatchOut_LF_LS_CA_P)*100
sd(EMatchOut_LF_LS_D_P-GMatchOut_LF_LS_D_P)*100
sd(EMatchOut_LF_LS_CP_P-GMatchOut_LF_LS_CP_P)*100

# All are significant 
t.test(EMatchOut_LF_LS_CA_P,GMatchOut_LF_LS_CA_P,paired=T)
t.test(EMatchOut_LF_LS_CP_P,GMatchOut_LF_LS_CP_P,paired=T)
t.test(EMatchOut_LF_LS_D_P,GMatchOut_LF_LS_D_P,paired=T)

# Option 2 
mean( (EMatchOut_LF_LS_CA_P/GTMatchOut_LF_LS_CA_P) - (GMatchOut_LF_LS_CA_P/GTMatchOut_LF_LS_CA_P) )*100
mean( (EMatchOut_LF_LS_CP_P/GTMatchOut_LF_LS_CP_P) - (GMatchOut_LF_LS_CP_P/GTMatchOut_LF_LS_CP_P) )*100
mean( (EMatchOut_LF_LS_D_P/GTMatchOut_LF_LS_D_P) - (GMatchOut_LF_LS_D_P/GTMatchOut_LF_LS_D_P) )*100

sd( (EMatchOut_LF_LS_CA_P/GTMatchOut_LF_LS_CA_P) - (GMatchOut_LF_LS_CA_P/GTMatchOut_LF_LS_CA_P) )*100
sd( (EMatchOut_LF_LS_CP_P/GTMatchOut_LF_LS_CP_P) - (GMatchOut_LF_LS_CP_P/GTMatchOut_LF_LS_CP_P) )*100
sd( (EMatchOut_LF_LS_D_P/GTMatchOut_LF_LS_D_P) - (GMatchOut_LF_LS_D_P/GTMatchOut_LF_LS_D_P) )*100 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSLF : Y1 Coefficients  SC ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BLSC= EMeanOut_LF_LS_CA_P
Coefs_GSC = GMeanOut_LF_LS_CA_P
Coefs_PBSC= GTMeanOut_LF_LS_CA_P
theta = c(17.1,.65,2.02,-1.20)

# Option 1 
CompareBL = (Coefs_BLSC - Coefs_PBSC)
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = (Coefs_GSC - Coefs_PBSC)
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

# Option 2
CompareBL = sweep(Coefs_BLSC,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GSC,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)
 

t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # GM
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # NOT
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # GM 

t.test(Coefs_BLSC[,1],Coefs_GSC[,1],paired=T) # 
t.test(Coefs_BLSC[,2],Coefs_GSC[,2],paired=T) #  
t.test(Coefs_BLSC[,3],Coefs_GSC[,3],paired=T) # NOT 
t.test(Coefs_BLSC[,4],Coefs_GSC[,4],paired=T) # 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSLF : Y1 Coefficients  W ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BLW= EMeanOut_LF_LS_D_P
Coefs_GW = GMeanOut_LF_LS_D_P
Coefs_PBW= GTMeanOut_LF_LS_D_P

# Option 1 
CompareBL = (Coefs_BLW - Coefs_PBW)
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = (Coefs_GW - Coefs_PBW)
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))
apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

# Option 2
CompareBL = sweep(Coefs_BLW,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GW,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # GM
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # GM
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # GM 


t.test(Coefs_BLW[,1],Coefs_GW[,1],paired=T) #  
t.test(Coefs_BLW[,2],Coefs_GW[,2],paired=T) # 
t.test(Coefs_BLW[,3],Coefs_GW[,3],paired=T) # 
t.test(Coefs_BLW[,4],Coefs_GW[,4],paired=T) # 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSLF : Y1 Coefficients  SI ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

Coefs_BLSI= EMeanOut_LF_LS_CP_P
Coefs_GSI = GMeanOut_LF_LS_CP_P
Coefs_PBSI= GTMeanOut_LF_LS_CP_P

# Option 2
CompareBL = sweep(Coefs_BLSI,2,theta,"-")
CompareBL = sweep(CompareBL, 2, theta,"/")

CompareG = sweep(Coefs_GSI,2,theta,"-")
CompareG = sweep(CompareG, 2, theta,"/")

Compare = ( abs(CompareBL) - abs(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)


t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(abs(CompareBL[,2]),abs(CompareG[,2]),paired=T) # GM
t.test(abs(CompareBL[,3]),abs(CompareG[,3]),paired=T) # GM
t.test(abs(CompareBL[,4]),abs(CompareG[,4]),paired=T) # GM 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSLF : Y2 Coefficients  SC ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

lm2 = lm(math~ prog + ses + female, data= ml )
eta = lm2$coefficients 

eta = apply(CoefsOut2_LF_LS_CA_P,2,mean)

Coefs_BLSC= EMeanOut2_LF_LS_CA_P[,1:6]
Coefs_GSC = GMeanOut2_LF_LS_CA_P[,1:6]
Coefs_PBSC= GTMeanOut2_LF_LS_CA_P[,1:6]

# Option 2
CompareBL = sweep(Coefs_BLSC,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")
CompareBL = abs(CompareBL)

CompareG = sweep(Coefs_GSC,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")
CompareG = abs(CompareG)

Compare = ((CompareBL) -(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(CompareBL[,2],CompareG[,2],paired=T) # BL
t.test(CompareBL[,3],CompareG[,3],paired=T) # GM 
t.test(CompareBL[,4],CompareG[,4],paired=T) # NOT SIGNIFICANT 
t.test(CompareBL[,5],CompareG[,5],paired=T) # BL
t.test(CompareBL[,6],CompareG[,6],paired=T) # GM   

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSLF : Y2 Coefficients  W ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

eta = apply(CoefsOut2_LF_LS_D_P,2,mean)

Coefs_BLW= EMeanOut2_LF_LS_D_P[,1:6]
Coefs_GW = GMeanOut2_LF_LS_D_P[,1:6]
Coefs_PBW= GTMeanOut2_LF_LS_D_P[,1:6]

# Option 2
CompareBL = sweep(Coefs_BLW,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")
CompareBL = abs(CompareBL)

CompareG = sweep(Coefs_GW,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")
CompareG = abs(CompareG)

Compare = ((CompareBL) -(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(CompareBL[,2],CompareG[,2],paired=T) # NOT SIGNIFICANT 
t.test(CompareBL[,3],CompareG[,3],paired=T) # NOT SIGNIFICANT 
t.test(CompareBL[,4],CompareG[,4],paired=T) # GM  
t.test(CompareBL[,5],CompareG[,5],paired=T) # BL
t.test(CompareBL[,6],CompareG[,6],paired=T) # BL   


t.test(Coefs_BLW[,1],Coefs_GW[,1],paired=T) # GM
t.test(Coefs_BLW[,2],Coefs_GW[,2],paired=T) # NOT
t.test(Coefs_BLW[,3],Coefs_GW[,3],paired=T) # NOT
t.test(Coefs_BLW[,4],Coefs_GW[,4],paired=T) # GM
t.test(Coefs_BLW[,5],Coefs_GW[,5],paired=T) # BL
t.test(Coefs_BLW[,6],Coefs_GW[,6],paired=T) # BL  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### LSLF : Y2 Coefficients  SI ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

eta = apply(CoefsOut2_LF_LS_CP_P,2,mean)

Coefs_BLSI= EMeanOut2_LF_LS_CP_P[,1:6]
Coefs_GSI = GMeanOut2_LF_LS_CP_P[,1:6]
Coefs_PBSI= GTMeanOut2_LF_LS_CP_P[,1:6]

# Option 2
CompareBL = sweep(Coefs_BLSI,2,eta,"-")
CompareBL = sweep(CompareBL, 2, eta,"/")
CompareBL = abs(CompareBL)

CompareG = sweep(Coefs_GSI,2,eta,"-")
CompareG = sweep(CompareG, 2, eta,"/")
CompareG = abs(CompareG)

Compare = ((CompareBL) -(CompareG))

apply(Compare*100,2,mean)
apply(Compare*100,2,sd)

t.test(CompareBL[,1],CompareG[,1],paired=T) # GM
t.test(CompareBL[,2],CompareG[,2],paired=T) # GM 
t.test(CompareBL[,3],CompareG[,3],paired=T) # GM
t.test(CompareBL[,4],CompareG[,4],paired=T) # GM  
t.test(CompareBL[,5],CompareG[,5],paired=T) # BL
t.test(CompareBL[,6],CompareG[,6],paired=T) # GM    


# All significant 
t.test(Coefs_BLSI[,1],Coefs_GSI[,1],paired=T) # 
t.test(Coefs_BLSI[,2],Coefs_GSI[,2],paired=T) # 
t.test(Coefs_BLSI[,3],Coefs_GSI[,3],paired=T) # 
t.test(Coefs_BLSI[,4],Coefs_GSI[,4],paired=T) # 
t.test(Coefs_BLSI[,5],Coefs_GSI[,5],paired=T) # 
t.test(Coefs_BLSI[,6],Coefs_GSI[,6],paired=T) #   


#%%%%%%%%%%%%%%%%%%%%%%%%%##
####     LSLF : RMSE     ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%##

RMSE_BLSC= sqrt(E_MSE_LF_LS_CA_P)
RMSE_GSC= sqrt(G_MSE_LF_LS_CA_P)
RMSE_PBSC= sqrt(GT_MSE_LF_LS_CA_P)

RMSE_BLW= sqrt(E_MSE_LF_LS_D_P)
RMSE_GW= sqrt(G_MSE_LF_LS_D_P)
RMSE_PBW= sqrt(GT_MSE_LF_LS_D_P)

RMSE_BLSI= sqrt(E_MSE_LF_LS_CP_P)
RMSE_GSI= sqrt(G_MSE_LF_LS_CP_P)
RMSE_PBSI= sqrt(GT_MSE_LF_LS_CP_P)

# All significant 
t.test(RMSE_BLSC,RMSE_GSC,paired=T) 
t.test(RMSE_BLW,RMSE_GW,paired=T)
t.test(RMSE_BLSI,RMSE_GSI,paired=T)

# Option 2 
mean( (RMSE_BLSC/RMSE_PBSC) - (RMSE_GSC/RMSE_PBSC) )*100
mean( (RMSE_BLW/RMSE_PBW) - (RMSE_GW/RMSE_PBW) )*100
mean( (RMSE_BLSI/RMSE_PBSI) - (RMSE_GSI/RMSE_PBSI) )*100

sd( (RMSE_BLSC/RMSE_PBSC) - (RMSE_GSC/RMSE_PBSC) )*100
sd( (RMSE_BLW/RMSE_PBW) - (RMSE_GW/RMSE_PBW) )*100
sd( (RMSE_BLSI/RMSE_PBSI) - (RMSE_GSI/RMSE_PBSI) )*100

##%%% %%%%%%%%%%%%%%%%%%%% ## 
###  Plots : Match Rate #### 
## %%%%%%%%%%%%%%%%%%%%%%% ## 

# Plot for paper
pdf("Paper_LFLS_Match.pdf",width = 11, height = 7)
par(mar=c(5,3,2,2)+0.1)
boxplot((GTMatchOut_LF_LS_CA_N-.2)/.8*100,(GTMatchOut_LF_LS_D_N-.2)/.8*100,(GTMatchOut_LF_LS_CP_N-.2)/.8*100,(GMatchOut_LF_LS_CA_P-.2)/.8*100,(GMatchOut_LF_LS_D_P-.2)/.8*100,(GMatchOut_LF_LS_CP_P-.2)/.8*100,(EMatchOut_LF_LS_CA_P-.2)/.8*100, (EMatchOut_LF_LS_D_P-.2)/.8*100,(EMatchOut_LF_LS_CP_P-.2)/.8*100,names=c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"),ylab="NonSeed Posterior Match Rate")
dev.off()

##%%% %%%%%%%%%%%%%%%%%%%% ## 
###      Plots : RMSE    #### 
## %%%%%%%%%%%%%%%%%%%%%%% ## 

pdf("Paper_MSE.pdf",width = 11, height = 7)
par(mar=c(5,3,2,2)+0.1)
boxplot(sqrt(GT_MSE_LF_LS_CA_N),sqrt(GT_MSE_LF_LS_D_N),sqrt(GT_MSE_LF_LS_CP_N),sqrt(G_MSE_LF_LS_CA_P), sqrt(G_MSE_LF_LS_D_P),sqrt(G_MSE_LF_LS_CP_P),sqrt(E_MSE_LF_LS_CA_P),sqrt(E_MSE_LF_LS_D_P),sqrt(E_MSE_LF_LS_CP_P), ylab= "Root MSE", names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"))
dev.off()

###### Plot: Parameter Estimates ########

install.packages("plotrix")
require(plotrix)

meanVec = c( mean(GTMeanOut_LF_LS_CA_N[,1]),mean(GTMeanOut_LF_LS_D_N[,1]),mean(GTMeanOut_LF_LS_CP_N[,1]),mean(GMeanOut_LF_LS_CA_P[,1]),mean(GMeanOut_LF_LS_D_P[,1]),mean(GMeanOut_LF_LS_CP_P[,1]),mean(EMeanOut_LF_LS_CA_P[,1]),mean(EMeanOut_LF_LS_D_P[,1]),mean(EMeanOut_LF_LS_CP_P[,1]))

li1 = c( quantile(GTMeanOut_LF_LS_CA_N[,1],c(0.025)),quantile(GTMeanOut_LF_LS_D_N[,1],c(0.025)),quantile(GTMeanOut_LF_LS_CP_N[,1],c(0.025)),quantile(GMeanOut_LF_LS_CA_P[,1],c(0.025)),quantile(GMeanOut_LF_LS_D_P[,1],c(0.025)),quantile(GMeanOut_LF_LS_CP_P[,1],c(0.025)),quantile(EMeanOut_LF_LS_CA_P[,1],c(0.025)),quantile(EMeanOut_LF_LS_D_P[,1],c(0.025)),quantile(EMeanOut_LF_LS_CP_P[,1],c(0.025)))

ui1 = c( quantile(GTMeanOut_LF_LS_CA_N[,1],c(0.975)),quantile(GTMeanOut_LF_LS_D_N[,1],c(0.975)),quantile(GTMeanOut_LF_LS_CP_N[,1],c(0.975)),quantile(GMeanOut_LF_LS_CA_P[,1],c(0.975)),quantile(GMeanOut_LF_LS_D_P[,1],c(0.975)),quantile(GMeanOut_LF_LS_CP_P[,1],c(0.975)),quantile(EMeanOut_LF_LS_CA_P[,1],c(0.975)),quantile(EMeanOut_LF_LS_D_P[,1],c(0.975)),quantile(EMeanOut_LF_LS_CP_P[,1],c(0.975)))

pdf('Paper_LF_LS_All_Intercept_CIs.pdf',width = 8, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeLS and axes 
 boxplot(  GTMeanOut_LF_LS_CA_N[,1], GTMeanOut_LF_LS_D_P[,1], GTMeanOut_LF_LS_CP_N[,1],GMeanOut_LF_LS_CA_P[,1], GMeanOut_LF_LS_D_P[,1], GMeanOut_LF_LS_CP_P[,1], EMeanOut_LF_LS_CA_P[,1], EMeanOut_LF_LS_D_P[,1], EMeanOut_LF_LS_CP_P[,1], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x= c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[1],lty = 2)
dev.off()

meanVec = c( mean(GTMeanOut_LF_LS_CA_N[,2]),mean(GTMeanOut_LF_LS_D_N[,2]),mean(GTMeanOut_LF_LS_CP_N[,2]),mean(GMeanOut_LF_LS_CA_P[,2]),mean(GMeanOut_LF_LS_D_P[,2]),mean(GMeanOut_LF_LS_CP_P[,2]),mean(EMeanOut_LF_LS_CA_P[,2]),mean(EMeanOut_LF_LS_D_P[,2]),mean(EMeanOut_LF_LS_CP_P[,2]))

li1 = c( quantile(GTMeanOut_LF_LS_CA_N[,2],c(0.025)),quantile(GTMeanOut_LF_LS_D_N[,2],c(0.025)),quantile(GTMeanOut_LF_LS_CP_N[,2],c(0.025)),quantile(GMeanOut_LF_LS_CA_P[,2],c(0.025)),quantile(GMeanOut_LF_LS_D_P[,2],c(0.025)),quantile(GMeanOut_LF_LS_CP_P[,2],c(0.025)),quantile(EMeanOut_LF_LS_CA_P[,2],c(0.025)),quantile(EMeanOut_LF_LS_D_P[,2],c(0.025)),quantile(EMeanOut_LF_LS_CP_P[,2],c(0.025)))

ui1 = c( quantile(GTMeanOut_LF_LS_CA_N[,2],c(0.975)),quantile(GTMeanOut_LF_LS_D_N[,2],c(0.975)),quantile(GTMeanOut_LF_LS_CP_N[,2],c(0.975)),quantile(GMeanOut_LF_LS_CA_P[,2],c(0.975)),quantile(GMeanOut_LF_LS_D_P[,2],c(0.975)),quantile(GMeanOut_LF_LS_CP_P[,2],c(0.975)),quantile(EMeanOut_LF_LS_CA_P[,2],c(0.975)),quantile(EMeanOut_LF_LS_D_P[,2],c(0.975)),quantile(EMeanOut_LF_LS_CP_P[,2],c(0.975)))

pdf('Paper_LF_LS_All_Math_CIs.pdf',width = 8, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeLS and axes 
 boxplot(  GTMeanOut_LF_LS_CA_N[,2], GTMeanOut_LF_LS_D_N[,2], GTMeanOut_LF_LS_CP_N[,2],GMeanOut_LF_LS_CA_P[,2], GMeanOut_LF_LS_D_P[,2], GMeanOut_LF_LS_CP_P[,2], EMeanOut_LF_LS_CA_P[,2], EMeanOut_LF_LS_D_P[,2], EMeanOut_LF_LS_CP_P[,2], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white",ylim = c(0.2,0.7), at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x=c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[2],lty = 2)
dev.off()

meanVec = c( mean(GTMeanOut_LF_LS_CA_N[,3]),mean(GTMeanOut_LF_LS_D_N[,3]),mean(GTMeanOut_LF_LS_CP_N[,3]),mean(GMeanOut_LF_LS_CA_P[,3]),mean(GMeanOut_LF_LS_D_P[,3]),mean(GMeanOut_LF_LS_CP_P[,3]),mean(EMeanOut_LF_LS_CA_P[,3]),mean(EMeanOut_LF_LS_D_P[,3]),mean(EMeanOut_LF_LS_CP_P[,3]))

li1 = c( quantile(GTMeanOut_LF_LS_CA_N[,3],c(0.025)),quantile(GTMeanOut_LF_LS_D_N[,3],c(0.025)),quantile(GTMeanOut_LF_LS_CP_N[,3],c(0.025)),quantile(GMeanOut_LF_LS_CA_P[,3],c(0.025)),quantile(GMeanOut_LF_LS_D_P[,3],c(0.025)),quantile(GMeanOut_LF_LS_CP_P[,3],c(0.025)),quantile(EMeanOut_LF_LS_CA_P[,3],c(0.025)),quantile(EMeanOut_LF_LS_D_P[,3],c(0.025)),quantile(EMeanOut_LF_LS_CP_P[,3],c(0.025)))

ui1 = c( quantile(GTMeanOut_LF_LS_CA_N[,3],c(0.975)),quantile(GTMeanOut_LF_LS_D_N[,3],c(0.975)),quantile(GTMeanOut_LF_LS_CP_N[,3],c(0.975)),quantile(GMeanOut_LF_LS_CA_P[,3],c(0.975)),quantile(GMeanOut_LF_LS_D_P[,3],c(0.975)),quantile(GMeanOut_LF_LS_CP_P[,3],c(0.975)),quantile(EMeanOut_LF_LS_CA_P[,3],c(0.975)),quantile(EMeanOut_LF_LS_D_P[,3],c(0.975)),quantile(EMeanOut_LF_LS_CP_P[,3],c(0.975)))

pdf('Paper_LF_LS_All_Academic_CIs.pdf',width = 8, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeLS and axes 
 boxplot(  GTMeanOut_LF_LS_CA_N[,3], GTMeanOut_LF_LS_D_N[,3], GTMeanOut_LF_LS_CP_N[,3],GMeanOut_LF_LS_CA_P[,3], GMeanOut_LF_LS_D_P[,3], GMeanOut_LF_LS_CP_P[,3], EMeanOut_LF_LS_CA_P[,3], EMeanOut_LF_LS_D_P[,3], EMeanOut_LF_LS_CP_P[,3], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x=c(2:4,6:8,10:12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[3],lty = 2)
dev.off()


meanVec = c( mean(GTMeanOut_LF_LS_CA_N[,4]),mean(GTMeanOut_LF_LS_D_N[,4]),mean(GTMeanOut_LF_LS_CP_N[,4]),mean(GMeanOut_LF_LS_CA_P[,4]),mean(GMeanOut_LF_LS_D_P[,4]),mean(GMeanOut_LF_LS_CP_P[,4]),mean(EMeanOut_LF_LS_CA_P[,4]),mean(EMeanOut_LF_LS_D_P[,4]),mean(EMeanOut_LF_LS_CP_P[,4]))

li1 = c( quantile(GTMeanOut_LF_LS_CA_N[,4],c(0.025)),quantile(GTMeanOut_LF_LS_D_N[,4],c(0.025)),quantile(GTMeanOut_LF_LS_CP_N[,4],c(0.025)),quantile(GMeanOut_LF_LS_CA_P[,4],c(0.025)),quantile(GMeanOut_LF_LS_D_P[,4],c(0.025)),quantile(GMeanOut_LF_LS_CP_P[,4],c(0.025)),quantile(EMeanOut_LF_LS_CA_P[,4],c(0.025)),quantile(EMeanOut_LF_LS_D_P[,4],c(0.025)),quantile(EMeanOut_LF_LS_CP_P[,4],c(0.025)))

ui1 = c( quantile(GTMeanOut_LF_LS_CA_N[,4],c(0.975)),quantile(GTMeanOut_LF_LS_D_N[,4],c(0.975)),quantile(GTMeanOut_LF_LS_CP_N[,4],c(0.975)),quantile(GMeanOut_LF_LS_CA_P[,4],c(0.975)),quantile(GMeanOut_LF_LS_D_P[,4],c(0.975)),quantile(GMeanOut_LF_LS_CP_P[,4],c(0.975)),quantile(EMeanOut_LF_LS_CA_P[,4],c(0.975)),quantile(EMeanOut_LF_LS_D_P[,4],c(0.975)),quantile(EMeanOut_LF_LS_CP_P[,4],c(0.975)))

pdf('Paper_LF_LS_All_Vocational_CIs.pdf',width = 9, height = 7)
par(mar=c(5,3,2,2)+0.1)
# Plot a blank boxplot graph to set the labeLS and axes 
 boxplot(  GTMeanOut_LF_LS_CA_N[,4], GTMeanOut_LF_LS_D_N[,4], GTMeanOut_LF_LS_CP_N[,4],GMeanOut_LF_LS_CA_P[,4], GMeanOut_LF_LS_D_P[,4], GMeanOut_LF_LS_CP_P[,4], EMeanOut_LF_LS_CA_P[,4], EMeanOut_LF_LS_D_P[,4], EMeanOut_LF_LS_CP_P[,4], names = c("PB","PB","PB","GM","GM","GM","BL.CA","BL.D","BL.CP"), ylab = "Posterior Means",border="white", at = c(2,3,4,6,7,8,10,11,12))
# add on the CIS
plotCI(x= c(2,3,4,6,7,8,10,11,12),y = meanVec, ui=ui1, li=li1,ylab = "Posterior Means", xlab = "",add=TRUE)
abline(h = GenCoefs[4],lty = 2)
dev.off()
