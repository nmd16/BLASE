# Evaluation Script for Multiple Runs

indices  = 1:50 
lastcoef = 5 

LambdaOut <- NULL
EDiffOut  <- NULL
GDiffOut  <- NULL 
EMatchOut <- NULL
GMatchOut <- NULL 
a = 0 
 
for( i in c(1:50)[-c(6,16,21,26,33,40,44)]){ 
  a = a + 1 
  # Load the Data 
  load(paste("DataIn_",i,".RData",sep=""))
  N = nrow(File1)   
  # Compute the Truth 
  truthModel = (lm(File1[,1]~File2[,1]+File1[,"prog"]))
  truth = truthModel$coefficients
  truth = truth[-c(1,lastcoef)] 
  # Load in the Gutman Output 
  load( paste("Run",i,"/GOut",perc.type1,"T1_0T2_",error.perc,"progError_Run",i,".RData",sep="") )
 
  # Read the Extension Output  
  Chain1 <- read.table( paste("Run",i,"/",perc.type1,"T1_0T2_",error.perc,"progError_Run",i,"ThetaOut.txt",sep="") )

  # Read the Lambda Output while we are here 
  Lambda1 <- read.table( paste("Run",i,"/",perc.type1,"T1_0T2_",error.perc,"progError_Run",i,"LambdaOut.txt",sep="") )
   GMean = apply( as.matrix(G.outlist$out.coefs[-c(1:newburn),-c(1,lastcoef)]),2,mean)
  EMean = apply( as.matrix(Chain1[-c(1:newburn),-c(1,lastcoef)]),2,mean)

  #GDiff = abs( GMean - truth )
  #EDiff = abs( EMean - truth )

  GDiff = ( GMean - truth)/abs(truth)
  EDiff = ( EMean - truth )/abs(truth)
  
  GDiffOut <- rbind(GDiffOut,GDiff)
  EDiffOut <- rbind(EDiffOut,EDiff)
  
  if( i ==1){
    # For the first run 
    par(pty="s") # make the axes square 
    plot(GDiffOut[1,]~EDiffOut[1,],col="white",pch = c("M","A","V"),xlim = c(-1,1), ylim= c(-1,1), xlab = "(1/|Matched|)*|Posterior Mean Extension - Matched|", ylab = "(1/|Matched|)*|Posterior Mean Gutman - Matched|")
abline( a = 0, b = 1)
text(EDiffOut[i,], GDiffOut[i,], labels = c('M','A','V'), adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 1, col = NULL)
  } else{ 
    # Add the points for the next run
    points( GDiffOut[a,]~EDiffOut[a,],col="white",pch = c("M","A","V"),xlim = c(-1,1), ylim= c(-1,1), xlab = "(1/|Matched|)*(Posterior Mean Extension - Matched)", ylab = "(1/|Matched|)*(Posterior Mean Gutman - Matched)")
    text(EDiffOut[a,], GDiffOut[a,], labels = c('M','A','V'), adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 1, col = NULL)

  }

  # Compute the Average Match Rate for each simulation 

  HOLDERG = apply( G.outlist$out.lambda[-c(1:newburn),],1, function(x) sum(x[1:N]== 1:N) )
  # Overall Match Rate
  GMatchOut <- c(GMatchOut,mean(HOLDERG/N))
  # Overall Match Rate for the Extension 
  EMatchOut <- c(EMatchOut,mean(as.matrix(Lambda1[-c(1:newburn),]/N)))	
} 

# Compute the AVERAGE difference

Avg.Diff = apply(GDiffOut - EDiffOut, 2, mean) 
print(Avg.Diff) 

# Create a matrix which determines the difference for each simulation
DiffAll = GDiffOut - EDiffOut
colnames(DiffAll) = c("Math","Academic","Vocational")
# Create a Boxplot of those differences 
jpeg('Boxplot.jpg')
boxplot(DiffAll,main = "Average Difference in (1/|Matched|)*|Posterior Mean - Matched|", ylab = "Gutman - Extension")
abline(h= 0)
dev.off()

# Determine which percentage of the simulations showed an increase 
apply( DiffAll, 2, function(x) length(which(x > 0 )))/num.sim

# Determine which simulations showed an increase on which field
which.improve.Math = which(DiffAll[,1] > 0)
which.improve.Ac = which(DiffAll[,2] > 0)
which.improve.Voc = which(DiffAll[,3] > 0)
which.worse.Math = which(DiffAll[,1] < 0)
which.worse.Ac = which(DiffAll[,2] < 0)
which.worse.Voc = which(DiffAll[,3] < 0)

# If we improve, the average improvement is 
mean(DiffAll[which.improve.Math,1])
mean(DiffAll[which.improve.Ac,2])
mean(DiffAll[which.improve.Voc,3])
# If we WORSEN, the average decrease is 
mean(DiffAll[which.worse.Math,1])
mean(DiffAll[which.worse.Ac,2])
mean(DiffAll[which.worse.Voc,3])

# Compute the average match rate for each simulation 
