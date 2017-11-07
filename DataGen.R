# Create a script to automatically generate the data. 

## %%%%%%%%%%%%%%%%%%%%%% ## 
## Load Required Packages ## 
## %%%%%%%%%%%%%%%%%%%%%% ## 

require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)

## %%%%%%%%%%%%%%%%%%%%%% ## 
##      Load the Data     ## 
## %%%%%%%%%%%%%%%%%%%%%% ## 

ml <- read.csv("ml.csv")

## %%%%%%%%%%%%%%%%%%%%%% ## 
##      Load Functions    ## 
## %%%%%%%%%%%%%%%%%%%%%% ## 

summarize_BlockMembership<-function( Binput , NonSeed, d, K,J){
  
  # Identify the block membership for all records in the input. 
  membership = identify.blocks( Binput ,d, J)
  
  # Look at the seeds; count the number in each block
  block.seeds = rep(0,K)
  for( i in NonSeed){
    k = membership[i]
    block.seeds[k] = block.seeds[k] + 1 
  }
  
  block.seeds = subset(block.seeds, block.seeds > 0)
  
  return( block.seeds )
  
}

identify.blocks<-function(B,d,J){
  # Get the block indicator
  block.cat <-list()
  if( is.list(B)==TRUE){
    for( f in 1:2){
      block.cat[[f]] = apply( as.matrix(B[[f]]), 1, function(x) category.block( x,d,J ) )
    }
  }
  if( is.list(B)==FALSE){
    block.cat = apply( as.matrix(B), 1, function(x) category.block( x,d,J ) )
  }
  return( block.cat )
}

category.block <-function(x,d,J){
  x = as.numeric(x)
  category=0
  for (j in 1:(J-1)) category=category+prod( d[(j+1):J])*(x[j]-1)
  category = category+x[J]
  return(category)
}

## %%%%%%%%%%%%%%%%%%%%%% ## 
##      Generate Sex      ## 
## %%%%%%%%%%%%%%%%%%%%%% ## 

female = sample( levels(ml$female), N, replace = T, prob = table(ml$female)/200)

## %%%%%%%%%%%%%%%%%%%%%% ## 
##      Generate CID      ## 
## %%%%%%%%%%%%%%%%%%%%%% ## 

cid = sample( 1:30, N, replace = T)

## %%%%%%%%%%%%%%%%%%%%%% ## 
##    Generate Schtyp     ## 
## %%%%%%%%%%%%%%%%%%%%%% ## 

ml.hold = ml[,c("schtyp","cid")]
#table( ml.hold )
# What we see from the table is that there are a few combinations of schtyp and cid which do not occur. We reintroduce at least some possibility of seeing this combination into the data. 
ml.hold = rbind( ml.hold, c("private",1))
ml.hold = rbind( ml.hold, c("private",2))
ml.hold = rbind( ml.hold, c("private",9))
ml.hold = rbind( ml.hold, c("private",14))
ml.hold = rbind( ml.hold, c("private",17))

#test = multinom(schtyp ~ as.factor(cid), data = ml.hold)
#dses <- data.frame(as.factor(cid))
#pred.probs = predict(test, newdata = dses, "probs")

#fitted.schtyp = NULL
#for(i in 1:N){
#  fitted.schtyp[i] = sample( levels(ml$schtyp), 1, replace=T, prob = c(1-pred.probs[i],pred.probs[i]))
#}

fitted.schtyp = sample( levels(ml$schtyp), N, replace=T, prob = table(ml.hold[,"schtyp"])/200)
new.B         = cbind(female,cid,fitted.schtyp)
colnames(new.B)[3] = "schtyp"

## %%%%%%%%%%%%%%%%%%%%%% ## 
##    Generate Program    ## 
## %%%%%%%%%%%%%%%%%%%%%% ## 

test = multinom(prog~schtyp, data = ml)
dses <- data.frame(schtyp=as.factor(new.B[,"schtyp"]))
pred.probs = predict(test, newdata = dses, "probs")

fitted.prog = NULL
for(i in 1:N){
  fitted.prog[i] = sample( levels(ml$prog), 1, replace=T, prob = c(pred.probs[i,]))
}

new.B = cbind(new.B,fitted.prog)
colnames(new.B)[4] = "prog"

## %%%%%%%%%%%%%%%%%%%%%% ## 
##      Generate SES      ## 
## %%%%%%%%%%%%%%%%%%%%%% ## 

test = multinom(ses~schtyp+prog, data = ml)
dses <- data.frame(schtyp=as.factor(new.B[,"schtyp"]),prog=new.B[,"prog"])
pred.probs = predict(test, newdata = dses, "probs")

fitted.ses = NULL
for(i in 1:N){
  fitted.ses[i] = sample( levels(ml$ses), 1, replace=T, prob = c(pred.probs[i,]))
}

#table(new.B)
new.B = cbind(new.B,fitted.ses)
colnames(new.B)[5] = "ses"

## %%%%%%%%%%%%%%%%%%%%%% ## 
##    Generate Honors     ## 
## %%%%%%%%%%%%%%%%%%%%%% ## 

test = multinom(honors~ses+prog, data = ml)
dses <- data.frame(ses=as.factor(new.B[,"ses"]),prog=new.B[,"prog"])
pred.probs = predict(test, newdata = dses, "probs")

fitted.honors = NULL
for(i in 1:N){
  fitted.honors[i] = sample( levels(ml$honors), 1, replace=T, prob = c(1-pred.probs[i],pred.probs[i]))
}

new.B = cbind(new.B,fitted.honors)
colnames(new.B)[6] = "honors"

## %%%%%%%%%%%%%%%%%%%%%% ## 
##     Store the BVS      ## 
## %%%%%%%%%%%%%%%%%%%%%% ## 

ml.Ready = apply(new.B,2,function(x) as.numeric(factor(x)) )

B.Indices = identify.blocks(ml.Ready,d=c(2,30,2,3,3,2),J=6)

new.B.plus = cbind(new.B, B.Indices)

## %%%%%%%%%%%%%%%%%%%%%% ## 
##   Delete Superfluous   ## 
## %%%%%%%%%%%%%%%%%%%%%% ## 

rm(dses,ml.hold,cid,female,pred.probs,test)

## %%%%%%%%%%%%%%%%%%%%%% ## 
##     Store the BVS      ## 
## %%%%%%%%%%%%%%%%%%%%%% ##

new.B.df = data.frame(new.B)
new.B.df[,"cid"] = as.numeric(levels(new.B.df[,"cid"]))[as.integer(new.B.df[,"cid"])]

new.B.df$honors<- relevel(new.B.df$honors, ref = "not enrolled")

levels(new.B.df$prog) <- levels(ml$prog)

#levels(ml$female)
new.B.df$female<- relevel(new.B.df$female, ref = "male")
#levels(new.B.df$female)

levels(new.B.df$ses)   <- levels(ml$ses)

levels(new.B.df$schtyp)<- levels(ml$schtyp)

B <- new.B

################################################

# Have NOT edited after this point 

HOLDER = new.B.df
HOLDER[,1] =  as.numeric(new.B.df[,1])
HOLDER[,3] =  as.numeric(new.B.df[,3])
HOLDER[,4] =  as.numeric(new.B.df[,4])
HOLDER[,5] =  as.numeric(new.B.df[,5])
HOLDER[,6] =  as.numeric(new.B.df[,6])

new.B = HOLDER

# Model Choices 
# Primary
Model3 = lm(read~math + prog , data= ml )

# Secondary
Model8 = lm(math~ honors +prog+ses+female, data= ml )

# Generate the Math Data 
Merged.Data         = cbind("math" = 0 ,new.B.df)
Pred.Math.Mat       = model.matrix(Model8,data=data.frame(Merged.Data))
Pred.Math           = rnorm( N, Pred.Math.Mat%*%Model8$coefficients, sd(ml$math)-3)
Merged.Data[,"math"]= Pred.Math

# Generate the Read Data 
Data.Frame.In          = data.frame(Merged.Data)
Data.Frame.In[,"math"] = Pred.Math
Data.Frame.In = cbind("read"=rep(0,5000),Data.Frame.In)
Pred.Read.Mat = model.matrix(Model3,data=Data.Frame.In)
Pred.Read     = rnorm( N, Pred.Read.Mat%*%Model3$coefficients, sd(ml$read)-4)

# %%%%%%%%%%%%%%%%%%%%%%%%%%
# Check to make sure we still have some correlations in place that we would expect. This is the FIXED True Data. 

if( cor(Pred.Math,Pred.Read) > .75 | cor(Pred.Math,Pred.Read) < 0.5){
  stop(paste("The correlations of", cor(Pred.Math,Pred.Read), "is not in the specified regions"))
}

Merged.Data = cbind(Pred.Read, Pred.Math,new.B)
colnames(Merged.Data)[1:2] = c("read","math")
Merged.DataFrame = cbind(Pred.Read,Pred.Math,new.B.df)
colnames(Merged.DataFrame)[1:2] = c("read","math")

#summary( lm(read~math + prog, data = Merged.DataFrame ) )
#summary(Model3)

#summary( lm( math ~ honors + prog + ses + female, data = Merged.DataFrame ) )
#summary(Model8)

## %%%%%%%%%%%%%%%%%%%%%% ## 
##     Format the Data    ## 
## %%%%%%%%%%%%%%%%%%%%%% ##

File1 = Merged.DataFrame[,-2]
File2 = Merged.DataFrame[,-1]

## %%%%%%%%%%%%%%%%%%%%%% ## 
##        Data Check      ## 
## %%%%%%%%%%%%%%%%%%%%%% ##

Merged.Data = as.matrix(Merged.Data)
B.Indices = identify.blocks(Merged.Data[,3:8],d=c(2,30,2,3,3,2),J=6)

over.10 = which(table(B.Indices) >10 ) # Which blocks have more than 11 records?
over.10.names = as.integer( names(over.10) )

#sum( table(B.Indices)[over.10]-10 )# How many over 10? 

type1seeds = rep(0,N)
a = 1 
for( i in over.10.names ){  
  in.block = which(B.Indices==i)
  need.to.sample = table(B.Indices)[over.10[a]]-10
  sampled.choices = sample(1:length(in.block), need.to.sample, replace=F)
  sampled.choices = in.block[sampled.choices]
  type1seeds[sampled.choices] = 1  
  a = a + 1
}

Merged.Data = cbind( Merged.Data, "Index" = B.Indices, "Seed" = type1seeds)
Merged.DataFrame = cbind(Merged.DataFrame,"Index" = B.Indices, "Seed" = type1seeds)

#table( summarize_BlockMembership( Merged.Data[,3:8], c(1:N)[-(which(type1seeds==1))], d=c(2,30,2,3,3,2), prod(c(2,30,2,3,3,2)),J=6) )


rm(Data.Frame.In, Pred.Math.Mat,Pred.Read.Mat,new.B.plus,Model3, Model8, Pred.Math, Pred.Read, a, i, in.block,need.to.sample,over.10,over.10.names, sampled.choices,HOLDER,new.B,new.B.df,ml,ml.Ready,B.Indices)

Seeds = type1seeds
type1seeds = which(Seeds==1)
d = apply(B, 2, function(x) length(unique(x)))

truthSim = list()
truthSim$coefsPrimary = lm(read~math + prog, data = Merged.DataFrame )$coefficients
truthSim$coefsSecondary = lm( math ~ honors + prog + ses + female, data = Merged.DataFrame )$coefficients

rm(fitted.honors,fitted.prog,fitted.schtyp,fitted.ses, Merged.Data, Merged.DataFrame,B)

## %%%%%%%%%%%%%%%%%%%%%% ## 
##     Save the image     ## 
## %%%%%%%%%%%%%%%%%%%%%% ##

save.image(paste("DataIn.RData",sep=""))



