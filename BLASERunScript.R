## %%%%%%%%%%%%%%%%%%% ##
##  Read User Inputs   ##
## %%%%%%%%%%%%%%%%%%% ## 

# Put the specific name of the input file you want here
source("Inputs.R") 

if( class(MCMCSeed) != "numeric"){
  MCMCSeed <-as.numeric(Sys.time())
  MCMCSeed = (MCMCSeed - floor(MCMCSeed))*1e8
  warning(paste("No random seed specified, seed of", MCMCSeed, "utilized."))
}

## %%%%%%%%%%%%%%%%%%% ## 
##    Libraries        ## 
## %%%%%%%%%%%%%%%%%%% ##  

library(NPBayesImpute)

## %%%%%%%%%%%%%%%%%%%%%% ## 
## Create the Data/ Store ##
## %%%%%%%%%%%%%%%%%%%%%% ## 

source("Functions.R")

## %%%%%%%%%%%%%%%%%%%%%% ## 
## Create the Data/ Store ##
## %%%%%%%%%%%%%%%%%%%%%% ## 

source("DataGen.R")

# Define the name
namepart = paste(perc.type1,"T1_0T2_",error.perc,"progError_Run",sep="")

# Set/Create the Working Directory
dir.create(file.path(mainDir=getwd(), subDir=namepart))
setwd(file.path(mainDir=getwd(), subDir=namepart) )

# Sample the Type 1 seeds
seeds.orig = sort(type1seeds)
perc.in = 100*length(seeds.orig)/N 
perc.needed = perc.type1 - perc.in 
if( perc.type1 != 0 ){
	seeds = sample( c(1:nrow(File1))[-seeds.orig],perc.needed/100*nrow(File1),replace=F)
	seeds = sort( c(seeds.orig, seeds) )
} else{ 
 	perc.type1 = perc.in
        seeds      = seeds.orig 
}

# %%%%%%%%%%%%%%%%%%%%%%%%% #  
# What are the type 1 seeds?  
# %%%%%%%%%%%%%%%%%%%%%%%%% #

matches <-matrix(NA,nrow = length(seeds), ncol = 2)
matches[,2] = matches[,1] = seeds

## %%%%%%%%%%%%%%%%%%%%% # 
##     Pointers         ##
## %%%%%%%%%%%%%%%%%%%%% # 

Y1.name = "read"
Y2.name = "math"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
####            Run PB            #####
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

method = "GOnly"

run_PB(MCMCSeed = MCMCSeed , File1,File2, method,its,burnin,gap, reps, thinning, printgap, Comment, coefs.primary, coefs.primary.all, coefs.secondary, col.Y1, col.Y2, secondary.option, secondary.option2,error.file, error.fields, blockingVars, matches, Y1.name, Y2.name,paste("PB",namepart,sep=""),threshold,complete.out,type2seeds= NULL)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
####     Introduce Block Faults       #### 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##   

Error.Options =  setdiff( 1:N, seeds)
in.error       = sample( Error.Options, N*error.perc/100, replace = F)
rm(Error.Options)

field.options = levels(File1[,error.fields])

for( i in in.error){
  
  if( error.file == 1 ){
    current.field   = File1[ i ,error.fields]
    field.choices   = setdiff( field.options,current.field )
    File1[i,"prog"] = sample(field.choices,1)
  } else{ 
    current.field   = File2[ i ,error.fields]
    field.choices   = setdiff( field.options,current.field )
    File2[i,"prog"] = sample(field.choices,1)
  }
    
}

rm(i,current.field,field.choices,field.options)

## %%%%%%%%%%%%%%%%%%%%%% ## 
####  Check the Inputs ####
## %%%%%%%%%%%%%%%%%%%%% ##

realmin =  1e-08

if( class(File1) != "data.frame" ){
  stop('File1 must be a data frame.')
}

if( class(File2) != "data.frame" ){
  stop('File2 must be a data frame.')
}

if( sum(matches !=0) > 0 ){
  if( class(matches) != "matrix" | ncol(matches) != 2 ){
    stop('matches must be a matrix holding in the type 1 seeds in which the first column denotes the record in File1, the second the corresponding record in File2.')
  }
}

if( exists("blockingVars") == "FALSE" ){
  stop('blockingVars is empty, ie there are no specified blocking variables.')
}

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
## Extract the Global Constants ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

n1 = nrow(File1)
n2 = nrow(File2)
N  = n1 + n2 # The possible number of unique individuals in the two files 

# Subtract 1 from the number of latent classes
HstarM1 = Hstar - 1 

# In the linked file, in which columns will Y1 and Y2 be stored? 
if( exists("which.Y1") == FALSE){
  which.Y1 = 1 
  which.Y2 = 2 
}

# In their respective files, in which columns are Y1 and Y2 stored 
col.Y1   = which(colnames(File1) == Y1.name)
col.Y2   = which(colnames(File2) == Y2.name) 

## %%%%%%%%%%%%%%%%%%%%%%% # 
## The Blocking Variables ##
## %%%%%%%%%%%%%%%%%%%%%%% #

# Determine how many blocking variables there are
J = length(blockingVars)

if( J == 0 ){
  stop('blockingVars is empty, ie there are no specified blocking variables.')
}

where.bvars.F1 = which(colnames(File1) %in% blockingVars)
where.bvars.F2 = which(colnames(File2) %in% blockingVars)

if( exists("where.bvars.F1") == "FALSE"){
  stop('The specified blocking variables are not present in File 1.')
}

if( length(where.bvars.F1) != J ){
  stop('Some of the specified blocking variables are not present in File 1.')
}

if( exists("where.bvars.F2") == "FALSE"){
  stop('The specified blocking variables are not present in File 2.')
}

if( length(where.bvars.F2) != J ){
  stop('Some of the specified blocking variables are not present in File 2.')
}

if( sum( blockingVars == names(File1)[where.bvars.F1]) != J ){
  # The order of the columns need to be changed.
  File1[,where.bvars.F1] = File1[,blockingVars]
  names(File1)[where.bvars.F1] <- blockingVars
  warning('File 1 columns re-ordered to match the input blocking variable ordering.')
}

if( sum( blockingVars == names(File2)[where.bvars.F2]) != J ){
  # The order of the columns need to be changed.
  File2[,where.bvars.F2] = File2[,blockingVars]
  names(File2)[where.bvars.F2] <- blockingVars
  warning('File 2 columns re-ordered to match the input blocking variable ordering.')
}

## %%%%%%%%%%%%% # 
##  Type 1 Seeds # 
## %%%%%%%%%%%%% # 

if( exists("error.file") == FALSE ){
  stop('error.file undefined. Please specify which file may have noisy blocking variables.')
}

if( sum(matches !=0) > 0  ){
  
  type1seeds = matches[,2]
  
  n.type1    = length(type1seeds)
} else{ 
  type1seeds <- NULL
  n.type1 = 0 
}

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
## Create Bhat, Bhat1, Bhat2     # 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 

# We need Bhat to be a data frame with each column a factor.

Bhat1          = File1[,where.bvars.F1]
Bhat2          = File2[,where.bvars.F2]
Bhat           = rbind(Bhat1,Bhat2)
rownames(Bhat) = 1:N

B     <- Bhat 

d = sapply(Bhat,function(x) length(unique(x)) )


## %%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
##  For the Dirichlet Process  ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# Is the number of blocking variables different from those used in the DP? 
if( exists("DPOffset") == "FALSE"){
  stop('DOffset not specified')
}

if( class(Bhat) != "data.frame"){
  stop('The blocking variables were not successfully stored in a data frame.')
}
  
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
##  Assigning Observed Blocks  ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% #

blocks <-NULL 

for( i in 1:N){
  blocks[i] = category.block(B[i,],d,J)
}

blocks1 = blocks[1:n1]
blocks2 = blocks[(n1+1):N]


if( sum(matches !=0) > 0 ){
  if( sum(blocks1[type1seeds] == blocks2[type1seeds])!= n.type1){
    stop('Some declared Type 1 Seed pairs are not in the same block.')
  }
}

K = length(unique(blocks))

observed.blocks = sort(blocks)

## %%%%%%%%%%%%%%%%%%% # 
## Block Restriction  ##
## %%%%%%%%%%%%%%%%%%% #

if( exists("ObsBlock_restriction") == FALSE ){
  warning('Block restriction not specified. By default, moves are restricted to observed blocks.')
  ObsBlock_restriction = "yes"
} 

if( ObsBlock_restriction == "yes"){
  # We restrict possible moves to observed blocks.
  # Convert from indices to blocks 
  Block_Ops      = blocks 
  Index_to_Block = Legal_Index = sort( unique(Block_Ops) )
  Block_Ops      = as.numeric(as.factor(Block_Ops))
  indices1       = blocks1
  indices2       = blocks2
  blocks1        = Block_Ops[1:n1]
  blocks2        = Block_Ops[(n1+1):N]
  K              = length(unique(blocks))
  
  if( max( max(blocks1),max(blocks2)) != K){
    stop('Observed blocking restriction error.')
  }
  
  
  if( sum(matches !=0) > 0  ){
    if( sum(blocks1[type1seeds] == blocks2[type1seeds])!= n.type1){
      stop('Some declared Type 1 Seed pairs are not in the same block post restriction.')
    }
  }
}

## %%%%%%%%%%%%%%%%%%%%% # 
## Records in the Block ##
## %%%%%%%%%%%%%%%%%%%%% #

n.RIB = matrix( NA, nrow= 2, ncol = K )
RIB   <- vector("list", length = 2 )
RIB[[1]] <- vector( "list", length = K )
RIB[[2]] <- vector( "list", length = K )

for( k in 1:K){
  RIB[[1]][[k]] = which(blocks1==k)
  RIB[[2]][[k]] = which(blocks2==k )
  
}
rm(k)

n.RIB[1,] = sapply(RIB[[1]],length)
n.RIB[2,] = sapply(RIB[[2]],length)
  
  
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
## Format File1 and File 2 for the regression ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

File1_data = File1
File2_data = File2 


coefs.all = union(coefs.primary,coefs.secondary)
f1 <- as.formula(paste("File1[,",col.Y1,"]~", paste(coefs.all, collapse="+")))
modelPart <- model.matrix( lm(File1[,col.Y1]~File2[,col.Y2] + prog,data=File1))

File1_data = cbind( File1[,col.Y1], modelPart[,-c(1:2)])
colnames(File1_data)[col.Y1] = Y1.name
modelPart <- model.matrix( lm(File1[,col.Y1]~File2[,col.Y2] + prog,data=File2))
modelPartPrimary = modelPart 
File2_data = cbind( File2[,col.Y2], modelPart[,-c(1:2)])
colnames(File2_data)[col.Y2] = Y2.name

# Create the Secondary Model Matrix 
f1 <- as.formula(paste("File2[,",col.Y2,"]~", paste(coefs.secondary, collapse="+")))
modelPartSecondary <- model.matrix( lm(f1,data=File2))
rm(f1)

coefs.all = union(coefs.primary.all,coefs.all)

# Identify which of these are primary
where.coefs.primary = which(coefs.all %in% coefs.primary.all)
# To count the columns, consult the levels associated with each. 
where.coefs.primary = where.coefs.primary[1]:(where.coefs.primary[1]+sum(d[coefs.primary])-length(coefs.primary))

# Identify which of these are secondary
where.coefs.secondary = which( coefs.all %in% coefs.secondary)
# To count the columns, consult the levels associated with each. 
where.coefs.secondary = where.coefs.secondary[1]:(1+sum(d[coefs.secondary])-length(coefs.secondary))

where.coefs = sort(union(where.coefs.primary,where.coefs.secondary))

p  = length(where.coefs.primary)+1
p2 = length(where.coefs.secondary)+1


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
## Store the Seed and Block Information ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

BlockRow   = Seed = Block= 0
File1_data = cbind( File1_data, BlockRow, Seed,Block)
File2_data = cbind( File2_data, BlockRow, Seed, Block)

# Exact Matches
File1_data[type1seeds,"Seed"] = 1
File2_data[type1seeds,"Seed"] = 1

possibleEthnic <- how_many_possible( B, N = N/2 , d, field = which(blockingVars == error.fields), J, type1seeds, Legal_Index,blocks)
  
possible.length.E = unlist( lapply( possibleEthnic, length ))
type2seeds = which(possible.length.E==1)
type2seeds = setdiff( type2seeds, type1seeds)

# Perfect Blocking 
if( length(type2seeds) > 0 ){
  File1_data[type2seeds,"Seed"] = 2
  File2_data[type2seeds,"Seed"] = 2
}


# Which records are now candidates to be moved?
if( error.file == 1){
  MoveCand_F2  = which(File1_data[,"Seed"]==0)
  n.error2     = length(MoveCand_F2)
} else{ 
  MoveCand_F2  = which(File2_data[,"Seed"]==0)
  n.error2     = length(MoveCand_F2)
}

File1_data[,"Block"] = blocks1
File2_data[,"Block"] = blocks2

count.matches   = rep(0,K)

if( sum(matches !=0) > 0 ){
  if( sum(File1_data[type1seeds,"Block"] == File2_data[type1seeds,"Block"])!= n.type1){
    stop('Some declared Type 1 Seed pairs are not in the same block post block storage.')
  }
  
  # How many matches in each block? 
  for( k in 1:K){
    count.matches[k] = sum(File1_data[type1seeds,"Block"]==k)
  }
  rm(k)
} 

## %%%%%%%%%%%%%%%%% ##
##  Format the Data  ## 
## %%%%%%%%%%%%%%%%% ##

record.holder  = records.in.block( File1_data, File2_data,K )
n.RIB          = record.holder$n.RIB
RIB            = record.holder$RIB
rm(record.holder)

lambda = 1:n1
Imputed = 0 
File1_data = cbind( File1_data,lambda,Imputed )
lambda = 1:n2
File2_data = cbind( File2_data, lambda,Imputed)
rm(lambda,Imputed)

FileHolder  = rbind(File1_data,File2_data)
FileHolder  = FileHolder[order(FileHolder[,"Block"]),]
Imp.Holder1 = make_imputed_holder( FileHolder,1 )
Imp.Holder2 = make_imputed_holder( FileHolder,1 )
Holder      = add_imputed_rows( n.RIB, Imp.Holder1, Imp.Holder2, RIB,File1_data,File2_data)
rm(FileHolder)
File1_data  = Holder$File1_data
File2_data  = Holder$File2_data

rm(Holder)

File1_data[,"lambda"] = c(1:nrow(File1_data))
File2_data[,"lambda"] = c(1:nrow(File2_data))

nrow.F1 <<- nrow(File1_data)
nrow.F2 <<- nrow(File2_data)
rownames(File1_data) = 1:nrow.F1  

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
## Save the formatted data sets  ## 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

save(File1_data,file ="File1Formatted.RData")
save(File2_data,file ="File2Formatted.RData")
  save.image(paste(namepart,"PreMethod.RData"))


# %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
#####       Run GM           ####
## %%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

n = nrow.F1
impute_list<-list()
impute_list[[1]] = 1
impute_list[[2]] =1 

record.holder  = G.records.in.block( File1_data, File2_data,K,type1seeds )
n.RIB          = record.holder$n.RIB
RIB            = record.holder$RIB
rm(record.holder)

PermSampling = ifelse( n.RIB[1,] > 1, ifelse(n.RIB[1,] < threshold, 1, 2), 0 )
SeedsOnly    = which(n.RIB==0)
needs.C      = which(PermSampling>0)

holder     = G.get_starting_permutation(K, RIB, File1_data, File2_data,nrow.F1,nrow.F2,type1seeds)
File1_data = holder$File1
File2_data = holder$File2
rm(holder) 

nonseed.F2 = which(File2_data[,"Seed"]==0)
n.nonseed.F2 = length(nonseed.F2)

C = diag(p)

# Run the script 
run_GM(MCMCSeed, its,burnin, thinning, gap, n.RIB, RIB, File1_data, File2_data,K,which.Y2, which.Y1,type1seeds, where.coefs.primary, where.coefs.secondary, p, p2, n.type1,needs.C,reps,PermSampling,secondary.option="dependent", n1, n2,col.Y1, col.Y2,paste("GM",namepart,sep=""))


## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ### 
#####        Run BLASE         ####
## %%%%%%%%%%%%%%%%%%%%%%%%%%% #### 

run_BLASE()








