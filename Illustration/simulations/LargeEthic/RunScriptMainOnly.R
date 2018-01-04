##%%%%%%%%%%%%%%%%%%% ##
## From the Args File ##
## %%%%%%%%%%%%%%%%%% ##

load("N:/Simulations/DataSet/AllDataOut.RData")

MCMCSeed = 1029394

## %%%%%%%%%%%%%%%%%%% ##
##  Set the Directory  ##
## %%%%%%%%%%%%%%%%%%% ##

setwd("N:/Simulations")

## %%%%%%%%%%%%%%%%%%% ##
##    Load the Data    ##
## %%%%%%%%%%%%%%%%%%% ##

#load("N:/Simulations/DataSet/DataFormatOut.RData")

#File2012[,"ethnic"]=File2011[,"ethnic"]

## %%%%%%%%%%%%%%%%%%% ##
##  Read User Inputs   ##
## %%%%%%%%%%%%%%%%%%% ## 

source("N:/Simulations/SimulationScripts/Inputs.R") 

# Change the priors! 
#a.beta.prior[,6]  = 0
#b.beta.prior[,6]  = 0 

#a.beta.prior[2,5] = 1300 
#b.beta.prior[2,5] = 7500

## %%%%%%%%%%%%%%%%%%% ## 
##     Scripts         ## 
## %%%%%%%%%%%%%%%%%%% ## 

# Loading the DP package

myrepo = getOption("repos")
myrepo["CRAN"] = "http://archive.linux.duke.edu/cran/"
options(repos=myrepo)
install.packages("NPBayesImpute")
library(NPBayesImpute)
rm(myrepo)

require(NPBayesImpute)
#require(parallel)

## %%%%%%%%%%%%%%%%%%%%%% ## 
## Create the Data/ Store ##
## %%%%%%%%%%%%%%%%%%%%%% ## 

source("N:/Simulations/SimulationScripts/SansParallelFunctions.R")
source("N:/Simulations/SimulationScripts/FunctionsDalzell.R")
source("N:/Simulations/SimulationScripts/FunctionsGOnly.R")
source("N:/Simulations/SimulationScripts/Functions_Needed.R")
source('N:/downloads/run_rlMCMC_fast.R')
source('N:/Simulations/SimulationScripts/SansParallelFunctionsNew.R')

propose.B <-function(J, d, File1_dataBlock, File2_dataBlock, Estar, Move.F1, Move.F2, current.B, Bhat1, Bhat2, labelsJ_E, Index_to_Block, Legal_Index,phi,z,File2_dataLambda, File1_dataLambda, gamma.samp,A.probs,BlockSeed_restriction,SeedsOnly,possible1,possible2){
  
  # Create an empty list 
  proposal.info <-list() 
  
  # Initiate the Length of the list 
  H = 0 
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
  #     Begin with File 1 
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
  
  Epart  = Estar[[1]]
  f      = 1 
  Lambda = File1_dataLambda
  B.part = current.B[[1]]
  
  for( i in Move.F1 ){
    
    lam = Lambda[i]
    # Build a Space to store the proposal 
    B.holder = B.part[i,]
    
    # Create the Proposal
    
    #Which field are we proposing to change? 
    which.changed  = Epart[i,"jstar"]
    # What is the original observed field value? 
    observed.field = Bhat1[i,which.changed]
    # What is the current field value? 
    current.field  = as.integer( B.holder[which.changed] )
    
    # The multinomial probabilities associated with this latent class, this field 
    initial.probs  = phi[1:d[which.changed],z[lam],which.changed] 
    
    together = c(current.field,observed.field)
    
    # What options do we have for moving? 
    field.options  = 1:d[which.changed]
    field.options  = field.options[ - together ] 
    num.fields     = length( field.options )
    
    if( Epart[i, which.changed]==0 ){
      
      # Option 3: We move from an error to observed 
      
      if( observed.field == current.field){
        print(i)
        print("WHOA! PROBLEM!")
      }
      
      B.holder[which.changed] = observed.field 
      
      # For the transition probability ratio, need prob moving to current from observed
      
      norm.probs  = initial.probs[-c(observed.field)] 
      norm.probs  = norm.probs/ sum(norm.probs) 
      desired     = which( c(1:d[which.changed])[-observed.field] == current.field)
      result.prob = norm.probs[desired] 
      
      prior.part = initial.probs[observed.field]/ initial.probs[current.field]
      
    } else{ 
      
      # Options 1 and 2 : We move to an error prone field       
      
      if(  num.fields == 1){
        
        result     = field.options
        #B.holder[which.changed] = result
        
        prior.part = initial.probs[result]/ initial.probs[current.field]
        
        result = 1 
        norm.probs = 1 
        
      } else{
        
        norm.probs = initial.probs[-together] 
        norm.probs = norm.probs/ sum(norm.probs) 
        result     = sample( num.fields,size=1,replace=T, prob= norm.probs )
        B.holder[which.changed]  = field.options[result]
        
        prior.part = initial.probs[field.options[result]]/ initial.probs[current.field]
        
      }      
      
    } 
    
    # Which block index are we proposing to move to? 
    to.block             = category.block(B.holder,J,d)
    
    # If this index belongs to a filled block, 
    if( to.block %in% Legal_Index & ObsBlock_restriction == "yes"){
      to.block            = which(Index_to_Block==to.block)
      if(to.block %in% SeedsOnly & BlockSeed_restriction=="yes"){
        stop("Whoopse! We propose a seed move")
        #break 
      }
      #to.block           = Index_to_Block[2,to.block] 
      
      from.block          = File1_dataBlock[i]
      
      if(  to.block != from.block ){  
        H = H + 1
        
        if( Epart[i,"Aratio"] ==1){
          
          A.ratio = prior.part 
          A.ratio = A.ratio * ( norm.probs[result] )^(-1)
          
          # For transition prob ratio, need prob moving to current from proposed
          B.changed = unlist(B.holder[which.changed])
          norm.probs  = initial.probs[-c(observed.field,B.changed)] 
          norm.probs  = norm.probs/ sum(norm.probs) 
          desired     = which( c(1:d[which.changed])[-c(observed.field,B.changed)] == current.field)
          result.prob = norm.probs[desired] 
          
          A.ratio = A.ratio * result.prob
          #cat("The A ratio is", A.ratio ,"\n")
          
        } else if (Epart[i,"Aratio"]==2){
          
          A.ratio = 1/num.fields
          A.ratio = A.ratio*( norm.probs[result] )^(-1)
          A.ratio = A.ratio* prior.part 
          #cat("The A ratio is", A.ratio ,"\n")
          
        } else if (Epart[i,"Aratio"]==3){
          
          A.ratio = num.fields
          A.ratio = A.ratio*result.prob
          A.ratio = A.ratio*prior.part
          #cat("The A ratio is", A.ratio ,"\n")
        }
        
        proposed.field      = B.holder[which.changed]
        
        propose.move        = matrix( c(f,i,from.block,to.block, A.ratio,proposed.field,which.changed),ncol=7, nrow=1)
        
        proposal.info[[H]]  = propose.move
      }
    }
  }
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
  #     And Now File 2 
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
  
  Epart   = Estar[[2]]
  f       = 2 
  Lambda  = File2_dataLambda
  B.part  = current.B[[2]]
  
  for( r in Move.F2 ){ 
    
    lam            = Lambda[r]
    # Build a Space to store the proposal 
    B.holder       = B.part[r,]
    
    # Create the Proposal
    which.changed  = Epart[r,"jstar"]
    observed.field = as.integer(Bhat2[r,which.changed])
    current.field  = as.integer( B.holder[which.changed] )
    #together       = c(current.field,observed.field)
    initial.probs  = phi[1:d[which.changed],z[lam],which.changed]
    
    # Which values are possible for the record? 
    #if( which.changed == 4){
    which.possible = possible1[[r]]
    #} else{ 
    #  which.possible = possible2[[r]]
    #}
    field.options  = 1:d[which.changed]
    field.options  = field.options[ - observed.field ] 
    field.options  = intersect(which.possible,field.options)
    num.fields     = length( field.options )
    
    if( num.fields >= 1){
      # if there is at least one legal option (which there should be in less we are dealing with a type 2 seed)
      if( Epart[r, which.changed]==0 ){
        
        # Option 3: We move from an error to observed 
        B.holder[which.changed] = observed.field 
        
        # For the transition probability ratio, need prob moving to current from observed
        
        norm.probs  = initial.probs[field.options] 
        norm.probs  = norm.probs/ sum(norm.probs) 
        desired = which( which.possible == observed.field) 
        desired = which ( which.possible[-desired] == current.field)
        result.prob = norm.probs[desired] 
        
        prior.part = initial.probs[observed.field]/ initial.probs[current.field]
        
      } else{ 
        
        # Options 1 and 2 : We move to an error prone field       
        
        if(  num.fields == 1){
          # There is exactly one option
          result     = field.options
          B.holder[which.changed] = result
          
          prior.part = initial.probs[result]/ initial.probs[current.field]
          
          result = 1 
          norm.probs = 1 
          
        } else{
          
          norm.probs = initial.probs[field.options] 
          norm.probs = norm.probs/ sum(norm.probs) 
          result     = sample( num.fields,size=1,replace=T, prob= norm.probs )
          B.holder[which.changed]  = field.options[result]
          
          prior.part = initial.probs[field.options[result]]/ initial.probs[current.field]
          
        }      
        
      } 
      
      # Which block index are we proposing to move to? 
      to.block             = category.block(B.holder,d,J)   
      
      if( (to.block%in% Legal_Index)==FALSE){
        stop( paste("Record", r, "was proposed to be moved to an illegal block. The record is in block", c(B.holder), "and is proposed to be moved to block ",to.block,"\n"))
      }
      
      # If this index belongs to a filled block, 
      if( ObsBlock_restriction == "no" | to.block %in% Legal_Index){
        to.block            = which(Index_to_Block==to.block)
        
        if(to.block %in% SeedsOnly & BlockSeed_restriction=="yes"){
          stop("In propose.B File 2, we propose a move to a seeds only block, even though the Block Seed Restriction is on.")
        }
        
        from.block          = File2_dataBlock[r]
        
        # if( from.block == to.block){
        #   stop( paste("Record",r,"was proposed to move to the block it is already in."))
        # }
        
        if(  to.block != from.block ){  
          H = H + 1
          
          if( Epart[r,"Aratio"] ==1){
            # E(s) = E* = 1 
            A.ratio = prior.part 
            A.ratio = A.ratio * ( norm.probs[result] )^(-1)
            
            # For transition prob ratio, need prob moving to current from proposed
            #B.changed = unlist(B.holder[which.changed])
            norm.probs  = initial.probs[-observed.field] 
            norm.probs  = norm.probs/ sum(norm.probs) 
            desired     = which( c(1:d[which.changed])[-c(observed.field)] == current.field)
            result.prob = norm.probs[desired] 
            
            A.ratio = A.ratio * result.prob
            #cat("The A ratio is", A.ratio ,"\n")
            
          } else if (Epart[r,"Aratio"]==2){
            # E(s) = 0, E* = 1 
            A.ratio = 1/num.fields
            A.ratio = A.ratio*( norm.probs[result] )^(-1)
            A.ratio = A.ratio* prior.part 
            #cat("The A ratio is", A.ratio ,"\n")
            
          } else if (Epart[r,"Aratio"]==3){
            
            A.ratio = num.fields
            A.ratio = A.ratio*result.prob
            A.ratio = A.ratio*prior.part
            #cat("The A ratio is", A.ratio ,"\n")
          }    
          
          proposed.field      = B.holder[which.changed]
          
          propose.move        = matrix( c(f,r,from.block,to.block, A.ratio,proposed.field,which.changed),ncol= 7, nrow=1)
          
          #print(class(propose.move))
          
          proposal.info[[H]]  = propose.move 
        }
      }
    }
  }
  
  newlist <-list( "no.proposed" = H, "move.holder"= proposal.info )
  
  return( newlist ) 
  
}

## %%%%%%%%%%%%%%%%%%%%%% ## 
## Create the Data/ Store ##
## %%%%%%%%%%%%%%%%%%%%%% ## 

#source("/dscrhome/nmd16/NewSims/DataGen.R")

# Define the name
namepart = "School"

type1seeds = type1Seeds
rm(type1Seeds)

# Set/Create the Working Directory
dir.create(file.path(mainDir=getwd(), subDir=namepart))
setwd(file.path(mainDir=getwd(), subDir=namepart) )

save.image("StartSpace.RData")


# %%%%%%%%%%%%%%%%%%%%%%%%% #  
# What are the type 1 seeds?  
# %%%%%%%%%%%%%%%%%%%%%%%%% #

matches <-matrix(NA,nrow = length(type1seeds), ncol = 2)
matches[,2] = matches[,1] = type1seeds

## %%%%%%%%%%%%%%%%%%%%% # 
##     Pointers         ##
## %%%%%%%%%%%%%%%%%%%%% # 

Y1.name = "mathscal"
Y2.name = "Cmathscal"

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
##     Format the Data for the Run        ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
##    Run the script : Gutman Test   ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

File2011 = Data[,c(11,1:3,9,8,23)]
File2012 = Data[,c(22,12:14,20,19,24)]
colnames(File2011) = colnames(File2012) = c("mathscal","birthday","birthmonth","birthyear","sex","ethnic","sch")

File2012[,1]=File2012[,1]-mean(File2012[,1])

# What are we running? Gutman Only (GOnly), a test run (GTest), or the full model? 
method = "GOnly"

run_rlMCMC_fast(MCMCSeed="test", File2011T,File2012T, method,its,burnin,DPits, DPburnin, gap, reps, thinning, DPthinning, printgap, Comment, coefs.primary, coefs.primary.all, coefs.secondary, col.Y1, col.Y2, secondary.option, secondary.option2, ObsBlock_restriction="yes", BlockSeed_restriction="no", error.file, error.fields, blockingVars, matches, a.beta.prior, b.beta.prior, Hstar, DPOffset, Y1.name, Y2.name,paste(namepart,"Test",sep=""),threshold,complete.out,type2seeds= NULL)


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
##     Introduce The Errors           ## 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##   

# Error.Options =  setdiff( 1:N, type1seeds)
# error.perc = 18
# in.error       = sample( Error.Options, N*error.perc/100, replace = F)
# rm(Error.Options)
# 
# field.options = levels(File2011[,error.fields])
# 
# fields.probs = table(File2011[,"ethnic"])/N
# 
# for( i in in.error){
#   
#   if( error.file == 1 ){
#     current.field   = File1[ i ,error.fields]
#     field.choices   = setdiff( field.options,current.field )
#     File1[i,"prog"] = sample(field.choices,1)
#   } else{ 
#     current.field   = File2[ i ,error.fields]
#     field.choices   = setdiff( field.options,current.field )
#     probs.here = field.probs[-current.field]
#     File2[i,"prog"] = sample(field.choices,1,prob)
#   }
#     
# }
# 
# rm(i,current.field,field.choices,field.options)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
##    Format for Data Runs     ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

  realmin =  1e-08
  
  ## %%%%%%%%%%%%%%%%% ## 
  ##  Check the Inputs ##
  ## %%%%%%%%%%%%%%%%% ##

File1 = File2011T
File2 = File2012T
  
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
  
MCMCSeed = 9191020
  if( class(MCMCSeed) != "numeric"){
    MCMCSeed <-as.numeric(Sys.time())
    MCMCSeed = (MCMCSeed - floor(MCMCSeed))*1e8
    warning(paste("No random seed specified, seed of", MCMCSeed, "utilized."))
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

set.seed(MCMCSeed)
  
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
  
  # Remember that to use the DP code from Quanli, we need Bhat to be a data frame with each column a factor.
  
  Bhat1          = File1[,where.bvars.F1]
  Bhat2          = File2[,where.bvars.F2]
  #Bhat1[,4]      = factor(Bhat1[,4])
  #Bhat2[,4]      = factor(Bhat2[,4])
  Bhat           = rbind(Bhat1,Bhat2)
  rownames(Bhat) = 1:N
  
  B     <- Bhat 
  #B[,4] <- factor(B[,4])
  
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

  B[,6]= as.numeric(as.factor(B[,6]))
  B[,6]=factor(B[,6])
  
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

  coefs.primary = c("sex","ethnic") 
  coefs.secondary = c("sex","ethnic")

  coefs.all = union(coefs.primary,coefs.secondary)
  modelPart <- model.matrix( lm( File2011[,1] ~ Cmathscal + ethnic+sex, data = File2011))
  modelPartSecondary <-  model.matrix( lm(formula = File2011[, 1] ~ sex+ethnic, data = File2012))

  setdiff(colnames(modelPartSecondary),colnames(modelPart))
  
  File1_data = cbind( File1[,col.Y1], modelPart[,-c(1:2)])
  colnames(File1_data)[col.Y1] = Y1.name

modelPart <- model.matrix( lm( File2011[,1] ~ Cmathscal + ethnic + sex, data = File2012))


  File2_data = cbind( File2[,col.Y2], modelPart[,-c(1:2)])
  colnames(File2_data)[col.Y2] = Y2.name

  coefs.all = union(coefs.primary.all,coefs.all)

  coefs.all =colnames(File2_data)
  
   where.coefs.primary = 1:4
  
   where.coefs.secondary = c(1:4)
  
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
  
  #if( method =="Extension"){
      
    ## Check for Type 2 Seeds
    #possibleSchool <- how_many_possible( B, Nhere = N/2, d, field = 6,J, type1seeds, Legal_Index,blocks)

how_many_possible <- function( B.new, Nhere, d, field, J, type1seeds, Legal_Index,blocks){
  
  # Checking for School possible block moves  
  SchoolLevel = B.new[,field]
  SchoolLevel = as.numeric(as.factor(SchoolLevel))
  
  # This part is fixed and fed in 
  check.gap.1   = category.block( c(1,1,1,1,1,1), d, J)
  second        = c(1,1,1,1,1,1)
  second[field] = 2
  check.gap.2   = category.block( second, d, J)
  check.gap     = check.gap.2 - check.gap.1
  school.levels = d[field]
  add.vec       = cumsum( rep(check.gap,school.levels-1))
  
  possible.out = vector( "list", length = Nhere )
  
  a = 1 
  
  for( i in 1:Nhere){
    
    if( i %in% type1seeds){
      possible.out[[a]] = 0 
    } else{ 
      k = i + Nhere
      # which index are we in? 
      x = blocks[k]
      # which level of school? 
      s = SchoolLevel[k] 
      
      up.move   = school.levels - s
      down.move = s -1 
      
      if(up.move > 0 ){
        up.check  = x + add.vec[1:up.move]
      } else{
        up.check = NULL
      }
      if(down.move > 0){
        down.check  = x - add.vec[1:down.move]
      }else{
        down.check = NULL 
      }
      
      possible.out[[a]] = s
      # Which of the up checks are legal?
      possible.up <- which( up.check %in% Legal_Index)
      if( length(possible.up) > 0 ){
        possible.out[[a]] = c(possible.out[[a]],s+possible.up)
        #which.up = s+ possible.up  
      }
      # Which of the down checks are legal?
      possible.down <- which( down.check %in% Legal_Index)
      if( length(possible.down) > 0 ){
        possible.out[[a]] = c(possible.out[[a]],s-possible.down)
        possible.out[[a]] = sort(possible.out[[a]])
      }
      
      # Now, we check to see which of these are legal indices
      
      possible.moves = c(up.check, down.check)
      possible.moves = intersect(possible.moves,Legal_Index)
      #possible.moves = length(possible.moves)
      # We need what schools are associated with these moves. 
      
    } 
    
    a = a + 1 
    
  }
  
  return(possible.out)
  
  
}
    
   error.fields = "ethnic"
   possibleEthnic <- how_many_possible( B, Nhere = N/2 , d, field = which(blockingVars == error.fields), J, type1seeds, Legal_Index,blocks)
    
    ## Check for Type 2 Seeds
    #possible.length.S = unlist( lapply( possibleSchool, length ))
    #seeds.school = which( possible.length.S ==1)
    possible.length.E = unlist( lapply( possibleEthnic, length ))
    T2check = lapply(possibleEthnic,function(x) 0 %in% x)
    type2seeds = which(T2check=="TRUE")
    #type2seeds = which(possible.length.S==1 & possibleSchool==0)
    #seeds.ethnic = which( possible.length.E ==1)
    #seeds.all = intersect( seeds.school, seeds.ethnic)
    type2seeds = setdiff( type2seeds, type1seeds)
    #seeds.ethnic = setdiff(seeds.ethnic,c(seeds.school,type1seeds))
    #seeds.school = setdiff(seeds.school,c(seeds.ethnic,type1seeds))
    # Any records for which possibleSchool = 0 is going to be declared a type 2 seeds
    #possible.length = unlist( lapply( possibleSchool, length ))
    #seeds.all = which( possible.length ==1)
    #type2seeds = setdiff( seeds.all, type1seeds)
  #}
  
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
  
  #possible <- how_many_possible( B, N = 10000, d, field = 6,J, type1seeds, Legal_Index,blocks)
  
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
  ## Save the formatted data sets  ## 
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
  
  save(File1_data,file ="File1Formatted.RData")
  save(File2_data,file ="File2Formatted.RData")
  save.image(paste(namepart,"PreMain.RData",sep=""))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
##    Run the script: Gutman   ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
  ####      Pre Run Steps       ####
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
  
  method = "GOnly"

  #namepart = "MainOnly"
  
  if( method == "GOnly"){ 
    
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
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
    #### Run the Gutman MCMC Only ####
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
    
    nonseed.F2 = which(File2_data[,"Seed"]==0)
    n.nonseed.F2 = length(nonseed.F2)
    
    # Run the script 
    gutman_mcmc_Speed("no",MCMCSeed, its,burnin, thinning, gap, n.RIB, RIB, File1_data, File2_data,K,which.Y2, which.Y1,type1seeds, where.coefs.primary, where.coefs.secondary, p, p2, n.type1,needs.C,reps,PermSampling,secondary.option="dependent", n1, n2,col.Y1, col.Y2,namepart)
    print("we completed the Gutman Run")
    
  }

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
##    Run the script: Extension  ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

  save.image(paste("PostGutman",repindex,".RData",sep=""))

  if( class(MCMCSeed) != "integer"){
    MCMCSeed <-as.numeric(Sys.time())
    MCMCSeed = (MCMCSeed - floor(MCMCSeed))*1e8
    warning(paste("No random seed specified, seed of", MCMCSeed, "utilized."))
  }

  Bhat1[,6]=  as.numeric(as.factor(B[1:(N/2),6]))
  Bhat2[,6]=  as.numeric(as.factor(B[(N/2+1):N,6]))

  Bhat1[,6]= B[1:n1,6]
  Bhat2[,6]= B[(n1+1):N,6]
nonseed.F2 = setdiff(1:(N/2), type1seeds)
  nonseed.F2 = nonseed.F2[ which(nonseed.F2 <= n2)]
  run_Extension()

  save.image(paste("PostDalzell",repindex,".RData",sep=""))

namepart = "MainOnlyRun2"







