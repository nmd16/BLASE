# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
####            Functions                 ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

get_starting_permutation<-function(K, RIB, matches,File1_data, File2_data,n1,n2 ){
  
  holder1 = rep( 0, n1 ) 
  holder2 = rep( 0, n2 )
  lambdaholder = rep( 0, n2 )
  
  for(k in 1:K){
    
    #Step 1: Look at the records in the block
    in.block   = list( RIB[[1]][[k]], RIB[[2]][[k]] )
    max.length = max( length(in.block[[1]]),length(in.block[[2]]) )
    seeded  = list( intersect(in.block[[1]],matches[,1] ), intersect(in.block[[2]], matches[,1] ) )
    
    # If the length is not 0, proceed. This means that we have some elements in the block.
    if( max.length >0 ){
      
      # Begin with the seeded 
      b                = 0
      if( length(seeded[[1]]) > 0 ){
        
        # We begin by looking at the matched data. 
        
        for( s in c( seeded[[1]])){
          b = b + 1
          #matchup            = which(matches[,1]==s) #Which row in the matches table?
          #matchup            = matches[ matchup,1 ]
          holder1[s]         = b 
          holder2[s]   = b 
          lambdaholder[ s ]  = File1_data[s,"lambda"]
        }
        
        in.block[[1]]        = setdiff( c(in.block[[1]]),c(seeded[[1]]))
        in.block[[2]]        = setdiff( c(in.block[[2]]),c(seeded[[2]])) 
        max.length           = max( length(in.block[[1]]),length(in.block[[2]]) )        
        
      }   
      
      if( max.length >0 ){
        # Fill in the information from file 1
        places1 = 1:length(in.block[[1]])
        places1 = places1 + b 
        holder1[in.block[[1]]] = places1 
        #Fill in the information from file 2
        places2 = permute.places( 1:length(in.block[[2]]) )
        #if(k ==16){
        #  print( in.block)
        #}
        places2 = places2 + b
        
        holder2[in.block[[2]]] = places2
        
        # Fill in the value for lambda2
        places2 = as.numeric(as.factor(places2))
        if( length(places2) > 1 ){ 
          lambdaholder[ in.block[[2]] ] = c(in.block[[1]])[ order(places1) ][ c(places2) ] 
          #lambda_F1[ order(BlockRow_F1) ][result]
        } else{
          lambdaholder[ in.block[[2]] ] = in.block[[1]]
        }
      }
      
    }
    
    if(length(unique(lambdaholder[ which(lambdaholder> 0) ])) != length( lambdaholder[ which(lambdaholder> 0)]) ) {
      print(k)
      stop( paste("In block", k, "we have assigned at least one individual label to multiple records."))
    }
  }
  
  File1_data[,"BlockRow"] = holder1 
  File2_data[,"BlockRow"]= holder2
  File2_data[,"lambda"]   = lambdaholder
  
  newlist <- list( "File1" = File1_data, "File2" = File2_data )
  
  return( newlist )
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
####          Main Function               ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# Purpose: 
# Runs the new model MCMC

# Saves output as .RData files, exports matched data files as .csv 

# Returns: 
# out.list - a list which contains the regression parameters samples from the MCMC
new_mcmc <- function(its, burnin,thinning,gap,seed,namepart,File1_data, File2_data,is.test,printgap, Legal_Index,Y2.name,n1,n2,p, p2, secondary.option, J,nonseed.F1,nonseed.F2,n.error1,n.error2,Bhat1,Bhat2,type1Seeds, PermSampling, SeedsOnly, RIB, n.RIB,Hstar,where.coefs.primary,where.coefs.secondary,where.coefs, Comment,Imp.Holder1,Imp.Holder2, J_E, K, realmin, reps, col.Y1, col.Y2, count.matches,d,labelsJ_E,n.type1,threshold,which.Y1,which.Y2,Index_to_Block, a.beta.prior, b.beta.prior,secondary.option2 = "dependent",error.files,DPoffset,dpmpm,model,Bhat,where.bvars.F1,where.bvars.F2,DPits,DPburnin,DPthinning){
  
  # %%%%%%%%%%%%%%%%%%%%%%%%% # 
  ##          Inputs         ##
  ## %%%%%%%%%%%%%%%%%%%%%%%  #
  
  # its         = how many iterations
  # burnin      = number of discarded iterations
  # thinning    = gap between saved post-burnin iterations
  # gap         = gap between saving .csv files for the matched data 
  # seed        = the random seed for the MCMC
  # namepart    = the name attached to this simulation for saving purposes 
  # File1_data  = The File 1 information; see CodeKey.R for column breakdown
  # File2_data  = The File 2 information
  # is.test     = do we fix the regression parameters in advance? 
  # printgap    = after how many iterations should we print out the iteration as a progress check
  # Legal_Index = which indices are we allowed to propose to move to
  # Y2.name     = what is the name of the continuous feature from File 2? 
  # n1          = the original number of records in File 1
  # n2          = the original number of records in File 2
  # p           = the number of parameters in the primary regression
  # p2          = the number of parameters in the secondary regression
  # secondary.option = "independent" if the secondary regression does not depend on the blocking variables, = "dependent" otherwise 
  # J           = the number of blocking variables
  # nonseed.F2  = a vector containg the rownumbers for all File 2 records which are not seeds
  # n.error1    = how many records can possibly contain errors in File1?
  # n.error2    = how many records can possibly contain errors in File 2? 
  # Bhat1       = the recorded blocking variables for File 1 (n1 by J matrix)
  # Bhat2       = the recorded blocking variables for File 2 (n2 by J matrix)
  # type1Seeds  = a vector containing the row labels for all type 1 seeds
  # maxDPcount  = the number of stored DP output lists 
  # dpName      = the name attached to these output lists 
  # PermSampling = vector takes a 1 if we use exact sampling, 2 if we use MH and 0 otherwise.
  # SeedsOnly   = vector containing all blocks that initially have only seeds
  # RIB         = list of RIB[[1]] and RIB[[2]], each lists of K list telling which records from file f are in each block
  # n.RIB       = 2 by K matrix where n.RIB[f,k] counts how many records from file f are in block k
  # Hstar       = the number of latent components in the DP 
  # where.coefs.primary = the columns in File2_data which hold information for the primary regression
  # where.coefs.secondary = the columns in File2_data which hold information for the secondary regression
  # where.coefs = a vector containing the columns involved in both the primary and secondary regressions
  
  # Added!!! 
  
  # Comment     = "YES" if we want a ton of information to print out
  # Imp.Holder1 = a list which stores skeleton imputations for a record from each block for File 1 
  # Imp.Holder2 = a list which stores skeleton imputations for a record from each block for File 2 
  # J_E         = a vector which tells us which column of the B matrix are allowed to be in error
  # K           = the number of possible blocks  
  # realmin     = 1e-08
  # reps        = the number of times in the Guthan MCMC that we propose to switch the position of two records
  # col.Y1      = pointer to which column in File1_data has the continuous response
  # col.Y2      = pointer to which column in File2_data has the continuous predictor 
  # count.matches = count of how many type 1 seeds are in each block 
  # d           = vector which tells us how many levels there are for each blocking variable
  # labelsJ_E   = the names of the columns in B which are allowed to be in error 
  # n.type1     = how many type 1 seeds? 
  # threshold   = if the number of records is at least the threshold, we cannot use exact sampling 
  # which.Y1    = 
  # which.Y2    = which column in the design matrix for the primary regression pertains to Y2 
  # Index_to_Block = a 2 x K matrix, in which the first row gives the index and the second gives the corresponding block label. 
  # a.beta.prior   = the hyperparamter for gamma
  # b.beta.prior   = hyperparamter for gamma 
  # nonseed.F1     = vector of which records in File 1 are not type 1 seeds 
  
  
  # %%%%%%%%%%%%%%%%%%%%%%%%% # 
  ##          Script         ##
  ## %%%%%%%%%%%%%%%%%%%%%%%  #
  
  # Set the Seed 
  set.seed(seed)
  
  # Identify the number of records in each block which are imputations 
  Imputed.Block = matrix( 0, nrow = K, ncol = 2) 
  # In row 1, column 1, we see the number of imputed records in file 1 for that block. 
  Imputed.Records = replicate(K, list())
  # Each row is a block, each column is a File
  
  for( k in 1:K){
    # File 1
    tester = (RIB[[1]][[k]] > n1)
    sum.tester = sum(tester)
    
    Imputed.Block[k,1] = sum.tester 
    
    if( sum.tester > 0 ){
      Imputed.Records[[k]][[1]] = which( RIB[[1]][[k]] > n1)
    }
    # File 2 
    tester = (RIB[[2]][[k]] > n2)
    sum.tester = sum(tester)
    
    Imputed.Block[k,2] = sum.tester 
    
    if( sum.tester > 0 ){
      Imputed.Records[[k]][[2]] = which( RIB[[2]][[k]] > n2)
    }
    
  }  
  
  ##  %%%%%%%%%%%%%%%%%%%  ##
  #### Initialize Gamma  ####
  ### %%%%%%%%%%%%%%%%%%  ###  
  
  gamma.samp = matrix( 0, nrow = 2, ncol = J)
  
  for( u in labelsJ_E){
    gamma.samp[1,u] = rbeta( 1, a.beta.prior[1,u], b.beta.prior[1,u] )
    gamma.samp[2,u] = rbeta( 1, a.beta.prior[2,u], b.beta.prior[2,u] )
  }
  
  ### %%%%%%%%%%%%%%%%%%% ##
  #### Initialize E, B, Bhat   ####
  ### %%%%%%%%%%%%%%%%%%% ##
  
  E1     = matrix(0,nrow=n1, ncol=J)
  jstar  = rep(0,n1)
  Aratio = rep(0,n1)
  E1     = cbind( E1, jstar,Aratio)
  
  E2 = matrix(0,nrow=n2, ncol=J)
  jstar = rep(0,n2)
  Aratio = rep(0,n2)
  E2 = cbind( E2, jstar,Aratio)
  
  current.B       = list( "File1" = Bhat1, "File2" = Bhat2)
  Bhat            = rbind(Bhat1,Bhat2)
  
  ### %%%%%%%%%%%%%%%%%%% ##
  ###    Utilities      ###
  ### %%%%%%%%%%%%%%%%%%% ##
  
  # Initialize an empty list 
  empty.blocks      <-c() 
  
  # How many rows in File 1? in File 2? 
  nrow.F1 = nrow(File1_data)
  nrow.F2 = nrow(File2_data)
  
  # Which records from File2 were imputed?
  if( nrow.F2 > n2){
    imputed.F2 = (n2+1):nrow.F2
    n.imputed.F2 = length(imputed.F2)
  } else{     
    imputed.F2   = NULL
    n.imputed.F2 = 0   
  }
  
  # Which records from File1 were imputed?
  if( nrow.F1 > n1){
    imputed.F1 = (n1+1):nrow.F1
    n.imputed.F1 = length(imputed.F1)
  } else{     
    imputed.F1   = NULL
    n.imputed.F1 = 0   
  }
  
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
  #### Adaptation Set Up: Data Storage  ###########
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##  
  
  theta.out              <- NULL
  reg2.out               <- NULL
  E.out                  <- NULL
  gamma.out              <- NULL
  empty.blocks.store     <-list()
  blocks.out             <- NULL
  lambda.out             <- NULL 
  move.counter.store     <- c()
  Y2.out                 <-c()
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ###########  Initialization  Theta #############
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  # We begin by using ONLY the seeds information 
  
  if( class(type1Seeds) != "NULL"){
    
    Type1SeedIndicator = 1
  
    Y1 = File1_data[type1Seeds,col.Y1]
    Y2 = File2_data[type1Seeds,col.Y2]
  
    # Create the model matrix based on the seeds
    Vtheta = cbind( 1, File2_data[type1Seeds, where.coefs.primary ] )
        
    Sigma1.Part = (n.type1 - p )/2
    
  } else{ 
    Type1SeedIndicator = 0
    
    Y1 = File1_data[,col.Y1]
    Y2 = File2_data[,col.Y2]
  
    # Create the model matrix based on the seeds
    Vtheta = cbind( 1, File2_data[, where.coefs.primary ] )
    
    Sigma1.Part = (nrow.F1 - p )/2
  }
  
  Vtheta       = as.matrix(Vtheta)
  
  VthetaTVtheta= crossprod(Vtheta)
  Gen          = zapsmall( eigen( VthetaTVtheta )$values )
  Gen          = ifelse( sum(Gen==0) > 0 , "TRUE","FALSE")
  invert.part  = invert_XTX( VthetaTVtheta , Gen )
  
  t.Vtheta     = t(Vtheta)
  Proj         = invert.part%*%t.Vtheta
  hat.beta     = Proj%*%Y1
  
  ResidSS      = (Y1 - Vtheta%*%hat.beta) 
  ResidSS      = t(ResidSS)%*%ResidSS
  ResidSS      = .5* ResidSS  

  sigmasq      = rigamma(1,Sigma1.Part,ResidSS) 
  sigma        = sqrt(sigmasq)
  var.beta     = sigmasq*invert.part
  beta.star    = mvrnorm( 1, hat.beta, var.beta ) 
  cat( "The inital beta draw is",beta.star,"\n")
  
  rm( Sigma1.Part,var.beta,ResidSS,hat.beta, Proj, t.Vtheta,Gen,invert.part,VthetaTVtheta)
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ##########  Initialization : Y2 Coefs  ##########
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  if( secondary.option == "independent"){
    
    if( class(type1Seeds) != "NULL"){
      Sigma2.Part = (n.type1 - p2 )/2
      barY2 = mean(Y2)
      Var.Part = sum( (Y2 - barY2)^2 )
    
      sigma2sq = rigamma(1,Sigma2.Part,.5*Var.Part)
      eta      = rnorm( 1, barY2, (1/n.type1)*sigma2sq)
    
      rm(Var.Part,barY2,Sigma2.Part)
    } else{
      Sigma2.Part = (n2 - p2 )/2
      barY2 = mean(Y2)
      Var.Part = sum( (Y2 - barY2)^2 )
    
      sigma2sq = rigamma(1,Sigma2.Part,.5*Var.Part)
      eta      = rnorm( 1, barY2, (1/n2)*sigma2sq)
    
      rm(Var.Part,barY2,Sigma2.Part)
    }
    
  } else{
    
    # The Y2s are dependent on the blocking variables 
    
    # Create the model matrix based on the seeds
    if( class(type1Seeds) != "NULL"){
       V2 = cbind( 1, File2_data[type1Seeds, where.coefs.secondary ] )
      Sigma2.Part = (n.type1 - p2 )/2
    } else{
       V2 = cbind( 1, File2_data[, where.coefs.secondary ] )
       Sigma2.Part = (n2 - p2 )/2
    }
     
    V2 = as.matrix(V2)
    
    V2TV2        = crossprod(V2)
    Gen          = zapsmall( eigen( V2TV2 )$values )
    Gen          = ifelse( sum(Gen==0) > 0 , "TRUE","FALSE")
    invert.part  = invert_XTX( V2TV2 , Gen )
    
    t.V2         = t(V2)
    Proj         = invert.part%*%t.V2
    hat.eta      = Proj%*%Y2
    est.mean.Y2  = V2%*%hat.eta
    
    ResidSS      = (Y2 - est.mean.Y2) 
    ResidSS      = t(ResidSS)%*%ResidSS
    ResidSS      = .5* ResidSS    
  
    sigma2sq     = rigamma(1,Sigma2.Part,ResidSS) 
    sigma2       = sqrt(sigma2sq)
    var.eta      = sigmasq*invert.part
    eta          = mvrnorm( 1, hat.eta, var.eta ) 
    cat( "The inital eta draw is",eta,"\n")  
    
    
    rm( Sigma2.Part,var.eta,ResidSS,hat.eta, Proj, t.V2,Gen,invert.part,V2TV2,est.mean.Y2) 
    
  } 
  
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ####### Initialization : Y2 Imputations #########
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  theta.small       = beta.star[ which.Y2 ] 
  var.Y2impute      = ( sigma2sq^(-1) + sigmasq^(-1)*(theta.small^2) )^(-1)
  sqrt.var.Y2impute = sqrt(var.Y2impute)
  
  if( nrow.F2 > n2){
    
    # Which Y1 values are matched to imputations?
    F1.matched.to.imputations = File2_data[imputed.F2,"lambda"]
    
    # Store the Y1 associated with these values
    Y1.matched.to.imputations = File1_data[F1.matched.to.imputations,col.Y1]
    # Store the parts of the design matrix which go along with these values. 
    V1 = cbind( 1, File1_data[F1.matched.to.imputations, where.coefs.primary[-1] ] )
    V1 = as.matrix(V1)  
    
    beta.impute = beta.star[ -which.Y2 ] 
    
    mean.Y2impute = V1%*%as.matrix(beta.impute) - Y1.matched.to.imputations
    mean.Y2impute = sigmasq^(-1)*theta.small*mean.Y2impute
    
    if( secondary.option != "independent"){
      V2.imp = cbind( 1, File2_data[imputed.F2, where.coefs.secondary ] )%*%eta
      
      mean.Y2impute = sigma2sq^(-1)*V2.imp - mean.Y2impute
    } else{
      mean.Y2impute = sigma2sq^(-1)*eta - mean.Y2impute
    }
    
    mean.Y2impute = var.Y2impute*mean.Y2impute 
    imputed.all   = rnorm( n.imputed.F2, mean.Y2impute, sqrt.var.Y2impute )
    
    File2_data[ imputed.F2, col.Y2 ] = imputed.all  
    
  }
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ####### Initialization : Y1 Imputations #########
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  if( nrow.F1 > n1){
    
    # Reorganize the Vtheta based on the lambda ordering
    Vtheta = cbind( 1, File1_data[, where.coefs.primary] )
    Vtheta = Vtheta[order(File2_data[,"lambda"]),]
    Vtheta = as.matrix(Vtheta)
    Vtheta = Vtheta[imputed.F1,]
    
    imputed.all = rnorm( n.imputed.F1, Vtheta %*% as.matrix(beta.star), sqrt(sigmasq) ) 
    
    File1_data[ imputed.F1, col.Y1 ] = imputed.all
    
    
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
    ##########  Initialization : Deletions  ##########
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
    
    # Delete everything that we no longer need.   
    
    rm( Vtheta,imputed.all,mean.Y2impute, beta.impute,theta.small,V1,Y1.matched.to.imputations,F1.matched.to.imputations,Y1,Y2)
  }
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ##########  Pre Sampler: Constants     ##########
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  # These remain constant 
  HstarM1     = Hstar-1 
  hvec        = 1:Hstar
  
  Empty.Blocks <- setdiff(1:K,unique(c(File1_data[,"Block"],File2_data[,"Block"])))
  
  # Which records do NOT only contain seeds?
  NotOnlySeeds = setdiff(1:K,c(SeedsOnly))  
  
  max.RIB = apply( n.RIB, 2, max)
  max.RIB = max.RIB - count.matches
  non.seed.RIB = max.RIB
  
  Y1          = File1_data[,col.Y1]
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
  #################               Run New Sampler       ######################
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
  
  require(NPBayesImpute)
  
  for( j in 1:J){
    if(class(Bhat[,j]) != "factor"){
      Bhat[,j] = as.factor(Bhat[,j])
    }
  }
  
  model <-CreateModel( Bhat, NULL, Hstar, 0, 0.25,0.25)
  model$Run(DPburnin,DPits,DPthinning)
  dpmpm <-model$snapshot
  
  
  for( s in 1:its){ 

    
    if(s%%printgap == 0){ print(s) }
    
    ######################
    ### Propose a new E ##
    ######################
    
    Estar = propose.E(E1,E2,labelsJ_E,nonseed.F1,nonseed.F2,gamma.samp)  
    
    ## %%%%%%%%%%%%% ##
    ## DP Parameters ## 
    ## %%%%%%%%%%%%% ## 
    
    model$Run(0,2,2)
    dpmpm <-model$snapshot
    
    #################################################################################
    ####  Sample the Permutation Conditional on the Current Blocking Structure ######
    #################################################################################
    
    to.sampleC = which(PermSampling>0)
    to.sampleC = setdiff( to.sampleC, Empty.Blocks) # This SHOULD be superfluous, but is an okay check    
    
    BlockRow.F1 = File1_data[,"BlockRow"]
    BlockRow    = File2_data[,"BlockRow"]
    # Take the current ordering of the File 2 data 
    F1.Lambda   = 1:nrow.F1
    F2.Lambda   = File2_data[,"lambda"]
    Y2          = File2_data[,col.Y2]
    # Store a copy in the original ordering
    Y2.in       = Y2
    D.Holder    = File2_data[,where.coefs.primary]
    # Ordered
    V.Holder   = cbind(1, D.Holder[order(F2.Lambda),] )
    New.Ordering = F2.Lambda
    Sigma1.Part = (nrow.F1 - p )/2
    Sigma2.Part = (nrow.F2 - p2 )/2 
    
    # Sample a Permutation C for every block which requires it
    for( k in to.sampleC){
      # Which records from File 1 are in the block?
      F1.records = RIB[[1]][[k]] 
      # Which records from File 2 are in the block?
      F2.records = RIB[[2]][[k]]
      # What is the current ordering of the records within the block?
      order.in   = BlockRow[F2.records]
      # How many pairs are in the block?
      pairs      = n.RIB[1,k]
      # How many type 1 seeds?
      if( class(type1Seeds) != "NULL"){
        num.matched = count.matches[k]
      }else{
        num.matched = 0 
      }
      # Store the File 1 block information
      F1.holder  = Y1[F1.records]
      # Store the File 2 block information
      F2.holder  = cbind(1,D.Holder[F2.records,])
      # Sample a new permutation 
      block.order.out=A.sample.perm(F1.holder,F2.holder,pairs,order.in,beta.star,sigma,S.ind=PermSampling[k],reps,k,num.matched)     
      #Update the order of the File 2 Data
      New.Ordering[F2.records] = F1.Lambda[F1.records][block.order.out]
      #Update the Block Row 
      BlockRow[F2.records]    = block.order.out 
    } 
    
    # Reorder based on the sampled C 
    Y2                = Y2.in[order(New.Ordering)]
    V.Holder[,col.Y2] = Y2
    F2.Lambda = New.Ordering 
    
    # Store in the File 2 data space 
    File2_data[,"lambda"]   = New.Ordering
    File2_data[,"BlockRow"] = BlockRow 
    
    ###########################
    ##  Propose Block Moves  ##
    ###########################
    
    Move.F1 = which( Estar[[1]][,"Aratio"]>0)
    Move.F2 = which( Estar[[2]][,"Aratio"]>0)
    
    
    move.holder =  propose.B( J, d, File1_data[,"Block"], File2_data[,"Block"], Estar, Move.F1, Move.F2, current.B, Bhat1, Bhat2, Blocking_WithError, Index_to_Block, Legal_Index,dpmpm$psi,dpmpm$z+1, File2_data[,"lambda"], File1_data[,"lambda"],gamma.samp,A.probs,BlockSeed_restriction,SeedsOnly )
    
    ###############################
    ## Accept/Reject Block Moves ##
    ###############################
    
    nrow.F1              = nrow(File1_data) #Technically, this should already be stored, check! 
    nrow.F2              = nrow(File2_data) #Technically, this should already be stored, check! 
    rownames(File1_data) = 1:nrow.F1
    rownames(File2_data) = 1:nrow.F2
    
    beta.impute    = beta.star[-2]
    theta.small    = beta.star[2]
    
    Move.Holder = move.holder$move.holder
    
    Accept.Count = 0     
    
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%###
    ####     Begin the Loop     ####
    
    # Step 1: Grab the move holder    
    
    for( u in 1:length(Move.Holder)){
      
      needtobreak= 0 
      #File2_C = File2_data
      #File1_C = File1_data
      
      if( Comment == "YES"){
        print(u)
      }
      
      ###    %%%%%%%%%%%%%%%%%%%%%%%%%%%    ###
      ## Sanity Check Points: Delete Later   ##
      ###    %%%%%%%%%%%%%%%%%%%%%%%%%%%    ###    
      
      if( class(Empty.Blocks)=="NULL" ){
        choices = c(1:K)
      }else{
        choices = c(1:K)[-Empty.Blocks]
      }    
      
      for( r in choices){
        checkout = File1_data[RIB[[1]][[r]],"BlockRow"]
        if( length(unique(checkout)) != length(checkout)){
          needtobreak = 1 
          stop( paste("For record", r, "in File 1, more than one record is assigned to the same block row."))
        }
        if( max(unique(checkout)) != length(checkout)){
          needtobreak = 1 
          stop( paste("For record", r, "in File 1, more than one record is assigned to the same block row."))
        }
      } 
      
      
      for( r in choices){
        checkout = File2_data[RIB[[2]][[r]],"BlockRow"]
        if( length(unique(checkout)) != length(checkout)){
          needtobreak = 1 
          stop( paste("For record", r, "in File 2, more than one record is assigned to the same block row."))
        }
        if( max(unique(checkout)) != length(checkout)){
          needtobreak = 1 
          stop( paste("For record", r, "in File 2, more than one record is assigned to the same block row."))
        }
      } 
      
      if(needtobreak==1){
        stop('We have a break issue ')
      }
      
      ###    %%%%%%%%%%%%%%%%%%%%%%%%%%%    ###
      ## Store the information for this move ##
      ###    %%%%%%%%%%%%%%%%%%%%%%%%%%%    ###  
      
      holder.all   = Move.Holder[[u]]
      holder.all   = unlist(holder.all)
      holder.file  = holder.all[1]
      holder.r     = holder.all[2]
      holder.from  = holder.all[3]
      holder.to    = holder.all[4] 
      holder.A     = holder.all[5] 
      holder.jstar = holder.all[7]
      
      # The first element of this tells us which file we are working with (1 or 2)
      # The second tells us which record within that file
      # The third tells us the block that we were in
      # The fourth tells us the block that we are moving to
      # The fifth gives us part of the acceptance probabilities 
      # The sixth, which we haven't called here, tells us the proposed B field that caused the change. 
      
      AcceptOption = c(0,0)
      Deleted      = c()
      Added        <- list()
      exact        = c(0,0,0,0) #Default: Not exact for C*, but MH        
      
      #### Check for Empty Blocks   #### 
      
      EmptyTo   = holder.to %in% Empty.Blocks 
      EmptyFrom = 0      
      
      
      ####################################################################################################
      ################                      File 1 Record Moves           ################################
      ####################################################################################################  
      
      
      if( holder.file == 1 ){ 
        
        ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
        ####Step 2: Creating the To Block: File 1 Record Move : Is the block empty?####  
        ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
        
        if( EmptyTo == "TRUE" ){
          
          if( Comment == "YES"){
            print("EmptyTo")
          }
          
          # Store an indicator that the block is empty 
          exact[1] = 3
          exact[3] = 2
          
          # Grab the Y1 information 
          
          F1.recordsT     = holder.r
          F2.recordsT     = -99 
          
          y.to            = Y1[holder.r]        
          
          # Grab the shell of a record which defines the block. 
          to.addin        = Imp.Holder2[[holder.to]]
          
          ##############################
          # Impute a record for file 2 #
          ##############################
          
          V1.part = c(1, to.addin[where.coefs.primary[-1]])
          if( secondary.option == "independent"){
            V2.part = eta
          } else{ 
            V2.part = c(1, to.addin[where.coefs.secondary])
            V2.part = V2.part%*%eta 
          } 
          
          mean.Y2impute = y.to - V1.part%*%beta.impute 
          mean.Y2impute = sigmasq^(-1)*theta.small*mean.Y2impute
          mean.Y2impute = sigma2sq^(-1)*V2.part + mean.Y2impute
          
          mean.Y2impute = var.Y2impute*mean.Y2impute 
          
          imputed.res2 = rnorm( 1, mean.Y2impute, sqrt.var.Y2impute)
          
          ############
          #  Output  #
          ############
          
          Design.To = c( 1, imputed.res2, V1.part[-1] )
          
          AcceptOption[1] = 2 # We imputed a record for file 2 
          Added$To        = imputed.res2 # We have imputed this value 
          
          orderC.1T       = 1
          orderC.2T       = 1 
          
          n.to.prop = 1 
          
          ######################################################################################################
          # Step 2: Creating the To Block: File 1 Record Move : Non Empty Block  
          ######################################################################################################
          
        } else{     
          
          F1.records  = RIB[[1]][[ holder.to ]] 
          F2.records  = RIB[[2]][[ holder.to ]]
          
          orderC.1T.orig    = BlockRow.F1[ F1.records ]
          orderC.2T.orig    = BlockRow[ F2.records ]
          
          # Store the Original Design To
          
          Design.To.Orig = D.Holder[F2.records,]
          if( class(Design.To.Orig) != "matrix"){
            Design.To.Orig = c(1,Design.To.Orig)
          } else{ 
            Design.To.Orig = cbind(1,Design.To.Orig)
          }
          
          y.to.orig     = Y1[F1.records]
          
          F1.records.to = F1.records
          F2.records.to = F2.records  
          
          n.to           = length(F1.records.to)
          num.matched.to = count.matches[holder.to]
          
          # Original To
          exact[1] = ifelse( n.to-num.matched.to > threshold, 0, ifelse( n.to-num.matched.to ==0 | n.to-num.matched.to ==1 , 2, 1) )
          
          ####################################################################################################
          # Step 2: Creating the To Block: File 1 Record Move : Replace an Imputation in the To Block
          ####################################################################################################
          
          imputed.inBlock = Imputed.Block[ holder.to, holder.file ]
          
          if( imputed.inBlock > 0){
            
            if( Comment == "YES"){
              print("ToImputation")
            }
            
            # In this case, the number of record pairs in the to block does not change, and the File 2 info does not change. 
            
            AcceptOption[1] = 11 # We replace an imputation from File 1 
            
            # We have imputed.inBlock imputed records in the block. Choose 1. 
            where.replacing       = sample( imputed.inBlock, 1 ) 
            # Of the "in block" records in the current ordering, which record is this? 
            where.replacing       = Imputed.Records[[holder.to]][[1]][where.replacing]
            Deleted$To            = F1.records[where.replacing]      
            
            
            # Keep track of the records that we propose to put in this block
            F1.recordsT                  = F1.records
            F1.recordsT[where.replacing] = holder.r
            F2.recordsT                  = F2.records 
            
            # The current File 1 information 
            y.to.orig             = y.to.orig[order(orderC.1T.orig)]
            y.to                  = y.to.orig 
            y.to[where.replacing] = Y1[ holder.r ]
            
            # Do we sample the proposed to block exactly? 
            exact[3] = exact[1] 
            n.to.prop = n.to 
            
            # The current File 2 information (unchanged) 
            
            if( n.to > 1){
              Design.To.Orig = Design.To.Orig[order(orderC.2T.orig),]
            }
            Design.To      = Design.To.Orig
            
            orderC.2T    = orderC.2T.orig
            orderC.1T    = orderC.1T.orig  
            
            
            ####################################################################################################
            # Step 2: Creating the To Block: File 1 Record Move : No records to replace in the to block 
            ####################################################################################################    
            
          } else{
            
            if( Comment == "YES"){
              print("NoImputationTo")
            }
            
            AcceptOption[1] = 2 # We impute a record for File 2 
            # Add the new record to the list of those in the block 
            F1.recordsT     = c( F1.records, holder.r )
            n.to.prop       = n.to + 1 
            F2.recordsT     = c(F2.records,-99)
            
            # Do we sample the proposed block exactly? 
            exact[3] = ifelse(  n.to.prop-num.matched.to < threshold, 1, 0 )
            # The proposed File 1 information  
            y.to.orig       = y.to.orig[order(orderC.1T.orig)]
            Y1.here         = Y1[ holder.r]
            y.to            = c(y.to.orig, Y1.here  )
            
            # Grab the shell of a record which defines the block. 
            to.addin         = Imp.Holder2[[holder.to]]
            
            #############################
            #Impute a record for file 2. 
            #############################
            
            V1.part = c(1, to.addin[where.coefs.primary[-1]]) 
            
            if( secondary.option == "independent"){
              V2.part = eta
            } else{ 
              V2.part = c(1, to.addin[where.coefs.secondary])
              V2.part = V2.part%*%eta 
            } 
            
            mean.Y2impute = Y1.here - V1.part%*%beta.impute 
            mean.Y2impute = sigmasq^(-1)*theta.small*mean.Y2impute
            mean.Y2impute = sigma2sq^(-1)*V2.part + mean.Y2impute
            
            mean.Y2impute = var.Y2impute*mean.Y2impute 
            
            imputed.res2 = rnorm( 1, mean.Y2impute, sqrt.var.Y2impute)
            
            if(Comment=="YES"){
              print(imputed.res2)
            }    
            
            # Before adjusting the ordering, put in the File 2 information
            
            # Impute a record from file 2 to balance 
            if( n.to > 1 ){
              Design.To.Orig = Design.To.Orig[order(orderC.2T.orig),]
            }
            Design.To      = rbind(Design.To.Orig,c(1,to.addin[where.coefs.primary]) )   
            
            # By default, new record placed on the bottom, so the order is now whatever it was with one more.
            orderC.1T = c(orderC.1T.orig, n.to.prop)
            orderC.2T = c(orderC.2T.orig, n.to.prop)
            
            Design.To[n.to.prop,which.Y2] =  imputed.res2   
            Added$To = imputed.res2
            
            colnames(Design.To)[2] = Y2.name   
            
          } 
        }
        # At the completion of this step, we have created the to block, as well as the likelihood from the to.block. We now need to do the same thing with the from.block 
        
        ################################################################################
        #  Step 3: Creating the From Block: File 1 Record Move: Create the design matrix  
        ################################################################################
        
        F1.records = RIB[[1]][[ holder.from ]] 
        F2.records = RIB[[2]][[ holder.from ]]
        
        orderC.1F.orig  = BlockRow.F1[ F1.records ]
        orderC.2F.orig  = BlockRow[ F2.records ]
        
        n.from = length(F1.records)
        num.matched.from = count.matches[holder.from]
        
        # Store the Original Design To
        Design.From.Orig = D.Holder[ F2.records, ] 
        if( n.from ==1){
          Design.From.Orig = c(1,Design.From.Orig)
        } else{ 
          Design.From.Orig = cbind(1,Design.From.Orig)
        }
        
        y.from.orig     = Y1[ F1.records ]   
        
        F1.records.From = F1.records
        F2.records.From = F2.records   
        
        if( holder.r %in% F1.records == FALSE ){
          stop( paste("We are proposing to move record",r,"in File 1, but the record is not in the from block",holder.from) )
        }
        
        # Original From
        exact[2] = ifelse( n.from-num.matched.from > threshold, 0, ifelse(  n.from-num.matched.from ==1, 2, 1) )
        
        ######################################################################################################
        # Step 3: Creating the From Block: File 1 Record Move : Records for us to replace 
        ######################################################################################################
        
        imputed.inBlock = Imputed.Block[ holder.from, ifelse(holder.file==1,2,1) ]
        
        if( imputed.inBlock > 0){ # Remember, now we are looking for imputations from the OPPOSITE file
          
          if( Comment == "YES"){
            print("FromImputation")
          }
          
          AcceptOption[2] = 22 # We replace an imputation from File 2 
          
          # Remove the moving record from the block 
          where.movingRecord = which(F1.records==holder.r)
          F1.recordsF        = F1.records[-where.movingRecord] 
          
          ######################################################
          # Was the moving record matched to an imputed record? 
          ######################################################
          
          lam         = F1.Lambda[ holder.r]
          which.File2 = which(F2.Lambda[F2.records] == lam )
          was.imputed = (File2_data[ F2.records[which.File2], "Imputed"] == 1)
          
          if( was.imputed == "TRUE" ){
            
            if( Comment == "YES"){
              print("WasImputed")
            }
            
            #################################
            ### Check For Empty Blocks    ###
            #################################
            
            if( n.from == 1){
              
              if( Comment == "YES"){
                print("From Block Now Empty")
              }
              
              #There was only one record from file 2 in the block, and it was imputed . 
              # We now declare the block empty. 
              EmptyFrom  = 1 
              exact[4] = 3
              Deleted$From = F2.records
              AcceptOption[2] = 22 
              y.from.orig     = y.from.orig[order(orderC.1F.orig)]
              one.record.only = 1
              
              #rm(F2.recordsF)
              
            } else{ 
              
              #################################
              ###  The Block is not empty #####
              #################################
              
              # The record that we are moving out was matched to an imputed record
              # Remove the corresponding row from the block 
              F2.recordsF            = F2.records[-which.File2]
              Deleted$From           = F2.records[which.File2]
              
              
              n.from.prop            = n.from - 1 
              
              # Remove the record that we have moved out 
              y.from                = y.from.orig[ -where.movingRecord ]
              # Reorganize 
              y.from.orig           = y.from.orig[order(orderC.1F.orig)]
              removed               = orderC.1F.orig[where.movingRecord]
              orderC.1F             = orderC.1F.orig[-where.movingRecord]
              orderC.1F             = ifelse( orderC.1F > removed, orderC.1F-1, orderC.1F)
              y.from                = y.from[ order(orderC.1F) ]
              
              #Design.From changes == check it!! 
              Design.From           = D.Holder[ F2.recordsF, ]
              removed               = orderC.2F.orig[which.File2]
              orderC.2F             = orderC.2F.orig[-which.File2]
              orderC.2F            = ifelse( orderC.2F >removed, orderC.2F - 1, orderC.2F)
              
              Design.From.Orig      = Design.From.Orig[order(orderC.2F.orig),]
              
              if( n.from.prop==1 ){
                Design.From = c(1,Design.From)
                one.record.only = 1 
                exact[4] = 2
              } else{
                Design.From = Design.From[order(orderC.2F),]
                Design.From = cbind(1,Design.From)
                if( n.from.prop-num.matched.from == 0){
                  exact[4] = 2
                } else{ 
                  exact[4] = ifelse( n.from.prop-num.matched.from < threshold, 1, 0 )
                }
              }    
              #n.from.prop = n.from.prop-num.matched.from
              
            }
            
            ######################################################################################################
            # Step 3: Creating the From Block: File 1 Record Move : Imputation to Replace 
            ######################################################################################################
            
          } else{
            
            ##################################################
            # Not matched to an imputed record, which imputed records are we replacing? 
            ##################################################
            
            # Record from F2 replaced by the non-imputed File2 previously matched to the moved File 1 record
            where.replacing = sample( imputed.inBlock, 1 )  
            where.replacing = Imputed.Records[[holder.from]][[2]][where.replacing]
            n.from.prop     =  n.from - 1 
            
            Deleted$From    = F2.records[where.replacing]
            
            # Which record from File 2 was matched to the one from File 1 that we moved out? 
            #which.File2        
            
            # Remove the record that we have moved out
            y.from               = y.from.orig[-where.movingRecord]
            y.from.orig          = y.from.orig[order(orderC.1F.orig)]
            
            # Remove the row corresponding to the moving record
            F2.recordsF     = F2.records[-where.replacing] 
            
            # Now the order has to change 
            removed         = orderC.1F.orig[where.movingRecord]
            orderC.1F       = orderC.1F.orig[-where.movingRecord]
            orderC.1F       = ifelse( orderC.1F>removed, orderC.1F - 1, orderC.1F)
            removed         = orderC.2F.orig[where.replacing]
            orderC.2F       = orderC.2F.orig[-where.replacing]
            orderC.2F       = ifelse( orderC.2F >removed, orderC.2F - 1, orderC.2F)
            
            # The File 1 information
            # CHANGED: The below used to have y.from.orig
            Design.From.Orig      = Design.From.Orig[order(orderC.2F.orig),]
            y.from      = y.from[ order(orderC.1F) ] 
            if( n.from.prop == 1 ){
              Design.From = File2_data[F2.recordsF,][where.coefs.primary]
              Design.From = c( 1, Design.From ) 
              one.record.only = 1
              exact[4] = 2 
            } else{ 
              Design.From = File2_data[F2.recordsF,where.coefs.primary]
              Design.From = cbind( 1, Design.From ) 
              Design.From = Design.From[ order(orderC.2F), ]
              if( n.from.prop-num.matched.from == 0){
                exact[4] = 2
              } else{ 
                exact[4] = ifelse( n.from.prop-num.matched.from < threshold, 1, 0 )
              }
            }
            
            
            #n.from.prop = n.from.prop-num.matched.from
            
          }
          
          ######################################################################################################
          # Step 3: Creating the From Block: File 1 Record Move : There are no imputations to replace 
          ######################################################################################################
          
        } else{   
          
          if( Comment == "YES"){
            print("NoImputationFrom")
          }
          
          # In this case, the number of total record pairs in the block does not change.
          n.from.prop = n.from
          exact[4]=exact[2]
          
          AcceptOption[2] = 1 # We impute a record from File 1 
          
          ######################################################################################
          # If no imputations to replace, impute a new value for the dangling Record from File 2 
          ######################################################################################
          
          # No imputations for File2. Impute a new record for File 1 to replace the dangling record in File 2 
          
          # The current record ordering: pre-removal
          y.from          = y.from.orig
          y.from.orig     = y.from.orig[order(orderC.1F.orig)]
          
          # Where was the record we are replacing?
          
          if( sum( F1.records==holder.r) ==0 ){
          stop( paste("We are proposing to move record",r,"in File 1, but the record is not in the from block",holder.from) )
          }
          where.replace              = replace.where = which( F1.records==holder.r )
          F1.recordsF                = F1.records
          F1.recordsF[where.replace] = -99
          y.from[where.replace]      = NA 
          F2.recordsF                = F2.records 
          
          # The File 2 information 
          
          if(  n.from > 1){
            check.length     = 0
            Design.From.Orig = Design.From.Orig[order(orderC.2F.orig),]
            Design.From      = Design.From.Orig
            
          } else{
            one.record.only  = 1
            Design.From      = Design.From.Orig
          }
          
          orderC.1F = orderC.1F.orig
          orderC.2F = orderC.2F.orig
          
          # The information we need for the imputation
          
          imputed.y = rnorm( n.from , Design.From%*%beta.star, sigmasq )
          
          y.from[where.replace] = Added$From = imputed.y[where.replace]
          
          if(Comment=="YES"){
            print(Added$From)
          }
          
          
          y.from = y.from[order(orderC.1F.orig)]
          
          #n.from.prop = n.from.prop-num.matched.from
          
        } 
        
        # At the completion of this step, we have created the from block.
        
        ########################################################################################################          
        #################################### File 2 Record Move  ###############################################
        ########################################################################################################
        
      } else{ 
        
        ##############################################################################################
        # Step 2: Creating the To Block: File 2 Record Move :  Empty To Block 
        ##############################################################################################
        
        #################################
        ###  If the to block is empty ###
        #################################
        
        if( EmptyTo == "TRUE" ){
          
          exact[1] = 3
          
          if( Comment == "YES"){
            print("EmptyTo")
          }
          
          AcceptOption[1] = 1 # We impute a record for file 1 
          Added$To        = 1 
          # Move in the target record from File 2
          to.addin     = Imp.Holder2[[holder.to]]
          
          to.addin     = c(1,Y2.in[holder.r],to.addin[where.coefs.primary][-1])  
          
          #Impute a record for file 1.  
          Design.To   = to.addin
          imputed.res = rnorm( 1 , Design.To%*%beta.star, sigma )
          F2.recordsT = holder.r
          F1.recordsT = -99    
          y.to        = imputed.res 
          
          orderC.1T       = 1
          orderC.2T       = 1 
          
          one.record.only = 1 
          exact[3]        = 2 
          
          n.to.prop       = 1 
          
          ##############################################################################################
          # Step 2: Creating the To Block: File 2 Record Move :  Non Empty To Block 
          ##############################################################################################
          
        } else{ 
          
          F1.records  = RIB[[1]][[ holder.to ]] 
          F2.records  = RIB[[2]][[ holder.to ]]
          
          orderC.1T.orig    = BlockRow.F1[ F1.records ]
          orderC.2T.orig    = BlockRow[ F2.records ]
          
          # Store the Original Design To
          Design.To.Orig = D.Holder[F2.records,]
          if( class(Design.To.Orig) != "matrix"){
            Design.To.Orig = c(1,Design.To.Orig)
          } else{ 
            Design.To.Orig = cbind(1,Design.To.Orig)
          }        
          
          y.to.orig     = Y1[F1.records]
          
          F1.records.to = F1.records
          F2.records.to = F2.records
          
          n.to           = length(F1.records.to)
          num.matched.to = count.matches[holder.to]
          
          # Original To
          exact[1] = ifelse( n.to-num.matched.to > threshold, 0, ifelse( n.to-num.matched.to ==0 | n.to-num.matched.to ==1 , 2, 1) )
          
          ##############################################################################################
          # Step 2: Creating the To Block: File 2 Record Move :  Imputations to Replace 
          ##############################################################################################
          
          # We are creating a to block, meaning that we are moving in a record for File 2. 
          
          imputed.inBlock = Imputed.Block[ holder.to, holder.file ]
          
          # Is there a record in the block from file 2 which is an imputation? 
          
          if( imputed.inBlock > 0){
            
            # In this case, the number of records in the block does not change. 
            exact[3] = exact[1] 
            
            if( Comment == "YES"){
              print("ImputationTo")
            }
            
            AcceptOption[1] = 22 # We replace an imputation from File 2 
            
            # Choose an imputed record to replace
            where.replacing       = sample( imputed.inBlock, 1 ) 
            # Of the "in block" records in the current ordering, which record is this? 
            where.replacing       = Imputed.Records[[ holder.to ]][[2]][ where.replacing ]
            Deleted$To            = F2.records[where.replacing]
            
            F1.recordsT                  = F1.records 
            F2.recordsT                  = F2.records
            F2.recordsT[where.replacing] = holder.r 
            orderC.1T                    = orderC.1T.orig
            orderC.2T                    = orderC.2T.orig           
            
            # The File 1 information (unchanged)  
            y.to.orig     = y.to.orig[order(orderC.1T)]
            y.to          = y.to.orig
            
            # Create V2.proposed 
            Imp.Holder.here = Imp.Holder2[[holder.to]]
            
            # The current File 2 information ( pre-replacment, pre-ordering ) 
            # Replace the chosen imputation with holder.r 
            if( n.to > 1 ){
              Design.To                           = Design.To.Orig 
              Design.To[where.replacing,which.Y2] = File2_data[holder.r,col.Y2] 
              # Now order according to the current permutation 
              Design.To.Orig = Design.To.Orig[order(orderC.2T.orig),]
              Design.To      = Design.To[order(orderC.2T.orig),]
            } else{
              Design.To           = Design.To.Orig
              Design.To[which.Y2] = File2_data[ holder.r, col.Y2 ] 
              one.record.only     = 1
              exact[3]            = 2
            }
            
            n.to.prop = n.to 
            
            
            ##############################################################################################
            # Step 2: Creating the To Block: File 2 Record Move :  No Imputations to Replace 
            ##############################################################################################
            
          } else{
            
            if( Comment == "YES"){
              print("NoImputationTo")
            }
            
            AcceptOption[1] = 1 # We impute a record for File 1 
            
            # Add the new record from File 2 to the list of those in the block 
            F2.recordsT   = c( F2.records, holder.r )
            F1.recordsT   = c( F1.records,-99 )
            n.to.prop     = n.to + 1 
            
            exact[3] = ifelse( n.to.prop-num.matched.to > threshold, 0, 1 )
            
            # The current File 1 information (pre-imputation) 
            y.to.orig   = y.to.orig[order(orderC.1T.orig)]
            y.to        = y.to.orig
            orderC.1T   = c(orderC.1T.orig,n.to.prop)
            
            # The proposed File 2 information
            Design.To   = Design.To.Orig
            if( n.to > 1){
              Design.To.Orig = Design.To.Orig[order(orderC.2T.orig),]
            }
            
            # Make a holder for moving in the new record
            add.in        = Imp.Holder2[[holder.to]] 
            add.in[col.Y2] = Y2.in[holder.r]
            to.add = c(1,add.in[where.coefs.primary])
            Design.To     = rbind(Design.To, to.add )   
            
            # We put in a new record, need a new ordering
            orderC.2T   = c(orderC.2T.orig,n.to.prop)
            
            
            Design.To   = Design.To[ order(orderC.2T),]
            
            # Impute a record from file 1 to balance     
            y.to             = c( y.to, NA )     
            
            imputed.y        = rnorm( n.to.prop, Design.To%*%beta.star, sigma ) 
            
            y.to[n.to.prop]  = Added$To = imputed.y[n.to.prop]  
            
            if(Comment=="YES"){
              print(Added$To)
            }
            
          } 
        }
        # At the completion of this step, we have created the to block. 
        
        ######################################################################################################
        # Step 3: Creating the From Block: File 2 Record Move : Matrices and Constants  
        ######################################################################################################
        
        F1.records = RIB[[1]][[ holder.from ]] 
        F2.records = RIB[[2]][[ holder.from ]]
        
        orderC.1F.orig  = BlockRow.F1[ F1.records ]
        orderC.2F.orig  = BlockRow[ F2.records ]
        
        # Store the Original Design From
        Design.From.Orig = D.Holder[ F2.records, ] 
        if( length(F2.records)==1){
          Design.From.Orig = c(1,Design.From.Orig)
        } else{ 
          Design.From.Orig = cbind(1,Design.From.Orig)
        }      
        
        y.from.orig     = Y1[ F1.records ]   
        
        F1.records.From = F1.records
        F2.records.From = F2.records   
        
        n.from           = length(F1.records.From)
        num.matched.from = count.matches[holder.from]
        
        # Original From
        exact[2] = ifelse( n.from-num.matched.from > threshold, 0, ifelse(  n.from-num.matched.from ==1, 2, 1) )
        
        
        if( holder.r %in% F2.records == FALSE ){
          stop( paste("We are proposing to move record",r,"in File 2, but the record is not in the from block",holder.from) )
        }
        
        
        
        ######################################################################################################
        # Step 3: Creating the From Block: File 2 Record Move : Check for Imputations  
        ######################################################################################################
        
        imputed.inBlock = Imputed.Block[ holder.from, ifelse(holder.file==1,2,1) ]
        
        if( imputed.inBlock > 0){ # Remember, now we are looking for imputations from the OPPOSITE file
          
          if( Comment == "YES"){
            print("ImputationFrom")
          }          
          
          AcceptOption[2] = 11 # We replace an imputation from File 1 
          
          # Remove the moving record from the block 
          where.movingRecord = which(F2.records==holder.r)
          F2.recordsF        = F2.records[-where.movingRecord] 
          
          ######################################################
          # Was the moving record matched to an imputed record? 
          ######################################################
          
          lam         = F2.Lambda[holder.r]
          which.File1 = which(F1.Lambda[F1.records] == lam ) # Has to have a match, so non-empty 
          was.imputed = (File1_data[lam, "Imputed"] == 1)
          
          if( was.imputed == "TRUE" ){            
            
            #################################
            ### Check For Empty Blocks    ###
            #################################
            
            if( n.from == 1){
              
              #There was only one record from file 1 in the block, and it was imputed. 
              #We now declare the block empty. 
              EmptyFrom       = 1 
              Deleted$From    = F1.records
              y.from.orig     = y.from.orig[order(orderC.1F.orig)]
              one.record.only = 1
              exact[4]        = 3
              
            } else{ 
              
              #################################
              ###  The Block is not empty #####
              #################################
              
              # The record that we are moving out was matched to an imputed record
              # Remove the corresponding row from the block 
              Deleted$From      = F1.records[ which.File1 ] 
              F1.recordsF       = F1.records[-which.File1] 
              
              
              y.from            = y.from.orig[ - which.File1 ]  
              y.from.orig       = y.from.orig[orderC.1F.orig]
              removed           = orderC.1F.orig[which.File1]
              orderC.1F         = orderC.1F.orig[-which.File1] 
              orderC.1F         = ifelse(orderC.1F > removed, orderC.1F-1,orderC.1F )
              y.from            = y.from[ order(orderC.1F) ]
              
              Design.From       = D.Holder[ F2.recordsF, ]
              removed           = orderC.2F.orig[which.File1]
              orderC.2F         = orderC.2F.orig[-which.File1] 
              orderC.2F         = ifelse(orderC.2F > removed, orderC.2F-1,orderC.2F )
              
              # The number of record pairs goes down by one.
              n.from.prop       = n.from - 1 
              Design.From.Orig  = Design.From.Orig[order(orderC.2F.orig),]
              
              if( n.from.prop ==1){
                Design.From     = c(1,Design.From)
                one.record.only =  1               
                exact[4] = 2
              } else{
                Design.From      = Design.From[order(orderC.2F),]
                Design.From      = cbind(1,Design.From)
                if( n.from.prop-num.matched.from == 0){
                  exact[4] = 2
                } else{ 
                  exact[4] = ifelse( n.from.prop-num.matched.from < threshold, 1, 0 )
                }
              }
            }
            
          } else{
            
            ##################################################
            # If not, which imputed records are we replacing? 
            ##################################################
            
            # Imputed record from F1 replaced with non-imputed F1 previously matched to the moved F2 record 
            where.replacing = sample( imputed.inBlock, 1 )  
            where.replacing = Imputed.Records[[holder.from]][[1]][where.replacing]
            Deleted$From    = F1.records[where.replacing]
            n.from.prop     = n.from - 1
            
            # Which record from File 1 was matched to the one from File 2 that we moved out? 
            #which.File1  
            y.from                = y.from.orig 
            y.from[ which.File1 ] = y.from.orig[where.replacing]
            y.from.orig           = y.from.orig[order(orderC.1F.orig)]
            y.from                = y.from[-where.replacing]
            
            F1.recordsF           = F1.records[-where.replacing]
            
            # Now the order has to change 
            removed         = orderC.1F.orig[where.replacing] 
            orderC.1F       = orderC.1F.orig[-where.replacing] 
            orderC.1F       = ifelse(orderC.1F > removed, orderC.1F-1,orderC.1F )
            removed         = orderC.2F.orig[where.movingRecord]
            orderC.2F       = orderC.2F.orig[-where.movingRecord]
            orderC.2F       = ifelse(orderC.2F > removed, orderC.2F-1,orderC.2F )
            y.from          = y.from.orig[ order(orderC.1F) ] 
            
            # Remove the row corresponding to the moving record
            Design.From.Orig = Design.From.Orig[order(orderC.2F.orig),]
            
            if( n.from.prop == 1 ){
              
              Design.From = D.Holder[F2.recordsF,]
              Design.From = c( 1, Design.From ) 
              one.record.only = 1 
              exact[4] = 2
              
              
            } else{ 
              Design.From = D.Holder[F2.recordsF,]
              Design.From = Design.From[ order(orderC.2F), ]
              Design.From = cbind( 1, Design.From ) 
              if( n.from.prop-num.matched.from == 0){
                exact[4] = 2
              } else{ 
                exact[4] = ifelse( n.from.prop-num.matched.from < threshold, 1, 0 )
              }
            }  
            
          }
          
        } else{
          
          if( Comment == "YES"){
            print("NoImputationFrom")
          }
          
          AcceptOption[2] = 2 # We impute a record for File 2  
          
          ##################################################
          # There are no imputations to replace 
          ##################################################
          
          # We have to impute a new record for File 2 to replace the dangling record in File 1 
          
          # In this case the number of records in the block stays the same
          exact[4] = exact[2] 
          
          # Remove the moving record from the block 
          
          where.impute  = which(F2.records==holder.r)
          F2.recordsF   = F2.records
          F2.recordsF[where.impute] = -99
          F1.recordsF   = F1.records 
          orderC.1F     = orderC.1F.orig
          orderC.2F     = orderC.2F.orig
          
          # The F1 record ordering (unchanged) 
          y.from        = Y1[F1.records][ order(orderC.1F) ]
          y.from.orig   = y.from.orig[order(orderC.1F.orig)]
          
          # Where was the record we are replacing?
          where.impute  = which( orderC.2F.orig== where.impute)
          
          # Grab the shell of a record which defines the block. 
          for.imputing  = Imp.Holder2[[holder.from]]
          
          #############################
          #Impute a record for file 2. 
          #############################
          to.addin = Imp.Holder2[[holder.from]]
          Y1.part  = y.from[where.impute]
          V1.part  = c(1, to.addin[where.coefs.primary[-1]]) 
          
          if( secondary.option == "independent"){
            V2.part = eta
          } else{ 
            V2.part = c(1, to.addin[where.coefs.secondary])
            V2.part = V2.part%*%eta 
          } 
          
          mean.Y2impute = Y1.part - V1.part%*%beta.impute 
          mean.Y2impute = sigmasq^(-1)*theta.small*mean.Y2impute
          mean.Y2impute = sigma2sq^(-1)*V2.part + mean.Y2impute
          
          mean.Y2impute = var.Y2impute*mean.Y2impute 
          
          imputed.reg2 = rnorm( 1, mean.Y2impute, sqrt.var.Y2impute)
          
          if(Comment=="YES"){
            print(imputed.reg2)
          }
          
          n.from.prop = n.from 
          
          # The File 2 information 
          
          if(  n.from > 1){
            check.length = 0
            Design.From  = D.Holder[F2.records,]
            Design.From  = Design.From[ order(orderC.2F), ]
            Design.From  = cbind(1,Design.From) 
            
            Design.From[where.impute,Y2.name] = Added$From = imputed.reg2
            
            #colnames(Design.From)=c("1",Y2.name, names(Design))    
            
            Design.From.Orig = Design.From.Orig[order(orderC.2F.orig),]
            
          } else{
            one.record.only= 1
            Design.From = File2_data[F2.records,][where.coefs.primary]
            Design.From = c(1,Design.From)          
            Design.From[where.impute] = Added$From = imputed.reg2
            #names(Design.From)=c("1",Cov_NonBl, Cov_Bl)    
            
          }  
          #n.from.prop = n.from.prop-num.matched.from
          
          
        }
        
        
      } 
      
      #####################################################################################
      ###                      Begin to Compute the Ratio                               ###
      #####################################################################################      
      
      log.A.ratio = log(holder.A) # The part of the ratio that does NOT depend on the MH sampled C part 
      
      #####################################################################################
      ### Here we are dealing with just the Acceptance Ratio part that has to do with C ###
      #####################################################################################
      
      exactALL = 0 
      
      c.log.A.ratio = 0           # The part of the ratio that DOES depend on the MH sampled C part 
      
      ### I need to compute the part of the acceptance probability for each of these for the initial accept/reject 
      
      ###################
      ### Original To ###
      ###################
      
      if( exact[1] == 0 ){
        
        # Before the move, the to block had ABOVE the threshold number of record pairs 
        
        like.part.Orig.To = dmvnorm( y.to.orig, Design.To.Orig%*%beta.star, sigmasq*diag(n.to),log=TRUE)
        log.A.ratio       = log.A.ratio - like.part.Orig.To
        
        # The Add On
        to.counter        = non.seed.RIB[holder.to]
        log.A.ratio       = log.A.ratio + log( 2*factorial(to.counter -2) )
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.to]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.to), ncol = p2 , byrow= TRUE)
          F2.part         = dmvnorm( Design.To.Orig[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.to),log=TRUE)
          log.A.ratio     = log.A.ratio - F2.part 
        }         
        
      } else if ( exact[1] == 1 ){
        
        # Before the move, the to block had BELOW the threshold number of record pairs
        
        C.part.Orig.To      = exactSampleC( y.to.orig, orderC.2T.orig, Design.To.Orig, n.to, beta.star, sigma, holder.to, count.matches[holder.to]) 
        
        to.counter          = non.seed.RIB[holder.to]
        Add.On              = log( factorial(to.counter) ) - C.part.Orig.To$Sum
        log.A.ratio         = log.A.ratio + Add.On         
        exactALL = exactALL + 1 
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.to]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.to), ncol = p2 , byrow= TRUE)
          Add.On          = dmvnorm( Design.To.Orig[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.to),log=TRUE)
          log.A.ratio     = log.A.ratio - Add.On 
        } 
        
        rm(Add.On) 
        
      } else if (exact[1] == 2 ){
        
        # Before the move, the to block had only one non.seed record pair. 
        
        like.part.Orig.To = dmvnorm( y.to.orig, Design.To.Orig%*%beta.star, sigmasq*diag(n.to),log=TRUE)
        log.A.ratio       = log.A.ratio - like.part.Orig.To
        
        exactALL = exactALL + 1
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.to]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.to), ncol = p2 , byrow= TRUE)
          if(n.to > 1){
            F2.part     = dmvnorm( Design.To.Orig[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.to),log=TRUE)
          } else{
            F2.part     = dnorm( Design.To.Orig[which.Y2], D.matrix.holder%*%eta, sigma2,log=TRUE)
          }
          log.A.ratio     = log.A.ratio - F2.part 
        } 
        
      } else if(exact[1]==3){
        # The To block was empty
        exactALL = exactALL + 1
      }
      
      # else the original to block was empty, and we have no contribution to the ratio 
      
      ###################
      ## Original From ##
      ###################
      
      if( exact[2] == 0 ){
        
        # Before the move, the from block had ABOVE the threshold number of record pairs 
        
        like.part.Orig.From = dmvnorm( y.from.orig, Design.From.Orig%*%beta.star, sigmasq*diag(n.from),log=TRUE)
        log.A.ratio         = log.A.ratio  - like.part.Orig.From
        
        # The Add On
        from.counter        = non.seed.RIB[holder.from]
        log.A.ratio         = log.A.ratio + log( 2*factorial(from.counter -2) )
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.from]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.from), ncol = p2 , byrow= TRUE)
          Add.On          = dmvnorm( Design.From.Orig[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.from),log=TRUE)
          log.A.ratio     = log.A.ratio - Add.On 
          rm(Add.On)
        } 
        
      } else if ( exact[2] == 1 ){
        
        # Before the move, the from block had BELOW the threshold number of record pairs 
        
        C.part.Orig.From     = exactSampleC( y.from.orig, orderC.2F.orig, Design.From.Orig, n.from , beta.star, sigma, holder.from, count.matches[holder.from]) 
        #orderC.2F.Orig.prob = C.part.Orig.From$Apart 
        #print(C.part.Orig.From)
        
        from.counter        = non.seed.RIB[holder.from]
        Add.On              = log( factorial(from.counter) ) - C.part.Orig.From$Sum
        log.A.ratio         = log.A.ratio + Add.On         
        exactALL = exactALL + 1 
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.from]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.from), ncol = p2 , byrow= TRUE)
          Add.On          = dmvnorm( Design.From.Orig[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.from),log=TRUE)
          log.A.ratio     = log.A.ratio - Add.On 
        } 
        
        rm(Add.On) 
        
      } else if (exact[2] == 2 ){
        
        # Before the move, the from block had only one record pair  
        
        like.part.Orig.From = dmvnorm( y.from.orig, Design.From.Orig%*%beta.star, sigmasq*diag(n.from),log=TRUE)
        log.A.ratio        = log.A.ratio - like.part.Orig.From
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.from]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.from), ncol = p2 , byrow= TRUE)
          if(n.from > 1){
            F2.part   = dmvnorm( Design.From.Orig[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.from),log=TRUE)
          } else{
            F2.part     = dnorm( Design.From.Orig[which.Y2], D.matrix.holder%*%eta, sigma2,log=TRUE)
          }
          log.A.ratio     = log.A.ratio - F2.part 
          rm(F2.part)
        } 
        
        exactALL = exactALL + 1 
        
      } 
      
      ###################
      ### Proposed To ###
      ###################
      
      if( exact[3] == 0 ){  
        
        # AFTER the move, the to block had ABOVE the threshold number of record pairs     
        
        if( num.matched.to > 0 ){
          switch.options = orderC.2T[orderC.2T > num.matched.to]
        } else{
          switch.options = orderC.2T
        } 
        
        chosen               = sample( switch.options, 2, replace= F )
        order.out.To         = orderC.2T
        order.out.To[chosen] = order.out.To[ chosen[2:1] ] 
        Proposed             = Design.To[ order.out.To, ] 
        
        like.part.Prop       = dmvnorm( y.to, Proposed%*%beta.star, sigmasq*diag(n.to.prop),log=TRUE)
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.to]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.to.prop), ncol = p2 , byrow= TRUE)
          Add.On          = dmvnorm( Design.To[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.to.prop),log=TRUE)
          log.A.ratio     = log.A.ratio + Add.On 
        } 
        
        c.log.A.ratio        = like.part.Prop
        
        # The Add On
        to.counter           = n.to.prop - count.matches[holder.to]
        log.A.ratio          = log.A.ratio - log( 2*factorial(to.counter -2) )
        
        to.order            = order.out.To
        
      } else if ( exact[3] == 1 ){  
        
        # After the move, the to block has BELOW the threshold number of record pairs
        
        C.part.Prop.To     = exactSampleC( y.to, orderC.2T, Design.To, n.to.prop, beta.star, sigma, holder.to, count.matches[holder.to]) 
        order.out.To       = C.part.Prop.To$order.out
        #orderC.2T.prob    = C.part.Prop.To$Apart
        #print(C.part.Prop.To)
        #c.log.A.ratio     = c.log.A.ratio + C.part.Prop.To$LikePart - C.part.Prop.To$Prob 
        to.counter         = n.to.prop - count.matches[holder.to]
        Add.On             = C.part.Prop.To$Sum - log( factorial(to.counter) )
        log.A.ratio        = log.A.ratio + Add.On  
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.to]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.to.prop), ncol = p2 , byrow= TRUE)
          Add.On          = dmvnorm( Design.To[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.to.prop),log=TRUE)
          log.A.ratio     = log.A.ratio + Add.On 
        }    
        
        to.order           = order.out.To
        exactALL           = exactALL + 1 
        
      } else if (exact[3] == 2 ){
        
        # After the move, the to block has only one record pair, or all seeds
        
        like.part.Prop.To = dmvnorm( y.to, Design.To%*%beta.star, sigmasq*diag(n.to.prop),log=TRUE)
        log.A.ratio       = log.A.ratio + like.part.Prop.To
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.to]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.to.prop), ncol = p2 , byrow= TRUE)
          if(n.to.prop > 1){
            Add.On   = dmvnorm( Design.To[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.to.prop),log=TRUE)
          } else{
            Add.On     = dnorm( Design.To[which.Y2], D.matrix.holder%*%eta, sigma2,log=TRUE)
          }
          log.A.ratio     = log.A.ratio + Add.On 
        }    
        
        to.order          = orderC.2T
        
        exactALL = exactALL + 1 
        
      } 
      
      
      ###################
      ## Proposed From ##
      ###################
      
      if( exact[4] == 0 ){
        
        # After the move, the from block has ABOVE the threshold number of record pairs
        
        if( num.matched.from > 0 ){
          switch.options = orderC.2F[orderC.2F > num.matched.from]
        } else{
          switch.options = orderC.2F
        } 
        
        chosen                 = sample( switch.options, 2, replace= F )
        order.out.From         = orderC.2F
        order.out.From[chosen] = order.out.From[ chosen[2:1] ] 
        Proposed               = Design.From[order.out.From,]
        
        like.part.Prop         = dmvnorm( y.from, Proposed%*%beta.star, sigmasq*diag(n.from.prop),log=TRUE)
        
        c.log.A.ratio          = like.part.Prop
        
        # The Add On
        from.counter           = n.from.prop - count.matches[holder.from]
        log.A.ratio            = log.A.ratio - log( 2*factorial(from.counter -2) )
        
        from.order             = order.out.From
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.from]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.from.prop), ncol = p2 , byrow= TRUE)
          Add.On          = dmvnorm( Design.From[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.from.prop),log=TRUE)
          log.A.ratio     = log.A.ratio + Add.On 
        }    
        
        
      } else if ( exact[4] == 1 ){
        
        # After the move, the from block has BELOW the threshold number of record pairs
        
        C.part.Prop.From   = exactSampleC( y.from, orderC.2F, Design.From, n.from.prop, beta.star, sigma, holder.from, count.matches[holder.from]) 
        order.out.From     = C.part.Prop.From$order.out
        #orderC.2F.prob    = C.part.Prop.From$Apart 
        #print(C.part.Prop.From)
        #c.log.A.ratio     = c.log.A.ratio + C.part.Prop.From$LikePart - C.part.Prop.From$Prob 
        from.counter       = n.from.prop - count.matches[holder.from]
        Add.On             = C.part.Prop.From$Sum - log( factorial(from.counter) )
        log.A.ratio        = log.A.ratio + Add.On         
        
        from.order        = order.out.From
        exactALL = exactALL + 1 
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.from]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.from.prop), ncol = p2 , byrow= TRUE)
          Add.On          = dmvnorm( Design.From[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.from.prop),log=TRUE)
          log.A.ratio     = log.A.ratio + Add.On 
        }          
        
      } else if (exact[4] == 2 ){
        
        # After the move, the to block has only one record pair 
        
        like.part.Prop.From = dmvnorm( y.from, Design.From%*%beta.star, sigmasq*diag(n.from.prop),log=TRUE)
        log.A.ratio         = log.A.ratio + like.part.Prop.From
        from.order          = orderC.2F
        
        if( secondary.option2 != "independent"){ 
          # The File 2 contribution 
          D.matrix.holder = c(1,Imp.Holder2[[holder.from]][where.coefs.secondary])
          D.matrix.holder = matrix( rep(D.matrix.holder, n.from.prop), ncol = p2 , byrow= TRUE)
          if(n.from.prop > 1){
            Add.On   = dmvnorm( Design.From[,which.Y2], D.matrix.holder%*%eta, sigma2sq*diag(n.from.prop),log=TRUE)
          } else{
            Add.On     = dnorm( Design.To[which.Y2], D.matrix.holder%*%eta, sigma2,log=TRUE)
          }
          log.A.ratio     = log.A.ratio + Add.On 
        }  
        
        exactALL = exactALL + 1 
        
      } 
      
      # else, the from block is empty 
      
      ######################
      ### Final A Ratio ####
      ######################
      
      Final.log.A.ratio = log.A.ratio + c.log.A.ratio 
      
      if(Comment=="YES"){    
        print(Final.log.A.ratio)
      }
      
      Sampled.Uniform = log(runif(1))
      
      #####################################################################################
      ###          After this, the Acceptance Ratio should be completed                 ###
      #####################################################################################
      
      # Now, whether we fail to accept or not, we are going to continue "walking" in C* space to see if we end up somewhere better now UNLESS we directly sampled all blocks
      # If we directly sampled all blocks, we only accept/reject once. 
      #print(exact)
      if( exactALL == 4){ 
        
        # We directly sampled all of the block permutatations 
        
        if( Sampled.Uniform < Final.log.A.ratio){ 
          if(Comment=="YES"){
            print("ACCEPT")
          }
          accept    = 1 
          #to.order = order.out.To
          #if( EmptyFrom == 0){ from.order = order.out.From}        
        } else{ accept = 0}
        
        
        
      } else{ 
        
        # We did NOT directly sample all of the block permutatations 
        
        if( Sampled.Uniform < Final.log.A.ratio){ 
          # We have accepted the move in E and B space. We continue to "walk" in C* space for reps iterations in the to block and from block independently,
          # excluding those blocks which have been directly sampled 
          accept =1 
          if( exact[4] == 0){
            from.order = GutmanMH( order.out.From ,y.from, Design.From[order.out.From,],n.from.prop, beta.star, sigma, reps, holder.from,count.matches )
          } 
          if( exact[3] == 0){
            to.order   = GutmanMH( order.out.To,y.to, Design.To[order.out.To,],n.to.prop, beta.star, sigma, reps, holder.to,count.matches )
          }
          
        } else{ 
          # We have FAILED to accepted the move in E and B space. We "walk" again in C* space for reps iterations, or until we accept, whatever comes first 
          accept = 0 
          for( r in 1:reps){     
            
            new.log.A.ratio  =0 
            
            if( exact[3] == 0 ){
              
              
              if( num.matched.to > 0 ){
                switch.options = orderC.2T[orderC.2T > num.matched.to] 
              } else{
                switch.options = orderC.2T
              } 
              
              chosen               = sample( switch.options, 2, replace= F )
              order.out.To         = orderC.2T
              order.out.To[chosen] = order.out.To[ chosen[2:1] ] 
              Proposed             = Design.To[ order.out.To, ] 
              
              like.part.Prop       = dmvnorm( y.to, Proposed%*%beta.star,  sigmasq*diag(length(y.to)), log = TRUE )
              
              new.log.A.ratio      = log.A.ratio + like.part.Prop
              
            } 
            
            if( exact[4] == 0 ){            
              # After the move, the from block has ABOVE the threshold number of record pairs
              
              if( num.matched.from > 0 ){
                switch.options = orderC.2F[orderC.2F > num.matched.from]
              } else{
                switch.options = orderC.2F
              } 
              
              chosen   = sample( switch.options, 2, replace= F )
              order.out.From = orderC.2F
              order.out.From[chosen] = order.out.From[ chosen[2:1] ] 
              Proposed  = Design.From[order.out.From,]
              
              like.part.Prop    = dmvnorm( y.from, Proposed%*%beta.star,  sigmasq*diag(n.from.prop), log = TRUE )
              
              new.log.A.ratio    = new.log.A.ratio + like.part.Prop
            }  
            
            if( log(runif(1)) < new.log.A.ratio){
              accept = 1 
              break
            }
          } 
          
          if( accept == 1 & exact[4] == 0 & r != reps){
            from.order = GutmanMH( order.out.From,y.from, Design.From[order.out.From,],n.from.prop, beta.star, sigma, reps-r, holder.from,count.matches )
          } 
          if( accept == 1 & exact[3] == 0 & r!=reps){
            to.order   = GutmanMH( order.out.To,y.to, Design.To[order.out.To,],n.to.prop, beta.star, sigma, reps-r, holder.to,count.matches )
          }
        } 
      } 
      
      #}
      
      #####################################################################################
      ###   We now have determined (1) whether we accept and (2) the permutation C      ###
      #####################################################################################      
      
      #####################################################################################
      ###                      If we do not accept                                      ###
      #####################################################################################
      
      if( accept == 0 ){
        
        if(Comment=="YES"){
          print( "Do Not Accept" ) 
        }
        
        # E(s+1) = E(s), though we do need to keep track of the jstar terms that we proposed 
        
        if( holder.file==1){
          E1[holder.r,"jstar"] = holder.jstar
        } else{ 
          E2[holder.r,"jstar"] = holder.jstar
        }
        
        # B1(s+1) = B1(s)  , B2(s+1) = B2(s)
        
        # C(s+1) = C(s) 
        
      } 
      
      # Note that what we get out of all of this is accept (do we accept the move in E,B space) and the new permutation in C space (from.order, to.order) 
      
      ########################################################################################################
      # If we accept : Now how to we format our output properly?  
      ########################################################################################################
      
      if( accept == 1 ){        
        
        Delete.Counter = 0 # Count the number of imputed records we delete
        Accept.Count   = Accept.Count + 1 # Keeps a running total of accepted moved this iteration
        
        if(holder.to %in% SeedsOnly){
          stop('We are proposing to move record to a block which contains only seeds, but the block seed restriction is on, which does not allow such a proposal.')
        }
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%% Update the Blocking information for the moving record %%#
        
        current.B[[holder.file]][ holder.r, holder.jstar ] = levels(Bhat[,holder.jstar])[holder.all[6]]
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%%     Save the block orders for filled blocks           %%#   
        
        
        if( EmptyFrom != 1){
          orderC.2F = from.order 
        }
        
        if(Comment=="YES"){
          print("Accepted")
        }
        
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%%     Change the Blocking for the moved record          %%#
        
        # All other changes will have to do with imputed records
        
        if( holder.file == 1){
          #File1_data[ holder.r, where.blocking[[1]] ] = Imp.Holder1[[holder.to]][where.blocking[[1]]] 
          File1_data[holder.r,"Block"] = holder.to 
          E1[holder.r,] = Estar[[1]][holder.r,]
        } else{ 
          File2_data[ holder.r, where.bvars.F2]  = Imp.Holder2[[holder.to]][where.bvars.F2] 
          File2_data[holder.r,"Block"] = holder.to 
          E2[holder.r,] = Estar[[2]][holder.r,]
          if( sum(File2_data[,Y2.name] ==0) > 0 ){
            stop("During the first check, one of the Y2 values is 0.")
          }
        } 
        
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%%     Have we added records to the to block?            %%#      
        
        
        if( class(Added$To) != "NULL"  ){          
          
          if( AcceptOption[1] == 1 ){
            # Added a record for File 1 
            added.in                  = Imp.Holder1[[holder.to]]
            added.in[col.Y1]          = Added$To 
            File1_data                = rbind( File1_data, added.in )
            nrow.F1                   = nrow.F1 + 1 
            File1_data[,"lambda"][nrow.F1]= nrow.F1
            F1.recordsT = ifelse( F1.recordsT ==-99, nrow.F1, F1.recordsT)
            File1_data[,"lambda"]= 1:nrow.F1 
          } else{ 
            # Added a record for File 2 
            added.in          = Imp.Holder2[[holder.to]]
            added.in[col.Y2]  = Added$To 
            nrow.F2           = nrow.F2 + 1
            File2_data        = rbind( File2_data, added.in) 
            F2.recordsT = ifelse( F2.recordsT ==-99, nrow.F2, F2.recordsT)
            
          }
          
        }
        
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%%     Have we added records to the from block?          %%#
        
        if( class(Added$From) != "NULL"  ){
          
          if( AcceptOption[2] == 1 ){
            # We have imputed a new file 1 record
            added.in                   = Imp.Holder1[[holder.from]]
            added.in[col.Y1]           = Added$From 
            File1_data                 = rbind( File1_data, added.in )
            nrow.F1                    = nrow.F1 + 1 
            File1_data[,"lambda"][nrow.F1]= nrow.F1  
            F1.recordsF = ifelse( F1.recordsF ==-99, nrow.F1, F1.recordsF)
            File1_data[,"lambda"]= 1:nrow.F1 
          } else{ 
            added.in                    = Imp.Holder2[[holder.from]]
            added.in[col.Y2] = Added$From
            nrow.F2                     = nrow.F2 + 1 
            File2_data                   = rbind( File2_data, added.in)
            F2.recordsF = ifelse( F2.recordsF ==-99, nrow.F2, F2.recordsF)
          }
        }  
        
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%%     If records were deleted from the from block       %%#  
        
        if( class(Deleted$From) != "NULL" ){ 
          
          Delete.Counter = Delete.Counter  + 1 
          
          if( AcceptOption[2] == 11 ){
            # We deleted records from File 1 
            DeletedFrom = Deleted$From
            File1_data  = File1_data[-DeletedFrom,]
            nrow.F1     = nrow.F1 - 1
            #File1_data[,"lambda"]= 1:nrow.F1
            
          } else{ 
            # We deleted records from File 2 
            DeletedFrom = Deleted$From
            File2_data = File2_data[-DeletedFrom, ] 
            nrow.F2    = nrow.F2 -1         
            
          } 
          
          # Note: we now have to define the in Block assignment for the imputed records, since they have shifted around. 
        }
        
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%%     If records were deleted from the TO block       %%#  
        
        
        if( class(Deleted$To) != "NULL" ){
          Delete.Counter = Delete.Counter  + 1 
          # Which file did we delete from? 
          if( AcceptOption[1] == 11 ){
            # If we deleted from File 1
            DeletedTo  = Deleted$To
            File1_data = File1_data[- DeletedTo,]
            nrow.F1    = nrow.F1 - 1             
            #File1_data[,"lambda"]= 1:nrow.F1
          } else{ 
            # If we deleted from File 2 
            DeletedTo  = Deleted$To
            File2_data = File2_data[- DeletedTo, ] 
            nrow.F2    = nrow.F2 -1 
          } 
          # Note: we now have to define the in Block assignment for the imputed records, since they have shifted around. 
        } 
        
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%%          Was the to block Empty?                       %%#    
        
        if( EmptyTo == "TRUE"){
          Empty.Blocks = Empty.Blocks[-which(Empty.Blocks==holder.to)]
        }
        
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%%          If the FROM block is Empty                   %%#   
        
        
        if( EmptyFrom == 1){
          
          #If the from block is empty       
          
          Empty.Blocks = union(Empty.Blocks,holder.from)
          #NotOnlySeeds = setdiff(1:K,c(SeedsOnly,Empty.Blocks))
          
          # Now, the correct rows should have been added/deleted. 
          # Update the block membership. 
          
          
          record.holder     = Reupdate_RIB( File1_data, File2_data,RIB,n.RIB, NotOnlySeeds )
          #record.holder     = update_RIB(File1_data,File2_data,K)
          n.RIB             = record.holder$n.RIB
          RIB               = record.holder$RIB
          rm(record.holder)
          
          max.RIB = apply( n.RIB, 2, max)
          max.RIB = max.RIB - count.matches
          non.seed.RIB = max.RIB
          PermSampling = ifelse( max.RIB > 1, ifelse(max.RIB < threshold, 1, 2), 0 )
          #SeedsOnly    = which(max.RIB==0)
          
          F1.recordsT = RIB[[1]][[holder.to]]
          F2.recordsT = RIB[[2]][[holder.to]]
          F1.recordsF = RIB[[1]][[holder.from]]
          F2.recordsF = RIB[[2]][[holder.from]]
          
          #For the records in the to block, update the block rows. 
          
          File1_data[,"BlockRow"][ F1.recordsT ] = orderC.1T
          File2_data[,"BlockRow"][ F2.recordsT ] = orderC.2T        
          
          lambda.order = File1_data[,"lambda"][ F1.recordsT ]
          lambda.order = lambda.order[ order(orderC.1T)]
          #lambda.order = lambda.order[ orderC.1T ]
          File2_data[,"lambda"][ F2.recordsT[order(orderC.2T)]]  = lambda.order
          #File2_data[,"lambda"][ F2.recordsT[ orderC.2T ] ]  = lambda.order
          
          # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
          #%%%          If the FROM block is NOT Empty               %%#  
          
        } else{ 
          
          # Now, the correct rows should have been added/deleted. Update the block membership.
          
          if( Delete.Counter > 0){
            # If we have deleted ANY records, run the full RIB
            
            #NotOnlySeeds = setdiff(1:K,c(SeedsOnly,Empty.Blocks))
            record.holder     = Reupdate_RIB( File1_data, File2_data,RIB, n.RIB, NotOnlySeeds )
            #record.holder     = update_RIB(File1_data,File2_data,K)
            n.RIB             = record.holder$n.RIB
            RIB               = record.holder$RIB
            rm(record.holder)
            
            max.RIB = apply( n.RIB, 2, max)
            max.RIB = max.RIB - count.matches
            non.seed.RIB = max.RIB
            PermSampling = ifelse( max.RIB > 1, ifelse(max.RIB < threshold, 1, 2), 0 )
            #SeedsOnly    = which(max.RIB==0)
            
            
          } else{ 
            
            for(k in c(holder.to,holder.from)){
              
              # First Pass through: Just the Seeds     
              data = File1_data[,"Block"]
              seedcheck = File1_data[,"Seed"]
              holder         = which(data==k & seedcheck==1 )
              RIB[[1]][[k]]  = holder
              holder         = which(data==k & seedcheck != 1 )
              RIB[[1]][[k]]  = c(RIB[[1]][[k]],holder) 
              n.RIB[1,k]     = length(RIB[[1]][[k]])
              data = File2_data[,"Block"]
              seedcheck = File2_data[,"Seed"]
              holder         = which(data==k & seedcheck==1 )
              RIB[[2]][[k]]  = holder
              holder         = which(data==k & seedcheck !=1 )
              RIB[[2]][[k]]  = c(RIB[[2]][[k]],holder) 
              n.RIB[2,k]     = length(RIB[[2]][[k]])
            }
            
            rm(holder) 
          }
          
          
          F1.recordsT = RIB[[1]][[holder.to]]
          F2.recordsT = RIB[[2]][[holder.to]] 
          
          File1_data[,"BlockRow"][F1.recordsT]   = orderC.1T
          File2_data[,"BlockRow"][F2.recordsT]   = orderC.2T
          
          #File1_data[,"lambda"]= 1:nrow.F1 
          
          lambda.order = File1_data[,"lambda"][ F1.recordsT ]
          lambda.order = lambda.order[ order(orderC.1T)]
          File2_data[,"lambda"][F2.recordsT[order(orderC.2T)]]= lambda.order
          
          F1.recordsF = RIB[[1]][[holder.from]]
          F2.recordsF = RIB[[2]][[holder.from]]    
          
          # For the records in the from block, update the block rows. 
          File1_data[,"BlockRow"][F1.recordsF]   = orderC.1F 
          File2_data[,"BlockRow"][F2.recordsF]   = orderC.2F        
          
          
          lambda.order = File1_data[,"lambda"][ F1.recordsF ]
          lambda.order = lambda.order[ order(orderC.1F)]
          File2_data[,"lambda"][ F2.recordsF[order(orderC.2F)] ]= lambda.order     
          
        }
        
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%% If any deletions, we need to rescale the lambdas     %%# 
        
        if(length(unique(File2_data[,"lambda"]))!= nrow.F2){
            stop('We have accepted a move and during the process somehow assigned two of the same individual labels to a record.')
        }
        
        if( Delete.Counter > 0 ){
          File1_data[,"lambda"] = order( order( File1_data[,"lambda"]))
          File2_data[,"lambda"] = order( order( File2_data[,"lambda"]))
        }
        
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%%                Sanity Checks                          %%# 
        
        if( sum(File2_data[F2.recordsT,"lambda"] %in% File1_data[F1.recordsT,"lambda"]) != length(F2.recordsT) ){
            stop('We have accepted a move, but within the TO block more than record has been assigned the same individual.' )
        }
        
        if( EmptyFrom != 1 & sum(File2_data[F2.recordsF,"lambda"] %in% File1_data[F1.recordsF,"lambda"]) != length(F2.recordsF) ){
            stop('We have accepted a move, but within the FROM block more than record has been assigned the same individual.' )
        }
        
        
        if( length(unique(File1_data[,"lambda"])) != length(unique(File2_data[,"lambda"]))){
          stop('We have accepted a move, but more than record has been assigned the same individual.' )
        }
        
        BRF1 = File1_data[,"BlockRow"]
        BRF2 = File2_data[,"BlockRow"]
        part1 = RIB[[1]][[holder.to]]
        part2 = RIB[[2]][[holder.to]]
        part1 = BRF1[part1]
        part2 = BRF2[part2]
        part3 = RIB[[1]][[holder.from]]
        part4 = RIB[[2]][[holder.from]]
        part3 = BRF1[part3]
        part4 = BRF2[part4]   
        # This is where we got stuck!!
        if( length(unique(part1)) != length(unique(part2)) | length(unique(part3)) != length(unique(part4))){
          stop('We are stopped at the parts check.')
          print("STUCK HERE")
          print(part1)
          print(part2)
          print(part3)
          print(part4)
          File1_data = File1_C
          File2_data = File2_C
          record.holder     = Reupdate_RIB( File1_data, File2_data ,RIB,n.RIB,NotOnlySeeds)
          n.RIB             = record.holder$n.RIB
          RIB               = record.holder$RIB
          break
        } 
        
        ################################
        #    Update Inputed.Block      #
        ################################
        
        #NotOnlySeeds = setdiff(1:K,c(SeedsOnly,Empty.Blocks))
        
        for( k in 1:K){
          # File 1
          tester = (RIB[[1]][[k]] > n1)
          sum.tester = sum(tester)
          
          Imputed.Block[k,1] = sum.tester 
          
          if( sum.tester > 0 ){
            Imputed.Records[[k]][[1]] = which( RIB[[1]][[k]] > n1)
          }
          # File 2 
          tester = (RIB[[2]][[k]] > n2)
          sum.tester = sum(tester)
          
          Imputed.Block[k,2] = sum.tester 
          
          if( sum.tester > 0 ){
            Imputed.Records[[k]][[2]] = which( RIB[[2]][[k]] > n2)
          }        
        }
        
        ################################
        #       Update Inputs          #
        ################################
        
        Y1          = File1_data[,col.Y1]
        Y2.in       = File2_data[,Y2.name]
        BlockRow.F1 = File1_data[,"BlockRow"]
        BlockRow    = File2_data[,"BlockRow"]
        D.Holder    = File2_data[,where.coefs.primary]
        F1.Lambda   = File1_data[,"lambda"]
        F2.Lambda   = File2_data[,"lambda"]
        
        nrow.F1              = nrow(File1_data) 
        nrow.F2              = nrow(File2_data) 
        rownames(File1_data) = 1:nrow.F1
        rownames(File2_data) = 1:nrow.F2
        
      } 
      
      
      if( sum(File2_data[,Y2.name] ==0) > 0 ){
        stop("During an acceptance step, one of the Y2 values is 0.")
      }
      
      
      if( EmptyFrom != 1){
        rm( Design.From)
        rm(orderC.2F) 
        rm(from.order)
      }
      #if( EmptyTo !="TRUE"){
      rm( Design.To )
      rm(to.order)
      rm(orderC.2T) 
      #}
      rm(Deleted)
      rm(Added) 
      
      
    } 
    
    record.holder     = Reupdate_RIB( File1_data, File2_data,RIB,n.RIB,NotOnlySeeds )
    #record.holder     = update_RIB(File1_data,File2_data,K)
    n.RIB             = record.holder$n.RIB
    RIB               = record.holder$RIB
    rm(record.holder) 
    
    max.RIB = apply( n.RIB, 2, max)
    max.RIB = max.RIB - count.matches
    PermSampling = ifelse( max.RIB > 1, ifelse(max.RIB < threshold, 1, 2), 0 )
    #SeedsOnly    = which(max.RIB==0)
    
    #################################
    ###  Update Gamma and A.probs ###
    #################################
    
    for( j in labelsJ_E ){
      if( error.files != 2){
        count = sum(E1[,j])
        zero.count = n.error1 -count 
        gamma.samp[1,j] = rbeta( 1, a.beta.prior[1,j] + count, b.beta.prior[1,j] + zero.count) 
      }
      if( error.files !=1){
        count = sum(E2[nonseed.F2,j])
        zero.count = n.error2 -count
        gamma.samp[2,j] = rbeta( 1, a.beta.prior[2,j] + count, b.beta.prior[2,j] + zero.count)  
      }
    }    
    
    #cat("The acceptance rate is", (Accept.Count)/length(Move.Holder) ,"for", length(Move.Holder)," proposed moves \n")  
    
    ##   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ##
    #### Current Blocks: Sample Theta and Reimpute ####
    
    imputed.F1 = n1:nrow.F1
    n.imputed.F1 = length(imputed.F1)
    
    imputed.F2 = n2:nrow.F2
    n.imputed.F2 = length(imputed.F2)
    
    # Reorder based on the sampled C
    Y2.in             = File2_data[,col.Y2]
    F2.Lambda         = File2_data[,"lambda"]
    Y1                = File1_data[,col.Y1]
    Y2                = Y2.in[order(F2.Lambda)]
    V.Holder          = File2_data[order(F2.Lambda),where.coefs.primary]
    V.Holder          = cbind(1,V.Holder)
    Sigma1.Part       = (nrow.F1 - p )/2
    Sigma2.Part       = (nrow.F2 - p2 )/2
    New.Ordering      = F2.Lambda
    
    V2 = cbind( 1, File2_data[, where.coefs.secondary ] )
    V2 = as.matrix(V2)
    
    V2TV2        = crossprod(V2)
    Gen          = zapsmall( eigen( V2TV2 )$values )
    Gen          = ifelse( sum(Gen==0) > 0 , "TRUE","FALSE")
    invert.part2  = invert_XTX( V2TV2 , Gen )
    
    t.V2         = t(V2)
    Proj2        = invert.part2%*%t.V2
    rm(V2TV2,Gen,t.V2)
    
    holder = Gupdate_and_Impute(Y1,Y2, Y2.in, V.Holder, New.Ordering,Sigma1.Part,Sigma2.Part,which.Y1,which.Y2, imputed.F2, imputed.F1,n.imputed.F2,n.imputed.F1,secondary.option,Proj2,V2,invert.part2)
    
    #Update based on the current imputations 
    Y1          = holder$Y1
    File1_data[imputed.F1,col.Y1] = Y1[imputed.F1] 
    Y2.in      = holder$Y2 # This is set back to the original ordering, but with updated imputations
    File2_data[imputed.F2,col.Y2] = Y2.in[imputed.F2] 
    beta.star  = holder$beta
    sigmasq    = holder$sigmasq
    sigma      = sqrt(sigmasq)
    eta        = holder$eta
    sigma2sq   = holder$sigma2sq
    theta      = holder$theta
    sqrt.var.Y2.impute = holder$impVar
    var.Y2impute = sqrt.var.Y2impute^2
    rm(holder)
    
    #D.Holder[,col.Y2] = Y2.in    
    
    #cat( "The beta draw is",beta.star,", the sigma draw is", sigma, "and the eta draw is", eta, "with variance", sigma2sq,"\n")    
    
    #####################################
    ###            STORAGE            ###
    #####################################
    
    if( s > burnin & s%%thinning==0){
      write( c(beta.star,sigmasq), file = paste(namepart,"ThetaOut.txt", sep=""), append = TRUE)
      write( c(eta,sigma2sq), file = paste(namepart,"EtaOut.txt", sep=""), append = TRUE)
      write( c(gamma.samp[error.file,]), file = paste(namepart,"GammaOut.txt", sep=""), append = TRUE)
      write( F2.Lambda[1:n2], file = paste(namepart,"LambdaOut.txt", sep=""), append = TRUE)
      write( Accept.Count, file = paste(namepart,"AcceptOut.txt", sep=""), append = TRUE)
      write( File2_data[1:n2,"Block"], file = paste(namepart,"AcceptOut.txt", sep=""), append = TRUE)
      # If we are interested in the completed data sets 
      if( complete.out == "yes"){
        write( F1.Lambda, file = paste(namepart,"F1Out.txt", sep=""), append = TRUE)
        write( Y2.in[order(F2.Lambda)], file = paste(namepart,"F2Out.txt", sep=""), append = TRUE)
      }
      if(error.files != 1 ){
        write( E2[,labelsJ_E], file = paste( namepart, "E2Out.txt",sep=""),append=TRUE)
      }
      if(error.files != 2 ){
        write( E1[,labelsJ_E], file = paste( namepart, "E1Out.txt",sep=""),append=TRUE)
      }
    }
    
  }
}



#%%%%%%%%%%%%%%%%%%%% #

# Purpose: Identifies whether a block requires a sampled permutation. If, for instance, there is one only non-seeded pair in the block, the permtuation is automatic

# Returns
# 1 if we do need to sample a permutation, 2 if not a matrix, and 0 otherwise

needs_C <-function(input,count.matches,input.type,n.RIB){
  
  if( input.type == "label"){
    block.count = n.RIB[1,input]
    
    if( block.count > 1){ #If there is more than one record pair
      
      if( (block.count - count.matches[input]) > 1 ){ # And at least two are not seeds
        return(1)
      } else { 
        return(0)    
      }
    } else {
      return(0)
    }
  }
  
  if( input.type == "Adapt"){
    
    if( block.count > 1){ #If there is more than one record pair
      
      if( (block.count - count.matches[input]) > 1 ){ # And at least two are not seeds
        return(1)
      } else { 
        return(0)    
      }
    } else {
      return(0)
    }
  }
  
  if( input.type == "records"){
    n.covariates = n.covariates 
    where.label1 = n.covariates + 2
    which.seeded= 0 
    if( input[1] !=0){
      which.seeded = intersect( seeds[,1], input[,where.label1])
    } 
    
    
    if( is.vector(input) == TRUE ){
      return( 0 )
    } else if( nrow(input) <= (length(which.seeded) +1) ){
      return( 0 )
    } else if(nrow(input)==1){
      return( 0 )
    } else if( input[1]==0){
      return(0)
    } else if(class(input) != "matrix"){
      return(2)
    } else{
      return(1)
    }
  }
  
}

#%%%%%%%%%%%%%%%%%%%% #

# Purpose: 
# Updates the regression parameters and re-imputes missing values 

# Returns: 
# list( "Y1" = Y1, "Y2" = Y2, "beta" = beta, "sigmasq" = sigmasq, "sigma" = sigma, "mu2" = mu2, "sigma2sq"=sigma2sq, "theta" =theta.small, "impVar" = sqrt(var.Y2impute))
# The re-imputed data and the regression parameters 

update_and_Impute <-function(Y1, Y2, Y2.in, Y2.info, New.Ordering, Sigma1.Part,Sigma2.Part,which.Y1,which.Y2,nrow.F1,nrow.F2){  
  
  imputed.F1 = n1:nrow.F1
  n.imputed.F1 = length(imputed.F1)
  
  imputed.F2 = n2:nrow.F2
  n.imputed.F2 = length(imputed.F2)
  
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ###########  Initialization  Theta #############
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###  
  
  # Create the model matrix 
  Vtheta = Y2.info 
  
  VthetaTVtheta= crossprod(Vtheta)
  Gen          = zapsmall( eigen( VthetaTVtheta )$values )
  Gen          = ifelse( sum(Gen==0) > 0 , "TRUE","FALSE")
  invert.part  = invert_XTX( VthetaTVtheta , Gen )
  
  t.Vtheta     = t(Vtheta)
  Proj         = invert.part%*%t.Vtheta
  hat.beta     = Proj%*%Y1
  
  ResidSS      = (Y1 - Vtheta%*%hat.beta) 
  ResidSS      = t(ResidSS)%*%ResidSS
  ResidSS      = .5* ResidSS
  
  sigmasq      = rigamma(1,Sigma1.Part,ResidSS) 
  sigma   = sqrt(sigmasq)
  #cat( "The variance is", sigmasq,"\n")
  var.beta     = sigmasq*invert.part
  #cat( "The primary variance is",diag(var.beta),"\n")
  #print( dim(var.beta) )
  beta         = mvrnorm( 1, hat.beta, var.beta ) 
  #cat( "The beta mean is",g1/(1+g1)*hat.beta + 1/(g1+1)*beta0,"\n")
  #cat( "The inital beta draw is",beta,"\n")
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ##########  Initialization : Y2 Coefs  ##########
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###  
  
  barY2    = mean(Y2)
  Var.Part = sum( (Y2 - barY2)^2 )
  
  sigma2sq = rigamma(1,Sigma2.Part,.5*Var.Part)
  mu2      = rnorm( 1, barY2, (1/n2)*sigma2sq) 
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ####### Initialization : Y2 Imputations #########
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  # Which Y1 values are matched to the imputed Y2 values?
  F1.matched.to.imputations = New.Ordering[imputed.F2]
  # Store the Y1 associated with these values
  Y1.matched.to.imputations = Y1[F1.matched.to.imputations]
  # Store the parts of the design matrix which go along with these values. 
  V1 = Vtheta[ F1.matched.to.imputations, -which.Y2]
  
  beta.impute = beta[ -which.Y2 ] 
  theta.small = beta[ which.Y2 ] 
  
  mean.Y2impute = V1%*%as.matrix(beta.impute) -  Y1.matched.to.imputations
  mean.Y2impute = sigmasq^(-1)*theta.small*mean.Y2impute
  mean.Y2impute = sigma2sq^(-1)*mu2 - mean.Y2impute
  
  var.Y2impute  = ( sigma2sq^(-1) + sigmasq^(-1)*(theta.small^2) )^(-1)
  mean.Y2impute = var.Y2impute*mean.Y2impute 
  #cat( "The imputation variance is", var.Y2impute,"\n")
  #print(summary(mean.Y2impute))
  imputed.all = rnorm( n.imputed.F2, mean.Y2impute, sqrt(var.Y2impute) )
  
  # Go to the original Y2 data (input ordering) and update the imputations 
  Y2.in[ imputed.F2 ]            = imputed.all 
  # Go to where the imputations are in the model matrix and impute
  Vtheta[F1.matched.to.imputations,which.Y2] = imputed.all 
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ####### Initialization : Y1 Imputations #########
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  Vtheta         = Vtheta[imputed.F1,]
  
  imputed.all    = rnorm( n.imputed.F1, Vtheta %*% as.matrix(beta), sqrt(sigmasq) ) 
  
  Y1[imputed.F1] = imputed.all  
  
  #################################################
  ## Step 6: Store and Release the Output 
  #################################################
  
  newlist <-list( "Y1" = Y1, "Y2" = Y2.in, "beta" = beta, "sigmasq" = sigmasq, "sigma" = sigma, "mu2" = mu2, "sigma2sq"=sigma2sq, "theta" =theta.small, "impVar" = sqrt(var.Y2impute))
  return( newlist ) 
  
}

#%%%%%%%%%%%%%%%%%%%% #

# Purpose: Updates two lists RIB[[1]] and RIB[[2]] which contain, respectively, a list of lists telling us which records from file 1 and file 2 respecitvely are in a given block. For instance, RIB[[1]][[24]] returns a vector of which records in File 1 belong to block 24

# Returns
# RIB[[1]] and RIB[[2]] 
# Also returns n.RIB, a matrix which counts the number of records from each file (row) in each block (column)

# THIS IS A SLOW FUNCTION

update_RIB<-function(File1, File2,K){
  n.RIB        = matrix( NA, nrow = 2, ncol = K )
  RIB          = vector( "list", length = 2 )
  RIB[[1]]     = vector( "list", length = K )
  RIB[[2]]     = vector( "list", length = K )
  
  for(k in 1:K){
    
    # First Pass through: Just the Seeds     
    data = File1[,"Block"]
    seedcheck = File1[,"Seed"]
    holder         = which(data==k & seedcheck==1 )
    RIB[[1]][[k]]  = holder
    holder         = which(data==k & seedcheck != 1 )
    RIB[[1]][[k]]  = c(RIB[[1]][[k]],holder) 
    n.RIB[1,k]     = length(RIB[[1]][[k]])
    data = File2[,"Block"]
    seedcheck = File2[,"Seed"]
    holder         = which(data==k & seedcheck==1 )
    RIB[[2]][[k]]  = holder
    holder         = which(data==k & seedcheck !=1 )
    RIB[[2]][[k]]  = c(RIB[[2]][[k]],holder) 
    n.RIB[2,k]     = length(RIB[[2]][[k]])
  }
  new.list <- list( "RIB" = RIB, "n.RIB" = n.RIB )
  return(new.list)
}

#%%%%%%%%%%%%%%%%%%%% #

# Purpose: 
# Exactly samples a new permutation for small blocks 

# Returns: 
# list( "order.out" = order.out, "Sum" = Sum.Part) 
# order.out = a vector giving the new ordering for the File 2 data in the block 
# Sum = a real number sum needed for the larger acceptance ratio

exactSampleC <-function( yData, orderC, Design.Matrix,num.pairs, beta, sigma,block,num.matched ){  

  if( num.matched > 0 ){
    # Remove the Type 1 Seeds
    switch.options = orderC[orderC > num.matched ]
    num.switch = num.pairs - num.matched 
    yData      = yData[-c(1:num.matched)] 
    order.out  = 1:num.matched
  } else{
    switch.options = orderC
    num.switch = num.pairs 
    order.out <- c()
  } 
  
  AllPossibleC = permutations( n =num.switch, r = num.switch, v= switch.options )
  
  Num.PossibleC = factorial(num.switch)
  
  likelihood.storage = rep(0, Num.PossibleC )   
  
  Variance.Component = sigma^2*diag(num.switch) 
  
  for( r in 1:Num.PossibleC ){ 
    orderC   = AllPossibleC[r,] 
    Current  = Design.Matrix[orderC,]
    likelihood.storage[r] = dmvnorm( yData, Current%*%beta, Variance.Component,log = FALSE)
    likelihood.storage[r] = ifelse( likelihood.storage[r] == 0, realmin, likelihood.storage[r])
  }
  
  Sum.Part = sum(likelihood.storage)
  C.Probabilities = likelihood.storage/Sum.Part
  SampledC = sample( Num.PossibleC, 1, replace= F, prob = C.Probabilities ) 
  
  order.out = c(order.out, AllPossibleC[ SampledC, ] ) 
  
  A.ratio.contribution = log( likelihood.storage[SampledC] )
  
  Sum.Part = log( Sum.Part )
  
  new.list <-list( "order.out" = order.out, "Sum" = Sum.Part) 
  
  return( new.list ) 
  
}

#%%%%%%%%%%%%%%%%%%%% #

# Purpose: 
# Performs the Gutman MH step - ie the block is large and it selects two records from File 2 and switches their position (repeats this step reps times)

# Returns: 
# order.out - a vector containing the selected order 

GutmanMH <-function( orderC,yData, Design.Matrix,num.pairs, beta, sqrt.sigma, reps, block,count.matches ){
  
  num.matched = count.matches[block]
  
  if( num.matched > 0 ){
    switch.options = orderC[-c(1:num.matched) ] 
  } else{
    switch.options = orderC
  } 
  
  chosen   = sample( switch.options, 2, replace= F )
  Current  = Design.Matrix[chosen,]
  Variance.Component = sqrt.sigma*diag(2)
  accepted = 0 
  
  for( r in 1:reps){ 
    
    Proposed = Design.Matrix[chosen[2:1],]
    
    like.part         = dmvnorm( yData[chosen], Proposed%*%beta, Variance.Component, log = TRUE)
    final.A.ratio     = like.part- dmvnorm( yData[chosen], Current%*%beta, Variance.Component, log = TRUE )
    #final.A.ratio     = log(like.part)
    
    if( log(runif(1)) < final.A.ratio ){
      Current = Proposed 
      orderC[chosen] = orderC[ chosen[2:1] ] 
      Design.Matrix[chosen,]= Design.Matrix[chosen[2:1],]
      accepted = 1 
    }
    
    chosen   = sample( switch.options, 2, replace= F )
    
  }
  
  order.out = orderC
  
  return(order.out)   
  
}

#%%%%%%%%%%%%%%%%%%%% #

# Purpose: 
# Performs the Gutman MH step within the new model- ie if the block is small, it exactly samples a permutation; if the block is large and it selects two records from File 2 and switches their position (repeats this step reps times)

# Returns: 
# File2_data, organized according to the new ordering 

A.sample.perm <- function(F1.holder,F2.holder, pairs, order.in, beta,sqrt.sigma,S.ind,reps,block,num.matched){  
  
  # If we get here, we know that we need to sample a permutation
  Design.Holder = F2.holder
  
  # S.ind takes a 1 if we use exact sampling, 2 if we use MH
  
  if( S.ind ==1 ){
    # Here we use exact sampling 
    order.out  = exactSampleC( F1.holder, 1:pairs, Design.Holder, pairs, beta, sqrt.sigma, block,num.matched)
    order.out  = order.out$order.out
  } else{ 
    
    # Order the Design Holder to match the original ordering 
    Design.Holder = Design.Holder[order.in,]
    
    orderC = 1:pairs
    
    if( num.matched > 0 ){
      switch.options = orderC[-c(1:num.matched) ] 
    } else{
      switch.options = orderC
    } 
    
    # Fix the variance component for the regression 
    Variance.Component = sqrt.sigma*diag(2)
    
    for( r in 1:reps){ 
      
      # Chose two non-type 1 seeds to swtich 
      chosen   = sample( switch.options, 2, replace= F )
      
      Current   = Design.Holder[chosen,]
      Proposed  = Design.Holder[chosen[2:1],]
      like.part = dmvnorm( F1.holder[chosen], Proposed%*%beta, Variance.Component, log = TRUE )
      like.part = like.part- dmvnorm( F1.holder[chosen], Current%*%beta, Variance.Component, log = TRUE)
      
      if( log(runif(1)) < like.part ){
        # if we accept, record the new ordering 
        orderC[chosen] = orderC[ chosen[2:1] ] 
        Design.Holder[chosen,]= Design.Holder[chosen[2:1],]
      }      
      
    }
    
    # Return the file ordering 
    order.out = orderC  
    
  }
  
  return(order.out) 
  
}


# A.sample.perm <-function(File1_data,File2_data,beta, sigma,RIB, n.RIB, PermSampling,reps,threshold,block,Y2.name){  
#   
#   F1.records = RIB[[1]][[block]] 
#   F2.records = RIB[[2]][[block]]
#   
#   pairs      = length(F1.records)
#   if(pairs==1){
#     break
#   }
#   
#   sqrt.sigma = sigma
# 
#   Design.Holder = File1_data[F1.records,]
#   Y1.part       = Design.Holder[,col.Y1]
#   Design.Matrix = Design.Holder[,where.coefs[-1]]
#   F1.order      = order( Design.Holder[,"BlockRow"] )
#   Design.Matrix = Design.Matrix[ F1.order, ]
#   Y1.part       = Y1.part[F1.order] 
#   Design.Holder = File2_data[F2.records,]
#   lambda.holder = Design.Holder[,"lambda"]
#   F2.order      = order( Design.Holder[,"BlockRow"] )
#   Design.Holder = Design.Holder[ F2.order, Y2.name ] 
#   Design.Matrix= cbind(1,Design.Holder, Design.Matrix ) 
#   lambda.holder = lambda.holder[F2.order]
#   
#   F1.records   = F1.records[ F1.order ] 
#   F2.records   = F2.records[ F2.order ] 
#     
#     if( PermSampling[block] ==1 ){
#       # Here we use exact sampling 
#       holder       = exactSampleC( Y1.part, 1:pairs, Design.Matrix, pairs, beta, sqrt.sigma, block, count.matches)
#       order.out    = holder$order.out
#       
#     } else{ 
#       
#       orderC = 1:pairs
#       
#       num.matched = count.matches[block]
#       
#       if( num.matched > 0 ){
#         switch.options = orderC[-c(1:num.matched) ] 
#       } else{
#         switch.options = orderC
#       } 
#       
#       chosen   = sample( switch.options, 2, replace= F )
#       Variance.Component = sqrt.sigma*diag(2)
#       accepted = 0 
#       
#       for( r in 1:reps){ 
#         
#         Current  = Design.Matrix[chosen,]
#         Proposed = Design.Matrix[chosen[2:1],]
#         #print(Proposed)
#         like.part         = dmvnorm( Y1.part[chosen], Proposed%*%beta, Variance.Component, log = TRUE )
#         like.part         = like.part- dmvnorm( Y1.part[chosen], Current%*%beta, Variance.Component, log = TRUE)
#         
#         if( log(runif(1)) < like.part ){
#           orderC[chosen] = orderC[ chosen[2:1] ] 
#           Design.Matrix[chosen,]= Design.Matrix[chosen[2:1],]
#           accepted = 1 
#         }
#         
#         chosen   = sample( switch.options, 2, replace= F )
#         
#       }
#       
#       order.out = orderC  
#     }
#     
#     File2_data[,"BlockRow"][ F2.records ]= order.out
#     File2_data[,"lambda"][ F2.records ]  = lambda.holder[ order.out ]
#     
#   return(File2_data) 
#   
# }


#%%%%%%%%%%%%%%%%%%%% #

# Purpose: 
# Same as update_RIB, but only for the blocks that require re-updating. Blocks with all seeds, for instance, do not need us to update the in-block records 

# Returns: 
# tow lists RIB[[1]] and RIB[[2]] which contain, respectively, a list of lists telling us which records from file 1 and file 2 respecitvely are in a given block. For instance, RIB[[1]][[24]] returns a vector of which records in File 1 belong to block 24
# n.RIB, a matrix which counts the number of records from each file (row) in each block (column)

Reupdate_RIB<-function(File1, File2,RIB,n.RIB,Possible){
  FileBlock <- list(c(File1[,"Block"]),File2[,"Block"])
  FileSeed <- list(c(File1[,"Seed"]),File2[,"Seed"])
  
  Possible <- sort.int(Possible)
  for (f in 1:2) {
    # First Pass through: Just the Seeds
    Index <- which(FileSeed[[f]]==1)
    Data <- FileBlock[[f]][Index]
    temp <- sort.int(Data, index.return=TRUE)
    Data1 <- temp[[1]]
    Index1 <- Index[temp[[2]]]
    
    Index <- which(FileSeed[[f]]!=1)
    Data <- FileBlock[[f]][Index]
    temp <- sort.int(Data, index.return=TRUE)
    Data0 <- temp[[1]]
    Index0 <- Index[temp[[2]]]
    
    for(k in Possible){
      RIB[[f]][[k]]  <- c(Index1[which(Data1 ==k)], Index0[which(Data0
                                                                 ==k)])
    }
    n.RIB[f,]     <- sapply(RIB[[f]],length)
  }
  new.list <- list( "RIB" = RIB, "n.RIB" = n.RIB )
  return(new.list)
}


