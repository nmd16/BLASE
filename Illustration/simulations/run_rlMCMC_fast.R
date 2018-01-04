run_rlMCMC_fast<-function( MCMCSeed, File1, File2, method, its, burnin, DPits, DPburnin,gap, reps, thinning, DPthinning, printgap, Comment, coefs.primary,coefs.primary.all,coefs.secondary,col.Y1, col.Y2, secondary.option, secondary.option2, ObsBlock_restriction, BlockSeed_restriction, error.file, error.fields, blockingVars,matches, a.beta.prior, b.beta.prior, Hstar, DPOffset,Y1.name,Y2.name,namepart,threshold,complete.out,type2seeds,blocksIn="FALSE",blocks=NULL,blocks1=NULL,blocks2=NULL,blockSizeIn= "FALSE", RIB = NULL, n.RIB = NULL){
  
  #library(rlMCMC)
  
  realmin =  1e-08
  
  #load("C:/Users/Nicole/Desktop/ForPackage/JustData.RData")
  
  ## %%%%%%%%%%%%%%%%% ## 
  ##  Check the Inputs ##
  ## %%%%%%%%%%%%%%%%% ##
  
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
  
  if( class(MCMCSeed) != "integer"){
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
  #d[6] = 192
  
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
  
  # If we fed in the observed blocks, great. 
  
  if(blocksIn == "FALSE"){
  
    blocks <-NULL 
    
    for( i in 1:N){
      blocks[i] = category.block(B[i,],d,J)
    }
    
    blocks1 = blocks[1:n1]
    blocks2 = blocks[(n1+1):N]
    
    K = length(unique(blocks))
    
  }
  
  
  if( sum(matches !=0) > 0 ){
    if( sum(blocks1[type1seeds] == blocks2[type1seeds])!= n.type1){
      stop('Some declared Type 1 Seed pairs are not in the same block.')
    }
  }
  
  observed.blocks = sort(blocks)
  
  ## %%%%%%%%%%%%%%%%%%% # 
  ## Block Restriction  ##
  ## %%%%%%%%%%%%%%%%%%% #
  
  if( exists("ObsBlock_restriction") == FALSE ){
    warning('Block restriction not specified. By default, moves are restricted to observed blocks.')
    ObsBlock_restriction = "yes"
  } 
  
  if( ObsBlock_restriction == "yes" & blocksIn ==FALSE){
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
  
  if(blocksSizeIn = "FALSE"){
  
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
  
  }
  
  
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
  ## Format File1 and File 2 for the regression ##
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
  
  File1_data = File1
  File2_data = File2 
  
  # Now, we need to make sure that ALL of the coefficients which are in either model end up in File 1 and File2 ( changing this changes a lot of coding, and we can do it later, 
  #but not now)
  
  coefs.all = union(coefs.primary,coefs.secondary)
  if(modelPrimaryIn ==FALSE){
    f1 <- as.formula(paste("File1[,",col.Y1,"]~", paste(coefs.all, collapse="+")))
    modelPart <- model.matrix( lm(f1,data=File2))
  } else{
    modelPart <-model.matrix(truthModel)
  } 
  
  File1_data = cbind( File1[,col.Y1], modelPart[,-c(1:2)])
  colnames(File1_data)[col.Y1] = Y1.name
  File2_data = cbind( File2[,col.Y2], modelPart[,-c(1:2)])
  colnames(File2_data)[1] = Y2.name
  col.Y2 = 1 
  
  if(modelPrimaryIn ==FALSE){
    # Create the Primary Model Matrix 
    f1 <- as.formula(paste("File1[,",col.Y1,"]~", paste(coefs.primary.all, collapse="+")))
    modelPartPrimary <- model.matrix( lm(f1,data=File2))
    # Create the Secondary Model Matrix 
    f1 <- as.formula(paste("File2[,",col.Y2,"]~", paste(coefs.secondary, collapse="+")))
    modelPartSecondary <- model.matrix( lm(f1,data=File2))
    rm(f1)
  } else{ 
    modelPartPrimary = modelPart
    modelPartSecondary <-model.matrix( truthSecondary )
  }
  

  
  # In the future, I think I want to work with vectors instead of matrices, and just replace the vectors that we need with the reorderings. This might make the code more efficient, but for now, leave it alone. 
  
  # Now, I have to know which of these columns are used for which regressions, and this somehow needs to be automated and stored 
  
  coefs.all = union(coefs.primary.all,coefs.all)
  
  # Identify which of these are primary
  where.coefs.primary = which(coefs.all %in% coefs.primary.all)
  # To count the columns, consult the levels associated with each. 
  #where.coefs.primary = where.coefs.primary[1]:(where.coefs.primary[1]+sum(d[coefs.primary])-length(coefs.primary))
  where.coefs.primary = 1:18
  
  # Identify which of these are secondary
  where.coefs.secondary = which( coefs.all %in% coefs.secondary)
  # To count the columns, consult the levels associated with each. 
  #where.coefs.secondary = where.coefs.secondary[1]:(1+sum(d[coefs.secondary])-length(coefs.secondary))
  File2_data = cbind(File2_data, File2_data[,8]*File2_data[,2])
  colnames(File2_data)[length(colnames(File2_data))] = "YearIndTRUE:sexM"
  File1_data = cbind(File1_data, File1_data[,8]*File1_data[,2])
  colnames(File1_data)[length(colnames(File1_data))] = "YearIndTRUE:sexM"
  where.coefs.secondary = c(8,2:7,19,14:18)
  
  where.coefs = sort(union(where.coefs.primary,where.coefs.secondary))
  
  p  = length(where.coefs.primary)+1
  p2 = length(where.coefs.secondary)+1
  
  
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
  ## Store the Seed and Block Information ##
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
  
  BlockRow   = Seed = Block= 0
  File1_data = cbind( File1_data, BlockRow, Seed,Block)
  File2_data = cbind( File2_data, BlockRow, Seed, Block)
  
  if( SeedsIn == FALSE){
    # Exact Matches
    File1_data[type1seeds,"Seed"] = 1
    File2_data[type1seeds,"Seed"] = 1
    
    # Perfect Blocking 
    if( length(type2seeds) > 0 ){
      File1_data[type2seeds,"Seed"] = 2
      File2_data[type2seeds,"Seed"] = 2
    }
  } else{
    File1_data[,"Seed"] = Seeds
    File2_data[,"Seed"] = Seeds
  }
  
  if( method == "Extension"){
  
    ## Check for Type 2 Seeds
    #possibleSchool <- how_many_possible( B, N = 10000, d, field = 6,J, type1Seeds, Legal_Index,blocks)
    
    possibleEthnic <- how_many_possible( B, N = n1 , d, field = 4, J, type1Seeds, Legal_Index,blocks)
    
    ## Check for Type 2 Seeds
    #possible.length.S = unlist( lapply( possibleSchool, length ))
    #seeds.school = which( possible.length.S ==1)
    possible.length.E = unlist( lapply( possibleEthnic, length ))
    type2seeds = which(possible.length.E==1)
    #seeds.ethnic = which( possible.length.E ==1)
    #seeds.all = intersect( seeds.school, seeds.ethnic)
    type2seeds = setdiff( type2seeds, type1seeds)
    #seeds.ethnic = setdiff(seeds.ethnic,c(seeds.school,type1seeds))
    #seeds.school = setdiff(seeds.school,c(seeds.ethnic,type1seeds))
    # Any records for which possibleSchool = 0 is going to be declared a type 2 seeds
    #possible.length = unlist( lapply( possibleSchool, length ))
    #seeds.all = which( possible.length ==1)
    #type2seeds = setdiff( seeds.all, type1seeds)
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
    
  if( sum(matches !=0) > 0 ){
    if( sum(File1_data[type1seeds,"Block"] == File2_data[type1seeds,"Block"])!= n.type1){
      stop('Some declared Type 1 Seed pairs are not in the same block post block storage.')
    }
  }
  
    if( count.matchesIn == FALSE){
  
      count.matches   = rep(0,K)
      
      # How many matches in each block? 
      for( k in 1:K){
        count.matches[k] = sum(File1_data[type1seeds,"Block"]==k)
      }
      rm(k)
      } 
    } else{ 
      count.matches = n.RIB.T1    
      rm(n.RIB.T1)
    }
  
  ## %%%%%%%%%%%%%%%%% ##
  ##  Format the Data  ## 
  ## %%%%%%%%%%%%%%%%% ##
  
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
  
  #possible <- how_many_possible( B, N = 10000, d, field = 6,J, type1Seeds, Legal_Index,blocks)
  
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
  ## Save the formatted data sets  ## 
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
  
  save(File1_data,file ="File1Formatted.RData")
  save(File2_data,file ="File2Formatted.RData")
  save.image("PreRunStep.RData")
  
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
  ## UNCHANGED PAST THIS POINT ##
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
  ####      Pre Run Steps       ####
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
  
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
  
  
  # Load the output 
  #load(paste( "GOut", namepart, ".RData", sep= ""))
  #GonlyChain1 = G.outlist
  
  #save.image( paste(namepart, ".RData",sep=""))
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
  #### Run the CHECK Gutman MCMC Only ####
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
  
  if( method == "GTest"){
    # This still needs to be updated, but it's only going go be used for testing, so it's not a big deal 
    File1_True = File1_data[1:n1,]
    File2_True = File1_data[1:n1,]
    File2_True[,col.Y2] = File2_data[1:n1,col.Y2]
    
    record.holder  = G.records.in.block( File1_True, File2_True,K,type1seeds )
    n.RIB          = record.holder$n.RIB
    RIB            = record.holder$RIB
    rm(record.holder)
    
    PermSampling = ifelse( n.RIB[1,] > 1, ifelse(n.RIB[1,] < threshold, 1, 2), 0 )
    SeedsOnly    = which(n.RIB==0)
    needs.C      = which(PermSampling>0)
    
    holder     = G.get_starting_permutation( K, RIB, File1_True, File2_True,n1,n2 )
    File1_True = holder$File1
    File2_True = holder$File2
    rm(holder) 
    
    # Run the script 
    gutman_mcmc_Speed("no",MCMCSeed, its,burnin, thinning, gap, n.RIB, RIB, File1_True, File2_True,K,which.Y2, which.Y1,type1seeds, where.coefs.primary, where.coefs.secondary, p, p2, n.type1,needs.C,reps,PermSampling,secondary.option="dependent",n1,n2,col.Y1,col.Y2,paste(namepart,"Test",sep=""))
    
  }
  
  ## %%%%%%%%%%%%% ##
  ##   Extension   ##
  ## %%%%%%%%%%%%% ##
  
  if( method == "Extension"){
    
    record.holder  = records.in.block( File1_data, File2_data,K )
    n.RIB          = record.holder$n.RIB
    RIB            = record.holder$RIB
    rm(record.holder)
    
    max.RIB = apply( n.RIB, 2, max)
    max.RIB = max.RIB - count.matches
    
    PermSampling = ifelse( max.RIB > 1, ifelse(max.RIB < threshold, 1, 2), 0 )
    SeedsOnly    = which(max.RIB==0)
    # This vector takes a 1 if we use exact sampling, 2 if we use MH and 0 otherwise.
    needs.C = which(PermSampling > 0)
    
    holder  = get_starting_permutation( K, RIB, matches,File1_data, File2_data,nrow.F1,nrow.F2)
    File1_data = holder$File1
    File2_data = holder$File2
    rm(holder) 
    
    if( exists("BlockSeed_restriction")==FALSE){
      BlockSeed_restriction = "yes"
    }
    
    if( BlockSeed_restriction != "no"){
      if( length(SeedsOnly > 0) ){
        Legal_Index = Legal_Index[-SeedsOnly]
      }
    }
    
    # %%%%%%%%%%%%%%%%%%%%%%%%% ##
    ####  Initialize the DP   ####
    # %%%%%%%%%%%%%%%%%%%%%%%%% ##
    
    require(NPBayesImpute)
    
    for( j in 1:J){
      if(class(Bhat[,j]) != "factor"){
        Bhat[,j] = as.factor(Bhat[,j])
      }
    }
    #Bhat[,4] = factor(Bhat[,4])
    
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
    ####               Run the Main Script             ####
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #### 
    
    nonseed.F2 = MoveCand_F2
    
    
    if( error.file == 1){
      n.error1 = length(nonseed.F2)
      n.error2 = 0 
    }else{ 
      n.error2 = length(nonseed.F2)
      n.error1 = 0 
    }
    
    labelsJ_E = which( blockingVars %in% error.fields)
    
    seeds.school= NULL
    seeds.ethnic = NULL
    T2seeds1    = seeds.ethnic
    T2seeds2    = seeds.school     
    
 # %%%%%%%%%%%%%%%%%%%%%%%%% # 
  ##          Script         ##
  ## %%%%%%%%%%%%%%%%%%%%%%%  #
  
  # Set the Seed 
  set.seed(MCMCSeed)
  
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
    #gamma.samp[1,u] = rbeta( 1, a.beta.prior[1,u], b.beta.prior[1,u] )
    gamma.samp[error.file,u] = rbeta( 1, a.beta.prior[2,u], b.beta.prior[2,u] )
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
  
  if( nrow.F2 > n2){
    
    # Which Y1 values are matched to imputations?
    F1.matched.to.imputations = File2_data[imputed.F2,"lambda"]
    
    # Store the Y1 associated with these values
    Y1.matched.to.imputations = File1_data[F1.matched.to.imputations,col.Y1]
    # Store the parts of the design matrix which go along with these values. 
    V1 = cbind( 1, File1_data[F1.matched.to.imputations, where.coefs.primary[-1] ] )
    V1 = as.matrix(V1)  
    
    beta.impute = beta.star[ -which.Y2 ] 
    theta.small = beta.star[ which.Y2 ] 
    
    mean.Y2impute = V1%*%as.matrix(beta.impute) - Y1.matched.to.imputations
    mean.Y2impute = sigmasq^(-1)*theta.small*mean.Y2impute
    
    if( secondary.option != "independent"){
      V2.imp = cbind( 1, File2_data[imputed.F2, where.coefs.secondary ] )%*%eta
      
      mean.Y2impute = sigma2sq^(-1)*V2.imp - mean.Y2impute
    } else{
      mean.Y2impute = sigma2sq^(-1)*eta - mean.Y2impute
    }
    
    var.Y2impute  = ( sigma2sq^(-1) + sigmasq^(-1)*(theta.small^2) )^(-1)
    mean.Y2impute = var.Y2impute*mean.Y2impute 
    sqrt.var.Y2impute = sqrt(var.Y2impute)
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
  if( BlockSeed_restriction=="yes"){
    NotOnlySeeds = setdiff(1:K,c(SeedsOnly)) 
  } else{
    NotOnlySeeds = 1:K
  }
  
  max.RIB = apply( n.RIB, 2, max)
  max.RIB = max.RIB - count.matches
  non.seed.RIB = max.RIB
  
  Y1          = File1_data[,col.Y1]
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
  #################               Run New Sampler       ######################
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
  
  #require(NPBayesImpute)
  
  #for( j in 1:J){
  #  if(class(Bhat[,j]) != "factor"){
  #    Bhat[,j] = as.factor(Bhat[,j])
  #  }
  #}
  
  model <-CreateModel( Bhat, NULL, Hstar, 0, 0.25,0.25)
  model$Run(DPburnin,DPits,DPthinning)
  dpmpm <-model$snapshot
  
  nrow.F1              = nrow(File1_data) 
  nrow.F2              = nrow(File2_data) 
  rownames(File1_data) = 1:nrow.F1
  rownames(File2_data) = 1:nrow.F2
  
  
  for( s in 1:its){ 
    
    Accept.Out = 0 
    
    if(s%%printgap == 0){ print(s) }
    
    ######################
    ### Propose a new E ##
    ######################
    
    Estar = propose.E(E1,E2,labelsJ_E,nonseed.F1,nonseed.F2,gamma.samp,T2seeds1,T2seeds2)  
    
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
    
    
    move.holder =  propose.B( J, d, File1_data[,"Block"], File2_data[,"Block"], Estar, Move.F1, Move.F2, current.B, Bhat1, Bhat2, Blocking_WithError, Index_to_Block, Legal_Index,dpmpm$psi,dpmpm$z+1, File2_data[,"lambda"], File1_data[,"lambda"],gamma.samp,A.probs,BlockSeed_restriction="no",SeedsOnly,possibleEthnic,possibleSchool )
    
    ###############################
    ## Accept/Reject Block Moves ##
    ###############################    
  
    beta.impute    = beta.star[-2]
    theta.small    = beta.star[2]
  
    Move.Holder = move.holder$move.holder    
  
    Accept.Count = 0     
    
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%###
    ####     Begin the Loop     ####
    
    # Step 1: Grab the move holder
    
    DeletedSet1 <- NULL
    DeletedSet2 <-NULL
    ptm <-proc.time()
    for( u in 1:length(Move.Holder)){
      
        
      needtobreak= 0 
      
      if( Comment == "YES"){
        print(u)
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
  
      # The first element: which file we are working with (1 or 2)
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
          stop( "Error: The moving records is not in the from block;File1Move")
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
            stop( "Error: The moving record is not in the from block; File1 Move")
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
          stop( "Error: The moving records is not in the from block;File2Move")
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
        
         
            Design.From.Orig = Design.From.Orig[order(orderC.2F.orig),]
            
          } else{
            one.record.only= 1
            Design.From = File2_data[F2.records,][where.coefs.primary]
            Design.From = c(1,Design.From)          
            Design.From[where.impute] = Added$From = imputed.reg2           
          }         
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
      
      } 
      
      # Note that what we get out of all of this is accept (do we accept the move in E,B space) and the new permutation in C space (from.order, to.order) 
      
      ########################################################################################################
      # If we accept : Now how to we format our output properly?  
      ########################################################################################################
      
      
      if( accept == 1 ){        
        
        Delete.Counter = 0 # Count the number of imputed records we delete
        Accept.Count   = Accept.Count + 1 # Keeps a running total of accepted moved this iteration
        
        if(holder.to %in% SeedsOnly & BlockSeed_restriction=="yes"){
          stop('Post the acceptance step, we accepting a move to a seed block even though the seed restrition is on.')
          
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
          File1_data[holder.r,"Block"] = holder.to 
          E1[holder.r,] = Estar[[1]][holder.r,]
          RIB[[1]][[holder.to]] =  c(RIB[[1]][[holder.to]],holder.r)
          RIB[[1]][[holder.from]] = setdiff( RIB[[1]][[holder.to]],holder.r)
        } else{ 
          File2_data[ holder.r, where.bvars.F2]  = Imp.Holder2[[holder.to]][where.bvars.F2] 
          File2_data[holder.r,"Block"] = holder.to 
          E2[holder.r,] = Estar[[2]][holder.r,]
          if( sum(File2_data[,Y2.name] ==0) > 0 ){
            stop("During the first check, one of the Y2 values is 0.")
          }
          
          RIB[[2]][[holder.to]]   = c(RIB[[2]][[holder.to]],holder.r)
          RIB[[2]][[holder.from]] = setdiff( RIB[[2]][[holder.from]],holder.r)
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
            RIB[[1]][[holder.to]] = F1.recordsT
            
            # Update the imputation status
            Imputed.Block[holder.to,1] = Imputed.Block[holder.to,1] + 1 
            Imputed.Records[[holder.to]][[1]] = which(RIB[[1]][[holder.to]]>n1)
            
          } else{ 
            # Added a record for File 2 
            added.in          = Imp.Holder2[[holder.to]]
            added.in[col.Y2]  = Added$To 
            nrow.F2           = nrow.F2 + 1
            File2_data        = rbind( File2_data, added.in) 
            F2.recordsT = ifelse( F2.recordsT ==-99, nrow.F2, F2.recordsT)
            RIB[[2]][[holder.to]] = F2.recordsT
            
            # Update the imputation status
            Imputed.Block[holder.to,2] = Imputed.Block[holder.to,2] + 1 
            Imputed.Records[[holder.to]][[2]] = which(F2.recordsT>n2)
            
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
            RIB[[1]][[holder.from]] = F1.recordsF
            
            # Update the imputation status
            Imputed.Block[holder.from,1] = Imputed.Block[holder.from,1] + 1 
            Imputed.Records[[holder.from]][[1]] = which(F1.recordsF>n1)
            
          } else{ 
            added.in                    = Imp.Holder2[[holder.from]]
            added.in[col.Y2] = Added$From
            nrow.F2                     = nrow.F2 + 1 
            File2_data                   = rbind( File2_data, added.in)
            F2.recordsF = ifelse( F2.recordsF ==-99, nrow.F2, F2.recordsF)
            RIB[[2]][[holder.from]] = F2.recordsF
            
            # Update the imputation status
            Imputed.Block[holder.from,2] = Imputed.Block[holder.from,2] + 1 
            Imputed.Records[[holder.from]][[2]] = which(F2.recordsF>n2)
          }
        }  
        
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%%     If records were deleted from the from block       %%#  
        
        if( class(Deleted$From) != "NULL" ){ 
          
          Delete.Counter = Delete.Counter  + 1 
          
          if( AcceptOption[2] == 11 ){
            # We deleted records from File 1 
            DeletedSet1 = c(DeletedSet1, Deleted$From )
            RIB[[1]][[holder.from]] = setdiff( RIB[[1]][[holder.from]],Deleted$From)
            # Update the imputation status
          Imputed.Block[holder.from,1] = Imputed.Block[holder.from,1] - 1 
          Imputed.Records[[holder.from]][[1]]=which(RIB[[1]][[holder.from]]>n1)
          
            
          } else{ 
            # We deleted records from File 2 
            DeletedSet2 = c(DeletedSet2, Deleted$From )
            RIB[[2]][[holder.from]] = setdiff(RIB[[2]][[holder.from]],Deleted$From)
            # Update the imputation status
            Imputed.Block[holder.from,2] = Imputed.Block[holder.from,2] - 1 
            Imputed.Records[[holder.from]][[2]] = which(RIB[[2]][[holder.from]]>n2)                        
          } 
        }
        
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
        #%%%     If records were deleted from the TO block       %%#  
        
        
        if( class(Deleted$To) != "NULL" ){
          Delete.Counter = Delete.Counter  + 1 
          # Which file did we delete from? 
          if( AcceptOption[1] == 11 ){
            # If we deleted from File 1
            DeletedSet1 = c(DeletedSet1, Deleted$To )
            RIB[[1]][[holder.to]] = setdiff( RIB[[1]][[holder.to]], Deleted$To)        
            
            # Update the imputation status
            Imputed.Block[holder.to,1] = Imputed.Block[holder.to,1] - 1 
            Imputed.Records[[holder.to]][[1]] = which(RIB[[1]][[holder.to]]>n1)
            
          } else{ 
            # If we deleted from File 2 
            DeletedSet2 = c(DeletedSet2, Deleted$To )
            RIB[[2]][[holder.to]] = setdiff( RIB[[2]][[holder.to]], Deleted$To)
            
            # Update the imputation status
            Imputed.Block[holder.to,2] = Imputed.Block[holder.to,2] - 1 
            Imputed.Records[[holder.to]][[2]]= which(RIB[[2]][[holder.to]]>n2)
          } 
          
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
          # The block is now empty
          PermSampling[holder.from] = 0 
          #For the records in the to block, update the block rows.       
          File1_data[,"BlockRow"][ F1.recordsT ] = orderC.1T
          File2_data[,"BlockRow"][ F2.recordsT ] = orderC.2T  
          
          lambda.order = File1_data[,"lambda"][ F1.recordsT ]
          lambda.order = lambda.order[ order(orderC.1T)]
          File2_data[,"lambda"][ F2.recordsT[order(orderC.2T)]]  = lambda.order
          
          max.RIB[holder.from] = 0 
          F1.recordsF = F2.recordsF = NULL
          
          # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
          #%%%          If the FROM block is NOT Empty               %%#  
          
        } else{ 
          
          
          File1_data[,"BlockRow"][F1.recordsT]   = orderC.1T
          File2_data[,"BlockRow"][F2.recordsT]   = orderC.2T
          
          
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
          stop("STUCK HERE")
          
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
      rm( Design.To )
      rm(to.order)
      rm(orderC.2T) 
      rm(Deleted)
      rm(Added) 
    
      
    }
    proc.time() -ptm
    
    ################################
    #       Update Inputs          #
    ################################

    if( length(DeletedSet2) > 0 ){
      BlockRow    = BlockRow[-DeletedSet2]
      D.Holder    = D.Holder[-DeletedSet2,] 
      File2_data = File2_data[-DeletedSet2,]
    }
    if( length(DeletedSet1)){
      BlockRow.F1 = BlockRow.F1[-DeletedSet1]      
      F1.Lambda   = F1.Lambda[-DeletedSet1]
      File1_data = File1_data[-DeletedSet1,]
    }

    nrow.F1              = nrow(File1_data) 
    nrow.F2              = nrow(File2_data) 
    rownames(File1_data) = 1:nrow.F1
    rownames(File2_data) = 1:nrow.F2
    
    File1_data[,"lambda"] = order( order( File1_data[,"lambda"]))
    File2_data[,"lambda"] = order( order( File2_data[,"lambda"]))
    
    
    record.holder     = Reupdate_RIB( File1_data, File2_data,RIB,n.RIB,NotOnlySeeds)
    n.RIB             = record.holder$n.RIB
    RIB               = record.holder$RIB
    rm(record.holder) 
    
    max.RIB = apply( n.RIB, 2, max)
    max.RIB = max.RIB - count.matches
    PermSampling = ifelse( max.RIB > 1, ifelse(max.RIB < threshold, 1, 2), 0 )
    #SeedsOnly    = which(max.RIB==0)
    
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
    
    #################################
    ###  Update Gamma and A.probs ###
    #################################
    
    for( j in labelsJ_E ){
      if( error.file != 2){
        count = sum(E1[,j])
        zero.count = n.error1 -count 
        gamma.samp[1,j] = rbeta( 1, a.beta.prior[1,j] + count, b.beta.prior[1,j] + zero.count) 
      }
      if( error.file !=1){
        count = sum(E2[nonseed.F2,j])
        zero.count = n.error2 -count
        gamma.samp[2,j] = rbeta( 1, a.beta.prior[2,j] + count, b.beta.prior[2,j] + zero.count)  
      }
    }    
    

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
    
    holder = Gupdate_and_Impute(Y1,Y2, Y2.in, V.Holder, New.Ordering,Sigma1.Part,Sigma2.Part,which.Y1,which.Y2, imputed.F2, imputed.F1,n.imputed.F2,n.imputed.F1,secondary.option="dependent",Proj2,V2,invert.part2)
    
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
 
    
    #####################################
    ###            STORAGE            ###
    #####################################    
  
    if( s > burnin & s%%thinning==0){
      cat( c(round(beta.star,4),round(sigmasq,4),""), file = paste(namepart,"ThetaOut.txt", sep=""), fill=FALSE, append = TRUE)
      cat( c(round(eta,4),round(sigma2sq,4),""), file = paste(namepart,"EtaOut.txt", sep=""), fill=FALSE, append = TRUE)
      write( c(gamma.samp[error.file,c(4)]), file = paste(namepart,"GammaOut1.txt", sep=""), append = TRUE)
      #write( c(gamma.samp[error.file,c(6)]), file = paste(namepart,"GammaOut2.txt", sep=""), append = TRUE)
      Lambda.out <- sum(F2.Lambda[1:n2]==1:n2)
      write( Lambda.out, file = paste(namepart,"LambdaOut.txt", sep=""), append = TRUE)
      Lambda.out <- Lambda.out/nrow.F2
      write( Lambda.out, file = paste(namepart,"AllLambdaPerc.txt", sep=""), append = TRUE)
      write( nrow.F2, file = paste(namepart,"SizeOut.txt",sep=""),append=T)
      write( Accept.Count, file = paste(namepart,"AcceptOut.txt", sep=""), append = TRUE)
      Accept.Count = 0 
      cat( c(current.B[[2]][-type1seeds,4],""), file = paste(namepart,"EthAll.txt", sep=""), fill=F, append = T)
      # If we are interested in the completed data sets 
      if( s > 9000 & complete.out == "yes"){
        combined = cbind(  File1_data[,1], Y2.in[order(F2.Lambda)],"F2Lam"= F2.Lambda)
        write.table( combined, file = paste(namepart,"Data",s,".csv",sep=""),row.names = F, col.names= F, sep=",")
      }
      
    }
    
  }
  }
}

propose.B <- function(J, d, File1_dataBlock, File2_dataBlock, Estar, Move.F1, Move.F2, current.B, Bhat1, Bhat2, labelsJ_E, Index_to_Block, Legal_Index,phi,z,File2_dataLambda, File1_dataLambda, gamma.samp,A.probs,BlockSeed_restriction,SeedsOnly,possible1,possible2){
  
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
    if( which.changed == 4){
      which.possible = possible1[[r]]
    } else{ 
      which.possible = possible2[[r]]
    }
    field.options  = 1:d[which.changed]
    field.options  = field.options[ - observed.field ] 
    field.options  = intersect(which.possible,field.options)
    num.fields     = length( field.options )    

    if( num.fields >= 1){
      
      if( Epart[r, which.changed]==0 ){
      
        # Option 3: We move from an error to observed       
      
        B.holder[which.changed] = observed.field 
        
        # For the transition probability ratio, need prob moving to current from observed
        
        norm.probs  = initial.probs[field.options] 
        norm.probs  = norm.probs/ sum(norm.probs) 
        desired = which( which.possible == observed.field) 
        desired = which ( which.possible[-desired] == current.field)
        # This is the issue here. 
        result.prob = norm.probs[desired] 
        
        prior.part = initial.probs[observed.field]/ initial.probs[current.field]
        
      } else{ 
        
        # Options 1 and 2 : We move to an error prone field       
        
        if(  num.fields == 1){
          
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
        stop( paste("Record", r, "was proposed to be moved to an illegal block."))
      }
      
      # If this index belongs to a filled block, 
      if( ObsBlock_restriction == "no" | to.block %in% Legal_Index){
        to.block            = which(Index_to_Block==to.block)
        
        if(to.block %in% SeedsOnly & BlockSeed_restriction=="yes"){
          stop("Whoopse! We propose a seed move")
        }
        
        from.block          = File2_dataBlock[r]        
          
        if(  to.block != from.block ){  
          H = H + 1
          
          if( Epart[r,"Aratio"] ==1){
            # E(s) = E*  = 1 
            
            A.ratio = prior.part 
            A.ratio = A.ratio * ( norm.probs[result] )^(-1)
            
            # For transition prob ratio, need prob moving to current from proposed
            B.changed = unlist(B.holder[which.changed])
            norm.probs  = initial.probs[-c(observed.field, B.changed)] 
            norm.probs  = norm.probs/ sum(norm.probs) 
            desired     = which( c(1:d[which.changed])[-c(observed.field,B.changed)] == current.field)
            result.prob = norm.probs[desired] 
            
            A.ratio = A.ratio * result.prob
            #cat("The A ratio is", A.ratio ,"\n")
            
          } else if (Epart[r,"Aratio"]==2){
            
            # E(s) = 0, E*  = 1 
            
            A.ratio = 1/num.fields
            A.ratio = A.ratio*( norm.probs[result] )^(-1)
            A.ratio = A.ratio* prior.part 
            #cat("The A ratio is", A.ratio ,"\n")
            
          } else if (Epart[r,"Aratio"]==3){
            
            # E(s) = 1, E*  = 0 
            
            A.ratio = num.fields
            A.ratio = A.ratio*result.prob
            A.ratio = A.ratio*prior.part
            #cat("The A ratio is", A.ratio ,"\n")
          }    
          
          proposed.field      = B.holder[which.changed]
          
          propose.move        = matrix( c(f,r,from.block,to.block, A.ratio,proposed.field,which.changed),ncol= 7, nrow=1)
          
        
          proposal.info[[H]]  = propose.move 
        }
      }
    }
  }
  
  newlist <-list( "no.proposed" = H, "move.holder"= proposal.info )
  
  return( newlist ) 
  
} 

how_many_possible<-function( B.new, N, d, field, J, type1seeds, Legal_Index,blocks){
  
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
  
  possible.out = vector( "list", length = N )
  
  a = 1 
  
  for( i in (N+1):(2*N)){
    
    if( i %in% type1seeds){
      possible.out[a] = 0 
    } else{ 
      # which index are we in? 
      x = blocks[i]
      # which level of school? 
      s = SchoolLevel[i] 
      
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


propose.E <-function(E1,E2,labelsJ_E,nonseed.F1,nonseed.F2,gamma.samp,T2seeds1,T2seeds2){
  
  E1[,"Aratio"] = 0 
  E2[,"Aratio"] = 0 
  
  # File 1 
  #   for( i in nonseed.F1){
  #     lastjstar = E1[i,"jstar"] 
  #     #print(lastjstar)
  #     # For this record, choose a jstar
  #     if( lastjstar==0){
  #       possiblejstar = labelsJ_E
  #       
  #       n.options      = length(possiblejstar)
  #       jstar = sample( 1:n.options, 1)
  #       jstar = possiblejstar[jstar]
  #       
  #     } else{
  #       
  #       if( length(setJ_E) == 1 ){
  #         #There is only one choice for an error; we have to remove the restriction that we can't propose on the same field again. 
  #         jstar = labelsJ_E
  #       } else{ 
  #         possiblejstar = setdiff( labelsJ_E, lastjstar) 
  #         
  #         n.options      = length(possiblejstar)
  #         jstar = sample( 1:n.options, 1)
  #         jstar = possiblejstar[jstar]
  #       }
  #     }
  #     
  #     # Record the current value 
  #     E1current = E1[i,jstar]
  #     # Sample E[[1]][i,jstar]
  #     prob.jstar  = gamma.samp[1,jstar]
  #     E1[i,jstar] = sample( c(0,1), 1, replace=FALSE, prob = c(1-prob.jstar,prob.jstar) )  
  #     E1[i,"jstar"] = jstar  
  #     
  #     # Determing the Form of the acceptance ratio
  #     if( E1current == E1[i,jstar] & E1current==1){
  #       E1[i,"Aratio"] = 1 
  #     } else if( E1current ==0 & E1[i,jstar]==1){
  #       E1[i,"Aratio"] = 2 
  #     } else if( E1current ==1 & E1[i,jstar]==0){
  #       E1[i,"Aratio"] = 3 
  #     } 
  #   }
  
  # File 2 
  for( i in nonseed.F2){
    
    if( i %in% T2seeds1){
      jstar = labelsJ_E[2]
    } else if( i %in% T2seeds2){
      jstar = labelsJ_E[1]
    } else{ 
      lastjstar = E2[i,"jstar"] 
      # For this record, choose a jstar
      if( lastjstar==0){
        possiblejstar = labelsJ_E
      } else{
        
        if( length(labelsJ_E) == 1 ){
          #There is only one choice for an error; we have to remove the restriction that we can't propose on the same field again. 
          jstar = labelsJ_E
          possiblejstar = labelsJ_E
          
        } else{ 
          possiblejstar = setdiff( labelsJ_E, lastjstar) 
          
          n.options      = length(possiblejstar)
          jstar = sample( 1:n.options, 1)
          jstar = possiblejstar[jstar]
        }
        
      }
      n.options      = length(possiblejstar)
      jstar = sample( 1:n.options, 1)
      jstar = possiblejstar[jstar]
    }
    
    
    # Record the current value 
    E2current = E2[i,jstar]
    # Sample E[[2]][i,jstar]
    prob.jstar  = gamma.samp[2,jstar]
    E2[i,jstar] = sample( c(0,1), 1, replace=FALSE, prob = c(1-prob.jstar,prob.jstar) ) 
    E2[i,"jstar"] = jstar  
    
    # Determing the Form of the acceptance ratio
    if( E2current == E2[i,jstar] & E2current==1){
      E2[i,"Aratio"] = 1 
    } else if( E2current ==0 & E2[i,jstar]==1){
      E2[i,"Aratio"] = 2 
    } else if( E2current ==1 & E2[i,jstar]==0){
      E2[i,"Aratio"] = 3 
    }
  } 
  
  output = list( "E1" = E1, "E2" = E2 )
  return( output )  
}

category.block <-function(x,d,p){
  x = as.numeric(x)
  category=0
  for (j in 1:(p-1)) category=category+prod( d[(j+1):p])*(x[j]-1)
  category = category+x[p]
  return(category)
}

ginv <-function (X, tol = sqrt(.Machine$double.eps)) 
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X)) 
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                                               t(Xsvd$u[, Positive, drop = FALSE]))
}
#<bytecode: 0x0000000010f13b08>
#  <environment: namespace:MASS>

# How to remove all BUT the functions
#remove( list = ls()[sapply(ls(),function(n){!is.function(get(n))})])