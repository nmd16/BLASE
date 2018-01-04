# Needed Functions 


run_rlMCMC<-function( MCMCSeed, File1, File2, method, its, burnin, DPits, DPburnin,gap, reps, thinning, DPthinning, printgap, Comment, coefs.primary,coefs.primary.all,coefs.secondary,col.Y1, col.Y2, secondary.option, secondary.option2, ObsBlock_restriction, BlockSeed_restriction, error.file, error.fields, blockingVars,matches, a.beta.prior, b.beta.prior, Hstar, DPOffset,Y1.name,Y2.name,namepart,threshold,complete.out,type2seeds){
  
  library(rlMCMC)
  
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
  n.RIB[2,] = sapply(RIB[[1]],length)
  
  
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
  ## Format File1 and File 2 for the regression ##
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
  
  File1_data = File1
  File2_data = File2 
  
  # Now, we need to make sure that ALL of the coefficients which are in either model end up in File 1 and File2 ( changing this changes a lot of coding, and we can do it later, 
  #but not now)
  
  coefs.all = union(coefs.primary,coefs.secondary)
  f1 <- as.formula(paste("File1[,",col.Y1,"]~", paste(coefs.all, collapse="+")))
  modelPart <- model.matrix( lm(f1,data=File2))
  
  File1_data = cbind( File1[,col.Y1], modelPart[,-1])
  colnames(File1_data)[col.Y1] = Y1.name
  File2_data = cbind( File2[,col.Y2], modelPart[,-1])
  colnames(File2_data)[col.Y2] = Y2.name
  
  # Create the Primary Model Matrix 
  f1 <- as.formula(paste("File1[,",col.Y1,"]~", paste(coefs.primary.all, collapse="+")))
  modelPartPrimary <- model.matrix( lm(f1,data=File2))
  
  # Create the Secondary Model Matrix 
  f1 <- as.formula(paste("File2[,",col.Y2,"]~", paste(coefs.secondary, collapse="+")))
  modelPartSecondary <- model.matrix( lm(f1,data=File2))
  rm(f1)
  
  # In the future, I think I want to work with vectors instead of matrices, and just replace the vectors that we need with the reorderings. This might make the code more efficient, but for now, leave it alone. 
  
  # Now, I have to know which of these columns are used for which regressions, and this somehow needs to be automated and stored 
  
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
    gutman_mcmc("no",MCMCSeed, its,burnin, thinning, gap, n.RIB, RIB, File1_data, File2_data,K,which.Y2, which.Y1,type1Seeds = type1seeds, where.coefs.primary, where.coefs.secondary, p, p2, n.type1,needs.C,reps,PermSampling,secondary.option="dependent", n1, n2,col.Y1, col.Y2,namepart)
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
    gutman_mcmc("no",MCMCSeed, its,burnin, thinning, gap, n.RIB, RIB, File1_True, File2_True,K,which.Y2, which.Y1,type1seeds, where.coefs.primary, where.coefs.secondary, p, p2, n.type1,needs.C,reps,PermSampling,secondary.option="dependent",n1,n2,col.Y1,col.Y2,paste(namepart,"Test",sep=""))
    
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
    
    #model <-CreateModel( Bhat, NULL, Hstar, 0, 0.25,0.25)
    #model$Run(1000,10000,2)
    #dpmpm <-model$snapshot
    
    
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
    
    labelsJ_E = which(blockingVars==error.fields)
    
  
    new_mcmc(its, burnin ,thinning,gap,MCMCSeed,namepart,File1_data, File2_data,is.test="no",printgap, Legal_Index, Y2.name ,n1, n2,p,p2,secondary.option = "dependent", J,nonseed.F1,nonseed.F2,n.error1,n.error2,Bhat1, Bhat2,type1seeds,PermSampling,SeedsOnly, RIB,n.RIB, Hstar, where.coefs.primary, where.coefs.secondary, where.coefs=NULL, Comment,Imp.Holder1,Imp.Holder2, J_E, K, realmin, reps, col.Y1, col.Y2, count.matches,d,labelsJ_E ,n.type1,threshold,which.Y1,which.Y2,Index_to_Block, a.beta.prior, b.beta.prior, secondary.option2,error.file,DPoffset, dpmpm,model,Bhat,where.bvars.F1, where.bvars.F2,DPits, DPburnin, DPthinning)
  }
  
} 



make_traceplots <-function(object,truth,name,comp.object,option){ 
  
  max.point = max( max(object), max(comp.object), truth )
  min.point = min( min(object), min(comp.object), truth )
  
  
  plot( object ,type="l",ylab="Parameter", ylim = c(min.point, max.point), main=paste("Trace Plot for ", name, sep = ""))
  
  #abline(h=1,col="blue",lwd=3) #Prior 
  abline(h=truth,col="green",lwd=2) # Target 
  
  mcmc.mean = mean(object)
  mcmc.sd   = sd(object)
  
  type.mean = ifelse( abs( truth- mcmc.mean)> .05, 1, 2)
  
  if(option=="compare"){
    comp.object.mean = mean(comp.object)
    comp.object.sd   = sd(comp.object)
    #lines(comp.object, type="l")
    abline(h = comp.object.mean + 1.96*comp.object.sd, col="magenta",lty=5  )
    abline(h = comp.object.mean - 1.96*comp.object.sd, col="magenta",lty=5  )
    abline(h = comp.object.mean, col="magenta", lwd = 2,lty= type.mean)
  }
  
  if(option=="COUNTER"){
    comp.object.mean = sum( apply( comp.object, 1, function(x) ifelse(sum(x) == 0, 0, 1)))
    cat("There are ",comp.object.mean," records in the wrong block in the starting data.")
  }    
  
  abline(h = mcmc.mean,col="blue",lwd=3, lty= type.mean) # Posterior Mean 
  abline(h = mcmc.mean + 1.96*mcmc.sd, col="blue",lty=5 ,lwd = 2  )
  abline(h = mcmc.mean - 1.96*mcmc.sd, col="blue",lty=5 , lwd = 2  ) 
  
  
}

make_traceplotsALL <-function(object,truth,name,comp.object,option,test.object){ 
  
  max.point = max( max(object), max(comp.object), max(test.object), truth )
  min.point = min( min(object), min(comp.object), min(test.object), truth )
  
  
  plot( object ,type="l",ylab="Parameter", ylim = c(min.point, max.point), main=paste("Trace Plot for ", name, sep = ""))
  
  #abline(h=1,col="blue",lwd=3) #Prior 
  abline(h=truth,col="green",lwd=2) # Target 
  
  mcmc.mean = mean(object)
  mcmc.sd   = sd(object)
  
  type.mean = ifelse( abs( truth- mcmc.mean)> .05 | abs( truth- mean(test.object))> .05 , 1, 2)
  
  if(option=="compare"){
    comp.object.mean = mean(comp.object)
    comp.object.sd   = sd(comp.object)
    #lines(comp.object, type="l")
    abline(h = comp.object.mean + 1.96*comp.object.sd, col="darkmagenta",lty=5,lwd=2  )
    abline(h = comp.object.mean - 1.96*comp.object.sd, col="darkmagenta",lty=5,lwd=2 )
    abline(h = comp.object.mean, col="darkmagenta", lwd = 2,lty= type.mean)
    
    test.object.mean = mean(test.object)
    test.object.sd   = sd(test.object)
    #lines(comp.object, type="l")
    abline(h = test.object.mean + 1.96*comp.object.sd, col="darkgreen",lty=5,lwd=2  )
    abline(h = test.object.mean - 1.96*comp.object.sd, col="darkgreen",lty=5 ,lwd=2 )
    abline(h = test.object.mean, col="darkgreen", lwd = 2,lty= type.mean)
  }
  
  if(option=="COUNTER"){
    comp.object.mean = sum( apply( comp.object, 1, function(x) ifelse(sum(x) == 0, 0, 1)))
    cat("There are ",comp.object.mean," records in the wrong block in the starting data.")
  }    
  
  abline(h = mcmc.mean,col="blue",lwd=3, lty= type.mean) # Posterior Mean 
  abline(h = mcmc.mean + 1.96*mcmc.sd, col="blue",lty=5 ,lwd = 2  )
  abline(h = mcmc.mean - 1.96*mcmc.sd, col="blue",lty=5 , lwd = 2  ) 
  
  
}

category.block<-function(x,d,J){
  x = as.numeric(x)
  category=0
  for (j in 1:(J-1)) category=category+prod( d[(j+1):J])*(x[j]-1)
  category = category+x[J]
  return(category)
}

permute.places<-function(vec1){
  n=length(vec1)
  vec1.ordering=sample(1: n , n ,replace=FALSE)
  vec2=vec1[order(vec1.ordering)];vec2
}

records.in.block <-function(File1_data, File2_data, num.blocks){
  
  # Storage: How many records in each block?  
  n.RIB = matrix( NA, nrow = 2, ncol = num.blocks )
  #Storage: Which records in each block? 
  RIB     = vector( "list", length = 2 )
  RIB[[1]]     = vector( "list", length = num.blocks )
  RIB[[2]]     = vector( "list", length = num.blocks )
  
  # Start with File 1 
  data = File1_data[,"Block"]
  for(f in 1:2){
    if( f ==2){
      data = File2_data[,"Block"]
    }
    for(k in 1:num.blocks){
      holder        = which(data==k)
      RIB[[f]][[k]] = holder
      n.RIB[f,k]    = length(holder)
    }
  }
  new.list <- list( "RIB" = RIB, "n.RIB" = n.RIB )
  return(new.list)
}

make_imputed_holder<-function( data, ToImpute ){
  # Step 1: Store only one row from each block 
  ImpHolder = data[!duplicated(data[,"Block"]),]
  # Step 2: Set all record-specific row entries to 0 
  ImpHolder[,c(ToImpute)]= ImpHolder[,"BlockRow"] = ImpHolder[,"Seed"] = ImpHolder[,"lambda"] = 0
  # Step 3: Set the binary imputed indicator to 1 
  ImpHolder[,"Imputed"] = 1
  # Step 4: Change the row names
  rownames(ImpHolder) = 1:nrow(ImpHolder)
  #Step 4: Transpose ImpHolder so that each record is a column
  ImpHolder = as.data.frame( t(ImpHolder) )
  #Step 5: Convert to a list 
  ImpHolder = as.list(ImpHolder)
  #Return the matrix 
  return(ImpHolder)
}

add_imputed_rows <-function( in_Bl,ImpHolder1, ImpHolder2,RIB,File1_data,File2_data){
  
  to.impute = in_Bl[1,]-in_Bl[2,]
  
  F1 = which( to.impute < 0 )
  F2 = which( to.impute > 0 )
  
  #Add in rows for each of the blocks that need it in File 1
  for(k in F1 ){  
    timeshere = abs(to.impute[k])
    holder = matrix( rep(ImpHolder1[[k]],timeshere),nrow = timeshere, byrow=T)
    colnames(holder)=colnames(File1_data)
    preLength = nrow(File1_data)
    File1_data = rbind( File1_data, holder)
    postLength = nrow(File1_data)
    RIB[[1]][[k]] = c( RIB[[1]][[k]] , (preLength+1):postLength )
  }
  
  #Add in rows for each of the blocks that need it in File 2
  for(k in F2 ){  
    timeshere = abs(to.impute[k])
    holder = matrix( rep(ImpHolder2[[k]],timeshere),nrow = timeshere, byrow=T)
    colnames(holder)=colnames(File2_data)
    preLength = nrow(File2_data)
    File2_data = rbind( File2_data, holder)
    postLength = nrow(File2_data)
    RIB[[2]][[k]] = c( RIB[[2]][[k]] , (preLength+1):postLength )
  }
  
  new.list <-list( "File1_data"= File1_data, "File2_data" = File2_data,"RIB"=RIB)
  
  return(new.list)
}