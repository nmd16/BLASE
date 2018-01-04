# Implementing the Speed Up 

# We need to feed in A and Ainv to this, as well as the seeds. 

Bindex = c(1:nrow(File2011))[-c(type1seeds)]
C = diag(9)


# A different way 

# Find the difference between the original Y2 and the new ordering
u = Y2.in - Y2 
# Use this to update the inverse 


Gupdate_and_Impute_Speed<-function(Y1, Y2, Y2.in, Y2.info, New.Ordering, Sigma1.Part,Sigma2.Part,which.Y1,which.Y2,imputed.F2,imputed.F1,n.imputed.F2,n.imputed.F1,secondary.option, Proj2,V2,invert.part2,C,Bindex,Ainv){  
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ###########          Theta          #############
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###  
  
  # Create the model matrix 
  B       = Y2.info[Bindex,]
  Vtheta  = Y2.info
  B       = crossprod(B)
  
  invert.part = Ainv - Ainv%*%B%*%solve(C + Ainv%*%B )%*%Ainv
  
  Proj         = tcrossprod(invert.part,Vtheta)
  hat.beta     = Proj%*%Y1
  
  ResidSS      = (Y1 - Vtheta%*%hat.beta) 
  ResidSS      = crossprod(ResidSS)
  ResidSS      = .5*ResidSS
  
  sigmasq      = rigamma(1,Sigma1.Part,ResidSS) 
  sigma        = sqrt(sigmasq)
  var.beta     = sigmasq*invert.part
  beta         = mvrnorm( 1, hat.beta, var.beta ) 
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ##########  Initialization : Y2 Coefs  ##########
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###  
  
  if( secondary.option == "independent"){
    
    # The Y2s are not dependent on the blocking variables  
    
    barY2 = mean(Y2)
    Var.Part = sum( (Y2 - barY2)^2 )
    
    sigma2sq = rigamma(1,Sigma2.Part,.5*Var.Part)
    eta      = rnorm( 1, barY2, (1/n.type1)*sigma2sq)
    
  } else{
    
    # The Y2s are dependent on the blocking variables 
    # Now, we need to remember that Y2 is ordered to match up with Y1
    hat.eta      = Proj2%*%Y2.in
    est.mean.Y2  = V2%*%hat.eta
    
    ResidSS      = (Y2.in - est.mean.Y2) 
    ResidSS      = crossprod(ResidSS)
    ResidSS      = .5* ResidSS  
    
    sigma2sq     = rigamma(1,Sigma2.Part,ResidSS) 
    sqrt.sigma2  = sqrt(sigma2sq)
    var.eta      = sigma2sq*invert.part2
    eta          = mvrnorm( 1, hat.eta, var.eta ) 
    
    V2.imp       = V2[ imputed.F2 , ]%*%eta 
  } 
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  #######        Y2 Imputations           #########
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
  
  if( secondary.option == "independent"){
    mean.Y2impute = sigma2sq^(-1)*eta - mean.Y2impute
  } else{
    mean.Y2impute = sigma2sq^(-1)*V2.imp - mean.Y2impute
  }
  
  var.Y2impute  = ( sigma2sq^(-1) + sigmasq^(-1)*(theta.small^2) )^(-1)
  mean.Y2impute = var.Y2impute*mean.Y2impute 
  
  imputed.all = rnorm( n.imputed.F2, mean.Y2impute, sqrt(var.Y2impute) )
  
  # Go to the original Y2 data (input ordering) and update the imputations 
  Y2.in[ imputed.F2 ]            = imputed.all 
  
  # Go to where the imputations are in the model matrix and impute
  Vtheta[F1.matched.to.imputations,which.Y2] = imputed.all 
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  #######        Y1 Imputations           #########
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  Vtheta         = Vtheta[imputed.F1,]
  
  imputed.all    = rnorm( n.imputed.F1, Vtheta %*% as.matrix(beta), sqrt(sigmasq) ) 
  
  Y1[imputed.F1] = imputed.all  
  
  #################################################
  ## Step 6: Store and Release the Output 
  #################################################
  
  newlist <-list( "Y1" = Y1, "Y2" = Y2.in, "beta" = beta, "sigmasq" = sigmasq, "sigma" = sigma, "eta" = eta, "sigma2sq"=sigma2sq, "theta" =theta.small, "impVar" = sqrt(var.Y2impute))
  return( newlist ) 
  
}

system.time(
gutman_mcmc(is.test = "no",MCMCSeed, its=20,burnin=0, thinning=1, gap=0, n.RIB, RIB, File1_data, File2_data,G.K,which.Y2=2, which.Y1=1,type1Seeds=type1seeds, where.coefs.primary, where.coefs.secondary, p, p2, n.type1,needs.C,reps,PermSampling,secondary.option="dependent", n1, n2,col.Y1, col.Y2,namepart="Nope")
)

# The inital beta draw is 37.99038 0.8878842 0.2253935 0.6440161 -0.7753142 -0.8061862 -0.5348364 -0.6807192 -1.692191 -0.9508519 -0.3925808 -0.1452956 0.2836036 -0.02160754 1.272814 1.066031 0.8044887 1.408828 -0.2208385 
#The inital eta draw is 366.4405 -9.632747 -0.5014186 3.566508 -5.644215 -3.369255 -4.242791 -2.273526 1.097845 -1.499896 2.334574 3.710713 3.525368 1.637654 
#user  system elapsed 
#41.49    0.09   43.68

system.time(
  gutman_mcmc_Speed(is.test = "no",MCMCSeed, its=20,burnin=0, thinning=1, gap=0, n.RIB, RIB, File1_data, File2_data,G.K,which.Y2=2, which.Y1=1, type1seeds, where.coefs.primary, where.coefs.secondary, p, p2, n.type1,needs.C,reps,PermSampling,secondary.option="dependent", n1, n2,col.Y1, col.Y2,namepart="Nope")
)

#The inital beta draw is 37.99038 0.8878842 0.2253935 0.6440161 -0.7753142 -0.8061862 -0.5348364 -0.6807192 -1.692191 -0.9508519 -0.3925808 -0.1452956 0.2836036 -0.02160754 1.272814 1.066031 0.8044887 1.408828 -0.2208385 
#The inital eta draw is 366.3891 -9.620824 -0.483852 3.55737 -5.566724 -3.23763 -4.251198 -2.298108 1.107696 -1.500013 2.341641 3.719649 3.525038 1.637209 
#user  system elapsed 
#41.49    0.07   43.03

# What do we need in? 
# The V2 matrix (which is also Vtheta)
# The Y2
# The Y1
# The Initial Ordering
# The Initial BlockRow
# Which records were imputed 

# We ignore the first part, because it happens only once. 

gutman_mcmc_Speed <-function(is.test,seed, its, burnin, thinning, gap, n.RIB, RIB, File1_data, File2_data,K,which.Y2, which.Y1,type1seeds, where.coefs.primary, where.coefs.secondary, p,p2, n.type1,needs.C,reps, PermSampling, secondary.option, n1, n2,col.Y1, col.Y2,namepart){
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ############   Constants Chunk  #################
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  # Initialize the Random Seed 
  set.seed(seed)
  
  s = 1 
  
  File2Records             = as.numeric( rownames(File2_data) )
  
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
  # Which records from File2 were imputed?
  if( nrow.F1 > n1){
    imputed.F1   = (n1+1):nrow.F1
    n.imputed.F1 = length(imputed.F1)
  } else{     
    imputed.F1   = NULL
    n.imputed.F1 = 0   
  }
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ############ Data Storage Chunk #################
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  # We need to store...
  counter               = 1
  # The primary regression coefficients
  out.coefs            <- NULL
  # The File 2 Lambdas
  out.lambda           <-NULL
  # The secondary regression coefficients
  out.coefs2          <- NULL
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ###########  Initialization  Theta #############
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  # We begin by using ONLY the seeds information 
  
  Y1 = File1_data[type1seeds,col.Y1]
  Y2 = File2_data[type1seeds,col.Y2]
  
  # Create the model matrix based on the seeds
  Vtheta = cbind( 1, File2_data[type1seeds, where.coefs.primary ] )
  Vtheta = as.matrix(Vtheta)
  A = crossprod(Vtheta)
  Ainv = solve(A)
  
  Proj         = tcrossprod(Ainv,Vtheta)
  hat.beta     = Proj%*%Y1
  
  ResidSS      = (Y1 - Vtheta%*%hat.beta) 
  ResidSS      = crossprod(ResidSS)
  ResidSS      = .5* ResidSS
  
  Sigma1.Part = (n.type1 - p )/2
  
  sigmasq      = rigamma(1,Sigma1.Part,ResidSS) 
  sqrt.sigma   = sqrt(sigmasq)
  var.beta     = sigmasq*Ainv
  beta         = mvrnorm( 1, hat.beta, var.beta ) 
  cat( "The inital beta draw is",beta,"\n")
  
  rm( Sigma1.Part,var.beta,ResidSS,hat.beta, Proj)
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ##########  Initialization : Y2 Coefs  ##########
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  if( secondary.option == "independent"){
    
    # The Y2s are not dependent on the blocking variables 
    
    Sigma2.Part = (n.type1 - p2 )/2
    
    # Again, we are just using the seeds right now. 
    
    barY2 = mean(Y2)
    Var.Part = sum( (Y2 - barY2)^2 )
    
    sigma2sq = rigamma(1,Sigma2.Part,.5*Var.Part)
    eta      = rnorm( 1, barY2, (1/n.type1)*sigma2sq)
    
    rm(Var.Part,barY2,Sigma2.Part)
    
  } else{
    
    # The Y2s are dependent on the blocking variables 
    
    # Create the model matrix based on the seeds
    V2 = cbind( 1, File2_data[type1seeds, where.coefs.secondary ] )
    V2 = as.matrix(V2)
    
    V2TV2        = crossprod(V2)
    Gen          = zapsmall( eigen( V2TV2 )$values )
    Gen          = ifelse( sum(Gen==0) > 0 , "TRUE","FALSE")
    invert.part  = invert_XTX( V2TV2 , Gen )
    
   #t.V2         = t(V2)
    Proj         = tcrossprod(invert.part,V2)
    hat.eta      = Proj%*%Y2
    est.mean.Y2  = V2%*%hat.eta
    
    ResidSS      = (Y2 - est.mean.Y2) 
    ResidSS      = crossprod(ResidSS)
    ResidSS      = .5* ResidSS
    
    Sigma2.Part = (n.type1 - p2 )/2
    
    sigma2sq     = rigamma(1,Sigma2.Part,ResidSS) 
    sqrt.sigma2  = sqrt(sigma2sq)
    var.eta      = sigmasq*invert.part
    eta          = mvrnorm( 1, hat.eta, var.eta ) 
    cat( "The inital eta draw is",eta,"\n")  
    
    
    rm( Sigma2.Part,var.eta,ResidSS,hat.eta, Proj,Gen,invert.part,V2TV2,est.mean.Y2) 
    
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
    
    beta.impute = beta[ -which.Y2 ] 
    theta.small = beta[ which.Y2 ] 
    
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
    imputed.all   = rnorm( n.imputed.F2, mean.Y2impute, sqrt(var.Y2impute) )
    
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
    
    imputed.all = rnorm( n.imputed.F1, Vtheta %*% as.matrix(beta), sqrt(sigmasq) ) 
    
    File1_data[ imputed.F1, col.Y1 ] = imputed.all
    
    
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
    ##########  Initialization : Deletions  ##########
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
    
    # Delete everything that we no longer need.   
  
    rm( Vtheta,imputed.all,mean.Y2impute, var.Y2impute,beta.impute,theta.small,V1,Y1.matched.to.imputations,F1.matched.to.imputations,Y1,Y2)
  }
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ############    Run the MCMC    #################
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%### 
  
  # Added for the speed up in matrix inversion 
  Bindex = c(1:nrow(File2_data))[-c(type1seeds)]
  C      = diag(p)
  
  Sigma1.Part = (nrow.F1 - p )/2
  Sigma2.Part = (nrow.F2 - p2 )/2
  
  V2 = cbind( 1, File2_data[, where.coefs.secondary ] )
  V2 = as.matrix(V2)
  
  V2TV2        = crossprod(V2)
  Gen          = zapsmall( eigen( V2TV2 )$values )
  Gen          = ifelse( sum(Gen==0) > 0 , "TRUE","FALSE")
  invert.part2 = invert_XTX( V2TV2 , Gen )
  
  #t.V2         = t(V2)
  Proj2        = tcrossprod(invert.part2,V2)#%*%t.V2
  
  rm(V2TV2,Gen)#,t.V2)
  
  BlockRow = File2_data[,"BlockRow"]
  
  # Take the current ordering of the File 2 data 
  F1.Lambda  = 1:nrow.F1
  F2.Lambda  = File2_data[,"lambda"]
  F2.Lambda[type1seeds] = type1seeds
  Y1         = File1_data[,col.Y1]
  Y2         = File2_data[,col.Y2]
  # Store a copy in the original ordering
  Y2.in      = Y2
  D.Holder   = File2_data[,where.coefs.primary]
  # Ordered
  V.Holder   = cbind(1, D.Holder[order(F2.Lambda),] )
  
  New.Ordering = F2.Lambda
  
  for( s in 1:its){     
    
    # Sample a Permutation C for every block which requires it
    for( k in needs.C){
      # Which records from File 1 are in the block?
      F1.records = RIB[[1]][[k]] 
      # Which records from File 2 are in the block?
      F2.records = RIB[[2]][[k]]
      # What is the current ordering of the records within the block?
      order.in   = BlockRow[F2.records]
      # How many pairs are in the block?
      pairs      = n.RIB[1,k]
      # Store the File 1 block information
      F1.holder  = Y1[F1.records]
      # Store the File 2 block information
      F2.holder  = D.Holder[F2.records,]
      # Sample a new permutation 
      block.order.out = SamplePerm(F1.holder,F2.holder,pairs,order.in,beta,sqrt.sigma,S.ind=PermSampling[k],reps,k)     
      #Update the order of the File 2 Data
      New.Ordering[F2.records] = F1.Lambda[F1.records][block.order.out]
      #Update the Block Row 
      BlockRow[F2.records]    = block.order.out 
    } 
    
    # Reorder based on the sampled C 
    Y2                  = Y2.in[order(New.Ordering)]
    V.Holder[,which.Y2] = Y2    
    
    holder = Gupdate_and_Impute_Speed(Y1, Y2, Y2.in,V.Holder, New.Ordering, Sigma1.Part,Sigma2.Part,which.Y1,which.Y2,imputed.F2,imputed.F1,n.imputed.F2,n.imputed.F1,secondary.option, Proj2,V2,invert.part2,C,Bindex,Ainv)
    
    #Update based on the current imputations 
    Y1         = holder$Y1
    Y2.in      = holder$Y2 # This is set back to the original ordering, but with updated imputations
    beta       = holder$beta
    sigmasq    = holder$sigmasq
    sigma      = sqrt(sigmasq)
    eta        = holder$eta
    sigma2sq   = holder$sigma2sq
    rm(holder)  
    
    D.Holder[,col.Y2] = Y2.in
    
    # Store the output 
    
    if( s > burnin & s%%thinning==0){
      
      # The regression coefficients
      out.coefs <- rbind( out.coefs, c(beta,sigmasq) )
      # The File 2 Lambdas
      out.lambda <-rbind( out.lambda,New.Ordering )
      # The mean and variance for the File 2 continuous variable 
      out.coefs2  <- rbind( out.coefs2, c(eta,sigma2sq) )
      
      G.outlist <- list( "out.coefs" = out.coefs, "out.lambda" = out.lambda, "out.coefs2" = out.coefs2)
      
      # Save the current output 
      save( G.outlist, file = paste( "GOut", namepart, ".RData", sep= ""))     
      
    }
    
    if(s%%printgap == 0){ print(s) }
    
  }
  
}