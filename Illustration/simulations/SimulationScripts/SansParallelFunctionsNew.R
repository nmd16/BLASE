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
  
  BlockRow.F2 = File2_data[,"BlockRow"]
  
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
      order.in   = BlockRow.F2[F2.records]
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
      BlockRow.F2[F2.records]    = block.order.out 
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
      #out.coefs <- rbind( out.coefs, c(beta,sigmasq) )
      # The File 2 Lambdas
      #out.lambda <-rbind( out.lambda,New.Ordering )
      # The mean and variance for the File 2 continuous variable 
      #out.coefs2  <- rbind( out.coefs2, c(eta,sigma2sq) )
      
      # We want to save these to text files. 
      
      write( c(beta,sigmasq), file = paste("G",namepart,"ThetaOut.txt", sep=""), append = TRUE)
      write( c(eta,sigma2sq), file = paste("G",namepart,"EtaOut.txt", sep=""), append = TRUE)
      Lambda.out <- sum(New.Ordering[1:n2]==1:n2)
      write( Lambda.out, file = paste("G",namepart,"LambdaOut.txt", sep=""), append = TRUE)
      
      #G.outlist <- list( "out.coefs" = out.coefs, "out.lambda" = out.lambda, "out.coefs2" = out.coefs2)
      
      # Save the current output 
      #save( G.outlist, file = paste( "GOut", namepart, ".RData", sep= ""))     
      
    }
    
    if(s%%printgap == 0){ print(s) }
    
  }
  
}

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
    
    # Now, because of the restriction, we could have the possibility of no valid moves. (1) If a record in the incorrect block, we may have no choice but to move it back to the observed (if we have only two levels). This means that we need to just skip that move, and proceed to the next. I am going to gamble that that is the issue. 
    
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


accept.moves.func <- function(Accept.Count,a.beta.prior, B,b.beta.prior,beta.impute,beta.star,Bhat,Bhat1,Bhat2,Block,Block_Ops,blockingVars,BlockRow,BlockRow.F1,BlockRow.F2,blocks,blocks.out,blocks1,blocks2,BlockSeed_restriction,burnin,coefs.all,coefs.primary,coefs.primary.all,coefs.secondary,col.Y1,col.Y2,complete.out,count.matches,current.B,d,D.Holder, DPburnin,DPits,dpmpm,DPOffset,DPthinning,E1,E2,Empty.Blocks,error.fields,error.file,error.perc,Estar,eta,F1.Lambda,F1.records,F2.Lambda,File1_data,File2_data,gamma.samp,gap,Hstar,HstarM1,hvec,i,Imp.Holder1,Imp.Holder2,Imputed.Block,imputed.F1,imputed.F2,Imputed.Records,in.error,Index_to_Block,InF1,InF2,its,J,K,labelsJ_E,lambda.out,Legal_Index, matches,max.RIB,MCMCSeed, MoveCand_F2,N,n.error1,n.error2,n.imputed.F1,n.imputed.F2,n.RIB,n.type1,n1,n2,needs_C, New.Ordering,non.seed.RIB,nonseed.F2,nrow.F1,nrow.F2,num.matched,ObsBlock_restriction,observed.blocks,p,p2,pairs,perc.in,perc.needed,perc.type1,realmin,reps,result,result.out,RIB,secondary.option,secondary.option2,SeedsOnly,sigma,Sigma1.Part,sigma2,Sigma2.Part,sigma2sq,sigmasq,sqrt.var.Y2impute,theta.small,threshold,to.sampleC,truthSim,type1seeds,type2seeds,V.Holder,V2,V2.imp,var.Y2impute,where.bvars.F1,where.bvars.F2,where.coefs,where.coefs.primary,where.coefs.secondary,which.Y1,which.Y2,Y1,Y1.name,Y2,Y2.in,Y2.name,Y2.out){
  
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
            where.replacing       = ImpF1Which[[holder.to]][where.replacing]
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
            where.replacing = ImpF2Which[[holder.from]][where.replacing]
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
            where.replacing       = ImpF2Which[[ holder.to ]][ where.replacing ]
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
            where.replacing = ImpF1Which[[holder.from]][where.replacing]
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
        Accepted.Moves.Holder = NULL      
      } 
      
      ########################################################################################################
      # If we accept :
      ########################################################################################################
            
      if( accept == 1 ){ 
        AddedTo = ifelse( class(Added$To) != "NULL", Added$To, 0 )  
        AddedFrom = ifelse( class(Added$From) != "NULL", Added$From, 0 ) 
        Accepted.Moves.Holder = c(holder.all,"Y1Imp" = AddedTo, "Y2Imp" = AddedFrom)  
      } 
    
      return(Accepted.Moves.Holder)
    
    }

# This is the edited on that works when are are not parallelizing 

run_Extension <-function(){
  
    record.holder  =  Reupdate_RIB( File1_data, File2_data,RIB ,n.RIB , 1:K)
    
    #record.holder  = records.in.block( File1_data, File2_data,K )
    n.RIB          = record.holder$n.RIB
    RIB            = record.holder$RIB
    rm(record.holder)
    
    max.RIB        = apply( n.RIB, 2, max)
    max.RIB        = max.RIB - count.matches
    
    PermSampling   = ifelse( max.RIB > 1, ifelse(max.RIB < threshold, 1, 2), 0 )
    SeedsOnly      = which(max.RIB==0)
    # This vector takes a 1 if we use exact sampling, 2 if we use MH and 0 otherwise.
    needs.C        = which(PermSampling > 0)
    
    holder     = get_starting_permutation( K, RIB, matches,File1_data, File2_data,nrow.F1,nrow.F2)
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
    
    C      = diag(p)
    
    
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
  Imputed.Block <- NULL 
  # Which records in a block from a given file are imputations?
  ImpF2Which      = mapply( function(x) which( x > n2), RIB[[2]])
  ImpF1Which      = mapply( function(x) which( x > n1), RIB[[1]])
  # Identify the number of records in each block which are imputations 
  Imputed.Block = cbind(Imputed.Block,unlist(mapply( length, ImpF1Which)))
  Imputed.Block = cbind(Imputed.Block,unlist(mapply( length, ImpF2Which)))
  
  ##  %%%%%%%%%%%%%%%%%%%  ##
  #### Initialize Gamma  ####
  ### %%%%%%%%%%%%%%%%%%  ###  
  
  gamma.samp = matrix( 0, nrow = 2, ncol = J)
    
   # a.beta.prior[,6]= c(0,0)
    #b.beta.prior[,6]=c(0,0)
    #a.beta.prior[,5]=c(2,2)
    #b.beta.prior[,5] =c(10,10)
  
  for( u in labelsJ_E){
    #gamma.samp[1,u] = rbeta( 1, a.beta.prior[1,u], b.beta.prior[1,u] )
    gamma.samp[error.file,u] = rbeta( 1, a.beta.prior[2,u], b.beta.prior[2,u] )
  }
    
  print(gamma.samp)
  
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
    
    Bhat1 = B[1:n1,]
    Bhat2 = B[(n1+1):N,]
  
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
  
  if( class(type1seeds) != "NULL"){
    
    Type1SeedIndicator = 1
    
    Y1 = File1_data[type1seeds,col.Y1]
    Y2 = File2_data[type1seeds,col.Y2]
    
    # Create the model matrix based on the seeds
    Vtheta = cbind( 1, File2_data[type1seeds, where.coefs.primary ] )
    
    # Create the model matrix based on the seeds
    Vtheta = as.matrix(Vtheta)
    A = crossprod(Vtheta)
    Ainv = solve(A)
    
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
    
    if( class(type1seeds) != "NULL"){
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
    if( class(type1seeds) != "NULL"){
      V2 = cbind( 1, File2_data[type1seeds, where.coefs.secondary ] )
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
    
    beta.impute   = beta.star[-which.Y2 ] 
    theta.small   = beta.star[ which.Y2 ] 
    
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
  
  for( j in 1:J){
    if(class(Bhat[,j]) != "factor"){
      Bhat[,j] = as.factor(Bhat[,j])
    }
  }
  
  model <-CreateModel( Bhat, NULL, Hstar, 0, 0.25,0.25)
  model$Run(DPburnin,DPits,DPthinning)
  dpmpm <-model$snapshot
  
  nrow.F1              = nrow(File1_data) 
  nrow.F2              = nrow(File2_data) 
  rownames(File1_data) = 1:nrow.F1
  rownames(File2_data) = 1:nrow.F2
 
     #################################################################################
    ####  Sample the Permutation Conditional on the Current Blocking Structure ######
    #################################################################################
    
    to.sampleC = which(PermSampling>0)
    to.sampleC = setdiff( to.sampleC, Empty.Blocks) # This SHOULD be superfluous, but is an okay check    
    
    BlockRow.F1 = File1_data[,"BlockRow"]
    BlockRow.F2    = File2_data[,"BlockRow"]
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
      order.in   = BlockRow.F2[F2.records]
      # How many pairs are in the block?
      pairs      = n.RIB[1,k]
      # How many type 1 seeds?
      if( class(type1seeds) != "NULL"){
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
      BlockRow.F2[F2.records]    = block.order.out 
    } 
    
    # Reorder based on the sampled C 
    Y2                = Y2.in[order(New.Ordering)]
    V.Holder[,col.Y2+1] = Y2
    F2.Lambda = New.Ordering 
    
    # Store in the File 2 data space 
    File2_data[,"lambda"]   = New.Ordering
    File2_data[,"BlockRow"] = BlockRow.F2
  
  sys <-proc.time() 
    
  RIB1Base  = mapply( function(x) x[which(x %in% type1seeds)], RIB[[1]])
  RIB2Base  = mapply( function(x) x[which(x %in% type1seeds)], RIB[[2]])
    
  nonseed.F2 = nonseed.F2[ which(nonseed.F2 <= n2)]
    
  #save.image("PreBLASE.RData")
  #nonseed.F1 = setdiff(nonseed.F1,type1seeds)
  for( s in 5:100){ 
    
   # File1Pre = File1_data
   # File2Pre = File2_data
    
    Accept.Out = 0 
    
    if(s%%printgap == 0){ print(s) }
    
    ######################
    ### Propose a new E ##
    ######################
    
    Estar = propose.E(E1,E2,labelsJ_E,nonseed.F1,nonseed.F2,gamma.samp,T2seeds1,T2seeds2)  
    
    ## %%%%%%%%%%%%% ##
    ## DP Parameters ## 
    ## %%%%%%%%%%%%% ##
    
    B = do.call(rbind, current.B )
    
    for( j in 1:J){
      if(class(B[,j])!="factor"){
        B[,j] = as.factor(B[,j])
      }
    }
    
    UpdateX(model,B)
    model$Run(0,2,2)
    dpmpm <-model$snapshot    
    
    ###########################
    ##  Propose Block Moves  ##
    ###########################
    
    Move.F1 = which( Estar[[1]][,"Aratio"]>0)
    Move.F2 = which( Estar[[2]][,"Aratio"]>0)  
    
    if( error.file==2 & length(Move.F2)> 0){
    
    move.holder =  propose.B( J, d, File1_data[,"Block"], File2_data[,"Block"], Estar, Move.F1, Move.F2, current.B, Bhat1, Bhat2, Blocking_WithError, Index_to_Block, Legal_Index,dpmpm$psi,dpmpm$z+1, File2_data[,"lambda"], File1_data[,"lambda"],gamma.samp,A.probs,BlockSeed_restriction="no",SeedsOnly,possibleEthnic,possibleSchool )
    
    
    if((move.holder$no.proposed)>0){
    ###############################
    ## Accept/Reject Block Moves ##
    ###############################    
  
    beta.impute    = beta.star[-2]
    theta.small    = beta.star[2]
  
    Move.Holder    = move.holder$move.holder    
  
    Accept.Count   = 0     
    
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%###
    ####      Run the Loop       ####
 
    #result.out = AcceptFunction(Move.Holder)    
    
    uSet = 1:length(Move.Holder)
  
  result.out = NULL 
  STARTSPACE = noquote( ls() )
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
          orderC.2T.orig    = BlockRow.F2[ F2.records ]
          
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
            where.replacing       = ImpF1Which[[holder.to]][where.replacing]
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
        orderC.2F.orig  = BlockRow.F2[ F2.records ]
        
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
            where.replacing = ImpF2Which[[holder.from]][where.replacing]
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
          orderC.2T.orig    = BlockRow.F2[ F2.records ]
          
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
            where.replacing       = ImpF2Which[[ holder.to ]][ where.replacing ]
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
        orderC.2F.orig  = BlockRow.F2[ F2.records ]
        
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
            where.replacing = ImpF1Which[[holder.from]][where.replacing]
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
        Accepted.Moves.Holder = NULL      
      } 
      
      ########################################################################################################
      # If we accept :
      ########################################################################################################
            
      if( accept == 1 ){ 
        AddedTo = ifelse( class(Added$To) != "NULL", Added$To, 0 )  
        AddedFrom = ifelse( class(Added$From) != "NULL", Added$From, 0 ) 
        Accepted.Moves.Holder = c(holder.all,"Y1Imp" = AddedTo, "Y2Imp" = AddedFrom)  
        result.out = rbind(result.out, Accepted.Moves.Holder)
      }
     
       rm( list = setdiff(ls(),c(STARTSPACE,"STARTSPACE") ) )

    
    } 
  
    if( length(result.out) > 0){
      Accept.Count = nrow(result.out)
    }
      # The first element: which file we are working with (1 or 2)
      # The second tells us which record within that file
      # The third tells us the block that we were in
      # The fourth tells us the block that we are moving to
      # The fifth gives us part of the acceptance probabilities 
      # The sixth, which we haven't called here, tells us the proposed B field that caused the change.
    
    if( Accept.Count > 0){
      # Update Current B 
      holder.jstar =  result.out[,7] 
      current.B[[result.out[1,1]]][ result.out[,2], holder.jstar[1] ] = levels(Bhat[,holder.jstar[1]])[result.out[,6]]
      # Update E 
      if( error.file == 1){      
        E1[result.out[,2],] = Estar[[1]][result.out[,2],]      
      } else{       
        E2[result.out[,2],] = Estar[[2]][result.out[,2],]    
      }  
  
      # Step 1: All at once adapt the block 
      File2_data[ result.out[,2], "Block" ] = result.out[,4]
      File2_data[ result.out[,2], c("lambda","BlockRow") ] = 0
      # Step 2: All at once adapt the block variables
      holder = Imp.Holder2[result.out[,4]]
      holder = do.call(rbind,holder)
      holder = holder[,where.bvars.F2]
      File2_data[ result.out[,2], where.bvars.F2]  = holder
      # Step 3: Remove the records from the blocks 
      # Split the records by FROM block 
      out.records.from = split(result.out[,2],result.out[,3])
      length.out.records.from = mapply(length, out.records.from)
      # Sort the unique FROM blocks
      unique.from = sort(unique(result.out[,3]))
      RIB[[error.file]][ unique.from ] = mapply( setdiff, RIB[[error.file]][unique.from], out.records.from,SIMPLIFY = FALSE )
      # Step 4: Put the new records in the block 
      # Split the records by TO block 
      out.records.to        = split(result.out[,2],result.out[,4])
      length.out.records.to = mapply(length, out.records.to)
      # Sort the unique TO blocks
      unique.to             = sort( unique(result.out[,4]) )
      RIB[[error.file]][ unique.to ] = mapply(append, RIB[[error.file]][ unique.to ], out.records.to, SIMPLIFY = FALSE)
      # Step 5: Make all records that were moved have a lambda value of 0 
      File2_data[result.out[,2],"lambda"] = 0
      
      # Step 6: Delete Imputations 
      # Which records in File 2 are imputations?
      ImpF2Which = mapply( function(x) x[which( x > n2)], RIB[[2]])
      # Which records in File 1 are imputations?
      ImpF1Which = mapply( function(x) x[which( x > n1)], RIB[[1]])
      # How many records in File 1 are imputations?
      l1 = mapply(length,ImpF1Which)
      # How many records in File 2 are imputations?
      l2 = mapply(length,ImpF2Which)
      # How many imputations do we need to remove from each block? 
      sample.from = mapply(min,l1[unique.from],length.out.records.from)
      sample.to   = mapply(min,l2[unique.to],length.out.records.to)
      
      # Remove records from File 2, and from RIB[[2]]
      imps.to.remove = mapply( function(x,y,z) if( y !=0 ){ImpF2Which[[z]][sample(unlist(x),y,replace=F)]}, l2[unique.to], sample.to,unique.to)
      RIB[[error.file]][ unique.to ] = mapply( setdiff, RIB[[error.file]][unique.to], imps.to.remove,SIMPLIFY = FALSE )
      
      if( length(unlist(imps.to.remove)) > 0 ){
         File2_data = File2_data[-unlist(imps.to.remove),]
      }
      
      # Remove records from File 1, and from RIB[[1]]
      imps.from.remove = mapply( function(x,y,z) if( y !=0 ){ImpF1Which[[z]][sample(unlist(x),y,replace=F)]}, l1[unique.from], sample.from,unique.from)
      RIB[[1]][ unique.from ] = mapply( setdiff, RIB[[1]][unique.from], imps.from.remove,SIMPLIFY = FALSE )
      
      if( length(unlist(imps.from.remove)) > 0 ){
        File1_data = File1_data[-unlist(imps.from.remove),]
      }
           
      # We need some kind of static count for the blocks in File 1, and then a dynamic count for the blocks in File 2. 
      # Static: 
      InF1         = mapply( function(x) length(which( x <= n1)), RIB[[1]])
      InF2         = mapply( function(x) length(which( x <= n2)), RIB[[2]])
      ImpF1        = mapply( function(x) length(which( x > n1)), RIB[[1]])
      ImpF2        = mapply( function(x) length(which( x > n2)), RIB[[2]])
      ImpF2Which   = mapply( function(x) which( x > n2), RIB[[2]])
      ImpF1Which   = mapply( function(x) which( x > n1), RIB[[1]])
      Empty.Blocks = which(InF1==0 & InF2==0)
  
     # Subtract the number of records in block from File 2 from the number in File 1 
      to.impute = InF1- InF2
      
      # How many imputations do we need in each file?
      F1 = which( to.impute < 0 )
      F2 = which( to.impute > 0 )
    
      impsneeded.F1 = ImpF1- abs(to.impute) 
      # Negative: We have fewer than we need; add some 
      to.add.F1     = intersect( which(impsneeded.F1< 0), F1)
      
      impsneeded.F2 = ImpF2 - to.impute
      to.add.F2     = intersect( which(impsneeded.F2< 0), F2)  
      
      # Determine which imputations we have stored from result.out 
       Y1Potential = subset(result.out, result.out[,"Y1Imp"] != 0 )   
       Y1Imps      = sort(unique(Y1Potential[,4]))
       Y2Potential = subset(result.out, result.out[,"Y2Imp"] != 0 )  
       Y2Imps      = sort(unique(Y2Potential[,3]))
      # Intersect the blocks which have imputations which the blocks that actually need them. 
       Y1Imps.none = setdiff(to.add.F1,Y1Imps)
       Y1Imps      = intersect(to.add.F1,Y1Imps)
       Y2Imps.none = setdiff(to.add.F2,Y2Imps)
       Y2Imps      = intersect(to.add.F2,Y2Imps)
       result.out1 = subset(result.out,result.out[,4] %in% Y1Imps)
       result.out2 = subset(result.out,result.out[,3] %in% Y2Imps)
      # For the blocks that have useful imputations, count how many we have. 
       Y1ImpsHave = split( result.out1[,"Y1Imp"], result.out1[,4])
       Y2ImpsHave = split( result.out2[,"Y2Imp"], result.out2[,3])
      
       if(length(Y1ImpsHave) > 0){
         count.Y1ImpsHave = mapply(length,Y1ImpsHave)
         
         # Are there any blocks in which I do not have enough? 
         Y1.check = abs(impsneeded.F1[Y1Imps]) - count.Y1ImpsHave
         
         toofew.Y1 = which(Y1.check > 0 )
         
         # Step 1: Determine which of the imputed response values we will keep. 
         sample.to.add.Y1 = mapply( function(x,y) sample( count.Y1ImpsHave[x], y, replace=F), 1:length(Y1Imps), abs(impsneeded.F1[Y1Imps]))
         
         # Step 2: Pull those response variables
         Add1 =  mapply( function(x,y) Y1ImpsHave[[x]][unlist(y)], 1:length(Y1Imps), sample.to.add.Y1 )        
         
         # Step 3: Find the body of the record that we need. Fill in the response variable 
         ToRBindF1 = mapply( function(x,y) rep(Imp.Holder1[[x]],y) , Y1Imps, abs(impsneeded.F1[Y1Imps]) )
         ToRBindF1 = matrix(unlist(ToRBindF1),ncol = ncol(File1_data),byrow =T)
         ToRBindF1[,col.Y1] = unlist(Add1)

         
         # Step 4: Determine how many more imputations are needed.  
         more.neededY1  = length(Y1Imps.none)+length(toofew.Y1)
         
         if( more.neededY1 > 0 ){
           # Determine which blocks need these imputations
           which.neededY1 = c(Y1Imps.none,toofew.Y1)
           # Grab the bcdies of these records 
           MoreRBindF1 = mapply( function(x,y) rep(Imp.Holder1[[x]],y) , which.neededY1, abs(impsneeded.F1[which.neededY1]) )
           MoreRBindF1 = matrix(unlist(MoreRBindF1),ncol = ncol(File1_data),byrow =T)
           # Find the within block medians for these blocks 
           Medians     = mapply( function(x,y) rep(median(File1_data[RIB[[1]][[x]],1]),y) , which.neededY1,abs(impsneeded.F1[which.neededY1]) )
           Medians = unlist(Medians)
           Medians[is.na(Medians)] = median(File1_data[,1])
           MoreRBindF1[,1] = unlist(Medians) 
           ToRBindF1 =rbind(ToRBindF1,MoreRBindF1)
         }
         
         # Step 5: Tack on all imputed records 
         colnames(ToRBindF1) = colnames(File1_data) 
         end.F1              = nrow(File1_data)+1
         File1_data          = rbind(File1_data, ToRBindF1)
       }
      
      if(length(Y1ImpsHave) ==0 & length(Y1Imps.none)>0){
        
          which.neededY1 = Y1Imps.none
          # Grab the bcdies of these records 
          MoreRBindF1 = mapply( function(x,y) rep(Imp.Holder1[[x]],y) , which.neededY1, abs(impsneeded.F1[which.neededY1]) )
          MoreRBindF1 = matrix(unlist(MoreRBindF1),ncol = ncol(File1_data),byrow =T)
          # Find the within block medians for these blocks 
          Medians     = mapply( function(x,y) rep(median(File1_data[RIB[[1]][[x]],1]),y) , which.neededY1,abs(impsneeded.F1[which.neededY1]) )
          Medians = unlist(Medians)
          Medians[is.na(Medians)] = median(File1_data[,1])
          MoreRBindF1[,1] = unlist(Medians) 
          ToRBindF1 = MoreRBindF1        
          # Step 5: Tack on all imputed records 
          colnames(ToRBindF1) = colnames(File1_data) 
          end.F1              = nrow(File1_data)+1
          File1_data          = rbind(File1_data, ToRBindF1)
      }
      
      if(length(Y2ImpsHave) > 0 ){
        count.Y2ImpsHave = mapply(length,Y2ImpsHave)
        
        # Are there any blocks in which I do not have enough? 
        Y2.check = abs(impsneeded.F2[Y2Imps]) - count.Y2ImpsHave
        
        toofew.Y2 = which(Y2.check > 0 )
        
        # Step 1: Determine which of the imputed response values we will keep.        
        sample.to.add.Y2 = mapply( function(x,y) sample( count.Y2ImpsHave[x], y, replace=F), 1:length(Y2Imps), abs(impsneeded.F2[Y2Imps]))
        
        # Step 2: Pull those response variables
        Add2 =  mapply( function(x,y) Y2ImpsHave[[x]][unlist(y)], 1:length(Y2Imps), sample.to.add.Y2 )
        
        
        # Step 3: Find the body of the record that we need. Fill in the response variable 
        ToRBindF2 = mapply( function(x,y) rep(Imp.Holder2[[x]],y) , Y2Imps, abs(impsneeded.F2[Y2Imps]) )
        ToRBindF2 = matrix(unlist(ToRBindF2),ncol = ncol(File2_data),byrow =T)
        ToRBindF2[,col.Y2] = unlist(Add2)
        
        # Step 4: Determine how many more imputations are needed.  
        more.neededY2  = length(Y2Imps.none)+length(toofew.Y2)
        
        
        if( more.neededY2 > 0 ){
          # Determine which blocks need these imputations
          which.neededY2 = c(Y2Imps.none,toofew.Y2)
          # Grab the bcdies of these records 
          MoreRBindF2 = mapply( function(x,y) rep(Imp.Holder2[[x]],y) , which.neededY2, abs(impsneeded.F2[which.neededY2]) )
          MoreRBindF2 = matrix(unlist(MoreRBindF2),ncol = ncol(File2_data),byrow =T)
          # Find the within block medians for these blocks 
          Medians     = mapply( function(x,y) rep(median(File2_data[RIB[[2]][[x]],1]),y) , which.neededY2,abs(impsneeded.F2[which.neededY2]) )
          Medians = unlist(Medians)
          # For any Medians that are NA, replace them with the global median
          Medians[is.na(Medians)] = median(File2_data[,1])      
          MoreRBindF2[,1] = unlist(Medians) 
          ToRBindF2 = rbind(ToRBindF2, MoreRBindF2)
        }
        
        # Step 5: Tack on all imputed records 
        colnames(ToRBindF2) = colnames(File2_data)
        end.F2              = nrow(File2_data)+1
        File2_data          = rbind(File2_data, ToRBindF2)
      }
      
      if(length(Y2ImpsHave) == 0 & length(Y2Imps.none) > 0 ){
        
          # Determine which blocks need these imputations
          which.neededY2 = Y2Imps.none
          # Grab the bcdies of these records 
          MoreRBindF2 = mapply( function(x,y) rep(Imp.Holder2[[x]],y) , which.neededY2, abs(impsneeded.F2[which.neededY2]) )
          MoreRBindF2 = matrix(unlist(MoreRBindF2),ncol = ncol(File2_data),byrow =T)
          # Find the within block medians for these blocks 
          Medians     = mapply( function(x,y) rep(median(File2_data[RIB[[2]][[x]],1]),y) , which.neededY2,abs(impsneeded.F2[which.neededY2]) )
          Medians = unlist(Medians)
          # For any Medians that are NA, replace them with the global median
          Medians[is.na(Medians)] = median(File2_data[,1])      
          MoreRBindF2[,1] = unlist(Medians) 
          ToRBindF2 = MoreRBindF2
            
        # Step 5: Tack on all imputed records 
        colnames(ToRBindF2) = colnames(File2_data)
        end.F2              = nrow(File2_data)+1
        File2_data          = rbind(File2_data, ToRBindF2)
      }

      nrow.F1             = nrow(File1_data)
      nrow.F2             = nrow(File2_data)
      
      if( nrow.F1 != nrow.F2){
        stop("The files do not have the same number of records")
      }
  
      # Step 6: update RIB 
      # We are going to try only updating this once,at the beginning
      #RIB1Base  = mapply( function(x) x[which(x %in% type1seeds)], RIB[[1]])
      #RIB2Base  = mapply( function(x) x[which(x %in% type1seeds)], RIB[[2]])
      blocks1   = sort(unique(File1_data[-type1seeds,"Block"]))
      blocks2   = sort(unique(File2_data[-type1seeds,"Block"]))
      InF1      = split( c(1:nrow.F1)[-type1seeds] , File1_data[-type1seeds,"Block"])
      InF2      = split( c(1:nrow.F2)[-type1seeds] , File2_data[-type1seeds,"Block"])
      RIB[[1]][blocks1] = mapply(union,RIB1Base[blocks1],InF1)
      RIB[[2]][blocks2] = mapply(union,RIB2Base[blocks2],InF2)
      n.RIB[1,] = mapply(length,RIB[[1]])
      n.RIB[2,] = mapply(length,RIB[[2]])
  
      # Step 7: In row 1, column 1, we see the number of imputed records in file 1 for that block. 
      ImpF2Which      = mapply( function(x) which( x > n2), RIB[[2]])
      ImpF1Which      = mapply( function(x) which( x > n1), RIB[[1]])
      # Step 8: Identify the number of records in each block which are imputations 
      Imputed.Block = NULL 
      Imputed.Block = cbind(Imputed.Block,unlist(mapply( length, ImpF1Which)))
      Imputed.Block = cbind(Imputed.Block,unlist(mapply( length, ImpF2Which)))
  
      # Step 9: Update Lambda
      File1_data[,'lambda'] = 1:nrow.F1
      # Note that we have not updated the File 2, so we are going to have to use the File 1 stuff. 
  
      # Sample C : 
  
       # Take the current ordering of the File 2 data 
       F1.Lambda   = 1:nrow.F1
       F2.Lambda   = File2_data[,"lambda"]
    
       max.RIB = apply( n.RIB, 2, max)
       max.RIB = max.RIB - count.matches
       non.seed.RIB = max.RIB
    
       PermSampling = ifelse( max.RIB > 1, ifelse(max.RIB < threshold, 1, 2), 0 )
       SeedsOnly    = which(max.RIB==0)
       # This vector takes a 1 if we use exact sampling, 2 if we use MH and 0 otherwise.
       to.sampleC   = which(PermSampling>0)
       to.sampleC   = setdiff( to.sampleC, Empty.Blocks) # This SHOULD be superfluous, but is an okay check
       
       # The only blocks which should have BlockRow issues are those that had blocks moved in or out. 
       BlockRowCheck = setdiff( sort(unique(result.out[,3:4])), Empty.Blocks)
       BlockRow.F1   = File1_data[,"BlockRow"]
       BlockRow.F2   = File2_data[,"BlockRow"]
       BlockRow.F2[result.out[,2]] = 0
     
       # For these blocks, we may need to fix the BlockRow. The records that were in the correct block already will have numbers that make sense. Everyone else will have 0s. We just need to do setdiff, and fill in the missing values. 
       HOLDER    = mapply( function(x,y) setdiff(1:x,BlockRow.F1[unlist(y)]), n.RIB[1,BlockRowCheck],RIB[[1]][BlockRowCheck])
       NEWHOLDER = mapply( function(x,y) x[which(BlockRow.F1[unlist(x)]>y)], RIB[[1]][BlockRowCheck], n.RIB[1,BlockRowCheck])
       BlockRow.F1[unlist(NEWHOLDER)] = 0
     
       # These are for the ones that are 0. 
       NEWHOLDER = mapply( function(x) x[which(BlockRow.F1[unlist(x)]==0)], RIB[[1]][BlockRowCheck])
       NEWHOLDER = unlist(NEWHOLDER)
       BlockRow.F1[NEWHOLDER] = unlist(HOLDER)
     
       # File 2 
       HOLDER    = mapply( function(x,y) setdiff(1:x,BlockRow.F2[unlist(y)]), n.RIB[2,BlockRowCheck],RIB[[2]][BlockRowCheck])
       NEWHOLDER = mapply( function(x,y) x[which(BlockRow.F2[unlist(x)]>y)], RIB[[2]][BlockRowCheck], n.RIB[2,BlockRowCheck])
       BlockRow.F2[unlist(NEWHOLDER)] = 0

      # These are for the ones that are 0. 
       NEWHOLDER            = mapply( function(x) x[which(BlockRow.F2[unlist(x)]==0)], RIB[[2]][BlockRowCheck])
       NEWHOLDER            = unlist(NEWHOLDER)
       #NEWHOLDER2 = unlist(mapply( function(x) x[sample(length(x))], HOLDER))
       BlockRow.F2[NEWHOLDER]  = unlist(mapply( function(x) x[sample(length(x))], HOLDER))
       # We also need to give these guys lambdas, since their lambda values will also be 0. 
       F2.LambdaNew         = mapply( function(x,y)  F1.Lambda[x][BlockRow.F2[y]], RIB[[1]], RIB[[2]])
     
       # Score! Now, we just have to put them back in order. This is a list which is block specific. We need them back in the actual order of the records. I guess I just unlist RIB[[2]] 
       F2.Lambda[unlist(RIB[[2]])]= unlist(F2.LambdaNew)
    }
    }
    }

     # Take the current ordering of the File 2 data 
     Y1          = File1_data[,col.Y1]
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
      order.in   = BlockRow.F2[F2.records]
      # How many pairs are in the block?
      pairs      = n.RIB[1,k]
      # How many type 1 seeds?
      if( class(type1seeds) != "NULL"){
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
      BlockRow.F2[F2.records]    = block.order.out 
    } 
    
    # Reorder based on the sampled C 
    Y2                  = Y2.in[order(New.Ordering)]
    V.Holder[,col.Y2+1] = Y2
    F2.Lambda           = New.Ordering 
    
    # Store in the File 2 data space 
    File2_data[,"lambda"]   = New.Ordering
    File2_data[,"BlockRow"] = BlockRow.F2 
    File1_data[,"BlockRow"] = BlockRow.F1

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
    
    imputed.F1   = n1:nrow.F1
    n.imputed.F1 = length(imputed.F1)
    
    imputed.F2   = n2:nrow.F2
    n.imputed.F2 = length(imputed.F2)    
    
    V2 = cbind( 1, File2_data[, where.coefs.secondary ] )
    V2 = as.matrix(V2)
    
    V2TV2        = crossprod(V2)
    Gen          = zapsmall( eigen( V2TV2 )$values )
    Gen          = ifelse( sum(Gen==0) > 0 , "TRUE","FALSE")
    invert.part2  = invert_XTX( V2TV2 , Gen )
    
    t.V2         = t(V2)
    Proj2        = invert.part2%*%t.V2
    rm(V2TV2,Gen,t.V2)
    
    # Added for the speed up in matrix inversion 
    Bindex = c(1:nrow(File2_data))[-c(type1seeds)]
    
    holder = Gupdate_and_Impute_Speed(Y1,Y2, Y2.in, V.Holder, New.Ordering,Sigma1.Part,Sigma2.Part,which.Y1,which.Y2, imputed.F2, imputed.F1,n.imputed.F2,n.imputed.F1,secondary.option="dependent",Proj2,V2,invert.part2,C,Bindex,Ainv)
    
  
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
      write( c(beta.star,sigmasq), file = paste(namepart,"ThetaOut.txt", sep=""), append = TRUE)
      write( c(eta,sigma2sq), file = paste(namepart,"EtaOut.txt", sep=""), append = TRUE)
      write( c(gamma.samp[error.file,c(5)]), file = paste(namepart,"GammaOut2.txt", sep=""), append = TRUE)
      #write( c(gamma.samp[error.file,c(6)]), file = paste(namepart,"GammaOut2.txt", sep=""), append = TRUE)
      Lambda.out <- sum(F2.Lambda[1:n2]==1:n2)
      write( Lambda.out, file = paste(namepart,"LambdaOut.txt", sep=""), append = TRUE)
      Lambda.out <-Lambda.out/nrow.F2
      write(Lambda.out, file = paste(namepart,"AllLambdaPerc.txt",sep=""),append=TRUE)
      write(nrow.F2,file=paste(namepart,"SizeOut.txt",sep=""),append=TRUE)      
      write( Accept.Count, file = paste(namepart,"AcceptOut.txt", sep=""), append = TRUE)
      Accept.Count = 0 
      # If we are interested in the completed data sets 
      if( s%%50 == 0){
        write.csv(File2_data,file=paste(namepart,"F2DataOut",s,".csv",sep=""),row.names=F)
      }
     
    }
    
  }
  time.final = proc.time()-sys
  }