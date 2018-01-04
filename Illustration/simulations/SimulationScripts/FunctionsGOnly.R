# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
####    Functions: Gutman Only            ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

invert_XTX <-function(M,Gen){  
  if(Gen=="FALSE"){  
    return(chol2inv(chol(M) ))
  } else {
    return(ginv(M) )
  }
  
}

rigamma <-function (n, alpha, beta) {
  if (alpha > 0 & beta > 0) 
    1/rgamma(n = n, alpha, beta)
  else stop("rigamma: invalid parameters\n")
}

mvrnorm<-function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE, EISPACK = FALSE) {
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  if (missing(EISPACK)) 
    EISPACK <- getOption("mvnorm_use_EISPACK", FALSE)
  eS <- eigen(Sigma, symmetric = TRUE, EISPACK = EISPACK)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1) 
    drop(X)
  else t(X)
}

permutations <-function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE){  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 
                                                                      0) 
  stop("bad value of n")
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 
        0) 
    stop("bad value of r")
  if (!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if ((r > n) & repeats.allowed == FALSE) 
    stop("r > n and repeats.allowed=FALSE")
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n) 
      stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  if (repeats.allowed) 
    sub <- function(n, r, v) {
      if (r == 1) 
        matrix(v, n, 1)
      else if (n == 1) 
        matrix(v, 1, r)
      else {
        inner <- Recall(n, r - 1, v)
        cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner), 
                                                  ncol = ncol(inner), nrow = nrow(inner) * n, 
                                                  byrow = TRUE))
      }
    }
  else sub <- function(n, r, v) {
    if (r == 1) 
      matrix(v, n, 1)
    else if (n == 1) 
      matrix(v, 1, r)
    else {
      X <- NULL
      for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n - 
                                                        1, r - 1, v[-i])))
      X
    }
  }
  sub(n, r, v[1:n])
}

dmvnorm <-function (x, mean, sigma, log = FALSE) 
{
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  if (missing(mean)) {
    mean <- rep(0, length = ncol(x))
  }
  if (missing(sigma)) {
    sigma <- diag(ncol(x))
  }
  if (NCOL(x) != NCOL(sigma)) {
    stop("x and sigma have non-conforming size")
  }
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  distval <- mahalanobis(x, center = mean, cov = sigma)
  logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
  logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
  if (log) 
    return(logretval)
  exp(logretval)
}

posterior_summaries<-function(object){
  apply(object,2,function(x) mean(x) + 1.96*sd(x))
  apply(object,2,function(x) mean(x) - 1.96*sd(x))
  cat("Posterior CI Lower: ",  apply(object,2,function(x) mean(x) - 1.96*sd(x)), "\n" )
  cat("Posterior Means: ",  c(apply(object,2,mean)), "\n")
  cat("Posterior CI Upper: ",  apply(object,2,function(x) mean(x) + 1.96*sd(x)), "\n" )
  cat("Posterior SE: ",  apply(object,2,sd),  "\n ") 
}

posterior_summaries_neat <-function (object,name.vector,truth) {
  
  length.out = length(name.vector)
  
  for( i in 1:length.out){
    cat("The posterior estimate for", name.vector[i] ,"is", mean(object[,i]), " with a CI of (", (mean(object[,i]) - 1.96 * sd(object[,i]) ) ,",",( mean(object[,i])+1.96 * sd(object[,i])),") with a SE of", sd(object[,i]),"\n")
    cat("The true value of", truth[i] , ifelse(truth[i] < (mean(object[,i])+1.96 * sd(object[,i])) & truth[i] > (mean(object[,i])-1.96 * sd(object[,i])),"is","is not" ) , "in the CI. \n")
  }
  
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# Purpose: 
# Runs the Gutman MCMC

# Returns: 
# Gout.list - a list which contains the regression parameters samples from the MCMC
# Saves out.put as .RData files

gutman_mcmc <-function(is.test,seed, its, burnin, thinning, gap, n.RIB, RIB, File1_data, File2_data,K,which.Y2, which.Y1,type1Seeds, where.coefs.primary, where.coefs.secondary, p,p2, n.type1,needs.C,reps, PermSampling, secondary.option, n1, n2,col.Y1, col.Y2,namepart){
  
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
  
  # We begin by using ONLY the seeds information, if there are seeds. Otherwise, we use everything. 
  
  if( class(type1Seeds) != "NULL"){
  
    Y1 = File1_data[type1Seeds,col.Y1]
    Y2 = File2_data[type1Seeds,col.Y2]
  
    # Create the model matrix based on the seeds
    Vtheta = cbind( 1, File2_data[type1Seeds, where.coefs.primary ] )
    
    Sigma1.Part = (n.type1 - p )/2
    
  } else{ 
      Y1 = File1_data[,col.Y1]
      Y2 = File2_data[order(File2_data[,"lambda"]),col.Y2]
  
    # Create the model matrix based on the seeds
    Vtheta = cbind( 1, File2_data[order(File2_data[,"lambda"]), where.coefs.primary ] )
    
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
  sqrt.sigma   = sqrt(sigmasq)
  #cat( "The variance is", sigmasq,"\n")
  var.beta     = sigmasq*invert.part
  #cat( "The primary variance is",diag(var.beta),"\n")
  #print( dim(var.beta) )
  beta         = mvrnorm( 1, hat.beta, var.beta ) 
  #cat( "The beta mean is",g1/(1+g1)*hat.beta + 1/(g1+1)*beta0,"\n")
  cat( "The inital beta draw is",beta,"\n")
  
  rm( Sigma1.Part,var.beta,ResidSS,hat.beta, Proj, t.Vtheta,Gen,invert.part,VthetaTVtheta)
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ##########  Initialization : Y2 Coefs  ##########
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  if( secondary.option == "independent"){
    
    # The Y2s are not dependent on the blocking variables 
    
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
    sqrt.sigma2  = sqrt(sigma2sq)
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
  
  Sigma1.Part = (nrow.F1 - p )/2
  Sigma2.Part = (nrow.F1 - p2 )/2
  
  V2 = cbind( 1, File2_data[, where.coefs.secondary ] )
  V2 = as.matrix(V2)
  
  V2TV2        = crossprod(V2)
  Gen          = zapsmall( eigen( V2TV2 )$values )
  Gen          = ifelse( sum(Gen==0) > 0 , "TRUE","FALSE")
  invert.part2  = invert_XTX( V2TV2 , Gen )
  
  t.V2         = t(V2)
  Proj2        = invert.part2%*%t.V2
  
  rm(V2TV2,Gen,t.V2)
  
  BlockRow = File2_data[,"BlockRow"]
  
  # Take the current ordering of the File 2 data 
  F1.Lambda  = 1:nrow.F1
  F2.Lambda  = File2_data[,"lambda"]
  if( class(type1Seeds) != "NULL"){
    F2.Lambda[type1Seeds] = type1Seeds
  }
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
    
    holder = Gupdate_and_Impute(Y1, Y2, Y2.in,V.Holder, New.Ordering, Sigma1.Part,Sigma2.Part,which.Y1,which.Y2,imputed.F2,imputed.F1,n.imputed.F2,n.imputed.F1,secondary.option, Proj2,V2,invert.part2)
    
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

# %%%%%%%%%%%%%%%%%%% # 

# Purpose: 
# Exactly samples a new permutation for small blocks 

# Returns: 
# order.out = a vector giving the new ordering for the File 2 data in the block 

G.exactSampleC <-function( yData, orderC, Design.Matrix,num.pairs, beta, sqrt.sigma,block){
  
  switch.options = orderC
  num.switch = num.pairs 
  order.out <- c()
  
  AllPossibleC = permutations( n =num.switch, r = num.switch, v= switch.options )
  
  Num.PossibleC = factorial(num.switch)
  
  likelihood.storage = rep(0, Num.PossibleC )   
  
  Variance.Component = sqrt.sigma^2*diag(num.switch) 
  
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
  
  #A.ratio.contribution = log( likelihood.storage[SampledC] )
  
  #new.list <-list( "order.out" = order.out, "Sum" = Sum.Part) 
  
  return( order.out ) 
  
}


#%%%%%%%%%%%%%%%%%%%% #

# Purpose: 
# Performs the Gutman MH step within the new model- ie if the block is small, it exactly samples a permutation; if the block is large and it selects two records from File 2 and switches their position (repeats this step reps times)

# Returns: 
# order.out = a vector giving the new ordering for the File 2 data in the block

SamplePerm <- function(F1.holder,F2.holder, pairs, order.in, beta,sqrt.sigma,S.ind,reps,block){  
  
  # If we get here, we know that we need to sample a permutation
  Design.Holder = cbind( 1, F2.holder)
  
  # S.ind takes a 1 if we use exact sampling, 2 if we use MH
  
  if( S.ind ==1 ){
    # Here we use exact sampling 
    order.out  = G.exactSampleC( F1.holder, 1:pairs, Design.Holder, pairs, beta, sqrt.sigma, block)
    
  } else{ 
    
    # Order the Design Holder to match the original ordering 
    Design.Holder = Design.Holder[order.in,]
    
    orderC = 1:pairs
    switch.options = orderC
  
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

#%%%%%%%%%%%%%%%%%%%% #

# Purpose: 
# Updates the regression parameters and re-imputes missing values 

# Returns: 
# list( "Y1" = Y1, "Y2" = Y2, "beta" = beta, "sigmasq" = sigmasq, "sigma" = sigma, "mu2" = mu2, "sigma2sq"=sigma2sq, "theta" =theta.small, "impVar" = sqrt(var.Y2impute))
# The re-imputed data and the regression parameters 

Gupdate_and_Impute<-function(Y1, Y2, Y2.in, Y2.info, New.Ordering, Sigma1.Part,Sigma2.Part,which.Y1,which.Y2,imputed.F2,imputed.F1,n.imputed.F2,n.imputed.F1,secondary.option, Proj2,V2,invert.part2){  
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ###########  Initialization  Theta #############
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###  
  
  # Create the model matrix 
  Vtheta       = Y2.info 
  
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
    ResidSS      = t(ResidSS)%*%ResidSS
    ResidSS      = .5* ResidSS  
    
    sigma2sq     = rigamma(1,Sigma2.Part,ResidSS) 
    sqrt.sigma2  = sqrt(sigma2sq)
    var.eta      = sigma2sq*invert.part2
    eta          = mvrnorm( 1, hat.eta, var.eta ) 
    
    V2.imp       = V2[ imputed.F2 , ]%*%eta 
  } 
  
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
  ####### Initialization : Y1 Imputations #########
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

#%%%%%%%%%%%%%%%%%%%% #

# Purpose: 
# Obtains the starting permutation when we are simulating data; this is NOT necessary if the data is being input 

# Returns: 
# holder$File1, holder$File2 (the data with the inital block rows)

G.get_starting_permutation <-function(K, RIB, File1_data, File2_data,n1,n2,type1seeds ){
  
  holder1 = rep( 0, n1 ) 
  holder2 = rep( 0, n2 )
  lambdaholder = rep( 0, n2 )
  
  for(k in 1:K){
    
    #Step 1: Look at the records in the block
    in.block   = list( RIB[[1]][[k]], RIB[[2]][[k]] )
    max.length = max( length(in.block[[1]]),length(in.block[[2]]) )
    
    # If the length is not 0, proceed. This means that we have some elements in the block.
    if( max.length >0 ){
      
      # Fill in the information from file 1
      places1 = 1:length(in.block[[1]])
      holder1[in.block[[1]]] = places1 
      #Fill in the information from file 2
      places2 = permute.places( 1:length(in.block[[2]]) )
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
    
    
    if(length(unique(lambdaholder[ which(lambdaholder> 0) ])) != length( lambdaholder[ which(lambdaholder> 0)]) ) {
      print(k)
      break 
    }
  }
  
  File1_data[,"BlockRow"] = holder1 
  File2_data[,"BlockRow"]= holder2
  File2_data[,"lambda"]   = lambdaholder
  File2_data[type1seeds,"lambda"] = File1_data[type1seeds,"lambda"]
  
  newlist <- list( "File1" = File1_data, "File2" = File2_data )
  
  return( newlist )
}

# %%%%%%%%%%%%%%%%%%%% # 

# Purpose: 
# Identify the records in each block; this excludes the Type 1 Seeds. For the Gutman algorithm, we do not care which blocks the type 1 seeds belong to 

# Returns: 
# RIB[[1]] and RIB[[2]] 
# Also returns n.RIB, a matrix which counts the number of records from each file (row) in each block (column)

G.records.in.block <- function(File1_data, File2_data, num.blocks,type1seeds){
  
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
      holder        = setdiff(holder,type1seeds)
      RIB[[f]][[k]] = holder
      n.RIB[f,k]    = length(holder)
    }
  }
  new.list <- list( "RIB" = RIB, "n.RIB" = n.RIB )
  return(new.list)
}

