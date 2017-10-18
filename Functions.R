#### %%%%%%%%%%%%%%%%%%%%%%%%% ## 
####  Functions: DataGen.R  #### 
#### %%%%%%%%%%%%%%%%%%%%%%%%% ## 

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

## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
#####  Functions: PB Required #### 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

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

permute.places<-function(vec1){
  n=length(vec1)
  vec1.ordering=sample(1: n , n ,replace=FALSE)
  vec2=vec1[order(vec1.ordering)];vec2
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

Gupdate_and_Impute_Speed <- function (Y1, Y2, Y2.in, Y2.info, New.Ordering, Sigma1.Part, Sigma2.Part, which.Y1, which.Y2, imputed.F2, imputed.F1, n.imputed.F2, n.imputed.F1, secondary.option, Proj2, V2,invert.part2, C, Bindex, Ainv) {
  
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
  ####                  Theta                       ####
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
  
  # Create the model matrix 
  B      = Y2.info[Bindex, ]
  Vtheta = Y2.info
  B      = crossprod(B)
  
  invert.part = Ainv -Ainv%*%B%*%solve(C + Ainv%*%B)%*%Ainv
  
  Proj = tcrossprod( invert.part,Vtheta)
  hat.beta = Proj%*%Y1

  ResidSS = (Y1 - Vtheta %*% hat.beta)
  ResidSS = crossprod(ResidSS)
  ResidSS = 0.5 * ResidSS
  
  sigmasq = rigamma(1, Sigma1.Part, ResidSS)
  sigma = sqrt(sigmasq)
  var.beta = sigmasq * invert.part
  beta = mvrnorm(1, hat.beta, var.beta)
    if (secondary.option == "independent") {
        barY2 = mean(Y2)
        Var.Part = sum((Y2 - barY2)^2)
        sigma2sq = rigamma(1, Sigma2.Part, 0.5 * Var.Part)
        eta = rnorm(1, barY2, (1/n.type1) * sigma2sq)
    }
    else {
        hat.eta = Proj2 %*% Y2.in
        est.mean.Y2 = V2 %*% hat.eta
        ResidSS = (Y2.in - est.mean.Y2)
        ResidSS = crossprod(ResidSS)
        ResidSS = 0.5 * ResidSS
        
        sigma2sq = rigamma(1, Sigma2.Part, ResidSS)
        sqrt.sigma2 = sqrt(sigma2sq)
        var.eta = sigma2sq * invert.part2
        eta = mvrnorm(1, hat.eta, var.eta)
        
        V2.imp = V2[imputed.F2, ] %*% eta
    }
  
    F1.matched.to.imputations = New.Ordering[imputed.F2]
    Y1.matched.to.imputations = Y1[F1.matched.to.imputations]
    V1 = Vtheta[F1.matched.to.imputations, -which.Y2]
  
    beta.impute = beta[-which.Y2]
    theta.small = beta[which.Y2]
  
    mean.Y2impute = V1 %*% as.matrix(beta.impute) - Y1.matched.to.imputations
    mean.Y2impute = sigmasq^(-1) * theta.small * mean.Y2impute
  
    if (secondary.option == "independent") {
        mean.Y2impute = sigma2sq^(-1) * eta - mean.Y2impute
    }
    else {
        mean.Y2impute = sigma2sq^(-1) * V2.imp - mean.Y2impute
    }
  
    var.Y2impute = (sigma2sq^(-1) + sigmasq^(-1) * (theta.small^2))^(-1)
    mean.Y2impute = var.Y2impute * mean.Y2impute
  
    imputed.all = rnorm(n.imputed.F2, mean.Y2impute, sqrt(var.Y2impute))
  
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
    ####                  Y1 Imptuations               ####
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
  
    Y2.in[imputed.F2] = imputed.all
    Vtheta[F1.matched.to.imputations, which.Y2] = imputed.all
    Vtheta = Vtheta[imputed.F1, ] 
  
    imputed.all = rnorm(n.imputed.F1, Vtheta %*% as.matrix(beta), 
        sqrt(sigmasq))
    Y1[imputed.F1] = imputed.all
  
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
    ####          Store and Release Output             ####
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
  
    newlist <- list("Y1" = Y1, "Y2" = Y2.in, "beta" = beta, "sigmasq" = sigmasq, 
        "sigma" = sigma, "eta" = eta, "sigma2sq" = sigma2sq, "theta" = theta.small, 
        "impVar" = sqrt(var.Y2impute))
    return(newlist)
}

run_GM <-function(seed, its, burnin, thinning, gap, n.RIB, RIB, File1_data, File2_data,K,which.Y2, which.Y1,type1Seeds, where.coefs.primary, where.coefs.secondary, p,p2, n.type1,needs.C,reps, PermSampling, secondary.option, n1, n2,col.Y1, col.Y2,namepart){
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  #%%%%%%%%   Constants Chunk   %%%%%%%%%%%%##
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
  ############ Data Storage Chunk  %%%%%%%%%%%%##
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
  ###########  Initialization  Theta  %%%%%%%%%%%%##
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  
  # We begin by using ONLY the seeds information, if there are seeds. Otherwise, we use everything. 
  
  if( class(type1Seeds) != "NULL"){
    

    Y1 = File1_data[type1Seeds,col.Y1]
    Y2 = File2_data[type1Seeds,col.Y2]
  
    # Create the model matrix based on the seeds
    Vtheta = cbind( 1, File2_data[type1Seeds, where.coefs.primary ] )
    Vtheta = as.matrix(Vtheta)
    A      = crossprod(Vtheta)
    Ainv   = solve(A)
    
    Sigma1.Part = (n.type1 - p )/2
    
  } else{ 
      Y1 = File1_data[,col.Y1]
      Y2 = File2_data[order(File2_data[,"lambda"]),col.Y2]
  
    # Create the model matrix based on the seeds
    Vtheta = cbind( 1, File2_data[order(File2_data[,"lambda"]), where.coefs.primary ] )
    Vtheta       = as.matrix(Vtheta)
    Sigma1.Part = (nrow.F1 - p )/2
    A      = crossprod(Vtheta)
    Ainv   = solve(A)
  }
  
  Proj = tcrossprod(Ainv,Vtheta)
  hat.beta     = Proj%*%Y1
  
  ResidSS      = (Y1 - Vtheta%*%hat.beta) 
  ResidSS      = crossprod(ResidSS)
  ResidSS      = .5* ResidSS
  
  sigmasq      = rigamma(1,Sigma1.Part,ResidSS) 
  sqrt.sigma   = sqrt(sigmasq)
  var.beta     = sigmasq*Ainv
  beta         = mvrnorm( 1, hat.beta, var.beta ) 
  #cat( "The beta mean is",g1/(1+g1)*hat.beta + 1/(g1+1)*beta0,"\n")
  cat( "The inital beta draw is",beta,"\n")
  
  rm( Sigma1.Part,var.beta,ResidSS,hat.beta, Proj)
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ##########  Initialization : Y2 Coefs   %%%%%%%##
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
    
    #t.V2         = t(V2)
    Proj         = tcrossprod(invert.part,V2)
    hat.eta      = Proj%*%Y2
    est.mean.Y2  = V2%*%hat.eta
    
    ResidSS      = (Y2 - est.mean.Y2) 
    ResidSS      = crossprod(ResidSS)
    ResidSS      = .5* ResidSS
    
  
    sigma2sq     = rigamma(1,Sigma2.Part,ResidSS) 
    sqrt.sigma2  = sqrt(sigma2sq)
    var.eta      = sigmasq*invert.part
    eta          = mvrnorm( 1, hat.eta, var.eta ) 
    cat( "The inital eta draw is",eta,"\n")  
    
    
    rm( Sigma2.Part,var.eta,ResidSS,hat.eta, Proj,Gen,invert.part,V2TV2,est.mean.Y2) 
    
  } 
  
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ####### Initialization : Y2 Imputations  %%%%%%##
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
  ####### Initialization : Y1 Imputations  %%%%%%##
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
    ##########  Initialization : Deletions   %%%%%%##
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
    
    # Delete everything that we no longer need.   
    
    rm( Vtheta,imputed.all,mean.Y2impute, var.Y2impute,beta.impute,theta.small,V1,Y1.matched.to.imputations,F1.matched.to.imputations,Y1,Y2)
  }
  
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
  ############    Run the MCMC     %%%%%%%%%%%%##
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%### 
  
  Bindex = c(1:nrow(File2_data))[-c(type1seeds)]
  C      = diag(p)
  
  Sigma1.Part = (nrow.F1 - p )/2
  Sigma2.Part = (nrow.F1 - p2 )/2
  
  V2 = cbind( 1, File2_data[, where.coefs.secondary ] )
  V2 = as.matrix(V2)
  
  V2TV2        = crossprod(V2)
  Gen          = zapsmall( eigen( V2TV2 )$values )
  Gen          = ifelse( sum(Gen==0) > 0 , "TRUE","FALSE")
  invert.part2  = invert_XTX( V2TV2 , Gen )
  
  #t.V2         = t(V2)
  Proj2        = tcrossprod(invert.part2,V2)
  
  rm(V2TV2,Gen)
  
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

      write( c(beta,sigmasq), file = paste(namepart,"ThetaOut.txt", sep=""), append = TRUE)
      write( c(eta,sigma2sq), file = paste(namepart,"EtaOut.txt", sep=""), append = TRUE)
      Lambda.out <- sum(New.Ordering[1:n2]==1:n2)
      write( Lambda.out, file = paste(namepart,"LambdaOut.txt", sep=""), append = TRUE)
      Lambda.out <-Lambda.out/nrow.F2
      write(Lambda.out, file = paste(namepart,"AllLambdaPerc.txt",sep=""),append=TRUE)

      
    }
    
    if(s%%printgap == 0){ print(s) }
    
  }
  
}

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

run_PB<-function( MCMCSeed, File1, File2, method, its, burnin, gap, reps, thinning, printgap, Comment, coefs.primary,coefs.primary.all,coefs.secondary,col.Y1, col.Y2, secondary.option, secondary.option2, error.file, error.fields, blockingVars,matches,Y1.name,Y2.name,namepart,threshold,complete.out,type2seeds){
  
  realmin =  1e-08
  
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

  coefs.all = union(coefs.primary.all,coefs.all)
  
  # Identify which of these are primary
  where.coefs.primary = which(coefs.all %in% coefs.primary.all)
  # To count the columns, consult the levels associated with each. 
  where.coefs.primary = where.coefs.primary[1]:(where.coefs.primary[1]+sum(d[coefs.primary])-length(coefs.primary))
  #where.coefs.primary = 1:8
  
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
  save.image("PB_PreRunStep.RData")  
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
  ####      Pre Run Steps       ####
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%% ###  
  
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
    run_GM("no",MCMCSeed, its,burnin, thinning, gap, n.RIB, RIB, File1_data, File2_data,K,which.Y2, which.Y1,type1seeds, where.coefs.primary, where.coefs.secondary, p, p2, n.type1,needs.C,reps,PermSampling,secondary.option, n1, n2,col.Y1, col.Y2,namepart)
    print("PB run completed")
    

}

## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
#####  Functions: GM Required #### 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

ginv<-function (X, tol = sqrt(.Machine$double.eps)) 
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


## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
#####  Functions: BLASE Required #### 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

run_BLASE <-function(){
  
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
    
    C = diag(p)
    
    
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
  
  if( class(type1seeds) != "NULL"){
    
    Type1SeedIndicator = 1
    
    Y1 = File1_data[type1seeds,col.Y1]
    Y2 = File2_data[type1seeds,col.Y2]
    
    # Create the model matrix based on the seeds
    Vtheta = cbind( 1, File2_data[type1seeds, where.coefs.primary ] )
    A    = crossprod(Vtheta)
    Ainv = solve(A)
    
    Sigma1.Part = (n.type1 - p )/2
    
  } else{ 
    Type1SeedIndicator = 0
    
    Y1 = File1_data[,col.Y1]
    Y2 = File2_data[,col.Y2]
    
    # Create the model matrix based on the seeds
    Vtheta = cbind( 1, File2_data[, where.coefs.primary ] )
    A    = crossprod(Vtheta)
    Ainv = solve(A)
    
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
  ResidSS      = crossprod(ResidSS)
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
    ResidSS      = crossprod(ResidSS)
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
      BlockRow[F2.records]    = block.order.out 
    } 
    
    # Reorder based on the sampled C 
    Y2                = Y2.in[order(New.Ordering)]
    V.Holder[,col.Y2+1] = Y2
    F2.Lambda = New.Ordering 
    
    # Store in the File 2 data space 
    File2_data[,"lambda"]   = New.Ordering
    File2_data[,"BlockRow"] = BlockRow 
 
    RIB1Base = mapply( function(x) x[which(x %in% type1seeds)], RIB[[1]])
    RIB2Base = mapply( function(x) x[which(x %in% type1seeds)], RIB[[2]])
  
  sys <-proc.time() 
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
    
    B = do.call( rbind, current.B)
    
    for( j in 1:J){
      if(class(B[,j]) != "factor"){
        B[,j] = as.factor(B[,j])
      }
    }
    
    #make changes to B, then call UpdateX() to update the data
    UpdateX(model,B)
    model$Run(0,2,2) 
    dpmpm <-model$snapshot    
    
    ###########################
    ##  Propose Block Moves  ##
    ###########################
    
    Move.F1 = which( Estar[[1]][,"Aratio"]>0)
    Move.F2 = which( Estar[[2]][,"Aratio"]>0)    
    
    if( error.file ==2 & length(Move.F2) > 0 ){
      
      move.holder =  propose.B( J, d, File1_data[,"Block"], File2_data[,"Block"], Estar, Move.F1, Move.F2, current.B, Bhat1, Bhat2, Blocking_WithError, Index_to_Block, Legal_Index,dpmpm$psi,dpmpm$z+1, File2_data[,"lambda"], File1_data[,"lambda"],gamma.samp,A.probs,BlockSeed_restriction="no",SeedsOnly,possibleEthnic,possibleSchool )
    
    if( move.holder$no.proposed > 0 ){
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
          # Step 2: Creating the k*: File 2 Record Move :  Empty k*
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
            
            SampleCkstar    = 1 
            
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
        
        c.log.A.ratio = 0     # The part of the ratio that DOES depend on the MH sampled C part 
        
        
        ###################
        ### Original To ###
        ###################
        
        if( exact[1] == 0 ){
          
          # Before the move, the to block had ABOVE the threshold number of record pairs 
          
          like.part.Orig.To = dmvnorm( y.to.orig, Design.To.Orig%*%beta.star, sigmasq*diag(n.to),log=TRUE)
          log.A.ratio       = log.A.ratio - like.part.Orig.To
          
          # The Add On
          to.counter        = non.seed.RIB[holder.to]
          log.A.ratio       = log.A.ratio + log( 2*factorial(to.counter-2) )
          
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
            
            #if( accept == 1 & exact[4] == 0 & r != reps){
             # from.order = GutmanMH( order.out.From,y.from, Design.From[order.out.From,],n.from.prop, beta.star, sigma, reps-r, holder.from,count.matches )
            #} 
            #if( accept == 1 & exact[3] == 0 & r!=reps){
            #  to.order   = GutmanMH( order.out.To,y.to, Design.To[order.out.To,],n.to.prop, beta.star, sigma, reps-r, holder.to,count.matches )
            #}
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
  
    if( length(result.out)>0){
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
      File2_data[ result.out[,2], "Block" ]                = result.out[,4]
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
      
      if( length(unlist(imps.to.remove)) > 0){
        File2_data = File2_data[-unlist(imps.to.remove),]
      }
      # Remove records from File 1, and from RIB[[1]]
      imps.from.remove = mapply( function(x,y,z) if( y !=0 ){ImpF1Which[[z]][sample(unlist(x),y,replace=F)]}, l1[unique.from], sample.from,unique.from)
      RIB[[1]][ unique.from ] = mapply( setdiff, RIB[[1]][unique.from], imps.from.remove,SIMPLIFY = FALSE )
      
      if( length(unlist(imps.to.remove)) > 0){
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
     
       if( length(Y1ImpsHave) > 0 ){
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
           ToRBindF1 = rbind(ToRBindF1, MoreRBindF1)
         }
         
          # Step 5: Tack on all imputed records 
          colnames(ToRBindF1) = colnames(File1_data) 
          end.F1              = nrow(File1_data)+1
          File1_data          = rbind(File1_data, ToRBindF1)
          nrow.F1             = nrow(File1_data)
       }
     
       if( length(Y2ImpsHave) > 0 ){
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
    

      nrow.F2             = nrow(File2_data)
      nrow.F1             = nrow(File1_data)
     
      # Step 6: update RIB 
      #RIB1Base  = mapply( function(x) x[which(x %in% type1seeds)], RIB[[1]])
      #RIB2Base  = mapply( function(x) x[which(x %in% type1seeds)], RIB[[2]])
      blocks1   = sort(unique(File1_data[-type1seeds,"Block"]))
      blocks2   = sort(unique(File2_data[-type1seeds,"Block"]))
      InF1      = split( c(1:nrow.F1)[-type1seeds] , File1_data[-type1seeds,"Block"])
      InF2      = split( c(1:nrow.F2)[-type1seeds] , File2_data[-type1seeds,"Block"])
      RIB[[1]][blocks1] = mapply(union,RIB1Base[blocks1],InF1)
      RIB[[2]][blocks2] = mapply(union,RIB2Base[blocks2],InF2)
      n.RIB[1,] = mapply(length,RIB[[1]])
      n.RIB[2,] = n.RIB[1,]
  
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
       BlockRow      = File2_data[,"BlockRow"]
     
       # For these blocks, we may need to fix the BlockRow. The records that were in the correct block already will have numbers that make sense. Everyone else will have 0s. We just need to do setdiff, and fill in the missing values. 
       HOLDER    = mapply( function(x,y) setdiff(1:x,BlockRow.F1[unlist(y)]), n.RIB[1,BlockRowCheck],RIB[[1]][BlockRowCheck])
       NEWHOLDER = mapply( function(x,y) x[which(BlockRow.F1[unlist(x)]>y)], RIB[[1]][BlockRowCheck], n.RIB[1,BlockRowCheck])
       BlockRow.F1[unlist(NEWHOLDER)] = 0
     
       # These are for the ones that are 0. 
       NEWHOLDER = mapply( function(x) x[which(BlockRow.F1[unlist(x)]==0)], RIB[[1]][BlockRowCheck])
       NEWHOLDER = unlist(NEWHOLDER)
       BlockRow.F1[NEWHOLDER] = unlist(HOLDER)
     
       # File 2 
       HOLDER    = mapply( function(x,y) setdiff(1:x,BlockRow[unlist(y)]), n.RIB[2,BlockRowCheck],RIB[[2]][BlockRowCheck])
       NEWHOLDER = mapply( function(x,y) x[which(BlockRow[unlist(x)]>y)], RIB[[2]][BlockRowCheck], n.RIB[2,BlockRowCheck])
       BlockRow[unlist(NEWHOLDER)] = 0
     
       # These are for the ones that are 0. 
       NEWHOLDER            = mapply( function(x) x[which(BlockRow[unlist(x)]==0)], RIB[[2]][BlockRowCheck])
       NEWHOLDER            = unlist(NEWHOLDER)
       BlockRow[NEWHOLDER]  = unlist(mapply( function(x) x[sample(length(x))], HOLDER))
       # We also need to give these guys lambdas, since their lambda values will also be 0. 
       F2.LambdaNew         = mapply( function(x,y)  F1.Lambda[x][BlockRow[y]], RIB[[1]], RIB[[2]])
     
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
      order.in   = BlockRow[F2.records]
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
      BlockRow[F2.records]    = block.order.out 
    } 
    
    # Reorder based on the sampled C 
    Y2                  = Y2.in[order(New.Ordering)]
    V.Holder[,col.Y2+1] = Y2
    F2.Lambda           = New.Ordering 
    
    # Store in the File 2 data space 
    File2_data[,"lambda"]   = New.Ordering
    File2_data[,"BlockRow"] = BlockRow 
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
      write( c(gamma.samp[error.file,c(4)]), file = paste(namepart,"GammaOut1.txt", sep=""), append = TRUE)
      #write( c(gamma.samp[error.file,c(6)]), file = paste(namepart,"GammaOut2.txt", sep=""), append = TRUE)
      Lambda.out <- sum(F2.Lambda[1:n2]==1:n2)
      write( Lambda.out, file = paste(namepart,"LambdaOut.txt", sep=""), append = TRUE)
      Lambda.out <-Lambda.out/nrow.F2
      write(Lambda.out, file = paste(namepart,"AllLambdaPerc.txt",sep=""),append=TRUE)
      write(nrow.F2,file=paste(namepart,"SizeOut.txt",sep=""),append=TRUE)      
      write( Accept.Count, file = paste(namepart,"AcceptOut.txt", sep=""), append = TRUE)
      Accept.Count = 0 
      #write( File2_data[1:n2,"Block"], file = paste(namepart,"AcceptOut.txt", sep=""), append = TRUE)
      # If we are interested in the completed data sets 
      if( complete.out == "yes"){
        write.table( F1.Lambda, file = paste(namepart,"F1Out.csv", row.names = F, col.names= F, sep=","), append = TRUE)
        write( Y2.in[order(F2.Lambda)], file = paste(namepart,"F2Out.txt", sep=""), append = TRUE)
      }
     
    }
    
  }
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

select.permutations<-function(length.in,star.record){
  vec.in      = 1:length.in
  vec.change  = vec.in[-star.record]
  out.choice  = matrix(vec.in,nrow=length.in,ncol=length.in,byrow=T)
  out.choice[vec.change,star.record] = vec.change
  for( i in 1:length.in){
    out.choice[vec.change[i],vec.change[i]]  = star.record  
  }
  return(out.choice)
}

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
        stop( paste("Record", r, "was proposed to be moved to an illegal block."))
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
         if( is.numeric(A.ratio)==FALSE){
          stop('We have an empty acceptance ratio during the propose B step.') 
	 } 
          #print(class(propose.move))
          
          proposal.info[[H]]  = propose.move 
        }
      }
    }
  }
  
  newlist <-list( "no.proposed" = H, "move.holder"= proposal.info )
  
  return( newlist ) 
  
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
