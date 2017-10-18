# User Inputs
# Simulation Settings: Seed: High, Fault: High, Prior: SC

# How many data points do we want? 
N = 5000

# Type 1 Seeds? 
perc.type1 = 60

# Error? 
error.perc = 40

# Prior Choice" 
prior.choice = "CA"

# Simulation Name
sim.name = "HS_HF_CA"

# How many iterations
its = 10000
# How many draws are discarded? 
burnin = 1
# How many draws in between each saved data set? 
gap = 25
# How many maximum "Swap steps" do we allow in the Gutman MCMC? 
reps = 20
# How do we thin our chain? 
thinning = 2
# How often do we print out the iteration numbers?
printgap = 1000
# Do you want MUCH simulation information to output as you run? 
Comment = "FALSE"

DPburnin   = 1000
DPits      = 10000
DPthinning = 2

# Do we restrict possible moves to the observed blocks? "yes" or "no", default "yes"
ObsBlock_restriction = "yes"

# Do we allow moves to blocks which only have seeds? "yes" or "no", default "yes"
BlockSeed_restriction = "no"
#Do we store the matched data sets? 
complete.out = "no" 

MCMCSeed = NULL 
## %%%%%%%%%%%%%%%%%%%%% # 
## Regression Constants ##
## %%%%%%%%%%%%%%%%%%%%% # 

# In File2_data, what are the catergorical variables in the primary regression? 
coefs.primary      = c("prog")
# In File2_data, what are the variables in the primary regression? 
coefs.primary.all  = c("math","prog")
# In File2_data, where are the coefficients for the secondary regression?
coefs.secondary    = c("female", "prog","ses","honors")

# Note to Self: Consider the possibility that the Files that are input may have things that are NOT used for linking. We need to have a way to create (1) Y1, (2) Y2, (3) the other common variables for the regression, (4) the variables not involved in the linking process which we don't need until the very end when we release the linked data copies. 

# In which column in File1_data will we find Y1?
col.Y1 = 1
# In which column in File2_data will we find Y2?
col.Y2 = 1

## %%%%%%%%%%%%%%%%%%%%% # 
## Matching Constants ##
## %%%%%%%%%%%%%%%%%%%%% # 

# Which file has error? Options = 1, 2, "both"
error.file  = 2

# Which field(s) are allowed to be in error (default NULL). If NULL, run Gutman only. 
error.fields = c("prog")

# Which variables are involved in the DP? 
# This should probably be names, and the package should locate them and convert them into numbers. 
# The ordering HAS to match the ordering of the variables in File 1 and File 2 
blockingVars = c("female","cid","schtyp","prog","ses","honors")
J            = length(blockingVars)

secondary.option = "dependent"
# Does the regression for y2 depend on B? 
secondary.option2 = "dependent"

## %%%%%%%%%%%%%%%%%%%%% # 
##     Prior Values     ##
## %%%%%%%%%%%%%%%%%%%%% # 

# a.beta.prior
a.beta.prior <- matrix(0,nrow = 2, ncol = J)
colnames(a.beta.prior) = blockingVars
a.beta.prior[,error.fields] = 90000 
# b.beta.prior 
b.beta.prior <- matrix(0,nrow = 2, ncol = J)
colnames(b.beta.prior) = blockingVars
b.beta.prior[,error.fields] = 10000

# Define the number of latent components to be used in the DP
Hstar = 100 

# Do we use all of the available blockign variables in the DP? If not, specify what we shift by. 
DPOffset = 0

threshold = 5  

