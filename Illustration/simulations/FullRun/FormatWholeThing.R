# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 
# Format the full data to be used in the small sample collection
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # 

# Load in the Full Data 
load("C:/Users/nmd16/Desktop/Summer2015/PreErrorIntro.RData")

# Set a Seed 
set.seed(9292015)

# %%%%%%%%%%%%%%%%%%%%%%%%%% # 
#  Create the school code  
# %%%%%%%%%%%%%%%%%%%%%%%%%% # 

Code2012 = paste(WholeThing$lea,WholeThing$schlcode,sep="")
Code2011 = paste(WholeThing$lea2011,WholeThing$schlcode2011,sep="")

# Now, we need to go ahead and create the rest based on these new values. 

File2011$schlcode = Code2011
File2012$schlcode = Code2012

# %%%%%%%%%%%%%%%%%%%%%% # 
# Ethnicity P : Remove 
# %%%%%%%%%%%%%%%%%%%%%% # 

# Remove the individuals who had ethnicity P. 
PEth = which(WholeThing[,"ethnic"]=="P")

File2012 = File2012[-PEth,]
File2011 = File2011[-PEth,]
rownames(File2012) = 1:nrow(File2012)
rownames(File2011) = 1:nrow(File2011)

# Remove the factor from consideration
File2012[,"ethnic"] = factor(File2012[,"ethnic"])
File2011[,"ethnic"] = factor(File2011[,"ethnic"])

# Ditch what we don't need
rm(DisagreePool,Month_Diff,RIB,Seeds,Singletons,bvars,i,inError,type1Seeds,Type1Only,Code2011,Code2012)

# %%%%%%%%%%%%%%%%%%%%%%%%%% # 
#  Make all other fields agree  
# %%%%%%%%%%%%%%%%%%%%%%%%%% # 

# Set all other fields correct.  
File2011Orig              = File2011
File2011[,"birthday"]     = File2012[,"birthday"]
File2011[,"birthmonth"]   = File2012[,"birthmonth"]
File2011[,"birthyear"]    = File2012[,"birthyear"]
File2011[,"ethnic"]       = File2012[,"ethnic"]
File2011[,"sex"]          = File2012[,"sex"]
File2011[,"schlcode"]     = File2012[,"schlcode"]
File2011[,"YearPost1996"] = File2012[,"YearPost1996"]

# The ratio of seeds to non seeds is 60% seeds, so let's keep that. 
File2011New = File2011
File2012New = File2012

File2012New[,"ethnic"] = relevel(File2012New[,"ethnic"],ref="W")
File2011New[,"ethnic"] = relevel(File2011New[,"ethnic"],ref="W")

# %%%%%%%%%%%%%%%%%%%%%%%%%% # 
#  Remove High Leverage Points  
# %%%%%%%%%%%%%%%%%%%%%%%%%% #

# Fit the Desired Model to this data 
Model1 = lm( File2012New[,1] ~ File2011New[,"mathscal"] + YearPost1996 + ethnic + sex, data = File2012New)

# Find all the points with extreme cooks distance. 
# See Document Model Diagnostics for Justification for the removal of these points 
CD = cooks.distance( Model1)
CDWeird = which(CD > 0.00080)

# Which of these have ethnicity I? 
ethnic.I = which(File2012New[,"ethnic"]=="I")
CDWeirdI = intersect( ethnic.I, CDWeird)

ModelNoCD = lm(File2012New[-CDWeirdI,1] ~ File2011New[-CDWeirdI,"mathscal"] + YearPost1996 + ethnic + sex, data = File2012New[-CDWeirdI,])

# Remove all points with strange ethnicity
CD = cooks.distance( ModelNoCD)
CDWeird = which(CD > 0.00080)

toremove = c(CDWeirdI,CDWeird)
ModelRemoved = lm(File2012New[-toremove,1] ~ File2011New[-toremove,"mathscal"] + YearPost1996 + ethnic + sex, data = File2012New[-toremove,])

truthModel = ModelRemoved
truthCoefs = truthModel$coefficients

# %%%%%%%%%%%%%%%%%%%%%%%%%% # 
# Convert to Factors   
# %%%%%%%%%%%%%%%%%%%%%%%%%% #

colnames(File2011New)
File2011New[,"birthday"] = factor(File2011New[,"birthday"])
File2012New[,"birthday"] = factor(File2012New[,"birthday"])

File2011New[,"birthmonth"] = factor(File2011New[,"birthmonth"])
File2012New[,"birthmonth"] = factor(File2012New[,"birthmonth"])

File2011New[,"birthyear"] = factor(File2011New[,"birthyear"])
File2012New[,"birthyear"] = factor(File2012New[,"birthyear"])

UniqueList = sort(unique(c(File2012New[,"schlcode"],File2011New[,"schlcode"])))
SchoolOut1 = apply( as.matrix(File2012[,"schlcode"]),1, function(x) which(UniqueList == x) )
File2012New[,"schlcode"] = SchoolOut1
SchoolOut1 = apply( as.matrix(File2011[,"schlcode"]),1, function(x) which(UniqueList == x) )
File2011New[,"schlcode"] = SchoolOut1

File2011New[,"schlcode"] = factor(File2011New[,"schlcode"])
File2012New[,"schlcode"] = factor(File2012New[,"schlcode"])

# %%%%%%%%%%%%%%%%%%%%%%%%%% # 
# Save and Delete 
# %%%%%%%%%%%%%%%%%%%%%%%%%% #

File2011 <- File2011New
File2012 <- File2012New

rm( File2011New, File2012New, File2011Orig, WholeThing, CD, CDWeird, CDWeirdI, Model1, ModelNoCD, ModelRemoved, ModelSeeds, N, PEth, SchoolOut1, UniqueList, ethnic.I, toremove)

save.image("WholePostFormat.RData")