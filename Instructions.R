# Code Use instructions 

# Step 1 : Set your working directory to be inside the code folder. 

# Step 2: To run a given simulation setting, choose folder corresponding to that setting. For example, to run the High Fault, High Seed setting with the Diffuse prior, choose folder "HS_HF_D". 

# Step 3: Within the folder, select the "Input" file and load it. 

# Step 4: Run the "BLASERunScript.R" 

# NOTE: The BLASERunScript runs all three methods, i.e., PG, GAZM and BLASE, creating the simulation results from the paper. 

# NOTE 2: The script exports output to .RData and .csv files. Make sure to set the working directory (setwd) to the folder where you want such output to be stored. 

# Step 5: To create the output from the paper, use the file BLASE_Eval_Paper.R. This will produce the plots/results presented in the paper and appendix. 