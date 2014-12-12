# FILE:     Data_Gen_pipeline.R
# PURPOSE:  Generated data simulated from a HMM with specified parameter vals
# OUTPUT:   .Rdata file of simulated data from a HMM with emission parameters
#           lambda.0, lambda.1
setwd(getwd())
source('Drugs_Functions.R')
set.seed(42)

# Get parameter values from batch script
#  NPat     - number of patients
#  nday     - number of days of study
#  lambda.0 - Log odds of observing drug use while in state 0
#  lambda.1 - Log OR of observing drug use in state 1 compared to
#             state 0
#  outfile  - Name of file in which simulated data is saved
args <- commandArgs(trailingOnly = TRUE)
Npat <- as.numeric(args[1])
nday <- as.numeric(args[2])
lambda.0 <- as.numeric(args[3])
lambda.1 <- as.numeric(args[4])
outfile <- args[5]

# Not run:
# For running this file in non-batch mode:
# Npat <-     156
# nday <-     84
# lambda.0 <- -1.5
# lambda.1 <- 4

# Initialize parameter values
#   beta.0   - Intercept for being in state 1 compared to state 0
#              (log odds of being in state 0, given in state 0 at
#              previous time point and in control group)
#   beta.1   - Log OR of being in state 1 given in state 1 at previous
#              time point
#   beta.2   - Treatment effect of being in state 1 compared to state
#              0 (log OR)

beta.0 <- -3.85
beta.1 <- 7.334
beta.2 <- 0 # For these simulations we do not consider any treatment effect
b.0 <- rnorm(Npat, 0, 2)

# Generate subject ID's, Day of trial, and Treatment Groups
# Save as a nday x 3 matrix, X.
Subject <- rep(x = 1:Npat, each = nday)
Day <- rep(x = 1:nday, times = Npat)
TRT1 <- c(rep(0,each=Npat*nday/2), rep(1,each=Npat*nday/2))

X <- as.matrix(data.frame(
  Subject = Subject,
  Day  = Day,
  TRT1 = TRT1
))

# Initialize an empty matrix with ncol(X) + 2 columns in which to store
# simulated data
D <- matrix(data = NA, nrow = 0, ncol = ncol(X)+2)

# For each subject, select their rows of data in X and simulate hidden states
# and observed use across all days. Append the simulated data to D.
for (i in unique(X[,1])) {
  X.i <- X[which(X[,1]  == i),]
  D.i <- Data.gen(X.i = X.i, nday = nday, Beta = beta.2, beta.0 = beta.0,
                  beta.1 = beta.1, b.0 = b.0[i], lambda.0 = lambda.0,
                  lambda.1 = lambda.1)
  D <- rbind(D, D.i)
}

# Convert the simulated data matrix to a data frame
Simulated_Data <- data.frame(D)

# Output the simulated data to 'outfile'.RData
save(Simulated_Data, file = paste(outfile,'.RData',sep=''))
