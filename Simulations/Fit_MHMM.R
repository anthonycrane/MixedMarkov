# FILE:    HMM_Model1_R2jags.R
# PURPOSE: Fits a mHMM using to data generated from Data_Gen_Pipeline.R
# OUTPUT:  .RData file of a fitted MHMM model

#setwd(getwd())
library(R2jags)
library(pROC)

# Get arguments from batch script
#    Sim_File - .RData file containing simulated data
#    model    - .txt file containing JAGS model
#    outfile  - what to name the saved fitted model
args <- commandArgs(trailingOnly = TRUE)
Sim_File <- args[1]
model <- args[2]
outfile <- args[3]
lambda1 <- as.numeric(args[4])
lambda2 <- as.numeric(args[5])
file <- args[6]

# load the simulated data
load(file = Sim_File)

# Number of patients, first day of observation, and total number of observations
Npat <- length(unique(Simulated_Data$Subject))
t_1v <- rep(1,Npat)
nday <- max(Simulated_Data$Day)

# Get the last day of observation for each subject
t_endv <- c()
for(pat in unique(Simulated_Data$Subject)){
  len <- nrow(Simulated_Data[Simulated_Data$Subject==pat,])
  t_endv <- append(t_endv, len)
}


# Define a Npat x (nday + 1) matrix to hold the responses and the treatment
# arm. The first colum of Y is the treatment arm, and columns 2:(nday+1) are
# the responses
Y <- matrix(nrow=Npat,ncol=nday+1)
Zt <- data.frame(cbind(Simulated_Data$TRT1, Simulated_Data$Z,
                       Simulated_Data$Subject))

# For each subject in the dataset, extract the treatment arm and responses
# from Simulated_Data, and enter this as a row of data in Y. Y[i,1] is the
# treatment arm for subject i, Y[i,2:(nday+1)] is the responses for subject
# i across all time points
for(i in 1:Npat){
  Y[i,1] <- Zt$X1[Zt$X3==i][1]
  Y[i,2:(nday+1)] <-Zt$X2[Zt$X3==i]
}

# Save a the matrix of observations for all subjects
# and the column vector of treatment for all subjects
Z <- Y[,2:(nday+1)]
TRT1 <- Y[,1]

# Define the data and initial parameters for JAGS
# data <- c('Npat', 't_1v', 't_endv', 'Z', 'TRT1')
data <- list(Npat=Npat, t_1v=t_1v, t_endv=t_endv, Z=Z)
inits <- function(){
  list(
    lambda = c(0,0),
    beta = c(0,0),
    b = rbind(rep(0,Npat)),
    p1 = 1
  )
}

jags.m <- jags.model(file = model, data = data, inits = inits, n.chains = 4)

params<- c('beta','lambda',paste("p[1:",Npat,",1:",nday,"]",sep=''))
samps <- coda.samples(jags.m, params, n.iter=5000)
summary <- summary(samps)

beta.0 <- summary$statistics['beta[1]',1]
beta.1 <- summary$statistics['beta[2]',1]
#beta.2 <- summary$statistics['beta[3]',1]
lambda.1 <- summary$statistics['lambda[1]',1]
lambda.2 <- summary$statistics['lambda[2]',1]
pMeans <- summary(samps)$statistics[-c(1:4),1]

fitted <- matrix(nrow=Npat,ncol=nday)
k <- 0
for(i in 1:nday){
  for(j in 1:Npat){
    k <- k+1
    fitted[j,i] <- pMeans[k]
  }
}

Simulated_Data$C.posterior <- NA
for(sub in unique(Simulated_Data$Subject)){
       Simulated_Data$C.posterior[Simulated_Data$Subject==sub] <- fitted[sub,]
     }

Simulated_Data$C.pred <- ifelse(Simulated_Data$C.posterior > 0.5, 1, 0)
tab <- table(Simulated_Data$C,Simulated_Data$C.pred)
acc <- sum(diag(tab))/sum(tab)

auc <- roc(C~C.posterior, data = Simulated_Data)
AUC <- auc$auc
by.sub <- with(Simulated_Data,
     by(list(C,C.pred), INDICES = Subject, function(x){
       sum(x[[1]]==x[[2]])/length(x[[1]])
     }
)
)

# plot(C~Day, data = Simulated_Data[Simulated_Data$Subject==145,],
#      ylim=c(-.5,1.5), pch=19, type='s')
# points(C.pred~Day, data = Simulated_Data[Simulated_Data$Subject==145,],
#        col='red',pch=19, type='s')
# points(Z~Day, data = Simulated_Data[Simulated_Data$Subject==145,],pch=19)
#
# plot(C~Day, data = Simulated_Data[Simulated_Data$Subject==14,],
#      ylim=c(-.5,1.5), pch=19, type='s')
# points(C.pred~Day, data = Simulated_Data[Simulated_Data$Subject==14,],
#        col='red',pch=19, type='s')
# points(Z~Day, data = Simulated_Data[Simulated_Data$Subject==14,],pch=19)

write(c(lambda1,lambda2, paste(acc), paste(AUC),
        paste(lambda.1), paste(lambda.2), paste(beta.0), paste(beta.1)),
	 file=file, append = TRUE,ncolumns = 8)

#save(out, file=paste(outfile,'.RData',sep=''))

