#!/bin/sh                                                                       
#$ -l mem=4G,time=10:: -S /bin/bash -N LAMBDAS -j y -cwd                          

R=/nfs/apps/R/3.1.1/bin/R
Rscript=/nfs/apps/R/3.1.1/bin/Rscript
R_LIBS_USER=/ifs/home/msph/software/library:~/local/R/titan:/ifs/home/msph/biostat/acs2244/local/R/hpc:$R_LIBS_USER

lambda1=`awk 'NR=='$SGE_TASK_ID'{print $1}' lambda.txt`
lambda2=`awk 'NR=='$SGE_TASK_ID'{print $2}' lambda.txt`

${Rscript} --vanilla  Simulate_Data.R 156 42 $lambda1 $lambda2 Simulated_Data_$SGE_TASK_ID

${Rscript} --vanilla  Fit_MHMM.R Simulated_Data_$SGE_TASK_ID.RData JAGS_Model.txt Fitted_Data_$SGE_TASK_ID $lambda1 $lambda2 ROC_$SGE_TASK_ID.txt

mv ROC_$SGE_TASK_ID.txt output/

rm Simulated_Data_$SGE_TASK_ID.RData