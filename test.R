#####################
# PLEASE READ CAREFULLY BEFORE EXECUTING!
#
# This code would normally be run on supercomputers. The entirety of the computation is expected to take 
# 20000 hours on a single core processor running at 2.5 GHz. 
#
# Please only run one condition at a time if running on a personal computer. This means only one choice 
# among 'designs' (single-case design), 'models' (simulation data model), 'ESMs' (effect size 
# measure/test statistic), 'ESs' (effect size applied), 'Ns' (number of measurements), 'methods' (missing 
# data handling method), and 'missings' (missing data proportion) should be kept, all others should to 
# be deleted.
#
# We strongly recommend setting parameter 'replications' (number of simulated datasets) to 1000 or lower.
# This would lead to lower accuracy for power/type I error estimate, however otherwise the computation 
# time would be too long.
#
# PLease ensure the following files are in the same folder (or in the working directory):
# 'RT.R' 'power_calc.R' 'NAP.cpp'
#
# Please ensure the following R packages are installed:
# 'foreach' 'parallel' 'doParallel' 'imputeTS' 'data.table' 'mice' 'Rcpp' 
#####################

# designs <- c("MBD", "RBD", "ABAB")
# models <- c("AR1", "normal", "uniform")
# ESMs <- c("MD", "NAP")
# ESs <- c(0, 1, 2)
# Ns <- c(40, 30, 20)
# methods <- c("full", "marker", "TS", "MI")
# missings <- c(0.1, 0.3, 0.5)

designs <- c("RBD", "ABAB")
models <- c("normal")
ESMs <- c("NAP")
ESs <- c(1)
Ns <- c(30)
methods <- c("marker")
missings <- c(0.3)

#####################

AR = 0.6
limit_phase = 3
reps_MBD = 4
num_MI = 10
number = 1
number_MC = 1000
alfa = 0.05
replications = 100000

######################

library(foreach)
library(parallel)
library(doParallel)

cores <- detectCores()
print(cores)

source("RT.R")
source("power_calc.R")

Result_table <- data.frame(
  design = character(),
  model = character(),
  ESM = character(),
  ES = numeric(),
  N = integer(),
  method = character(),
  missing = numeric(),
  Reps = integer(),
  nMC = integer(),
  nCP = integer(),
  timer = numeric(),
  power = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(designs))
{
  for (j in 1:length(models)) 
  {
    for (k in 1:length(ESMs)) 
    {
      for (l in 1:length(ESs)) 
      {
        for (m in 1:length(Ns)) 
        {
          for(n in 1:length(methods))
          {
            if(methods[n]=="full")
            {
              cl <- makeCluster(max(cores - 2, 1))
              registerDoParallel(cl) 
              
              start <- Sys.time()
              set.seed(1000)
              
              output <- foreach(it=1:replications, .inorder=FALSE, .combine='c') %dopar%
              {
                result <- Calculate_conditional_power_random(
                  design=designs[i], 
                  model=models[j], 
                  number=number, 
                  ESM=ESMs[k], 
                  limit_phase=limit_phase,
                  reps_MBD=reps_MBD,
                  num_MI=num_MI,
                  direction="+", 
                  alfa=alfa, 
                  N=Ns[m], 
                  ES=ESs[l], 
                  AR=AR,
                  number_MC=number_MC,
                  method=methods[n]
                )
              }
              
              power <- mean(output)
              end <- Sys.time()
              
              outlist <- data.frame(
                design = designs[i], 
                model = models[j], 
                ESM = ESMs[k], 
                ES = ESs[l], 
                N = Ns[m], 
                method = methods[n],
                missing = 0,
                Reps = replications,
                nMC = number_MC,
                nCP = number,
                timer = as.numeric(difftime(end, start, units = "secs")),
                power = power
              )
              Result_table <- rbind(Result_table, outlist)
              
              print(power)
            
            } else
            {
              for(o in 1:length(missings))
              {
                cl <- makeCluster(max(cores - 2, 1))
                registerDoParallel(cl)
                
                start <- Sys.time()
                set.seed(1000)
                
                output <- foreach(it=1:replications, .inorder=FALSE, .combine='c') %dopar%
                {
                  result <- Calculate_conditional_power_random(
                    design=designs[i], 
                    model=models[j], 
                    number=number, 
                    ESM=ESMs[k], 
                    limit_phase=limit_phase, 
                    reps_MBD=reps_MBD,
                    num_MI=num_MI,
                    direction="+", 
                    alfa=alfa, 
                    N=Ns[m], 
                    ES=ESs[l], 
                    AR=AR,
                    number_MC=number_MC,
                    method=methods[n],
                    missing=missings[o]
                  )
                }
                
                power <- mean(output)
                end <- Sys.time()
                
                outlist <- data.frame(
                  design = designs[i], 
                  model = models[j], 
                  ESM = ESMs[k], 
                  ES = ESs[l], 
                  N = Ns[m], 
                  method = methods[n],
                  missing = missings[o],
                  Reps = replications,
                  nMC = number_MC,
                  nCP = number,
                  timer = as.numeric(difftime(end, start, units = "secs")),
                  power = power
                )
                Result_table <- rbind(Result_table, outlist)
                
                print(power)
                
              }
            }
          }
        }
      }
    }
  }
}

print(Result_table)
write.table(Result_table, "Results", append = TRUE, row.names = FALSE, col.names = FALSE)
