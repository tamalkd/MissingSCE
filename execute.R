#########################
# This code is only intended for executing the simulation on supercomputers.
# For running the simulations on a personal computer, please refer to 'test.R'.
#########################

args <- commandArgs(TRUE)

print(args)

design <- as.character(args[1])
model <- as.character(args[2])
ESM <- as.character(args[3])
ES <- as.integer(args[4])
N <- as.integer(args[5])
method <- as.character(args[6])
missing <- as.numeric(args[7])

options(warn=-1)

#####################
       
AR=0.6
limit_phase=3
reps_MBD=4
num_MI=10
number=1
number_MC=1000
alfa=0.05
replications=100000

#######################

strt <- Sys.time()

# library(foreach)
# library(parallel)
# library(doParallel)
# cores <- detectCores()
# print(cores)
# cl <- makeCluster(cores - 1)
# registerDoParallel(cl) 

source("RT.R")
source("power_calc.R")

set.seed(1000)

#######################

output <- numeric()
for(it in 1:replications)
{
  output[it] <- Calculate_conditional_power_random(
    design=design, 
    model=model, 
    number=number, 
    ESM=ESM, 
    limit_phase=limit_phase,
    reps_MBD=reps_MBD,
    num_MI=num_MI,
    direction="+", 
    alfa=alfa, 
    N=N, 
    ES=ES, 
    AR=AR,
    number_MC=number_MC,
    method=method,
    missing=missing
  )
}

power <- mean(output)
end <- Sys.time()

Result_table <- data.frame(
  design = design, 
  model = model, 
  ESM = ESM, 
  ES = ES, 
  N = N, 
  method = method,
  missing = missing,
  Reps = replications,
  nMC = number_MC,
  nCP = number,
  timer = difftime(end, strt, units = "secs"),
  power = power,
  stringsAsFactors = FALSE
)

print(Result_table)
write.table(Result_table, "Results", append = TRUE, row.names = FALSE, col.names = FALSE)
