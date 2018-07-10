#####################

designs <- c("MBD", "RBD", "ABAB")
models <- c("AR1", "normal", "uniform")
ESMs <- c("MD", "NAP")
ESs <- c(0, 1, 2)
Ns <- c(40, 30, 20)
methods <- c("full", "marker", "TS", "MI")
missings <- c(0.1, 0.3, 0.5)

#####################

AR=0.6
limit_phase=3
reps_MBD=4
num_MI=10
replications=1000
number=1
number_MC=1000
alfa=0.05

#######################

library(foreach)
# library(parallel)
# library(doParallel)
# cores <- detectCores()
# print(cores)

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

start <- Sys.time()

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
              #cl <- makeCluster(cores - 2)
              #registerDoParallel(cl) 
              
              set.seed(1000)
              
              output <- foreach(it=1:replications, .inorder=FALSE, .combine='c') %do%
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
              start <- Sys.time()
              
            } else
            {
              for(o in 1:length(missings))
              {
                #cl <- makeCluster(cores-2)
                #registerDoParallel(cl) 
                
                set.seed(1000)
                
                output <- foreach(it=1:replications, .inorder=FALSE, .combine='c') %do%
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
                start <- Sys.time()
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
