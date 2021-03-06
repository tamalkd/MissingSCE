#########################
# This code is only intended for generating the simulation conditions for 
# running the simulations on supercomputers.
# For running the simulations on a personal computer, please refer to 'test.R'.
#########################

### All simulation conditions

designs <- c("MBD", "RBD", "ABAB")          # SCE design types
models <- c("AR1", "normal", "uniform")     # Data models
ESMs <- c("MD", "NAP")                      # Test statistics
ESs <- c(0, 1, 2)                           # Effect sizes
Ns <- c(40, 30, 20)                         # Number of measurements
methods <- c("full", "marker", "TS", "MI")  # Missing data handling methods
missings <- c(0.1, 0.3, 0.5)                # Proportion of missing data 

### Generate all possible combinations of simulation conditions

outfile <- "start.csv"

Result_table <- data.frame(
  design = character(),
  model = character(),
  ESM = character(),
  ES = numeric(),
  N = integer(),
  method = character(),
  missing = numeric(),
  stringsAsFactors = FALSE
)

print(Sys.time())

for (i in 1:length(designs))
{
  for (m in 1:length(Ns)) 
  {
    for (k in 1:length(ESMs)) 
    {
      for (j in 1:length(models)) 
      {
        for (l in 1:length(ESs)) 
        {
          for(n in 1:length(methods))
          {
            if(methods[n]=="full")
            {
              outlist <- data.frame(
                design = designs[i], 
                model = models[j], 
                ESM = ESMs[k], 
                ES = ESs[l], 
                N = Ns[m], 
                method = methods[n],
                missing = 0
              )
              Result_table <- rbind(Result_table, outlist)
            } else
            {
              for(o in 1:length(missings))
              {
                outlist <- data.frame(
                  design = designs[i], 
                  model = models[j], 
                  ESM = ESMs[k], 
                  ES = ESs[l], 
                  N = Ns[m], 
                  method = methods[n],
                  missing = missings[o]
                )
                Result_table <- rbind(Result_table, outlist)
              }
            }
          }
        }
      }
    }
  }
}

### Publish table containing all possible combinations

print(Result_table)
write.csv(Result_table, outfile, row.names = FALSE)
