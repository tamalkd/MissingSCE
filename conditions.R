#########################
# This code is only intended for generating the simulation conditions for 
# running the simulations on supercomputers.
# For running the simulations on a personal computer, please refer to 'test.R'.
#########################

designs <- c("MBD", "RBD", "ABAB")
models <- c("AR1", "normal", "uniform")
ESMs <- c("MD", "NAP")
ESs <- c(0, 1, 2)
Ns <- c(40, 30, 20)
methods <- c("full", "marker", "TS", "MI")
missings <- c(0.1, 0.3, 0.5)

#######################

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

print(Result_table)
write.csv(Result_table, outfile, row.names = FALSE)
