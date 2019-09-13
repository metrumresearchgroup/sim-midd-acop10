library(tidyverse)

# data <- readRDS('opgpost.RDS')
# 
# names(data)[1:13] <- paste0("THETA", 1:13)
# 
# write_csv(data, "post.csv")

revar(mod)

library(mrgsolve)
 
mod <- mread("opg.mod")

dose <- ev(amt = 100)

out <- mrgsim(mod, dose)

plot(out)



