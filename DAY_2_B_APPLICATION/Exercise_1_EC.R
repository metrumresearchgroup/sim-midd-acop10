# Packages
library(mrgsolve)
library(tidyverse)

# Exercise 1:  Deterministic Simulation

# Question
# What is the maximum single dose of OPG that can be administered to these 
# patients without Cmax exceeding 5,000 ng/mL?

#OPG PK model
#Need to add AUC compartment
mod <- mread("models/opg_pk_AUC.mod")
mod

#Checkout the random effects structure
revar(mod)

#set OMEGA and SIGMA matrices to zero to remove IIV and RUV
mod <- mod %>% zero_re()
revar(mod)

#read patient weights and parameter values
cohort <- read.table("data/SAD_params.csv",header=T,sep=",")
cohort

#simple function to simulate a cohort at a specified dose and calculate cmax 
#for each patient
simCohort <- function(dose)  {
  sc <- expand.ev(ID=1:8, amt=1, cmt=1)
  sc <- sc %>% left_join(cohort,by="ID")
  sc <- sc %>% mutate(amt=dose*WT) %>% dplyr::select(-WT)
  res <- mod %>% Req(PKDV,AUC) %>% mrgsim_df(data=sc,end=168,obsonly=TRUE)
  out <- res %>% group_by(ID) %>% summarize(CMAX=max(PKDV)) %>% mutate(DOSE=dose)
  return(out)
}

#simulate doses from 1-20 mg/kg
cmax <- map_df(1:20,simCohort)

#calculate maximum cmax in these patients for each dose
results <- cmax %>% group_by(DOSE) %>% summarize(MAX=max(CMAX))
results

#plot maximum Cmax by dose

qplot(DOSE,MAX,data = results, geom = c("line", "point"), colour=I("red3")) + 
  ylab("Maximum Concenration (ng/mL)") + 
  xlab("Dose (mg/kg)") +
  geom_hline(yintercept=5000,linetype='dashed')

