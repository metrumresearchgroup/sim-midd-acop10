# Packages
library(mrgsolve)
library(tidyverse)

# Question
# What is the maximum single dose of OPG that can be administered to these patients without Cmax exceeding 5,000 ng/mL?

#OPG PK model
#Need to add AUC compartment
mod <- mread("opg_pk_AUC.mod")
mod

#Checkout the random effects structure
revar(mod)

#set OMEGA and SIGMA matrices to zero to remove IIV and RUV
mod <- mod %>% zero_re()
revar(mod)

#read patient weights and parameter values
cohort <- read.table("SAD_params.csv",header=T,sep=",")
cohort

#simple function to simulate a cohort at a specified dose and calculate cmax for each patient
simCohort <- function(dose)  {
  sc <- expand.ev(ID=1:8, amt=1, cmt=1)
  sc <- sc %>% left_join(cohort)
  sc <- sc %>% mutate(amt=dose*WT) %>% dplyr::select(-WT)
  des <- tgrid(end=-1,add=seq(0,168,1))
  res <- mod %>% Req(PKDV,AUC) %>% mrgsim(data=sc,tgrid=des,obsonly=TRUE) %>% as.data.frame() 
  out <- res %>% group_by(ID) %>% summarize(CMAX=max(PKDV)) %>% mutate(DOSE=dose)
  return(out)
}

#simulate doses from 1-20 mg/kg
cmax <- lapply(1:20,simCohort) %>% bind_rows()

#calculate maximum cmax for each dose
results <- cmax %>% group_by(DOSE) %>% summarize(MAX=max(CMAX))
results

ggplot(data=results) + geom_line(aes(x=DOSE,y=MAX),color='red',size=1.5) + ylab("Maximcum Concenration (ng/mL)") + 
                       xlab("Dose (mg/kg)") + geom_hline(yintercept=5000,linetype='dashed')



