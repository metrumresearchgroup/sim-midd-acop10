# Packages
library(mrgsolve)
library(tidyverse)

# Exercise 1:  Deterministic Simulation

# Question
# What would these patient’s exposures (AUC0-168, Cmax) be 
# if they received a single dose of 1 mg/kg OPG in the next cohort?
#OPG PK model
mod <- mread("models/opg_pk.mod")
mod

#Checkout the random effects structure
revar(mod)

#set OMEGA and SIGMA matrices to zero to remove IIV and RUV
mod <- mod %>% zero_re()
revar(mod)

#read patient weights and parameter values
cohort <- read.table("data/SAD_params.csv",header=T,sep=",")
cohort

# Input / template data set
#Here is a simple dosing data set to simulate 8 patients receiving 1 mg/kg SC
sc <- expand.ev(ID=1:8, amt=1, cmt=1)
sc

#merge in patient weights and parameter values
sc <- sc %>% left_join(cohort,by="ID")
sc

#calculate AMT in mg/kg using patient weights
sc <- sc %>% mutate(amt=1*WT)

#remove WT, so that it doesn't impact fixed effects in model
sc <- sc %>% dplyr::select(-WT)
sc

#Simulate profiles for each patient
out <- mod %>% Req(PKDV,AUC) %>% mrgsim_df(data=sc,end=168,obsonly=TRUE)

head(out)

#plot profiles for each patient
ggplot(data=out) + 
  geom_line(aes(x=time,y=PKDV,group=ID,color=factor(ID)),size=1) + 
  ylab("Concentration (ng/mL)") + xlab("Time (hours)")

#gather AUC and CMAX
summ <- out %>% group_by(ID) %>% summarize(CMAX=max(PKDV),AUC=max(AUC))
summ

summ$CMAX
summ$AUC
