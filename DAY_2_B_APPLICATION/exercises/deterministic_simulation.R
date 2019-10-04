# Packages
library(mrgsolve)
library(tidyverse)

# Question
# What would these patient’s exposures (AUC, Cmax) be if they receive a single dose of 1 mg/kg OPG in the next cohort?

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

# Input / template data set
#Here is a simple dosing data set to simulate 8 patients receiving 1 mg/kg SC
sc <- expand.ev(ID=1:8, amt=1, cmt=1)
sc

#merge in patient weights and parameter values
sc <- sc %>% left_join(cohort)
sc

#calculate AMT in mg/kg using patient weights
sc <- sc %>% mutate(amt=1*WT)
#remove WT, so that it doesn't impact fixed effects in model
sc <- sc %>% dplyr::select(-WT)
sc

#We will get the observation design for the simulation through a `tgrid` object
des <- tgrid(end=-1,add=seq(0,168,1))
des

#Simulate profiles for each patient
out <- mod %>% Req(PKDV,AUC) %>% mrgsim(data=sc,tgrid=des,obsonly=TRUE) %>% as.data.frame()
head(out)

#plot profiles for each patient
ggplot(data=out) + geom_line(aes(x=time,y=PKDV,group=ID,color=factor(ID)),size=1.5) + ylab("Concentration (ng/mL)") + xlab("Time (hours)")

#gather AUC and CMAX
summ <- out %>% group_by(ID) %>% summarize(CMAX=max(PKDV),AUC=max(AUC))
summ

summ$CMAX
summ$AUC
