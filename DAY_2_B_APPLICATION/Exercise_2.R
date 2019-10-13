# Packages
library(mrgsolve)
library(tidyverse)

# Exercise 2:  Stochastic Simulation

# Question
# What dose is required to maintain steady-state trough concentrations above 90 ng/mL in 80% of patients, given weekly dosing?

#OPG PK model
mod <- mread("models/opg_pk.mod")
mod

#Checkout the random effects structure
revar(mod)

#Set SIGMA to zero to remove observation error
mod <- mod %>% zero_re(sigma)
revar(mod)

#load database of patient weights and sample 1000 values to use in dose amount calculation
demo <- read.table("data/demographics.csv",header=T,sep=",") %>% sample_n(1000)

# Input / template data set
#dosing for 4 weeks, so addl=3
# 1000 patients with doses from 1-12 mg/kg
sc <- expand.ev(ID=1, amt=0, cmt=1, addl=3, ii=168, DOSE=1:12, WT=demo$WT)
#set dose amount using patient weight
sc <- sc %>% mutate(amt=DOSE*WT)
sc

#We will get the observation design for the simulation through a `tgrid` object
#collect concentrations over 4 week study
des <- tgrid(end=-1,add=seq(0,672,1))
des

#Simulate profiles for all patients
out <- mod %>% carry_out(DOSE) %>% mrgsim(data=sc,tgrid=des,obsonly=TRUE)
head(out)

#summarize and plot
summ <- out %>% group_by(DOSE,time) %>% summarize(Q20=quantile(PKDV,prob=c(0.20)),  # 20th percentile to visualize minimum for 80% of patients
                                                  Q50=quantile(PKDV,prob=c(0.50)),
                                                  Q90=quantile(PKDV,prob=c(0.90)))  # 90th percentile to visualize maximum for 90% of patients
ggplot(data=summ) + geom_ribbon(aes(x=time,ymin=Q20,ymax=Q90),fill='cyan',alpha=0.5) + 
                    geom_line(aes(x=time,y=Q50)) + facet_wrap(~DOSE) +
                    geom_hline(yintercept=90,linetype="dashed") +
                    geom_hline(yintercept=5000,linetype="dashed") 

#function to calculate fraction of patients above a threshold

#summarize and plot
#calculate what fraction of trough concentrations are above 90 ng/mL at 672 hours
fraction_above <- function(x,boundary) { length(x[x>boundary])/length(x) }
above <- out %>% filter(time==672) %>% group_by(DOSE) %>% summarize(FRAC=fraction_above(PKDV,90))
above
ggplot(data=above) + geom_line(aes(x=DOSE,y=FRAC),color='red',size=1.5) + geom_hline(yintercept=0.8,linetype="dashed") +
                     xlab("Dose (mg/kg") + ylab("Fraction Above 90 ng/mL")


