# Packages
library(mrgsolve)
library(tidyverse)
library(parallel)

# Question
# What dose would provide a 90% probability of achieving a 40% reduction in NTX with Q2W dosing? 

#OPG PK/PD model
mod <- mread("opg_PCFB.mod")
mod

#Checkout the random effects structure
revar(mod)

#set OMEGA and SIGMA matrices to zero to remove IIV and RUV
mod <- zero_re(mod)
revar(mod)

# Load the simulated posterior
post <- read.table("post_POC.csv",header=T,sep=",") %>% sample_n(1000)

#load database of patient weights and sample 1000 values to use in dose amount calculation
demo <- read.table("demographics.csv",header=T,sep=",") %>% sample_n(1000)

#Q2W doses
# Input / template data set
sc <- expand.ev(ID=1, amt=0, cmt=1, dose=1:10, ii=336, addl=5)
# tgrid object
des <- tgrid(end=-1,add=seq(0,2016,12))

#Simulation function
sim <- function(i,data,des) {
  #sample WTs and set dosing amount
  data[,"amt"] <- demo$WT[i]*data[,"dose"]
  #set parameter value
  mod <- mod %>% param(slice(post,i)) 
  #simulate and add replicate number to data frame
  mod %>% carry_out(evid,dose) %>% mrgsim(data=data,tgrid=des,obsonly=TRUE) %>% mutate(irep=i,dose=dose) %>% filter(evid==0)
}

options(mc.cores=8)
mcRNG()
set.seed(9999)
out <- mclapply(1:1000, sim, data=sc, des=des) %>% bind_rows
head(out)

summ <- out %>% group_by(dose,time) %>% summarize(Q10=quantile(PCFB,prob=c(0.10)),
                                                  Q50=quantile(PCFB,prob=c(0.50)),
                                                  Q90=quantile(PCFB,prob=c(0.90)))

#Plot the distribution of percent change from baseline
ggplot(data=summ) + geom_ribbon(aes(x=time,ymin=Q10,ymax=Q90),fill='cyan',alpha=0.5) + geom_line(aes(x=time,y=Q50)) + 
  facet_wrap(~dose) + geom_hline(yintercept=-50,linetype="dashed") 

fraction_below <- function(x,boundary) { length(x[x<boundary])/length(x) }
below <- out %>% filter(time==2016) %>% group_by(dose) %>% summarize(FRAC=fraction_below(PCFB,-40))
below
ggplot(data=below) + geom_line(aes(x=dose,y=FRAC),color='red',size=1.5) + geom_hline(yintercept=0.9,linetype="dashed") +
  xlab("Dose (mg/kg)") + ylab("Probability of Achieving 40% NTX Reduction")



#Q3W doses

# Input / template data set
sc <- expand.ev(ID=1, amt=0, cmt=1, dose=seq(10,100,10), ii=504, addl=3)
# tgrid object
des <- tgrid(end=-1,add=seq(0,2016,12))

out <- mclapply(1:1000, sim, data=sc, des=des) %>% bind_rows
head(out)

summ <- out %>% group_by(dose,time) %>% summarize(Q10=quantile(PCFB,prob=c(0.10)),
                                                  Q50=quantile(PCFB,prob=c(0.50)),
                                                  Q90=quantile(PCFB,prob=c(0.90)))

#Plot the distribution of percent change from baseline
ggplot(data=summ) + geom_ribbon(aes(x=time,ymin=Q10,ymax=Q90),fill='cyan',alpha=0.5) + geom_line(aes(x=time,y=Q50)) + 
  facet_wrap(~dose) + geom_hline(yintercept=-50,linetype="dashed") 

below <- out %>% filter(time==2016) %>% group_by(dose) %>% summarize(FRAC=fraction_below(PCFB,-40))
ggplot(data=below) + geom_line(aes(x=dose,y=FRAC),color='red',size=1.5) + geom_hline(yintercept=0.9,linetype="dashed") +
  xlab("Dose (mg/kg)") + ylab("Probability of Achieving 40% NTX Reduction")
