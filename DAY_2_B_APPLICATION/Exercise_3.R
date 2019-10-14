# Packages
library(mrgsolve)
library(tidyverse)
library(furrr)

options(future.fork.enable=TRUE)
plan(multiprocess,workers=4L)
opt <- future_options(seed = TRUE)

# Exercise 3:  Uncertainty Simulation

# Question
# What is the probability that a 5 mg/kg dose would achieve the target of 
# 40% mean reduction in NTX, if administered every 2 weeks for 3 months?

#OPG PK/PD model
mod <- mread("models/opg_pkpd.mod")
mod

#Checkout the random effects structure
revar(mod)

#set OMEGA and SIGMA matrices to zero to remove IIV and RUV
mod <- zero_re(mod)
revar(mod)

# Load the simulated posterior
post <- read.table("data/post_POC.csv",header=T,sep=",") %>% sample_n(1000)

#load database of patient weights and calculate mean to use in dose amount calculation
wt <- read.table("data/demographics.csv",header=T,sep=",") %>% summarize(WT=mean(WT))

# Input / template data set
#Here is a simple dosing data set to simulate 1 patient
#3 months of Q2W dosing is 6 total doses, so addl=5
sc <- expand.ev(ID=1, amt=0, cmt=1, dose=5, ii=336, addl=5)

#We will get the observation design for the simulation through a `tgrid` object
#Interested in 2 weeks following final dose, 2016 hours
des <- tgrid(end  = 2016)

des

# Replicate simulation with uncertainty
#When we do replicate simulation, it almost always pays off
#to construct a function that will carry out one replicate.

#Arguments to the function are
#- `i` the current simulation replicate
#- `data` the dosing data set
#- `des` the observation design

sim <- function(i,data,des) {
  #sample WTs and set dosing amount
  data[,"amt"] <- wt$WT*data[,"dose"]
  #set parameter value
  mod <- mod %>% param(slice(post,i)) 
  #simulate and add replicate number to data frame
  mod %>% 
    mrgsim(data=data,tgrid=des,obsonly=TRUE) %>% 
    mutate(irep=i) 
}

#simulate by calling the function 1000 times
set.seed(9999)
out <- future_map_dfr(1:1000, sim, data=sc, des=des, .options=opt)

head(out)

#summarize and plot
summ <- 
  out %>% 
  group_by(time) %>% 
  summarize(
    Q05=quantile(PCFB,prob=c(0.05)),
    Q50=quantile(PCFB,prob=c(0.50)),
    Q95=quantile(PCFB,prob=c(0.95))
  )

#Plot the distribution of percent change from baseline
ggplot(data=summ) + 
  geom_ribbon(aes(x=time,ymin=Q05,ymax=Q95),fill='cyan',alpha=0.5) + 
  geom_line(aes(x=time,y=Q50)) + 
  geom_hline(yintercept=-40,linetype="dashed") +
  ylab("Change in NTX (%)")

ggplot(data=(out %>% filter(time==2016))) + 
  geom_histogram(aes(x=PCFB),fill='cyan',bins=50) + 
  geom_vline(xintercept=-40,linetype='dashed') +
  xlab("Change in NTX at 3 Months (%)")

#calculate the fraction below the -40% threshold
out %>% filter(time==2016) %>% summarize(FRAC=mean(PCFB < -40))
