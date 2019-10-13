# Packages
library(mrgsolve)
library(tidyverse)
library(parallel)

# Exercise 4:  Probabilistic Simulation

# Question
# If a 400 mg Q2W dosing regimen is taken into Phase 3,  what is the probability 
# that an individual patient will achieve a 40% reduction in NTX at trough?

mod <- mread("models/opg_pkpd.mod")
mod

#Checkout the random effects structure
revar(mod)

#zero out sigma and collapse omega matrices
mod <- mod %>% zero_re(sigma)
mod <- mod %>% mrgsolve:::collapse_omega() 
revar(mod)

# Load the simulated posterior
post <- read.table("data/post.csv",header=T,sep=",") %>% sample_n(1000)

omegas <- as_bmat(post,"OMEGA")
omegas[[1]]

#load database of patient weights and sample 1000 values to use in dose amount calculation
demo <- read.table("data/demographics.csv",header=T,sep=",") %>% sample_n(1000)

# Input / template data set
#Here is a simple dosing data set to simulate 1 patient
sc <- expand.ev(ID=1, amt=400, cmt=1, ii=336, ss=1)
sc

#We will get the observation design for the simulation through a `tgrid` object
des <- tgrid(end=-1,add=seq(0,336,6))
des

# Replicate simulation with uncertainty
#When we do replicate simulation, it almost always pays off
#to construct a function that will carry out one replicate.

#Arguments to the function are
#- `i` the current simulation replicate
#- `data` the dosing data set
#- `des` the observation design

sim <- function(i,data,des) {
  #sample WT
  data <- data %>% mutate(WT=demo$WT[i])
  #set parameter value
  mod <- mod %>% param(slice(post,i)) %>% omat(omegas[[i]])
  #simulate and add replicate number to data frame
  mod %>% Req(PCFB) %>% mrgsim(data=data,tgrid=des,obsonly=TRUE,recsort=3) %>% mutate(irep=i) 
}

#simulate by calling the function 1000 times
options(mc.cores=4)
mcRNG()
set.seed(9999)
out <- mclapply(1:1000, sim, data=sc, des=des) %>% bind_rows
head(out)

#summarize and plot
summ <- out %>% group_by(time) %>% summarize(Q05=quantile(PCFB,prob=c(0.05)),
                                             Q50=quantile(PCFB,prob=c(0.50)),
                                             Q95=quantile(PCFB,prob=c(0.95)))
ggplot(data=summ) + geom_ribbon(aes(x=time,ymin=Q05,ymax=Q95),fill='cyan',alpha=0.5) + geom_line(aes(x=time,y=Q50)) +
                    geom_hline(yintercept=-40,linetype="dashed") + xlab("Time (hours)") + ylab("Change in NTX (%)") 


ggplot(data=(out %>% filter(time==336))) + geom_histogram(aes(x=PCFB),fill='cyan',bins=50) + 
                                           geom_vline(xintercept=-40,linetype='dashed') +
                                           xlab("Change in NTX at Trough (%)")

#calculate the fraction below threshold
fraction_below <- function(x,boundary) { length(x[x<boundary])/length(x) }
out %>% filter(time==336) %>% summarize(FRAC=fraction_below(PCFB,-40))


