# Packages
library(mrgsolve)
library(tidyverse)
library(parallel)

# Exercise 5:  Population Simulation

# Question
# What is the probability that at least a 43% mean reduction in NTX  would be observed at steady-state in the OPG-treated study arm?
# 600 patients receiving 400 mg Q2W

mod <- mread("models/opg_pkpd.mod")
mod

#Checkout the random effects structure
revar(mod)

# collapse omega and sigma matrices
mod <- mod %>% mrgsolve:::collapse_omega() 
mod <- mod %>% mrgsolve:::collapse_sigma() 
revar(mod)

# Load the simulated posterior
post <- read.table("data/post.csv",header=T,sep=",") %>% sample_n(1000)

omegas <- as_bmat(post,"OMEGA")
omegas[[1]]
sigmas <- as_bmat(post,"SIGMA")
sigmas[[1]]

#load database of patient weights and sample 600 values
demo <- read.table("data/demographics.csv",header=T,sep=",") %>% sample_n(600)

# Input / template data set
#Here is a simple dosing data set to simulate 1 patient
sc <- expand.ev(ID=1, amt=400, cmt=1, ii=336, ss=1, WT=demo$WT)
sc

#We will get the observation design for the simulation through a `tgrid` object
des <- tgrid(end=-1,add=c(0,336))
des

# Replicate simulation with uncertainty
#When we do replicate simulation, it almost always pays off
#to construct a function that will carry out one replicate.

#Arguments to the function are
#- `i` the current simulation replicate
#- `data` the dosing data set
#- `des` the observation design

sim.mean <- function(i,data,des) {
  #set parameter values
  mod <- mod %>% param(slice(post,i)) %>% omat(omegas[[i]]) %>% smat(sigmas[[i]])
  #simulate one study arm
  res <- mod %>% Req(PCFB) %>% mrgsim(data=data,tgrid=des,obsonly=TRUE,recsort=3)
  #calculate mean and add replicated number to data frame
  out <- res %>% filter(time==336) %>% summarize(MEAN=mean(PCFB)) %>% mutate(irep=i) %>% as.data.frame()
  return(out)
}

#simulate by calling the function 1000 times
options(mc.cores=4)
mcRNG()
set.seed(9999)
out <- mclapply(1:1000, sim.mean, data=sc, des=des) %>% bind_rows
head(out)

ggplot(data=out) + geom_histogram(aes(x=MEAN),fill='cyan',bins=50) + 
                   geom_vline(xintercept=-53,linetype='dashed') +
                   xlab("Mean Change in NTX at Trough (%)")

#calculate fraction of studies below threshold
fraction_below <- function(x,boundary) { length(x[x<=boundary])/length(x) }
out %>% summarize(FRAC=fraction_below(MEAN,-53))



######## Sensitivity Analysis ############

# pull THETA values out of posteriors and merge with mean change in NTX by replicate number
data <- post %>% dplyr::select(THETA1,THETA2,THETA3,THETA4,THETA5,THETA6,THETA7,THETA8,THETA9,THETA10,THETA11,THETA12,THETA13) %>% 
                 mutate(irep=row_number()) %>% left_join(out)

# plot posteriors vs mean change in NTX
ggplot(data=data,aes(x=THETA1,y=MEAN)) + geom_point() + geom_smooth() + xlab("Clearance") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA2,y=MEAN)) + geom_point() + geom_smooth() + xlab("Central Volume") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA3,y=MEAN)) + geom_point() + geom_smooth() + xlab("Peripheral Volume 1") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA4,y=MEAN)) + geom_point() + geom_smooth() + xlab("Peripheral Volume 2") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA5,y=MEAN)) + geom_point() + geom_smooth() + xlab("Distribution Clearance 1") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA6,y=MEAN)) + geom_point() + geom_smooth() + xlab("Distribution Clearance 2") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA7,y=MEAN)) + geom_point() + geom_smooth() + xlab("Absorption Rate") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA8,y=MEAN)) + geom_point() + geom_smooth() + xlab("Maximum Velocity") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA9,y=MEAN)) + geom_point() + geom_smooth() + xlab("Michaelis Constant") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA10,y=MEAN)) + geom_point() + geom_smooth() + xlab("Bioavailability") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA11,y=MEAN)) + geom_point() + geom_smooth() + xlab("NTX Synthesis Rate") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA12,y=MEAN)) + geom_point() + geom_smooth() + xlab("NTX Elimination Rate") + ylab("Mean Change in NTX (%)")
ggplot(data=data,aes(x=THETA13,y=MEAN)) + geom_point() + geom_smooth() + xlab("IC50") + ylab("Mean Change in NTX (%)")