# Packages
library(mrgsolve)
library(tidyverse)
library(furrr)

options(future.fork.enable=TRUE)
plan(multiprocess,workers=4L)
opt <- future_options(seed = TRUE)

# Exercise 5:  Population Simulation

# Question
# Using your competitor’s fracture model, what is the probability that fracture risk 
# will be reduced by 4.6% in the OPG-treated arm relative to placebo?
# 600 patients receiving 400 mg Q2W

mod <- mread("models/opg_fracture.mod")
mod

#Checkout the random effects structure
revar(mod)

# collapse omega and sigma matrices
mod <- mod %>% mrgsolve:::collapse_omega() 
mod <- mod %>% zero_re(sigma)  #don't want "observation error" in NTX to impact fracture probabilities
revar(mod)

# Load the simulated posterior
post <- read.table("data/post.csv",header=T,sep=",") %>% sample_n(1000)

omegas <- as_bmat(post,"OMEGA")
omegas[[1]]
sigmas <- as_bmat(post,"SIGMA")
sigmas[[1]]

#load database of patient weights and sample 600 values
demo <- read.table("data/demographics.csv",header=T,sep=",")

# Input / template data set
#Here is a simple dosing data set to simulate 1 patient
sc <- expand.ev(ID=1:600, amt=400, cmt=1, ii=336, ss=1)
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
  #grab 600 patients demographics and merge into data
  dem600 <- demo %>% sample_n(600) %>% mutate(ID=row_number())
  data <- data %>% left_join(dem600,by="ID")
  #set parameter values
  mod <- mod %>% param(slice(post,i)) %>% omat(omegas[[i]]) 
  #simulate one study arm
  res <- 
    mod %>% 
    Req(FRAC) %>% 
    mrgsim(data=data,tgrid=des,obsonly=TRUE,recsort=3,ss_n=100, ss_fixed=TRUE)
  
  #calculate mean and add replicated number to data frame
  out <- res %>% filter(time==336) %>% summarize(FRAC=mean(FRAC),irep=i)
  return(out)
}

#simulate by calling the function 1000 times
set.seed(9999)
out <- future_map_dfr(1:1000, sim.mean, data=sc, des=des, .options=opt)
head(out)

#simulate placebo fracture rate
sc <- sc %>% mutate(amt=0)
placebo <- future_map_dfr(1:1000, sim.mean, data=sc, des=des, .options=opt)
head(placebo)

#merge results and calculate fracture rate vs placebo for each study
res <- out %>% mutate(PLACEBO=placebo$FRAC,DIFF=FRAC-PLACEBO)
head(res)

#calculate fraction of studies with fracture rate at least 4.6% better than placebo
res %>% summarize(FRAC=mean(DIFF < -0.046))

