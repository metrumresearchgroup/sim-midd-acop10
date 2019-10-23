#' Packages
library(mrgsolve)
library(tidyverse)
library(furrr)

options(future.fork.enable=TRUE)
plan(multiprocess,workers=4L)
opt <- future_options(seed = TRUE)

#' Exercise 4:  Probabilistic Simulation

#' Question
#' If two doses are taken into Phase 3,  one for patients < 70 kg and one for 
#' patients >= 70 kg, what doses are required for patients to have a 90% 
#' probability of achieving the target?

#' Read model
mod <- mread("models/opg_pkpd.mod")
mod

#' Checkout the random effects structure
revar(mod)

#' zero out sigma and collapse omega matrices
mod <- mod %>% zero_re(sigma)
mod <- mod %>% mrgsolve:::collapse_omega() 
revar(mod)

#' Load the simulated posterior
post <- read.table("data/post.csv",header=T,sep=",") %>% sample_n(1000)

omegas <- as_bmat(post,"OMEGA")
omegas[[1]]

#' Patients < 70 kg

#' load database of patient weights and sample 1000 values to use in dose amount 
#' calculation
demo <- read.table("data/demographics.csv",header=T,sep=",") %>% filter(WT < 70) %>% sample_n(1000)

#' Input / template data set
#' Here is a simple dosing data set to simulate 1 patient
sc <- expand.ev(ID=1, amt=seq(300,800,50), cmt=1, ii=336, ss=1)
sc <- sc %>% mutate(dose=amt)
sc

#' We will get the observation design for the simulation through a `tgrid` object
des <- tgrid(0,336,6)
des

#' Replicate simulation with uncertainty
#' When we do replicate simulation, it almost always pays off
#' to construct a function that will carry out one replicate.

#' Arguments to the function are
#' - `i` the current simulation replicate
#' - `data` the dosing data set
#' - `des` the observation design

sim <- function(i,data,des) {
  #sample WT
  data <- data %>% mutate(WT=demo$WT[i])
  #set parameter value
  mod <- mod %>% param(slice(post,i)) %>% omat(omegas[[i]])
  #simulate and add replicate number to data frame
  mod %>% 
    carry_out(dose) %>% 
    mrgsim(data=data,tgrid=des,obsonly=TRUE,recsort=3,
           ss_n = 50, ss_fixed=TRUE) %>% 
    mutate(irep=i) 
}

#' simulate by calling the function 1000 times
set.seed(1234)
out <- future_map_dfr(1:1000, sim, data=sc, des=des,.options=opt)
head(out)

#' summarize and plot
summ <- 
  out %>% 
  group_by(dose,time) %>% 
  summarize(
    Q20=quantile(PCFB,prob=c(0.20)),
    Q50=quantile(PCFB,prob=c(0.50)),
    Q80=quantile(PCFB,prob=c(0.80))
  )

ggplot(data=summ) + 
  geom_ribbon(aes(x=time,ymin=Q20,ymax=Q80),fill='cyan',alpha=0.5) + 
  geom_line(aes(x=time,y=Q50)) +
  geom_hline(yintercept=-40,linetype="dashed") +
  xlab("Time (hours)") + ylab("Change in NTX (%)") + 
  facet_wrap(~dose)

below <- out %>% filter(time==336) %>% group_by(dose) %>% summarize(FRAC=mean(PCFB < -40))

below

ggplot(data=below) + 
  geom_line(aes(x=dose,y=FRAC),color='red',size=1.5) + 
  geom_hline(yintercept=0.8,linetype='dashed') +
  xlab("Dose (mg)") + ylab("Probability dNTX < -40%") + 
  ggtitle("Patients < 70 kg")

#' WT >= 70

#' load database of patient weights and sample 1000 values to use in dose amount 
#' calculation
demo <- read.table("data/demographics.csv",header=T,sep=",") %>% filter(WT >= 70) %>% sample_n(1000)

set.seed(5678)
out <- future_map_dfr(1:1000, sim, data=sc, des=des, .options=opt)
head(out)

#' summarize and plot
summ <- 
  out %>% 
  group_by(dose,time) %>% 
  summarize(
    Q20=quantile(PCFB,prob=c(0.20)),
    Q50=quantile(PCFB,prob=c(0.50)),
    Q80=quantile(PCFB,prob=c(0.80))
  )

ggplot(data=summ) + 
  geom_ribbon(aes(x=time,ymin=Q20,ymax=Q80),fill='cyan',alpha=0.5) + 
  geom_line(aes(x=time,y=Q50)) +
  geom_hline(yintercept=-40,linetype="dashed") + 
  xlab("Time (hours)") + ylab("Change in NTX (%)") + 
  facet_wrap(~dose)

#' calculate the fraction of patients below the threshold
below <- out %>% filter(time==336) %>% group_by(dose) %>% summarize(FRAC=mean(PCFB < -40))
below

ggplot(data=below) + 
  geom_line(aes(x=dose,y=FRAC),color='red',size=1.5) + 
  geom_hline(yintercept=0.8,linetype='dashed') +
  xlab("Dose (mg)") + ylab("Probability dNTX < -40%") + 
  ggtitle("Patients >= 70 kg")
