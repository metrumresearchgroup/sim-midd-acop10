#' Packages
library(mrgsolve)
library(tidyverse)
library(furrr)

options(future.fork.enable=TRUE)
plan(multiprocess,workers=4L)
opt <- future_options(seed = TRUE)

#' Exercise 3:  Simulation with Uncertainty

#' Question
#' What dose would provide a 90% probability of achieving a 40% reduction in NTX 
#' with Q2W dosing?  With Q3W dosing?

#' OPG PK/PD model
mod <- mread("models/opg_pkpd.mod")
mod

#' Checkout the random effects structure
revar(mod)

#' set OMEGA and SIGMA matrices to zero to remove IIV and RUV
mod <- zero_re(mod)
revar(mod)

#' Load the simulated posterior
post <- read.table("data/post_POC.csv",header=T,sep=",") %>% sample_n(1000)

#' load database of patient weights and calculate mean to use in dose amount calculation
wt <- read.table("data/demographics.csv",header=T,sep=",") %>% summarize(WT=mean(WT))

#' Q2W doses
#' Input / template data set
sc <- expand.ev(ID=1, amt=0, cmt=1, dose=1:10, ii=336, addl=5)

#' tgrid object
des <- tgrid(0,2016,12)

#' Simulation function
sim <- function(i,data,des) {
  #sample WTs and set dosing amount
  data[,"amt"] <- wt$WT*data[,"dose"]
  #set parameter value
  mod <- mod %>% param(slice(post,i)) 
  #simulate and add replicate number to data frame
  mod %>% 
    carry_out(dose) %>% 
    mrgsim(data=data,tgrid=des,obsonly=TRUE) %>% 
    mutate(irep=i,dose=dose) 
}

set.seed(9999)
out <- future_map_dfr(1:1000, sim, data=sc, des=des, .options=opt)
head(out)

summ <- 
  out %>% 
  group_by(dose,time) %>% 
  summarize(
    Q10=quantile(PCFB,prob=c(0.10)),
    Q50=quantile(PCFB,prob=c(0.50)),
    Q90=quantile(PCFB,prob=c(0.90))
  )

#' Plot the distribution of percent change from baseline
ggplot(data=summ) + 
  geom_ribbon(aes(x=time,ymin=Q10,ymax=Q90),fill='cyan',alpha=0.5) + 
  geom_line(aes(x=time,y=Q50)) + 
  facet_wrap(~dose) + 
  geom_hline(yintercept=-50,linetype="dashed") 

below <- 
  out %>% 
  filter(time==2016) %>% 
  group_by(dose) %>% 
  summarize(FRAC=mean(PCFB < -40))

below

ggplot(data=below) + 
  geom_line(aes(x=dose,y=FRAC),color='red',size=1.5) + 
  geom_hline(yintercept=0.9,linetype="dashed") +
  xlab("Dose (mg/kg)") + 
  ylab("Probability of Achieving 40% NTX Reduction")



#' Q3W doses

#' Input / template data set
sc <- expand.ev(ID=1, amt=0, cmt=1, dose=seq(10,100,10), ii=168*3, addl=3)

#' tgrid object
des <- tgrid(0,2016,12)

out <- future_map_dfr(1:1000, sim, data=sc, des=des, .options=opt)
head(out)

summ <- 
  out %>% 
  group_by(dose,time) %>% 
  summarize(
    Q10=quantile(PCFB,prob=c(0.10)),
    Q50=quantile(PCFB,prob=c(0.50)),
    Q90=quantile(PCFB,prob=c(0.90))
  )

#' Plot the distribution of percent change from baseline
ggplot(data=summ) + 
  geom_ribbon(aes(x=time,ymin=Q10,ymax=Q90),fill='cyan',alpha=0.5) + 
  geom_line(aes(x=time,y=Q50)) + 
  facet_wrap(~dose) + 
  geom_hline(yintercept=-50,linetype="dashed") 

below <- out %>% filter(time==2016) %>% group_by(dose) %>% summarize(FRAC=mean(PCFB < -40))

ggplot(data=below) + 
  geom_line(aes(x=dose,y=FRAC),color='red',size=1.5) + 
  geom_hline(yintercept=0.9,linetype="dashed") +
  xlab("Dose (mg/kg)") + 
  ylab("Probability of Achieving 40% NTX Reduction")

