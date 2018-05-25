#predict outcomes of competition between TT and CH for different temperatures and dilution rates
library(rootSolve)
library(deSolve)
library(reshape2)
library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

###########################################################################
#competition

#differential equations with one consumer
dNRdt = function(t, state, parameters) {
	with(
		as.list(c(state,parameters)),{
			dN1 = N1 * ((u1 * R / (K1 + R)) - d) #consumer1
			dR = d * (S - R) - Q1 * 1000* ( dN1 + d * N1)  #resource
			return(list(c(dN1, dR)))			
		}
	)
}
	
	
#times
time = seq(from = 0, to = 500, by = 1) 

TT_ks_umax<-cbind(read.csv("data-processed/TT_ks_umax.csv"), Species="TT")
CH_ks_umax<-cbind(read.csv("data-processed/CH_ks_umax.csv"), Species="CH")
AC_ks_umax<-cbind(read.csv("data-processed/AC_ks_umax.csv"), Species="AC")
CS_ks_umax<-cbind(read.csv("data-processed/CS_ks_umax.csv"), Species="CS")

all_ks_umax<-rbind(TT_ks_umax, CS_ks_umax, AC_ks_umax, CH_ks_umax)

QCH<-28.74/((130000-5000)*1000)
QTT<-28.74/((90000-5000)*1000)

#parameters in experiment
Nlimiting<-4.5
Dexp<-0.1

state = c(N1 = 1E3, R = Nlimiting)

testks<-all_ks_umax %>%
  filter(term=="ks") %>%
  filter(Species=="TT") %>%
  filter(Temperature==28) %>% 
  .$estimate
testumax<-all_ks_umax %>%
  filter(term=="umax") %>%
  filter(Species=="TT") %>%
  filter(Temperature==28) %>% 
  .$estimate

testparameters = c(u1 = testumax, K1 = testks, Q1 =  QTT, d = Dexp, S = Nlimiting)


#input N concentration in medium and original number of cells 
out = ode(y = state, times = time, func = dNRdt, parms = testparameters)
out.df = as.data.frame(out) # required by ggplot: data object must be a data frame
out.m = melt(out.df, id.vars='time') # this makes plotting easier by puting all variables in a single column
with(subset(out.m, out.m$variable!="R"), plot(value  ~time, type="l", xlim=c(0, 60)))#plots in base R

with(subset(out.m, out.m$variable=="R"), plot(value  ~time, col=as.numeric(out.m$variable), xlim=c(0, 60), main="R, 15°C"))#plots in base R


#predict for 25°C
parameters25 = c(u1 = umaxTT25, K1 = KSTT25, Q1 =  QTT, u2 = umaxCH25, K2 = KSCH25, Q2 = QCH, d = Dexp, S = Nlimiting)

#input N concentration in medium and original number of cells 
out = ode(y = state, times = time, func = dNNRdt, parms = parameters25)

out.df = as.data.frame(out) # required by ggplot: data object must be a data frame
out.m = melt(out.df, id.vars='time') # this makes plotting easier by puting all variables in a single column
with(subset(out.m, out.m$variable!="R"), plot(log(value, 10)  ~time, col=as.numeric(out.m$variable), xlim=c(0, 60), pch=as.numeric(out.m$variable), main="25°C", ylim=c(0, 5)))#plots in base R


with(subset(out.m, out.m$variable=="R"), plot(value  ~time, col=as.numeric(out.m$variable), xlim=c(0, 60), main="R, 25°C"))#plots in base R
