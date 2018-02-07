#Goal: plot all temperatures for a given species
#need to refit after grouping by replicate

#read in data
Rtemp_all<-read.csv("data-processed/Rtemp_all_manualedit.csv")
Rtemp_all$P<- Rtemp_all$cell_density
#fix errors in temperature treatment names
Rtemp_all$temperature[Rtemp_all$temperature=="1"]<-"3"
Rtemp_all$temperature[Rtemp_all$temperature=="30"]<-"31"
Rtemp_all$P<- Rtemp_all$cell_density
Rtemp_all$temperature<-as.numeric(Rtemp_all$temperature)


unique(Rtemp_all$replicate)

par(mfcol=c(4,6))
par(mar=c(2,2,0.5, 0.5), oma=c(0,0,0,0))
resultsr<-1:4
resultsk<-1:4


#for species 1, temp=3, use:
Parameters <- c(r = 4, K = 15000)
LowerBound <- c(r = 0.01, K = 10000)
UpperBound <- c(r = 5, K = 20000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,75), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[2] 
                    & Rtemp_all$species==1 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
  }
results31<-data.frame(fitr=resultsr, fitk=resultsk, species=1, rep=1:5)

#for species 2, temp=3, use
Parameters <- c(r = 0.5, K = 15000)
LowerBound <- c(r = 0, K = 100)
UpperBound <- c(r = 1, K = 16000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[2] 
                    & Rtemp_all$species==2 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
#save data, convert r to a negative value
results32<-data.frame(fitr=resultsr*(-1), fitk=resultsk, species=2, rep=1:5)

#for species 3, temp=3, use
Parameters <- c(r = 0.5, K = 15000)
LowerBound <- c(r = 0, K = 100)
UpperBound <- c(r = 1, K = 16000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[2] 
                    & Rtemp_all$species==3 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

#save data, convert r to a negative value
results33<-data.frame(fitr=resultsr*(-1), fitk=resultsk, species=3, rep=1:5)

#for species 4, temp=3, use
Parameters <- c(r = 2, K = 15000)
LowerBound <- c(r = 0, K = 100)
UpperBound <- c(r = 4, K = 18000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[2] 
                    & Rtemp_all$species==4 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results34<-data.frame(fitr=resultsr, fitk=resultsk, species=4, rep=1:5)

#combine all results for temp=3
results3<-rbind(results31, results32, results33, results34) %>%
  rename(r=fitr, k=fitk) %>%
  mutate(temp=3)

#for species 1, temp=10, use:
Parameters <- c(r = 0.5, K = 25000)
LowerBound <- c(r = 0.13, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,75), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[1] 
                    & Rtemp_all$species==1 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results101<-data.frame(fitr=resultsr, fitk=resultsk, species=1, rep=1:5)

#for species 2, temp=10, use:
Parameters <- c(r = 0.1, K = 25000)
LowerBound <- c(r = 0.09, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[1] 
                    & Rtemp_all$species==2 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results102<-data.frame(fitr=resultsr, fitk=resultsk, species=2, rep=1:5)

#for species 3, temp=10, use:
Parameters <- c(r = 0.1, K = 25000)
LowerBound <- c(r = 0.08, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[1] 
                    & Rtemp_all$species==3 & Rtemp_all$replicate==j & time_since_innoc_days<50)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results103<-data.frame(fitr=resultsr, fitk=resultsk, species=3, rep=1:5)

#for species 4, temp=10, use:
Parameters <- c(r = 0.1, K = 25000)
LowerBound <- c(r = 0.08, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[1] 
                    & Rtemp_all$species==4 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results104<-data.frame(fitr=resultsr, fitk=resultsk, species=4, rep=1:5)

#combine all results for temp=10
results10<-rbind(results101, results102, results103, results104) %>%
  rename(r=fitr, k=fitk) %>%
  mutate(temp=10)

#for species 1, temp=17, use:
Parameters <- c(r = 0.5, K = 25000)
LowerBound <- c(r = 0.1, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[3] 
                    & Rtemp_all$species==1 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results171<-data.frame(fitr=resultsr, fitk=resultsk, species=1, rep=1:5)

#for species 1, temp=17, use:
Parameters <- c(r = 0.5, K = 25000)
LowerBound <- c(r = 0.1, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound


for(j in 1:length(unique(Rtemp_all$replicate))){
  plot(c(1,50), c(0,30000), type="n")
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[3] 
                    & Rtemp_all$species==2 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results172<-data.frame(fitr=resultsr, fitk=resultsk, species=2, rep=1:5)

#for species 3, temp=17, use:
Parameters <- c(r = 0.5, K = 25000)
LowerBound <- c(r = 0.1, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[3] 
                    & Rtemp_all$species==3 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results173<-data.frame(fitr=resultsr, fitk=resultsk, species=3, rep=1:5)

#for species 4, temp=17, use:
Parameters <- c(r = 0.5, K = 25000)
LowerBound <- c(r = 0.1, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[3] 
                    & Rtemp_all$species==4 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results174<-data.frame(fitr=resultsr, fitk=resultsk, species=4, rep=1:5)

#combine all results for temp=17
results17<-rbind(results171, results172, results173, results174) %>%
  rename(r=fitr, k=fitk) %>%
  mutate(temp=17)


#for species 1, temp=24, use:
Parameters <- c(r = 10, K = 18000)
LowerBound <- c(r = 0.8, K = 17000) 
UpperBound <- c(r = 10, K = 21000) 
ParamScaling <- 0.0001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[4] 
                    & Rtemp_all$species==1 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results241<-data.frame(fitr=resultsr, fitk=resultsk, species=1, rep=1:5)
results241<-filter(results241, rep!=4) #remove rep 4 because it didn't have initial sample

#for species 2, temp=24, use:
Parameters <- c(r = 0.5, K = 5000)
LowerBound <- c(r = 0.5, K = 3000)
UpperBound <- c(r = 2, K = 10000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[4] 
                    & Rtemp_all$species==2 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results242<-data.frame(fitr=resultsr, fitk=resultsk, species=2, rep=1:5)

#for species 3, temp=24, use:
Parameters <- c(r = 0.5, K = 25000)
LowerBound <- c(r = 0.1, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[4] 
                    & Rtemp_all$species==3 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results243<-data.frame(fitr=resultsr, fitk=resultsk, species=3, rep=1:5)

#for species 4, temp=24, use:
Parameters <- c(r = 0.5, K = 25000)
LowerBound <- c(r = 0.1, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[4] 
                    & Rtemp_all$species==4 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results244<-data.frame(fitr=resultsr, fitk=resultsk, species=4, rep=1:5)

#combine all results for temp=24
results24<-rbind(results241, results242, results243, results244) %>%
  rename(r=fitr, k=fitk) %>%
  mutate(temp=24)


#for species 1, temp=31, use: (not the all r converged but point is negative growth)
Parameters <- c(r = 0.5, K = 300)
LowerBound <- c(r = 0.06, K = 0.00001)
UpperBound <- c(r = 0.7, K = 1300) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[5] 
                    & Rtemp_all$species==1 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results311<-data.frame(fitr=resultsr*(-1), fitk=resultsk, species=1, rep=1:5)

#for species 2, temp=31, use: (not the all r converged but point is negative growth)
Parameters <- c(r = 0.5, K = 300)
LowerBound <- c(r = 0.009, K = 50)
UpperBound <- c(r = 0.4, K = 1300) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[5] 
                    & Rtemp_all$species==2 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results312<-data.frame(fitr=resultsr*(-1), fitk=resultsk, species=2, rep=1:5)

#for species 3, temp=31, use: (not the all r converged but point is negative growth)
Parameters <- c(r = 0.5, K = 300)
LowerBound <- c(r = 0.06, K = 0.00001)
UpperBound <- c(r = 0.5, K = 1300) 
ParamScaling <- 0.001 / UpperBound


plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[5] 
                    & Rtemp_all$species==3 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results313<-data.frame(fitr=resultsr*(-1), fitk=resultsk, species=3, rep=1:5)

#for species 4, temp=31, use: (not the all r converged but point is negative growth)
Parameters <- c(r = 0.5, K = 300)
LowerBound <- c(r = 0.06, K = 0.00001)
UpperBound <- c(r = 1, K = 3000) 
ParamScaling <- 0.001 / UpperBound


plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[5] 
                    & Rtemp_all$species==4 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results314<-data.frame(fitr=resultsr*(1), fitk=resultsk, species=4, rep=1:5)

#combine all results for temp=31
results31<-rbind(results311, results312, results313, results314) %>%
  rename(r=fitr, k=fitk) %>%
  mutate(temp=31)


#for species 1, temp=38, use:
Parameters <- c(r = 0.05, K = 300)
LowerBound <- c(r = 0.01, K = 0.00011)
UpperBound <- c(r = 1, K = 8000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[6] 
                    & Rtemp_all$species==1 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results381<-data.frame(fitr=resultsr*(-1), fitk=resultsk, species=1, rep=1:5)

#for species 2, temp=38, use: (not converged but captures negative growth)
Parameters <- c(r = 0.05, K = 300)
LowerBound <- c(r = 0.01, K = 1)
UpperBound <- c(r = 1, K = 8000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[6] 
                    & Rtemp_all$species==2 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results382<-data.frame(fitr=resultsr*(-1), fitk=resultsk, species=2, rep=1:5)

#for species 3, temp=38, use: (not converged but captures negative growth)
Parameters <- c(r = 0.05, K = 300)
LowerBound <- c(r = 0.01, K = 100)
UpperBound <- c(r = 1, K = 5000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[6] 
                    & Rtemp_all$species==3 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results383<-data.frame(fitr=resultsr*(-1), fitk=resultsk, species=3, rep=1:5)


#for species 4, temp=38, use: (not converged but captures negative growth)
Parameters <- c(r = 0.5, K = 300)
LowerBound <- c(r = 0.06, K = 0.00001)
UpperBound <- c(r = 1, K = 3000) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,30000), type="n")
for(j in 1:length(unique(Rtemp_all$replicate))){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[6] 
                    & Rtemp_all$species==4 & Rtemp_all$replicate==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results384<-data.frame(fitr=resultsr*(-1), fitk=resultsk, species=4, rep=1:5)

#combine all results for temp=38
results38<-rbind(results381, results382, results383, results384) %>%
  rename(r=fitr, k=fitk) %>%
  mutate(temp=38)


#
#
#
#

cell_results<-rbind(results3, results10, results17, results24, results31, results38)
write.csv(cell_results, file="data-processed/logistic_growth_fits_rep_r-star.csv")

cell_results<-read.csv("data-processed/logistic_growth_fits_rep_r-star.csv")

cell_results %>%
  mutate(r_corr=r+0.1) %>% #correct for an added death rate of 0.1 per day
  ggplot(aes(x = temp, y = r_corr, color=as.factor(species))) + geom_point(size = 2) +
  theme_bw() + facet_grid(~species) +
  geom_hline(yintercept = 0) + geom_smooth(method = 'loess')
ggsave("growth_rate_with_temp_2017.pdf", width=7, height=2.5)
