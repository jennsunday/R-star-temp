#Goal: plot all temperatures for a given species

#figure out outliers at temp=31 - messing up curve fits at crucial part of curve...

#to do: investigate spurious outliers and refit models.
#to do: use this model to fit decay in N

par(mfcol=c(4,6))
par(mar=c(2,2,0.5, 0.5), oma=c(0,0,0,0))
resultsr<-1:4
resultsk<-1:4


#for species 1-4, temp 1, use:
Parameters <- c(r = 0.5, K = 15000)
LowerBound <- c(r = 0, K = 100)
UpperBound <- c(r = 1, K = 16000) 
ParamScaling <- 0.001 / UpperBound


for(i in 1:4){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[2] 
                    & Rtemp_all$species==i)
  plotsinglefit(curvedata)
  resultsr[i]<-controlfit(curvedata)$r[1]
  resultsk[i]<-controlfit(curvedata)$K[1]
}

results3<-data.frame(r=resultsr*(c(1, -1, -1, 1)), K=resultsk, temp=rep(3, 4), species=1:4)

#for species i, temp 2, use:
Parameters <- c(r = 0.5, K = 25000)
LowerBound <- c(r = 0.1, K = 10 ^ 2)
UpperBound <- c(r = 20, K = 40000) 
ParamScaling <- 0.001 / UpperBound

#plot all temperatures for a given species (just changing speices # by hand for now)
for(i in 1:4){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[1] & Rtemp_all$species==i)
  plotsinglefit(curvedata)
  resultsr[i]<-controlfit(curvedata)$r[1]
  resultsk[i]<-controlfit(curvedata)$K[1]
}
results10<-data.frame(r=resultsr, K=resultsk, temp=rep(10, 4), species=1:4)

#for species i, temp 3, use:
Parameters <- c(r = 0.5, K = 25000)
LowerBound <- c(r = 0.1, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound

#plot all temperatures for a given species (just changing speices # by hand for now)
for(i in 1:4){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[3] 
                    & Rtemp_all$species==i & time_since_innoc_days<40)
  plotsinglefit(curvedata)
  resultsr[i]<-controlfit(curvedata)$r[1]
}
results17<-data.frame(r=resultsr, K=resultsk, temp=rep(17, 4), species=1:4)

#for species 1, temp 4, use: (it wants to be a very low r and really high k, completely off of the data)
Parameters <- c(r = 10, K = 20000)
LowerBound <- c(r = 0.8, K = 17000) 
UpperBound <- c(r = 10, K = 20000) 
ParamScaling <- 0.0001 / UpperBound

#plot all temperatures for a given species (just changing speices # by hand for now)
for(i in 1){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[4] 
                    & Rtemp_all$species==i)
  plotsinglefit(curvedata)
  resultsr[i]<-controlfit(curvedata)$r[1]
  resultsk[i]<-controlfit(curvedata)$K[1]
}

#for species 2, temp 4, use: (it wants to be a really high r - I just let it go)
Parameters <- c(r = 0.5, K = 5000)
LowerBound <- c(r = 0.5, K = 3000)
UpperBound <- c(r = 10, K = 10000) 
ParamScaling <- 0.001 / UpperBound

#plot all temperatures for a given species (just changing speices # by hand for now)
for(i in 2){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[4] 
                    & Rtemp_all$species==i)
  plotsinglefit(curvedata)
  resultsr[i]<-controlfit(curvedata)$r[1]
  resultsk[i]<-controlfit(curvedata)$K[1]
}

#for species 3-4, temp 4, use:
Parameters <- c(r = 0.5, K = 25000)
LowerBound <- c(r = 0.1, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 40000) 
ParamScaling <- 0.001 / UpperBound

#plot all temperatures for a given species (just changing speices # by hand for now)
for(i in 3:4){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[4] 
                    & Rtemp_all$species==i)
  plotsinglefit(curvedata)
  resultsr[i]<-controlfit(curvedata)$r[1]
  resultsk[i]<-controlfit(curvedata)$K[1]
}
results24<-data.frame(r=resultsr, K=resultsk, temp=rep(24, 4), species=1:4)


#for species 1-3, temp 5, use:
Parameters <- c(r = 0.5, K = 1000)
LowerBound <- c(r = 0, K = 0.1)
UpperBound <- c(r = 1, K = 1500) 
ParamScaling <- 0.001 / UpperBound

#plot all temperatures for a given species (just changing speices # by hand for now)
for(i in 1:3){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[5] 
                    & Rtemp_all$species==i)
  plotsinglefit(curvedata)
  resultsr[i]<-controlfit(curvedata)$r[1]
  resultsk[i]<-controlfit(curvedata)$K[1]
}

#for species 1-3, temp 5, use:
Parameters <- c(r = 0.5, K = 1000)
LowerBound <- c(r = 0, K = 0.1)
UpperBound <- c(r = 5, K = 3000) 
ParamScaling <- 0.001 / UpperBound

for(i in 4){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[5] 
                    & Rtemp_all$species==i)
  plotsinglefit(curvedata)
  resultsr[i]<-controlfit(curvedata)$r[1]
  resultsk[i]<-controlfit(curvedata)$K[1]
}
results31<-data.frame(r=resultsr*(c(-1, -1, -1, 1)), K=resultsk, temp=rep(31, 4), species=1:4)

#for species 4, temp 6, use:
Parameters <- c(r = 0.05, K = 2000)
LowerBound <- c(r = 0, K = 1)
UpperBound <- c(r = 1, K = 8000) 
ParamScaling <- 0.001 / UpperBound

#plot all temperatures for a given species (just changing speices # by hand for now)
for(i in 1:4){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[6] 
                    & Rtemp_all$species==i)
  plotsinglefit(curvedata)
  resultsr[i]<-controlfit(curvedata)$r[1]
  resultsk[i]<-controlfit(curvedata)$K[1]
}

results38<-data.frame(r=(-resultsr), K=resultsk, temp=rep(38, 4), species=1:4)

#
#
#
#

results<-rbind(results3, results10, results17, results24, results31, results38)
write.csv(results, file="data-processed/logistic_growth_fits_r-star.csv")

cell_results<-read.csv("data-processed/logistic_growth_fits_r-star.csv")

par(mfrow=c(1,1))
with(results, plot(r~temp, type="n", las=2))
for(i in 1:4){
with(subset(results, cell_results$species==i), points(r~temp, col=i))
with(subset(results, cell_results$species==i), lines(r~temp, col=i))  
}
abline(0,0)
mtext("r", 2, 3, las=1)
mtext("Temperature (°C)", 1, 3)
legend("topleft", c("TT", "CS", "AC", "Chlamy"), lty=1, col=1:4, pch=1, bty="n")

par(mfrow=c(1,1))
with(results, plot(K~temp, type="n", las=1))
for(i in 1:4){
  with(subset(results, results$species==i), points(K~temp, col=i))
  with(subset(results, results$species==i), lines(K~temp, col=i))  
}
abline(0,0)
mtext("K", 2, 3, las=1)
mtext("Temperature (°C)", 1, 3)
legend("topright", c("TT", "CS", "AC", "Chlamy"), lty=1, col=1:4, pch=1, bty="n")
