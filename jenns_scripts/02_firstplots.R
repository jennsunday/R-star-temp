Rtemp_all<-read_csv("data-processed/Rtemp_all.csv")
names(Rtemp_all)

####add data from mix experiment#####
mixnomix_all<-read_csv("data-processed/mixnomix_all.csv")
mixnomix_all$treatmentspecies<-mixnomix_all$species
mixnomix_all$treatment<-ifelse(mixnomix_all$treatmentspecies %in% c("mx1", "mx2", "mx3", "mx4"), "mx", "nm")
mixnomix_all$treatment<-as.factor(mixnomix_all$treatment)
mixnomix_all$species<-ifelse(mixnomix_all$treatmentspecies %in% c("mx1", "nm1"), 1, 
                             ifelse(mixnomix_all$treatmentspecies %in% c("mx2", "nm2"), 2, 
                                    ifelse(mixnomix_all$treatmentspecies %in% c("mx3", "nm3"), 3,
                                           ifelse(mixnomix_all$treatmentspecies %in% c("mx4", "nm4"), 4, NA))))

names(mixnomix_all)
mixonly<-subset(mixnomix_all, mixnomix_all$treatment=="mx") #just add the mixed data
mixonly<-mixonly[, c(1:16)] #drop extra treatment columns

Rtemp_all<-rbind(Rtemp_all, mixonly)
#####################################

#plot log cell density through time
par(mfrow=c(2,2))

sp1data<-subset(Rtemp_all, Rtemp_all$species==1)
with(sp1data, plot(log(cell_density)~time_since_innoc_days, col=as.numeric(temperature), main=species[1]))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))

sp2data<-subset(Rtemp_all, Rtemp_all$species==2)
with(sp2data, plot(log(cell_density)~time_since_innoc_days, col=as.numeric(temperature), main=species[1]))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))

sp3data<-subset(Rtemp_all, Rtemp_all$species==3)
with(sp3data, plot(log(cell_density)~time_since_innoc_days, col=as.numeric(temperature), main=species[1]))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))

sp4data<-subset(Rtemp_all, Rtemp_all$species==4)
with(sp4data, plot(log(cell_density)~time_since_innoc_days, col=as.numeric(temperature), main=species[1]))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))


#plot cell density through time
par(mfrow=c(2,2))

sp1data<-subset(Rtemp_all, Rtemp_all$species==1)
with(sp1data, plot((cell_density)~time_since_innoc_days, col=as.numeric(temperature), main=species[1]))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))

sp2data<-subset(Rtemp_all, Rtemp_all$species==2)
with(sp2data, plot((cell_density)~time_since_innoc_days, col=as.numeric(temperature), main=species[1]))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))

sp3data<-subset(Rtemp_all, Rtemp_all$species==3)
with(sp3data, plot((cell_density)~time_since_innoc_days, col=as.numeric(temperature), main=species[1]))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))

sp4data<-subset(Rtemp_all, Rtemp_all$species==4)
with(sp4data, plot((cell_density)~time_since_innoc_days, col=as.numeric(temperature), main=species[1]))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))

names(Rtemp_all)



#plot cell volume through time
par(mfrow=c(2,2))

sp1data<-subset(Rtemp_all, Rtemp_all$species==1)
with(sp1data, plot(cell_volume~time_since_innoc_days, col=as.numeric(temperature), main=species[1], ylim=c(0,2000)))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))

sp2data<-subset(Rtemp_all, Rtemp_all$species==2)
with(sp2data, plot(cell_volume~time_since_innoc_days, col=as.numeric(temperature), main=species[1], ylim=c(0,2000)))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))

sp3data<-subset(Rtemp_all, Rtemp_all$species==3)
with(sp3data, plot(cell_volume~time_since_innoc_days, col=as.numeric(temperature), main=species[1], ylim=c(0,2000)))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))

sp4data<-subset(Rtemp_all, Rtemp_all$species==4)
with(sp4data, plot(cell_volume~time_since_innoc_days, col=as.numeric(temperature), main=species[1], ylim=c(0,2000)))
legend("topleft", pch=1, unique(sp1data$temperature), col=as.numeric(unique(sp1data$temperature)))


