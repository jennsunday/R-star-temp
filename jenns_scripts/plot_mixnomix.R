mixnomix_all<-read_csv("data-processed/mixnomix_all.csv")

mixnomix_all$treatmentspecies<-mixnomix_all$species
mixnomix_all$treatment<-ifelse(mixnomix_all$treatmentspecies %in% c("mx1", "mx2", "mx3", "mx4"), "mx", "nm")
mixnomix_all$treatment<-as.factor(mixnomix_all$treatment)
mixnomix_all$species<-ifelse(mixnomix_all$treatmentspecies %in% c("mx1", "nm1"), 1, 
                             ifelse(mixnomix_all$treatmentspecies %in% c("mx2", "nm2"), 2, 
                                    ifelse(mixnomix_all$treatmentspecies %in% c("mx3", "nm3"), 3,
                                           ifelse(mixnomix_all$treatmentspecies %in% c("mx4", "nm4"), 4, NA))))
head(mixnomix_all)

mix10<-subset(mixnomix_all, mixnomix_all$temperature==10)
mix17<-subset(mixnomix_all, mixnomix_all$temperature==17)
mix24<-subset(mixnomix_all, mixnomix_all$temperature==24)

pdf("figures/mix_nomix.pdf", 6,4)
par(mfrow=c(1,3))
with(mix10, plot(cell_density~as.numeric(treatment), col=species, ylim=c(5000, 45000), main=temperature[1], xlab="", xaxt="n"))
for(i in 1:4){
  with(subset(mix10, mix10$species==i), abline(lm(cell_density~as.numeric(treatment)), col=i))
 }
axis(1, at=c(1,2), c("mix", "no mix"))

with(mix17, plot(cell_density~as.numeric(treatment), col=species, ylim=c(5000, 45000), main=temperature[1], xlab="", xaxt="n"))
for(i in 1:4){
   with(subset(mix17, mix17$species==i), abline(lm(cell_density~as.numeric(treatment)), col=i))
 }
axis(1, at=c(1,2), c("mix", "no mix"))

with(mix24, plot(cell_density~as.numeric(treatment), col=species, ylim=c(5000, 45000), main=temperature[1], xlab="", xaxt="n"))
for(i in 1:4){
  with(subset(mix24, mix24$species==i), abline(lm(cell_density~as.numeric(treatment)), col=i))
}
axis(1, at=c(1,2), c("mix", "no mix"))
dev.off()