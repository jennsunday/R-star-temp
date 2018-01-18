
TT_ks_umax<-cbind(read_csv("data-processed/TT_ks_umax_boot.csv"), Species="TT")
CS_ks_umax<-cbind(read_csv("data-processed/CS_ks_umax_boot.csv"), Species="CS")
AC_ks_umax<-cbind(read_csv("data-processed/AC_ks_umax_boot.csv"), Species="AC")
CH_ks_umax<-cbind(read_csv("data-processed/CH_ks_umax_boot.csv"), Species="CH")

all_ks_umax<-rbind(TT_ks_umax, CS_ks_umax, AC_ks_umax, CH_ks_umax)


all_ks_umax %>%
  #filter(Species=="CH") %>%
  #group_by(term) %>%
  ggplot(aes(x = Temperature, y = estimate, color = factor(Species))) + geom_point(size = 4) +
  theme_bw() + facet_wrap(~term, scales = "free") + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)
ggsave("figures/all_ks_umax.png")

# estimate R-star ------------
#R* = mKs / umax - m
#
#set m at 0.1
all_umax<- all_ks_umax %>%
  filter(term=="umax") 
all_ks<- all_ks_umax %>%
  filter(term=="ks") 
m<-0.1
all_umax$rstar<-m * all_ks$estimate / (all_umax$estimate - m)
all_rstar<- all_umax %>% 
  select(Temperature, run, Species, rstar) 
write_csv(all_rstar, "data-processed/all_rstar.csv")

all_rstar%>% 
  ggplot(aes(x = Temperature, y = rstar, color=factor(Species))) + geom_point(size = 4) +
  facet_wrap(~Species, scales = "free")





#next redo this by pulling 1000 iterations from umax and k estimates with norm dist with se
#write a function to draw random numbers with a boundary of 0:
rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

TT_rstar_norm_summ<-data.frame(median=1:6, lci=1:6, uci=1:6, Temperature=TT_ks$Temperature)
for(i in 1:length(TT_umax$estimate)){
  rstar_norm<-(rtnorm(n=100, mean=TT_ks$estimate[i], sd=TT_ks$std.error[i], a=0, b=Inf)*m)/
    (rnorm(100, mean=TT_umax$estimate[i], sd=TT_umax$std.error[i])-m)
  TT_rstar_norm_summ[i,c(1:3)]<-quantile(rstar_norm, c(0.5, 0.05, 0.95))
}
write_csv(TT_rstar_norm_summ, "data-processed/TT_rstar.csv")


all_rstar %>%
  ggplot(aes(x = Temperature, y = median, color = factor(Species))) + geom_point(size = 4) +
  theme_bw() + geom_line() + labs(y="R-star") +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.2)
ggsave("figures/all_rstar.png")
