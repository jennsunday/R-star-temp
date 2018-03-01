library(ggplot2)
library(broom)
library(tidyverse)


#Read in standard curves

stcurve.2606<-read.csv("nitrate_data/standard_curves/standard_curve_260617.csv")
stcurve.2706<-read.csv("nitrate_data/standard_curves/standard_curve_270617.csv")
stcurve.2806<-read.csv("nitrate_data/standard_curves/standard_curve_280617.csv")
stcurve.1107<-read.csv("nitrate_data/standard_curves/standard_curve_110717.csv")
stcurve.1207<-read.csv("nitrate_data/standard_curves/standard_curve_120717.csv")
stcurve.1407<-read.csv("nitrate_data/standard_curves/standard_curve_140717.csv") #this one has a different strucutre

#remove abherrant point because there was a bubble on the cuvette
stcurve.2806$absorbance.1 [stcurve.2806$sample=="st.3." & stcurve.2806$rep=="B"]<- 0.044



head(stcurve.2606)

new_1<-stcurve.2606 %>%
  select(sample, absorbance.1, absorbance.2, nitrate, rep) %>%
  mutate(date.measured="17.06.26", day=1)
new_2<-stcurve.2706 %>%
  select(sample, absorbance.1, absorbance.2, nitrate, rep) %>%
  mutate(date.measured="17.06.27", day=2)
new_3<-stcurve.2806 %>%
  select(sample, absorbance.1, absorbance.2, nitrate, rep) %>%
  mutate(date.measured="17.06.28", day=3) %>%
  filter(rep!="B") #remembering that I immediately redid rep B on 17/06/28 as rep B2, therefore ok to remove rep B
new_4<-stcurve.1107 %>%
  select(sample, absorbance.1, absorbance.2, nitrate, rep) %>%
  mutate(date.measured="17.07.11", day=4)
new_5<-stcurve.1207 %>%
  select(sample, absorbance.1, absorbance.2, nitrate, rep) %>%
  mutate(date.measured="17.07.12", day=5)
new_6<-stcurve.1407 %>%
  select(sample, absorbance.1, absorbance.2, nitrate, rep) %>%
  mutate(date.measured="17.07.14", day=6)


all_standard_curves<-rbind(new_1, new_2, new_3, new_4, new_5, new_6)

#get 2 absorbance measures from each samples into a single column labelled as "take" 1 and 2
all_standard_curves_vert<-rbind(all_standard_curves, all_standard_curves) %>%
  mutate(absorbance=c(all_standard_curves$absorbance.1, all_standard_curves$absorbance.2), 
         take=rep(1:2, length(all_standard_curves$absorbance.2))) %>%
  select(-absorbance.1, -absorbance.2) %>%
  filter(!sample %in% c("blank", "Blank"))  %>%
  filter(!rep %in% c("B2_500_a", "B2_500_b", "B2_cuvette")) 

       
names(all_standard_curves_vert)
as.factor(all_standard_curves_vert$day)

#adventures in plotting
all_standard_curves_vert %>%
  filter(absorbance<0.12)%>%
  ggplot(aes(y=nitrate, x=absorbance)) +
  geom_point() + stat_smooth(method="lm") + facet_wrap(~date.measured) + theme_bw()
ggsave("standard_curves.pdf")

all_standard_curves_vert %>%
  filter(absorbance<0.12)%>%
  ggplot(aes(y=nitrate, x=absorbance)) +
  geom_point() + stat_smooth(method="lm")


#Fit model
all_st_curves_model<-all_standard_curves_vert %>%
  filter(absorbance<0.12)%>%
  group_by(date.measured, rep)%>%
  do(tidy(lm(nitrate~absorbance, data=.)))

predict(all_st_curves_model, newdata=0.5)

#fit a different model for different dates
for (i in 1:length(unique(all_standard_curves_vert$date.measured))){
temp<-all_standard_curves_vert %>%
  filter(absorbance<0.12)%>%
  filter(date.measured==unique(all_standard_curves_vert$date.measured)[i])
mod<-lm(nitrate~absorbance, data=temp)
assign(paste("mod", unique(all_standard_curves_vert$date.measured)[i], sep="_"), mod)
}

#these are the models:
mod_17.06.26
mod_17.06.27
mod_17.06.28
mod_17.07.11
mod_17.07.12
mod_17.07.14

all_st_curves_model  %>% 
  filter(term=="absorbance") %>%
  ggplot(aes(y=estimate, x=date.measured, col=rep)) + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)

