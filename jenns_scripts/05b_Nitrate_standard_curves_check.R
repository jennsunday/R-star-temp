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
  mutate(date.measured="17/06/26", day=1)
new_2<-stcurve.2706 %>%
  select(sample, absorbance.1, absorbance.2, nitrate, rep) %>%
  mutate(date.measured="17/06/27", day=2)
new_3<-stcurve.2806 %>%
  select(sample, absorbance.1, absorbance.2, nitrate, rep) %>%
  mutate(date.measured="17/06/28", day=3)
new_4<-stcurve.1107 %>%
  select(sample, absorbance.1, absorbance.2, nitrate, rep) %>%
  mutate(date.measured="17/07/11", day=4)
new_5<-stcurve.1207 %>%
  select(sample, absorbance.1, absorbance.2, nitrate, rep) %>%
  mutate(date.measured="17/07/12", day=5)
new_6<-stcurve.1407 %>%
  select(sample, absorbance.1, absorbance.2, nitrate, rep) %>%
  mutate(date.measured="17/07/14", day=6)


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

all_standard_curves_vert %>%
  filter(absorbance<0.12)%>%
  ggplot(aes(y=nitrate, x=absorbance, col=rep)) +
  geom_point() + stat_smooth(method="lm") + facet_wrap(~date.measured)


all_standard_curves_vert %>%
  filter(absorbance<0.12)%>%
  ggplot(aes(y=nitrate, x=absorbance)) +
  geom_point() + stat_smooth(method="lm")

all_st_curves_model<-all_standard_curves_vert %>%
  filter(absorbance<0.12)%>%
  group_by(date.measured)%>%
  do(tidy(lm(nitrate~absorbance, data=.)))

all_st_curves_model  %>% 
  filter(term=="absorbance") %>%
  ggplot(aes(y=estimate, x=date.measured)) + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)




#plot standard curves from different dates on top of each other
with(stcurve.2606, plot(absorbance.1~nitrate,  xlim=c(0,5), ylim=c(-0.02,0.1)))
with(stcurve.2606, points(absorbance.2~nitrate, col=rep))
abline(lm(meanabs~nitrate, data=stcurve.2606))
summary(lm(meanabs~nitrate, data=stcurve.2606))

with(stcurve.2706, points(absorbance.1~nitrate, col=rep, pch=2))
with(stcurve.2706, points(absorbance.2~nitrate, col=rep, pch=2))
abline(lm(absorbance.1~nitrate, data=stcurve.2706))
summary(lm(meanabs~nitrate, data=stcurve.2706))

#remove erroneous reads - clearly out of line, especially after cadmium bottle was replaced
stcurve.2806.corr<-subset(stcurve.2806, rep=="C" | 
                            rep=="A" & !nitrate %in% c(4,5,6) | 
                            rep=="B" & !nitrate %in% c(4,7))
with(stcurve.2806.corr, points(meanabs~nitrate, col=rep, pch=3, ylim=c(0, 0.11)))

with(subset(stcurve.2806, stcurve.2806$rep %in% c("B2_500_a", "B2_500_b", "B2_cuvette")), 
            points(meanabs~nitrate, col=rep, pch=3))

with(stcurve.1107, points(meanabs~nitrate, col=rep, pch=1))
with(stcurve.1207, points(meanabs~nitrate, col=rep, pch=2))
with(stcurve.1407, points(meanabs~nitrate, col=1, pch=2))

legend("topleft", levels(stcurve.1107$rep), levels(as.factor((stcurve.1107$rep))), pch=1, col=1:4)


#define standard curve equations for each date
st.2606<-lm(meanabs~nitrate, data=stcurve.2606)
st.2706<-lm(meanabs~nitrate, data=stcurve.2706)

st.2806<-lm(meanabs~nitrate, data=stcurve.2806.corr)
st.1107<-lm(meanabs~nitrate, data=stcurve.1107)
st.1207<-lm(meanabs~nitrate, data=stcurve.1207)
st.1407<-lm(meanabs~nitrate, data=stcurve.1407)

#look at the lines
abline(st.2806, col=3)
abline(st.1107, col=1)
abline(st.1207, col=4)
abline(st.1407, col=4)

#plot slopes of this curve as a function of date
slopes<-c(coef(st.2606)[1], coef(st.2706)[1], coef(st.2806)[1], coef(st.1107)[1], coef(st.1207)[1], coef(st.1407)[1])
dates<-c(26-26, 27-26, 28-26, 41-26, 42-26, 44-26)
plot(slopes~dates)

#plot intercepts of this curve as a function of date
ints<-c(coef(st.2606)[2], coef(st.2706)[2], coef(st.2806)[2], coef(st.1107)[2], coef(st.1207)[2], coef(st.1407)[2])
dates<-c(26-26, 27-26, 28-26, 41-26, 42-26, 44-26)
plot(ints~dates)

#July 5,6,7 are the measurement dates that don't have a standard curve
#these would be 9,10,11

