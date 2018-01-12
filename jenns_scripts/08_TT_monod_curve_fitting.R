library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

TT <- read_csv("data/monod_data/TT_15_08_24.csv")

#plot raw data 
TT %>% 
	mutate(day = Hours.since.Innoc/24) %>% 
	#mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
	ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
	facet_wrap( ~ Temperature, scales = "free") + geom_line() + geom_abline(slope = 1, intercept = 5)


#filter out data where K seems to have been reached beyond the variability expected by noise (quite arbitrary)
TTfiltered<- TT %>% 
	mutate(day = Hours.since.Innoc/24) %>% 
	filter( Temperature == 16 & N.Treatment == 8 |
				 	Temperature == 16 & N.Treatment == 7 |
				 	Temperature == 16 & N.Treatment == 6 & day<6.5 | 
				 	Temperature == 16 & N.Treatment == 5 & day<6.5 |
				 	Temperature == 16 & N.Treatment == 4 & day<5 |
				 	Temperature == 16 & N.Treatment == 3 |
				 	Temperature == 16 & N.Treatment == 2 |
				 	Temperature == 16 & N.Treatment == 1 |
				 	Temperature == 16 & N.Treatment == 0 |
						Temperature == 19 & N.Treatment == 8  |
						Temperature == 19 & N.Treatment == 7  |
						Temperature == 19 & N.Treatment == 6 | 
						Temperature == 19 & N.Treatment == 5 |
						Temperature == 19 & N.Treatment == 4 |
						Temperature == 19 & N.Treatment == 3 & day<5 |
						Temperature == 19 & N.Treatment == 2 |
						Temperature == 19 & N.Treatment == 1 |
						Temperature == 19 & N.Treatment == 0 |
					Temperature == 22 & N.Treatment == 8 | 
					Temperature == 22 & N.Treatment == 7 & day<4.1 |
					Temperature == 22 & N.Treatment == 6 & day<4.1 | 
					Temperature == 22 & N.Treatment == 5 & day<4.5 |
					Temperature == 22 & N.Treatment == 4 & day<4.5 |
					Temperature == 22 & N.Treatment == 3 & day<4.5 |
					Temperature == 22 & N.Treatment == 2 |
					Temperature == 22 & N.Treatment == 1 |
					Temperature == 22 & N.Treatment == 0 | 		
						  Temperature == 25 & N.Treatment == 8 |
						  Temperature == 25 & N.Treatment == 7 |
						  Temperature == 25 & N.Treatment == 6 & day<4.5 | 
						  Temperature == 25 & N.Treatment == 5 & day<4.5 |
						  Temperature == 25 & N.Treatment == 4 |
						  Temperature == 25 & N.Treatment == 3 |
						  Temperature == 25 & N.Treatment == 2 |
						  Temperature == 25 & N.Treatment == 1 |
						  Temperature == 25 & N.Treatment == 0 | 		
								Temperature == 28 & N.Treatment == 8 & day<4.1 |
								Temperature == 28 & N.Treatment == 7 & day<4.1 |
								Temperature == 28 & N.Treatment == 6 & day<4.1 | 
								Temperature == 28 & N.Treatment == 5 & day<4.5 |
								Temperature == 28 & N.Treatment == 4 |
								Temperature == 28 & N.Treatment == 3 |
								Temperature == 28 & N.Treatment == 2 |
								Temperature == 28 & N.Treatment == 1 |
								Temperature == 28 & N.Treatment == 0 |
									Temperature == 13|
									Temperature == 10) 

#plot filtered data
TTfiltered %>% 
	mutate(day = Hours.since.Innoc/24) %>% 
	#mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
	ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
	facet_wrap( ~ Temperature, scales = "free") + geom_line() + geom_abline(slope = 1, intercept = 5)


#take filtered data, add actual N concentrations, plot monod curves
TTfilteredN<- TTfiltered %>% 	
	mutate(N.Treatment = str_replace(N.Treatment, "1", "11")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "2", "22")) %>% 
	mutate(N.Treatment = str_replace(N.Treatment, "3", "33")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "4", "55")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "5", "110")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "6", "220")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "7", "330")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "8", "440")) %>%  
	mutate(N.Treatment = str_replace(N.Treatment, "1105", "55")) %>% 
	mutate(N.Treatment = as.numeric(N.Treatment))


TTfilteredN %>% 
	group_by(Temperature, N.Treatment) %>% 
	do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
	geom_line() + theme_bw() + facet_wrap( ~ Temperature)
ggsave("figures/TT_monod_curves.png")

TT_r <- TTfilteredN %>% 
	group_by(Temperature, N.Treatment) %>% 
	do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) 

write_csv(TT_r, "data-processed/fitted_r_TT_from_2015.csv")


TTfilteredN %>% 
	group_by(Temperature, N.Treatment) %>% 
	do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ggplot(aes(x = Temperature, y = estimate, color = factor(Temperature))) + geom_point(size = 4) +
	geom_line() + theme_bw() + facet_wrap( ~ N.Treatment)
ggsave("figures/TT_TPC_by_nitrate_curves.png")

#set up the jacknifing
for (j in c(13, 16, 19, 22, 25, 28)){
  for (k in c(0, 11, 22, 33, 55, 110, 220, 440)){
    test<-subset(TTfilteredN, TTfilteredN$Temperature==j & TTfilteredN$N.Treatment==k)
    boot_r<-1:100
    for(i in 1:100){
      mod<-sample_n(test, length(test$day)-1) %>%
        do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
                    data= .,  start=list(a=0.01),
                    control = nls.control(maxiter=100, minFactor=1/204800000))))
      boot_r[i]<-mod$estimate
      name<-data.frame(a=unique(boot_r), Temperature=test$Temperature[1], N.Treatment=test$N.Treatment[1])
      assign(paste("boot_r_unique", test$Temperature[1], test$N.Treatment[1], sep ="_"), name)
    }
  }
}

# bind all of the objects that start with boot_r_unique_!
allboots<-do.call("rbind", mget(ls(pattern="boot_r_unique")))

#plot the jacknifed data - monod
allboots %>% 
  group_by(Temperature, N.Treatment) %>% 
  ggplot(aes(x = N.Treatment, y = a, color = factor(N.Treatment))) + geom_point(size = 2) +
  geom_line() + theme_bw() + facet_wrap( ~ Temperature)
ggsave("figures/TT_monod_jacknifed.png")