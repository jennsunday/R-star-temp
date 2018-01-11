library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

TT <- read_csv("data-raw/TT_2015.csv")

#plot raw data 
TTfiltered %>% 
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
	filter(Temperature == 28) %>% 
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
ggsave("figures/monod_curves.png")

TT_r <- TTfilteredN %>% 
	group_by(Temperature, N.Treatment) %>% 
	do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) 

write_csv(TT_r, "data-processed/fitted_r_TT_2018_nls.csv")


TTfilteredN %>% 
	group_by(Temperature, N.Treatment) %>% 
	do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ggplot(aes(x = Temperature, y = estimate, color = factor(Temperature))) + geom_point(size = 4) +
	geom_line() + theme_bw() + facet_wrap( ~ N.Treatment)
ggsave("figures/TPC_by_nitrate_curves.png")
