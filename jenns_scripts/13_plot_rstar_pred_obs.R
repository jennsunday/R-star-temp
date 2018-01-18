
all_rstar_2105<-read_csv("data-processed/all_rstar.csv")
obs_rstar_2017<-read_csv("data-processed/logistic_N_decay_fits_r-star.csv")

with(obs_rstar_2017, plot(K~temp, type="n", las=1))
for(i in 1:4){
  with(subset(results, results$species==i), points(K~temp, col=i))
  with(subset(results, results$species==i), lines(K~temp, col=i))  
  with(subset(results, results$species==i), points(K_flat~temp, col=i, pch=2))
  with(subset(results, results$species==i), segments(temp, K_flat+K_flat_SE, temp, K_flat-K_flat_SE, col=i))
}

obs_rstar_2017 %>%
  rename(N_equib=K) %>%
  ggplot(aes(x = temp, y = K_flat, color=factor(species))) + geom_point(size = 4) + geom_line() +
  geom_errorbar(aes(ymin=K_flat-K_flat_SE, ymax=K_flat+K_flat_SE), width=.2) + 
  geom_point(data = all_rstar_2105, aes(x = Temperature, y = rstar, color=factor(Species)))