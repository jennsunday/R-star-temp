#Code by Jennifer Sunday

#packages####

#load data####
data<-read.csv(file="./data/gapminder_2007.csv")

#load functions####
source(file = "./scripts/functions.R")

## Hi Jenn!!!!!!

#process data####
cv(data$gdpPercap)

#figures####
pdf(file="./figures/Figure 1.pdf", width = 8, height = 6)
plot(data$gdpPercap~data$lifeExp)
dev.off()

#save.workspace####
save.image(file = "./workspace/gapminer.RData")
load(file = "./workspace/gapminer.RData")

save.image(file = "./workspace/gapminer.RData")
load(file = "./workspace/gapminer.RData")

