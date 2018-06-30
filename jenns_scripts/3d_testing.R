library(ggplot2)
library(plotly)
library(tidyverse)

#scatterplot
indirect_r<-read_csv("data-processed/indirect_r.csv")
indirect_r_AC<-filter(indirect_r, Species=="AC")
p <- plot_ly(indirect_r_AC, x = ~Temperature, y = ~N.Treatment, z = ~estimate, 
             color = ~Temperature, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Temperature'),
                      yaxis = list(title = 'Growth rate'),
                      zaxis = list(title = 'Resource')))
p


#smooth plot of model
#write function
growth_function<-function (N,Temp,a,b,z,w,ks) {
  (a*exp(b*Temp)*(1-((Temp-z)/(w/2))^2)) * ((N )/(ks + N)) 
}

predict_grid<-expand.grid(N=seq(0, 10, 0.1), Temp=seq(8, 33, 1))

NB_fit_AC<-read_csv("data-processed/all_params_fit_NB_fixed.csv") %>% filter(Species=="AC")


predict_grid_full <- predict_grid %>%
  mutate(r_pred=growth_function(predict_grid$N,
                                predict_grid$Temp,
                                filter(NB_fit_TT, term=="a")$estimate,
                                filter(NB_fit_TT, term=="b")$estimate,
                                filter(NB_fit_TT, term=="z")$estimate,
                                filter(NB_fit_TT, term=="w")$estimate,
                                filter(NB_fit_TT, term=="ks")$estimate))
head(predict_grid_full)

predict_grid_full %>%
  ggplot(aes(x=Temp, y=N, z=r_pred)) + geom_raster(aes(fill = r_pred)) + 
  scale_fill_gradientn(colours = c("#00007F", "red","yellow")) +
  stat_contour(binwidth=0.05)
ggplotly()

z<-predict_grid_full$r_pred
dim(z)<-c(length(unique(predict_grid$N)), length(unique(predict_grid$Temp)))

p <- plot_ly(predict_grid_full, x = ~Temp, y = ~N, z= ~ r_pred) %>% 
  add_surface(z = ~z)

p<-plot_ly(predict_grid_full, x = ~Temp, y = ~N, z = ~r_pred, type = 'mesh3d') 
l <- plotly_build(p)
l$data[[1]]$color <- toRGB("grey")
l


p<-plot_ly(predict_grid_full, x = ~Temp, y = ~N, z = ~r_pred, type = 'heatmap') 

p <- plot_ly() %>% 
  add_surface(z = ~z)

p <- plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~z)




library(plotly)

axx <- list(
  nticks = 4,
  range = c(-25,75)
)

axy <- list(
  nticks = 4,
  range = c(-25,75)
)

axz <- list(
  nticks = 4,
  range = c(0,50)
)
x <- 70*(runif(70, 0, 1))
y <- 55*(runif(70, 0, 1))
z <- 40*(runif(70, 0, 1))

p <- plot_ly(x = ~x, y = ~y, z = ~z, type = 'mesh3d') 

fdejong <- function (x, y) {
  return (x^2 + y^2)
}

x <- seq(-10, 10, length= 30)
y <- x
z <- outer(x, y, fdejong)
z[is.na(z)] <- 1

require(lattice)
wireframe(z, drape=T, col.regions=rainbow(100))
