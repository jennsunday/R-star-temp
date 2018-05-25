pt3<-ggplot(sp1temp3, aes(date, absorbance.1)) + coord_cartesian(ylim = c(0, 0.15)) +
  geom_point(data=sp1temp3, colour='blue') +
  geom_point(data = sp1summ3, aes(y = mean), colour = 'blue', size = 3) + 
  geom_line(data = sp1summ3, aes(y = mean), colour = 'blue') +
  geom_point(data=sp2temp3, colour='blue') +
  geom_point(data = sp2summ3, aes(y = mean), colour = 'blue', size = 3) + 
  geom_line(data = sp2summ3, aes(y = mean), colour = 'blue') +
  geom_point(data=sp3temp3, colour='blue') +
  geom_point(data = sp3summ3, aes(y = mean), colour = 'blue', size = 3) + 
  geom_line(data = sp3summ3, aes(y = mean), colour = 'blue') +
  geom_point(data=sp4temp3, colour='blue') +
  geom_point(data = sp4summ3, aes(y = mean), colour = 'blue', size = 3) + 
  geom_line(data = sp4summ3, aes(y = mean), colour = 'blue')

pt10<-ggplot(sp1temp10, aes(date, absorbance.1)) + coord_cartesian(ylim = c(0, 0.15)) +
  geom_point(data=sp1temp10, colour='turquoise') +
  geom_point(data = sp1summ10, aes(y = mean), colour = 'turquoise', size = 3) + 
  geom_line(data = sp1summ10, aes(y = mean), colour = 'turquoise') +
  geom_point(data=sp2temp10, colour='turquoise') +
  geom_point(data = sp2summ10, aes(y = mean), colour = 'turquoise', size = 3) + 
  geom_line(data = sp2summ10, aes(y = mean), colour = 'turquoise') +
  geom_point(data=sp3temp10, colour='turquoise') +
  geom_point(data = sp3summ10, aes(y = mean), colour = 'turquoise', size = 3) + 
  geom_line(data = sp3summ10, aes(y = mean), colour = 'turquoise') +
  geom_point(data=sp4temp10, colour='turquoise') +
  geom_point(data = sp4summ10, aes(y = mean), colour = 'turquoise', size = 3) + 
  geom_line(data = sp4summ10, aes(y = mean), colour = 'turquoise')

pt17<-ggplot(sp1temp17, aes(date, absorbance.1)) + coord_cartesian(ylim = c(0, 0.15)) +
  geom_point(data=sp1temp17, colour='green') +
  geom_point(data = sp1summ17, aes(y = mean), colour = 'green', size = 3) + 
  geom_line(data = sp1summ17, aes(y = mean), colour = 'green') +
  geom_point(data=sp2temp17, colour='green') +
  geom_point(data = sp2summ17, aes(y = mean), colour = 'green', size = 3) + 
  geom_line(data = sp2summ17, aes(y = mean), colour = 'green') +
  geom_point(data=sp3temp17, colour='green') +
  geom_point(data = sp3summ17, aes(y = mean), colour = 'green', size = 3) + 
  geom_line(data = sp3summ17, aes(y = mean), colour = 'green') +
  geom_point(data=sp4temp17, colour='green') +
  geom_point(data = sp4summ17, aes(y = mean), colour = 'green', size = 3) + 
  geom_line(data = sp4summ17, aes(y = mean), colour = 'green')

pt24<-ggplot(sp1temp24, aes(date, absorbance.1)) + coord_cartesian(ylim = c(0, 0.15)) +
  geom_point(data=sp1temp24, colour='purple') +
  geom_point(data = sp1summ24, aes(y = mean), colour = 'purple', size = 3) + 
  geom_line(data = sp1summ24, aes(y = mean), colour = 'purple') +
  geom_point(data=sp2temp24, colour='purple') +
  geom_point(data = sp2summ24, aes(y = mean), colour = 'purple', size = 3) + 
  geom_line(data = sp2summ24, aes(y = mean), colour = 'purple') +
  geom_point(data=sp3temp24, colour='purple') +
  geom_point(data = sp3summ24, aes(y = mean), colour = 'purple', size = 3) + 
  geom_line(data = sp3summ24, aes(y = mean), colour = 'purple') +
  geom_point(data=sp4temp24, colour='purple') +
  geom_point(data = sp4summ24, aes(y = mean), colour = 'purple', size = 3) + 
  geom_line(data = sp4summ24, aes(y = mean), colour = 'purple')

pt31<-ggplot(sp1temp31, aes(date, absorbance.1)) + coord_cartesian(ylim = c(0, 0.15)) +
  geom_point(data=sp1temp31, colour='orange') +
  geom_point(data = sp1summ31, aes(y = mean), colour = 'orange', size = 3) + 
  geom_line(data = sp1summ31, aes(y = mean), colour = 'orange') +
  geom_point(data=sp2temp31, colour='orange') +
  geom_point(data = sp2summ31, aes(y = mean), colour = 'orange', size = 3) + 
  geom_line(data = sp2summ31, aes(y = mean), colour = 'orange') +
  geom_point(data=sp3temp31, colour='orange') +
  geom_point(data = sp3summ31, aes(y = mean), colour = 'orange', size = 3) + 
  geom_line(data = sp3summ31, aes(y = mean), colour = 'orange') +
  geom_point(data=sp4temp31, colour='orange') +
  geom_point(data = sp4summ31, aes(y = mean), colour = 'orange', size = 3) + 
  geom_line(data = sp4summ31, aes(y = mean), colour = 'orange')

pt38<-ggplot(sp1temp38, aes(date, absorbance.1)) + coord_cartesian(ylim = c(0, 0.15)) +
  geom_point(data=sp1temp38, colour='red') +
  geom_point(data = sp1summ38, aes(y = mean), colour = 'red', size = 3) + 
  geom_line(data = sp1summ38, aes(y = mean), colour = 'red') +
  geom_point(data=sp2temp38, colour='red') +
  geom_point(data = sp2summ38, aes(y = mean), colour = 'red', size = 3) + 
  geom_line(data = sp2summ38, aes(y = mean), colour = 'red') +
  geom_point(data=sp3temp38, colour='red') +
  geom_point(data = sp3summ38, aes(y = mean), colour = 'red', size = 3) + 
  geom_line(data = sp3summ38, aes(y = mean), colour = 'red') +
  geom_point(data=sp4temp38, colour='red') +
  geom_point(data = sp4summ38, aes(y = mean), colour = 'red', size = 3) + 
  geom_line(data = sp4summ38, aes(y = mean), colour = 'red')


grid.arrange(pt3, pt10, pt17, pt24, pt31, pt38, ncol = 3)
