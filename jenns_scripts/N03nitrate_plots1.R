maxy<-7.5

###sp.1
psp1t3<-ggplot(sp1temp3, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy), xlim=c(0, 70)) +
  geom_point(data=sp1temp3, colour='blue') +
  geom_point(data = sp1summ3, aes(y = mean), colour = 'blue', size = 3) + 
  geom_line(data = sp1summ3, aes(y = mean), colour = 'blue')
psp1t10<-ggplot(sp1temp10, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy), xlim=c(0, 70)) +
  geom_point(data=sp1temp10, colour='turquoise') +
  geom_point(data = sp1summ10, aes(y = mean), colour = 'turquoise', size = 3) + 
  geom_line(data = sp1summ10, aes(y = mean), colour = 'turquoise')
psp1t17<-ggplot(sp1temp17, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp1temp17, colour='green') +
  geom_point(data = sp1summ17, aes(y = mean), colour = 'green', size = 3) +
  geom_line(data = sp1summ17, aes(y = mean), colour = 'green')
psp1t24<-ggplot(sp1temp24, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy), xlim=c(0, 70)) +
  geom_point(data=sp1temp24, colour='yellow') +
  geom_point(data = sp1summ24, aes(y = mean), colour = 'yellow', size = 3) +
  geom_line(data = sp1summ24, aes(y = mean), colour = 'yellow')
psp1t31<-ggplot(sp1temp31, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy), xlim=c(0, 70)) +
  geom_point(data=sp1temp31, colour='orange') +
  geom_point(data = sp1summ31, aes(y = mean), colour = 'orange', size = 3) +
  geom_line(data = sp1summ31, aes(y = mean), colour = 'orange')
psp1t38<-ggplot(sp1temp38, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy), xlim=c(0, 70)) +
  geom_point(data=sp1temp38, colour='red') +
  geom_point(data = sp1summ38, aes(y = mean), colour = 'red', size = 3) + 
  geom_line(data = sp1summ38, aes(y = mean), colour = 'red')

quartz()
grid.arrange(psp1t3, psp1t10, psp1t17, psp1t24, psp1t31, psp1t38, ncol = 3)

#set initial data for all species from mean and sd of total
#figure out what happened to May 3 readings - many high
#suspect a bunch of "outlying" high points are events of cells getting through syringe.

###sp.2
psp2t3<-ggplot(sp2temp3, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy), xlim=c(0, 70)) +
  geom_point(data=sp2temp3, colour='blue') +
  geom_point(data = sp2summ3, aes(y = mean), colour = 'blue', size = 3) + 
  geom_line(data = sp2summ3, aes(y = mean), colour = 'blue')
psp2t10<-ggplot(sp2temp10, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy), xlim=c(0, 70)) +
  geom_point(data=sp2temp10, colour='turquoise') +
  geom_point(data = sp2summ10, aes(y = mean), colour = 'turquoise', size = 3) + 
  geom_line(data = sp2summ10, aes(y = mean), colour = 'turquoise')
psp2t17<-ggplot(sp2temp17, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy), xlim=c(0, 70)) +
  geom_point(data=sp2temp17, colour='green') +
  geom_point(data = sp2summ17, aes(y = mean), colour = 'green', size = 3) +
  geom_line(data = sp2summ17, aes(y = mean), colour = 'green')
psp2t24<-ggplot(sp2temp24, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy), xlim=c(0, 70)) +
  geom_point(data=sp2temp24, colour='yellow') +
  geom_point(data = sp2summ24, aes(y = mean), colour = 'yellow', size = 3) +
  geom_line(data = sp2summ24, aes(y = mean), colour = 'yellow')
psp2t31<-ggplot(sp2temp31, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy), xlim=c(0, 70)) +
  geom_point(data=sp2temp31, colour='orange') +
  geom_point(data = sp2summ31, aes(y = mean), colour = 'orange', size = 3) +
  geom_line(data = sp2summ31, aes(y = mean), colour = 'orange')
psp2t38<-ggplot(sp2temp38, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy), xlim=c(0, 70)) +
  geom_point(data=sp2temp38, colour='red') +
  geom_point(data = sp2summ38, aes(y = mean), colour = 'red', size = 3) + 
  geom_line(data = sp2summ38, aes(y = mean), colour = 'red')

quartz()
grid.arrange(psp2t3, psp2t10, psp2t17, psp2t24, psp2t31, psp2t38, ncol = 3)



####sp. 3
psp3t3<-ggplot(sp3temp3, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp3temp3, colour='blue') +
  geom_point(data = sp3summ3, aes(y = mean), colour = 'blue', size = 3) + 
  geom_line(data = sp3summ3, aes(y = mean), colour = 'blue')
psp3t10<-ggplot(sp3temp10, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp3temp10, colour='turquoise') +
  geom_point(data = sp3summ10, aes(y = mean), colour = 'turquoise', size = 3) + 
  geom_line(data = sp3summ10, aes(y = mean), colour = 'turquoise')
psp3t17<-ggplot(sp3temp17, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp3temp17, colour='green') +
  geom_point(data = sp3summ17, aes(y = mean), colour = 'green', size = 3) +
  geom_line(data = sp3summ17, aes(y = mean), colour = 'green')
psp3t24<-ggplot(sp3temp24, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp3temp24, colour='yellow') +
  geom_point(data = sp3summ24, aes(y = mean), colour = 'yellow', size = 3) +
  geom_line(data = sp3summ24, aes(y = mean), colour = 'yellow')
psp3t31<-ggplot(sp3temp31, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp3temp31, colour='orange') +
  geom_point(data = sp3summ31, aes(y = mean), colour = 'orange', size = 3) +
  geom_line(data = sp3summ31, aes(y = mean), colour = 'orange')
psp3t38<-ggplot(sp3temp38, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp3temp38, colour='red') +
  geom_point(data = sp3summ38, aes(y = mean), colour = 'red', size = 3) + 
  geom_line(data = sp3summ38, aes(y = mean), colour = 'red')

quartz()
grid.arrange(psp3t3, psp3t10, psp3t17, psp3t24, psp3t31, psp3t38, ncol = 3)

####sp. 4
psp4t3<-ggplot(sp4temp3, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp4temp3, colour='blue') +
  geom_point(data = sp4summ3, aes(y = mean), colour = 'blue', size = 3) + 
  geom_line(data = sp4summ3, aes(y = mean), colour = 'blue')
psp4t10<-ggplot(sp4temp10, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp4temp10, colour='turquoise') +
  geom_point(data = sp4summ10, aes(y = mean), colour = 'turquoise', size = 3) + 
  geom_line(data = sp4summ10, aes(y = mean), colour = 'turquoise')
psp4t17<-ggplot(sp4temp17, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp4temp17, colour='green') +
  geom_point(data = sp4summ17, aes(y = mean), colour = 'green', size = 3) +
  geom_line(data = sp4summ17, aes(y = mean), colour = 'green')
psp4t24<-ggplot(sp4temp24, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp4temp24, colour='yellow') +
  geom_point(data = sp4summ24, aes(y = mean), colour = 'yellow', size = 3) +
  geom_line(data = sp4summ24, aes(y = mean), colour = 'yellow')
psp4t31<-ggplot(sp4temp31, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp4temp31, colour='orange') +
  geom_point(data = sp4summ31, aes(y = mean), colour = 'orange', size = 3) +
  geom_line(data = sp4summ31, aes(y = mean), colour = 'orange')
psp4t38<-ggplot(sp4temp38, aes(day, nitrate)) + coord_cartesian(ylim = c(0, maxy)) +
  geom_point(data=sp4temp38, colour='red') +
  geom_point(data = sp4summ38, aes(y = mean), colour = 'red', size = 3) + 
  geom_line(data = sp4summ38, aes(y = mean), colour = 'red')

quartz()
grid.arrange(psp4t3, psp4t10, psp4t17, psp4t24, psp4t31, psp4t38, ncol = 3)

