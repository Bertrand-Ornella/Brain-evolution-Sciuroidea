#Make Zachos time temperature
library(ggplot2) #plots

#directory
setwd("~/Desktop/Squirrel_June_8_2020/Code")

#Import squirrel data
Zachos.data<-read.csv("Final_time_temperature5.csv", header=T)

ggplot(Zachos.data, aes(x= Age, y = d18O5pt))+
  theme_classic()+
  geom_point(shape=16, size=0.5, color="turquoise3")+
  #geom_line(color="grey", size=0.1)+
  scale_y_reverse()+
  scale_x_reverse(breaks=c(70, 60, 50, 40, 30, 20, 10, 0))
