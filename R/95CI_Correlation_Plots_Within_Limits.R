library(ggplot2)

mydata <- read.csv(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/Combined_EC50_v3.csv",sep=""),header=T)

mydata.cor <- mydata[which(mydata$SARS2.RBD.mam > 320 & mydata$SARS2.RBD.mam < 40960),]
mydata.cor.lm <- lm(SARS2.RBD.mam ~ FRNT, data = mydata.cor)
newx = seq(0,2048,by = 0.05)
conf_interval <- as.data.frame(predict(mydata.cor.lm, newdata=data.frame(FRNT=newx), interval="confidence",
                         level = 0.95))
conf_interval$x <- newx

ggplot(mydata, aes(x=FRNT, y = SARS2.RBD.mam)) +
  geom_point(size = 3) +
  geom_line(data = conf_interval, aes(x=x,y=fit)) +
  geom_line(data = conf_interval, aes(x=x,y=lwr),color="darkgray") +
  geom_line(data = conf_interval, aes(x=x,y=upr), color="darkgray") +
  theme_classic() +
  ggtitle(" ")  +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust=0.5),
        axis.text = element_text(size=16, face="bold"),
        axis.title = element_text(size=18, face="bold"),
        axis.line = element_line(size = 1.25)) +
  xlab("Reciprocal FRNT50") +
  ylab("Reciprocal EC50") +
  scale_x_continuous(trans="log2") +
  scale_y_continuous(trans="log2",breaks = c(256,256*2^2,256*2^4,256*2^6,256*2^8)) + 
  coord_cartesian(xlim = c(16,2048), ylim = c(256,65536))
 

