library(ggplot2)
library(tidyverse)
library(pheatmap)

setwd("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/")
mydata <- read.csv(paste("EC50_All_Mayo_Data_Patient_Days.csv",sep=""),header=TRUE)

rownames(mydata) <- mydata[c(1)]$X
mydata <- as.data.frame(mydata[c(2:13)])

dt2 <- mydata %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)
  
p<-pheatmap(mat = log2(mydata), legend_breaks = 9:15, legend_labels = paste("1:",2^(9:15),sep=""))

ggsave("EC50_All_Mayo_Data_Patient_Days.tiff", p)
ggsave("EC50_All_Mayo_Data_Patient_Days.eps",p)