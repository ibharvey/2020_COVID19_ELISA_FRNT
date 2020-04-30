library(ggplot2)
library(tidyverse)

setwd("~/Box/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/")
mydata <- read.csv(paste("EC50_Correlations_Pearson_Values.csv",sep=""),header=TRUE)

rownames(mydata) <- mydata[c(1)]$X
mydata <- as.data.frame(mydata[c(2:11)])

dt2 <- mydata %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)
#tiff(filename = "EC50_Pearson_Correlations.tiff", width = 7000, height = 8000)
p <- ggplot(data = dt2, aes(x = rowname, y=colname,fill=value)) +
  geom_tile()
p
#dev.off()
ggsave("EC50_Pearson_Correlations.tiff",p)

mydata <- read.csv(paste("EC50_Correlations_PValues_Values.csv",sep=""),header=TRUE)

rownames(mydata) <- mydata[c(1)]$X
mydata <- as.data.frame(mydata[c(2:11)])

dt2 <- mydata %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)

#tiff(filename = "EC50_Pvalue_Correlations.tiff", width = 700, height = 800)
q <- ggplot(data = dt2, aes(x = rowname, y=colname,fill=value)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("green","yellow","white"), trans = "log10")
q
#dev.off()
ggsave("EC50_Pvalue_Correlations.tiff",q)
