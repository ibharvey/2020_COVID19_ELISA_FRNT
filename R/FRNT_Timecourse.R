
library(ggplot2)

# Working on FRNT timecourse plot

mydata <- read.csv(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/",
                         "FRNT_Timecourse.csv",sep=""),header=TRUE)

rownames(mydata) <- mydata[,1]
mydata <- mydata[c(2:19)]

pdata <- mydata %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname, na.rm=T)

pdata$rowname <- as.numeric(pdata$rowname)
pdata$group_num <- as.numeric(as.factor(pdata$colname))

g <- ggplot(data = pdata, aes(x=rowname, y=value, group = colname)) + 
  geom_line() +
  geom_point(aes(shape = group_num, size = 1.5)) +
  scale_shape_identity() +
  scale_y_continuous(trans='log2') +
  theme_classic() +
  ggtitle("FRNT Timecourse for COVID19+ Patients") +
  xlab("Days Post Symptom Onset") +
  ylab("Reciprocal FRNT EC50 (1/[Serum])") +
  scale_size(guide = 'none') +
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16, face = "bold"),
        axis.line = element_line(size = 1))
ggsave("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/FRNT_Timecourse.eps",g)
g