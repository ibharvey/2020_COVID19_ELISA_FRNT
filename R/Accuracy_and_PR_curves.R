library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(caret)
library(tidyverse)

dilution <- "320"

######### Negs
mynegs <- read.csv(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/",
                         "Chicago_duplicates_", dilution, ".csv",sep=""),header=TRUE)
rownames(mynegs) <- c("KL7","AM21","SW9","PZ24","JS33","AH48","CR35","RW26","VK50","CG8","BL67","BG20","JG64","JS19","CK34","JT6")
mesnegs <- read.csv(paste("/home/ian/Dropbox/Fremont/Bioinformatics/Code/Github/2020_COVID19_ELISA_FRNT/data/",
                          dilution, "_all_UChicago.csv",sep=""),header=TRUE)
rownames(mesnegs) <- mesnegs$X
mesnegs <- mesnegs[c(2:11)]
names(mesnegs) <- c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA","SARS2.NP.RNA")
mesnegs <- mesnegs[rownames(mynegs),]
mynegs <- mynegs[names(mesnegs)]


######### MayoClinic
mydata <- read.csv(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/",
                         dilution, ".csv",sep=""),header=F)
names(mydata) <- c("FRNT","ChkP62E1","BSA","HA","MERS.bRBD","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer","SARS2.ORF8","SARS2.NP.RNA","SARS2.NP.FL","SARS2.ORF7a")
#mydata <- mydata[c("FRNT","ChkP62E1","HA","MERS.bRBD","SARS2.ORF8","SARS2.NP.RNA","SARS2.NP.FL","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer")]

mydata2 <- read.csv(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/AliConvalescent/", dilution, ".csv",sep=""),header=T)

mydata_cat <- bind_rows(mydata,mydata2)
mydata_cat <- mydata_cat[c("FRNT","ChkP62E1","HA","MERS.bRBD","SARS2.ORF8","SARS2.NP.RNA","SARS2.NP.FL","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer", "SARS1.Trimer","MERS.Trimer","SARS1.mRBD","SARS2.ORF7a","BSA")]


mesdata <- read.csv(paste("/home/ian/Dropbox/Fremont/Bioinformatics/Code/Github/2020_COVID19_ELISA_FRNT/data/",
                          dilution, "_all_MayoClinic.csv",sep=""),header=TRUE)
rownames(mesdata) <- mesdata$X
mesdata <- mesdata[c(2:11)]
names(mesdata) <- c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA","SARS2.NP.RNA")


mesdata2 <- read.csv(paste("/home/ian/Dropbox/Fremont/Bioinformatics/Code/Github/2020_COVID19_ELISA_FRNT/data/",
                           dilution, "_all_Convalescent.csv",sep=""),header=TRUE)
mesdata2 <- mesdata2[c(2:11)]
names(mesdata2) <- c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA","SARS2.NP.RNA")


mesdata2$FRNT <- mydata2$FRNT
mesdata$FRNT <- mydata$FRNT

mesdata <- mesdata[c("FRNT","ChkP62E1","HA","MERS.bRBD","SARS2.ORF8","SARS2.NP.RNA","SARS2.NP.FL",
                     "SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer")]

mesdata_cat <- bind_rows(mesdata, mesdata2)

set.seed(123)
output <- data.frame(matrix(ncol = 8, nrow = 0))
for(j in 40:1000)
{
  for(k in 1:100)
  {
    mydata.cutoff <- mydata_cat
    mydata.cutoff$FRNT <- ifelse(mydata.cutoff$FRNT > j, 1, 0)
    mydata.cutoff$FRNT <- as.factor(mydata.cutoff$FRNT)
    train.samples <- mydata.cutoff$FRNT %>%
      createDataPartition(p = 0.6, list=F)
    train.data <- mydata.cutoff[train.samples,]
    test.data <- mydata.cutoff[-train.samples,]
    modelFit <- pcaNNet(x=train.data[,c("SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer")], y=train.data$FRNT, size = 1)
    test.predict <- predict(modelFit, test.data[,c("SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer")])
    tbl <- as.matrix(table(test.data$FRNT,round(test.predict[,2])))
    n = sum(tbl) # number of instances
    nc = nrow(tbl) # number of classes
    diag = diag(tbl) # number of correctly classified instances per class 
    rowsums = apply(tbl, 1, sum) # number of instances per class
    colsums = apply(tbl, 2, sum) # number of predictions per class
    p = rowsums / n # distribution of instances over the actual classes
    q = colsums / n # distribution of instances over the predicted classes
    accuracy <- sum(diag) / n
    precision <- diag / colsums
    recall <- diag / rowsums
    f1 = 2 * precision * recall / (precision + recall)
    out <- data.frame( Dilution = f,
                       Combination = co,
                       CoIndex = i,
                       Cutoff = j,
                       Accuracy = mean(accuracy,na.rm=T),
                       Precision = mean(precision,na.rm=T),
                       Recall = mean(recall,na.rm=T),
                       F1 = mean(f1,na.rm=T)
    )
    output <- rbind(output, out)
  }
  print(j)
}

gagg <- aggregate( Accuracy ~ Cutoff, output, mean )

g <- ggplot() + 
  geom_point(data = gagg, aes(x = Cutoff, y = Accuracy), color = "green") +
  ylim(0,1) + xlim(40,1000)
g