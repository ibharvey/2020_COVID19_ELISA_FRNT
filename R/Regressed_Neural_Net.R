library(tidyverse)
library(caret)
library(glmnet)

setwd("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/")


lambda <- 10^seq(-3,3,length=100)
antigens <- c("SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer","SARS2.NP.RNA","SARS2.NP.FL")
output <- data.frame(matrix(ncol = 8, nrow = 0))
names(output) <- c("Dilution","Combination","CoIndex","Cutoff","Accuracy","Precision","Recall","F1")
for (f in c(320,1280))
{
  # Load in the MayoClinic Dataset
  mydata <- read.csv(paste(f,".csv",sep=""),header=F)
  mdf <- data.frame(mydata)
  heading <- c("FRNT","Chk.P62E1","BSA","HA","MERS.bRBD","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer","SARS2.ORF8","SARS2.NP.RNA","SARS2.NP.FL","SARS2.ORF7a")
  # First dilution
  colnames(mdf) <- heading
  
  # Load in the WashU Dataset
  myconv <- read.csv(paste("../../AliConvalescent/",f,".csv",sep=""),header=T)
  cdf <- data.frame(myconv)
  
  # For each n in (6 choose n) 
  for (co in 2)
  {
    a_set <- combn(1:6,co)
    # For each combination of antigens possible
    for (i in 1:length(a_set[c(1),]))
    {
      lambda <- 10^seq(-3,3, length = 100)
      # Extract the relevant antigens for this nnet
      this_antigen_set <- c("FRNT",antigens[a_set[,i]])
      mydata.whole <- rbind(mdf[this_antigen_set], cdf[this_antigen_set])
      for (j in 40:1000)
      {
        set.seed(123)
        mydata.cutoff <- mydata.whole
        mydata.cutoff$FRNT <- ifelse(mydata.cutoff$FRNT > j, 1, 0)
        mydata.cutoff$FRNT <- as.factor(mydata.cutoff$FRNT)
        train.samples <- mydata.cutoff$FRNT %>%
          createDataPartition(p = 0.6, list=F)
        train.data <- mydata.cutoff[train.samples,]
        test.data <- mydata.cutoff[-train.samples,]
        lasso <- train(
          FRNT~., data = train.data, method = "glmnet",
          trControl = trainControl("cv", number = 10),
          tuneGrid = expand.grid(alpha = 1, lambda = lambda)
        )
        #out1 <- coef(lasso$finalModel, lasso$bestTune$lambda)
        predictions <- lasso %>% predict(test.data)
        tbl <- as.matrix(table(test.data$FRNT,predictions))
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
        print(paste(f,j, toString(this_antigen_set)))
        print(out)
        print('----------------------------------')
      }
      
      
    }
  }
}