library(ggplot2)
library(tidyverse)
library(pheatmap)

setwd("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/")
mydata <- read.csv(paste("EC50_All_Mayo_Data_Patient_Days.csv",sep=""),header=TRUE)

rownames(mydata) <- mydata[c(1)]$X
mydata <- as.data.frame(mydata[c("Chk.P62E1","BSA","HA","MERS.bRBD","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD",
                                 "SARS2.Trimer","SARS2.ORF8","SARS2.NP.FL","SARS2.ORF7a")])

colnames(mydata) <- c("CHIKV p62-E1","BSA","Influenza HA","MERS-CoV bRBD","SARS-CoV bRBD","SARS-CoV-2 bRBD",
                      "SARS-CoV-2 mRBD","SARS-CoV-2 S","SARS-CoV-2 ORF8","SARS-CoV-2 NP","SARS-CoV-2 ORF7a")

dt2 <- mydata %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)
  
p<-pheatmap(mat = log2(mydata), legend_breaks = 9:15, legend_labels = paste("1:",2^(9:15),sep=""), 
            fontsize_row = 9, fontsize_col = 14)

ggsave("EC50_All_Mayo_Data_Patient_Days_v2.tiff", p)
ggsave("EC50_All_Mayo_Data_Patient_Days_v2.eps",p)