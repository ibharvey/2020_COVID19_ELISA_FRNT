library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(corrplot)
library(abind)
library(pheatmap)
setwd("/home/ian/Dropbox/Fremont/Bioinformatics/ELISA/Prelim_Screen/")

temp1 <- read.csv("big_antigen_screen_1.csv",header=T)[2:37]
temp2 <- read.csv("big_antigen_screen_2.csv",header=T)[2:37]
temp3 <- read.csv("big_antigen_screen_3.csv",header=T)[2:37]
temp4 <- read.csv("big_antigen_screen_4.csv",header=T)[2:37]

# for(i in 2:length(temp1))
# {
#   g<- ggplot() +
#     geom_line(data = temp1, aes(x = X, y = temp1[,i]), linetype = "dotted", color = "gray") +
#     geom_line(data = temp2, aes(x = X, y = temp2[,i]), linetype = "dotdash", color = "gray") +
#     geom_line(data = temp3, aes(x = X, y = temp3[,i]), linetype = "dashed") +
#     geom_line(data = temp4, aes(x = X, y = temp4[,i]), linetype = "twodash") +
#     theme_classic() +
#     scale_x_continuous(trans = "log2") +
#     ylim(0,3)
#   g
# }

high_conc <- rbind(cbind(temp1[1:12]/temp1$Influenza.HA[1], 
                         temp1[13:24]/temp1$Influenza.HA.1[1], 
                         temp1[25:36]/temp1$Influenza.HA.2[1])[1,],
                   cbind(temp2[1:12]/temp2$Influenza.HA[1], 
                         temp2[13:24]/temp2$Influenza.HA.1[1],
                         temp2[25:36]/temp2$Influenza.HA.2[1])[1,],
                   cbind(temp3[1:12]/temp3$Influenza.HA[1], 
                         temp3[13:24]/temp3$Influenza.HA.1[1],
                         temp3[25:36]/temp3$Influenza.HA.2[1])[1,],
                   cbind(temp4[1:12]/temp4$Influenza.HA[1], 
                         temp4[13:24]/temp4$Influenza.HA.1[1],
                         temp4[25:36]/temp4$Influenza.HA.2[1])[1,]
)

mydata.screen <- high_conc[c("Influenza.HA", "SARS2.bRBD","SARS2.mRBD","SARS2.Trimer",
                             "NSP1","NSP3","NSP7.8","NSP8","NSP9","NSP10","NSP15","NSP16","ORF8.pure",
                             "SARS1.bRBD","MERS.bRBD", "ORF7a","NP.RNA","NP.FL")]

pheatmap(mydata.screen, labels_row = c("Negative 1","Negative 2","Positive 1","Positive 2"),
         labels_col = c("Influenza HA","SARS2 bRBD","SARS2 mRBD","SARS2 Spike Trimer","SARS2 NSP1",
                        "SARS2 NSP3","SARS2 NSP7+8","SARS2 NSP8","SARS2 NSP9","SARS2 NSP10","SARS2 NSP15",
                        "SARS2 NSP16","SARS2 ORF8","SARS1 bRBD","MERS bRBD", "SARS2 ORF7a","SARS2 NP-RNA",
                        "SARS2 NP-FL"),
         fontsize = 18)
