library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(corrplot)
library(scales)

add_PRA <- function(tbl)
{
  #tbl[[1]] <- tbl[[1]] + 16 # If assuming common CoV serum does not neutralize SARS-CoV-2
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
  out <- data.frame( Cutoff = it,
                     Accuracy =accuracy,
                     PPV = precision[2],
                     NPV = precision[1],
                     Sensitivity =recall[2],
                     Specificity =recall[1],
                     #F1F = f1[1],
                     F1 = f1[2]
  )
  out
}

dilution <- "1280"

######### Negs
mynegs <- read.csv(paste("~/Dropbox/WashU/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/",
                         "Chicago_duplicates_", dilution, ".csv",sep=""),header=TRUE)
rownames(mynegs) <- c("KL7","AM21","SW9","PZ24","JS33","AH48","CR35","RW26","VK50","CG8","BL67","BG20","JG64","JS19","CK34","JT6")
mesnegs <- read.csv(paste("~/Dropbox/WashU/Fremont/Bioinformatics/Code/Github/2020_COVID19_ELISA_FRNT/data/",
                          dilution, "_all_UChicago.csv",sep=""),header=TRUE)
rownames(mesnegs) <- mesnegs$X
mesnegs <- mesnegs[c(2:11)]
names(mesnegs) <- c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA","SARS2.NP.RNA")
mesnegs <- mesnegs[rownames(mynegs),]
mynegs <- mynegs[names(mesnegs)]


######### MayoClinic
mydata <- read.csv(paste("~/Dropbox/WashU/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/",
                          dilution, ".csv",sep=""),header=F)
names(mydata) <- c("FRNT","ChkP62E1","BSA","HA","MERS.bRBD","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer","SARS2.ORF8","SARS2.NP.RNA","SARS2.NP.FL","SARS2.ORF7a")
#mydata <- mydata[c("FRNT","ChkP62E1","HA","MERS.bRBD","SARS2.ORF8","SARS2.NP.RNA","SARS2.NP.FL","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer")]

mydata2 <- read.csv(paste("~/Dropbox/WashU/Fremont/Bioinformatics/ELISA/AliConvalescent/", dilution, ".csv",sep=""),header=T)

mydata_cat <- bind_rows(mydata,mydata2)
mydata_cat <- mydata_cat[c("FRNT","ChkP62E1","HA","MERS.bRBD","SARS2.ORF8","SARS2.NP.RNA","SARS2.NP.FL","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer", "SARS1.Trimer","MERS.Trimer","SARS1.mRBD","SARS2.ORF7a","BSA")]


mesdata <- read.csv(paste("~/Dropbox/WashU/Fremont/Bioinformatics/Code/Github/2020_COVID19_ELISA_FRNT/data/",
                          dilution, "_all_MayoClinic.csv",sep=""),header=TRUE)
rownames(mesdata) <- mesdata$X
mesdata <- mesdata[c(2:11)]
names(mesdata) <- c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA","SARS2.NP.RNA")


mesdata2 <- read.csv(paste("~/Dropbox/WashU/Fremont/Bioinformatics/Code/Github/2020_COVID19_ELISA_FRNT/data/",
                           dilution, "_all_Convalescent.csv",sep=""),header=TRUE)
mesdata2 <- mesdata2[c(2:11)]
names(mesdata2) <- c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA","SARS2.NP.RNA")


############################ Accuracy of Single-point ELISAs ###########################


ggplot(data = mydata, aes(x = FRNT, y=SARS1.bRBD)) +
  geom_point(na.rm=T) +
  scale_x_continuous(trans="log2") +
  theme_classic() +
  ggtitle(paste("SARS1.bRBD vs FRNT50 [1:", dilution, "]", sep = ""))  +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab("FRNT50") + ylab("Abs 450") + 
  geom_abline(intercept = 3.2, slope=0,linetype="solid",color="dark gray",size=2) + 
  geom_abline(intercept = 1, slope=0,linetype="dotted",color="dark gray",size=2)


############################ Add FRNT info to Mesoscale data ############################


mesdata2$FRNT <- mydata2$FRNT
mesdata$FRNT <- mydata$FRNT

mesdata_cat <- bind_rows(mesdata, mesdata2)


########## Mesoscale

qplot <- ggplot(data = mesdata_cat, aes(x = FRNT, y=SARS2.Trimer)) +
  scale_colour_gradient(trans='log2', limits=c(32,4096), na.value = 'black', name='FRNT') +
  geom_point(na.rm=T, size = 2) + 
  # S
  geom_hline(yintercept = 1e7, linetype="solid",color="dark gray",size=2) +
  geom_hline(yintercept = 0.5e6, linetype="dotted",color="dark gray",size=2)+
  # # NP
  # geom_hline(yintercept = 1e7, linetype="solid",color="dark gray",size=2) +
  # geom_hline(yintercept = 0.5e6, linetype="dotted",color="dark gray",size=2)+
  # # mRBD opt
  # geom_hline(yintercept = 3e6, linetype="solid",color="dark gray",size=2) +
  # geom_hline(yintercept = 7e4, linetype="dotted",color="dark gray",size=2)+
  # # S2 bRBD
  # geom_hline(yintercept = 1e6, linetype="solid",color="dark gray",size=2) +
  # geom_hline(yintercept = 1e5, linetype="dotted",color="dark gray",size=2)+
  # # S1 bRBD and ORF8
  # geom_hline(yintercept = 2.5e5, linetype="solid",color="dark gray",size=2) +
  # geom_hline(yintercept = 5e4, linetype="dotted",color="dark gray",size=2)+
  # # Influenza H.A
  # geom_hline(yintercept = 9e6, linetype="solid",color="dark gray",size=2) +
  # geom_hline(yintercept = 4e6, linetype="dotted",color="dark gray",size=2)+
  scale_x_continuous(trans="log2") +
  scale_y_continuous(trans="log2") +
  theme_classic() +
  ggtitle("SARS-CoV-2 S",
    subtitle=paste("Mesoscale [1:", dilution, "]", sep = ""))  +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab("FRNT50") + ylab("Electrochemoluminescence") 

# Patient 2
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mesdata_cat[c(11,20),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), color = "dark red", size = 2)
#qplot <- qplot + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 2"), nudge_x = 0.65, nudge_y = -0.1)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mesdata_cat[c(20,21),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), data = df3, color = "dark red", size = 2)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mesdata_cat[c(21,22),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), data = df3, color = "dark red", size = 2, arrow = arrow(length = unit(0.2, "inches")))


#Patient 8
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mesdata_cat[c(6,31),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), color = "dark green", size = 2)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mesdata_cat[c(31,32),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), color = "dark green", size = 2, arrow = arrow(length = unit(0.2, "inches")))

#Patient 3
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mesdata_cat[c(13,26),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), color = "dark orange", size = 2)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mesdata_cat[c(26,39),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), color = "dark orange", size = 2, arrow = arrow(length = unit(0.2, "inches")))

qplot
ggsave(paste("~/Dropbox/WashU/Fremont/Bioinformatics/ELISA/Meso_320dil_FRNT_Trimer_122_datapoints.tiff",sep = ""), 
       plot = qplot,
       width = 23.230/1.5,
       height = 14.949/1.5,
       units = "cm")

results.weak <- data.frame(matrix(ncol = 8, nrow = 0))
results.strong <- data.frame(matrix(ncol = 8, nrow = 0))
for(it in 20:1000)
{
  results.weak <- rbind(results.weak,
                         add_PRA(as.matrix(table(mesdata_cat$FRNT > it,mesdata_cat$SARS2.Trimer>0.5e6))))
  results.strong <- rbind(results.strong,
                           add_PRA(as.matrix(table(mesdata_cat$FRNT > it,mesdata_cat$SARS2.Trimer>1e7))))
}

#Accuracy plot based on cutoff
plot_out3 <- ggplot() + 
  scale_size(guide = 'none') +
  geom_line(data = results.weak, aes(x = Cutoff, y = Accuracy), linetype = "dotted", size = 1) +
  geom_line(data = results.strong, aes(x = Cutoff, y = Accuracy), linetype = "solid", size = 1) +
  theme_classic() +
  ggtitle("Prediction Accuracy of Influenza SARS2.Trimer",
    subtitle=paste("Mesoscale [1:",dilution,"]",sep="")) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab("Lower Limit Cutoff 1/[FRNT50]") +
  scale_y_continuous(labels=function(x) paste0(x*100,"%"), limits = c(0,1)) +
  scale_x_continuous(trans="log2") +
  coord_cartesian(xlim = c(20,1000), ylim = c(.5,1))

plot_out3
ggsave(paste("~/Dropbox/WashU/Fremont/Bioinformatics/ELISA/Accuracy_Meso_320dil_FRNT_Trimer_122_datapoints.tiff",sep = ""), 
       plot = plot_out3,
       width = 23.230/1.5,
       height = 14.949/1.5,
       units = "cm")

# ELISA

qplot <- ggplot(data = mydata_cat, aes(x = FRNT, y=SARS2.Trimer)) +
  geom_point(na.rm=T, size = 2) + 
  # S and NP
  geom_hline(yintercept = 3.2, linetype="solid",color="dark gray",size=2) +
  geom_hline(yintercept = 0.5, linetype="dotted",color="dark gray",size=2)+
  # # mRBD opt
  # geom_hline(yintercept = 2.0, linetype="solid",color="dark gray",size=2) +
  # geom_hline(yintercept = 0.2, linetype="dotted",color="dark gray",size=2)+
  # # S2 bRBD
  # geom_hline(yintercept = 2.0, linetype="solid",color="dark gray",size=2) +
  # geom_hline(yintercept = 0.5, linetype="dotted",color="dark gray",size=2)+
  # # S1 bRBD and ORF8
  # geom_hline(yintercept = 1.0, linetype="solid",color="dark gray",size=2) +
  # geom_hline(yintercept = 0.3, linetype="dotted",color="dark gray",size=2)+
  # # Influenza H.A
  # geom_hline(yintercept = 3.4, linetype="solid",color="dark gray",size=2) +
  # geom_hline(yintercept = 2.7, linetype="dotted",color="dark gray",size=2)+
  scale_x_continuous(trans="log2") +
  scale_y_continuous(trans="log2") +
  theme_classic() +
  ggtitle("SARS-CoV-2 S",
    subtitle=paste("ELISA [1:", dilution, "]", sep = ""))  +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab("FRNT50") + ylab("Absorbance (450nm)") 

# Patient 2
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata_cat[c(11,20),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), color = "dark red", size = 2)
#qplot <- qplot + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 2"), nudge_x = 0.65, nudge_y = -0.1)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata_cat[c(20,21),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), data = df3, color = "dark red", size = 2)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata_cat[c(21,22),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), data = df3, color = "dark red", size = 2, arrow = arrow(length = unit(0.2, "inches")))


#Patient 8
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata_cat[c(6,31),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), color = "dark green", size = 2)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata_cat[c(31,32),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), color = "dark green", size = 2, arrow = arrow(length = unit(0.2, "inches")))

#Patient 3
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata_cat[c(13,26),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), color = "dark orange", size = 2)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata_cat[c(26,39),])))[c("FRNT","SARS2.Trimer")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.Trimer1, xend = FRNT2, yend = SARS2.Trimer2), color = "dark orange", size = 2, arrow = arrow(length = unit(0.2, "inches")))

qplot


ggsave(paste("~/Dropbox/WashU/Fremont/Bioinformatics/ELISA/ELISA_320dil_FRNT_Trimer_122_datapoints.tiff",sep = ""), 
       plot = qplot,
       width = 23.230/1.5,
       height = 14.949/1.5,
       units = "cm")

results.weak <- data.frame(matrix(ncol = 8, nrow = 0))
results.strong <- data.frame(matrix(ncol = 8, nrow = 0))
for(it in 20:1000)
{
  results.weak <- rbind(results.weak,
                         add_PRA(as.matrix(table(mydata_cat$FRNT > it,mydata_cat$SARS2.Trimer>0.5))))
  results.strong <- rbind(results.strong,
                           add_PRA(as.matrix(table(mydata_cat$FRNT > it,mydata_cat$SARS2.Trimer>3.2))))
}

#Accuracy plot based on cutoff
plot_out3 <- ggplot() + 
  scale_size(guide = 'none') +
  geom_line(data = results.weak, aes(x = Cutoff, y = Accuracy), linetype = "dotted", size = 1) +
  geom_line(data = results.strong, aes(x = Cutoff, y = Accuracy), linetype = "solid", size = 1) +
  theme_classic() +
  ggtitle("Prediction Accuracy of SARS-CoV-2 S",
    subtitle=paste("ELISA [1:",dilution,"]",sep="")) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab("Lower Limit Cutoff 1/[FRNT50]") +
  scale_y_continuous(labels=function(x) paste0(x*100,"%"), limits = c(0,1)) +
  scale_x_continuous(trans="log2") +
  coord_cartesian(xlim = c(20,1000), ylim = c(.5,1))

plot_out3
ggsave(paste("~/Dropbox/WashU/Fremont/Bioinformatics/ELISA/Accuracy_ELISA_320dil_FRNT_Trimer_122_datapoints.tiff",sep = ""), 
       plot = plot_out3,
       width = 23.230/1.5,
       height = 14.949/1.5,
       units = "cm")

