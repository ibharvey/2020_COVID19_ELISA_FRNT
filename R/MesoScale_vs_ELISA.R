library(ggplot2)
library(grid)
library(gridExtra)

dilution <- "320"

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

p <- ggplot(data = )