library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(corrplot)

dilution <- "320"
dilution.mes <- "1280"

######### Negs
mynegs <- read.csv(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/",
                         "Chicago_duplicates_", dilution, ".csv",sep=""),header=TRUE)
mynegs <- mynegs[c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA")]

rownames(mynegs) <- c("KL7","AM21","SW9","PZ24","JS33","AH48","CR35","RW26","VK50","CG8","BL67","BG20","JG64","JS19","CK34","JT6")
mesnegs <- read.csv(paste("~/Dropbox/Fremont/Bioinformatics/Code/Github/2020_COVID19_ELISA_FRNT/data/",
                          dilution.mes, "_all_UChicago.csv",sep=""),header=TRUE)
rownames(mesnegs) <- mesnegs$X
mesnegs <- mesnegs[c(2:11)]
names(mesnegs) <- c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA","SARS2.NP.RNA")
# Reorder rows to match mynegs dataframe
mesnegs <- mesnegs[rownames(mynegs),]
# Reorder cols to match mynegs dataframe
mesnegs <- mesnegs[names(mynegs)]


######### MayoClinic
mydata <- read.csv(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/",
                          dilution, ".csv",sep=""),header=F)
names(mydata) <- c("FRNT","ChkP62E1","BSA","HA","MERS.bRBD","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer","SARS2.ORF8","SARS2.NP.RNA","SARS2.NP.FL","SARS2.ORF7a")
# Reorder cols to match mynegs
mydata <- mydata[c("FRNT","SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA")]

mydata2 <- read.csv(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/AliConvalescent/", dilution, ".csv",sep=""),header=T)


mydata_cat <- bind_rows(mydata,mydata2)
mydata_cat <- mydata_cat[c("FRNT","ChkP62E1","HA","MERS.bRBD","SARS2.ORF8","SARS2.NP.RNA","SARS2.NP.FL","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer", "SARS1.Trimer","MERS.Trimer","SARS1.mRBD","SARS2.ORF7a","BSA")]


mesdata <- read.csv(paste("~/Dropbox/Fremont/Bioinformatics/Code/Github/2020_COVID19_ELISA_FRNT/data/",
                          dilution.mes, "_all_MayoClinic.csv",sep=""),header=TRUE)
rownames(mesdata) <- mesdata$X
mesdata <- mesdata[c(2:11)]
names(mesdata) <- c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA","SARS2.NP.RNA")


mesdata2 <- read.csv(paste("~/Dropbox/Fremont/Bioinformatics/Code/Github/2020_COVID19_ELISA_FRNT/data/",
                           dilution.mes, "_all_Convalescent.csv",sep=""),header=TRUE)
mesdata2 <- mesdata2[c(2:11)]
names(mesdata2) <- c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA","SARS2.NP.RNA")

############################ Box Plot of Mesoscale data #################################

box_cols <- c("ChkP62E1", "HA","MERS.bRBD","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD",
              "SARS2.Trimer","SARS2.NP.FL","SARS2.ORF8")
box_data <- mesdata[,box_cols]
colnames(box_data) <- c("CHIKV p62-E1","Influenza HA","MERS-CoV bRBD","SARS-CoV bRBD",
                        "SARS-CoV-2 bRBD","SARS-CoV-2 mRBD","SARS-CoV-2 S",
                        "SARS-CoV-2 NP","SARS-CoV-2 ORF8")
box_data2 <- mesdata2[,box_cols]
colnames(box_data2) <- c("CHIKV p62-E1","Influenza HA","MERS-CoV bRBD","SARS-CoV bRBD",
                        "SARS-CoV-2 bRBD","SARS-CoV-2 mRBD","SARS-CoV-2 S",
                        "SARS-CoV-2 NP","SARS-CoV-2 ORF8")
box_negs <- mesnegs[,box_cols]
colnames(box_negs) <- c("CHIKV p62-E1","Influenza HA","MERS-CoV bRBD","SARS-CoV bRBD",
                        "SARS-CoV-2 bRBD","SARS-CoV-2 mRBD","SARS-CoV-2 S",
                        "SARS-CoV-2 NP","SARS-CoV-2 ORF8")

dbox <- box_data %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname)

dbox <- cbind(dbox, "SARS2-CoV-2 #1")
names(dbox)[4] <- "sample"

nbox <- box_negs %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname)

nbox <- cbind(nbox, "Common CoV")
names(nbox)[4] <- "sample"

cbox <- box_data2 %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname)

cbox <- cbind(cbox, "SARS2-CoV-2 #2")
names(cbox)[4] <- "sample"

box_results <- rbind(dbox,cbox,nbox)


q <- ggplot(data = box_results, aes(x = colname, y=value, fill=sample)) +
  scale_y_continuous(trans='log10') +
  geom_boxplot(outlier.size = 2.0) + 
  ggtitle(paste("Mesoscale 1:",dilution," Seroreactivity",sep="")) +
  theme_classic()+
  theme(plot.title = element_text(size = 16,hjust=0.5),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle=45,hjust=1)
        ) + 
  xlab("") + 
  ylab("Electrochemiluminescence")
q
ggsave(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/Mesoscale_Box_Plot_v3_",dilution,".eps",sep=""),q)


############################ Box Plot of ELISA data #################################

box_cols <- c("ChkP62E1", "HA","MERS.bRBD","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer","SARS2.NP.FL","SARS2.ORF8")
box_data <- mydata[,box_cols]
colnames(box_data) <- c("CHIKV p62-E1","Influenza HA","MERS-CoV bRBD","SARS-CoV bRBD",
                        "SARS-CoV-2 bRBD","SARS-CoV-2 mRBD","SARS-CoV-2 S",
                        "SARS-CoV-2 NP","SARS-CoV-2 ORF8")
box_data2 <- mydata2[,c("ChkP62E1", "HA","SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer","SARS2.NP.FL","SARS2.ORF8")]
colnames(box_data2) <- c("CHIKV p62-E1","Influenza HA","SARS-CoV bRBD",
                        "SARS-CoV-2 bRBD","SARS-CoV-2 mRBD","SARS-CoV-2 S",
                        "SARS-CoV-2 NP","SARS-CoV-2 ORF8")
box_negs <- mynegs[,box_cols]
colnames(box_negs) <- c("CHIKV p62-E1","Influenza HA","MERS-CoV bRBD","SARS-CoV bRBD",
                        "SARS-CoV-2 bRBD","SARS-CoV-2 mRBD","SARS-CoV-2 S",
                        "SARS-CoV-2 NP","SARS-CoV-2 ORF8")

  
dbox <- box_data %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname)

dbox <- cbind(dbox, "SARS2-CoV-2 #1")
names(dbox)[4] <- "sample"

nbox <- box_negs %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname)

nbox <- cbind(nbox, "Common CoV")
names(nbox)[4] <- "sample"

cbox <- box_data2 %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname)

cbox <- cbind(cbox, "SARS2-CoV-2 #2")
names(cbox)[4] <- "sample"

box_results <- rbind(dbox,cbox,nbox)
# Converting to 10^x and then scaling back down so the formatting of ELISA and Mesoscale plots are equivalent
box_results$value <- 10^box_results$value

q <- ggplot(data = box_results, aes(x = colname, y=value, fill=sample)) +
  scale_y_continuous(trans='log10') +
  geom_boxplot(outlier.size = 2.0) + 
  ggtitle(paste("ELISA 1:",dilution," Seroreactivity",sep="")) +
  theme_classic()+
  theme(plot.title = element_text(size = 16,hjust=0.5),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle=45,hjust=1)
  ) + 
  xlab("") + 
  ylab("Absorbance (450nm)")
q
ggsave(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/ELISA_Box_Plot_v3_",dilution,".eps",sep=""),q)

############################ Accuracy of Single-point ELISAs ###########################


ggplot(data = mydata, aes(x = FRNT, y=SARS2.mRBD)) +
  geom_point(na.rm=T) +
  scale_x_continuous(trans="log2") +
  theme_classic() +
  ggtitle(paste("SARS2.mRBDr vs FRNT50 [1:", dilution, "]", sep = ""))  +
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

mesdata <- mesdata[c("FRNT","ChkP62E1","HA","MERS.bRBD","SARS2.ORF8","SARS2.NP.RNA","SARS2.NP.FL",
                     "SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer")]

mesdata_cat <- bind_rows(mesdata, mesdata2)

mesdata_catcor <- mesdata_cat[c("FRNT","ChkP62E1","HA","SARS2.ORF8","SARS2.NP.FL",
                                "SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer")]

mesdata_catcor$FRNT <- log2(mesdata_catcor$FRNT)

mydata_catcor <- mydata_cat[c("FRNT","ChkP62E1","HA","SARS2.ORF8","SARS2.NP.FL",
                              "SARS1.bRBD","SARS2.bRBD","SARS2.mRBD","SARS2.Trimer")]
colnames(mesdata_catcor) <- c("log2(FRNT50)","CHIKV p62-E1","Influenza HA","SARS-CoV-2 ORF8",
                              "SARS-CoV-2 NP","SARS-CoV bRBD","SARS-CoV-2 bRBD","SARS-CoV-2 mRBD","SARS-CoV-2 S")

mydata_catcor$FRNT <- log2(mydata_catcor$FRNT)
colnames(mydata_catcor) <- c("log2(FRNT50)","CHIKV p62-E1","Influenza HA","SARS-CoV-2 ORF8",
                             "SARS-CoV-2 NP","SARS-CoV bRBD","SARS-CoV-2 bRBD","SARS-CoV-2 mRBD","SARS-CoV-2 S")

res <- cor(mydata_catcor,y=mesdata_catcor)

tiff(filename="~/Dropbox/Fremont/Bioinformatics/ELISA/Corrplot_Meso_ELISA.tiff")
#ELISA on rows, Mesoscale on columns
corrplot(res,tl.col = "black", method = "color", 
         addCoef.col = "black",
         tl.cex = 1.6, cl.cex = 1.4, cl.align.text = "l",
         number.cex = 1)
dev.off()

##################################################################

######################################################

mes_set <- c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA","SARS2.NP.RNA")

cdf <- mesdata_cat[c("SARS2.NP.FL", "SARS2.Trimer", "SARS2.mRBD", "SARS2.ORF8")]
ddf <- mesdata_cat[c("SARS2.NP.FL", "SARS2.Trimer", "SARS2.mRBD", "SARS2.ORF8")]
ndf <- mesnegs[c("SARS2.NP.FL", "SARS2.Trimer", "SARS2.mRBD", "SARS2.ORF8")]
# 

# cdf <- mesdata2[c("SARS1.bRBD","SARS2.mRBD")]
# ddf <- mesdata[c("SARS1.bRBD","SARS2.mRBD")]
# ndf <- mesnegs[c("SARS1.bRBD","SARS2.mRBD")]

pca_df <- prcomp(cdf, retx=TRUE, center=TRUE, scale=TRUE)
nega_points <- as.data.frame(predict(pca_df, ndf))
mayo_points <- as.data.frame(predict(pca_df, ddf))
mayo_points$FRNT <- mesdata_cat$FRNT
eigs <- pca_df$sdev^2
provar <- eigs / sum(eigs)

my_points <- as.data.frame(predict(pca_df, cdf))
my_points$FRNT <- mesdata_cat$FRNT

plot_out2 <- ggplot()
plot_out2 <- plot_out2 + scale_size(guide = 'none')

# Diagonal cutoffs front-slash
# # S & NP
# plot_out2 <- plot_out2 + geom_abline(intercept = -1.9, slope=-1,linetype="dotted",color="dark gray",size=2)
# plot_out2 <- plot_out2 + geom_abline(intercept = -1.4, slope=-1,linetype="dashed",color="dark gray", size=2)
# plot_out2 <- plot_out2 + geom_abline(intercept = 1, slope=-1,linetype="solid",color="dark gray", size=2)
# # S, mRBD, NP
# plot_out2 <- plot_out2 + geom_abline(intercept = -2.8, slope=-1,linetype="dotted",color="dark gray",size=2)
# plot_out2 <- plot_out2 + geom_abline(intercept = -2, slope=-1,linetype="dashed",color="dark gray", size=2)
# plot_out2 <- plot_out2 + geom_abline(intercept = 0.6, slope=-1,linetype="solid",color="dark gray", size=2)
# S, mRBD, NP, ORF8
plot_out2 <- plot_out2 + geom_abline(intercept = -2.4, slope=-1,linetype="dotted",color="dark gray",size=2)
plot_out2 <- plot_out2 + geom_abline(intercept = -1.6, slope=-1,linetype="dashed",color="dark gray", size=2)
plot_out2 <- plot_out2 + geom_abline(intercept = 0.6, slope=-1,linetype="solid",color="dark gray", size=2)
plot_out2 <- plot_out2 + geom_point(data = my_points, aes(x=PC1, y=PC2, size=3, colour=FRNT))
plot_out2 <- plot_out2 + geom_point(data = mayo_points, aes(x=PC1, y=PC2, size=3, colour=FRNT))
plot_out2 <- plot_out2 + geom_point(data = nega_points, aes(x=PC1, y=PC2, size=1.5), color="red")
plot_out2 <- plot_out2 + geom_point(data = nega_points, aes(x=PC1, y=PC2, size=3), color="red")


# plot_out2 <- plot_out2 + geom_point(data = my_points, aes(x=PC1, y=PC2, size=3, colour=FRNT, alpha = 0.3))
# #plot_out2 <- plot_out2 + geom_point(data = mayo_points, aes(x=PC1, y=PC2, size=3, colour=FRNT, alpha = 0.3))
# plot_out2 <- plot_out2 + geom_point(data = nega_points, aes(x=PC1, y=PC2, size=1.5), color="red", alpha = 0.0)
# # Patient 2
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(11,20),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(data = df3, aes(x = PC11, y = PC21, xend = PC12, yend = PC22), color = "dark red", size = 5)
# #plot_out2 <- plot_out2 + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 2"), nudge_x = 0.65, nudge_y = -0.1)
# 
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(20,21),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(aes(x = PC11, y = PC21, xend = PC12, yend = PC22), data = df3, color = "dark red", size = 5)
# 
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(21,22),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(aes(x = PC11, y = PC21, xend = PC12, yend = PC22), data = df3, color = "dark red", size = 5, arrow = arrow(length = unit(0.5, "inches")))
# 
# 
# #Patient 8
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(6,31),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(data = df3, aes(x = PC11, y = PC21, xend = PC12, yend = PC22), color = "dark green", size = 5)
# 
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(31,32),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(data = df3, aes(x = PC11, y = PC21, xend = PC12, yend = PC22), color = "dark green", size = 5, arrow = arrow(length = unit(0.5, "inches")))
# 
# #Patient 3
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(13,26),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(data = df3, aes(x = PC11, y = PC21, xend = PC12, yend = PC22), color = "dark orange", size = 5)
# 
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(26,39),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(data = df3, aes(x = PC11, y = PC21, xend = PC12, yend = PC22), color = "dark orange", size = 5, arrow = arrow(length = unit(0.5, "inches")))



plot_out2 <- plot_out2 + scale_colour_gradient(trans='log2', limits=c(32,4096), na.value = 'black', name='FRNT')
plot_out2 <- plot_out2 + ggtitle(paste("PCA of MSD [1:", dilution, "] vs FRNT", sep = ""),
                                 subtitle = "SARS-CoV-2 S, mRBD, NP, ORF8") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)) + 
  xlab(paste("PC1 [", trimws(format(round(provar[1]*100,2))), "%]", sep="")) + 
  ylab(paste("PC2 [", trimws(format(round(provar[2]*100,2))), "%]", sep=""))
ggsave(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/PCA_MSD_320_FRNT_with_arrow.tiff",sep = ""), 
       plot = plot_out2,
       width = 23.230,
       height = 14.949,
       units = "cm")
plot_out2

############################### PCA Linear FRNT Cutoff Precision #######################################

#Precision/Recall/Accuracy
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

results.weak <- data.frame(matrix(ncol = 8, nrow = 0))
results.moderate <- data.frame(matrix(ncol = 8, nrow = 0))
results.strong <- data.frame(matrix(ncol = 8, nrow = 0))
for(it in 20:1000)
{
  # Diagonal cutoffs - back-slash
  results.weak <- rbind(results.weak,
                         add_PRA(as.matrix(table(mayo_points$FRNT > it,mayo_points$PC1 * -1.0 - 2.4 < mayo_points$PC2))))
  results.moderate <- rbind(results.moderate,
                       add_PRA(as.matrix(table(mayo_points$FRNT > it,mayo_points$PC1 * -1.0 - 1.6 < mayo_points$PC2))))
  results.strong <- rbind(results.strong,
                           add_PRA(as.matrix(table(mayo_points$FRNT > it,mayo_points$PC1 * -1.0 + 0.6 < mayo_points$PC2))))
}

plot_out2 <- ggplot()
#plot_out2 <- plot_out2 + geom_point(data = df_out, aes(x=PC1, y=PC2, colour=group1, size=1.5))
#  + geom_point(data = df_out, alpha = 0.3, color="white",aes(x=PC1, y=PC2,size=1.5))
plot_out2 <- plot_out2 + scale_size(guide = 'none')
plot_out2 <- plot_out2 + geom_point(data = results.weak, aes(x=Sensitivity, y=Precision, size=1.5, colour=Cutoff))
plot_out2 <- plot_out2 + scale_colour_gradient(trans='log2', limits=c(16,1024), na.value = 'black', name='FRNT')
plot_out2 <- plot_out2 + ggtitle(paste("PR Curve - Weak Cutoff", sep = ""))  +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab(paste("Recall", sep="")) + 
  ylab(paste("Precision", sep="")) +
  scale_y_continuous(labels=function(x) paste0(x*100,"%"), limits = c(0.25,1)) + 
  scale_x_continuous(labels=function(x) paste0(x*100,"%"), limits = c(0.25,1))
#ggsave(paste("20200423_PCA_",dilution,"S2S1RBD_vs_FRNT_chicago_conv.tiff",sep = ""), plot = plot_out2)
plot_out2

#Accuracy plot based on cutoff
plot_out3 <- ggplot() + 
  scale_size(guide = 'none') +
  geom_line(data = results.weak, aes(x = Cutoff, y = Accuracy), linetype = "dotted", size = 1) +
  geom_line(data = results.moderate  , aes(x = Cutoff, y = Accuracy), linetype = "dashed", size = 1) +
  geom_line(data = results.strong  , aes(x = Cutoff, y = Accuracy), linetype = "solid", size = 1) +
  theme_classic() +
  ggtitle(paste("Prediction Accuracy of MSD PCA [1:", dilution, "]", sep = ""),
          subtitle = "SARS-CoV-2 mRBD & S")  +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlim(20,1000) + xlab("Lower Limit Cutoff 1/[FRNT50]") +
  scale_y_continuous(labels=function(x) paste0(x*100,"%"), limits = c(0,1)) +
  scale_x_continuous(trans="log2") +
  coord_cartesian(xlim = c(20,1000), ylim = c(.5,1))

ggsave(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/Accuracy_MSD_FRNT_Prediction.tiff",sep = ""), 
       plot = plot_out3,
       width = 23.230/1.6,
       height = 14.949/1.1,
       units = "cm")
plot_out3


##################### PCA ELISA #########################

my_set <- c("SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD","SARS2.ORF8","SARS2.NP.FL","ChkP62E1","HA","SARS2.NP.RNA")

cdf <- mydata_cat[c("SARS2.NP.FL","SARS2.Trimer", "SARS2.mRBD", "SARS2.ORF8")]
ddf <- mydata_cat[c("SARS2.NP.FL","SARS2.Trimer", "SARS2.mRBD", "SARS2.ORF8")]
ndf <- mynegs[c("SARS2.NP.FL","SARS2.Trimer", "SARS2.mRBD", "SARS2.ORF8")]

# cdf <- mydata2[c("SARS1.bRBD","SARS2.mRBD")]
# ddf <- mydata[c("SARS1.bRBD","SARS2.mRBD")]
# ndf <- mynegs[c("SARS1.bRBD","SARS2.mRBD")]

# cdf <- mydata2[c("SARS1.bRBD","SARS2.mRBD","SARS2.bRBD","SARS2.Trimer","SARS2.NP.RNA","SARS2.NP.FL")]
# ddf <- mydata[c("SARS1.bRBD","SARS2.mRBD","SARS2.bRBD","SARS2.Trimer","SARS2.NP.RNA","SARS2.NP.FL")]
# ndf <- mynegs[c("SARS1.bRBD","SARS2.mRBD","SARS2.bRBD","SARS2.Trimer","SARS2.NP.RNA","SARS2.NP.FL")]

pca_df <- prcomp(cdf, retx=TRUE, center=TRUE, scale=TRUE)
nega_points <- as.data.frame(predict(pca_df, ndf))
mayo_points <- as.data.frame(predict(pca_df, ddf))
mayo_points$FRNT <- mesdata_cat$FRNT
eigs <- pca_df$sdev^2
provar <- eigs / sum(eigs)

my_points <- as.data.frame(predict(pca_df, cdf))
my_points$FRNT <- mesdata_cat$FRNT

plot_out2 <- ggplot()

# #Vertical cutoffs
# plot_out2 <- plot_out2 + geom_vline(xintercept = -2.6, linetype = "dotted", color="dark gray",size = 2)
# plot_out2 <- plot_out2 + geom_vline(xintercept = -0.5, linetype = "dashed", color="dark gray",size = 2)
# plot_out2 <- plot_out2 + geom_vline(xintercept = 0.5, linetype = "solid", color="dark gray",size = 2)
#Forward slash cutoffs
# # S and NP
# plot_out2 <- plot_out2 + geom_abline(intercept = 3.4, slope=1,linetype="dotted",color="dark gray", size=2)
# plot_out2 <- plot_out2 + geom_abline(intercept = 0.6, slope=1,linetype="dashed",color="dark gray",size=2)
# plot_out2 <- plot_out2 + geom_abline(intercept = -0.5, slope=1,linetype="solid",color="dark gray",size=2)
# # S, NP, mRBD
# plot_out2 <- plot_out2 + geom_abline(intercept = -3, slope=-1,linetype="dotted",color="dark gray", size=2)
# plot_out2 <- plot_out2 + geom_abline(intercept = -1, slope=-1,linetype="dashed",color="dark gray",size=2)
# plot_out2 <- plot_out2 + geom_abline(intercept = 0.5, slope=-1,linetype="solid",color="dark gray",size=2)
# S, NP, mRBD, ORF8
plot_out2 <- plot_out2 + geom_abline(intercept = -2.5, slope=-1,linetype="dotted",color="dark gray", size=2)
plot_out2 <- plot_out2 + geom_abline(intercept = -1, slope=-1,linetype="dashed",color="dark gray",size=2)
plot_out2 <- plot_out2 + geom_abline(intercept = 0.5, slope=-1,linetype="solid",color="dark gray",size=2)

plot_out2 <- plot_out2 + scale_size(guide = 'none')
plot_out2 <- plot_out2 + geom_point(data = my_points, aes(x=PC1, y=PC2, size=3, colour=FRNT))
#plot_out2 <- plot_out2 + geom_point(data = mayo_points, aes(x=PC1, y=PC2, size=3, colour=FRNT))
plot_out2 <- plot_out2 + geom_point(data = nega_points, aes(x=PC1, y=PC2, size=1.5), color="red")
plot_out2 <- plot_out2 + geom_point(data = nega_points, aes(x=PC1, y=PC2, size=3), color="red")


# plot_out2 <- plot_out2 + geom_point(data = my_points, aes(x=PC1, y=PC2, size=3, colour=FRNT, alpha = 0.3))
# plot_out2 <- plot_out2 + geom_point(data = mayo_points, aes(x=PC1, y=PC2, size=3, colour=FRNT, alpha = 0.3))
# plot_out2 <- plot_out2 + geom_point(data = nega_points, aes(x=PC1, y=PC2, size=1.5), color="red", alpha = 0.0)
# #Patient 8
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(6,31),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(data = df3, aes(x = PC11, y = PC21, xend = PC12, yend = PC22), color = "dark green", size = 5)
# 
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(31,32),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(data = df3, aes(x = PC11, y = PC21, xend = PC12, yend = PC22), color = "dark green", size = 5, arrow = arrow(length = unit(0.5, "inches")))
# 
# # Patient 2
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(11,20),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(data = df3, aes(x = PC11, y = PC21, xend = PC12, yend = PC22), color = "dark red", size = 5)
# #plot_out2 <- plot_out2 + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 2"), nudge_x = 0.65, nudge_y = -0.1)
# 
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(20,21),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(aes(x = PC11, y = PC21, xend = PC12, yend = PC22), data = df3, color = "dark red", size = 5)
# 
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(21,22),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(aes(x = PC11, y = PC21, xend = PC12, yend = PC22), data = df3, color = "dark red", size = 5, arrow = arrow(length = unit(0.5, "inches")))
# 
# #Patient 3
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(13,26),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(data = df3, aes(x = PC11, y = PC21, xend = PC12, yend = PC22), color = "dark orange", size = 5)
# 
# df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mayo_points[c(26,39),])))[c("PC1","PC2")])))
# plot_out2 <- plot_out2 + geom_segment(data = df3, aes(x = PC11, y = PC21, xend = PC12, yend = PC22), color = "dark orange", size = 5, arrow = arrow(length = unit(0.5, "inches")))



plot_out2 <- plot_out2 + scale_colour_gradient(trans='log2', limits=c(32,4096), na.value = 'black', name='FRNT')
plot_out2 <- plot_out2 + ggtitle(paste("PCA of ELISA [1:", dilution, "] vs FRNT", sep = ""),
                                 subtitle = "SARS-CoV-2 S, mRBD, NP, ORF8") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18,hjust=0.5),
        plot.subtitle = element_text(size = 14, hjust=0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)) + 
  xlab(paste("PC1 [", trimws(format(round(provar[1]*100,2))), "%]", sep="")) + 
  ylab(paste("PC2 [", trimws(format(round(provar[2]*100,2))), "%]", sep=""))
ggsave(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/PCA_ELISA_320_FRNT_without_arrow.tiff",sep = ""), 
       plot = plot_out2,
       width = 23.230,
       height = 14.949,
       units = "cm")
plot_out2

############################   Prediction PR Curves - ELISA   #####################################

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


results.weak <- data.frame(matrix(ncol = 8, nrow = 0))
results.moderate <- data.frame(matrix(ncol = 8, nrow = 0))
results.strong <- data.frame(matrix(ncol = 8, nrow = 0))
for(it in 20:1000)
{
  # # Forward-slash cutoffs
  # results.weak <- rbind(results.weak,
  #                        add_PRA(as.matrix(table(mayo_points$FRNT > it,mayo_points$PC1 * 1.0 + 3.4 > mayo_points$PC2))))
  # results.moderate <- rbind(results.moderate,
  #                      add_PRA(as.matrix(table(mayo_points$FRNT > it,mayo_points$PC1 * 1.0 + 0.6 > mayo_points$PC2))))
  # results.strong <- rbind(results.strong,
  #                       add_PRA(as.matrix(table(mayo_points$FRNT > it,mayo_points$PC1 * 1.0 + -0.5 > mayo_points$PC2))))
  
  # Back-slash cutoffs
  results.weak <- rbind(results.weak,
                        add_PRA(as.matrix(table(mayo_points$FRNT > it,mayo_points$PC1 * -1.0 - 2.5 < mayo_points$PC2))))
  results.moderate <- rbind(results.moderate,
                            add_PRA(as.matrix(table(mayo_points$FRNT > it,mayo_points$PC1 * -1.0 - 1.0 < mayo_points$PC2))))
  results.strong <- rbind(results.strong,
                          add_PRA(as.matrix(table(mayo_points$FRNT > it,mayo_points$PC1 * -1.0 + 0.5 < mayo_points$PC2))))
  
  # # Vertical cutoffs
  # results.weak <- rbind(results.weak,
  #                        add_PRA(as.matrix(table(mayo_points$FRNT > it, mayo_points$PC1 > -2.6))))
  # results.moderate <- rbind(results.moderate,
  #                      add_PRA(as.matrix(table(mayo_points$FRNT > it, mayo_points$PC1 > -0.5))))
  # results.strong <- rbind(results.strong,
  #                          add_PRA(as.matrix(table(mayo_points$FRNT > it, mayo_points$PC1 >  0.5))))
}


############

plot_out2 <- ggplot()
#plot_out2 <- plot_out2 + geom_point(data = df_out, aes(x=PC1, y=PC2, colour=group1, size=1.5))
#  + geom_point(data = df_out, alpha = 0.3, color="white",aes(x=PC1, y=PC2,size=1.5))
plot_out2 <- plot_out2 + scale_size(guide = 'none')
plot_out2 <- plot_out2 + geom_point(data = results.weak, aes(x=RecallT, y=PrecisionT, size=1.5, colour=Cutoff))
plot_out2 <- plot_out2 + scale_colour_gradient(trans='log2', limits=c(16,1024), na.value = 'black', name='FRNT')
plot_out2 <- plot_out2 + ggtitle(paste("PR Curve - Weak Cutoff", sep = ""))  +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab(paste("Recall", sep="")) + 
  ylab(paste("Precision", sep="")) +
  scale_y_continuous(labels=function(x) paste0(x*100,"%"), limits = c(0.25,1)) + 
  scale_x_continuous(labels=function(x) paste0(x*100,"%"), limits = c(0.25,1))
#ggsave(paste("20200423_PCA_",dilution,"S2S1RBD_vs_FRNT_chicago_conv.tiff",sep = ""), plot = plot_out2)
plot_out2

#Accuracy plot based on cutoff
plot_out3 <- ggplot() + 
  scale_size(guide = 'none') +
  geom_line(data = results.weak, aes(x = Cutoff, y = Accuracy), linetype = "dotted", size = 1) +
  geom_line(data = results.moderate  , aes(x = Cutoff, y = Accuracy), linetype = "dashed", size = 1) +
  geom_line(data = results.strong, aes(x = Cutoff, y = Accuracy), linetype = "solid", size = 1) +
  theme_classic() +
  ggtitle(paste("Prediction Accuracy of ELISA PCA [1:", dilution, "]", sep = ""),
          subtitle = "SARS-CoV-2 mRBD & S")  +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab("Lower Limit Cutoff 1/[FRNT50]") +
  scale_y_continuous(labels=function(x) paste0(x*100,"%"), limits = c(0,1)) +
  scale_x_continuous(trans="log2") +
  coord_cartesian(xlim = c(20,1000), ylim = c(.5,1))

ggsave(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/Accuracy_ELISA_FRNT_Prediction.tiff",sep = ""), 
       plot = plot_out3,
       width = 23.230/1.6,
       height = 14.949/1.1,
       units = "cm")
plot_out3
