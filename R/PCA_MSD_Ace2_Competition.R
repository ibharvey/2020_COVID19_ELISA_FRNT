library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(corrplot)
library(abind)
setwd("~/Dropbox/Fremont/Bioinformatics/ELISA/Ace2_Competition/")

dilution = 80

temp1 <- read.csv(paste("Mesoscale_Ace2_Competition_1/",dilution,"_MSD_1.csv",sep=""),header=T)
temp2 <- read.csv("Mesoscale_Ace2_Competition_2/2BK46AS018_2020-06-16-210443.txt_Rinput.csv",header=T)
temp_combine <- as.data.frame(rowMeans(abind(temp1,temp2,along=3), dims=2))
mydata.neg <- temp_combine[c(81:96),]
mydata.conv <- temp_combine[c(1:80),]


temp3 <- read.csv("Mesoscale_Ace2_Competition_2/2BK46AL049_2020-06-16-210905.txt_Rinput.csv",header=T)
temp3.a <- temp3[c(1:42),]
temp3.b <- temp3[c(49:90),]
mydata.mayo <- as.data.frame(rowMeans(abind(temp3.a,temp3.b,along=3),dims = 2))

################## Correlations ##########################
cor_data <- rbind(mydata.conv, mydata.mayo)[c("FRNT","SARS2.Trimer","SARS2.mRBD","SARS2.bRBD","SARS1.bRBD","MERS.bRBD")]
cor_data$FRNT <- log(cor_data$FRNT)
colnames(cor_data) <- c("log2(FRNT50)","SARS2-CoV-2 bRBD",
                        "SARS-CoV-2 mRBD","SARS-CoV-2 bRBD",
                        "SARS-CoV bRBD","MERS-CoV bRBD")
res <- cor(cor_data)
# 
# #ELISA on rows, Mesoscale on columns
tiff(filename="~/Dropbox/Fremont/Bioinformatics/ELISA/Ace2_Competition/Corrplot_AceMeso.tiff")
corrplot(res,tl.col = "black", method = "color", 
         addCoef.col = "black",
         tl.cex = 1.6, cl.cex = 1.4, cl.align.text = "l",
         number.cex = 1.2, type = "upper")
dev.off()

mydata.cm <- rbind(mydata.mayo, mydata.conv)

# With assumption of negative FRNT in common coronavirus
# mydata.neg$FRNT <- 20
# mydata.cm <- rbind(mydata.mayo, mydata.conv, mydata.neg)

gplot <- ggplot(data = mydata.cm, aes(x = FRNT, y=SARS2.bRBD)) +
  geom_point(na.rm=T, size = 3) + 
  # # S2 S
  # geom_abline(intercept = 1.5e6, slope=0,linetype="solid",color="dark gray",size=2) + 
  # geom_abline(intercept = 6e6, slope=0,linetype="dotted",color="dark gray",size=2)+
  # S2 bRBD
  geom_abline(intercept = 1e4, slope=0,linetype="solid",color="dark gray",size=2) + 
  geom_abline(intercept = 5e4, slope=0,linetype="dotted",color="dark gray",size=2)+
  scale_x_continuous(trans="log2") +
  scale_y_continuous(trans="log2", breaks = c(5e4, 5e5, 5e6)) +
  theme_classic() +
  ggtitle(paste("SARS2-CoV-2 bRBD vs FRNT50 [1:", dilution, "]", sep = ""))  +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab("FRNT50") + ylab("Electrochemoluminescence") 

# Patient 2
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(11,20),])))[c("FRNT","SARS2.bRBD")])))
gplot <- gplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), color = "dark red", size = 3)
#gplot <- gplot + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 2"), nudge_x = 0.65, nudge_y = -0.1)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(20,21),])))[c("FRNT","SARS2.bRBD")])))
gplot <- gplot + geom_segment(aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), data = df3, color = "dark red", size = 3)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(21,22),])))[c("FRNT","SARS2.bRBD")])))
gplot <- gplot + geom_segment(aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), data = df3, color = "dark red", size = 3, arrow = arrow(length = unit(0.3, "inches")))


#Patient 8
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(6,31),])))[c("FRNT","SARS2.bRBD")])))
gplot <- gplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), color = "dark green", size = 3)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(31,32),])))[c("FRNT","SARS2.bRBD")])))
gplot <- gplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), color = "dark green", size = 3, arrow = arrow(length = unit(0.3, "inches")))

#Patient 3
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(13,26),])))[c("FRNT","SARS2.bRBD")])))
gplot <- gplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), color = "dark orange", size = 3)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(26,39),])))[c("FRNT","SARS2.bRBD")])))
gplot <- gplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), color = "dark orange", size = 3, arrow = arrow(length = unit(0.3, "inches")))


################# Plot linear
qplot <- ggplot(data = mydata.cm, aes(x = FRNT, y=SARS2.bRBD)) +
  geom_point(na.rm=T, size = 2) + 
  # # S2 S
  # geom_abline(intercept = 1.5e6, slope=0,linetype="solid",color="dark gray",size=2) + 
  # geom_abline(intercept = 6e6, slope=0,linetype="dotted",color="dark gray",size=2)+
  geom_abline(intercept = 1e4, slope=0,linetype="solid",color="dark gray",size=2) + 
  geom_abline(intercept = 5e4, slope=0,linetype="dotted",color="dark gray",size=2)+
  scale_x_continuous(trans="log2") +
  theme_classic() +
  ggtitle(paste("SARS2-CoV-2 bRBD vs FRNT50 [1:", dilution, "]", sep = ""))  +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab("FRNT50") + ylab("Electrochemoluminescence") 

# Patient 2
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(11,20),])))[c("FRNT","SARS2.bRBD")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), color = "dark red", size = 2)
#qplot <- qplot + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 2"), nudge_x = 0.65, nudge_y = -0.1)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(20,21),])))[c("FRNT","SARS2.bRBD")])))
qplot <- qplot + geom_segment(aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), data = df3, color = "dark red", size = 2)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(21,22),])))[c("FRNT","SARS2.bRBD")])))
qplot <- qplot + geom_segment(aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), data = df3, color = "dark red", size = 2, arrow = arrow(length = unit(0.2, "inches")))


#Patient 8
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(6,31),])))[c("FRNT","SARS2.bRBD")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), color = "dark green", size = 2)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(31,32),])))[c("FRNT","SARS2.bRBD")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), color = "dark green", size = 2, arrow = arrow(length = unit(0.2, "inches")))

#Patient 3
df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(13,26),])))[c("FRNT","SARS2.bRBD")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), color = "dark orange", size = 2)

df3 <- as.data.frame(t(unlist(do.call(rbind.data.frame, list(c(mydata.cm[c(26,39),])))[c("FRNT","SARS2.bRBD")])))
qplot <- qplot + geom_segment(data = df3, aes(x = FRNT1, y = SARS2.bRBD1, xend = FRNT2, yend = SARS2.bRBD2), color = "dark orange", size = 2, arrow = arrow(length = unit(0.2, "inches")))

qplot
ggsave(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/Ace2_Competition/Log2_FRNT_Trimer_122_datapoints.tiff",sep = ""), 
       plot = gplot,
       width = 23.230/1.5,
       height = 14.949/1.5,
       units = "cm")

ggsave(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/Ace2_Competition/Linear_FRNT_Trimer_122_datapoints.tiff",sep = ""), 
       plot = qplot,
       width = 23.230,
       height = 14.949*1.2,
       units = "cm")

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
results.strong <- data.frame(matrix(ncol = 8, nrow = 0))
for(it in 20:1000)
{
  results.weak <- rbind(results.weak,
                         add_PRA(as.matrix(table(mydata.cm$FRNT > it,mydata.cm$SARS2.bRBD<5e4))))
  results.strong <- rbind(results.strong,
                           add_PRA(as.matrix(table(mydata.cm$FRNT > it,mydata.cm$SARS2.bRBD<1e4))))
}


#Accuracy plot based on cutoff
plot_out3 <- ggplot() + 
  scale_size(guide = 'none') +
  geom_line(data = results.weak, aes(x = Cutoff, y = Accuracy), linetype = "dotted", size = 1) +
  geom_line(data = results.strong, aes(x = Cutoff, y = Accuracy), linetype = "solid", size = 1) +
  theme_classic() +
  ggtitle(paste("Prediction Accuracy of Ace2 Competition"),
          subtitle = paste("SARS2-CoV-2 bRBD ", "[1:", dilution, "]",sep = "")) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust=0.5),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab("Lower Limit Cutoff 1/[FRNT50]") +
  scale_y_continuous(labels=function(x) paste0(x*100,"%"), limits = c(0,1)) +
  scale_x_continuous(trans="log2") +
  coord_cartesian(xlim = c(20,1000), ylim = c(.5,1))

plot_out3
ggsave(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/Ace2_Competition/Accuracy_FRNT_Trimer_122_datapoints.tiff",sep = ""), 
       plot = plot_out3,
       width = 23.230/1.5,
       height = 14.949/1.5,
       units = "cm")

# ########################## PCA ###############################
# 
# df.conv <- mydata.conv[c("SARS2.Trimer","SARS2.bRBD","SARS2.Trimer")]
# pca_df <- prcomp(df.conv, retx=TRUE, center=TRUE, scale=TRUE)
# 
# eigs <- pca_df$sdev^2
# provar <- eigs / sum(eigs)
# 
# my_points.conv <- as.data.frame(predict(pca_df, df.conv))
# my_points.conv$FRNT <- mydata.conv$FRNT
# df.neg <- mydata.neg[c("SARS2.Trimer","SARS2.bRBD","SARS2.Trimer")]
# my_points.neg <- as.data.frame(predict(pca_df, df.neg))
# df.mayo <- mydata.mayo[c("SARS2.Trimer","SARS2.bRBD","SARS2.Trimer")]
# my_points.mayo <- as.data.frame(predict(pca_df, df.mayo))
# my_points.mayo$FRNT <- mydata.mayo$FRNT
# 
# plot_out2 <- ggplot()
# plot_out2 <- plot_out2 + scale_size(guide = 'none')
# 
# plot_out2 <- plot_out2 + geom_point(data = my_points.neg, aes(x=PC1, y=PC2, size=1.5), color="red")
# plot_out2 <- plot_out2 + geom_point(data = my_points.conv, aes(x=PC1, y=PC2, size=3, colour=FRNT))
# plot_out2 <- plot_out2 + geom_point(data = my_points.mayo, aes(x=PC1, y=PC2, size=3, colour=FRNT))
# plot_out2 <- plot_out2 + geom_point(data = my_points.neg, aes(x=PC1, y=PC2, size=3), color="red")
# 
# plot_out2 <- plot_out2 + geom_abline(intercept = -6.25, slope=2,linetype="dotted",color="dark gray",size=2)
# plot_out2 <- plot_out2 + geom_abline(intercept = -3, slope=2,linetype="dashed",color="dark gray",size=2)
# plot_out2 <- plot_out2 + geom_abline(intercept = 1, slope=2,linetype="solid",color="dark gray",size=2)
# 
# plot_out2 <- plot_out2 + scale_colour_gradient(trans='log2', limits=c(32,4096), na.value = 'black', name='FRNT')
# plot_out2 <- plot_out2 + ggtitle(paste("PCA of MSD Ace2-Comp [1:", dilution, "] vs FRNT", sep = ""),
#                                  subtitle = "S2.bRBD S2.mRBD S2.Trimer") +
#   theme_minimal() +
#   theme(plot.title = element_text(size = 18, hjust = 0.5),
#         plot.subtitle = element_text(size = 14, hjust = 0.5),
#         axis.title = element_text(size = 18),
#         axis.text = element_text(size = 16)) + 
#   xlab(paste("PC1 [", trimws(format(round(provar[1]*100,2))), "%]", sep="")) + 
#   ylab(paste("PC2 [", trimws(format(round(provar[2]*100,2))), "%]", sep=""))
# ggsave(paste("~/Dropbox/Fremont/Bioinformatics/ELISA/PCA_MSD_320_FRNT_arrow.tiff",sep = ""), 
#        plot = plot_out2,
#        width = 23.230,
#        height = 14.949,
#        units = "cm")
# plot_out2
# 
# 
# results.weak <- data.frame(matrix(ncol = 8, nrow = 0))
# results.far <- data.frame(matrix(ncol = 8, nrow = 0))
# results.strong <- data.frame(matrix(ncol = 8, nrow = 0))
# for(it in 20:1000)
# {
#   results.weak <- rbind(results.weak,
#                          add_PRA(as.matrix(table(my_points.mayo$FRNT > it,my_points.mayo$PC1 * 2.0 - 6.25 < my_points.mayo$PC2))))
#   results.far <- rbind(results.far, 
#                        add_PRA(as.matrix(table(my_points.mayo$FRNT > it,my_points.mayo$PC1 * 2.0 - 3 < my_points.mayo$PC2))))
#   results.strong <- rbind(results.strong,
#                            add_PRA(as.matrix(table(my_points.mayo$FRNT > it,my_points.mayo$PC1 * 2.0 + 1 < my_points.mayo$PC2))))
#   
# }
# 
# 
# #Accuracy plot based on cutoff
# plot_out3 <- ggplot() + 
#   scale_size(guide = 'none') +
#   geom_line(data = results.weak, aes(x = Cutoff, y = Accuracy), linetype = "dotted", size = 1) +
#   geom_line(data = results.far, aes(x = Cutoff, y = Accuracy), linetype = "dashed", size = 1) +
#   geom_line(data = results.strong, aes(x = Cutoff, y = Accuracy), linetype = "solid", size = 1) +
#   theme_classic() +
#   ggtitle(paste("Prediction Accuracy of Ace2 Competition [1:", dilution, "]", sep = ""),
#           subtitle = "S2.bRBD S2.mRBD S2.Trimer")  +
#   theme(plot.title = element_text(size = 18, hjust = 0.5),
#         plot.subtitle = element_text(size = 14, hjust=0.5),
#         axis.text = element_text(size=16),
#         axis.title = element_text(size=18)) +
#   xlab("Lower Limit Cutoff 1/[FRNT50]") +
#   scale_y_continuous(labels=function(x) paste0(x*100,"%"), limits = c(0,1)) +
#   scale_x_continuous(trans="log2") +
#   coord_cartesian(xlim = c(20,1000), ylim = c(.5,1))
# 
# plot_out3


