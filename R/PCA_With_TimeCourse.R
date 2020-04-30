library(ggplot2)
library(grid)
library(gridExtra)
setwd("~/Box/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/")

dilution <- "320"

mydata <- read.csv(paste(dilution, ".csv",sep=""),header=FALSE)
df <- data.frame(mydata)
df2 <- df[c(6,8)]
pca_df2 <- prcomp(df2, retx=TRUE, center=TRUE, scale=TRUE)
#Flip the PCA2 values
#pca_df2$x[,2] <- -pca_df2$x[,2]
eigs <- pca_df2$sdev^2
provar <- eigs / sum(eigs)
df_out <- as.data.frame(pca_df2$x)
df_out$group1 <- df$V1


plot_out <- ggplot()

mynegs <- read.csv(paste("Chicago_duplicates_", dilution, ".csv",sep=""),header=TRUE)
ndf2 <- mynegs[c(5,7)]
colnames(ndf2) <- c("V6","V8")
neg_points <- as.data.frame(predict(pca_df2, ndf2))



# Flip PCA2
#neg_points$PC2 <- -neg_points$PC2




plot_out <- plot_out + geom_point(data = neg_points, aes(x=PC1, y=PC2, size=1.5), color="red")


plot_out <- plot_out + geom_point(data = df_out, aes(x=PC1, y=PC2, colour=group1, size=1.5)) + scale_size(guide = 'none')

plot_out <- plot_out + scale_colour_gradient(trans='log2', limits=c(64,2048), na.value = 'black', name='FRNT')

plot_out <- plot_out + ggtitle(paste("PCA of bS1-RBD mS2-RBD ELISA [1:", dilution, "] vs FRNT", sep = "")) +
  theme(plot.title = element_text(size = 16)) + 
  xlab(paste("PC1 [", trimws(format(round(provar[1]*100,2))), "%]", sep="")) + 
  ylab(paste("PC2 [", trimws(format(round(provar[2]*100,2))), "%]", sep=""))

ggsave(paste("20200423_PCA_",dilution,"S1_S2_RBD_vs_FRNT.tiff",sep = ""), plot = plot_out)
plot_out2 <- plot_out
plot_out

#Patient 2
df3 <- do.call(rbind.data.frame, list(c(pca_df2$x[c(11,20),])))
colnames(df3) <- c("x1","x2","y1","y2")
plot_out <- plot_out + geom_segment(data = df3, aes(x = x1, y = y1, xend = x2, yend = y2), color = "dark red", size = 1.5)
plot_out <- plot_out + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 2"), nudge_x = -0.6, nudge_y = -0.05)

df3 <- do.call(rbind.data.frame, list(c(pca_df2$x[c(20,21),])))
colnames(df3) <- c("x1","x2","y1","y2")
plot_out <- plot_out + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df3, color = "dark red", size = 1.5)

df3 <- do.call(rbind.data.frame, list(c(pca_df2$x[c(21,22),])))
colnames(df3) <- c("x1","x2","y1","y2")
plot_out <- plot_out + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df3, color = "dark red", size = 1.5, arrow = arrow(length = unit(0.15, "inches")))


#Patient 8
df3 <- do.call(rbind.data.frame, list(c(pca_df2$x[c(6,31),])))
colnames(df3) <- c("x1","x2","y1","y2")
plot_out <- plot_out + geom_segment(data = df3, aes(x = x1, y = y1, xend = x2, yend = y2), color = "dark green", size = 1.5)
plot_out <- plot_out + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 8"), nudge_x = -1.5, nudge_y = 0.4)

df3 <- do.call(rbind.data.frame, list(c(pca_df2$x[c(31,32),])))
colnames(df3) <- c("x1","x2","y1","y2")
plot_out <- plot_out + geom_segment(data = df3, aes(x = x1, y = y1, xend = x2, yend = y2), color = "dark green", size = 1.5, arrow = arrow(length = unit(0.15, "inches")))

#Patient 18
df3 <- do.call(rbind.data.frame, list(c(pca_df2$x[c(23,34),])))
colnames(df3) <- c("x1","x2","y1","y2")
plot_out <- plot_out + geom_segment(data = df3, aes(x = x1, y = y1, xend = x2, yend = y2), color = "dark orange", size = 1.5)
plot_out <- plot_out + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 18"), nudge_x = -.15, nudge_y = 0.7)

df3 <- do.call(rbind.data.frame, list(c(pca_df2$x[c(34,40),])))
colnames(df3) <- c("x1","x2","y1","y2")
plot_out <- plot_out + geom_segment(data = df3, aes(x = x1, y = y1, xend = x2, yend = y2), color = "dark orange", size = 1.5, arrow = arrow(length = unit(0.15, "inches")))

ggsave(paste("20200423_PCA_",dilution,"S1_S2_RBD_vs_FRNT_labelled.tiff",sep = ""), plot = plot_out)
plot_out3 <- plot_out
plot_out

mynegs <- read.csv(paste("Chicago_duplicates_", dilution, ".csv",sep=""),header=TRUE)
ndf2 <- mynegs[c(5,7)]
colnames(ndf2) <- c("V6","V8")
neg_points <- as.data.frame(predict(pca_df2, ndf2))



# Flip PCA2
#neg_points$PC2 <- -neg_points$PC2




plot_out <- plot_out + geom_point(data = neg_points, aes(x=PC1, y=PC2, size=1.5), color="red")
ggsave(paste("20200423_PCA_",dilution,"S1_S2_RBD_vs_FRNT_labelled_chicago.tiff",sep = ""), plot = plot_out)

plot_out2 <- plot_out2 + geom_point(data = neg_points, aes(x=PC1, y=PC2, size=1.5), color="red")
ggsave(paste("20200423_PCA_",dilution,"S1_S2_RBD_vs_FRNT_chicago.tiff",sep = ""), plot = plot_out2)
plot_out2

p <- ggplot() + geom_point(data = df, aes(x=V6, y=V7, colour=V1, size = 2)) +
  scale_size(guide = 'none') +
  scale_colour_gradient(trans='log2', limits=c(64,2048), na.value = 'black', name='FRNT') +
  ggtitle(paste("S1-RBD S2-RBD ELISA [1:", dilution, "] vs FRNT", sep = "")) +
  theme(plot.title = element_text(size = 16)) + 
  xlab("SARS1 bRBD Absorbance (450nm)") + 
  ylab("SARS2 bRBD Absorbance (450nm)")

ggsave("20200428_S1S2_RBD_ELISA_640_vs_FRNT.tiff",p)

plot_out3