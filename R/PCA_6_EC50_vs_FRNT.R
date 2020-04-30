library(ggplot2)
library(grid)
library(gridExtra)
setwd("~/Box/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/")

mydata <- read.csv(paste("FRNT_and_EC50_duplicates.csv",sep=""),header=TRUE)
df <- data.frame(mydata)
df2 <- df[c(4:7,9:10)]
pca_df2 <- prcomp(df2)
eigs <- pca_df2$sdev^2
provar <- eigs / sum(eigs)
df_out <- as.data.frame(pca_df2$x)
df_out$group <- df$FRNT

plot_out <- ggplot(data = df_out, aes(PC1,PC2))
plot_out <- plot_out + geom_point(data = df_out, aes(x=PC1, y=PC2, colour=group,size=1.5))  + scale_size(guide = 'none')
plot_out <- plot_out + scale_colour_gradient(trans='log2', limits=c(64,2048), na.value = 'black', name='FRNT')

plot_out <- plot_out + ggtitle(paste("PCA of Six Correlated ELISA EC50s vs FRNT", sep = "")) +
  theme(plot.title = element_text(size = 16)) + 
  xlab(paste("PC1 [", trimws(format(round(provar[1]*100,2))), "%]", sep="")) + 
  ylab(paste("PC2 [", trimws(format(round(provar[2]*100,2))), "%]", sep=""))

plot_out2 <- plot_out
ggsave("PCA_6_sigEC50_vs_FRNT.tiff", plot_out)

#Patient 2
df3 <- do.call(rbind.data.frame, list(c(pca_df2$x[c(11,20),])))
colnames(df3) <- c("x1","x2","y1","y2")
plot_out <- plot_out + geom_segment(data = df3, aes(x = x1, y = y1, xend = x2, yend = y2), color = "dark red", size = 1.5)
plot_out <- plot_out + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 2"), nudge_x = 0.2, nudge_y = 2500)

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
plot_out <- plot_out + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 8"), nudge_x =18000, nudge_y = -8000)

df3 <- do.call(rbind.data.frame, list(c(pca_df2$x[c(31,32),])))
colnames(df3) <- c("x1","x2","y1","y2")
plot_out <- plot_out + geom_segment(data = df3, aes(x = x1, y = y1, xend = x2, yend = y2), color = "dark green", size = 1.5, arrow = arrow(length = unit(0.15, "inches")))
plot_out
#Patient 18
df3 <- do.call(rbind.data.frame, list(c(pca_df2$x[c(23,34),])))
colnames(df3) <- c("x1","x2","y1","y2")
plot_out <- plot_out + geom_segment(data = df3, aes(x = x1, y = y1, xend = x2, yend = y2), color = "dark orange", size = 1.5)
plot_out <- plot_out + geom_text(data = df3, aes(x = x1, y = y1, label = "Patient 18"), nudge_x = 2000, nudge_y = -11000)

df3 <- do.call(rbind.data.frame, list(c(pca_df2$x[c(34,40),])))
colnames(df3) <- c("x1","x2","y1","y2")
plot_out <- plot_out + geom_segment(data = df3, aes(x = x1, y = y1, xend = x2, yend = y2), color = "dark orange", size = 1.5, arrow = arrow(length = unit(0.15, "inches")))

plot_out
ggsave(paste("20200428_PCA_Six_sigEC50_vs_FRNT_labelled.tiff",sep = ""), plot = plot_out)


k <- kmeans(df2,2,nstart=100,iter.max=10000)
a <- data.frame(matrix(ncol = 2, nrow = 42))
for(v in 1:42) {a[[v,k$cluster[[v]]]] <- df[[v,1,1]]}
cTemp <- colMeans(a,na.rm=T)
bigCluster <- ifelse(cTemp[[1]] > cTemp[[2]], 1, 2)
df_out$kclusterN <- ifelse(k$cluster == bigCluster,"green","red")

plot_out2 <- ggplot(data = df_out, aes(PC1,PC2))
plot_out2 <- plot_out2 + geom_point(data = df_out, colour=df_out$kclusterN, aes(x=PC1, y=PC2, size=1.5, fill=group), pch=21) + scale_size(guide = 'none')

plot_out2 <- plot_out2 + scale_fill_gradient(trans='log2', limits=c(64,2048), na.value = 'black', name='FRNT')

plot_out2 <- plot_out2 + ggtitle(paste("PCA of Six Correlated ELISA EC50s vs FRNT", sep = "")) +
  xlab(paste("PC1 [", trimws(format(round(provar[1]*100,2))), "%]", sep="")) + 
  ylab(paste("PC2 [", trimws(format(round(provar[2]*100,2))), "%]", sep=""))

plot_out2

ggsave("20200428_PCA_Six_sigEC50_vs_FRNT_kclustered.tiff",plot_out2)
