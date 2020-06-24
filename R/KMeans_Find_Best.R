library(ggplot2)
library(tidyverse)
library(grid)
library(gridExtra)
setwd("~/Dropbox/Fremont/Bioinformatics/ELISA/MayoCOVID19/20200421_ELISA_Processing/")

d <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(d) <- paste(320*2^(0:7),sep="")

# for each number of antigens in a combo
totIter = 1
for (co in 2:6)
{
  # for each antigen combination
  r6cn <- combn(1:6,co)
  for (i in 1:length(r6cn[c(1),]))
  {
    # for each titer
    for (iter in 0:7 )
    {
      dilution <- toString(320*2^iter)
      
      mydata <- read.csv(paste(dilution, ".csv",sep=""),header=FALSE)
      df <- data.frame(mydata)
      colnames(df) <- c("frnt","chk","bsa","ha","mers","s1brbd","s2brbd","s2mrbd","s2spike","orf8","np-rna","np-fl","orf7a")
      # Preselect only the antigens that independently correlate with FRNT values
      df2 <- df[c(6:9,11:12)]
      # Iterate through subsets
      df3 <- df2[r6cn[,c(i)]]
      k<-kmeans(df3,2,nstart=100,iter.max = 10000)
      a <- data.frame(matrix(ncol = 2,nrow = 42))
      for(v in 1:42)
      {
        a[[v,k$cluster[[v]]]] <- df[[v,1,1]]
      }
      #This is the absolute difference of the mean FRNT value within the two clusters
      # for a specific titer
      d[totIter,iter+1] <- abs(diff(colMeans(a, na.rm = TRUE)))[[1]]
      if(d[totIter,iter+1] > 700 & d[totIter,iter+1] < 800.0) print(paste(c(d[totIter,iter+1],dilution,colnames(df3))))
    }
    totIter = totIter + 1
  }
}


d_plot <- d %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

d_plot$colname <- as.numeric(d_plot$colname)
d_plot$rowname <- as.numeric(d_plot$rowname)
d_plot$value <- as.numeric(d_plot$value)

q <- ggplot(data = d_plot, aes(x = reorder(colname,X=as.numeric(colname)), y=value)) +
  geom_jitter(color = ifelse(d_plot$value > 800, "black","gray"), 
             show.legend = FALSE, 
             size = ifelse(d_plot$value > 800, 3,2),
             width = 0.2) +
  geom_boxplot(outlier.size = 2.0, 
               outlier.shape = NA,
               fill="gray") +
  geom_hline(yintercept = 774.25, linetype = "dashed", color = "gray") +
  ggtitle("ΔAverage FRNT of ELISA k-means clusters") + 
  xlab("Reciprocal ELISA Sera Dilution") +
  ylab("Absolute difference of cluster FRNT values") +
  theme_classic() + 
  theme(axis.title = element_text(size = 18,hjust=0.5),
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.text = element_text(size = 16))
  
ggsave("Average_FRNT_of_ELISA_kmeans_clusters.tiff", q)
q


