library(gplots)

# First, set the path to the working directory. You will need to change this.
setwd("/Users/blancha/Desktop/R_workshop")

# Read tab-separated value into data frame. Set fill to TRUE to add NAs when the rows are imcomplete
data <- read.table("deseq/cyto_vs_total.txt", fill=TRUE)

# Keep only the 1st 30 genes in the file
data <- data[1:30,]

# Keep only the counts columns
data <- subset(data, select=c(Mean_cyto, Mean_total))

# Set a logarithmic scale to better distinguish the samples
data <- log(data)

# Draw the heatmap of the first 30 genes in the file.
distance <- dist(data)
heatmap.2(as.matrix(data), trace="none", margin=c(16, 10), dendrogram="row")