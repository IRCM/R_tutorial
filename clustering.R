# Set path to R workshop directory. You will need to change this.
setwd("/Users/blancha/Desktop/R_workshop")

# Read tab-separated file into data frame
data <- read.table("deseq/normalized_counts.txt")

# Calculate distance between samples.
distance <- dist(t(data))

# Compute clustering
cluster <- hclust(distance)

# Plot clustering
plot(cluster)


