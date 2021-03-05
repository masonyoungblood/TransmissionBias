library(data.table)
library(dtwclust)

load("data.RData")

#log transform each syllable individually
data$meanfreq <- sapply(1:nrow(data), function(x){log(data$meanfreq[[x]])})

#z-score normalization across entire dataset
m <- mean(unlist(data$meanfreq))
sd <- sd(unlist(data$meanfreq))
data$meanfreq <- sapply(1:nrow(data), function(x){(data$meanfreq[[x]]-m)/sd})

#set window size
window <- round(mean(lengths(data$meanfreq))*0.10) #10

#calculate distances
distances <- proxy::dist(data$meanfreq, method = "dtw_basic", window.size = window, normalize = TRUE)

#strip everything from distances
attr(distances, "dimnames") <- NULL
attr(distances, "method") <- NULL
attr(distances, "call") <- NULL
class(distances) <- "matrix"

#replace infinite distances with the maximum value
distances[which(distances == Inf)] <- max(distances[-which(distances == Inf)])

#convert to format compatible with hclust
dtw_dist <- stats::as.dist(distances)

#run clustering
clustering <- fastcluster::hclust(dtw_dist, method = "average")

#save clustering
save(clustering, file = "clustering.RData")

#run hybrid tree cut
hybrid_cut <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = distances, method = "hybrid", minClusterSize = 1, deepSplit = 3)

#save hybrid tree cut
save(hybrid_cut, file = "hybrid_cut.RData")
