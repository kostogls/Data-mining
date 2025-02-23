# load libraries
library(ggplot2)
library(dplyr)
library(ggalt)
library(ggforce)
library("factoextra")
#install.packages("ggalt")
#install.packages("ggforce")
#install.packages("factoextra")

# load dataset
wines <- as.matrix(read.table("wine.data", header = FALSE, sep = ","))
# circular data set, uncomment to use
# it only needs normalization, not pca

#circle <-  function(x = x, y = y, r = radius, n = n.faces){
#t <- seq(from = 0, to = 2 * pi, length = n + 1)[-1]
#t <- cbind(x = x + r * sin(t), y = y+ r * cos(t))
#t <- rbind(t, t[1,])
#return(t)
#}

#dataset<- data.frame(csr(circle(0, 0, 100, 30), 1000))

# it only needs normalization, not pca
#Uncomment to use
#size <- 2000             
#set.seed(1) 
#x <- rnorm(size)          
#y <-rnorm(size) 
#dataset <-data.frame(x,y)

#fviz_nbclust(dataset, kmeans, method='silhouette')
fviz_nbclust(wines, kmeans, method='silhouette')
i = 1

#Features
column_names = list("Index", "Alcohol", "Malic acid", "Ash", "Alcalinity of ash", "Magnesium", "Total phenols", "Flavanoids", 
                    "Nonflavanoid phenols", "Proanthocyanins", "Color intensity", "Hue", "OD280/OD315 of diluted wines", 
                    "Proline")

for (index in colnames(wines)) {
  colnames(wines)[i] = column_names[i]
  i = i + 1
}

# Select rows and columns
n = dim(wines)[1]
d = dim(wines)[2]

#for random
#n = dim(dataset)[1]
#d = dim(dataset)[2]

# Preprocessing(normalization)
for (feauture in 2:d) {
  wines[, feauture] = wines[, feauture] / max(wines[, feauture])
}
## only perform it when normalization is desired
#Preprocessing(normalization)
#for (feauture in 1:d) {
#  dataset[, feauture] = dataset[, feauture] / max(dataset[, feauture])
#}


# delete first row which is index
wines_p=wines[,2:14]

# perform pca
res.pca <- prcomp(wines_p,center=TRUE)
fviz_eig(res.pca)
fviz_pca_ind(res.pca)
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
ind.coord = as.matrix(ind.coord)
ind.coord = ind.coord[,1:2]


## function that calculates euclidean distance
Distance <- function(X, c) {
  dists = apply(X, 1, function(point)
    sapply(1:nrow(c), function(dim)
      dist(rbind(point, c[dim, ]))))
  return(t(dists))
}

k_means <- function(data, k, max_iter=100){
  
  
  # Select rows and columns
  n = dim(data)[1]
  d = dim(data)[2]
  
  
  ## select 3 random points for initial means
  ids = sample(1:n, k)
  centroids = data[ids, ]
  #dist = matrix(0, nrow = 3, ncol = 1)
  ## initialize matrices
  point_centroid = matrix(0, nrow = n, ncol = 1)
  mean_centroid = matrix(0, nrow = k, ncol = d+1)
  dist = matrix(0, nrow = n, ncol = d)
  
  times = 0
  epsilon=10e-5
  
  while (times < max_iter) {
    ## find the distance of every point to the 3 clusters
    dist=Distance(data, centroids)
    ## store the minimum distance
    point_c <- apply(dist, 1, which.min)
    point_centroid = as.matrix(point_c)
    ## find the new centroids
    for (i in 1:n) {
      for (j in 1:d) {
        mean_centroid[point_centroid[i,1], j] = mean_centroid[point_centroid[i,1], j] + data[i,j] 
      }
      mean_centroid[point_centroid[i,1], d+1] = mean_centroid[point_centroid[i,1], d+1] + 1
    }
    for (i in 1:k){
      mean_centroid[i,] = mean_centroid[i,] / mean_centroid[i, d+1]
    }
    # store mean centroid to an array A to delete the last column
    A = mean_centroid
    A<-A[,-(d+1)]
    ## find the change in the centroids
    delta = sum((A - centroids) ^ 2)
    ## if the change is bigger than the error, then continue the loop else break it
    if (delta > epsilon) {
      ## use new centroids for next iteration
      centroids = A
    } else {
      return(point_centroid)
      break
    }
    
    times = times + 1
    return(point_centroid)
  }
}

# tried to implement Jaccard similarity metric
# a are the predictions and b the labels
jaccard <- function(a) {
  count = 0
  element = which.max(a[1:59, ])
  for (i in 1:59) {
    if (a[i] == element) {
      count = count + 1
    }
  }
  print(count)
  
  element = which.max(a[60:130, ])
  for (i in 60:130) {
    if (a[i] == element) {
      count = count + 1
    }
  }
  print(count)
  
  element = which.max(a[131:171, ])
  for (i in 131:171) {
    if (a[i] == element) {
      count = count + 1
    }
  }
  print(count)
  
  count = count/171
  return(count)
  
}
#kmeans for random
#km = k_means(dataset,3)
## k means without pca
#km = k_means(wines_p,3)
# jac = jaccard(wines[,1], c(km))
# print(jac)
## with pca
km = k_means(ind.coord,3)
#print(jaccard(km))
km = as.character(km) 
x <- ind.coord[,1]
y <- ind.coord[,2]
#x <- dataset[,1]
#y <- dataset[,2]
df <- data.frame(x,y,km)
## plot the k means results
ggplot(df,aes(x,y,colour=km))+geom_point(aes(color = km))
# plot with ellipse
ggplot(df,aes(x=x,y=y))+
  geom_mark_ellipse(expand = 0,aes(fill=km))+
  theme_bw()+
  geom_point()
  
dev.off()

