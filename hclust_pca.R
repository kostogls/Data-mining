# load libraries
#install.packages("factoextra")
#install.packages("igraph")
library("factoextra")
library("igraph")


# read wine data
wines <- as.matrix(read.table("wine.data", header = FALSE, sep = ","))

# circular data, uncomment to use  
# it only needs normalization, not pca
# the way of generating circular data was found online
#circle <-  function(x = x, y = y, r = radius, n = n.faces){
#  t <- seq(from = 0, to = 2 * pi, length = n + 1)[-1]
#  t <- cbind(x = x + r * sin(t), y = y+ r * cos(t))
#  t <- rbind(t, t[1,])
#  return(t)
#}
#dataset<- data.frame(csr(circle(0, 0, 100, 30), 1000))

#size <- 2000             
#set.seed(1) 
#x <- rnorm(size)          
#y <-rnorm(size) 
#dataset <-data.frame(x,y)

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

##uncomment below if circular data is used
#n = dim(dataset)[1]
#d = dim(dataset)[2]

# only perform it whene normalization is desired. Else don't run it.
# Preprocessing (normalization)
for (feauture in 2:d) {
  wines[, feauture] = wines[, feauture] / max(wines[, feauture])
}

## uncomment below if circular data is used
#for (feauture in 1:d) {
#dataset[, feauture] = dataset[, feauture] / max(dataset[, feauture])
#}

# delete first row which is index
wines_p = wines[,2:14]

# perform pca
res.pca <- prcomp(wines_p,center=TRUE)
fviz_eig(res.pca)
fviz_pca_ind(res.pca)
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
ind.coord = as.matrix(ind.coord)
# matrix that holds the results of pca
ind.coord = ind.coord[,1:2]

# hclustering function  
hc <-function (data)
{
  # calculate the distance matrix
  dist_matrix <- as.matrix(dist(data))
  print(dist_matrix)
  # change zeros with NA to avoid them
  dist_matrix[dist_matrix==0] <- NA
  
  #In order to find pairs of a tree, treeGroup is made a negative matrix from 1 to nrow(dist_matrix)
  # and then treeGroup include a positive number after combine two attributes that are the smallest distance if the group used previously
  tree_Group <- -(1:nrow(dist_matrix))
  
  # put merged data in variable mergeData
  mergeData <- matrix(0, nrow = nrow(dist_matrix) - 1, ncol = 2)
  mergerows <- matrix(0,nrow = nrow(dist_matrix) - 1,ncol = 2)
  
  # minvalue means height in hc
  minValue <- rep(0,nrow(dist_matrix)-1)
  
  for (i in 1:(nrow(dist_matrix)-1))
  {
    # In order to find index minimum value in DistMatrix, it used which function.
    idx <- which(dist_matrix == min(dist_matrix,na.rm = TRUE), arr.ind = TRUE)
    
    minValue[i] <- min(dist_matrix,na.rm = TRUE)
    #browser()
    
    mPair <- tree_Group[idx[1,]]
    mergeData[i,] <- mPair[order(mPair)]
    #matrix mergerows shows where merges between attributes
    mergerows[i,] <- idx[1,]
    
    #In order to merge data bottom up, it bounds current pair and all previous groups they belong to
    # initially, tree_group compose of negative numbers but if column changed positive, it means that it has been clustered
    gIdx <- c(idx, which(tree_Group %in% tree_Group[idx[1,tree_Group[idx[1,]]>0]]))
    grouping <- c()
    
    if (length(gIdx) >0) {
      tree_Group[ c(idx[1,],gIdx) ] <- i
    }
    else {
      tree_Group[ c(idx[1,]) ] <- i
    }
    
    # To compute distance by each linkage, function 'apply' is used
    rDist <- apply(dist_matrix[idx[1,],],2,max)
    
    dist_matrix[min(idx),] <- rDist
    dist_matrix[,min(idx)] <- rDist
    
    # To move next the smallest distance, megered data is excluded so data is NA
    dist_matrix[min(idx),min(idx)] <- NA
    dist_matrix[max(idx),] <- NA
    dist_matrix[,max(idx)] <- NA
  }
  
  # making similar as function hclust in r, mergeData is converted to vector
  vc <- c()
  vc <- as.vector(abs(mergeData))
  vc <- unique(vc)
  
  # By using structure, this function returns merge data, height, and order
  structure (list(merge = mergeData, height = minValue, order = vc,labels = rownames(data)),class="hclust")
}

# pass wines_p to hc function
h=hc(wines_p)

# if pca is used uncomment the below row
## h=hc(ind.coord)

# if random data is used uncomment below
#h=hc(dataset)

# dendogram
fviz_dend(h, k = 3,                 
          cex = 0.5,                 
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800"),
          color_labels_by_k = TRUE,  
          ggtheme = theme_gray())
# phylogenic dendogram
fviz_dend(h, cex = 0.8, lwd = 0.8, k = 3, 
          rect = TRUE, 
          k_colors = "jco", 
          rect_border = "jco", 
          rect_fill = TRUE,
          type = "phylogenic",
          repel = TRUE)
fviz_dend(h, cex = 0.8, k=3, 
          rect = TRUE,  
          k_colors = "jco",
          rect_border = "jco", 
          rect_fill = TRUE, 
          horiz = TRUE)

dev.off()


