rm(list = ls()) # Clean variable
memory.limit(150000)

##### Current path and new folder setting* #####
ProjectName = "KMean"
Sampletype = "Iris"
#ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
Save.Path = paste0(getwd(),"/",Version)
## Create new folder
if (!dir.exists(Save.Path)){
  dir.create(Save.Path)
}

##### Load Packages #####
Package.set <- c("tidyverse","ggpubr","factoextra")
## Check whether the installation of those packages is required from basic
for (i in 1:length(Package.set)) {
  if (!requireNamespace(Package.set[i], quietly = TRUE)){
    install.packages(Package.set[i])
  }
}
## Load Packages
lapply(Package.set, library, character.only = TRUE)
rm(Package.set,i)


##### Data preparation #####
data("iris")
df <- iris
head(df, 3)

##### K-means clustering calculation example #####
# Compute k-means with k = 3
set.seed(123)
res.km <- kmeans(scale(df[, -5]), 3, nstart = 25)
# K-means clusters showing the group of each individuals
res.km$cluster

##### Plot k-means #####
##********************** Using the factoextra R package **********************##

fviz_cluster(res.km, data = df[, -5],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
) -> p1
p1

##********************** Using the ggpubr R package **********************##
##### Compute PCA and extract individual coordinates #####
# Dimension reduction using PCA
res.pca <- prcomp(df[, -5],  scale = TRUE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(res.km$cluster)
# Add Species groups from the original data sett
ind.coord$Species <- df$Species
# Data inspection
head(ind.coord)


# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)


##### Visualize k-means clusters #####
ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Species", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4) -> p2
p2

##### Export PDF #####
pdf(
  file = paste0(getwd(), "/",Version,"/", Sys.Date(), "_KMeans.pdf"),
  width = 7, height = 7
)
p1
p2
graphics.off()



