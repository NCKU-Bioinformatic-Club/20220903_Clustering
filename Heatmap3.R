##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


# install.packages("heatmap3")
if(!require("heatmap3")) install.packages("heatmap3")

# Load heatmap3 package
library(heatmap3)

# View help files
?heatmap3

# The examples of heatmap3 and other functions
example(heatmap3)
example(showLegend)
example(showAnn)
example(colByValue)


##### Loading data #####
# wget https://github.com/slzhao/heatmap3/releases/download/example/allSample_edgeR_result.csv
# wget https://github.com/slzhao/heatmap3/releases/download/example/BRCA_30Samples.csv
# wget https://github.com/slzhao/heatmap3/releases/download/example/BRCA_30Samples_clinic.csv

#Prepare expression data
counts <- read.csv("https://github.com/slzhao/heatmap3/releases/download/example/BRCA_30Samples.csv",header=T,row.names=1)
#Prepare column side annotation
clinic <- read.csv("https://github.com/slzhao/heatmap3/releases/download/example/BRCA_30Samples_clinic.csv",header=T,row.names=1)
#Prepare row side color bar annotation
edgeR_result <- read.csv("https://github.com/slzhao/heatmap3/releases/download/example/allSample_edgeR_result.csv",header=T,row.names=1)

temp1 <- (edgeR_result$logFC)
temp2 <- -log10(edgeR_result$FDR)
temp1 <- colByValue(as.matrix(temp1),range=c(-4,4),col=colorRampPalette(c('chartreuse4','white','firebrick'))(1024))
temp2 <- colByValue(as.matrix(temp2),range=c(0,5),col=heat.colors(1024))
colGene <- cbind(temp1,temp2)
row.names(colGene)<-row.names(edgeR_result)
colnames(colGene)<-c("log2FC","-Log10P")

#Generate Figure1
#counts, colGene and clinic were read throught the csv file
##Assume counts has counts information for each gene, colGene has the colors for each gene, clinic has the clinic information for each sample
temp<-apply(counts,1,sd)
selectedGenes<-rev(order(temp))[1:500]
heatmap3(counts[selectedGenes,],labRow="",margin=c(7,0),RowSideColors=colGene[selectedGenes,],ColSideCut=0.85,ColSideAnn=clinic,ColSideFun=function(x) showAnn(x),ColSideWidth=1.2,balanceColor=T)

#Generate Figure2
heatmap3(counts,topN=c(500,3000,nrow(counts)),Rowv=NA,labRow="",margin=c(7,0),RowSideColors=colGene,ColSideCut=0.85,ColSideAnn=clinic,ColSideFun=function(x) showAnn(x),ColSideWidth=1.2,balanceColor=T)
