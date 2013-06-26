####################################
### Plot OTU Networks using the network package
###
### Jarrett Byrnes
### 4/5/2013
####################################

######
#Load relevant methods
######

source("./getNetworkAIC.R")
allInfo <- read.csv("./seq_otu_abundance.csv")

#####################
###graphics
#####################

library(RColorBrewer)

#siteCols = rainbow(9)[-1]
siteCols <- brewer.pal(8, "RdYlGn")

netplot2(allInfo, names(allInfo)[12], 20:27,
         siteCol=siteCols,
         main=names(allInfo[12])) #lowest AICnet, exclude red from site colors

legend(5,10, "Hello!")

clustering <- names(allInfo)[12]#lowest AICall
#clustering <- names(allInfo)[2]#lowest AICall

#VIEW BY NUMBER OF PLOTS FOR EACH OTU
set.seed(2002)
netplot2(allInfo, clustering, 20:27,
         siteCol = rep("black", 8),
         main=clustering,
         site.cex=3, 
         setCol="occurance",
         edge.col="grey",
         occur.cols = siteCols)

legend(-130, -128, fill=siteCols, legend=paste(1:8, "plots", sep=" "))

#VIEW BY TOTAL ABUNDANCE FOR EACH OTU
set.seed(197)
siteCols <- brewer.pal(8, "Spectral")
netplot2(allInfo, clustering, 20:27,
         siteCol = siteCols,
         main=clustering,
         site.cex=4, 
         edge.col="grey",
         seqCol="black",
         setSize ="abundance",
         size.scale=3,
         site.sides=4) 

legend(x=- 66, y= -42, fill=siteCols, legend=gsub("X", "", names(allInfo)[20:27]))

#need max size
maxAbund <- max(rowSums(getOTUAbundMat(allInfo, clustering)[,2:9]))
legend(x= 123, y= -42, col="black", pch=19, legend=round(seq(1,maxAbund,length.out=8)), pt.cex=seq(0.1, 3, length.out=8))


#see who is who
plot(1:8, rep(1,8), col=rainbow(9)[-1], pch=19)
text(1:8, rep(1,8)+0.1, labels=names(allInfo)[20:27])

netPlot(allInfo, names(allInfo)[12], 20:27,high.spacing=0.0009)

#####################
###what about an adjacency matrix using connections?
#####################

head(allInfo)
otu_12 <- levels(allInfo$clusterDist_0.12)

getEdgeWeight <- function(x,y, acol="clusterDist_0.12") {
  otu1 <- colSums(allInfo[which(allInfo[[acol]] ==x),20:27])>0
  otu2 <- colSums(allInfo[which(allInfo[[acol]] ==y),20:27])>0
  sum(colSums(rbind(otu1, otu2))>1)
}

# #getEdgeWeight("OTU111", "OTU111")
# getEdgeWeight <- Vectorize(getEdgeWeight, c("x", "y"))

# #outer(otu_12[1:20], otu_12[1:20], FUN="getEdgeWeight")
# site_adjmat <- outer(otu_12, otu_12, FUN="getEdgeWeight")
# rownames(site_adjmat) <- colnames(site_adjmat) <- otu_12

# write.csv(site_adjmat, "./site_adjmat_12.csv", row.names=T)

site_adjmat<-read.csv("./site_adjmat_12.csv", row.names=1)
site_network <- network(as.matrix(site_adjmat), directed=F, loops=F)
plot(site_network, edge.lwd=site_adjmat)