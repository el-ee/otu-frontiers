####################################
### Methods to calculate LLs of
### OTU structures
### And then to plot them using the network package
###
### Jarrett Byrnes
### Last updated: 4/11/2013
###
### Changelog
### 4/11/2013 - added poisson LL for LikColDensity
####################################


library(bipartite)
library(network)
library(plyr)
library(parallel)

## Get the LL based on all sequences being assumed to have
## equal abundances within a group
LikColDensity <- function(acol, useLogLik=T, dist="binom"){
  if(sum(acol)==0) return(0) #defining 0*log(0) as 0 as in Allesina and Pascual 2009
  s <- sum(acol)
  if(dist=="binom") LL <- dbinom(acol, size=s, prob=1/length(acol), log=useLogLik)
  if(dist=="pois") LL <- dpois(acol, s/length(acol), log=useLogLik)
  if(useLogLik) return(sum(LL))
  prod(LL)
}  

## Get the LL based on network structure
LikColNet <- function(acol, useLogLik=T){
  if(sum(acol)==0) return(0) #defining 0*log(0) as 0 as in Allesina and Pascual 2009
  acol <- as.numeric(acol>0)
  s <- length(acol)
  l <- sum(acol)
  if(s==l) return(0) #perfect match, and if we're going with logLik, we'll get log(0) problems
  
  p <- l/(s)
  ret <- l*log(p) + (s-l)*log(1-p) 
  if(!useLogLik) ret <- p^l*(1-p)^(s-l)
  return(ret)
}


## Take a OTU column and calculate the LL of the grouping structure
getLogLik <- function(adf, otuCol, ids=4:6, getNet=T, binary=T, mc.cores=2, dist="binom"){
  otuCol <-adf[[otuCol]]
  liks <- mclapply(unique(otuCol), function(x){
    
    reducedDF <- adf[which(otuCol==x),]
    LL <- 0
    if(getNet) LL <- LL + sum(apply(reducedDF[,ids], 2, LikColNet))
    if(!binary) LL <- LL + sum(apply(reducedDF[,ids], 2, LikColDensity, dist=dist))
    
    LL
      
  }, mc.cores=mc.cores)
  
  liks <- as.vector(simplify2array(liks))
  
  sum(liks)
  
}

# loc <- matrix(c(1,0,0,
                # 1,0,0,
                # 0,1,1,
                # 0,1,1,
                # 1,0,1), ncol=3, byrow=T)

# names(loc) = c("one", "two", "three")
# mdf <- data.frame(low=1:5, med = c(1,1,2,2,3), high=c(1,1,2,2,2))

# mdf <- cbind(mdf, loc)

### Uses the bipartite library to make a bipirtate graph
### although I'm note a fan of how it looks
netPlot <- function(adf, otuCol, ids=4:6, ...){
  groupWeb <- ddply(adf, otuCol, function(x) colSums(x[,ids]))
  rownames(groupWeb) <- groupWeb[,1]
  groupWeb <- groupWeb[,-1]
  plotweb(t(groupWeb), high.lablength=0,...)
}



#netPlot(mdf, "low")
#netPlot(mdf, "med")
#netPlot(mdf, "high")


### Uses the network library to make a network graph
### with the sites as the blue nodes.  Can accept
### other arguments to plot.network
getOTUAbundMat <- function(adf, otuCol, ids=20:27) ddply(adf, otuCol, function(x) colSums(x[,ids]))

netplot2 <- function(adf, otuCol, ids=4:6, edge.scale=20, 
                     setCol="site",
                     setSize="same",
                     seqCol = "red",
                     siteCol = rep("grey", length(ids)),
                     site.cex=1,
                      edge.col="black",
                     size.scale=2,
                     site.sides=50,
                     occur.cols=rainbow(length(ids)),
                     ...){
  groupWeb <- getOTUAbundMat(adf, otuCol)
  groupWeb <- groupWeb[,-1]
  zmat1 <- matrix(rep(0, length(ids)^2), ncol=length(ids))
  g <- rbind(zmat1, as.matrix(groupWeb))
  zempty <-  matrix(rep(0, nrow(g)*nrow(groupWeb)), nrow=nrow(g))
  g <- cbind(g, zempty)
  n_plots <- rowSums(g>0)
  abund <- rowSums(g)
  g2 <- network(as.matrix(g))
  
  vertex.col <- c(siteCol, rep(seqCol, nrow(groupWeb)))
  vertex.cex <- c(rep(site.cex, length(ids)), rep(1, nrow(groupWeb))) 
  vertex.sides <- c(rep(site.sides,length(ids)), rep(50, nrow(groupWeb)))
  
  if(setCol=="occurance"){
    
    vertex.col <- c(siteCol, occur.cols[n_plots])
    vertex.cex <- c(rep(site.cex, length(ids)), rep(1, nrow(groupWeb)))
    vertex.sides <- c(3:(3+length(ids)), rep(50, nrow(groupWeb)))
    
    
  }
  
  if(setSize == "abundance") vertex.cex <- c(rep(site.cex, length(ids)), size.scale*abund/max(abund))
    
  #note, I still don't like how edges work here...need to figure out a better scheme
  plot(g2, vertex.col=vertex.col,
       vertex.cex = vertex.cex,
       vertex.sides = vertex.sides,
      # edge.lwd=g/max(g)*edge.scale, 
       usearrows=FALSE,
       edge.col=edge.col,
       ...)
}

