####################################
### Methods to calculate LLs of
### OTU structures
### And then to plot them using the network package
###
### Jarrett Byrnes
### 4/5/2013
####################################


library(bipartite)
library(network)
library(plyr)
library(parallel)

## Get the LL based on all sequences being assumed to have
## equal abundances within a group
LikColDensity <- function(acol, useLogLik=T){
  if(sum(acol)==0) return(0) #defining 0*log(0) as 0 as in Allesina and Pascual 2009
  s <- sum(acol)
  LL <- dbinom(acol, size=s, prob=1/length(acol), log=useLogLik)
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
getLogLik <- function(adf, otuCol, ids=4:6, getNet=T, binary=T, mc.cores=2){
  otuCol <-adf[[otuCol]]
  liks <- mclapply(unique(otuCol), function(x){
    
    reducedDF <- adf[which(otuCol==x),]
    LL <- 0
    if(getNet) LL <- LL + sum(apply(reducedDF[,ids], 2, LikColNet))
    if(!binary) LL <- LL + sum(apply(reducedDF[,ids], 2, LikColDensity))
    
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


### Uses the network library to make a network graph
### with the sites as the blue nodes.  Can accept
### other arguments to plot.network
netplot2 <- function(adf, otuCol, ids=4:6, edge.scale=20, ...){
  groupWeb <- ddply(adf, otuCol, function(x) colSums(x[,ids]))
  groupWeb <- groupWeb[,-1]
  zmat1 <- matrix(rep(0, length(ids)^2), ncol=length(ids))
  g <- rbind(zmat1, as.matrix(groupWeb))
  zempty <-  matrix(rep(0, nrow(g)*nrow(groupWeb)), nrow=nrow(g))
  g <- cbind(g, zempty)
  
  g2 <- network(as.matrix(g))
  
  #note, I still don't like how edges work here...need to figure out a better scheme
  plot(g2, vertex.col=c(rep("blue", length(ids)), rep("red", nrow(groupWeb))),
       edge.lwd=g/max(g)*edge.scale, ...)
}


#netPlot(mdf, "low")
#netPlot(mdf, "med")
#netPlot(mdf, "high")
