####################################
### Calculate AICs of OTU structures
### And then to plot them using the network package
###
### Jarrett Byrnes
### 4/5/2013
####################################

######
#Load relevant methods
######

source("./getNetworkAIC.R")
allInfo <- read.csv("./seq_otu_abundance.csv")

######
#Some work with AICs
######

#just to get a handle on things
LLMax <- getLogLik(allInfo, names(allInfo)[19], 20:27, mc.cores=2)
LLAMax <- getLogLik(allInfo, names(allInfo)[18], 20:27, mc.cores=2)
LLMin <- 0 #by definition
LLAMin <- getLogLik(allInfo, names(allInfo)[2], 20:27, mc.cores=6)

kmax <- length(unique(allInfo[,19]))
kamax <- length(unique(allInfo[,18]))
kamin <- length(unique(allInfo[,2]))
kmin <- length(unique(allInfo[,1]))


#Get the network LLs.  Note, the first one is 0 by definition
LLvec <-c(0, sapply(2:19, function(x) 
  getLogLik(allInfo, names(allInfo)[x], 20:27, binary=T, mc.cores=6)))

LLDvec <-c(0, sapply(2:19, function(x) 
  getLogLik(allInfo, names(allInfo)[x], 20:27, binary=F, network=F, mc.cores=6)))

#get the number og groups
kvec <- sapply(1:19, function(x) length(unique(allInfo[,x])))

#create a data frame with LLs
scores <- data.frame(group = names(allInfo)[1:19], K=kvec, LLnet = LLvec, LLdens = LLDvec)

#calculate AIC values
scores$AICnet <- with(scores, -2*LLnet + 2*K*6 + 2*nrow(allInfo))
scores$AICall <- with(scores, -2*LLnet - 2*LLdens+ 4*K*6 + 2*nrow(allInfo))
scores$group[which(scores$AIC==min(scores$AIC))]


write.csv(scores, "aic_scores.csv", row.names=F)

#####################
###graphics
#####################


netplot2(allInfo, names(allInfo)[12], 20:27) #lowest AICnet
netplot2(allInfo, names(allInfo)[9], 20:27) #lowest AICall

netPlot(allInfo, names(allInfo)[12], 20:27,high.spacing=0.0009)