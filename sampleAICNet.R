####################################
### Demo of calculating LLs of
### group structures
### And plotting them using the network package
###
### Jarrett Byrnes
### Last updated: 4/11/2013
####################################


loc <- matrix(c(50,0,0,
                45,0,0,
                0,100,1,
                0,112,7,
                0,12,110), ncol=3, byrow=T)

names(loc) = c("one", "two", "three")
mdf <- data.frame(low=1:5, med = c(1,1,2,2,3), high=c(1,1,2,2,2))

mdf <- cbind(mdf, loc)

par(mfrow=c(1,3))
netplot2(mdf, "low", vertex.cex=5, main="Separate");box()
netplot2(mdf, "med", vertex.cex=5, main="Moderate Grouping");box()
netplot2(mdf, "high", vertex.cex=5, main="Highly Grouped");box()
par(mfrow=c(1,1))

aicdf <- data.frame(
	k=c(5,3,2),
	LLNet = sapply(names(mdf)[1:3], function(x) getLogLik(mdf, x, ids=4:6)),
	LLBinomNet = sapply(names(mdf)[1:3], function(x) getLogLik(mdf, x, ids=4:6, binary=F, dist="binom")),
	LLPoisNet = sapply(names(mdf)[1:3], function(x) getLogLik(mdf, x, ids=4:6, binary=F, dist="pois")))

aicdf <- within(aicdf, {
	AICnet <- -2*LLNet+2*k*3
	AICbinom<- -2*LLBinomNet+4*k*3
	AICpois <- -2*LLPoisNet+2*k*3
	
})

aicdf
