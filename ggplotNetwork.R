library(reshape2)
library(plyr)
library(ggplot2)

allInfo <- read.csv("./seq_otu_abundance.csv")

ggnetplot <- function(adf=allInfo, otu=names(allInfo)[12], ids=20:27, colorvar="Plot"){

  meltdf <- melt(adf, id.vars=otu,
               measure.vars=names(allInfo)[20:27],
               variable.name="Plot",
               value.name="Abundance")
  names(meltdf)[1] <- "Cluster"

  dat <- ddply(meltdf, .(Plot, otu), Abundance=sum(Abundance))
  dat <- dat[which(dat$Abundance > 0),]
  levels(dat$Plot) <- gsub("X", "", levels(dat$Plot))
  
#  dat$Plot <- factor(dat$Plot, levels=sort(levels(dat$Plot))) #deal with bad sorting
  dat$Plot <- factor(dat$Plot, levels=c("3c", "7c", "1l", "5l", "2h", "9h", "6x", "8x")) #deal with bad sorting

  dat$xtop <- as.numeric(dat$Cluster)
  dat$xtop <- dat$xtop/max(dat$xtop)*8
  dat$xbottom <- as.numeric(dat$Plot)
  dat$y <- rep(1, nrow(dat))
  dat$yend <- rep(0, nrow(dat))
  
  if(colorvar=="Cluster") {
  	netplot <- ggplot(data=dat, 
                  aes(x=xtop,
                      xend=xbottom, 
                      y=y, 
                      yend=yend,
                      alpha=Abundance, color=Cluster)) 
    }
    
   if(colorvar=="Plot") {
  	netplot <- ggplot(data=dat, 
                  aes(x=xtop,
                      xend=xbottom, 
                      y=y, 
                      yend=yend,
                      alpha=Abundance, color=Plot)) 
    }
                      
  netplot <- netplot + 
  	geom_segment() + 
  	theme_bw() + 
  	ylab("") + 
  	xlab("") +
  	annotate("text", x=1:8, y=-0.012, label=levels(dat$Plot)) + 
  	theme(axis.ticks=element_blank(), axis.text=element_blank())

 return(netplot)
}
 
###Plot one 
 a<- ggnetplot(otu=names(allInfo)[12])
 a + ggtitle(names(allInfo)[12])
 ggsave(names(allInfo)[19])

###Color by OTU
 b<- ggnetplot(otu=names(allInfo)[12], colorvar="Cluster")
 b + ggtitle(names(allInfo)[12])


#### ALL OF THEM - note, do not run unless you want to lock up your computer!
sapply(names(allInfo)[1:19], function(x){
  a<- ggnetplot(otu=x)+ ggtitle(x)
  ggsave(paste("./networks/", x, "_bi_plot.jpg", sep=""), a)
  
  
})


sapply(names(allInfo)[1:19], function(x){
  a<- ggnetplot(otu=x, colorvar="Cluster") + ggtitle(x)
  ggsave(paste("./networks/", x, "_bi_cluster.jpg", sep=""), a)
  
})