#############################################
# Code to bring together abundances of different sequences 
# into a single data file and merge them with OTU information
#
# Jarrett Byrnes
# 3/28/13
#############################################

library(reshape2)

#A function to take a read in fasta file that returns a 
#data frame with sequences and OTU memberships
get_plotseq <- function (fastaSet){
  #use mclapply to scroll over all lines to make this shorter
  #in order to later transition to mclapply for speed
  densityTab <- mclapply(1:length(fastaSet), function(x) {
    #get OTU membership for a line
    b <- strsplit(names(fastaSet[x]), ":")[[1]]
    pl <- strsplit(b[2], "\\|")[[1]][1]
    cnt <- strsplit(b[3], "\\|")[[1]][1]
    class(cnt) <-"numeric"
    data.frame(Sequence=as.character(fastaSet[[x]]), pl=pl, Count = cnt)
  }, mc.cores=2)
  
  densityTab <- ldply(densityTab)
  #densityTab$Count <- as.numeric(as.character(densityTab$Count))
  
  densityTab <- dcast(densityTab, Sequence ~ pl, value.var="Count", 
                      fun.aggregate=sum,
                      fill=0)
  
  
  densityTab
}

#use one site to get the abundances for all sequences
aSite <-  readDNAStringSet("./clusterDist_0.12/concordance.fasta", "fasta")
sequence_abundance <- get_plotseq(aSite)
write.csv(sequence_abundance, "sequence_abundance.csv")

#merge this with info about the otu-sequence assignments
OTU_MAT <- read.csv("OTU_MAT.csv")
allInfo <- merge(OTU_MAT, sequence_abundance)
write.csv(allInfo, "seq_otu_abundance.csv", row.names=F)
