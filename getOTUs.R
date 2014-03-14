#############################################
# Code to bring together OTUs from different levels of clustering
# into a single data file
#
# Jarrett Byrnes
# 3/28/13
#############################################

library(Biostrings)
library(plyr)
library(multicore)

#a <- readDNAStringSet("./clusterDist_0.12/concordance.fasta", "fasta")

#sq <- as.character(a[[1]])

#the common file name for all .fasta files in your directories
FILE_NAME = "R1.fasta"

#A function to take a read in fasta file that returns a 
#data frame with sequences and OTU memberships
get_otuseq <- function (fastaSet){
  #use mclapply to scroll over all lines to make this shorter
  #in order to later transition to mclapply for speed
  SOTU <- mclapply(1:length(fastaSet), function(x) {
    #get OTU membership for a line
    b <- strsplit(names(fastaSet[x]), ":")[[1]]
    otu<-b[length(b)]
  
    data.frame(S=as.character(fastaSet[[x]]), OTU=otu)
  }, mc.cores=2)

  SOTU <- ldply(SOTU)

  names(SOTU) <- c("S", "OTU")
  SOTU$S <- as.character(SOTU$S)

  SOTU
}

#What are the directories
dirs <- dir()

# Lost on what is happening with this next line, exactly; previous command by itself fetches directories in current folder which is how i'm runnig this, possibly incorrectly? but...  
# dirs <- dirs[-grep("\\.R", dirs)]

###Now let's figure out how many unique sequences we have
# E: added separator as /
baseline <- readDNAStringSet(paste(dirs[1], FILE_NAME, sep="/"), "fasta")
#which are the unique sequences across all files
uqIdx <- which(!duplicated(as.vector(baseline))) 
seq_rows <- length(uqIdx)

concordance_array <- array(rep(NA,2*seq_rows*length(dirs)), 
                           c(seq_rows,#200,#
                             2,
                             length(dirs)))

# E: added separator as /
for(i in 1:length(dirs)){
  a_concordance <- readDNAStringSet(paste(dirs[i], FILE_NAME, sep="/"), "fasta")
  uqIdxI <- which(!duplicated(as.vector(a_concordance)))#[1:200] 
  concordance_array[,,i] <- as.matrix(get_otuseq(a_concordance[uqIdxI]) )  
}



#get a matrix of rows which match all sequences to the sequences seen in the first file
concordance_array[,1,1] -> base
match(base, concordance_array[,1,2])

#get the row numbers of matching sequences in the concordance array
idxMat <- apply(concordance_array, 3, function(x) match(base, x))
idxMat <- cbind(1:length(base), idxMat)


#now put together a matrix of OTU identities...
OTU_MAT <- sapply(1:length(dirs),  function(idx){
  concordance_array[idxMat[,idx],2,idx]
})


OTU_MAT<-cbind(base, OTU_MAT)
OTU_MAT<-as.data.frame(OTU_MAT)
names(OTU_MAT) <- c("Sequence", dirs)

write.csv(OTU_MAT, "OTU_MAT.csv", row.names=F)

