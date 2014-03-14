# otu-frontiers

Code for the analysis of sequence and OTU occurance data to derive optimal taxonomic structures for bacterial communities.

# THE FILES

Trying to make some sense of how to do anything with this code, by organizing the files.

Best guess at arranging them in some order of how they might be **used**


## First this one, reads in fasta, produces otu csv

### getOTUs.R

Puts together multiple .fasta files from different directories

* INPUT: Multiple directories each of which contain a file, `concordance.fasta`

* OUTPUT: `OTU_MAT.csv`


## Second this one

### getAbundance.R

Don't understand what is special about the one file this one takes in. 

* INPUT: `./clusterDist_0.12/concordance.fasta`
* INPUT: `OTU_MAT.csv`

* OUTPUT: `seq_otu_abundance.csv`


## Then these are things you can do with that csv file

### analyze_otu_networks.R

* INPUT: `./seq_otu_abundance.csv`

* OUTPUT: 

### plotOTUNetworks.R

* INPUT: `/seq_otu_abundance.csv`

* OUTPUT: a plot


### ggplotNetwork.R

* INPUT: `./seq_otu_abundance.csv`

* OUTPUT: a plot



## Other files

### getNetworkAIC.R

Helper file; used by `analyze_otu_networks.R` and `sampleAICNet.R`



### sampleAICNet.R

Confused.