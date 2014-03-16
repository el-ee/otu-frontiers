# otu-frontiers

Code for the analysis of sequence and OTU occurance data to derive optimal taxonomic structures for bacterial communities.

See original if you want to do anything with this.

I am just messing it up right now, as I struggle to figure out what is happening!

Paper describing technique: 

* http://journal.frontiersin.org/Journal/10.3389/fmicb.2013.00342/full  


# THE FILES

FunFrame Pipeline used to get files to the right starting point: 

* http://faculty.www.umb.edu/jennifer.bowen/software/FunFrame.zip
* See http://bioinformatics.oxfordjournals.org/content/29/9/1212.long
 

## Data organization


## 1. First use this file, which reads in this fasta data by traversing directories & produces otu csv

### getOTUs.R

Puts together multiple .fasta files from different directories??

* INPUT: Multiple directories each of which contain a `fasta` formatted file; specify common file name in variable `FILE_NAME`. THat is, for the example above, you would set `FILE_NAME=data.fasta`

* OUTPUT: `OTU_MAT.csv` This file is then used by the other files


## 2. After that, you can use this one to get abundance data?

### getAbundance.R

? Don't really understand which file here

* INPUT: `./clusterDist_0.12/concordance.fasta`
* INPUT: `OTU_MAT.csv`

* OUTPUT: `seq_otu_abundance.csv`


## Then these are things you can do with that csv file about abundance

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
