# otu-frontiers

Code for the analysis of sequence and OTU occurance data to derive optimal taxonomic structures for bacterial communities.

# THE FILES

Trying to make some sense of how to do anything with this code, by organizing the files.

Best guess at arranging them in some order of how they might be **used**


## Data organization
Best I can tell, data should be organized in a series of directories, one directory for each sample, fasta file for that sample inside its directory

like

```
/S01/data.fasta
/S02/data.fasta
/S03/data.fasta
...
/Sxx/data.fasta
```

You should run these R files from the parent directory that contains all of the S* directories. 

## 1. First use this file, which reads in this fasta data by traversing directories & produces otu csv

### getOTUs.R

Puts together multiple .fasta files from different directories

* INPUT: Multiple directories each of which contain a `fasta` formatted file; specify common file name in variable `FILE_NAME`. THat is, for the example above, you would set `FILE_NAME=data.fasta`

* OUTPUT: `OTU_MAT.csv` This file is then used by the other files


## 2. After that, you can use this one to get abundance data?

### getAbundance.R

? Don't understand what is special about the one file this one takes in. 

? Can it just be any one of the random data files?

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