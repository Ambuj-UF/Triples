#Triples
Calculates resolution of triples from the set of phylogeny trees

###Requires package ape

##Installation instruction
1. In R command line run install.packages(‘ape’) to install ape library
2. Download the zip and install by running 

```
R CMD INSTALL Triples
```

You can use datasets available in data directory to test this package.

##Loading Triples library

After installation type following in R command line

```
library(Triples)
triples("prev.tre", "authority.txt", outgroup = "Dasypus_novemcinctus", output="output.txt")
```