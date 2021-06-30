# Running basic EGAD
```
load("bin/run_GBA.Rdata")
load("data/GO.human.Rdata")
source("bin/helper_functions.r")

exprs_file = "data/GSE12946_expression_FPKM.parsed"
exprs = read.table(file=exprs_file, header=F)
geneids = exprs[,1]
exprs = exprs[,-1]
exprs = as.matrix(exprs) 
rownames(exprs) = geneids
```
 
### 0. Setup and get expression experiment 
``` 
load("bin/run_GBA.Rdata")
load("data/GO.human.Rdata")
source("bin/helper_functions.r")

exprs_file = "data/GSE12946_expression_FPKM.parsed"
exprs = read.table(file=exprs_file, header=F)
geneids = exprs[,1]
exprs = exprs[,-1]
exprs = as.matrix(exprs) 
rownames(exprs) = geneids
``` 

### 1. Calculate network edges using correlations and rank standardize
Note, spearman correlations are recommended. 
``` 
network = EGAD::build_coexp_network(exprs[1:1000,]) 
```
 
### 2. Hard threshold for top connections  
Hard threshold by taking the top 1% of connections:  
```
network_top = threshold_network_top_genes(network, 0.01) # code in helper_functions.r
```
### 3. Binary network using only connections >0.95 
```
network_bin = network[network>0.95] * 1  
```
  
### 4. Testing for network functionality 
``` 
GBA = run_GBA(network, GO.labels)
```

### 5. Testing for network assortativity



