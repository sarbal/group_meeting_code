# Notes 

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

### 5. Checking node properties through node degrees 
```
nd_net <- node_degree(network)
plot_distribution(nd_net, xlab="Node degrees")
``` 

### 6. Testing for network assortativity
```
A_net <- assortativity(network)
``` 

### 6. Testing for multifunctionality of annotations 
```
multifunc_assessment <- calculate_multifunc(annots) 
auc_mf <- auc_multifunc(annots, multifunc_assessment[,4]) 
hist <- plot_distribution(auc_mf, xlab="AUROC", med=FALSE, avg=FALSE)
```

### 7. Graphing networks using igraph
```
sub_net <- network[1:100,1:100] # just taking a subset again

## Removing self-connections
diag(sub_net) <-  0

## Getting the top triangle 
upper <- row(sub_net) < col(sub_net)
pairs <- which(upper, arr.ind = T )

## Getting labels and weights 
gene_names <- rownames(sub_net)
weights <- sub_net[upper]
pairs <- data.frame( p1 = gene_names[pairs[,1]], p2 = gene_names[pairs[,2]] , weights = weights )

## Build graph with igraph  
inet <- igraph::graph_from_data_frame(pairs, directed=F)

### Add weights to edges, and play around with other properties 
igraph::E(inet)$weight <- weights
igraph::E(inet)$width <- (weights^2 * 10)
igraph::E(inet)$edge.color <- "gray80"
igraph::E(inet)$color <- viridis(101)[ round(weights * 100) + 1  ]

o <- match(igraph::V(inet)$name, clust_net$clusters$genes)

### Add colors and remove labels from vertices (nodes) 
igraph::V(inet)$color <- as.character(clust_net$clusters$colors)[o]
igraph::V(inet)$label <- ""
igraph::V(inet)$size <- 4

### Tidy up network to only show weights of given threshold  
inet_sub <-  igraph::delete_edges(inet, igraph::E(inet)[weights > threshold])

plot(inet_sub )
```
