## Other ways to "build" gene-gene networks to incorporate additional gene-gene expression relationships
Standard: 
```
network = EGAD::build_coexp_network(exprs[1:1000,])
```
Including negative correlations (take the absolute values of the correlations)
```
networkabs  = EGAD::build_coexp_network(exprs[1:1000,], flag="absrank") 
```
Including non-linear relationships using mutual information
   * https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-328   
   * https://en.wikipedia.org/wiki/Mutual_information  
```
require(parmigene)
mi <- parmigene::knnmi.all(exprs[1:1000,]) 
# Rank the MI data -> makes this similar to the networks obtained above (rank standardizes) 
n = dim(mi)[1]
mi.rank = matrix(rank(mi, na.last="keep",ties.method="average"), nrow=n, ncol=n)
mi.rank = mi.rank/max(mi.rank, na.rm=T)
```

### Using node degrees to compare 
```
nd_mi = node_degree(mi.rank)
nd_co = node_degree(network)
nd_coa = node_degree(networkabs)
```
Scatter plots of the node degrees (standardized, so that the node degree is between 0 and 1)
```
### Plot between the MI and Coexp
plot(nd_mi/n, nd_co/n, pch=19)
### Plot between Coexp and absolute coexp 
plot(nd_mi/n, nd_coa/n, pch=19)
### Plots of the distributions of node degrees 
hco = hist(nd_co/n, breaks=c(10:80)/100 , col=make_transparent("orange"), border=NA) 
hist(nd_mi/n, breaks=hco$breaks, add=T, col=make_transparent("blue"), border=NA)
hist(nd_coa/n, breaks=hco$breaks, add=T, col=make_transparent("green"), border=NA)
```

Most uses for MI nets are to also build GRNs -> gene regulatory networks. These again try to search for the most likely "real" connections. There are a variety of algorithms that work with this, I've selected a few, and you can play around with their parameters. For the most part, coexpression is similar enough that simpler methods usually outperform them ie. like just thresholding for strong connections. More of this in my 2015 paper if you're interested. Some of the ones I tested are below. 
```
## Run ARACNE-a (Algorithm for the Reconstruction of Accurate Cellular Networks, additive model)
grn.a.a <- aracne.a(mi, 0.05) 

## Run ARACNE-m (Algorithm for the Reconstruction of Accurate Cellular Networks, multiplicative model)
grn.a.m <- aracne.m(mi,tau=0.15) 

## Run CLR (Context Likelihood or Relatedness Network)
grn.clr <- clr(mi) 

## Run MRNET (Maximum Relevance Minimum Redundancy)
grn.mrnet <- mrnet(mi)
```

Tidying up the networks to be able to compare them to the others. 
```
diag(grn.a.m) = max(grn.a.m, na.rm=T) 
diag(grn.a.a) = max(grn.a.a, na.rm=T)
diag(grn.clr) = max(grn.clr, na.rm=T) 
diag(grn.mrnet) = max(grn.mrnet, na.rm=T) 

grn.a.a.rank = matrix(rank(grn.a.a, na.last="keep",ties.method="average"), nrow=n, ncol=n)
grn.a.a.rank = grn.a.a.rank/max(grn.a.a.rank, na.rm=T)
# heatmap.3(grn.a.a.rank)
grn.a.m.rank = matrix(rank(grn.a.m, na.last="keep",ties.method="average"), nrow=n, ncol=n)
grn.a.m.rank = grn.a.m.rank/max(grn.a.m.rank, na.rm=T)
#heatmap.3(grn.a.m.rank)
grn.clr.rank = matrix(rank(grn.clr, na.last="keep",ties.method="average"), nrow=n, ncol=n)
grn.clr.rank = grn.clr.rank/max(grn.clr.rank, na.rm=T)
#heatmap.3(grn.clr.rank)
grn.mrnet.rank = matrix(rank(grn.mrnet, na.last="keep",ties.method="average"), nrow=n, ncol=n)
grn.mrnet.rank = grn.mrnet.rank/max(grn.mrnet.rank, na.rm=T)
# heatmap.3(grn.mrnet.rank)
```
Calculating node degrees
```
nd_grna = node_degree(grn.a.a.rank)
nd_grnm = node_degree(grn.a.m.rank)
nd_grnc = node_degree(grn.clr.rank)
nd_grnn = node_degree(grn.mrnet.rank)
```
Plots 
```
plot(nd_mi/n, nd_co/n, pch=19)
plot(nd_grna/n, nd_co/n, pch=19)
plot(nd_grnm/n, nd_co/n, pch=19)
plot(nd_grnc/n, nd_coa/n, pch=19)
plot(nd_grna/n, nd_coa/n, pch=19)
plot(nd_grnm/n, nd_coa/n, pch=19)
plot(nd_grnn/n, nd_coa/n, pch=19)
plot(nd_grnn/n, nd_mi/n, pch=19)
```

 
## Clustering and plotting 
```
## Calculate the euclidean distances between genes using their coexpression values
dist_net <- dist(network) 
## Cluster (defaults)
tree_net <- hclust( dist_net, method="average" )
## Extract dendrogram
dend_net <- as.dendrogram(tree_net)
## Get "order" of genes in dendrogram
colInd <- order.dendrogram(dend_net)

## Plot "heatmap" based on this ordering
image( network[colInd,colInd], col=heat.colors(100))
## Plot dendrogram
plot(dend_net)
### There are lots of other ways to visualise this with base R and other libs: 
### library(phylogram) -> https://cran.r-project.org/web/packages/phylogram/vignettes/phylogram-vignette.html 
### library(ape) -> https://rdrr.io/cran/ape/man/plot.phylo.html
```

Some other example plots using the above libraries 
```
require(ape)
require(phylogram)

plot(as.phylo(tree_net), type="radial") #  "phylogram", "cladogram", "fan", "unrooted", "radial" 
plot(as.phylo(tree_net), type="unrooted")
plot(as.phylo(tree_net), type="fan")
plot(as.phylo(tree_net), type="cladogram")

plot(as.phylo(tree_net), type="fan", tip.col = cutree( tree_net, h=10) )

plot(as.phylo(tree_net), type="unrooted", tip.col = cutree( tree_net, k=10) )

# Color the clusters from the cut tree function 
cluster_cols = plasma(15)[cutree( tree_net, k=10)]

plot(as.phylo(tree_net), type="unrooted", tip.col =  cluster_cols)
```

Heatmap plots, using the above clustering.
```
col_map = viridis(100) 
heatmap.2( network, density.info = "none", trace = "none",
           col = col_map,
           Rowv = dend_net, Colv = dend_net, 
           RowSideColors = cluster_cols,
           ColSideColors =  cluster_cols)

nd_cols = magma(101)[ round(nd_co/n,2)*100 ]
```

Heatmap plots, using the above clustering. Color clusters using other properties, eg node degrees. And color rows diff to columns. 
```
heatmap.2( network, density.info = "none", trace = "none",
           col = col_map,
           Rowv = dend_net, Colv = dend_net,
           RowSideColors = cluster_cols  ,
           ColSideColors =   nd_cols )
```


### Dynamic tree cutting to extract clusters 
```
filt_min = 6 
medK = 0.5
flag_plot = FALSE; flag_med = TRUE; flag_dist = TRUE;
frac = 0.995; deep_split = 2; min_cs = 2;

## Get distance matrix from co-expression
temp <- network

## Make as a binary network based on median
if( flag_med == TRUE) {
  temp[temp > medK] <- 1
  temp[temp <= medK] <- 0
  diag(temp) <- 0
}

gene_names <- rownames(temp)
```

```
# Re-calculate distances between genes for distance matrix
if( flag_dist == TRUE) {
  dist_temp <- dist(temp)
} else {
  dist_temp <- as.dist(temp)
}

# Cluster genes using distance matrix
clust_tree <- hclust(dist_temp ) # use this instead of default
clust_dend <- as.dendrogram(clust_tree)
```
 

```
# Extract clusters/modules
if( flag_dist == TRUE ){
  max_h = max(clust_tree$height)
} else {
  max_h = 1
}

unmerged_modules <- dynamicTreeCut::cutreeDynamic(dendro = clust_tree,
                                                  distM = as.matrix(dist_temp),
                                                  deepSplit = deep_split,
                                                  cutHeight = frac * max_h,
                                                  minClusterSize = min_cs,
                                                  pamRespectsDendro = FALSE)
### Play around with the parameters here 
unmerged_modules <- dynamicTreeCut::cutreeDynamic(dendro = clust_tree,
                                                  distM = as.matrix(dist_temp),
                                                  deepSplit = 4, 
                                                  cutHeight = frac * max_h,
                                                  minClusterSize = 50,
                                                  pamRespectsDendro = FALSE)

```

### 
```
n_max <- max(unmerged_modules)
n_l <- sum(unmerged_modules==0)
merged_modules <- unmerged_modules
if( n_l > 0) { merged_modules[unmerged_modules==0] <- (1:n_l ) + n_max }
merged_modules <- merged_modules[clust_tree$order]

# Re-label clusters ids
i.prev <- ""
ki <- 1
ji <- 0
remerged_modules <- as.numeric(merged_modules) * 0
for (ii in as.numeric(merged_modules) ) {
  if (ii == i.prev) {
    remerged_modules[ki] <- ji
  } else {
    i.prev <- ii
    ji <- ji + 1
    remerged_modules[ki] <- ji
    
  }
  ki <- ki + 1
}
```

```
# Generate colors for modules
# Total number of clusters
nsclust <- as.numeric(remerged_modules) + 1
merged_colors <- viridis::magma(max(nsclust))[nsclust]
m <- match( gene_names, gene_names[clust_tree$order] )

## Plot with new clusters 
heatmap.2( network, density.info = "none", trace = "none",
             col = col_map,
             Rowv = clust_dend, Colv = clust_dend,
             RowSideColors = merged_colors[m],
             ColSideColors = merged_colors[m] )
 

plot(as.phylo(clust_tree), type="unrooted", tip.col =  merged_colors[m])
```

```
# Compare to previous clusters
plot(as.phylo(clust_tree), type="c", tip.col =  cluster_cols )
plot(as.phylo(tree_net), type="c", tip.col =  merged_colors[m] )

#plot(as.phylo(clust_tree), type="c", tip.col =  nd_cols )
#plot(as.phylo(tree_net), type="c", tip.col =  nd_cols )
```




