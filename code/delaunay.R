---
title: "Spatial Domain Detection"
author: "Cindy Fang"
date: "2023-10-06"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
suppressPackageStartupMessages({
    library(Voyager)
    library(SFEData)
    library(patchwork)
    library(SpatialFeatureExperiment)
    library(SingleCellExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(bluster)
    library(stringr)
    library(scuttle)
    library(BiocSingular)
    library(scater)
    library(rjson)
    library(Matrix)
    library(DropletUtils)
    library(vroom)
    library(sf)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(schex)
    library(ggforce)
    library(interp)
    library(igraph)
    library(rlist)
})
source(here("code", "xenium_helpers.R"))
```

```{r}
data_dir <- "output-XETG00089__0005434__Region_1__20230831__172339"

# Read in the three regions
br6471_p <- readRDS(here("data", data_dir, "Br6471_Post_SFE.RDS"))

cols_use <- names(colData(br6471_p))[str_detect(
    names(colData(br6471_p)), "_percent$")]
cols_use <- cols_use[!grepl("anti", cols_use)]

for (n in cols_use) {
    br6471_p <- get_neg_ctrl_outliers(n, br6471_p)
}
inds_keep <- br6471_p$nucleus_area <=200 & !br6471_p$is_blank_outlier &
    !br6471_p$is_depr_outlier &!br6471_p$is_negProbe_outlier &
    !br6471_p$is_negCodeword_outlier & br6471_p$nCounts > 0

br6471_p <- br6471_p[,inds_keep]
```


```{r}
br6471_p <- scry::nullResiduals(br6471_p, assay="counts", fam="poisson", type="pearson")
br6471_p <- scater::runPCA(br6471_p, ncomponents=50, 
                      ntop = 1000,
                      exprs_values = "poisson_pearson_residuals",
                      scale = TRUE, name = "GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

plotReducedDim(br6471_p, ncomponents=4, colour_by="total_counts", dimred="GLM-PCA")
Voyager::ElbowPlot(br6471_p, reduction="GLM-PCA", ndims=50)
```

```{r}
set.seed()
g25 <- scran::buildSNNGraph(br6471_p, k=25, use.dimred = "GLM-PCA")
lou25 <- igraph::cluster_louvain(g25)
br6471_p$louvain25 <- paste0("Lou", lou25$membership)


br6471_p$isLou4 <- br6471_p$louvain25=="Lou4"
clusts.p <- plotSpatialFeature(br6471_p, "louvain25", colGeometryName = "cellSeg")
clusts.p

br6471_p$isLou15 <- br6471_p$louvain25=="Lou15"
#plotSpatialFeature(br6471_p, "isLou15", colGeometryName = "cellSeg") 
```

```{r}
delaunay <- function(sfe, cluster, seed){
    set.seed(seed)
    celltype.coords <- spatialCoords(sfe)[colData(sfe)[[cluster]],]
    
    tri <- interp::tri.mesh(x=celltype.coords[,1],
                        y=celltype.coords[,2])
    
    # save the arcs as a data frame for easier usage later 
    tri$arcs <- as.data.frame(tri$arcs)
    rownames(tri$arcs) <- paste(tri$arcs$from, tri$arcs$to, sep="-")
    return(tri)
    
}

plotDelaunay <- function(tri, sfe){
    # make the data frame for the line segments of the triangulation
    from.x <- tri$x[tri$arcs$from]
    from.y <- tri$y[tri$arcs$from]

    to.x <- tri$x[tri$arcs$to]
    to.y <- tri$y[tri$arcs$to]
    

    seg.df <- as.data.frame(cbind(from.x, from.y, to.x,to.y))
    
    # compute coordinates for placing the scale bar
    x_min <- min(spatialCoords(sfe)[,1])
    y_min <- min(spatialCoords(sfe)[,2])
    scale.df <- data.frame(x=x_min + 100, xend=x_min + 600,
                           y=y_min+100, yend=y_min+100, text="500um")
    
    p <- plotGeometry(sfe, type="cellSeg")+
        geom_segment(data=seg.df, aes(x=from.x,xend = to.x, y=from.y,yend = to.y))+
        geom_segment(data=scale.df, aes(x=x, xend=xend,y=y,yend=yend))+
        geom_text(data=scale.df, aes(label=text, x=xend, y=yend-300))
    
    return(p)
    }


```

```{r}
delaunayAdjacencyMatrix <- function(sfe, tri){
    # get the x and y indices of the spatial coordinates that are in the triangulation
    x.inds <- which(spatialCoords(sfe)[,1]%in%tri$x)
    y.inds <- which(spatialCoords(sfe)[,2]%in%tri$y)
    
    # find the intersecting x and y indices, these are the coordinates that are actually 
    # in the triangulation
    celltype.inds <- intersect(x.inds, y.inds)
    
    # get the ids of the cells in the triangulation 
    cell_ids <- sfe$cell_id[celltype.inds]
    
    # get the adjacency matrix for the triangulation
    adj <- get.adjacency(graph.edgelist(as.matrix(interp::arcs(tri)), directed=FALSE))
    
    # add the cell ids to the column and rownames
    dimnames(adj) <- list(cell_ids, cell_ids)
    
    return(adj)
}

# adjs <- c()
# for (i in 1:length(tris)){
#     adj <- delaunayAdjacencyMatrix(sfe=br6471_p, tri=tris[[i]])
#     adjs <- list.append(adjs, adj)
# }
# 
# gs <- sapply(adjs, graph_from_adjacency_matrix)
# joint <- do.call("union", gs)
# 
# joint.adj <- get.adjacency(joint)
```


```{r}
getEdgesSpatialDistance <- function(tri){
    edges <- tri$arcs
    # get the x,y coordinates for each edge
    from.x <- tri$x[edges$from]
    from.y <- tri$y[edges$from]
    
    to.x <- tri$x[edges$to]
    to.y <- tri$x[edges$to]
    
    # compute the edge length
    edge.lengths <- sqrt((from.x-to.x)^2 + (from.y-to.y)^2)
    
    tri$arcs <- cbind(edges, edge.lengths)
    
    return(tri)
}
tri <- getEdgesSpatialDistance(tri)

# prune_global_long_edges <- function(sfe, adj){
#     
# }
```

```{r}
#pdf("delaunayPlots-k25.pdf")
tri.ps <- c()
tris <- list()

for(i in 1:length(unique(br6471_p$louvain25))){
    clust <- unique(br6471_p$louvain25)[[i]]
    print(clust)
    colData(br6471_p)[clust] <- br6471_p$louvain25==clust
    
    # do the triangulation and plot it on the tissue
    tri <- delaunay(br6471_p, clust, seed=i)
    tris <- rlist::list.append(tris, tri)
    
    p <- plotDelaunay(tri, br6471_p)
    # get the adjacency matrix
    #adj <- delaunayAdjacencyMatrix(br6471_p, tri)
    
    # get the edge lengths
    #gdist <- getEdgesSpatialDistance(br6471_p, adj)
    
    # plot the histogram of the edge lengths
    # weight.df <- data.frame(weight=E(gdist)$weight) %>%
    #     filter(weight < 1000)
    # 
    # h <- ggplot(weight.df, aes(x=weight))+
    #     #geom_histogram(bins=50)+
    #     geom_density()+
    #     xlab("Edge lengths")
    # print(p + h)
    
}
#dev.off()
```

```{r}
computeGlobalEdgeConstraint <- function(tri, vertex){
    # compute the global edge constraint for a single vertex detailed in https://www.sciencedirect.com/science/article/pii/S0098300411004419#s0010
    edges <- tri$arcs
    # compute the global edge length mean and sd
    global_mean <- mean(edges$edge.lengths)
    global_variation <- sd(edges$edge.lengths)
    
    # find the neighbours of the given vertex
    nn <- rbind(edges[which(edges$from ==vertex),],
                edges[which(edges$to ==vertex),])
    print(nn)
    
    # compute the mean of neighbouring edge lengths
    local_mean <- mean(nn$edge.lengths)
    alpha <- global_mean / local_mean
    
    global_distance_constraint <- global_mean + alpha * global_variation
    
    return(global_distance_constraint)
}

gc <- computeGlobalEdgeConstraint(tri, vertex=1)
gc
```

```{r}
identifyGlobalLongEdges <- function(tri, vertex, global_constraint){
    
    edges <- tri$arcs
    
    # find the neighbours of the given vertex
    nn <- rbind(edges[which(edges$from == vertex),],
                edges[which(edges$to == vertex),])
    
    nnLong <- nn[which(nn$edge.lengths > global_constraint),]
    
    # indicate which edges are long and store it in the triangulation
    longInds <- which(rownames(edges) %in% rownames(nnLong))
    
    return(longInds)
}

pdf("delaunaypruned-k25.pdf")
# find the global constraint for each node
for (j in 1:length(tris)){
    tri <- tris[[j]]
    summary(tri)
    longInds <- c()
    for (i in 1:tri$n){
        gc <- computeGlobalEdgeConstraint(tri, i)
        print(gc)
        long <- identifyGlobalLongEdges(tri, i, gc)
        longInds <- c(longInds, long)
    # print("global constraint")
    # print(gc)
    }
    p <- plotDelaunay(tri, br6471_p)
    globalPrunedtri <- tri
    globalPrunedtri$arcs <- globalPrunedtri$arcs[-longInds,]

    pruned.p <- plotDelaunay(globalPrunedtri, br6471_p)
    print(p + pruned.p)
}
dev.off()
```


