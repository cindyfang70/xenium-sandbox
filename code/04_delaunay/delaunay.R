suppressPackageStartupMessages({
    library(Voyager)
    library(patchwork)
    library(SpatialFeatureExperiment)
    library(SingleCellExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(stringr)
    library(scuttle)
    library(BiocSingular)
    #library(scater)
    #library(rjson)
    library(Matrix)
    #library(vroom)
    library(sf)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(ggforce)
    library(interp)
    library(igraph)
    library(rlist)
})
#source(here("code", "cindy", "01_createSCE", "xenium_helpers.R"))

###############################################################
# helper functions for making delaunay triangulations.        #
###############################################################
delaunay <- function(sfe, cluster, seed){
    set.seed(seed)
    celltype.coords <- spatialCoords(sfe)[colData(sfe)[[cluster]],]
    
    tri <- interp::tri.mesh(x=celltype.coords[,1],
                            y=celltype.coords[,2])
    
    # save the arcs as a data frame for easier usage later 
    tri$arcs <- as.data.frame(tri$arcs)
    rownames(tri$arcs) <- paste(tri$arcs$from, tri$arcs$to, sep="-")
    
    # save the cluster that the triangulation is for in the edges df
    tri$arcs$clust <- cluster
    return(tri)
    
}

plotDelaunay <- function(tri, sfef){
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
        geom_text(data=scale.df, aes(label=text, x=xend, y=yend-300))+
        ggtitle(tri$arcs$clust[[1]])
    
    return(p)
}

delaunayAdjacencyMatrix <- function(tri){
    # get the x and y indices of the spatial coordinates that are in the triangulation
    #x.inds <- which(spatialCoords(sfe)[,1]%in%tri$x)
    #y.inds <- which(spatialCoords(sfe)[,2]%in%tri$y)
    
    # find the intersecting x and y indices, these are the coordinates that are actually 
    # in the triangulation
    #celltype.inds <- intersect(x.inds, y.inds)
    
    # get the ids of the cells in the triangulation 
    #cell_ids <- sfe$cell_id[celltype.inds]
    
    # get the adjacency matrix for the triangulation
    arcs <- tri$arcs[,-3]
    adj <- get.adjacency(graph.edgelist(as.matrix(arcs), directed=FALSE))
    
    # add the cell ids to the column and rownames
    dimnames(adj) <- list(arcs$to, arcs$to)
    
    return(adj)
}

getEdgesSpatialDistance <- function(tri){
    # compute the length of each edge (euclidean distance between the connected vertices)
    edges <- tri$arcs
    # get the x,y coordinates for each edge
    from.x <- tri$x[edges$from]
    from.y <- tri$y[edges$from]
    
    to.x <- tri$x[edges$to]
    to.y <- tri$y[edges$to]
    
    # compute the edge length
    edge.lengths <- sqrt((from.x-to.x)^2 + (from.y-to.y)^2)
    
    tri$arcs <- cbind(edges, edge.lengths)
    
    return(tri)
}

computeGlobalEdgeConstraint <- function(tri, vertex){
    # compute the global edge constraint for a single vertex detailed in https://www.sciencedirect.com/science/article/pii/S0098300411004419#s0010
    edges <- tri$arcs
    # compute the global edge length mean and sd
    global_mean <- mean(edges$edge.lengths)
    global_variation <- sd(edges$edge.lengths)
    
    # find the neighbours of the given vertex
    nn <- rbind(edges[which(edges$from ==vertex),],
                edges[which(edges$to ==vertex),])
    #print(nn)
    
    # compute the mean of neighbouring edge lengths
    local_mean <- mean(nn$edge.lengths)
    alpha <- global_mean / local_mean
    
    global_distance_constraint <- global_mean + alpha * global_variation
    
    return(global_distance_constraint)
}

identifyGlobalLongEdges <- function(tri, vertex, global_constraint){
    
    edges <- tri$arcs
    
    # find the neighbours of the given vertex
    nn <- rbind(edges[which(edges$from == vertex),],
                edges[which(edges$to == vertex),])
    
    # identify the neighbours with long edges
    nnLong <- nn[which(nn$edge.lengths > global_constraint),]
    
    # find the indices of long edges
    longInds <- which(rownames(edges) %in% rownames(nnLong))
    
    return(longInds)
}

plotEdgeLengthHistogram <- function(tri, title="Edge lengths"){
    edges <- tri$arcs
    edges <- edges %>%
        filter(edge.lengths <= 1000)
    
    p <- ggplot(edges, aes(x=edge.lengths))+
        geom_density()+
        ggtitle(title)
    
    return(p)
}

plotEdgeLengthViolin <- function(tris, removeLong=TRUE){
    # take in a list of triangulation objects
    allEdges <- list()
    for (i in 1:length(tris)){
        tri <- tris[[i]]
        edges <- tri$arcs
        head(edges)
        allEdges <- list.append(allEdges, edges)
    }
    
    allEdges <- do.call(rbind, allEdges)
    print(head(allEdges))
    
    if (removeLong){
        allEdges <- allEdges %>%
            as.data.frame() %>%
            filter(edge.lengths<=1000)
    }
    
    # make the cluster label a factor so they don't get rearranged
    allEdges$clust <- factor(allEdges$clust, 
                             sfelevels=unique(allEdges$clust))
    
    p <- ggplot(allEdges, aes(x=clust, y=edge.lengths))+
        geom_violin()
    
    return(p)
}
