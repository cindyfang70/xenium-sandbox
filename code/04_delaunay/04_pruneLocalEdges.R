suppressPackageStartupMessages({
    library(Voyager)
    library(SpatialFeatureExperiment)
    library(SingleCellExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(stringr)
    library(BiocSingular)
    library(Matrix)
    library(sf)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(ggforce)
    library(interp)
    library(igraph)
    library(rlist)
    library(alphahull)
    library(sp)
})

sfe <- readRDS("Br2743_Mid_SFE_filt.RDS")
tris <- readRDS("Br2743_Mid_5548-delaunay-lou25.RDS")
source(here("code","04_delaunay","delaunay.R"))
tri <- tris[[1]]
g <- igraph::graph.edgelist(arcs(tri), directed=FALSE)


v <- 1



#  2-Order_Mean(Pi) is the mean length of all edges that belong to a path of 
# two or fewer edges starting from Pi in Gi.

# compute the path lengths:
computePathLength <- function(tri, path){
    # compute the length of the path from the first vertex in the given path
    # to the last
    edges <- tri$arcs
    total.length <- 0
    for (i in 2:length(path)){
        from <- path[i-1]
        to <- path[i]
        
        # get the edge names, could be either from-to or to-from
        edge.name <- paste(from, to, sep="-") 
        edge.name <- c(edge.name, paste(to, from, sep="-"))
        
        # get the length of the edge
        el <- edges[rownames(edges) %in% edge.name,]
        length <- el$edge.lengths
        
        # add the edge length to the total
        total.length <- total.length + length
    }

    return(total.length)
}

# Mean_Variation(Pi) is the mean value of all Local_Variation(Qi) for 
# Qi belonging to a path of two or fewer edges starting at Pi.


# Local_Variation(Qi) be the standard deviation of the length of the edges 
# directly incident to Qi.
computeLocalVariation <- function(tri, vertex){
    edges <- tri$arcs
    # find the neighbours of the given vertex
    nn <- rbind(edges[which(edges$from ==vertex),],
                edges[which(edges$to ==vertex),])
    
    # compute the mean of neighbouring edge lengths
    local_variation <- sd(nn$edge.lengths)
    return(local_variation)
}



computeLocalEdgeConstraint <- function(tri, g, vertex, beta=2){
    
    # find all simple paths of length 2 originating from vertex
    paths <- all_simple_paths(g, v, V(g), cutoff=2)
    
    # find the mean of the edges involved in the paths (two order mean)
    path.edge.lengths <- sum(unlist(lapply(paths, computePathLength, tri=tri)))
    paths.num.edges <- sum(unlist(lapply(paths, length)))
    
    two_order_mean <- path.edge.lengths/paths.num.edges
    
    # find the nodes belonging to a path of two or fewer edges starting at v
    two.order.neighbours <- unlist(paths)
    two.order.neighbours <- unique(two.order.neighbours[two.order.neighbours!=v])
    
    # compute the mean sd of the edges of the two order neighbours
    two.order.vars <- lapply(two.order.neighbours, computeLocalVariation, tri=tri)
    mean_variation <- mean(unlist(two.order.vars))
    
    local_distance_constraint <- two_order_mean + beta*mean_variation
    return(local_distance_constraint)
}

computeLocalEdgeConstraint(tri, g, v)
