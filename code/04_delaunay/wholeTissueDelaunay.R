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
#
# Build the delaunay triangulation for the entire tissue and prune the edges

# read in the data

saveCellIds <- function(slide_number){
    
    slide_directory <- here("data", paste0("slide-", slide_number, "/"))
    sfe.paths <- system(paste0("ls ", slide_directory), intern=TRUE)
    sfe.paths <- sfe.paths[grepl("*filt.RDS", sfe.paths)]
    sfe.paths <- paste0(slide_directory, sfe.paths)
    sfes <- lapply(sfe.paths, readRDS)
    
    # save the cell ids in the filtered datasets, so that preprocessing doesn't 
    # need to be done in Python
    for (sfe in sfes){
        cell_ids <- sfe$cell_id
        write.csv(cell_ids, paste0(slide_directory, unique(sfe$region_id),"_cell_ids", ".csv"))
    
    }
    return(sfes) # return the sfes for downstream stuff
}
sfes <- saveCellIds("5434")

delaunayAllRegions <- function(sfes){
    tris <- list()
    for (i in 1:length(sfes)){
        sfe <- sfes[[i]]
        coords <- spatialCoords(sfe)
        tri <- interp::tri.mesh(x=coords[,1],
                                y=coords[,2])
        
        # save the arcs as a data frame for easier usage later 
        tri$arcs <- as.data.frame(tri$arcs)
        rownames(tri$arcs) <- paste(tri$arcs$from, tri$arcs$to, sep="-")
        
        tri <- getEdgesSpatialDistance(tri)
        
        # prune the global long edges
        longIndsGlobal <- c()
        for (i in 1:tri$n){
            gc <- computeGlobalEdgeConstraint(tri, i)
            globalLong <- identifyLongEdges(tri, i, gc)
            longIndsGlobal <- c(longIndsGlobal, globalLong)
        }
        globalPrunedtri <- tri
        globalPrunedtri$arcs <- globalPrunedtri$arcs[-longIndsGlobal,]
        tris <- rlist::list.append(tris, globalPrunedtri)
    }
    return(tris)
}

tris <- delaunayAllRegions(sfes)

# after pruning the global long edges, get the weighted adjacency matrix

getAdj <- function(tri){
    g <- graph_from_edgelist(arcs(tri))
    E(g)$weight <- tri$arcs$edge.lengths
    
    adj <- as_adjacency_matrix(g, attr="weight")
    
    return(adj)
}

adjs <- lapply(tris, getAdj)


slide_number <- "5434"
for (i in 1:length(adjs)){
    slide_directory <- here("data", paste0("slide-", slide_number, "/"))
    writeMM(adjs[[i]], file=paste0(slide_directory, unique(sfes[[i]]$region_id), "globalPrunedTriAdj.mtx"))
}
