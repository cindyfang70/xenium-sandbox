suppressPackageStartupMessages({
    library(Voyager)
    library(patchwork)
    library(SpatialFeatureExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(bluster)
    library(stringr)
    library(BiocSingular)
    library(scater)
    library(rjson)
    library(Matrix)
    library(vroom)
    library(sf)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(ggforce)
    library(interp)
    library(igraph)
    library(rlist)
    library(gridExtra)
})
###########################################################
# Perform Delaunay triangulation on each of the clusters  #
# in the tissue, prune the global long edges, and plot.   #
###########################################################
source(here("code", "cindy", "04_delaunay", "delaunay.R"))
args <- commandArgs(trailingOnly = TRUE)
print(args[[1]])
sfe <- readRDS(args[[1]])
clusterName <- args[[2]] # get the name of the clustering output to do the triangulation on


# sfe <- readRDS("Br2743_Mid_SFE_filt.RDS")
# tris <- readRDS("Br2743_Mid_5548-delaunay-lou25.RDS")
# source(here("code","04_delaunay","delaunay.R"))
# clusterName <- "lou25"

tris <- list()
gs <- list()
for(i in 1:length(unique(colData(sfe)[[clusterName]]))){
    clust <- unique(colData(sfe)[[clusterName]])[[i]]
    print(clust)
    colData(sfe)[clust] <- colData(sfe)[[clusterName]]==clust
    
    # do the triangulation 
    tri <- delaunay(sfe, clust, seed=i)
    summary(tri)
    tris <- rlist::list.append(tris, tri)
    g <- graph_from_edgelist(arcs(tri), directed=FALSE)
    gs <- list.append(gs, g)
}

fname <- paste(sfe$region_id[[1]], "delaunay", clusterName, sep="-")
pdfname <- paste0(fname, ".pdf")

pdf(here("plots", "cindy", "04_delaunay", pdfname))

plotSpatialFeature(sfe, clusterName, colGeometryName="cellSeg",
                   linewidth=0)
#print(length(tris))
tris.withLengths <- list()
globalPrunedTris <- list()
localPrunedTris <- list()
for (j in 1:length(tris)){
    tri <- tris[[j]]
    g <- gs[[j]]
    print(summary(tri))
    # compute the length of the edges
    tri <- getEdgesSpatialDistance(tri)
    tris.withLengths <- list.append(tris.withLengths, tri)
    
    longIndsGlobal <- c()
    longIndsLocal <- c()
    for (i in 1:tri$n){
        # compute the global constraint and find the long edges
        gc <- computeGlobalEdgeConstraint(tri, i)
        globalLong <- identifyLongEdges(tri, i, gc)
        longIndsGlobal <- c(longIndsGlobal, globalLong)
        
        lc <- computeLocalEdgeConstraint(tri, g, i, beta=2)
        print(lc)
        localLong <- identifyLongEdges(tri, i, lc)
        print(localLong)
        longIndsLocal <- c(longIndsLocal, localLong)
        
    }
    # plot the unpruned triangulation
    p <- plotDelaunay(tri, sfe)
    print(head(tri$arcs))
    hist <- plotEdgeLengthHistogram(tri)
    
    # plot the pruned tri
    globalPrunedtri <- tri
    globalPrunedtri$arcs <- globalPrunedtri$arcs[-longIndsGlobal,]
    print(head(globalPrunedtri$arcs))
    pruned.p <- plotDelaunay(globalPrunedtri, sfe)
    pruned.hist <- plotEdgeLengthHistogram(globalPrunedtri, title = "Edge lengths (global pruned)")
 
    globalPrunedTris <- list.append(globalPrunedTris, globalPrunedtri)
    
    # prune the local edges and plot
    localPrunedTri <- globalPrunedtri
    localPrunedTri$arcs <- localPrunedTri$arcs[-longIndsLocal,]
    #print(head(localPrunedTri$arcs))
    localPruned.p <- plotDelaunay(localPrunedTri, sfe)
    localPruned.hist <- plotEdgeLengthHistogram(localPrunedTri, title="Edge lengths (local and global pruned)")
    localPrunedTris <- list.append(localPrunedTris, localPrunedTri)
    
    # plot
    do.call(gridExtra::grid.arrange, c(list(p, pruned.p, localPruned.p, 
                                                        hist, pruned.hist, localPruned.hist), 
                                                   ncol=3))
}
plotEdgeLengthViolin(tris.withLengths)    
dev.off()
saveRDS(globalPrunedTris, here("processed-data", "cindy", "04_delaunay", paste0(fname, ".RDS")))
saveRDS(localPrunedTris, here("processed-data", "cindy", "04_delaunay", paste0(fname, "localPruned", ".RDS")))
