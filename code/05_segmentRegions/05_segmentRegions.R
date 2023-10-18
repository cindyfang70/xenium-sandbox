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
})
###########################################################
# Find the connected components within the pruned         #
# triangulation graph, and compute the alpha shape        #
###########################################################

source(here("code","cindy","04_delaunay","delaunay.R"))
sfe <- readRDS(args[[1]])
tri <- readRDS(args[[2]])

sfe <- readRDS("code/05_segmentRegions/Br2743_Mid_SFE_filt.RDS")
tris <- readRDS("code/05_segmentRegions/Br2743_Mid_5548-delaunay-lou25.RDS")
source(here("code","04_delaunay","delaunay.R"))

plist <- list()
for (i in 1:length(tri)){
    tri <- tris[[i]]
    g <- igraph::graph.edgelist(arcs(tri)) # build igraph
    
    # find biggest component
    comp <- components(g) # components() takes a graph as input, not adjacency matrix
    max.inds <- which.max(comp$csize)
    max.comp.vertices <- igraph::groups(comp)[[max.inds]]
    
    # compute the alpha shape using the maximal component
    comp.x <- tri$x[max.comp.vertices]
    comp.y <- tri$y[max.comp.vertices]
    
    alpha.shape <- ashape(comp.x, comp.y, alpha=100)
    alpha.shape$delvor.obj$tri.obj$arcs <- as.data.frame(
        arcs(alpha.shape$delvor.obj$tri.obj))
    
    # plot the alpha shape on the tissue
    alpha.edges <- as.data.frame(alpha.shape$edges)
    from.x <- alpha.edges$x1
    from.y <- alpha.edges$y1
    
    to.x <- alpha.edges$x2
    to.y <- alpha.edges$y2
    
    seg.df <- as.data.frame(cbind(from.x, from.y, to.x,to.y))
    


    p <- plotGeometry(sfe, type="cellSeg")+
        geom_segment(data=seg.df, aes(x=from.x,xend = to.x, y=from.y,yend = to.y))+
        geom_segment(data=scale.df, aes(x=x, xend=xend,y=y,yend=yend))+
        geom_text(data=scale.df, aes(label=text, x=xend, y=yend-300))
    
    plist <- list.append(plist, p)

}

fname <- paste(sfe$region_id[[1]], clusterName, "alphashape", sep="-")
pdfname <- paste0(fname, ".pdf")

pdf(here("plots", "cindy", "05_segmentRegions", pdfname))
do.call(gridExtra::grid.arrange, plist)
dev.off()
