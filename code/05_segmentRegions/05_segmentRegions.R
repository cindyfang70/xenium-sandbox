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
tri <- tris[[1]]
el <- igraph::graph.edgelist(arcs(tri))
# adj <- get.adjacency(el)
#dimnames(adj) <- list(V(el), V(el))

comp <- components(el) # components() takes a graph as input, not adjacency matrix
max.inds <- which.max(comp$csize)
max.comp.vertices <- groups(comp)[[max.inds]]

tri.comp <- tri

comp.x <- tri$x[max.comp.vertices]
comp.y <- tri$y[max.comp.vertices]

alpha.shape <- ashape(comp.x, comp.y, alpha=100)
alpha.shape$delvor.obj$tri.obj$arcs <- as.data.frame(arcs(alpha.shape$delvor.obj$tri.obj))
#library(Voyager)

alpha.edges <- as.data.frame(alpha.shape$edges)
from.x <- alpha.edges$x1
from.y <- alpha.edges$y1

to.x <- alpha.edges$x2
to.y <- alpha.edges$y2

seg.df <- as.data.frame(cbind(from.x, from.y, to.x,to.y))

# x_min <- min(spatialCoords(sfe)[,1])
# y_min <- min(spatialCoords(sfe)[,2])
# scale.df <- data.frame(x=x_min + 100, xend=x_min + 600,
#                        y=y_min+100, yend=y_min+100, text="500um")

p <- plotGeometry(sfe, type="cellSeg")+
    geom_segment(data=seg.df, aes(x=from.x,xend = to.x, y=from.y,yend = to.y))
    #geom_segment(data=scale.df, aes(x=x, xend=xend,y=y,yend=yend))+
    #geom_text(data=scale.df, aes(label=text, x=xend, y=yend-300))

comps <- list()
for (i in 1:length(tri)){
    
    # build an igraph object
    g <- graph.edgelist(arcs(tri))
    
    # find the components
    comp <- components(g)
    
    # which vertices are in each component
    groups(comp)
    
    # we only care about components that are sufficiently large. so we want
    # to identify the "large enough" components and then subset the triangulation
    # object to contain only those vertices in the component. Then, we can plot them.
    
    
}

