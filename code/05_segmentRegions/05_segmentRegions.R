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
###########################################################
# Find the connected components within the pruned         #
# triangulation graph, and compute the alpha shape        #
###########################################################
args <- commandArgs(trailingOnly=TRUE)
source(here("code","cindy","04_delaunay","delaunay.R"))
sfe <- readRDS(args[[1]])
tris <- readRDS(args[[2]])

sfe <- readRDS("Br2743_Mid_SFE_filt.RDS")
tris <- readRDS("Br2743_Mid_5548-delaunay-lou25.RDS")
source(here("code","04_delaunay","delaunay.R"))

# compute the coordinates for the scale bar
x_min <- min(spatialCoords(sfe)[,1])
y_min <- min(spatialCoords(sfe)[,2])
scale.df <- data.frame(x=x_min + 100, xend=x_min + 600,
                       y=y_min+100, yend=y_min+100, text="500um")
plist <- list()
alpha.shapes <- list()
for (i in 1:length(tris)){
    tri <- tris[[i]]
    #clustName <- unique(tri$arcs$clust) # get the name of the cluster
    clustName <- sprintf("clust%s", i)
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
    
    # area of the alpha shape
    ashape <- alpha.shape
    bds <-ashape$x[alpha.shape$alpha.extremes,]   # matrix of coordinates in ashapebds <- 
    bds <- rbind(bds, bds[1,])          # close the ring
    ashape <- Polygon(bds)   
    print(ashape@area)
    
    # plot the alpha shape on the tissue
    alpha.edges <- as.data.frame(alpha.shape$edges)
    from.x <- alpha.edges$x1
    from.y <- alpha.edges$y1
    
    to.x <- alpha.edges$x2
    to.y <- alpha.edges$y2
    
    seg.df <- as.data.frame(cbind(from.x, from.y, to.x,to.y))
    
    
    sfe <- logNormCounts(sfe)

    p <- plotSpatialFeature(sfe, "MOBP")+
        geom_segment(data=seg.df, aes(x=from.x,xend = to.x, y=from.y,yend = to.y))+
        geom_segment(data=scale.df, aes(x=x, xend=xend,y=y,yend=yend))+
        geom_text(data=scale.df, aes(label=text, x=xend, y=yend-300))+
        ggtitle(clustName)
    
    # save the alpha shapes to plot them all on the same tissue
    seg.df$clust <- clustName # need to change this later to reflect variable cluster names
    alpha.shapes <- list.append(alpha.shapes, seg.df)
    plist <- list.append(plist, p)

}

all.alpha.shapes <- do.call(rbind, alpha.shapes[c(2,4,6,12)])
all.alpha.shapes <- as.data.frame(all.alpha.shapes)

fname <- paste(sfe$region_id[[1]], "alphashape", sep="-")
pdfname <- paste0(fname, ".pdf")

pdf(here("plots", "cindy", "05_segmentRegions", pdfname))
do.call(gridExtra::grid.arrange, plist)
plotGeometry(sfe, type="cellSeg")+
    geom_segment(data=all.alpha.shapes, linewidth=1, aes(x=from.x,xend = to.x, y=from.y,yend = to.y,
                                            colour=as.factor(clust)))+
    geom_segment(data=scale.df, aes(x=x, xend=xend,y=y,yend=yend))+
    geom_text(data=scale.df, aes(label=text, x=xend, y=yend-300))
dev.off()

# see if the coords are in the shape
coords <- colGeometries(sfe)$centroid
bds <- st_as_sf(as.data.frame(bds), coords=c("V1", "V2") )
sf::st_intersects(coords, st_as_sf(Polygon(ashape)))

sp::over(coords, ashape)
sf::st_intersection(coords, st_as_sf(ashape))