suppressPackageStartupMessages({
    library(Voyager)
    library(SpatialFeatureExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(stringr)
    library(BiocSingular)
    #library(scater)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(gridExtra)
    library(Banksy)
    library(scuttle)
    library(escheR)
    library(RColorBrewer)
    library(data.table)
    library(pals)
    library(alphahull)
})
#------------------------------------------------------------#
# Run banksy on all 6 tissue regions at the ame time
# tutorial: https://github.com/prabhakarlab/Banksy/tree/bioc
#------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
source(here("code", "01_createSCE", "xenium_helpers.R"))
source(here("code", "04_delaunay", "delaunay.R"))

# read in the SPE
sfename <- here("processed-data", "cindy", "all-tissues-spe-with-banksy-cell-types.RDS")
spe <- readRDS(sfename)

tab <- prop.table(table(spe$clust_M1_lam0_k50_res1, spe$preds_from_bayesspace_NMF_k200),2)
tab

sub.spe <- spe[,which(spe$region_id == "Br6471_Post_5434")]
plist <- c()
for (i in 1:length(unique(sub.spe$preds_from_bayesspace_NMF_k200))){
    layer <- unique(sub.spe$preds_from_bayesspace_NMF_k200)[[i]]
    sub.spe.layer <- sub.spe[,which(sub.spe$preds_from_bayesspace_NMF_k200==layer)]
    
    tri <- interp::tri.mesh(x=spatialCoords(sub.spe.layer)[,1],
                            y=spatialCoords(sub.spe.layer)[,2])
    g <- graph_from_edgelist(arcs(tri), directed=FALSE)
    
    tri$arcs <- as.data.frame(tri$arcs)
    rownames(tri$arcs) <- paste(tri$arcs$from, tri$arcs$to, sep="-")
    print(summary(tri))
    
    tri <- getEdgesSpatialDistance(tri)
    
    longIndsGlobal <- c()
    longIndsLocal <- c()
    for (j in 1:tri$n){
        # compute the global constraint and find the long edges
        gc <- computeGlobalEdgeConstraint(tri, j)
        globalLong <- identifyLongEdges(tri, j, gc)
        longIndsGlobal <- c(longIndsGlobal, globalLong)
        #print(gc)
        
        lc <- computeLocalEdgeConstraint(tri, g, j, beta=2)
        #print(lc)
        localLong <- identifyLongEdges(tri, j, lc)
        #print(localLong)
        longIndsLocal <- c(longIndsLocal, localLong)
        
    }
    print(longIndsGlobal)
    print(longIndsLocal)
    
    globalPrunedtri <- tri
    globalPrunedtri$arcs <- globalPrunedtri$arcs[-longIndsGlobal,]
    
    print(summary(globalPrunedtri))
    
    pruned.g <- igraph::graph_from_edgelist(arcs(globalPrunedtri))
    comps <- components(pruned.g) 
    
    comp.vertices <- igraph::groups(comps)[[which.max(comps$csize)]] # only use the first component for now
    comp.x <- globalPrunedtri$x[comp.vertices]
    comp.y <- globalPrunedtri$y[comp.vertices]
    
    alpha.shape <- ashape(x=comp.x, y=comp.y, alpha=100)
    print(alpha.shape)

    seg.df <- as.data.frame(alpha.shape$edges)[,3:6]
    colnames(seg.df) <- c("from.x", "from.y", "to.x", "to.y")

    p <- make_escheR(sub.spe, y_reverse=FALSE) %>%
        add_fill("clust_M1_lam0_k50_res1")+
        scale_fill_discrete()+
        geom_segment(data=seg.df, aes(x=from.x,xend = to.x, y=from.y,yend = to.y))+
        ggtitle(layer)
    plist[[i]] <- p


}

pdf(here("plots", "cindy", "06_cellTyping", "alpha-shapes", "Br6471_post_5434-layer-alpha-shapes.pdf"),
    height=20, width=20)
do.call(gridExtra::grid.arrange, c(plist, ncol=3))
dev.off()
