suppressPackageStartupMessages({
    library(Voyager)
    #library(SFEData)
    #library(patchwork)
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
    library(sf)
    library(BiocParallel)
    library(dplyr)
    library(here)
    #library(cmdArgs)
})
###########################################################
# Cluster cells in the tissue using GLM-PCA normalization #
# and then SNN graph-based clustering                     #
###########################################################
args <- commandArgs(trailingOnly=TRUE)

# Read in the data
sfe <- readRDS(args[[1]])
numNeighbours <- args[[2]] # k parameter for clustering

print(sfe)

# compute null residuals and then perform PCA on them (GLM-PCA)
sfe <- scry::nullResiduals(sfe, assay="counts", fam="poisson", type="pearson")
print('null residuals done')
sfe <- scater::runPCA(sfe, ncomponents=50, 
                           ntop = 1000,
                           exprs_values = "poisson_pearson_residuals",
                           scale = TRUE, name = "GLM-PCA",
                           BSPARAM = BiocSingular::RandomParam())
print('pca done')

# cluster using louvain
set.seed(10115023)
g <- scran::buildSNNGraph(sfe, k=numNeighbours, use.dimred = "GLM-PCA")
lou <- igraph::cluster_louvain(g)

# add the clusters back into the SFE object and save it
clusterName <- paste0("lou",numNeighbours)
colData(sfe)[[clusterName]] <- paste0("Lou", lou$membership)
saveRDS(sfe, args[[1]])

pdfname <- paste(sfe$region_id[[1]], "louvain", "k", numNeighbours, sep="-")
pdfname <- paste0(pdfname, ".pdf")
print(pdfname)

# plot the PCA plots and the clusters on the tissue
pdf(here("plots", "cindy", "03_clustering", pdfname))
plotReducedDim(br6471_p, ncomponents=4, colour_by="total_counts", dimred="GLM-PCA")
ElbowPlot(br6471_p, reduction="GLM-PCA", ndims=50)
plotSpatialFeature(br6471_p, "isLou15", colGeometryName = "cellSeg") 
dev.off()

