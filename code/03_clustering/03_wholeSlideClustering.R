suppressPackageStartupMessages({
    library(Voyager)
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
})

#########################################################
# Whole slide clustering in the tissue using GLM-PCA    #
# normalization and then SNN graph-based clustering     #
#########################################################
args <- commandArgs(trailingOnly=TRUE)

k <- args[[length(args)]] 
args <- args[-length(args)]
print(args)

sfes <- lapply(args, readRDS)
sfe <- do.call(cbind, sfes)


print(sfe)

sfe <- scry::nullResiduals(sfe, assay="counts", fam="poisson", type="pearson")
print('null residuals done')
sfe <- scater::runPCA(sfe, ncomponents=50, 
                      ntop = 1000,
                      exprs_values = "poisson_pearson_residuals",
                      scale = TRUE, name = "WholeSlide-GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

set.seed(10222023)
g <- scran::buildSNNGraph(sfe, k=as.numeric(k), use.dimred = "WholeSlide-GLM-PCA")
lou <- igraph::cluster_louvain(g)

# add the clusters back into the SFE object and save it
clusterName <- paste0("louWholeSlide",numNeighbours)
colData(sfe)[[clusterName]] <- paste0("LouWholeSlide", lou$membership)

region_id <- unique(sfe$region_id)
slide <- unlist(strsplit(region_id, split="_"))[[3]]

sfeName <- paste0("slide", region_id, "-clustSFE.RDS")

saveRDS(sfe, here("processed-data", "cindy", paste0("slide-", slide), sfeName))


