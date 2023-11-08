suppressPackageStartupMessages({
    library(Voyager)
    library(patchwork)
    library(SpatialFeatureExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(bluster)
    library(stringr)
    library(BiocSingular)
    #library(scater)
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
    library(Banksy)
})
#------------------------------------------------------------#
# Trying out Banksy for spatial domain detection:
# tutorial: https://github.com/prabhakarlab/Banksy/tree/bioc
#------------------------------------------------------------#
sfe <- readRDS("data/slide-5434/Br8667_Mid_SFE_filt.RDS")
sfe <- computeLibraryFactors(sfe)
aname <- "normcounts"
assay(sfe, aname) <- normalizeCounts(sfe, log = FALSE)


# Compute neighbourhood matrices
lambda <- c(0, 0.2)
k_geom <- c(15, 30)

sfe <- Banksy::computeBanksy(sfe, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)

set.seed(1000)
sfe <- Banksy::runBanksyPCA(sfe, use_agf = TRUE, lambda = lambda)
sfe <- Banksy::runBanksyUMAP(sfe, use_agf = TRUE, lambda = lambda)
sfe <- Banksy::clusterBanksy(sfe, use_agf = TRUE, lambda = lambda, resolution = 1.2)

# connect the clusters as suggested in the tutorial
sfe <- Banksy::connectClusters(sfe)


cnames <- colnames(colData(sfe))
cnames <- cnames[grep("^clust", cnames)]
colData(sfe) <- cbind(colData(sfe), spatialCoords(sfe))

plotSpatialFeature(sfe, "clust_M1_lam0.2_k50_res1.2", colGeometryName = "cellSeg")


