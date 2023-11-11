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
})

#-------------------------------------------------------------------------------
# Run Banksy on the Visium samples
# TODO: try using scry::nullResiduals instead of normalizeCounts, then follow
# the canonical Banksy workflow
#-------------------------------------------------------------------------------

spe <- readRDS(here("processed-data", "cindy", "visium", "Br8667_mid-visium-SPE.RDS"))

# library size normalization
spe <- computeLibraryFactors(spe)
aname <- "normcounts"
assay(spe, aname) <- normalizeCounts(spe, log = FALSE)


# Compute neighbourhood matrices
lambda <- c(0, 0.9)
k_geom <- c(15, 30)


# run Banksy
spe <- Banksy::computeBanksy(spe, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)

set.seed(1000)
spe <- Banksy::runBanksyPCA(spe, use_agf = TRUE, lambda = lambda)
spe <- Banksy::runBanksyUMAP(spe, use_agf = TRUE, lambda = lambda)
spe <- Banksy::clusterBanksy(spe, use_agf = TRUE, lambda = lambda, resolution = 1.2)

# connect the clusters as suggested in the tutorial
spe <- Banksy::connectClusters(spe)

# Plot using escheR
print(head(colData(sfe)["clust_M1_lam0.9_k50_res1.2"]))
colourCount = nlevels(colData(sfe)[["clust_M1_lam0.9_k50_res1.2"]])
print(colourCount)
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

p1 <- make_escheR(sfe) %>%
    add_fill(var="clust_M1_lam0.9_k50_res1.2")+
    scale_fill_manual(values=getPalette(colourCount))