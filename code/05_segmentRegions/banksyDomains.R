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
#------------------------------------------------------------#
# Trying out Banksy for spatial domain detection:
# tutorial: https://github.com/prabhakarlab/Banksy/tree/bioc
#------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
print(args[[1]])
sfe <- readRDS(args[[1]])

#sfe <- readRDS("data/slide-5434/Br6522_Post_SFE_filt.RDS")
sfe <- computeLibraryFactors(sfe)
aname <- "normcounts"
assay(sfe, aname) <- normalizeCounts(sfe, log = FALSE)


# Compute neighbourhood matrices
lambda <- c(0, 0.9)
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

p1 <- make_escheR(sfe) %>%
    add_fill("clust_M1_lam0.9_k50_res1.2") +
    scale_fill_brewer("Set3")
    

#p1 <- plotSpatialFeature(sfe, "clust_M1_lam0.9_k50_res1.2", colGeometryName = "cellSeg")

fname <- paste(sfe$region_id[[1]], "Banksy", "lambda", 
               lambda[[2]], "res", 1.2, sep="-")
pdfname <- paste0(fname, ".pdf")
pdf(here("plots", "cindy", "05_segmentRegions", pdfname))
p1
dev.off()

