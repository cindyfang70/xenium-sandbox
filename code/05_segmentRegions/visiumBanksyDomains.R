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

args <- commandArgs(trailingOnly=TRUE)
spe_path <- args[[1]]
spe <- readRDS(spe_path)
colData(spe) <- colData(spe)[!grepl("^clust", colnames(colData(spe)))]

# spe_path <- here("processed-data", "cindy", "visium", "Br8667_mid-visium-SPE.RDS")
# spe <- readRDS(spe_path)

runBanksy <- function(spe, spe_path){
    
    # library size normalization
    spe <- computeLibraryFactors(spe)
    aname <- "normcounts"
    assay(spe, aname) <- normalizeCounts(spe, log = FALSE)
    
    
    # Compute neighbourhood matrices
    lambda <- c(0, 0.9)
    k_geom <- c(15, 30)
    spe <- Banksy::computeBanksy(spe, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
    
    set.seed(1000)
    spe <- Banksy::runBanksyPCA(spe, use_agf = TRUE, lambda = lambda)
    spe <- Banksy::runBanksyUMAP(spe, use_agf = TRUE, lambda = lambda)
    spe <- Banksy::clusterBanksy(spe, use_agf = TRUE, lambda = lambda, resolution = 1.2)
    
    # connect the clusters as suggested in the tutorial
    spe <- Banksy::connectClusters(spe)
    

    
    saveRDS(spe, spe_path)
    return(spe)
}

spe <- runBanksy(spe, spe_path)

# Plot using escheR

colourCount = nlevels(colData(spe)[["clust_M1_lam0.9_k50_res1.2"]])
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

p1 <- make_escheR(spe) %>%
        add_fill(var="clust_M1_lam0.9_k50_res1.2")+
        scale_fill_manual(values=getPalette(colourCount))+
        ggtitle(unique(spe$sample_id))

pdfname <- paste(spe$sample_id[[1]], "Banksy", "lambda", 
                 0.9, "res", 1.2, sep="-")

pdf(here("plots","cindy", "05_segmentRegions","banksy", "visium", pdfname))
p1
dev.off()