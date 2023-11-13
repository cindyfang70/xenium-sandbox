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

#------------------------------------------------------------------------------
# pseudobulking the domains in xenium 
# and seeing how they compare to the domains that you can get from visium. 
# 
# One would hope that similar domains would cluster together in PC space 
# We want to try and quantify any difference we 
# see between the two platforms (assuming an intersecting set of genes)
# Do pseudobulked domains cluster together? do they not?
#------------------------------------------------------------------------------

source(here("code", "01_createSCE", "xenium_helpers.R"))
# read in data
sfe <- readRDS(here("processed-data", "cindy", "slide-5434", 
                    "Br8667_Mid_SFE_filt.RDS"))

spe <- readRDS(here("Br8667_mid-visium-SPE.RDS"))

# pseudobulk across the domain labels
sfe.summed <- aggregateAcrossCells(sfe, id=colData(sfe)[,c("clust_M1_lam0.9_k50_res1.2")])
spe.summed <- aggregateAcrossCells(spe, id=colData(spe)[,c("clust_M1_lam0.9_k50_res1.2")])

# convert the SFE to SPE so that the two datasets can be combined
sfe.summed <- SFEtoSPE(sfe.summed)

# subset the visium data to the same genes that are in xenium
spe.summed <- spe.summed[which(rowData(spe.summed)$gene_name %in% rownames(sfe.summed)),]
# remove the imgData
imgData(spe.summed) <- NULL
spe.all <- cbind(sfe.summed, spe.summed)
