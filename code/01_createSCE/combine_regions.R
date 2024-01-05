suppressPackageStartupMessages({
    #library(Voyager)
    library(SpatialFeatureExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(stringr)
    library(dplyr)
    library(here)
    library(escheR)
    library(RColorBrewer)
})


#-------------------------------------------------------------------------------
# Combine the region-specific SFE objects into one large SFE.
#-------------------------------------------------------------------------------
source(here("code", "cindy", "01_createSCE", "xenium_helpers.R"))

args <- commandArgs(trailingOnly=TRUE)
print(args)

sfes <- lapply(args, readRDS)

# convert to SPE because cbind is super annoying for SFEs
spes <- lapply(sfes, SFEtoSPE)

sfe.all <- do.call(cbind, spes)
print(sfe.all)
sfe.all$donor_id <- unlist(lapply(strsplit(sfe.all$region_id, 
                                           split="_"), "[", 1))

print(rowData(sfes[[1]])$ID)
rowData(sfe.all)$Symbol <- rownames(sfe.all)
rowData(sfe.all)$ID <- rowData(sfes[[1]])$ID

print(rowData(sfe.all))
saveRDS(sfe.all, here("processed-data", "cindy", "sfe-all.RDS"))