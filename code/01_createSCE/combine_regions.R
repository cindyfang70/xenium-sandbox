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
# remove the rowData because they don't match
# sfes <- lapply(sfes, function(x){rowData(x)$means <- NULL
# rowData(x)$vars <- NULL
# rowData(x)$cv2 <- NULL
# return(x)})

# convert to SPE because cbind is super annoying for SFEs
sfes <- lapply(sfes, SFEtoSPE)

sfe.all <- do.call(cbind, sfes)
sfe.all$Sample <- "sample01"
print(sfe.all)
sfe.all$donor_id <- lapply(strsplit(sfe.all$region_id, split="_"), "[", 1)[[1]]

rowData(sfe.all)$Symbol <- rownames(sfe.all)
rowData(sfe.all)$ID <- rowData(sfes[[1]])$ID

saveRDS(sfe.all, here("processed-data", "cindy", "sfe-all.RDS"))