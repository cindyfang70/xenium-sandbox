suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SpatialExperiment)
    library(stringr)
    library(scuttle)
    library(scater)
    library(dplyr)
    library(here)
})

#-------------------------------------------------------------------------------
# Read in the Visum data from /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
# and subset the SPE object to the samples that we have Xenium replicates for
#-------------------------------------------------------------------------------

# Define the samples for which we want to get the Visium data
sample_info <- c("Br6471_post", "Br6522_post", "Br8667_mid", "Br2743_mid")

load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_raw_final.Rdata")


for (i in 1:length(sample_info)){
    sample_id <- sample_info[[i]]
    spe <- spe_raw[,spe_raw$sample_id == sample_id]
    spe
    saveRDS(spe, here("processed-data", "cindy", "visium", paste0(sample_id, "-visium-SPE.RDS")))
}