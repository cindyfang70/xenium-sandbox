suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SpatialExperiment)
    library(stringr)
    library(scuttle)
    library(scater)
    library(dplyr)
    library(here)
    library(spatialLIBD)
})

#-------------------------------------------------------------------------------
# Read in the Visum data from /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
# and subset the SPE object to the samples that we have Xenium replicates for
#-------------------------------------------------------------------------------
ehub <- ExperimentHub::ExperimentHub()
raw_sce <- fetch_data(type = "spatialDLPFC_Visium", eh = ehub)
sce <- raw_sce

subject_positions <- c("Br6471_Post", "Br6522_Post", "Br8667_Mid", "Br2743_Mid")
sce$position <- substr(sce$position, 1,3)
sce$position[grepl("pos", sce$position)] <- "Post"
sce$position <- str_to_title(sce$position)

sce$subject_position <- paste(sce$subject, sce$position, sep="_")


# Use BayesSpace_harmony_09 as the layer annotations.
# This is what's shown in the Visium paper and plotted here:
# https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/03_BayesSpace/03_BayesSpace_big_plots.R

# subset the big SPE to the ones we have Xenium replicates for
spes <- list()
for (id in subject_positions){
    spe <- sce[,sce$subject_position==id]
    spes <- rlist::list.append(spes, spe)
    fname <- paste0(id, "-Visium-SPE.RDS")
    saveRDS(spe, here("processed-data", "cindy", "visium", fname))
}