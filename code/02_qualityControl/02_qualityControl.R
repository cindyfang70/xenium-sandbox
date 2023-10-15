suppressPackageStartupMessages({
    library(Voyager)
    library(SFEData)
    library(patchwork)
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
    library(DropletUtils)
    library(vroom)
    library(sf)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(schex)
})


################################################
# Filter out low quality cells for each region #
################################################


##### Slide 5434 #####
data_dir <- "output-XETG00089__0005434__Region_1__20230831__172339"

# Read in the three regions
br6471_p <- readRDS(here("data", data_dir, "Br6471_Post_SFE.RDS"))
br6522_p <- readRDS(here("data", data_dir, "Br6522_Post_SFE.RDS"))
br8667_m <- readRDS(here("data", data_dir, "Br8667_Mid_SFE.RDS"))


filterCells <- function(sfe){
    cols_use <- names(colData(sfe))[str_detect(
        names(colData(sfe)), "_percent$")]
    cols_use <- cols_use[!grepl("anti", cols_use)]
    
    for (n in cols_use) {
        sfe <- get_neg_ctrl_outliers(n, sfe)
    }
    inds_keep <- sfe$nucleus_area <=200 & !sfe$is_blank_outlier &
        !sfe$is_depr_outlier &!sfe$is_negProbe_outlier &
        !sfe$is_negCodeword_outlier & sfe$nCounts > 0
    
    sfe <- sfe[,inds_keep]
    return(sfe)
}
br6471_p_filt <- filterCells(br6471_p)
br6522_p_filt <- filterCells(br6522_p)
br8667_m_filt <- filterCells(br8667_m)


saveRDS(br6471_p_filt, here("data", data_dir, "Br6471_Post_SFE_filt.RDS"))
saveRDS(br6522_p_filt, here("data", data_dir, "Br6522_Post_SFE_filt.RDS"))
saveRDS(br6522_p_filt , here("data", data_dir, "Br8667_Mid_SFE_filt.RDS"))


##### Slide 5548 #####

data_dir <- "output-XETG00089__0005548__Region_1__20230831__172339"

# Read in the three regions
br6471_p <- readRDS(here("data", data_dir, "Br6471_Post_SFE.RDS"))
br2743_m <- readRDS(here("data", data_dir, "Br2743_Mid_SFE.RDS"))
br8667_m <- readRDS(here("data", data_dir, "Br8667_Mid_SFE.RDS"))

br6471_p_filt <- filterCells(br6471_p)
br2743_m_filt <- filterCells(br2743_m)
br8667_m_filt <- filterCells(br8667_m)

saveRDS(br6471_p_filt, here("data", data_dir, "Br6471_Post_SFE_filt.RDS"))
saveRDS(br2743_m_filt, here("data", data_dir, "Br2743_Mid_SFE_filt.RDS"))
saveRDS(br8667_m_filt, here("data", data_dir, "Br8667_Mid_SFE_filt.RDS"))

