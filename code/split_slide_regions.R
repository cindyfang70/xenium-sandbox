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
    library(arrow)
    library(sf)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(BiocParallel)
    library(dplyr)
    library(here)
})
# load helper functions
source(here("code/xenium_helpers.R"))


#### Slide Number 5548 ####
data_dir <- "output-XETG00089__0005548__Region_1__20230831__172339/"
sfe <- readXenium(data_dir)

# Assign one of three region names to each cell, first create a dummy coldata variable.
sfe$region_id <- "Sample"

# Br2743_Mid
br2743_inds <- getRegionInds("Br2743_Mid", data_dir, sfe)
sfe$region_id[br2743_inds] <- "Br2743_Mid_5548"

# Br6471_Post
br6471_inds <- getRegionInds("Br6471_Post",data_dir, sfe)
sfe$region_id[br6471_inds] <- "Br6471_Post_5548"

# Br8667_Mid
br8667_inds <- getRegionInds("Br8667_Mid", data_dir, sfe)
sfe$region_id[br8667_inds] <- "Br8667_Mid_5548"

# Check for any cells that weren't assigned to a region
sum(sfe$region_id == "Sample") # 6 are not assigned. in the morphology image, it's hard to tell which cells belong to which region..

# Create the individual SFE objects for each region and save it as RDS
br2743 <- sfe[,br2743_inds]
saveRDS(br2743, here("data", data_dir, "Br2743_Mid_SFE.RDS"))

br6471 <- sfe[,br6471_inds]
saveRDS(br6471, here("data", data_dir, "Br6471_Post_SFE.RDS"))

br8667 <- sfe[,br8667_inds]
saveRDS(br8667, here("data", data_dir, "Br8667_Mid_SFE.RDS"))

# Save the big SFE 
saveRDS(sfe, here("data", data_dir, "xenium-0005548-SFE.RDS"))


#### Slide Number 5434 ####
data_dir <- "output-XETG00089__0005434__Region_1__20230831__172339"
sfe <- readXenium(data_dir)

# Assign one of three region names to each cell, first create a dummy coldata variable.
sfe$region_id <- "Sample"

# Br6522_Post
br6522_inds <- getRegionInds("Br6522_Post", data_dir, sfe)
sfe$region_id[br6522_inds] <- "Br6522_Post_5434"

# Br6471_Post
br6471_inds <- getRegionInds("Br6471_Post", data_dir, sfe)
sfe$region_id[br6471_inds] <- "Br6471_Post_5434"

# Br8667_Mid
br8667_inds <- getRegionInds("Br8667_Mid", data_dir, sfe)
sfe$region_id[br8667_inds] <- "Br8667_Mid_5434"

# Check for any cells that weren't assigned to a region
sum(sfe$region_id == "Sample") # 5 are not assigned. in the morphology image, it's hard to tell which cells belong to which region..

# Create the individual SFE objects for each region and save it as RDS
br6522 <- sfe[,br6522_inds]
saveRDS(br6522, here("data", data_dir, "Br6522_Post_SFE.RDS"))

br6471 <- sfe[,br6471_inds]
saveRDS(br6471, here("data", data_dir, "Br6471_Post_SFE.RDS"))

br8667 <- sfe[,br8667_inds]
saveRDS(br8667, here("data", data_dir, "Br8667_Mid_SFE.RDS"))

# Save the big SFE 
saveRDS(sfe, here("data", data_dir, "xenium-0005434-SFE.RDS"))
