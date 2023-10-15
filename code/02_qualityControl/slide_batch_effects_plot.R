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
})

######################################
# Make some PCA plots to look at batch effects between slides#
######################################

source(here("code", "xenium_helpers.R"))

# Read in the data
data_dir_5434 <- "output-XETG00089__0005434__Region_1__20230831__172339"
data_dir_5548 <- "output-XETG00089__0005548__Region_1__20230831__172339"

sfe1 <- readRDS(here("data", data_dir_5434, "Br6471_Post_SFE.RDS"))
sfe2 <- readRDS(here("data", data_dir_5548, "Br6471_Post_SFE.RDS"))

# Add slide number to the SFEs
sfe1$sample_id <- "slide_5434"
sfe2$sample_id <- "slide_5548"

# Convert to SPE and combine 
spe1 <- SFEtoSPE(sfe1)
spe2 <- SFEtoSPE(sfe2)

spe <- cbind(spe1, spe2)


#### Quick QC
cols_use <- names(colData(spe))[str_detect(names(colData(spe)), "_percent$")]
cols_use <- cols_use[!grepl("anti", cols_use)] # no antisense controls detected, so avoid dividing by zero


cols_use2 <- names(colData(spe))[str_detect(names(colData(spe)), "_detected$")]
cols_use2 <- cols_use2[!grepl("anti", cols_use)] 


for (n in cols_use) {
    spe <- get_neg_ctrl_outliers(n, spe)
}


inds_keep <- spe$nCounts > 0  & !spe$is_blank_outlier & !spe$is_negCodeword_outlier & !spe$is_negProbe_outlier & !spe$is_depr_outlier
spe <- spe[,inds_keep] # 8342 cells were discarded

### normalize and PCA
spe <- logNormCounts(spe)
spe <- scater::runPCA(spe, ncomponents=30, scale=TRUE,
                      exprs_values="logcounts")

p<-plotReducedDim(spe, ncomponents=4, dimred="PCA", colour_by="sample_id")

pdf("pca_plot.pdf")
print(p)
dev.off()

spe <- schex::make_hexbin(spe, nbins = 25, 
                        dimension_reduction = "PCA", use_dims=c(1,2))

spe$sample_id_num <- ifelse(spe$sample_id=="slide_5434", 1, 0)
p <- schex::plot_hexbin_meta(spe, col="sample_id",
                             action="prop")

###### Now let's try this for each slide #####

sfe5434 <- readRDS(here("data", data_dir_5434, "xenium-0005434-SFE.RDS"))
#sfe5548 <- readRDS(here("data", data_dir_5548, "xenium-0005548-SFE.RDS"))

cols_use <- names(colData(sfe5434))[str_detect(names(colData(sfe5434)), "_percent$")]
cols_use <- cols_use[!grepl("anti", cols_use)] 

spe5434 <- SFEtoSPE(sfe5434)
for (n in cols_use) {
    spe5434 <- get_neg_ctrl_outliers(n, spe5434)
}


inds_keep <- spe5434$nCounts > 0  & !spe5434$is_blank_outlier & 
    !spe5434$is_negCodeword_outlier & !spe5434$is_negProbe_outlier & 
    !spe5434$is_depr_outlier

spe5434 <- spe5434[,inds_keep]

spe5434 <- logNormCounts(spe5434)
spe5434 <- scater::runPCA(spe5434, ncomponents=30, scale=TRUE,
                      exprs_values="logcounts")

spe5434 <- schex::make_hexbin(spe5434, nbins = 25, 
                          dimension_reduction = "PCA", use_dims=c(1,2))


p1 <- schex::plot_hexbin_meta(spe5434, col="region_id",
                             action="prop", no=1)

p2 <- schex::plot_hexbin_meta(spe5434, col="region_id",
                              action="prop", no=2)

p3 <- schex::plot_hexbin_meta(spe5434, col="region_id",
                              action="prop", no=3)
p1+p2+p3

sfe5548 <- readRDS(here("data", data_dir_5548, "xenium-0005548-SFE.RDS"))
spe5548 <- SFEtoSPE(sfe5548)

for (n in cols_use) {
    spe5548 <- get_neg_ctrl_outliers(n, spe5548)
}


inds_keep <- spe5548$nCounts > 0  & !spe5548$is_blank_outlier & 
    !spe5548$is_negCodeword_outlier & !spe5548$is_negProbe_outlier & 
    !spe5548$is_depr_outlier

spe5548 <- spe5548[,inds_keep]

spe5548 <- logNormCounts(spe5548)
spe5548 <- scater::runPCA(spe5548, ncomponents=30, scale=TRUE,
                          exprs_values="logcounts")

spe5548 <- schex::make_hexbin(spe5548, nbins = 25, 
                              dimension_reduction = "PCA", use_dims=c(1,2))


pl1 <- schex::plot_hexbin_meta(spe5548, col="region_id",
                              action="prop", no=1)

pl2 <- schex::plot_hexbin_meta(spe5548, col="region_id",
                              action="prop", no=2)

pl3 <- schex::plot_hexbin_meta(spe5548, col="region_id",
                              action="prop", no=3)


pdf(here("plots", "xenium_plots", "batch_effect_plots.PDF"),
    height=5, width=12)
print(p1+p2+p3)
schex::plot_hexbin_meta(spe5434, col="nucleus_area",
                        action="mean") + 
    schex::plot_hexbin_meta(spe5434, col="total_counts",
                        action="mean")
print(pl1 + pl2 + pl3)
schex::plot_hexbin_meta(spe5548, col="nucleus_area",
                            action="mean") +
    schex::plot_hexbin_meta(spe5548, col="total_counts",
                        action="mean")
dev.off()
