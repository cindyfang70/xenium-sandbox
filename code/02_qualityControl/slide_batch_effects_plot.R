suppressPackageStartupMessages({
    library(Voyager)
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
# Make some PCA plots to look at batch effects between slides #
######################################

source(here("code", "01_createSCE", "xenium_helpers.R"))

# Read in the data
data_dir_5434 <- "slide-5434"
data_dir_5548 <- "slide-5548"

sfe1 <- readRDS(here("processed-data", "cindy", data_dir_5434, "slide5434-filt_clustSFE.RDS"))
sfe2 <- readRDS(here("processed-data", "cindy", data_dir_5548, "slide5548-filt_clustSFE.RDS"))

sfe1 <- sfe1[,sfe1$region_id != "Sample"]
sfe2 <- sfe2[,sfe2$region_id != "Sample"]

# Plot UMAP
set.seed(3557)
sfe1 <- scater::runUMAP(sfe1, dimred="GLM-PCA")
plotUMAP(sfe1, colour_by="region_id")+
    ggtitle(data_dir_5434)

sfe2 <- scater::runUMAP(sfe2, dimred="GLM-PCA")
plotUMAP(sfe2, colour_by="region_id")+
    ggtitle(data_dir_5548)

# Make hexbins
sfe1 <- schex::make_hexbin(sfe1, dimension_reduction="UMAP", nbins=40)
schex::plot_hexbin_meta(sfe1, col="region_id", action="majority")

sfe2 <- schex::make_hexbin(sfe2, dimension_reduction="UMAP", nbins=40)
schex::plot_hexbin_meta(sfe2, col="region_id", action="majority")

br8667 <- SFEtoSPE(sfe1[,grepl("Br8667", sfe1$region_id)])


br8667_2 <- SFEtoSPE(sfe2[,grepl("Br8667", sfe2$region_id)])
colData(br8667_2) <- colData(br8667_2)[,colnames(colData(br8667))]


br8667_all <- cbind(br8667, br8667_2)

br8667_all <- scry::nullResiduals(br8667_all, assay="counts", fam="poisson", type="pearson")
br8667_all <- scater::runPCA(br8667_all, ncomponents=50, 
                                    ntop = 1000,
                                    exprs_values = "poisson_pearson_residuals",
                                    scale = TRUE, name = "GLM-PCA",
                                    BSPARAM = BiocSingular::RandomParam())

br8667_all <- scater::runUMAP(br8667_all, dimred="GLM-PCA")

plotUMAP(br8667_all, colour_by="region_id")
br8667_all <- schex::make_hexbin(br8667_all, dimension_reduction="UMAP", nbins=40)

schex::plot_hexbin_meta(br8667_all, col="region_id", action="majority")+
    ggtitle("Br8667_Mid")


br6471_1 <- SFEtoSPE(sfe1[,grepl("Br6471",sfe1$region_id)])
br6471_2 <- SFEtoSPE(sfe2[,grepl("Br6471",sfe2$region_id)])

colData(br6471_2) <- colData(br6471_2)[,colnames(colData(br6471_1))]
br6471_all <- cbind(br6471_1, br6471_2)

br6471_all <- scry::nullResiduals(br6471_all, assay="counts", fam="poisson", type="pearson")
br6471_all <-  scater::runPCA(br6471_all, ncomponents=50, 
                              ntop = 1000,
                              exprs_values = "poisson_pearson_residuals",
                              scale = TRUE, name = "GLM-PCA",
                              BSPARAM = BiocSingular::RandomParam())

br6471_all <- scater::runUMAP(br6471_all, dimred="GLM-PCA")

plotUMAP(br6471_all, colour_by="region_id")
br6471_all <- schex::make_hexbin(br6471_all, dimension_reduction="UMAP", nbins=40)

schex::plot_hexbin_meta(br6471_all, col="region_id", action="majority")+
    ggtitle("Br6471_Post")


#plotUMAP(sfe1, colour_by="region_id")


## TODO: run multisample banksy and plot labels on UMAP

