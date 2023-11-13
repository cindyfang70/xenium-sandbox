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

# Now subset the colData for both to have the same
colData(sfe.summed) <- colData(sfe.summed)[c("sample_id", "clust_M1_lam0.9_k50_res1.2")]
colData(spe.summed) <- colData(spe.summed)[c("sample_id", "clust_M1_lam0.9_k50_res1.2")]

# rename the spatial coordinates to be the same 
col <- spatialCoords(spe.summed)[,1]
row <- spatialCoords(spe.summed)[,2]
spatialCoords(spe.summed) <- cbind(row=row, col=col)

col_xen <- spatialCoords(sfe.summed)[,2]
row_xen <- spatialCoords(sfe.summed)[,1]
spatialCoords(sfe.summed) <- cbind(row=row_xen, col=col_xen)

metadata(spe.summed) <- NULL
spe.all <- cbind(sfe.summed, spe.summed)

# create new objects to remove unnecessary information from the individual SPEs
new.spe <- SpatialExperiment(assay=list(counts=counts(spe.summed)), 
                                    spatialCoords=spatialCoords(spe.summed))
colData(new.spe) <- colData(spe.summed)

new.sfe <- SpatialExperiment(assays=list(counts=counts(sfe.summed)),
                             spatialCoords=spatialCoords(sfe.summed))
colData(new.sfe) <- colData(sfe.summed)

# add modality to the colData
new.sfe$platform <- "Xenium"
new.spe$platform <- "Visium"

spes.all <- cbind(new.sfe, new.spe)

spes.all <- scry::nullResiduals(spes.all, assay="counts", fam="poisson", type="pearson")
spes.all <- scater::runPCA(spes.all, ncomponents=10, 
                      ntop = 1000,
                      exprs_values = "poisson_pearson_residuals",
                      scale = TRUE, name = "GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())
plotReducedDim(spes.all,ncomponents=4, colour_by="platform", dimred="GLM-PCA")
