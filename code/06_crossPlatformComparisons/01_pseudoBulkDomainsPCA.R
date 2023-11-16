suppressPackageStartupMessages({
    library(Voyager)
    library(SpatialFeatureExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(stringr)
    library(dplyr)
    library(here)
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

source(here("code", "cindy", "01_createSCE", "xenium_helpers.R"))
# read in data
args <- commandArgs(trailingOnly=TRUE)

sfe <- readRDS(args[[1]])
spe <- readRDS(args[[2]])
# sfe <- readRDS(here("processed-data", "cindy", "slide-5434", 
#                     "Br8667_Mid_SFE_filt.RDS"))
# 
# spe <- readRDS(here("Br8667_mid-visium-SPE.RDS"))


# pseudobulk across the domain labels
sfe.summed <- aggregateAcrossCells(sfe, id=colData(sfe)[,c("clust_M1_lam0.9_k50_res1.2")])
spe.summed <- aggregateAcrossCells(spe, id=colData(spe)[,c("BayesSpace_harmony_06")])

# convert the SFE to SPE so that the two datasets can be combined
sfe.summed <- SFEtoSPE(sfe.summed)

# subset the visium data to the same genes that are in xenium
spe.summed <- spe.summed[which(rowData(spe.summed)$gene_name %in% rownames(sfe.summed)),]
# rename the coldata column denoting the sample to match what's in the sfe
spe.summed$sample_id <- spe.summed$subject_position 

# Now subset the colData for both to have the same
colData(sfe.summed) <- colData(sfe.summed)[c("sample_id", "clust_M1_lam0.9_k50_res1.2")]
colData(spe.summed) <- colData(spe.summed)[c("sample_id", "BayesSpace_harmony_06")]

# rename the spatial coordinates to be the same 
col <- spatialCoords(spe.summed)[,1]
row <- spatialCoords(spe.summed)[,2]
spatialCoords(spe.summed) <- cbind(row=row, col=col)

col_xen <- spatialCoords(sfe.summed)[,2]
row_xen <- spatialCoords(sfe.summed)[,1]
spatialCoords(sfe.summed) <- cbind(row=row_xen, col=col_xen)

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

new.sfe$clust <- paste0(new.sfe$platform, new.sfe$clust_M1_lam0.9_k50_res1.2)
new.spe$clust <- paste0(new.spe$platform, new.spe$BayesSpace_harmony_06)
spes.all <- cbind(new.sfe, new.spe)

#spes.all <- logNormCounts(spes.all)

mat <- counts(spes.all)
colnames(mat) <- spes.all$clust

m <- cor(mat, method="pearson")
ComplexHeatmap::Heatmap(m)

m <- m[!grepl("Visium", rownames(m)),!grepl("Xenium", colnames(m))]


pdfname <- here("plots", "cindy", "06_crossPlatformComparisons0", 
                paste0(unique(sfe$region_id), 
                "Visium-Xenium-Correlation-Heatmap.pdf"))

pdf(pdfname)
ComplexHeatmap::Heatmap(m, name="Pearson Cor")
dev.off()

