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
    library(spatialLIBD)
})

#-------------------------------------------------------------------------------
# Trying out spatial registration following this tutorial:
# https://bioconductor.org/packages/devel/data/experiment/vignettes/spatialLIBD/inst/doc/guide_to_spatial_registration.html
# Overview:
# 1. Perform gene set enrichment analysis between spatial features 
# (ex. anatomical features, histological layers) on reference spatial data set. 
# Or access existing statistics.
# 2. Perform gene set enrichment analysis between features 
# (ex. new annotations, data-driven clusters) on new query data set.
# 3. Correlate the t-statistics between the reference and query features.
# 4. Annotate new spatial features with the most strongly associated reference feature.\
# 5. Plot correlation heat map to observe patterns between the two data sets.
#-------------------------------------------------------------------------------

# Step 1: Access existing statistics for the Visium DLPFC dataset with 
# manual annotations

layer_modeling_results <- fetch_data(type = "modeling_results")


# Step 2:  Perform gene set enrichment analysis between Banksy layers on the query
# Xenium dataset. For now, use just one of the datasets as an example.

sfe <- readRDS(here("processed-data/", "cindy", "slide-5434", "Br6471_Post_SFE_filt-with-banksy-domains.RDS"))
