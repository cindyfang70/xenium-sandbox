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
# Xenium dataset. In our case, it's all of the regions cbind-ed into a large SFE
# so that we can pseudobulk across the donor.

sfe <- readRDS(here("processed-data/", "cindy", "slide-5434", "slide5434-filt_clustSFE.RDS"))
sfe <- sfe[,!grepl("Sample", sfe$region_id)]
colData(sfe) <- colData(sfe)[,unique(colnames(colData(sfe)))]
length(unique(sfe$clust_M1_lam0.9_k50_res0.4))
# [1] 9

table(sfe$clust_M1_lam0.9_k50_res0.4)
# 1     2     3     4     5     6     7     8     9    10    11 
# 43364 24864 15938 15914 14650 12467  9346  8259  3193  2563  1376 

sfe$donor_id <- unlist(lapply(strsplit(sfe$region_id, split="_"), '[', 1))
sfe$clust_M1_lam0.9_k50_res0.4 <- as.character(sfe$clust_M1_lam0.9_k50_res0.4)
sfe$clust_M1_lam0.9_k50_res0.4 <- paste0("Banksy", sfe$clust_M1_lam0.9_k50_res0.4)

# Perform the spatial registration
sfe_modeling_results <- registration_wrapper(
    sce = sfe,
    var_registration = "clust_M1_lam0.9_k50_res0.4",
    var_sample_id = "donor_id",
    gene_ensembl = "ID",
    gene_name = "Symbol"
)
# 2023-11-19 13:54:31.528475 make pseudobulk object
# 2023-11-19 13:54:33.190859 dropping 1 pseudo-bulked samples that are below 'min_ncells'.
# 2023-11-19 13:54:33.208148 drop lowly expressed genes
# 2023-11-19 13:54:33.2521 normalize expression
# 2023-11-19 13:54:33.2952 create model matrix
# 2023-11-19 13:54:33.324226 run duplicateCorrelation()
# 2023-11-19 13:54:33.43173 The estimated correlation is: 0.209745825848492
# 2023-11-19 13:54:33.432894 computing enrichment statistics
# 2023-11-19 13:54:33.816696 extract and reformat enrichment results
# 2023-11-19 13:54:33.824632 running the baseline pairwise model
# 2023-11-19 13:54:33.829013 computing pairwise statistics
# 2023-11-19 13:54:33.877568 computing F-statistics

# extract the t stats
registration_t_stats <- sfe_modeling_results$enrichment[, grep("^t_stat", colnames(sfe_modeling_results$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

#Step 3: correlate query statistics with layer reference
cor_layer <- layer_stat_cor(
    stats = registration_t_stats,
    modeling_results = layer_modeling_results,
    model_type = "enrichment",
    top_n = 100
)
lp <- layer_stat_cor_plot(cor_layer, max = max(cor_layer))


clustName <- "clust_M1_lam0.9_k50_res0.4"
# plot layers with correlation
colourCount = length(unique((colData(sfe)[[clustName]])))
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

plist <- list()
for (i in 1:length(unique(sfe$region_id))){
    region <- unique(sfe$region_id)[[i]]
    sub.sfe <- sfe[,sfe$region_id==region]
    
    p <- make_escheR(sub.sfe, y_reverse=FALSE) %>%
        add_fill(var=clustName)+
        scale_fill_manual(values=getPalette(colourCount))+
        ggtitle(region)
    plist[[i]] <- p
}

plist[[4]] <- lp
do.call(gridExtra::grid.arrange, c(plist,  ncol=3))
