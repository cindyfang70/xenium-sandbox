suppressPackageStartupMessages({
    library(Voyager)
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
    library(sf)
    library(BiocParallel)
    library(dplyr)
    library(here)
})

#########################################################
# Whole slide clustering in the tissue using GLM-PCA    #
# normalization and then SNN graph-based clustering     #
#########################################################
# Read in the data
args <- commandArgs(trailingOnly=TRUE)
sfe <- readRDS(args[[1]])
config <- read.table(args[[2]], header=TRUE)
k <- args[[length(args)]] 
region_id <- config$Slide[[1]]
sub.sfes.paths <- paste0(here("processed-data", "cindy", paste0("slide-",region_id)), 
                   "/", config$SampleName, "_SFE_filt.RDS")

sub.sfes <- lapply(sub.sfes.paths, readRDS)

# QC the whole slide
cols_use <- names(colData(sfe))[str_detect(
    names(colData(sfe)), "_percent$")]
cols_use <- cols_use[!grepl("anti", cols_use)]

get_neg_ctrl_outliers <- function(col, spe) {
    inds <- colData(spe)$nCounts > 0 & colData(spe)[[col]] > 0
    df <- colData(spe)[inds,]
    outlier_inds <- isOutlier(df[[col]], type = "higher")
    outliers <- rownames(df)[outlier_inds]
    col2 <- str_remove(col, "^subsets_")
    col2 <- str_remove(col2, "_percent$")
    new_colname <- paste("is", col2, "outlier", sep = "_")
    colData(spe)[[new_colname]] <- colnames(spe) %in% outliers
    return(spe)
}
    
for (n in cols_use) {
    sfe <- get_neg_ctrl_outliers(n, sfe)
}
inds_keep <- sfe$nucleus_area <=200 & !sfe$is_blank_outlier &
        !sfe$is_depr_outlier &!sfe$is_negProbe_outlier &
        !sfe$is_negCodeword_outlier & sfe$nCounts > 0
    
sfe <- sfe[,inds_keep]
    
# remove the negative controls
sfe <- sfe[which(rowData(sfe)$Type=="Gene Expression"), ]

# GLM-PCA
sfe <- scry::nullResiduals(sfe, assay="counts", fam="poisson", type="pearson")
print('null residuals done')
sfe <- scater::runPCA(sfe, ncomponents=50, 
                      ntop = 1000,
                      exprs_values = "poisson_pearson_residuals",
                      scale = TRUE, name = "GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

# graph based clustering
set.seed(10222023)
g <- scran::buildSNNGraph(sfe, k=as.numeric(k), use.dimred = "GLM-PCA")
lou <- igraph::cluster_louvain(g)

# add the clusters back into the SFE object and save it
clusterName <- paste0("louWholeSlide",k)
colData(sfe)[[clusterName]] <- paste0("LouWholeSlide", lou$membership)


clustSfeName <- paste0("slide", region_id, "-filt_clustSFE.RDS")
saveRDS(sfe, here("processed-data", "cindy", paste0("slide-", region_id), clustSfeName))

# now put the cluster labels into the individual tissue SFE objects
for (i in 1:length(sub.sfes)){
    sub.sfe <- sub.sfes[[i]]
    whole.subset <- sfe[,which(sfe$cell_id %in% sub.sfe$cell_id)]
    print(any(!(whole.subset$cell_id == sub.sfe$cell_id)))
    colData(sub.sfe)[[clusterName]] <- colData(whole.subset)[[clusterName]]
    
    saveRDS(sub.sfe, sub.sfes.paths[[i]])
}


