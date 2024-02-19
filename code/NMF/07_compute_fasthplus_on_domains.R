if (!require("fasthplus")) devtools::install_github(repo="ntdyjack/fasthplus", ref = "main")
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
    library(escheR)
    library(RColorBrewer)
    library(data.table)
    library(fasthplus)
})
#-------------------------------------------------------------------------------#
# Compute fast h plus metric on the Banksy cell type clusters to evaluate clustering
# results
# reference: https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/06_fasthplus/01_fasthplus.R
#-------------------------------------------------------------------------------#
# read in the SPE object
source(here("code", "cindy", "01_createSCE", "xenium_helpers.R"))
banksy_sfe <- readRDS(here('processed-data', "cindy", "all-tissues-sfe-with-Banksy-spatial-domains.RDS"))

sfe_5434_nmf <- readRDS(here("processed-data", "cindy", "slide-5434","slide-5434-spe-with-nmf-k200-bayesspace.RDS"))
sfe_5548_nmf <- readRDS(here("processed-data", "cindy", "slide-5548", "slide-5548-spe-with-nmf-k200-bayesspace.RDS"))

sfe_5434_nmf <- SFEtoSPE(sfe_5434_nmf)
sfe_5548_nmf <- SFEtoSPE(sfe_5548_nmf)

# keep only the colData columns that are in both
colData(sfe_5548_nmf) <- colData(sfe_5548_nmf)[colnames(colData(sfe_5434_nmf))]

nmf_sfe <- cbind(sfe_5434_nmf, sfe_5548_nmf)

# make unique cell ids
nmf_sfe$cell_sample_id <- paste(nmf_sfe$cell_id, nmf_sfe$region_id, sep="-")
banksy_sfe$cell_sample_id <- paste(banksy_sfe$cell_id, banksy_sfe$region_id, sep="-")

# make the ordering of cells match in the two separate spe objects
spe <- nmf_sfe[,match(banksy_sfe$cell_sample_id, nmf_sfe$cell_sample_id)]
sum(spe$cell_sample_id == banksy_sfe$cell_sample_id) == ncol(spe) # sanity check

reducedDim(spe, "UMAP") <- reducedDim(banksy_sfe, "UMAP_M1_lam0.9")

spe$banksy_domains <- banksy_sfe$clust_M1_lam0.9_k50_res0.45

# proportion of each nmf label in each banksy domain
df <- colData(spe)[,c("banksy_domains", "preds_from_bayesspace_NMF_k200")]
prop.table(table(df), 1)

spe$nmf_preds <- spe$preds_from_bayesspace_NMF_k200

lambda <- 0.9
res <- 0.45


clustName <- "banksy_domains"
k <- length(unique(colData(spe)[[clustName]]))

# hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations

find_t <- function(L, proportion = 0.05) {
    initial_t <- floor(length(L) * proportion)
    smallest_cluster_size <- min(table(L))
    n_labels <- length(unique(L))
    ifelse(smallest_cluster_size > (initial_t / n_labels), initial_t, smallest_cluster_size * n_labels)
}

initial_t <- find_t(L = colData(spe)[[clustName]], proportion = 0.01)

cluster_prop <- table(colData(spe)[[clustName]]) / ncol(spe)
bad_clusters <- which(cluster_prop < 0.01 / k)
if (length(bad_clusters) > 0) {
    message("For k: ", k, " we are dropping small clusters: ", paste(names(bad_clusters), collapse = ", "))
    spe <- spe[, !colData(spe)[[clustName]] %in% as.integer(names(bad_clusters))]
    updated_t <- find_t(colData(spe)[[clustName]], 0.01)
    message("initial t: ", initial_t, "; updated t: ", updated_t)
}else{
    updated_t <- initial_t
}


set.seed(1211)
fasthplus <- hpb(D = reducedDims(spe)$"UMAP", L = colData(spe)[[clustName]], t = updated_t, r = 30)
results <- data.frame(k = k, fasthplus = fasthplus)
print(results)
write.table(results, file = here::here("processed-data", "cindy", "NMF", "fasthplus_results_banksy_domains.csv"), append = TRUE)
