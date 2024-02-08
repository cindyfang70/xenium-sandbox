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
    library(data.table)
    library(pals)
    library(alphahull)
})
#------------------------------------------------------------#
# Run banksy on all 6 tissue regions at the ame time
# tutorial: https://github.com/prabhakarlab/Banksy/tree/bioc
#------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
source(here("code", "01_createSCE", "xenium_helpers.R"))

# read in the SPE
sfename <- here("processed-data", "cindy", "all-tissues-spe-with-banksy-cell-types.RDS")
spe <- readRDS(sfename)

tab <- prop.table(table(spe$clust_M1_lam0_k50_res1, spe$preds_from_bayesspace_NMF_k200),2)
tab

sub.spe <- spe[,which(spe$region_id == "Br6471_Post_5434")]
for (i in 1:length(unique(sub.spe$preds_from_bayesspace_NMF_k200))){
    layer <- unique(sub.spe$preds_from_bayesspace_NMF_k200)[[i]]
    sub.spe.layer <- sub.spe[,which(sub.spe$preds_from_bayesspace_NMF_k200==layer)]
    
    alpha.shape <- ashape(x=spatialCoords(sub.spe.layer)[,1], y=spatialCoords(sub.spe.layer)[,2], alpha=100)
    print(alpha.shape)
    
    
}
