suppressPackageStartupMessages({
    library(SpatialFeatureExperiment)
    library(SingleCellExperiment)
    library(Banksy)
    library(here)
    library(escheR)
    library(ggplot2)
})

# Run the Banksy multisample workflow:
# https://github.com/prabhakarlab/Banksy/blob/bioc/vignettes/multi-sample.Rmd

source(here("code","cindy", "01_createSCE","xenium_helpers.R"))

# sfe.paths1 <- system("ls processed-data/cindy/slide-5434", intern=TRUE)
# sfe.paths1 <- sfe.paths1[grepl("*filt.RDS", sfe.paths1)]
# sfe.paths1 <- paste0("processed-data/cindy/slide-5434/", sfe.paths1)

args <- commandArgs(trailingOnly = TRUE)
print(args)
# get the SFE paths
sfe.paths <- args[1:3]

# get the Banksy params
lambda <- as.numeric(args[[4]])
k_geom <- as.numeric(args[[5]])
res <- as.numeric(args[[6]])

# get the slide number
slide <- args[[7]]

# read in the individual sfes and 
aname <- "normcounts"
spe_list <- lapply(sfe.paths, function(x){
    x <- readRDS(x)
    x <- SFEtoSPE(x) # note: this function only keeps counts assay
    assay(x, aname) <- normalizeCounts(x, log = FALSE)
    return(x) # return SPE for easier downstream manipulations
})

use_agf <- FALSE

# Compute neighbourhood matrices individually
compute_agf <- FALSE
spe_list <- lapply(spe_list, computeBanksy, 
                   assay_name = aname,
                   compute_agf = compute_agf, 
                   k_geom = k_geom)
# cbind the sfes together
spe_joint <- do.call(cbind, spe_list) 


# PCA 
spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda,
                    group = "region_id",
                    seed = 0858)

# UMAP
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, 
                          lambda = lambda, 
                          seed = 0922)

# Cluster

spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, 
                           lambda = lambda, resolution = res, seed = 0923)

cnm <- sprintf("clust_M0_lam%s_k50_res%s", lambda,
              res)


# Split the samples back up into their own SPEs
spe_list <- lapply(unique(spe_joint$region_id), function(x) 
    spe_joint[, spe_joint$region_id == x])

print(length(spe_list))

plist <- list()
for (i in seq_along(spe_list)){
    spe <- spe_list[[i]]
    p <- make_escheR(spe, y_reverse=FALSE) |>
        add_fill(cnm) +
        scale_fill_discrete()+
        ggtitle(unique(spe$region_id))
    plist[[i]] <- p
}

pdf(here("plots", "cindy", "05_segmentRegions", "banksy",
         sprintf("slide%s-multisample-banksy-plots-res%s.pdf", slide, res)), height=15, width=20)
do.call(gridExtra::grid.arrange, c(plist, ncol=2))
dev.off()

saveRDS(spe_joint, here("processed-data", "cindy", 
                        sprintf("slide%s-all-samples-spe-with-banksy.RDS",slide)))

###############################################
# set.seed(1015)
# spe_joint <- scry::nullResiduals(spe_joint, assay="counts", fam="poisson", type="pearson")
# spe_joint <- scater::runPCA(spe_joint, ncomponents=50, 
#                             ntop = 1000,
#                             exprs_values = "poisson_pearson_residuals",
#                             scale = TRUE, name = "GLM-PCA",
#                             BSPARAM = BiocSingular::RandomParam())
# 
# spe_joint <- scater::runUMAP(spe_joint, dimred="GLM-PCA")
# 
# pdf(here("plots", "05_segmentRegions", "banksy", 
#          "multisample-slide-5434.pdf"), height=15, width=20)
# do.call(gridExtra::grid.arrange, c(plist, ncol=2))
# 
# 
# plotReducedDim(spe_joint,ncomponents=4, colour_by="region_id", dimred="GLM-PCA")
# 
# plotUMAP(spe_joint, colour_by="region_id") + plotUMAP(spe_joint, colour_by=cnm)
# 
# 
# spe_joint <- schex::make_hexbin(spe_joint, nbins = 40, 
#                            dimension_reduction = "UMAP", use_dims=c(1,2))
# schex::plot_hexbin_meta(spe_joint, col="region_id", action="majority")
# # schex::plot_hexbin_meta_plus(spe_joint,
# #                          col1=cnm,
# #                          col2="region_id", action="prop")
# schex::plot_hexbin_meta(spe_joint, col="total", action="median")
# schex::plot_hexbin_meta(spe_joint, col="subsets_any_neg_sum", action="median")
# schex::plot_hexbin_meta(spe_joint, col="nucleus_area", action="median")
# dev.off()