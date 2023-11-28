library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(Banksy)

# Run the Banksy multisample workflow:
# https://github.com/prabhakarlab/Banksy/blob/bioc/vignettes/multi-sample.Rmd

source(here("code","01_createSCE","xenium_helpers.R"))
    
sfe.paths1 <- system("ls processed-data/cindy/slide-5434", intern=TRUE)
sfe.paths1 <- sfe.paths1[grepl("*filt.RDS", sfe.paths1)]
sfe.paths1 <- paste0("processed-data/cindy/slide-5434/", sfe.paths1)

# read in the individual sfes and 
#aname <- "normcounts"
spe_list <- lapply(sfe.paths1, function(x){
    x <- readRDS(x)
    x <- SFEtoSPE(x) # note: this function only keeps counts assay
    #assay(x, aname) <- normalizeCounts(x, log = FALSE)
    x <- scry::nullResiduals(x,
                             assay="counts",
                             fam="poisson",
                             type="pearson") # try GLM PCA
    return(x) # return SPE for easier downstream manipulations
})


lambda <- 0.9
use_agf <- FALSE
k_geom <- 6

# Compute neighbourhood matrices individually
compute_agf <- FALSE
spe_list <- lapply(spe_list, computeBanksy, 
                   assay_name = "poisson_pearson_residuals", 
                   compute_agf = compute_agf, 
                   k_geom = k_geom)
# cbind the sfes together
spe_joint <- do.call(cbind, spe_list) 

assay(spe_joint, "joint_poisson_pearson_residuals") <- 
    assay(scry::nullResiduals(spe_joint, assay="counts",
                           fam="poisson",
                           type="pearson"), "poisson_pearson_residuals")

# PCA 
spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda,
                    group = "region_id",
                    assay_name="joint_poisson_pearson_residuals",
                    seed = 0858)

# UMAP
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, 
                          lambda = lambda, 
                          seed = 0922)

# Cluster
res <- 0.23
spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, 
                           lambda = lambda, resolution = res, seed = 0923)

cnm <- sprintf("clust_M0_lam%s_k50_res%s", lambda,
              res)


# colnames(colData(spe_joint))[grepl(cnm, colnames(colData(spe_joint)))] <- 
#     paste0("glmpca", cnm)

# cnm <- paste0("glmpca", cnm)

# Split the samples back up into their own SPEs
spe_list <- lapply(unique(spe_joint$region_id), function(x) 
    spe_joint[, spe_joint$region_id == x])

plist <- list()
for (i in seq_along(spe_list)){
    spe <- spe_list[[i]]
    p <- make_escheR(spe, y_reverse=FALSE) |>
        add_fill(cnm) +
        scale_fill_discrete()+
        ggtitle(unique(spe$region_id))
    plist[[i]] <- p
}


###############################################
set.seed(1015)
spe_joint <- scry::nullResiduals(spe_joint, assay="counts", fam="poisson", type="pearson")
spe_joint <- scater::runPCA(spe_joint, ncomponents=50, 
                            ntop = 1000,
                            exprs_values = "poisson_pearson_residuals",
                            scale = TRUE, name = "GLM-PCA",
                            BSPARAM = BiocSingular::RandomParam())

spe_joint <- scater::runUMAP(spe_joint, dimred="GLM-PCA")

pdf(here("plots", "05_segmentRegions", "banksy", 
         "multisample-slide-5434.pdf"), height=15, width=20)
do.call(gridExtra::grid.arrange, c(plist, ncol=2))


plotReducedDim(spe_joint,ncomponents=4, colour_by="region_id", dimred="GLM-PCA")

plotUMAP(spe_joint, colour_by="region_id") + plotUMAP(spe_joint, colour_by=cnm)


spe_joint <- schex::make_hexbin(spe_joint, nbins = 40, 
                           dimension_reduction = "UMAP", use_dims=c(1,2))
schex::plot_hexbin_meta(spe_joint, col="region_id", action="majority")
# schex::plot_hexbin_meta_plus(spe_joint,
#                          col1=cnm,
#                          col2="region_id", action="prop")
schex::plot_hexbin_meta(spe_joint, col="total", action="median")
schex::plot_hexbin_meta(spe_joint, col="subsets_any_neg_sum", action="median")
schex::plot_hexbin_meta(spe_joint, col="nucleus_area", action="median")
dev.off()