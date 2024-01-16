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
    library(harmony)
})
#------------------------------------------------------------#
# Run banksy on one whole slide at a time
# tutorial: https://github.com/prabhakarlab/Banksy/tree/bioc
#------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)

args <- c("processed-data/cindy/slide-5434/slide5434-filt_clustSFE.RDS", 
          "0.9", "50", "0.45")
sfe <- readRDS(args[[1]])
lambda <- as.numeric(args[[2]])
k <- args[[3]]
res <- as.numeric(args[[4]])
harmony_lam <- as.numeric(args[[5]])

slide <- unlist(strsplit(sfe$region_id[[1]], split="_"))[[3]]

colData(sfe) <- colData(sfe)[,!grepl("clust", colnames(colData(sfe)))]

sfe <- computeLibraryFactors(sfe)
aname <- "normcounts"
assay(sfe, aname) <- normalizeCounts(sfe, log = FALSE)


# Compute neighbourhood matrices
k_geom <- c(15, 30)

sfe <- Banksy::computeBanksy(sfe, assay_name = aname, compute_agf = TRUE,
                             k_geom = k_geom)

# run PCA
# lambda is in [0,1], higher values of lambda puts more weight on spatial
# information
set.seed(1000)
sfe <- Banksy::runBanksyPCA(sfe, use_agf = TRUE, lambda = lambda)
print("PCA done")


# Run Harmony batch effect correction on the Banksy matrix
set.seed(1106)
reducedDimName <- sprintf("PCA_M1_lam%s", lambda)
harmony_embedding <- RunHarmony(
    data_mat = reducedDim(sfe, reducedDimName),
    meta_data = colData(sfe),
    vars_use = c("region_id"),
    do_pca = FALSE,
    verbose = TRUE,
    lambda=harmony_lam,
    max_iter=100,
    kmeans_init_nstart=20, 
    kmeans_init_iter_max=5000
)

reducedDim(sfe, "Harmony_BANKSY") <- harmony_embedding

# run UMAP on both the raw Banksy matrix and the Harmony corrected one
sfe <- Banksy::runBanksyUMAP(sfe, use_agf = TRUE, lambda = lambda) # raw
sfe <- runBanksyUMAP(sfe, dimred = "Harmony_BANKSY") # Harmony corrected

# Visualize the UMAPs annotated by subject ID:
umapName <- sprintf("UMAP_M1_lam%s", lambda)
pdf(here("plots", "cindy", "05_segmentRegions", "banksy", 
     sprintf("slide%s-wholeSlide-Harmony-corrected-lambda%s-harmony-lam%s.pdf", slide, lambda, harmony_lam)))
cowplot::plot_grid(
    scater::plotReducedDim(sfe, umapName, 
                   point_size = 0.1,
                   point_alpha = 0.5,
                   color_by = "region_id") +
        theme(legend.position = "none"),
    scater::plotReducedDim(sfe, "UMAP_Harmony_BANKSY", 
                   point_size = 0.1,
                   point_alpha = 0.5,
                   color_by = "region_id") +
        theme(legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))),
    nrow = 1,
    rel_widths = c(1, 1.2)
)
dev.off()

# for leiden, higher resolution = more clusters, 
# lower resolution = fewer clusters
sfe <- Banksy::clusterBanksy(sfe, use_agf = TRUE, lambda = lambda, 
                             resolution = res, dimred = "Harmony_BANKSY")
print("clustering done")

# save the Banksy results
#saveRDS(sfe, args[[1]])

# default for the leiden algorithm in clusterBanksy is k_neighbours=50
clustName <- sprintf("clust_Harmony_BANKSY_k%s_res%s",  50, res)


# plot using escheR 
colourCount = nlevels(colData(sfe)[[clustName]])
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

plist <- list()
for (i in 1:length(unique(sfe$region_id))){
    region <- unique(sfe$region_id)[[i]]
    sub.sfe <- sfe[,sfe$region_id==region]
    
    p <- make_escheR(sub.sfe, y_reverse=FALSE) %>%
            add_fill(var=clustName)+
            scale_fill_manual(values=getPalette(colourCount))
    plist[[i]] <- p
}

pdfname <- paste0(sprintf("banksy-wholeslide-harmony-%s-res%s-harmony-lambda%s",
                          slide, res, harmony_lam), ".pdf")

pdf(here("plots", "cindy", "05_segmentRegions", "banksy", pdfname), height=15, width=20)
    do.call(gridExtra::grid.arrange, c(plist, ncol=3))
dev.off()
