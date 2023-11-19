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
#------------------------------------------------------------#
# Run banksy on one whole slide at a time
# tutorial: https://github.com/prabhakarlab/Banksy/tree/bioc
#------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
sfe <- readRDS(args[[1]])
lambda <- as.numeric(args[[2]])
k <- args[[3]]
res <- as.numeric(args[[4]])


print(lambda)
print(res)


sfe <- computeLibraryFactors(sfe)
aname <- "normcounts"
assay(sfe, aname) <- normalizeCounts(sfe, log = FALSE)


# Compute neighbourhood matrices
k_geom <- c(15, 30)

sfe <- Banksy::computeBanksy(sfe, assay_name = aname, compute_agf = TRUE,
                             k_geom = k_geom)

# run PCA and UMAP, then Leiden clustering to find spatial domains
# lambda is in [0,1], higher values of lambda puts more weight on spatial
# information
set.seed(1000)
sfe <- Banksy::runBanksyPCA(sfe, use_agf = TRUE, lambda = lambda)
print("PCA done")
sfe <- Banksy::runBanksyUMAP(sfe, use_agf = TRUE, lambda = lambda)
print("UMAP done")

# for leiden, higher resolution = more clusters, 
# lower resolution = fewer clusters
sfe <- Banksy::clusterBanksy(sfe, use_agf = TRUE, lambda = lambda, 
                             resolution = res)
print("clustering done")

# save the Banksy results
saveRDS(sfe, args[[1]])

# default for the leiden algorithm in clusterBanksy is k_neighbours=50
clustName <- sprintf("clust_M1_lam%s_k%s_res%s", lambda, 50, res)


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
slide <- unlist(strsplit(sfe$region_id[[1]], split="_"))[[3]]
pdfname <- paste0(sprintf("banksy-wholeslide-%s", slide), ".pdf")

pdf(here("plots", "cindy", "05_segmentRegions", "banksy", pdfname), height=15, width=20)
    do.call(gridExtra::grid.arrange, c(plist, ncol=3))
dev.off()
