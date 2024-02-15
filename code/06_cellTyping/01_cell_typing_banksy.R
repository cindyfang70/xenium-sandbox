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
})
#------------------------------------------------------------#
# Run banksy on all 6 tissue regions at the ame time
# tutorial: https://github.com/prabhakarlab/Banksy/tree/bioc
#------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
source(here("code", "01_createSCE", "xenium_helpers.R"))

args <- c("processed-data/cindy/slide-5434/slide-5434-spe-with-nmf-k200-bayesspace.RDS",
          "processed-data/cindy/slide-5548/slide-5548-spe-with-nmf-k200-bayesspace.RDS",
          0, 50, 1.0)
sfe1 <- readRDS(args[[1]])
sfe2 <- readRDS(args[[2]])
lambda <- as.numeric(args[[3]])
k <- args[[4]]
res <- as.numeric(args[[5]])


print(lambda)
print(res)

spe1 <- SFEtoSPE(sfe1)
spe2 <- SFEtoSPE(sfe2)

colData(spe1) <- colData(spe1)[,!grepl("clust", colnames(colData(spe1)))]
spe1 <- spe1[,!grepl("Sample", spe1$region_id)]

colData(spe2) <- colData(spe2)[,!grepl("clust", colnames(colData(spe2)))]
spe2 <- spe2[,!grepl("Sample", spe2$region_id)]

spe <- cbind(spe1, spe2)

locs <- spatialCoords(spe)
locs <- cbind(locs, sample_id = factor(spe$sample_id))
locs_dt <- data.table(locs)

colnames(locs_dt) <- c("sdimx", "sdimy", "group")
locs_dt[, sdimx := sdimx - min(sdimx), by = group]

global_max <- max(locs_dt$sdimx) * 1.5
locs_dt[, sdimx := sdimx + group * global_max]

locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(spe)

spatialCoords(spe) <- locs

sfe <- spe

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
sfename <- here("processed-data", "cindy", "all-tissues-spe-with-banksy-cell-types.RDS")
saveRDS(sfe, sfename)

# default for the leiden algorithm in clusterBanksy is k_neighbours=50
clustName <- sprintf("clust_M1_lam%s_k%s_res%s", lambda, 50, res)


# plot using escheR 
colourCount = nlevels(colData(sfe)[[clustName]])
getPalette = colorRampPalette(brewer.pal(colourCount, "Set3"))

myPal <- brewer.pal(12, "Set3")

myPal <- c(myPal, "#83b9c9","#71ab91", "#71ab51","#43b9c9", "#cf6390","#abb37b","#def485",
           "#583d3f", "#89a61e","#0ff23d", "#9a3ef6", "#1284d8","#8984d8", "#7d056a", "#5f0cab", "#935821" )
names(myPal) <- levels(colData(sfe)[[clustName]])

# myPalette <- pals::polychrome(colourCount)
# names(myPalette) <- levels(colData(sfe)[[clustName]])

plist <- list()
for (i in 1:length(unique(sfe$region_id))){
    region <- unique(sfe$region_id)[[i]]
    sub.sfe <- sfe[,sfe$region_id==region]
    
    p <- make_escheR(sub.sfe, y_reverse=FALSE) %>%
        add_fill(var=clustName)+
        #scale_fill_manual(values=getPalette(colourCount))
        scale_fill_manual(values=myPal)
    plist[[i]] <- p
}

pdfname <- paste0(sprintf("celltype-banksy-bothslides-res%s-lambda%s", res, lambda), ".pdf")

pdf(here("plots", "cindy", "06_cellTyping", "banksy", pdfname), height=15, width=20)
do.call(gridExtra::grid.arrange, c(plist, ncol=3))
for (i in 1:length(plist)){
    print(plist[[i]])
}
dev.off()

