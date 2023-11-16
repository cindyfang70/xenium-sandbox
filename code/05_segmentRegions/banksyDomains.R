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
# Trying out Banksy for spatial domain detection:
# tutorial: https://github.com/prabhakarlab/Banksy/tree/bioc
#------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
lambda <- args[[2]]
#k <- args[[3]]
res <- args[[3]]

# read in the data
sfe <- readRDS(args[[1]])
# delete this line after the first run
colData(sfe) <- colData(sfe)[,!grepl("^clust", colnames(colData(sfe)))]


sfe <- computeLibraryFactors(sfe)
aname <- "normcounts"
assay(sfe, aname) <- normalizeCounts(sfe, log = FALSE)


# Compute neighbourhood matrices
#lambda <- c(0, 0.9)
k_geom <- c(15, 30)

sfe <- Banksy::computeBanksy(sfe, assay_name = aname, compute_agf = TRUE,
                             k_geom = k_geom)

# run PCA and UMAP, then Leiden clustering to find spatial domains
# lambda is in [0,1], higher values of lambda puts more weight on spatial
# information
set.seed(1000)
sfe <- Banksy::runBanksyPCA(sfe, use_agf = TRUE, lambda = lambda)
sfe <- Banksy::runBanksyUMAP(sfe, use_agf = TRUE, lambda = lambda)
sfe <- Banksy::clusterBanksy(sfe, use_agf = TRUE, lambda = lambda, 
                             resolution = res)

# save the Banksy results
sfeName <- str_replace(args[[1]], "\\.RDS", "-with-banksy-domains.RDS")
saveRDS(sfe, sfeName)

# default for the leiden algorithm in clusterBanksy is k_neighbours=50
clustName <- sprintf("clust_M1_lam%s_k%s_res%s", lambda, 50, res)


# plot using escheR 
colourCount = nlevels(colData(sfe)[[clustName]])
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

p1 <- make_escheR(sfe) %>%
    add_fill(var=clustName)+
    scale_fill_manual(values=getPalette(colourCount))
    

fname <- paste(sfe$region_id[[1]], "Banksy", "lambda", 
               lambda, "res", res, sep="-")
pdfname <- paste0(fname, ".pdf")
pdf(here("plots", "cindy", "05_segmentRegions", pdfname))
p1
dev.off()

