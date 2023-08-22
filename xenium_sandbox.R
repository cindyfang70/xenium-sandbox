library(Voyager)
library(SFEData)
library(patchwork)
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(ggplot2)
library(bluster)
library(stringr)
library(scuttle)
library(BiocSingular)
library(scater)

sfe <- JanesickBreastData(dataset = "rep2")
colnames(sfe) <- seq_len(ncol(sfe))
plotGeometry(sfe, "cellSeg")

n_panel <- 313 # this is part of the experimental design
colData(sfe)$nCounts_normed <- sfe$nCounts/n_panel
colData(sfe)$nGenes_normed <- sfe$nGenes/n_panel
colData(sfe)$prop_nuc <- sfe$nucleus_area / sfe$cell_area

plotColDataHistogram(sfe, c("nCounts_normed", "nGenes_normed"))
plotColDataBin2D(sfe, "cell_area", "prop_nuc") + 
    scale_fill_viridis_c()

plotColDataBin2D(sfe, "nucleus_area", "prop_nuc") + 
    scale_fill_viridis_c()

is_blank <- str_detect(rownames(sfe), "^BLANK_")
sum(is_blank)

is_neg <- str_detect(rownames(sfe), "^NegControlProbe")
sum(is_neg)

is_neg2 <- str_detect(rownames(sfe), "^NegControlCodeword")
sum(is_neg2)

is_anti <- str_detect(rownames(sfe), "^antisense")
sum(is_anti)

is_any_neg <- is_blank | is_neg | is_neg2 | is_anti

sfe <- addPerCellQCMetrics(sfe, subsets = list(blank = is_blank,
                                               negProbe = is_neg,
                                               negCodeword = is_neg2,
                                               anti = is_anti,
                                               any_neg = is_any_neg))


cols_use <- names(colData(sfe))[str_detect(names(colData(sfe)), "_percent$")]
plotColDataHistogram(sfe, cols_use, bins = 100, ncol = 3)+
    scale_x_log10()

get_neg_ctrl_outliers <- function(col, sfe) {
    inds <- colData(sfe)$nCounts > 0 & colData(sfe)[[col]] > 0
    df <- colData(sfe)[inds,]
    outlier_inds <- isOutlier(df[[col]], type = "higher")
    outliers <- rownames(df)[outlier_inds]
    col2 <- str_remove(col, "^subsets_")
    col2 <- str_remove(col2, "_percent$")
    new_colname <- paste("is", col2, "outlier", sep = "_")
    colData(sfe)[[new_colname]] <- colnames(sfe) %in% outliers
    sfe
}

cols_use <- names(colData(sfe))[str_detect(names(colData(sfe)), "_percent$")]
for (n in cols_use) {
    sfe <- get_neg_ctrl_outliers(n, sfe)
}

inds_keep <- sfe$nCounts > 0 & sfe$nucleus_area < 400 & !sfe$is_anti_outlier &
    !sfe$is_blank_outlier & !sfe$is_negCodeword_outlier & !sfe$is_negProbe_outlier
(sfe <- sfe[,inds_keep])

rowData(sfe)$means <- rowMeans(counts(sfe))
rowData(sfe)$vars <- rowVars(counts(sfe))
rowData(sfe)$is_neg <- is_any_neg

sfe <- logNormCounts(sfe)


sfe$frac_non_zero <- nexprs(sfe, byrow=FALSE)/colSums(counts(sfe))
sfe$frac_zero <- 1-sfe$frac_non_zero


set.seed(0)
# try nonspatial PCA
system.time(
    sfe <- runPCA(sfe, ncomponents = 20, subset_row = !is_blank,
                  exprs_values = "logcounts",
                  scale = TRUE, BSPARAM = IrlbaParam())
)

sfe <- runTSNE(sfe, dimred="PCA")
plotTSNE(sfe, colour_by="frac_zero")+
    aes(alpha=0.2)
plotSpatialFeature(sfe, "frac_zero")+aes(alpha=0.5)

# try MULTISPATI PCA
system.time(
    colGraph(sfe, "knn5") <- findSpatialNeighbors(sfe, method = "knearneigh", 
                                                  dist_type = "idw", k = 5, 
                                                  style = "W")
)
system.time({
    sfe <- runMultivariate(sfe, "multispati", colGraphName = "knn5", nfposi = 20,
                           nfnega = 20)
})

sfe <- runTSNE(sfe, dimred="multispati")
plotTSNE(sfe, colour_by="frac_zero")+
    aes(alpha=0.2)
