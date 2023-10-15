suppressPackageStartupMessages({
    #library(Voyager)
    #library(SFEData)
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
    library(rjson)
    library(Matrix)
    library(DropletUtils)
    library(vroom)
    #library(arrow)
    library(sf)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(BayesSpace)
})

args <- commandArgs(trailingOnly = TRUE)

#data_dir <- "output-XETG00089__0005434__Region_1__20230831__172339"
#sfe <- readRDS("data/output-XETG00089__0005434__Region_1__20230831__172339/BR6471_Post-SFE.RDS")

sfe <- readRDS(args[[1]])

colData(sfe)$nCounts <- colSums(counts(sfe))
colData(sfe)$nGenes <- colSums(counts(sfe) > 0)

# Add some QC metrics for the negative control probes
is_blank <- str_detect(rownames(sfe), "^BLANK_")
is_neg <- str_detect(rownames(sfe), "^NegControlProbe")
is_neg2 <- str_detect(rownames(sfe), "^NegControlCodeword")
is_anti <- str_detect(rownames(sfe), "^antisense")
is_depr <- str_detect(rownames(sfe), "^DeprecatedCodeword")

is_any_neg <- is_blank | is_neg | is_neg2 | is_anti | is_depr
rowData(sfe)$is_neg <- is_any_neg

n_panel <- nrow(sfe) - sum(is_any_neg)

sfe <- addPerCellQCMetrics(sfe, subsets = list(blank = is_blank,
                                               negProbe = is_neg,
                                               negCodeword = is_neg2,
                                               anti = is_anti,
                                               depr = is_depr,
                                               any_neg = is_any_neg))

#plotColData(sfe, x="nucleus_area", y="nCounts", bin=100)
#plotColData(sfe, x="nucleus_area", y="subsets_any_neg_percent", bin=100)
#hist(sfe$nucleus_area)

# Do QC based on negative control probes
cols_use <- names(colData(sfe))[str_detect(names(colData(sfe)), "_percent$")]
cols_use <- cols_use[!grepl("anti", cols_use)] # no antisense controls detected, so avoid dividing by zero

#plotColDataHistogram(sfe, cols_use, bins = 100, ncol = 3)+
    #scale_x_log10()+
    #annotation_logticks(sides = "b")

cols_use2 <- names(colData(sfe))[str_detect(names(colData(sfe)), "_detected$")]
cols_use2 <- cols_use2[!grepl("anti", cols_use)] 

get_neg_ctrl_outliers <- function(col, sfe) {
    inds <- colData(sfe)$nCounts > 0 & colData(sfe)[[col]] > 0
    df <- colData(sfe)[inds,]
    #print(rownames(df))
    outlier_inds <- isOutlier(df[[col]], type = "higher")
    #print(outlier_inds)
    #print(outlier_inds)
    outliers <- rownames(df)[outlier_inds]
    #print(outliers)
    col2 <- str_remove(col, "^subsets_")
    col2 <- str_remove(col2, "_percent$")
    new_colname <- paste("is", col2, "outlier", sep = "_")
    colData(sfe)[[new_colname]] <- colnames(sfe) %in% outliers
    sfe
}

for (n in cols_use) {
    sfe <- get_neg_ctrl_outliers(n, sfe)
}

#plotColData(sfe, x="nucleus_area", y="subsets_anti_percent", bin=100)

inds_keep <- sfe$nCounts > 0  & !sfe$is_blank_outlier & !sfe$is_negCodeword_outlier & 
    !sfe$is_negProbe_outlier & !sfe$is_depr_outlier
sum(!inds_keep) # [1] 3072
sfe <- sfe[,inds_keep] # 3072 cells discarded

# Compute some gene-level summary stats
rowData(sfe)$means <- rowMeans(counts(sfe))
rowData(sfe)$vars <- rowVars(counts(sfe))
rowData(sfe)$cv2 <- rowData(sfe)$vars/rowData(sfe)$means^2

# Plot the mean-var relationship
#plotMeanVar(mean_emp=rowData(sfe)$means,
            #var_emp=rowData(sfe)$vars,
            #plotTitle="Xenium DLPFC")

# Try GLM-PCA approach, compute null residuals and run PCA on them
sfe <- scry::nullResiduals(sfe, assay="counts", fam="poisson", type="pearson")
sfe <- scater::runPCA(sfe, ncomponents=50, 
                      ntop = 1000,
                      exprs_values = "poisson_pearson_residuals",
                      scale = TRUE, name = "GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

# look at where the scree plot dips to determine # PCs
#ElbowPlot(sfe, reduction="GLM-PCA", ndims=50)

# Normalize and run PCA
sfe <- scuttle::logNormCounts(sfe)
sfe <- scater::runPCA(sfe, ncomponents=30, scale=TRUE,
                      exprs_values="logcounts")
#ElbowPlot(sfe, reduction="PCA", ndims=30)


# Convert the sfe to spe
sce <- SingleCellExperiment(counts(sfe))
colData(sce) <- colData(sfe)

# Get the x and y centroids for this region, and add them to the coldata
sce$x_centroid <- spatialCoords(sfe)[,1]
sce$y_centroid <- spatialCoords(sfe)[,2]

reducedDim(sce, "PCA") <- reducedDim(sfe, "PCA")

# Convert to SPE
spe <- toSpatialExperiment(sce, spatialCoordsNames=c("x_centroid", "y_centroid"))

spe$col <- spatialCoords(spe)[,2]
spe$row <- spatialCoords(spe)[,1]
# try clustering with BayesSpace
library(BayesSpace)
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

set.seed(149)
Sys.time()
spe <- spatialCluster(spe, use.dimred = "PCA", q = 10 ,nrep=5000)
Sys.time()

saveRDS(spe, "bayesSpace-BR_6471_post.RDS")
