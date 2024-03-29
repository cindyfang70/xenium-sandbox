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
library(rjson)
library(Matrix)
library(DropletUtils)
library(vroom)
library(arrow)
library(sf)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(BiocParallel)

# Reference: https://github.com/pachterlab/SFEData/blob/main/inst/scripts/make-data.R

#fromJSON(file="data/Human_Brain/Xenium_V1_FFPE_Human_Brain_Healthy_With_Addon_gene_panel.json")

# Read in the data
readXenium <- function(counts_path, cell_info_path, cell_poly_path, nuc_poly_path){
    sce <- read10xCounts(counts_path)
    counts(sce) <- as(realize(counts(sce)), "dgCMatrix")
    
    cell_info <- vroom(cell_info_path)
    
    cell_schema <- schema(cell_id=string(),
                          vertex_x=float64(),
                          vertex_y=float64())
    
    cell_poly <- open_dataset(cell_poly_path,
                              schema=cell_schema) %>%
        collect()
    nuc_poly <- open_dataset(nuc_poly_path,
                             schema=cell_schema) %>%
        collect()
    
    names(cell_poly)[1] <- "ID"
    names(nuc_poly)[1] <- "ID"
    
    cells_sf <- df2sf(cell_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
    nuc_sf <- df2sf(nuc_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
    
    all(st_is_valid(cells_sf))
    all(st_is_valid(nuc_sf))
    
    ind_invalid <- !st_is_valid(nuc_sf)
    nuc_sf[ind_invalid,] <- nngeo::st_remove_holes(st_buffer(nuc_sf[ind_invalid,], 0))
    
    colData(sce) <- cbind(colData(sce), cell_info[,-1])
    spe <- toSpatialExperiment(sce, spatialCoordsNames = c("x_centroid", "y_centroid"))
    sfe <- toSpatialFeatureExperiment(spe)
    
    cellSeg(sfe, withDimnames = FALSE) <- cells_sf
    nucSeg(sfe, withDimnames = FALSE) <- nuc_sf
    
    
    # Add some QC Metrics
    colData(sfe)$nCounts <- colSums(counts(sfe))
    colData(sfe)$nGenes <- colSums(counts(sfe) > 0)
    
    rowData(sfe)$means <- rowMeans(counts(sfe))
    rowData(sfe)$vars <- rowVars(counts(sfe))
    rowData(sfe)$cv2 <- rowData(sfe)$vars/rowData(sfe)$means^2
    
    colnames(sfe) <- seq_len(ncol(sfe))
    
    # map ensembl ids to gene symbols
    # rowData(sfe)$SYMBOL <- mapIds(org.Hs.eg.db, keys = rownames(sfe), keytype = "ENSEMBL", column="SYMBOL")
    # rowData(sfe)$ENSEMBL <- rownames(sfe)
    # rownames(sfe) <- uniquifyFeatureNames(ID=rownames(sfe),
    #                             names=rowData(sfe)$SYMBOL)
    # 
    return(sfe)
}


sfe.healthy <- readXenium("data/Human_Brain/cell_feature_matrix.h5", "data/Human_Brain/cells.csv.gz",
                          "data/Human_Brain/cell_boundaries.parquet",
                          "data/Human_Brain/nucleus_boundaries.parquet")

rownames(sfe.healthy) <- uniquifyFeatureNames(ID=rownames(sfe),  names=rowData(sfe)$Symbol)

plotSpatialFeature(sfe.healthy, "nCounts", colGeometryName = "cellSeg")

n_panel <- 319

colData(sfe.healthy)$nCounts_normed <- sfe.healthy$nCounts/n_panel
colData(sfe.healthy)$nGenes_normed <- sfe.healthy$nGenes/n_panel

plotColDataHistogram(sfe.healthy, c("nCounts_normed", "nGenes_normed"))
plotSpatialFeature(sfe.healthy, "nGenes", colGeometryName = "cellSeg")
plotColData(sfe.healthy, x="nCounts", y="nGenes", bins = 150)
plotColDataHistogram(sfe.healthy, c("cell_area", "nucleus_area"), scales = "free_y")

plotSpatialFeature(sfe.healthy, "cell_area", colGeometryName = "cellSeg")
plotSpatialFeature(sfe.healthy, "nucleus_area", colGeometryName = "nucSeg")

plotColData(sfe.healthy, x="cell_area", y="nucleus_area", bins = 100)
colData(sfe.healthy)$prop_nuc <- sfe.healthy$nucleus_area / sfe.healthy$cell_area

plotColDataHistogram(sfe.healthy, "prop_nuc")
plotSpatialFeature(sfe.healthy, "prop_nuc", colGeometryName = "cellSeg")
plotColData(sfe.healthy, x="cell_area", y="prop_nuc")

plotColData(sfe.healthy, x="nucleus_area", y="prop_nuc", bins = 100)

is_blank <- str_detect(rownames(sfe.healthy), "^BLANK_")
sum(is_blank)

is_neg <- str_detect(rownames(sfe.healthy), "^NegControlProbe")
sum(is_neg)

is_neg2 <- str_detect(rownames(sfe.healthy), "^NegControlCodeword")
sum(is_neg2)

is_anti <- str_detect(rownames(sfe.healthy), "^antisense")
sum(is_anti)

is_any_neg <- is_blank | is_neg | is_neg2 | is_anti

sfe.healthy <- addPerCellQCMetrics(sfe.healthy, subsets = list(blank = is_blank,
                                               negProbe = is_neg,
                                               negCodeword = is_neg2,
                                               anti = is_anti,
                                               any_neg = is_any_neg))

cols_use <- names(colData(sfe.healthy))[str_detect(names(colData(sfe.healthy)), "_percent$")]
cols_use <- cols_use[-4]
plotColDataHistogram(sfe.healthy, cols_use, bins = 100, ncol = 3)+
    scale_x_log10()+
    annotation_logticks(sides = "b")

cols_use2 <- names(colData(sfe.healthy))[str_detect(names(colData(sfe.healthy)), "_detected$")]
cols_use2 <- cols_use2[-4]
plotColDataHistogram(sfe.healthy, cols_use2, bins = 20, ncol = 3) +
    # Avoid decimal breaks on x axis unless there're too few breaks
    scale_x_continuous(breaks = scales::breaks_extended(Q = c(1,2,5)))


get_neg_ctrl_outliers <- function(col, sfe) {
    inds <- colData(sfe)$nCounts > 0 & colData(sfe)[[col]] > 0
    df <- colData(sfe)[inds,]
    print(rownames(df))
    outlier_inds <- isOutlier(df[[col]], type = "higher")
    #print(outlier_inds)
    #print(outlier_inds)
    outliers <- rownames(df)[outlier_inds]
    #print(outliers)
    col2 <- str_remove(col, "^subsets_")
    col2 <- str_remove(col2, "_percent$")
    new_colname <- paste("is", col2, "outlier", sep = "_")
    colData(sfe)[[new_colname]] <- colnames(sfe) %in% outliers
    print(sfe)
    sfe
}


cols_use <- names(colData(sfe.healthy))[str_detect(names(colData(sfe.healthy)), "_percent$")]
cols_use <- cols_use[-4]
for (n in cols_use) {
    sfe.healthy <- get_neg_ctrl_outliers(n, sfe.healthy)
}

plotSpatialFeature(sfe.healthy, "is_blank_outlier", colGeometryName = "cellSeg")

plotColData(sfe.healthy, y = "is_blank_outlier", x = "cell_area", 
            point_fun = function(...) list()) 

plotSpatialFeature(sfe.healthy, "is_negProbe_outlier", colGeometryName = "cellSeg")

plotColData(sfe.healthy, y = "is_negProbe_outlier", x = "cell_area", 
            point_fun = function(...) list()) 

plotSpatialFeature(sfe.healthy, "is_negCodeword_outlier", colGeometryName = "cellSeg")

plotColData(sfe.healthy, y = "is_negCodeword_outlier", x = "cell_area", 
            point_fun = function(...) list()) 

plotSpatialFeature(sfe.healthy, "is_any_neg_outlier", colGeometryName = "cellSeg")

plotColData(sfe.healthy, y = "is_any_neg_outlier", x = "cell_area", 
            point_fun = function(...) list()) 

inds_keep <- sfe.healthy$nCounts > 0 & sfe.healthy$nucleus_area < 200 &
    !sfe.healthy$is_blank_outlier & !sfe.healthy$is_negCodeword_outlier & !sfe.healthy$is_negProbe_outlier
sfe.healthy <- sfe.healthy[,inds_keep] # only 74 cells were discarded

plotColDataHistogram(sfe.healthy, cols_use2, bins = 20, ncol = 3) +
    # Avoid decimal breaks on x axis unless there're too few breaks
    scale_x_continuous(breaks = scales::breaks_extended(3, Q = c(1,2,5)))

plotColDataHistogram(sfe.healthy, cols_use, bins = 20, ncol = 3) +
    # Avoid decimal breaks on x axis unless there're too few breaks
    scale_x_continuous(breaks = scales::breaks_extended(3, Q = c(1,2,5)))

rowData(sfe.healthy)$means <- rowMeans(counts(sfe.healthy))
rowData(sfe.healthy)$vars <- rowVars(counts(sfe.healthy))

rowData(sfe.healthy)$is_neg <- is_any_neg
plotRowData(sfe.healthy, x = "means", y = "is_neg") +
    scale_y_log10() +
    annotation_logticks(sides = "b")

plotRowData(sfe.healthy, x="means", y="vars", color_by = "is_neg") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    scale_x_log10() + scale_y_log10() +
    annotation_logticks() +
    coord_equal() +
    labs(color = "Negative control")

sfe.healthy <- logNormCounts(sfe.healthy)
plotSpatialFeature(sfe.healthy, "SST", colGeometryName = "cellSeg")
plotSpatialFeature(sfe.healthy, "MOBP", colGeometryName = "cellSeg")

# Spatial autocorrelation
system.time(
    colGraph(sfe.healthy, "knn5") <- findSpatialNeighbors(sfe.healthy, method = "knearneigh", 
                                                  dist_type = "idw", k = 5, 
                                                  style = "W")
)

sfe.healthy <- colDataMoransI(sfe.healthy, c("nCounts", "nGenes", "cell_area", "nucleus_area"),
                      colGraphName = "knn5")
colFeatureData(sfe.healthy)[c("nCounts", "nGenes", "cell_area", "nucleus_area"),]

sfe.healthy <- colDataUnivariate(sfe.healthy, type = "localmoran", 
                         features = c("nCounts", "nGenes", "cell_area", 
                                      "nucleus_area"),
                         colGraphName = "knn5", BPPARAM = MulticoreParam(2))

plotLocalResult(sfe.healthy, "localmoran",
                features = c("nCounts", "nGenes", "cell_area", "nucleus_area"),
                colGeometryName = "centroids", scattermore = TRUE,
                divergent = TRUE, diverge_center = 0, pointsize = 1)

sfe.healthy <- colDataUnivariate(sfe.healthy, "moran.plot", "nCounts", colGraphName = "knn5")

p1 <- moranPlot(sfe.healthy, "nCounts", binned = TRUE, plot_influential = FALSE) 
p2 <- moranPlot(sfe.healthy, "nCounts", binned = TRUE)
p1 / p2 + plot_layout(guides = "collect")

# Moran's I 

system.time(
    sfe.healthy <- runMoransI(sfe.healthy, colGraphName = "knn5", BPPARAM = MulticoreParam(2))
)
rowData(sfe.healthy)$is_neg <- is_any_neg
plotRowData(sfe.healthy, x = "moran_sample01", y = "is_neg")

# what's the negative control probe with large positive moran's I?
max_ind <- which.max(rowData(sfe.healthy)$moran_sample01[is_any_neg])
top_neg <- rownames(sfe.healthy)[is_any_neg][max_ind]

plotSpatialFeature(sfe.healthy, top_neg, colGeometryName = "centroids",
                   scattermore = TRUE, pointsize = 1)

# Genes with top Moran's I
top_moran <- rownames(sfe.healthy)[order(rowData(sfe.healthy)$moran_sample01, decreasing = TRUE)[1:6]]
plotSpatialFeature(sfe.healthy, top_moran, colGeometryName = "centroids",
                   scattermore = TRUE, ncol = 2, pointsize = 0.5)

plotRowData(sfe.healthy, x = "means", y = "moran_sample01")

# Let's do a comparison of regular and spatial PCA 
set.seed(0)
sfe.healthy <- runPCA(sfe.healthy, ncomponents = 30, scale = TRUE, BSPARAM = IrlbaParam())
ElbowPlot(sfe.healthy, ndims = 30)
plotDimLoadings(sfe.healthy, dims = 1:6)

spatialReducedDim(sfe.healthy, "PCA", 6, colGeometryName = "centroids", divergent = TRUE,
                  diverge_center = 0, ncol = 2, scattermore = TRUE, pointsize = 0.5)

colData(sfe.healthy)$clust.nonspat <- clusterRows(reducedDim(sfe.healthy, "PCA")[,1:15],
                                    BLUSPARAM = SNNGraphParam(
                                        cluster.fun = "leiden",
                                        cluster.args = list(
                                            resolution_parameter = 0.5,
                                            objective_function = "modularity")))

plotPCA(sfe.healthy, ncomponents = 4, colour_by = "clust.nonspat", rasterise = FALSE)
plotSpatialFeature(sfe.healthy, "clust.nonspat", colGeometryName = "cellSeg")

markers <- scran::findMarkers(sfe.healthy, groups = colData(sfe.healthy)$clust.nonspat,
                       test.type = "wilcox", pval.type = "all", direction = "up")

genes_use <- vapply(markers, function(x) rownames(x)[1], FUN.VALUE = character(1))
plotExpression(sfe.healthy, genes_use, x = "clust.nonspat", point_fun = function(...) list())

genes_use2 <- unique(unlist(lapply(markers, function(x) rownames(x)[1:5])))
plotGroupedHeatmap(sfe.healthy, genes_use2, group = "clust.nonspat", colour = scales::viridis_pal()(100))

plotSpatialFeature(sfe.healthy, genes_use, colGeometryName = "centroids", ncol = 3,
                   pointsize = 0.3, scattermore = TRUE)

# spatial PCA
system.time({
    sfe.healthy <- runMultivariate(sfe.healthy, "multispati", colGraphName = "knn5", nfposi = 20,
                           nfnega = 20)
})
ElbowPlot(sfe.healthy, nfnega = 20, reduction = "multispati")







inds <- order(rowSums(logcounts(sfe.healthy)), decreasing = TRUE)[1:50]
mat <- logcounts(sfe.healthy)[inds,]
nf <- 20
g <- colGraph(sfe.healthy, "knn5")

listw<- g
x <- t(mat)
W <- listw2sparse(listw)

x <- sweep(x, 2, colMeans(x))
if (is(x, "dgeMatrix")) x <- as.matrix(x)
covar <- t(x) %*% (W + t(W)) %*% x / (2*nrow(x))

res <- eigs_sym(covar, k = 2*nf, which = "BE")
