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
library(AnnotationDbi)
library(org.Hs.eg.db)
library(BiocParallel)
library(dplyr)

# Reference: https://github.com/pachterlab/SFEData/blob/main/inst/scripts/make-data.R

# Read in Xenium data and create an SFE object
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
    
    is_blank <- str_detect(rownames(sfe), "^BLANK_")
    is_neg <- str_detect(rownames(sfe), "^NegControlProbe")
    is_neg2 <- str_detect(rownames(sfe), "^NegControlCodeword")
    is_anti <- str_detect(rownames(sfe), "^antisense")
    
    
    is_any_neg <- is_blank | is_neg | is_neg2 | is_anti
    rowData(sfe)$is_neg <- is_any_neg
    
    n_panel <- nrow(sfe) - sum(is_any_neg)
    #print(n_panel)
    
    colData(sfe)$nCounts_normed <- sfe$nCounts/n_panel
    colData(sfe)$nGenes_normed <- sfe$nGenes/n_panel
    colData(sfe)$prop_nuc <- sfe$nucleus_area / sfe$cell_area
    
    
    
    sfe <- addPerCellQCMetrics(sfe, subsets = list(blank = is_blank,
                                                                   negProbe = is_neg,
                                                                   negCodeword = is_neg2,
                                                                   anti = is_anti,
                                                                   any_neg = is_any_neg))
    
    rowData(sfe)$means <- rowMeans(counts(sfe))
    rowData(sfe)$vars <- rowVars(counts(sfe))
    rowData(sfe)$cv2 <- rowData(sfe)$vars/rowData(sfe)$means^2
    
    # Add cell ids and make gene names unique
    colnames(sfe) <- seq_len(ncol(sfe))
    rownames(sfe) <- uniquifyFeatureNames(ID=rownames(sfe),  names=rowData(sfe)$Symbol)
    
    return(sfe)
}

# Determine which cells are outliers based on percentage of counts from a negative control probe
get_neg_ctrl_outliers <- function(col, sfe) {
    # Only consider the cells which have non-zero percentage of counts from neg control
    inds <- colData(sfe)$nCounts > 0 & colData(sfe)[[col]] > 0
    df <- colData(sfe)[inds,]
    print(rownames(df))
    outlier_inds <- isOutlier(df[[col]], type = "higher")
    outliers <- rownames(df)[outlier_inds]
    col2 <- str_remove(col, "^subsets_")
    col2 <- str_remove(col2, "_percent$")
    new_colname <- paste("is", col2, "outlier", sep = "_")
    colData(sfe)[[new_colname]] <- colnames(sfe) %in% outliers
    print(sfe)
    sfe
}
