library(Voyager)
#library(SFEData)
#library(patchwork)
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
library(here)

# Reference: https://github.com/pachterlab/SFEData/blob/main/inst/scripts/make-data.R

# Read in Xenium data and create an SFE object
readXenium <- function(dir_name){
    # Find the files
    counts_path <- here("data", dir_name, "cell_feature_matrix.h5")
    cell_info_path <- here("data", dir_name, "cells.csv.gz")
    cell_poly_path <- here("data", dir_name, "cell_boundaries.parquet")
    nuc_poly_path <- here("data", dir_name, "nucleus_boundaries.parquet")
    
    # Read in the data
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
    
    colData(sce) <- cbind(colData(sce), cell_info)
    print(cell_info)
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
    is_depr <- str_detect(rownames(sfe), "^DeprecatedCodeword")
    
    is_any_neg <- is_blank | is_neg | is_neg2 | is_anti | is_depr
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
    depr = is_depr,
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


plotMeanVar <- function(mean_emp, var_emp, plotTitle){
    model = lm(var_emp ~ 1*mean_emp + I(mean_emp^2) + 0, tibble(mean_emp, var_emp))
    phi = 1/coef(model)["I(mean_emp^2)"]
    
    
    mean_var_tb <- tibble(mean_emp = mean_emp,
                          var_emp = var_emp,
                          nbinomial = mean_emp + mean_emp^2 * 1/phi)
    print(mean_var_tb)
    p <- mean_var_tb %>%
        ggplot(aes(x = mean_emp, y = var_emp)) + 
        geom_point(alpha = 0.3) + 
        geom_line(data = mean_var_tb %>% dplyr::select(mean_emp, var_emp, nbinomial) %>% 
                      tidyr::pivot_longer(cols=-mean_emp, names_to = "model", values_to = "var_value")%>%
                      filter(model %in% c("nbinomial")),
                  aes(x = mean_emp, y = var_value, colour=model)) +
        scale_x_log10() + scale_y_log10() +
        labs(x = "Log of mean expression",
             y = "Log of variance") +
        geom_abline(slope = 1, intercept = 0, color = "black")+
        ggtitle(plotTitle)+
        theme_bw()
    
    return(p)
}

getRegionInds <- function(region, data_dir, sfe){
    # Read in the cell ids of a particular region
    region_ids_fname <- paste0(region, "_cell_ids.csv")
    region_ids <- readr::read_csv2(here("data", data_dir, region_ids_fname))
    
    # Extract just the cell ids
    region_ids <- region_ids %>%
        tidyr::separate(everything(), sep="\\,", into=c("cell_id", NA, NA, NA))
    region_ids <- region_ids[-c(1,2),]
    
    # Get the sfe indices for the cells in the region
    region_inds <- which(sfe$cell_id %in% region_ids$cell_id)
    
    return(region_inds)
    
}

SFEtoSPE <- function(sfe){
    # First convert to an SCE
    sce <- SingleCellExperiment(assays=list(counts=counts(sfe)))
    colData(sce) <- colData(sfe)
    
    # Get the x and y centroids for this region, and add them to the coldata
    sce$x_centroid <- spatialCoords(sfe)[,1]
    sce$y_centroid <- spatialCoords(sfe)[,2]
    
    # Convert to SPE
    spe <- toSpatialExperiment(sce, spatialCoordsNames=c("x_centroid", "y_centroid"))
    
    return(spe)
}

get_neg_ctrl_outliers <- function(col, spe) {
    inds <- colData(spe)$nCounts > 0 & colData(spe)[[col]] > 0
    df <- colData(spe)[inds,]
    outlier_inds <- isOutlier(df[[col]], type = "higher")
    outliers <- rownames(df)[outlier_inds]
    col2 <- str_remove(col, "^subsets_")
    col2 <- str_remove(col2, "_percent$")
    new_colname <- paste("is", col2, "outlier", sep = "_")
    colData(spe)[[new_colname]] <- colnames(spe) %in% outliers
    return(spe)
}

filterCells <- function(sfe){
    cols_use <- names(colData(sfe))[str_detect(
        names(colData(sfe)), "_percent$")]
    cols_use <- cols_use[!grepl("anti", cols_use)]
    
    for (n in cols_use) {
        sfe <- get_neg_ctrl_outliers(n, sfe)
    }
    inds_keep <- sfe$nucleus_area <=200 & !sfe$is_blank_outlier &
        !sfe$is_depr_outlier &!sfe$is_negProbe_outlier &
        !sfe$is_negCodeword_outlier & sfe$nCounts > 0
    
    sfe <- sfe[,inds_keep]
    
    # remove the negative controls
    sfe <- sfe[which(rowData(sfe)$Type=="Gene Expression"), ]
    return(sfe)
}
