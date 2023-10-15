suppressPackageStartupMessages({
    library(Voyager)
    library(SFEData)
    library(patchwork)
    library(SpatialFeatureExperiment)
    library(SingleCellExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(stringr)
    library(scuttle)
    library(BiocSingular)
    library(scater)
    library(Matrix)
    library(DropletUtils)
    library(sf)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(BiocSingular)
    library(gridExtra)
})

######################################################
# Identify lowly expressed genes and make some plots #
######################################################

##### Slide 5434 #####
data_dir <- "output-XETG00089__0005434__Region_1__20230831__172339"

# Read in the three regions
br6471_p <- readRDS(here("data", data_dir, "Br6471_Post_SFE.RDS"))
br6522_p <- readRDS(here("data", data_dir, "Br6522_Post_SFE.RDS"))
br8667_m <- readRDS(here("data", data_dir, "Br8667_Mid_SFE.RDS"))

findLowlyExpressedGenes <- function(sfe){
    # Find the number of times each gene has non-zero expression in a cell
    sfe <- sfe[,sfe$nCounts > 0]
    sfe <- logNormCounts(sfe)
    geneCounts <- rowSums(counts(sfe) > 0)
    
    # Remove the negative controls and sort from low to high
    geneCounts <- sort(geneCounts[!grepl("Codeword|Neg|BLANK", names(geneCounts))])
    
    fiveLowest <- head(geneCounts, 5)
    print(fiveLowest)
    
    p <- plotSpatialFeature(sfe, features=names(fiveLowest), 
                            colGeometryName="cellSeg", exprs_values="logcounts")
    return(p)
}

pdf(here("plots", "xenium_plots", "lowlyExpressedGenes_5434.pdf"))
findLowlyExpressedGenes(br6471_p)
findLowlyExpressedGenes(br6522_p)
findLowlyExpressedGenes(br8667_m)
dev.off()

##### Slide 5548 #####
data_dir <- "output-XETG00089__0005548__Region_1__20230831__172339"

# Read in the three regions
br6471_p <- readRDS(here("data", data_dir, "Br6471_Post_SFE.RDS"))
br2743_m <- readRDS(here("data", data_dir, "Br2743_Mid_SFE.RDS"))
br8667_m <- readRDS(here("data", data_dir, "Br8667_Mid_SFE.RDS"))

pdf(here("plots", "xenium_plots", "lowlyExpressedGenes_5548.pdf"))
findLowlyExpressedGenes(br6471_p)
findLowlyExpressedGenes(br2743_m)
findLowlyExpressedGenes(br8667_m)
dev.off()
