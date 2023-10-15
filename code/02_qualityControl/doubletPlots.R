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
    library(scDblFinder)
    library(gridExtra)
})

###################################################
# Identify doublets and make some related plots   #
###################################################
findDbls <- function(sfe){
    # Computes doublet densities and doublet calling on an SFE object, and stores
    # the results in colData.
    dbl.dens <- computeDoubletDensity(sfe, dims=25) #set seed?
    sfe$doubletScore <- dbl.dens
    
    dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                     method="griffiths", returnType="call")
    sfe$isDoublet <- dbl.calls
    return(sfe)
}

makeDblPlots <- function(sfe){
    # make doublets-related plots for an SFE object that has doublet 
    # densities and doublet calls stored in colData
    
    # First find the number of doublets in each dataset
    ndbls <- sprintf("Number of doublets: %s", sum(sfe$isDoublet == "doublet"))
    ntotal <- sprintf("Total number of cells %s", dim(sfe)[[2]])
    
    ndbls.lab <- data.frame(x=c(10, 10), y=c(200,150), lab=c(ndbls, ntotal))
    
    p1 <- plotSpatialFeature(sfe, "isDoublet", colGeometryName="cellSeg")+
        ggtitle(unique(sfe$region_id))
    
    p2 <- plotColData(sfe,  x="doubletScore", y="nucleus_area", colour_by="isDoublet")+
        geom_smooth(method="lm")+
        ggtitle(unique(sfe$region_id))+
        geom_text(data=ndbls.lab, aes(x=x, y=y,label=lab), position="identity")
    
    p3 <- plotColData(sfe, x="doubletScore", y="subsets_any_neg_percent", colour_by="isDoublet")+
        geom_smooth(method="lm")+
        ggtitle(unique(sfe$region_id))
    
    p4 <- plotColData(sfe, x="doubletScore", y="unassigned_codeword_counts", colour_by="isDoublet")+
        geom_smooth(method="lm")+
        ggtitle(unique(sfe$region_id))
    
    return(list(p1, p2, p3, p4))
}

##### Slide 5434 #####
data_dir <- "output-XETG00089__0005434__Region_1__20230831__172339"

# Read in the three regions
br6471_p <- readRDS(here("data", data_dir, "Br6471_Post_SFE.RDS"))
br6522_p <- readRDS(here("data", data_dir, "Br6522_Post_SFE.RDS"))
br8667_m <- readRDS(here("data", data_dir, "Br8667_Mid_SFE.RDS"))

# Doublet detection and plotting for each of the three regions
br6471_p <- br6471_p[,br6471_p$nCounts > 0] 
br6471_p <- findDbls(br6471_p)
br6471_p_plots <- makeDblPlots(br6471_p)

br6522_p <- br6522_p[,br6522_p$nCounts > 0]
br6522_p <- findDbls(br6522_p)
br6522_p_plots <- makeDblPlots(br6522_p)

br8667_m <- br8667_m[,br8667_m$nCounts > 0]
br8667_m <- findDbls(br8667_m)
br8667_m_plots <- makeDblPlots(br8667_m)


pdf(here("plots", "xenium_plots", "doublet_plots_5434.pdf"), height=8, width=15)
do.call("grid.arrange", c(br6471_p_plots, ncol=2))
do.call("grid.arrange", c(br6522_p_plots, ncol=2))
do.call("grid.arrange", c(br8667_m_plots, ncol=2))
dev.off()

##### Slide 5548 #####
data_dir <- "output-XETG00089__0005548__Region_1__20230831__172339"

# Read in the three regions
br6471_p <- readRDS(here("data", data_dir, "Br6471_Post_SFE.RDS"))
br2743_m <- readRDS(here("data", data_dir, "Br2743_Mid_SFE.RDS"))
br8667_m <- readRDS(here("data", data_dir, "Br8667_Mid_SFE.RDS"))

br6471_p <- br6471_p[,br6471_p$nCounts > 0] 
br6471_p <- findDbls(br6471_p)
br6471_p_plots <- makeDblPlots(br6471_p)

br2743_m <- br2743_m[,br2743_m$nCounts > 0]
br2743_m <- findDbls(br2743_m)
br2743_m_plots <- makeDblPlots(br2743_m)

br8667_m <- br8667_m[,br8667_m$nCounts > 0]
br8667_m <- findDbls(br8667_m)
br8667_m_plots <- makeDblPlots(br8667_m)

pdf(here("plots", "xenium_plots", "doublet_plots_5548.pdf"), height=8, width=15)
do.call("grid.arrange", c(br6471_p_plots, ncol=2))
do.call("grid.arrange", c(br2743_m_plots, ncol=2))
do.call("grid.arrange", c(br8667_m_plots, ncol=2))
dev.off()

# 
# plotSpatialFeature(br6471_p, c("FASLG", "TOP2A", "NKG7",
#                                "CD3G", "IL7R"), colGeometryName="cellSeg", exprs_values="counts")

