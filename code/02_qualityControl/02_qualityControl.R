suppressPackageStartupMessages({
    library(Voyager)
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
    library(sf)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(schex)
})


################################################
# Filter out low quality cells for each region #
################################################
source(here("code", "01_createSCE", "xenium_helpers.R"))

##### Slide 5434 #####
data_dir <- "slide-5434"

# Read in the three regions
br6471_p <- readRDS(here("processed-data", "cindy", data_dir, "Br6471_Post_SFE.RDS"))
br6522_p <- readRDS(here("processed-data", "cindy", data_dir, "Br6522_Post_SFE.RDS"))
br8667_m <- readRDS(here("processed-data", "cindy", data_dir, "Br8667_Mid_SFE.RDS"))


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
    
    sfe$keep <- inds_keep
    sfe <- sfe[,inds_keep]
    
    # remove the negative controls
    sfe <- sfe[which(rowData(sfe)$Type=="Gene Expression"), ]
    return(sfe)
}
br6471_p_filt <- filterCells(br6471_p)
br6522_p_filt <- filterCells(br6522_p)
br8667_m_filt <- filterCells(br8667_m)

cols_use <- names(colData(br8667_m))[str_detect(names(colData(sfe)), "_percent$")]
cols_use <- cols_use[!grepl("subsets_anti_percent", cols_use)]

pdf(here("plots", "02_qualityControl", "post_filtering_histograms-5434.pdf"), height=5, width=12)
plotColDataHistogram(br6471_p_filt, cols_use, bins = 100, ncol = 3)+ 
    scale_x_log10() +
    annotation_logticks(sides = "b")+
    ggtitle(unique(br6471_p$region_id))+
    xlab("percentage of total counts")
plotColDataHistogram(br6522_p_filt, cols_use, bins = 100, ncol = 3)+ 
    scale_x_log10() +
    annotation_logticks(sides = "b")+
    ggtitle(unique(br6522_p$region_id))+
    xlab("percentage of total counts")
plotColDataHistogram(br8667_m_filt, cols_use, bins = 100, ncol = 3)+ 
    scale_x_log10() +
    annotation_logticks(sides = "b")+
    ggtitle(unique(br8667_m$region_id))+
    xlab("percentage of total counts")
dev.off()

# p1 <- plotColData(br6471_p_filt, y="subsets_any_neg_sum",
#             colour_by="keep",
#             show_median=TRUE,
#             point_alpha=0.3)+
#     ggtitle(unique(br6471_p_filt$region_id)) 
# 
# p2 <- plotColData(br6522_p_filt, y="subsets_any_neg_sum",
#             colour_by="keep",
#             show_median=TRUE,
#             point_alpha=0.3)+
#     ggtitle(unique(br6522_p_filt$region_id))
# 
# p3 <- plotColData(br8667_m_filt, y="subsets_any_neg_sum",
#             colour_by="keep",
#             show_median=TRUE,
#             point_alpha=0.3)+
#     ggtitle(unique(br8667_m_filt$region_id))
# 
# pdf("qcplots.pdf", height=8, width=20)
# p1+p2+p3
# dev.off()

saveRDS(br6471_p_filt, here("data", data_dir, "Br6471_Post_SFE_filt.RDS"))
saveRDS(br6522_p_filt, here("data", data_dir, "Br6522_Post_SFE_filt.RDS"))
saveRDS(br8667_m_filt , here("data", data_dir, "Br8667_Mid_SFE_filt.RDS"))


##### Slide 5548 ##### 

data_dir <- "slide-5548"

# Read in the three regions
br6471_p <- readRDS(here("data", data_dir, "Br6471_Post_SFE.RDS"))
br2743_m <- readRDS(here("data", data_dir, "Br2743_Mid_SFE.RDS"))
br8667_m <- readRDS(here("data", data_dir, "Br8667_Mid_SFE.RDS"))

br6471_p_filt <- filterCells(br6471_p)
br2743_m_filt <- filterCells(br2743_m)
br8667_m_filt <- filterCells(br8667_m)

saveRDS(br6471_p_filt, here("data", data_dir, "Br6471_Post_SFE_filt.RDS"))
saveRDS(br2743_m_filt, here("data", data_dir, "Br2743_Mid_SFE_filt.RDS"))
saveRDS(br8667_m_filt, here("data", data_dir, "Br8667_Mid_SFE_filt.RDS"))

