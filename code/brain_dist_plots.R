suppressPackageStartupMessages({
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
    library(here)
})
source(here("code/xenium_helpers.R"))

sfe.healthy <- readXenium("data/Human_Brain/cell_feature_matrix.h5", "data/Human_Brain/cells.csv.gz",
                          "data/Human_Brain/cell_boundaries.parquet",
                          "data/Human_Brain/nucleus_boundaries.parquet")


# cols_use <- names(colData(sfe.healthy))[str_detect(names(colData(sfe.healthy)), "_percent$")]
# cols_use <- cols_use[-4]
# 
# for (n in cols_use) {
#     sfe.healthy <- get_neg_ctrl_outliers(n, sfe.healthy)
# }
# 
# inds_keep <- sfe.healthy$nCounts > 0 & sfe.healthy$nucleus_area < 200 &
#     !sfe.healthy$is_blank_outlier & !sfe.healthy$is_negCodeword_outlier & !sfe.healthy$is_negProbe_outlier
# sfe.healthy <- sfe.healthy[,inds_keep] # only 74 cells were discarded
# 
# rowData(sfe.healthy)$means <- rowMeans(counts(sfe.healthy))
# rowData(sfe.healthy)$vars <- rowVars(counts(sfe.healthy))
# 
# plotRowData(sfe.healthy, x="means", y="vars", color_by = "is_neg") +
#     geom_abline(slope = 1, intercept = 0, color = "red") +
#     scale_x_log10() + scale_y_log10() +
#     annotation_logticks() +
#     coord_equal() +
#     labs(color = "Negative control")

plotMeanVar <- function(mean_emp, var_emp, plotTitle){
    model = lm(var_emp ~ 1*mean_emp + I(mean_emp^2) + 0, tibble(mean_emp, var_emp))
    phi = 1/coef(model)["I(mean_emp^2)"]
    
    
    mean_var_tb <- tibble(mean_emp = mean_emp,
                          var_emp = var_emp,
                          nbinomial = mean_emp + mean_emp^2 * 1/phi)
    print(mean_var_tb)
    p4 <- mean_var_tb %>%
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
    
    return(p4)
}

p1 <- plotMeanVar(rowData(sfe.healthy)$means, rowData(sfe.healthy)$vars, "Xenium Human Brain Healthy")

# Glioblastoma
glio.sfe <- readXenium("data/Human_Brain_Glioblastoma/cell_feature_matrix.h5", 
                       "data/Human_Brain_Glioblastoma/cells.csv.gz",
                        "data/Human_Brain_Glioblastoma/cell_boundaries.parquet",
                        "data/Human_Brain_Glioblastoma/nucleus_boundaries.parquet")

#plotColData(glio.sfe, x="cell_area", y="nucleus_area", bins = 100)
#plotSpatialFeature(glio.sfe, "nGenes", colGeometryName = "cellSeg")

#plotSpatialFeature(glio.sfe, "cell_area", colGeometryName = "cellSeg")
#plotSpatialFeature(glio.sfe, "nucleus_area", colGeometryName = "nucSeg")

#cols_use <- names(colData(glio.sfe))[str_detect(names(colData(glio.sfe)), "_percent$")]
#cols_use <- cols_use[-4]

# for (n in cols_use) {
#     glio.sfe <- get_neg_ctrl_outliers(n, glio.sfe)
# }
# 
# inds_keep <- glio.sfe$nCounts > 0 & 
#     !glio.sfe$is_blank_outlier & !glio.sfe$is_negCodeword_outlier & !glio.sfe$is_negProbe_outlier
# glio.sfe <- glio.sfe[,inds_keep] #201 cells were discarded

rowData(glio.sfe)$means <- rowMeans(counts(glio.sfe))
rowData(glio.sfe)$vars <- rowVars(counts(glio.sfe))

plotRowData(glio.sfe, x="means", y="vars", color_by = "is_neg") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    scale_x_log10() + scale_y_log10() +
    annotation_logticks() +
    coord_equal() +
    labs(color = "Negative control")

p2 <- plotMeanVar(rowData(glio.sfe)$means, rowData(glio.sfe)$vars, "Xenium Human Brain Glioblastoma")


# Alzheimers

alz.sfe <- readXenium("data/Human_Brain_Alzheimers/cell_feature_matrix.h5", 
                       "data/Human_Brain_Alzheimers/cells.csv.gz",
                       "data/Human_Brain_Alzheimers/cell_boundaries.parquet",
                       "data/Human_Brain_Alzheimers/nucleus_boundaries.parquet")

# plotColData(alz.sfe, x="cell_area", y="nucleus_area", bins = 100)
# plotSpatialFeature(alz.sfe, "nGenes", colGeometryName = "cellSeg")
# 
# plotSpatialFeature(alz.sfe, "cell_area", colGeometryName = "cellSeg")
# plotSpatialFeature(alz.sfe, "nucleus_area", colGeometryName = "nucSeg")
# 
# cols_use <- names(colData(alz.sfe))[str_detect(names(colData(alz.sfe)), "_percent$")]
# #cols_use <- cols_use[-4]
# 
# for (n in cols_use) {
#     alz.sfe <- get_neg_ctrl_outliers(n, alz.sfe)
# }
# 
# inds_keep <- alz.sfe$nCounts > 0 & 
#     !alz.sfe$is_blank_outlier & !alz.sfe$is_negCodeword_outlier & !alz.sfe$is_negProbe_outlier
# alz.sfe <- alz.sfe[,inds_keep] #109 cells were discarded

rowData(alz.sfe)$means <- rowMeans(counts(alz.sfe))
rowData(alz.sfe)$vars <- rowVars(counts(alz.sfe))

p3 <- plotMeanVar(rowData(alz.sfe)$means, rowData(alz.sfe)$vars, "Xenium Human Brain Alzheimer's")

# plotRowData(alz.sfe, x="means", y="vars", color_by = "is_neg") +
#     geom_abline(slope = 1, intercept = 0, color = "red") +
#     scale_x_log10() + scale_y_log10() +
#     annotation_logticks() +
#     coord_equal() +
#     labs(color = "Negative control")

p <- p1 + p2 + p3
ggsave(here("./plots/mean_var_plots/xenium_human_brain_all.png"), p, width=20, height=10)



mousebrain_vizgen <- readRDS("data/Vizgen/mousebrain-slide1-rep1-meanvar-df.RDS")
mean_emp <- mousebrain_vizgen$rowmeans
var_emp <- mousebrain_vizgen$rowvars
r1 <- plotMeanVar(mean_emp, var_emp, "Vizgen Mouse Brain Slice 1 Rep 1")

vizgen.s1.r2 <- readRDS("data/Vizgen/mousebrain-slice1-rep2-meanvar-df.RDS")
mean_emp <- vizgen.s1.r2$rowmeans
var_emp <- vizgen.s1.r2$rowvars
r2 <- plotMeanVar(mean_emp, var_emp, "Vizgen Mouse Brain Slice 1 Rep 2")



vizgen.s1.r3 <- readRDS("data/Vizgen/mousebrain-slice1-rep3-meanvar-df.RDS")
mean_emp <- vizgen.s1.r3$rowmeans
var_emp <- vizgen.s1.r3$rowvars
r3 <- plotMeanVar(mean_emp, var_emp, "Vizgen Mouse Brain Slice 1 Rep 3")

ggsave(here("./plots/mean_var_plots/vizgen_mouse_brain_slice1.png"), r1+r2+r3, width=20, height=10)

cosmx <- HeNSCLCData()
cosmx.mean_emp <- rowMeans(counts(cosmx))
cosmx.var_emp <- rowVars(counts(cosmx))
p.cm <- plotMeanVar(cosmx.mean_emp, cosmx.var_emp, "CosMx Human Liver")

ggsave(here("./plots/mean_var_plots/cosmx_human_liver.png"), p.cm)
