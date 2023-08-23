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
source(here("./code/xenium-helpers.R"))

sfe.healthy <- readXenium("data/Human_Brain/cell_feature_matrix.h5", "data/Human_Brain/cells.csv.gz",
                          "data/Human_Brain/cell_boundaries.parquet",
                          "data/Human_Brain/nucleus_boundaries.parquet")


cols_use <- names(colData(sfe.healthy))[str_detect(names(colData(sfe.healthy)), "_percent$")]
cols_use <- cols_use[-4]

for (n in cols_use) {
    sfe.healthy <- get_neg_ctrl_outliers(n, sfe.healthy)
}

inds_keep <- sfe.healthy$nCounts > 0 & sfe.healthy$nucleus_area < 200 &
    !sfe.healthy$is_blank_outlier & !sfe.healthy$is_negCodeword_outlier & !sfe.healthy$is_negProbe_outlier
sfe.healthy <- sfe.healthy[,inds_keep] # only 74 cells were discarded

rowData(sfe.healthy)$means <- rowMeans(counts(sfe.healthy))
rowData(sfe.healthy)$vars <- rowVars(counts(sfe.healthy))

plotRowData(sfe.healthy, x="means", y="vars", color_by = "is_neg") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    scale_x_log10() + scale_y_log10() +
    annotation_logticks() +
    coord_equal() +
    labs(color = "Negative control")

mean_emp <- rowData(sfe.healthy)$means
var_emp <- rowData(sfe.healthy)$vars

# Negative Binomial
# Estimate the overall size/dispersion parameters
model = lm(var_emp ~ 1*mean_emp + I(mean_emp^2) + 0, tibble(mean_emp, var_emp))
phi = 1/coef(model)["I(mean_emp^2)"]


mean_var_tb = tibble(mean_emp = mean_emp,
                     var_emp = var_emp,
                     nbinomial = mean_emp + mean_emp^2 * 1/phi) %>%
    tidyr::pivot_longer(cols = -mean_emp, names_to = "model", values_to = "var_value")


p = mean_var_tb %>%
    filter(model %in% c("var_emp")) %>%
    ggplot(aes(x = mean_emp, y = var_value)) + 
    geom_point(alpha = 0.3) + 
    geom_line(data = mean_var_tb %>% filter(model %in% c("nbinomial")),
              aes(x = mean_emp, y = var_value, color = model)) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Log of mean expression",
         y = "Log of variance") +
    geom_abline(slope = 1, intercept = 0, color = "red")+
    theme_bw()
