suppressPackageStartupMessages({
    library(Voyager)
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

######################################
# Make some QC plots for each region #
######################################


##### Slide 5434 #####
data_dir <- "slide-5434"

# Read in the three regions
br6471_p <- readRDS(here("processed-data", "cindy", data_dir, "Br6471_Post_SFE.RDS"))
br6522_p <- readRDS(here("processed-data", "cindy", data_dir, "Br6522_Post_SFE.RDS"))
br8667_m <- readRDS(here("processed-data", "cindy", data_dir, "Br8667_Mid_SFE.RDS"))

sfe.all <- readRDS(here("processed-data", "cindy", data_dir, "xenium-0005434-SFE.RDS"))

sfe.all <- sfe.all[,which(sfe.all$region_id != "Sample")]

br6471_p$total_gene_counts <- br6471_p$total_counts - br6471_p$subsets_any_neg_sum
br6522_p$total_gene_counts <- br6522_p$total_counts - br6522_p$subsets_any_neg_sum
br8667_m $total_gene_counts <- br8667_m$total_counts - br8667_m$subsets_any_neg_sum

h1 <- plotColData(sfe.all, x="region_id", y="nucleus_area", 
                  colour_by="region_id",
                  show_median=TRUE,
                  point_alpha=0.3)+
    ggtitle("Nucleus Area for Slide 5434")

h2 <- plotColData(sfe.all, x="region_id", y="subsets_any_neg_sum",
                  colour_by="region_id",
                  show_median=TRUE,
                  point_alpha=0.3)+
    ggtitle("Number of negative control counts for slide 5434")

p1 <- plotColData(br6471_p, x="nucleus_area", y="total_counts", bin=100)+
    ggtitle(unique(br6471_p$region_id)[[1]])

p2 <- plotColData(br6522_p, x="nucleus_area", y="total_counts", bin=100)+
    ggtitle(unique(br6522_p$region_id)[[1]])

p3 <- plotColData(br8667_m, x="nucleus_area", y="total_counts", bin=100)+
    ggtitle(unique(br8667_m$region_id)[[1]])



br6471_p$nucleus_area_200 <- br6471_p$nucleus_area >= 200
br6471_p$nucleus_area_50 <- br6471_p$nucleus_area >= 50

br6522_p$nucleus_area_200 <- br6522_p$nucleus_area >= 200
br8667_m$nucleus_area_200 <- br8667_m$nucleus_area >= 200

# n1 <- plotSpatialFeature(br6471_p, "nucleus_area_100", colGeometryName = "cellSeg")+
#     ggtitle(unique(br6471_p$region_id))
# 
# n2 <- plotSpatialFeature(br6522_p, "nucleus_area_100", colGeometryName = "cellSeg")+
#     ggtitle(unique(br6522_p$region_id))
# 
# n3 <- plotSpatialFeature(br8667_m, "nucleus_area_100", colGeometryName = "cellSeg")+
#     ggtitle(unique(br8667_m$region_id))

n1 <- make_escheR(br6471_p) |>
    add_fill("nucleus_area_200", point_size = 0.5)+
    scale_fill_manual(values=c("grey", "red"))+
    ggtitle(unique(br6471_p$region_id))

n2 <- make_escheR(br6522_p) |>
    add_fill("nucleus_area_200", point_size = 0.5)+
    scale_fill_manual(values=c("grey", "red"))+
    ggtitle(unique(br6522_p$region_id))

n3 <- make_escheR(br8667_m) |>
    add_fill("nucleus_area_200", point_size = 0.5)+
    scale_fill_manual(values=c("grey", "red"))+
    ggtitle(unique(br8667_m$region_id))


np1 <- plotColData(br6471_p, x="nucleus_area", y="subsets_any_neg_percent", bin=100)+
    ggtitle(unique(br6471_p$region_id)[[1]])+
    geom_smooth(method="lm")

np2 <- plotColData(br6522_p, x="nucleus_area", y="subsets_any_neg_percent", bin=100)+
    ggtitle(unique(br6522_p$region_id)[[1]])+
    geom_smooth(method="lm")

np3 <- plotColData(br8667_m, x="nucleus_area", y="subsets_any_neg_percent", bin=100)+
    ggtitle(unique(br8667_m$region_id)[[1]])+
    geom_smooth(method="lm")

# br6471_p <- addPerCellQCMetrics(br6471_p, subsets=list(
#     unassigned_codeword="unassigned_codeword_counts"))

c1 <- plotColData(br6471_p, x="nucleus_area", y="unassigned_codeword_counts", bin=100)+
    ggtitle(unique(br6471_p$region_id)[[1]])+
    geom_smooth(method="lm")

c2 <- plotColData(br6522_p, x="nucleus_area", y="unassigned_codeword_counts", bin=100)+
    ggtitle(unique(br6522_p$region_id)[[1]])+
    geom_smooth(method="lm")

c3 <- plotColData(br8667_m, x="nucleus_area", y="unassigned_codeword_counts", bin=100)+
    ggtitle(unique(br8667_m$region_id)[[1]])+
    geom_smooth(method="lm")

nt1 <- plotColData(br6471_p, x="total_gene_counts", y="subsets_any_neg_percent", bin=100)+
    ggtitle(unique(br6471_p$region_id)[[1]])+
    geom_smooth(method="lm")

nt2 <- plotColData(br6522_p, x="total_gene_counts", y="subsets_any_neg_percent", bin=100)+
    ggtitle(unique(br6522_p$region_id)[[1]])+
    geom_smooth(method="lm")

nt3 <- plotColData(br6522_p, x="total_gene_counts", y="subsets_any_neg_percent", bin=100)+
    ggtitle(unique(br6522_p$region_id)[[1]])+
    geom_smooth(method="lm")

pdf(here("plots", "02_qualityControl", "slide5434_QCPlots.pdf"), height=5, width=12)
print(h1)
print(h2)
print(p1+p2+p3)
print(n1+n2+n3)
print(np1+np2+np3)
print(c1+c2+c3)
print(nt1 + nt2 + nt3)
dev.off()


##### Slide 5548 #####



data_dir <- "slide-5548"

# Read in the three regions
br6471_p <- readRDS(here("processed-data", "cindy", data_dir, "Br6471_Post_SFE.RDS"))
br2743_m <- readRDS(here("processed-data", "cindy", data_dir, "Br2743_Mid_SFE.RDS"))
br8667_m <- readRDS(here("processed-data", "cindy", data_dir, "Br8667_Mid_SFE.RDS"))

sfe.all <- readRDS(here("processed-data", "cindy", data_dir, "xenium-0005548-SFE.RDS"))

br6471_p$total_gene_counts <- br6471_p$total_counts - br6471_p$subsets_any_neg_sum
br2743_m$total_gene_counts <- br2743_m$total_counts - br2743_m$subsets_any_neg_sum
br8667_m $total_gene_counts <- br8667_m$total_counts - br8667_m$subsets_any_neg_sum

h1 <- plotColData(sfe.all, x="region_id", y="nucleus_area", 
                  colour_by="region_id",
                  show_median=TRUE,
                  point_alpha=0.3)+
    ggtitle("Nucleus Area for Slide 5548")

h2 <- plotColData(sfe.all, x="region_id", y="subsets_any_neg_sum",
                  colour_by="region_id",
                  show_median=TRUE,
                  point_alpha=0.3)+
    ggtitle("Number of negative control counts for slide 5548")

p1 <- plotColData(br6471_p, x="nucleus_area", y="total_counts", bin=100)+
    ggtitle(unique(br6471_p$region_id)[[1]])

p2 <- plotColData(br2743_m, x="nucleus_area", y="total_counts", bin=100)+
    ggtitle(unique(br2743_m$region_id)[[1]])

p3 <- plotColData(br8667_m, x="nucleus_area", y="total_counts", bin=100)+
    ggtitle(unique(br8667_m$region_id)[[1]])


br6471_p$nucleus_area_200 <- br6471_p$nucleus_area >= 200
br6471_p$nucleus_area_50 <- br6471_p$nucleus_area >= 50

br2743_m$nucleus_area_200 <- br2743_m$nucleus_area >= 200
br8667_m$nucleus_area_200 <- br8667_m$nucleus_area >= 200

# n1 <- plotSpatialFeature(br6471_p, "nucleus_area_100", colGeometryName = "cellSeg")+
#     ggtitle(unique(br6471_p$region_id))
# 
# n2 <- plotSpatialFeature(br2743_m, "nucleus_area_100", colGeometryName = "cellSeg")+
#     ggtitle(unique(br2743_m$region_id))
# 
# n3 <- plotSpatialFeature(br8667_m, "nucleus_area_100", colGeometryName = "cellSeg")+
#     ggtitle(unique(br8667_m$region_id))

n1 <- make_escheR(br6471_p) |>
    add_fill("nucleus_area_200")+
    scale_fill_manual(values=c("grey", "red"))+
    ggtitle(unique(br6471_p$region_id))

n2 <- make_escheR(br2743_m) |>
    add_fill("nucleus_area_200")+
    scale_fill_manual(values=c("grey", "red"))+
    ggtitle(unique(br2743_m$region_id))

n3 <- make_escheR(br8667_m) |>
    add_fill("nucleus_area_200")+
    scale_fill_manual(values=c("grey", "red"))+
    ggtitle(unique(br8667_m$region_id))




np1 <- plotColData(br6471_p, x="nucleus_area", y="subsets_any_neg_percent", bin=100)+
    ggtitle(unique(br6471_p$region_id)[[1]])

np2 <- plotColData(br2743_m, x="nucleus_area", y="subsets_any_neg_percent", bin=100)+
    ggtitle(unique(br2743_m$region_id)[[1]])

np3 <- plotColData(br8667_m, x="nucleus_area", y="subsets_any_neg_percent", bin=100)+
    ggtitle(unique(br8667_m$region_id)[[1]])

#br6471_p <- addPerCellQCMetrics(br6471_p, subsets=list(
    #unassigned_codeword="unassigned_codeword_counts"))

c1 <- plotColData(br6471_p, x="nucleus_area", y="unassigned_codeword_counts", bin=100)+
    ggtitle(unique(br6471_p$region_id)[[1]])+
    geom_smooth(method="lm")

c2 <- plotColData(br2743_m, x="nucleus_area", y="unassigned_codeword_counts", bin=100)+
    ggtitle(unique(br2743_m$region_id)[[1]])+
    geom_smooth(method="lm")

c3 <- plotColData(br8667_m, x="nucleus_area", y="unassigned_codeword_counts", bin=100)+
    ggtitle(unique(br8667_m$region_id)[[1]])+
    geom_smooth(method="lm")

nt1 <- plotColData(br6471_p,  x="total_gene_counts", y="subsets_any_neg_percent", bin=100)+
    ggtitle(unique(br6471_p$region_id)[[1]])+
    geom_smooth(method="lm")

nt2 <- plotColData(br2743_m,  x="total_gene_counts", y="subsets_any_neg_percent", bin=100)+
    ggtitle(unique(br2743_m$region_id)[[1]])+
    geom_smooth(method="lm")

nt3 <- plotColData(br8667_m, x="total_gene_counts", y="subsets_any_neg_percent", bin = 100)+
    ggtitle(unique(br8667_m$region_id)[[1]])+
    geom_smooth(method="lm")


pdf(here("plots", "02_qualityControl", "slide5548_QCPlots.pdf"), height=5, width=12)
print(h1)
print(h2)
print(p1+p2+p3)
print(n1+n2+n3)
print(np1+np2+np3)
print(c1+c2+c3)
print(nt1 + nt2 + nt3)
dev.off()
