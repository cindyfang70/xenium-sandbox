library(escheR)
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(dplyr)
sfe <- readRDS("processed-data/cindy/slide-5434/slide5434-all-samples-spe-with-banksy.RDS")


sfe <- sfe[,which(sfe$region_id == "Br6471_Post_5434")]
sfe$counts_MOBP <- counts(sfe)[which(rownames(sfe)=="MOBP"),]
p <- make_escheR(sfe, y_reverse=FALSE) %>%
    add_fill("counts_MOBP", point_size=4.5) %>%
    add_ground("clust_M1_lam0.9_k50_res0.6", show.legend=FALSE)+
    scale_fill_gradient(low="white", high="black")

pdf("escheR_example.pdf", width=25, height=25)
print(p)
dev.off()

