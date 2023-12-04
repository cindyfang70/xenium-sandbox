library(brglm2)


k <- 20
model <- readRDS(here("processed-data", "cindy", "NMF", sprintf("visium-nmf-model-k%s", k)))
factors <- t(model$h)
colnames(factors) <- paste0("NMF", 1:dim(factors)[[2]])

# get the manually annotated Visium data
ehub <- ExperimentHub::ExperimentHub()
vis_anno <- fetch_data(type = "spe", eh = ehub)
